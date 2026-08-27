#!/usr/bin/env python3
"""Export the HTML slide deck to a validated 16:9 PDF."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from pathlib import Path


PRINT_STYLE = r"""
  <style id="pdf-export-style">
    @page {
      size: 1600px 900px;
      margin: 0;
    }

    @media print {
      html,
      body {
        width: 1600px;
        margin: 0 !important;
        padding: 0 !important;
        background: #fff !important;
      }

      .deck {
        display: block !important;
        width: 1600px !important;
        margin: 0 !important;
        gap: 0 !important;
      }

      .slide {
        width: 1600px !important;
        height: 900px !important;
        margin: 0 !important;
        break-after: page;
        page-break-after: always;
        box-shadow: none !important;
      }

      .slide:last-child {
        break-after: auto;
        page-break-after: auto;
      }
    }
  </style>
"""


def find_executable(candidates: list[str | Path | None]) -> Path | None:
    for candidate in candidates:
        if not candidate:
            continue
        resolved = shutil.which(str(candidate))
        if resolved:
            return Path(resolved).resolve()
        path = Path(candidate).expanduser()
        if path.is_file():
            return path.resolve()
    return None


def default_poppler_binary(name: str) -> Path:
    return (
        Path.home()
        / ".cache"
        / "codex-runtimes"
        / "codex-primary-runtime"
        / "dependencies"
        / "native"
        / "poppler"
        / "Library"
        / "bin"
        / f"{name}.exe"
    )


def count_slides(source: str) -> int:
    pattern = re.compile(
        r"<section\b[^>]*\bclass\s*=\s*(['\"])[^'\"]*\bslide\b[^'\"]*\1",
        re.IGNORECASE,
    )
    return len(pattern.findall(source))


def make_export_html(source: str) -> str:
    if "</head>" not in source:
        raise ValueError("The HTML file does not contain </head>.")
    return source.replace("</head>", f"{PRINT_STYLE}\n</head>", 1)


def browser_arguments(
    profile: Path,
    candidate_pdf: Path,
    temporary_html: Path,
) -> list[str]:
    return [
        "--headless",
        "--no-sandbox",
        "--disable-gpu",
        "--disable-gpu-compositing",
        "--disable-features=UseSkiaRenderer,Vulkan",
        "--use-angle=swiftshader",
        "--disable-extensions",
        "--allow-file-access-from-files",
        "--no-pdf-header-footer",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=15000",
        f"--user-data-dir={profile}",
        f"--print-to-pdf={candidate_pdf}",
        temporary_html.as_uri(),
    ]


def validate_page_count(pdf: Path, expected: int, pdfinfo: Path | None) -> None:
    if not pdfinfo:
        print("Warning: pdfinfo was not found; page-count validation was skipped.")
        return
    result = subprocess.run(
        [str(pdfinfo), str(pdf)],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    match = re.search(r"^Pages:\s+(\d+)", result.stdout, re.MULTILINE)
    if not match:
        raise RuntimeError("pdfinfo did not report a page count.")
    actual = int(match.group(1))
    if actual != expected:
        raise RuntimeError(
            f"Page-count mismatch: HTML has {expected} slides, PDF has {actual} pages."
        )


def validate_mathjax(pdf: Path, pdftotext: Path | None, workdir: Path) -> None:
    if not pdftotext:
        print("Warning: pdftotext was not found; TeX validation was skipped.")
        return
    text_file = workdir / "pdf-text.txt"
    subprocess.run(
        [str(pdftotext), "-enc", "UTF-8", str(pdf), str(text_file)],
        check=True,
        capture_output=True,
    )
    text = text_file.read_text(encoding="utf-8", errors="replace")
    if re.search(r"\\\(|\\\[|\\alpha|\\Psi|\\begin\{", text):
        raise RuntimeError(
            "Unrendered TeX was detected. Check the MathJax network connection and export again."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export an HTML slide deck to a validated 16:9 PDF."
    )
    parser.add_argument(
        "html",
        type=Path,
        help="Input HTML file path.",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Output PDF file path.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    html_path = args.html.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    if not html_path.is_file():
        raise FileNotFoundError(f"HTML file not found: {html_path}")
    edge = find_executable(
        [
            "msedge.exe",
            Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"),
            Path(r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"),
            "chrome.exe",
            Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
        ]
    )
    if not edge:
        raise FileNotFoundError("Microsoft Edge or Google Chrome was not found.")

    pdfinfo = find_executable(["pdfinfo.exe", default_poppler_binary("pdfinfo")])
    pdftotext = find_executable(
        ["pdftotext.exe", default_poppler_binary("pdftotext")]
    )

    source = html_path.read_text(encoding="utf-8")
    slide_count = count_slides(source)
    if slide_count == 0:
        raise ValueError('No <section class="slide ..."> elements were found.')

    temporary_html = html_path.with_name(f"_pdf-export-{uuid.uuid4().hex}.html")
    export_html = make_export_html(source)

    try:
        temporary_html.write_text(export_html, encoding="utf-8", newline="")
        with tempfile.TemporaryDirectory(prefix="slides-pdf-") as temporary_root:
            workdir = Path(temporary_root)
            profile = workdir / "browser-profile"
            profile.mkdir()
            candidate_pdf = workdir / "candidate.pdf"

            process = subprocess.run(
                [str(edge), *browser_arguments(profile, candidate_pdf, temporary_html)],
                check=False,
                creationflags=(
                    subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
                ),
            )
            if process.returncode != 0:
                raise RuntimeError(
                    f"Browser export failed with exit code {process.returncode}."
                )
            if not candidate_pdf.is_file():
                raise RuntimeError("Browser finished without creating the PDF.")

            validate_page_count(candidate_pdf, slide_count, pdfinfo)
            validate_mathjax(candidate_pdf, pdftotext, workdir)

            output_path.parent.mkdir(parents=True, exist_ok=True)
            os.replace(candidate_pdf, output_path)
    finally:
        temporary_html.unlink(missing_ok=True)

    size_mib = output_path.stat().st_size / (1024 * 1024)
    print(f"Created: {output_path}")
    print(f"Slides:  {slide_count}")
    print(f"Size:    {size_mib:.2f} MiB")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)

