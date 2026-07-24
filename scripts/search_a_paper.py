import html
import json
import os
from pathlib import Path
import re
import sys
import urllib.parse
import urllib.request

ROOT_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HEADERS = {"User-Agent": "doi-metadata-script/0.1"}


def get_crossref_metadata(doi: str) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)["message"]


def clean_title(title: str) -> str:
    title = re.sub(r"<[^>]+>", " ", title)
    title = re.sub(r"\s+", " ", html.unescape(title)).strip()
    return re.sub(r"\b([A-Z])\s+(\d+)\s+([a-z])\b", r"\2\1\3", title)


def print_head(doi: str, metadata: dict) -> None:
    date_parts = (metadata.get("issued") or {}).get("date-parts") or [[]]
    journal = (metadata.get("short-container-title") or [""])[0]
    volume = metadata.get("volume") or ""
    number = metadata.get("article-number") or ""
    page = metadata.get("page") or ""
    year = (date_parts[0] or [""])[0]
    print(f"doi: {metadata.get('DOI') or doi}")
    print(f"title: {clean_title((metadata.get('title') or [''])[0])}")
    print(f"Journal: {journal}, Vol. {volume}, Nol. {number}, P. {page}, Y. {year}")


def get_md_paths(raw_dir: str | Path = os.path.join(ROOT_PATH, "raw")) -> list[Path]:
    paths = Path(raw_dir).glob("*.md")
    return sorted(paths)


def get_json_paths(raw_dir: str | Path = os.path.join(ROOT_PATH, "raw")) -> list[Path]:
    paths = Path(raw_dir).glob("*.json")
    return sorted(paths)


def find_recorded(doi: str, json_data: dict) -> str | None:
    if (json_data.get("DOI") or "").strip().lower() != doi: return
    filename = json_data["filename"]
    print(f"recorded: {filename}")
    for path in sorted(Path(ROOT_PATH, "raw").glob(f"{filename}.*")):
        if path.is_file(): print(f"samename: {path.name}")
    return filename


def find_json_reference(doi: str, json_data: dict) -> None:
    filename = json_data["filename"]
    for ref in json_data.get("reference") or []:
        ref_doi = (ref.get("DOI") or "").strip().lower()
        if ref_doi == doi:
            print(f"json reference: {filename}: {ref.get('key') or ref_doi}")


def find_md_secondary(doi: str, md_path: Path) -> None:
    lines = md_path.read_text(encoding="utf-8").splitlines()
    start, end = lines.index("### Secondary Citations") + 1, lines.index("### Gaps")
    for line_number, line in enumerate(lines[start:end], start + 1):
        match = re.match(r"\[doi\]:\s*(.+)$", line.strip(), re.I)
        if not match: continue
        line_dois = [item.strip().lower() for item in match.group(1).split(";")]
        if doi in line_dois: print(f"md secondary: {md_path.name}: {line_number}: {line.strip()}")


def main(doi: str) -> None:
    doi = doi.strip().lower().removeprefix("https://doi.org/").removeprefix("doi:").strip()
    metadata = get_crossref_metadata(doi)
    print_head(doi, metadata)

    for json_path in get_json_paths():
        json_data = json.loads(json_path.read_text(encoding="utf-8"))
        find_recorded(doi, json_data)
        find_json_reference(doi, json_data)

    for md_path in get_md_paths():
        find_md_secondary(doi, md_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python scripts/search_a_doi.py DOI")
    main(sys.argv[1])
