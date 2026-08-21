import html
import json
from pathlib import Path
import re
import sys
import urllib.parse
import urllib.request

PAPERS_PATH = Path(__file__).resolve().parent.parent
RAW_PATH = PAPERS_PATH.joinpath("raw")
HEADERS = {"User-Agent": "doi-metadata-script/0.1"}


def load_crossref_metadata(doi: str) -> dict:
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    request = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)["message"]


def clean_title(title_str: str) -> str:
    title_str = re.sub(r"<[^>]+>", " ", title_str)
    title_str = re.sub(r"\s+", " ", html.unescape(title_str)).strip()
    return re.sub(r"\b([A-Z])\s+(\d+)\s+([a-z])\b", r"\2\1\3", title_str)


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


def get_raw_md_paths(raw_dir: str | Path = RAW_PATH) -> list[Path]:
    paths = Path(raw_dir).glob("*.md")
    return sorted(paths)


def get_raw_json_paths(raw_dir: str | Path = RAW_PATH) -> list[Path]:
    paths = Path(raw_dir).glob("*.json")
    return sorted(paths)


def find_recorded(doi: str, raw_json_data: dict) -> str | None:
    if (raw_json_data.get("DOI") or "").strip().lower() != doi: return
    raw_stem = raw_json_data["stem"]
    print(f"recorded: {raw_stem}")
    for path in sorted(RAW_PATH.glob(f"{raw_stem}.*")):
        if path.is_file(): print(f"samename: {path.name}")
    return raw_stem


def split_doi_line(line: str) -> list[str]:
    match = re.match(r"\[doi\]:\s*(.+)$", line.strip(), re.I)
    if not match: return []
    dois = []
    for item in match.group(1).split(";"):
        dois.append(re.sub(r"^\[[^\]]+\]\s*", "", item.strip()).lower())
    return dois


def find_md_secondary(doi: str, raw_md_path: Path) -> None:
    lines = raw_md_path.read_text(encoding="utf-8").splitlines()
    start, end = lines.index("### Secondary Citations") + 1, lines.index("### Gaps")
    for line_number, line in enumerate(lines[start:end], start + 1):
        if doi in split_doi_line(line):
            print(f"md secondary: {raw_md_path.name}: {line_number}: {line.strip()}")


def main(doi: str) -> None:
    doi = doi.strip().lower().removeprefix("https://doi.org/").removeprefix("doi:").strip()
    metadata = load_crossref_metadata(doi)
    print_head(doi, metadata)

    for raw_json_path in get_raw_json_paths():
        raw_json_data = json.loads(raw_json_path.read_text(encoding="utf-8"))
        find_recorded(doi, raw_json_data)

    for raw_md_path in get_raw_md_paths():
        find_md_secondary(doi, raw_md_path)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python papers/scripts/search_a_paper.py DOI")
    main(sys.argv[1])
