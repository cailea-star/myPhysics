import json
import os
from pathlib import Path
import re
import sys

ROOT_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def get_md_paths(raw_dir: str | Path = os.path.join(ROOT_PATH, "raw")) -> list[Path]:
    paths = Path(raw_dir).glob("*.md")
    return sorted(paths)


def get_json_paths(raw_dir: str | Path = os.path.join(ROOT_PATH, "raw")) -> list[Path]:
    paths = Path(raw_dir).glob("*.json")
    return sorted(paths)


def find_recorded(doi: str, fulljson: str) -> bool:
    data = json.loads(fulljson)
    return (data.get("DOI") or "").strip().lower() == doi


def find_json_reference(doi: str, fulljson: str) -> list[str]:
    data = json.loads(fulljson)
    hits = []
    for ref in data.get("reference") or []:
        ref_doi = (ref.get("DOI") or "").strip().lower()
        if ref_doi == doi:
            hits.append(ref.get("key") or ref_doi)
    return hits


def find_md_secondary(doi: str, fullmd: str) -> list[str]:
    hits = []
    for line_number, line in enumerate(fullmd.splitlines(), 1):
        match = re.match(r"\[doi\]:\s*(.+)$", line.strip(), re.I)
        if not match: continue
        line_dois = [item.strip().lower() for item in match.group(1).split(";")]
        if doi in line_dois:
            hits.append(f"{line_number}: {line.strip()}")
    return hits


def main(doi: str) -> None:
    doi = doi.strip().lower().removeprefix("https://doi.org/").removeprefix("doi:").strip()
    for path in get_json_paths():
        fulljson = path.read_text(encoding="utf-8")
        if find_recorded(doi, fulljson):
            print(f"recorded: {path}")
        for hit in find_json_reference(doi, fulljson):
            print(f"json reference: {path}: {hit}")

    for path in get_md_paths():
        fullmd = path.read_text(encoding="utf-8")
        for hit in find_md_secondary(doi, fullmd):
            print(f"md secondary: {path}: {hit}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python scripts/search_a_doi.py DOI")
    main(sys.argv[1])
