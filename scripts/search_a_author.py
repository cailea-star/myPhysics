import json
import os
from pathlib import Path
import re
import sys

ROOT_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def get_json_paths(raw_dir: str | Path = os.path.join(ROOT_PATH, "raw")) -> list[Path]:
    paths = Path(raw_dir).glob("*.json")
    return sorted(paths)


def find_author(given: str, family: str, fulljson: str) -> list[str]:
    target = re.sub(r"\W+", "", f"{given} {family}").lower()
    data = json.loads(fulljson)
    hits = []
    for author in data.get("author") or []:
        name = " ".join(filter(None, [author.get("given"), author.get("family")]))
        if re.sub(r"\W+", "", name).lower() == target:
            hits.append(name)
    return hits


def find_corresponding(given: str, family: str, fulljson: str) -> list[str]:
    target = re.sub(r"\W+", "", f"{given} {family}").lower()
    data = json.loads(fulljson)
    hits = []
    for author in data.get("author-corresponding-openalex") or []:
        name = " ".join(filter(None, [author.get("given"), author.get("family")]))
        if re.sub(r"\W+", "", name).lower() == target:
            hits.append(name)
    return hits


def main(given: str, family: str) -> None:
    for path in get_json_paths():
        fulljson = path.read_text(encoding="utf-8")
        for hit in find_author(given, family, fulljson):
            print(f"author: {path}: {hit}")
        for hit in find_corresponding(given, family, fulljson):
            print(f"corresponding: {path}: {hit}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/search_a_author.py GIVEN FAMILY")
    main(sys.argv[1], sys.argv[2])
