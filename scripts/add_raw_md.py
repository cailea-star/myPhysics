import json
import sys
import os
from pathlib import Path


def replace_source_json(text, json_data):
    return text.replace("[json_filename.json]", json_data["_json_filename"])


def replace_doi(text, json_data):
    return text.replace("[doi_number]", json_data.get("DOI") or "")


def replace_first_author(text, json_data):
    author = (json_data.get("author") or [{}])[0]
    return text.replace("[given_name]", author.get("given") or "").replace("[family_name]", author.get("family") or "")


def replace_affiliations(text, json_data):
    first_author = (json_data.get("author") or [{}])[0]
    corresponding = (json_data.get("author-corresponding-openalex") or [{}])[0]
    affiliations = corresponding.get("raw_affiliation_strings") or first_author.get("affiliation") or []
    names = [a.get("name") if isinstance(a, dict) else str(a) for a in affiliations]
    return text.replace("[affiliations]", ", ".join(filter(None, names)))


def replace_title(text, json_data):
    return text.replace("[title_name]", (json_data.get("title") or [""])[0])


def replace_abstract(text, json_data):
    return text.replace("[abstract]", json_data.get("abstract") or "")


def main(json_path):
    json_path = Path(json_path)
    md_path = json_path.with_suffix(".md")
    if json_path.suffix.lower() != ".json":
        raise SystemExit("input must be a .json file")
    if md_path.exists():
        raise SystemExit(f"{md_path} already exists")

    json_data = json.loads(json_path.read_text(encoding="utf-8"))
    json_data["_json_filename"] = json_path.name

    text = Path(os.path.join(os.path.dirname(__file__), "add_raw_md.md")).read_text(encoding="utf-8")
    text = replace_source_json(text, json_data)
    text = replace_doi(text, json_data)
    text = replace_first_author(text, json_data)
    text = replace_affiliations(text, json_data)
    text = replace_title(text, json_data)
    text = replace_abstract(text, json_data)
    md_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python scripts/add_raw_md.py raw/paper.json")
    main(sys.argv[1])
