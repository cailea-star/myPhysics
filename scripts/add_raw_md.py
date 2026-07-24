import json
import sys
import re
from pathlib import Path


def remove_html_comments(raw_template_str):
    return re.sub(r"<!--.*?-->", "", raw_template_str, flags=re.DOTALL)


def replace_source_json(raw_template_str, raw_json_data):
    return raw_template_str.replace("[this_source_json]", raw_json_data["_json_filename"])


def replace_doi(raw_template_str, raw_json_data):
    return raw_template_str.replace("[this_doi]", raw_json_data.get("DOI") or "")


def replace_first_author(raw_template_str, raw_json_data):
    author = (raw_json_data.get("author") or [{}])[0]
    return raw_template_str.replace("[this_first_name]", author.get("given") or "").replace("[this_family_name]", author.get("family") or "")


def replace_affiliations(raw_template_str, raw_json_data):
    first_author = (raw_json_data.get("author") or [{}])[0]
    corresponding = (raw_json_data.get("author-corresponding-openalex") or [{}])[0]
    affiliations = corresponding.get("raw_affiliation_strings") or first_author.get("affiliation") or []
    names = [a.get("name") if isinstance(a, dict) else str(a) for a in affiliations]
    return raw_template_str.replace("[this_affiliations]", ", ".join(filter(None, names)))


def replace_title(raw_template_str, raw_json_data):
    return raw_template_str.replace("[title_name]", (raw_json_data.get("title") or [""])[0])


def replace_abstract(raw_template_str, raw_json_data):
    return raw_template_str.replace("[abstract]", raw_json_data.get("abstract") or "")


def main(raw_json_path):
    raw_json_path = Path(raw_json_path)
    raw_md_path = raw_json_path.with_suffix(".md")
    if raw_json_path.suffix.lower() != ".json": raise SystemExit("input must be a .json file")
    if raw_md_path.exists(): raise SystemExit(f"{raw_md_path} already exists")

    raw_json_data = json.loads(raw_json_path.read_text(encoding="utf-8"))
    raw_json_data["_json_filename"] = raw_json_path.name

    raw_template_path = Path(__file__).resolve().with_name("add_raw_md.md")
    raw_template_str = raw_template_path.read_text(encoding="utf-8")
    raw_template_str = remove_html_comments(raw_template_str)
    raw_template_str = replace_source_json(raw_template_str, raw_json_data)
    raw_template_str = replace_doi(raw_template_str, raw_json_data)
    raw_template_str = replace_first_author(raw_template_str, raw_json_data)
    raw_template_str = replace_affiliations(raw_template_str, raw_json_data)
    raw_template_str = replace_title(raw_template_str, raw_json_data)
    raw_template_str = replace_abstract(raw_template_str, raw_json_data)
    raw_md_path.write_text(raw_template_str, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python scripts/add_raw_md.py raw/paper.json")
    main(sys.argv[1])
