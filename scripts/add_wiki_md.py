import json
import re
import sys
from pathlib import Path

ROOT_PATH = Path(__file__).resolve().parent.parent
WIKI_PATH = ROOT_PATH.joinpath("wiki")
TAGS_PATH = ROOT_PATH.joinpath("vocab", "tags.json")


def check_tag(tag: str) -> dict:
    tag_data = {item["tag"]: item for item in json.loads(TAGS_PATH.read_text(encoding="utf-8"))}
    if tag not in tag_data: raise SystemExit(f"invalid tag: {tag}")
    return tag_data[tag]


def remove_html_comments(wiki_template_str):
    return re.sub(r"<!--.*?-->", "", wiki_template_str, flags=re.DOTALL)


def replace_tag_name(wiki_template_str, tag):
    return wiki_template_str.replace("[this_tag_name]", tag)


def replace_tag_types(wiki_template_str, tag_types):
    return wiki_template_str.replace("[this_tag_types]", ", ".join(tag_types))


def main(tag):
    matched_tag = check_tag(tag)

    wiki_md_path = WIKI_PATH.joinpath(f"{tag}.md")
    if wiki_md_path.exists(): raise SystemExit(f"{wiki_md_path} already exists")

    tag_types = matched_tag["types"]
    wiki_template_filename = "add_wiki_method.md" if "method" in tag_types else "add_wiki_topic.md"
    wiki_template_path = Path(__file__).resolve().with_name(wiki_template_filename)
    wiki_template_str = wiki_template_path.read_text(encoding="utf-8")
    wiki_template_str = remove_html_comments(wiki_template_str)
    wiki_template_str = replace_tag_name(wiki_template_str, tag)
    wiki_template_str = replace_tag_types(wiki_template_str, tag_types)

    wiki_md_path.parent.mkdir(exist_ok=True)
    wiki_md_path.write_text(wiki_template_str, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python scripts/add_wiki_md.py TAG")
    main(sys.argv[1])
