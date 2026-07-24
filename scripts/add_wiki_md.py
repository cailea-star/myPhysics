import re
import sys
from pathlib import Path

ROOT_PATH = Path(__file__).resolve().parent.parent
WIKI_PATH = ROOT_PATH.joinpath("wiki")

def remove_html_comments(wiki_template_str):
    return re.sub(r"<!--.*?-->", "", wiki_template_str, flags=re.DOTALL)


def replace_tags_set(wiki_template_str, tag):
    return wiki_template_str.replace("[this_tags_set]", tag)


def replace_tag_types(wiki_template_str, tag_type):
    return wiki_template_str.replace("[this_tag_types]", tag_type)


def replace_tag_name(wiki_template_str, tag):
    return wiki_template_str.replace("[tag_name]", tag)


def main(tag, tag_type):
    wiki_md_path = WIKI_PATH.joinpath(f"{tag}.md")
    if wiki_md_path.exists(): raise SystemExit(f"{wiki_md_path} already exists")

    wiki_template_filename = "add_wiki_method.md" if tag_type == "method" else "add_wiki_topic.md"
    wiki_template_path = Path(__file__).resolve().with_name(wiki_template_filename)
    wiki_template_str = wiki_template_path.read_text(encoding="utf-8")
    wiki_template_str = remove_html_comments(wiki_template_str)
    wiki_template_str = replace_tags_set(wiki_template_str, tag)
    wiki_template_str = replace_tag_types(wiki_template_str, tag_type)
    wiki_template_str = replace_tag_name(wiki_template_str, tag)

    wiki_md_path.parent.mkdir(exist_ok=True)
    wiki_md_path.write_text(wiki_template_str, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/add_wiki_md.py TAG TAG_TYPE")
    main(sys.argv[1], sys.argv[2])
