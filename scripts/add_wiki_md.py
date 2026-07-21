import re
import sys
import os
from pathlib import Path


def remove_html_comments(text):
    return re.sub(r"<!--.*?-->", "", text, flags=re.DOTALL)


def replace_tags_set(text, tag):
    return text.replace("[this_tags_set]", tag)


def replace_tag_types(text, tag_type):
    return text.replace("[this_tag_types]", tag_type)


def replace_tag_name(text, tag):
    return text.replace("[tag_name]", tag)


def main(tag, tag_type):
    if not re.fullmatch(r"[a-z0-9_]+", tag): raise SystemExit("tag must contain only lowercase letters, digits, and underscores")
    md_path = Path(os.path.join(os.path.dirname(__file__), "..", "wiki", f"{tag}.md"))
    if md_path.exists(): raise SystemExit(f"{md_path} already exists")

    template = "add_wiki_method.md" if tag_type == "method" else "add_wiki_topic.md"
    text = Path(os.path.join(os.path.dirname(__file__), template)).read_text(encoding="utf-8")
    text = remove_html_comments(text)
    text = replace_tags_set(text, tag)
    text = replace_tag_types(text, tag_type)
    text = replace_tag_name(text, tag)

    md_path.parent.mkdir(exist_ok=True)
    md_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/add_wiki_md.py TAG TAG_TYPE")
    main(sys.argv[1], sys.argv[2])
