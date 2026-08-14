import json
from pathlib import Path
import re
import sys


ROOT_PATH = Path(__file__).resolve().parent.parent
RAW_PATH = ROOT_PATH.joinpath("raw")
TAGS_PATH = ROOT_PATH.joinpath("vocab", "tags.json")

TAGS_LINE = re.compile(r"^([ \t]*\[tags\]:[ \t]*)([^\r\n]*)(\r?\n?)$")


def check_tag(tag: str) -> dict:
    tag_data = {item["tag"]: item for item in json.loads(TAGS_PATH.read_text(encoding="utf-8"))}
    if tag not in tag_data: raise SystemExit(f"invalid tag: {tag}")
    return tag_data[tag]


def get_raw_md_paths(raw_dir: str | Path = RAW_PATH) -> list[Path]:
    """Return sorted raw/*.md paths."""
    return sorted(Path(raw_dir).glob("*.md"))


def find_tags_lines(raw_md_lines_str_list: list[str]) -> list[int]:
    """Return zero-based indices of exact `[tags]:` lines."""
    return [i for i, line in enumerate(raw_md_lines_str_list) if TAGS_LINE.match(line)]


def replace_tag(raw_md_line_str: str, old_tag: str, new_tag: str) -> str:
    """Replace an exact tag, preserving layout and deduplicating merges.

    >>> replace_tag("[tags]: old, keep, new\\n", "old", "new")
    '[tags]: new, keep\\n'
    """
    match = TAGS_LINE.match(raw_md_line_str)
    if not match: return raw_md_line_str
    prefix, value, ending = match.groups()
    tags = [tag.strip() for tag in value.split(",")]
    if old_tag not in tags: return raw_md_line_str
    tags = dict.fromkeys(new_tag if tag == old_tag else tag for tag in tags)
    return prefix + ", ".join(tags) + ending


def main(old_tag: str, new_tag: str) -> list[Path]:
    """Update matched files and return changed paths."""
    old_tag, new_tag = old_tag.strip(), new_tag.strip()
    check_tag(old_tag)
    check_tag(new_tag)

    changed = []
    for raw_md_path in get_raw_md_paths():
        with raw_md_path.open(encoding="utf-8", newline="") as file:
            raw_md_str_old = file.read()
        raw_md_lines = raw_md_str_old.splitlines(keepends=True)
        for tagline_number in find_tags_lines(raw_md_lines):
            raw_md_lines[tagline_number] = replace_tag(raw_md_lines[tagline_number], old_tag, new_tag)
        raw_md_str_new = "".join(raw_md_lines)
        if raw_md_str_new != raw_md_str_old:
            with raw_md_path.open("w", encoding="utf-8", newline="") as file:
                file.write(raw_md_str_new)
            changed.append(raw_md_path)
    return changed


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/rename_raw_tag.py OLD NEW")
    paths = main(sys.argv[1], sys.argv[2])
    for path in paths:
        print(path)
    print(f"renamed {sys.argv[1]} -> {sys.argv[2]} in {len(paths)} file(s).")
