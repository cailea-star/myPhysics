from pathlib import Path
import re
import sys


ROOT_PATH = Path(__file__).resolve().parent.parent
RAW_PATH = ROOT_PATH / "raw"
TAGS_LINE = re.compile(r"^([ \t]*\[tags\]:[ \t]*)([^\r\n]*)(\r?\n?)$")


def get_md_paths(raw_dir: Path = RAW_PATH) -> list[Path]:
    """Return sorted raw/*.md paths."""
    return sorted(raw_dir.glob("*.md"))


def match_tags_lines(lines: list[str]) -> list[int]:
    """Return zero-based indices of exact `[tags]:` lines."""
    return [i for i, line in enumerate(lines) if TAGS_LINE.match(line)]


def replace_tag(line: str, old: str, new: str) -> str:
    """Replace an exact tag, preserving layout and deduplicating merges.

    >>> replace_tag("[tags]: old, keep, new\\n", "old", "new")
    '[tags]: new, keep\\n'
    """
    match = TAGS_LINE.match(line)
    if not match: return line
    prefix, value, ending = match.groups()
    tags = [tag.strip() for tag in value.split(",")]
    if old not in tags: return line
    tags = dict.fromkeys(new if tag == old else tag for tag in tags)
    return prefix + ", ".join(tags) + ending


def main(old: str, new: str) -> list[Path]:
    """Update matched files and return changed paths."""
    old, new = old.strip(), new.strip()
    if not old or not new or old == new or "," in old + new:
        raise SystemExit("tags must be distinct, non-empty, and contain no comma")

    changed = []
    for path in get_md_paths():
        with path.open(encoding="utf-8", newline="") as file:
            original = file.read()
        lines = original.splitlines(keepends=True)
        for i in match_tags_lines(lines):
            lines[i] = replace_tag(lines[i], old, new)
        text = "".join(lines)
        if text != original:
            with path.open("w", encoding="utf-8", newline="") as file:
                file.write(text)
            changed.append(path)
    return changed


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/rename_raw_tag.py OLD NEW")
    paths = main(sys.argv[1], sys.argv[2])
    for path in paths:
        print(path)
    print(f"renamed {sys.argv[1]} -> {sys.argv[2]} in {len(paths)} file(s).")
