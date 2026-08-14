import json
from pathlib import Path
import re
import sys

PAPERS_PATH = Path(__file__).resolve().parent.parent
RAW_PATH = PAPERS_PATH.joinpath("raw")
TMP_PATH = PAPERS_PATH.joinpath("tmp")
TAGS_PATH = PAPERS_PATH.joinpath("vocab", "tags.json")


def check_tag(tag: str) -> dict:
    tag_data = {item["tag"]: item for item in json.loads(TAGS_PATH.read_text(encoding="utf-8"))}
    if tag not in tag_data: raise SystemExit(f"invalid tag: {tag}")
    return tag_data[tag]


# Resolve quotations from raw markdown files;

def find_section(raw_md_str: str, sectionname: str) -> str:
    sectionname = sectionname.strip()
    escaped_sectionname = re.escape(sectionname)
    header = re.search(rf"(?m)^(#+)\s+{escaped_sectionname}\s*$", raw_md_str)
    if not header: return ""
    body_start = header.end()
    header_marks = header.group(1)
    header_level = len(header_marks)
    section_tail = raw_md_str[body_start:]
    next_header = re.search(rf"(?m)^#{{1,{header_level}}}\s+", section_tail)
    if next_header: body_end = body_start + next_header.start()
    else: body_end = len(raw_md_str)
    section_body = raw_md_str[body_start:body_end]
    return section_body.strip()


def split_section_quotations(raw_md_section_str: str) -> list[dict]:
    quotations = []

    # split the section into blocks by "##### " header, don't include the first block which is the section header
    blocks = re.split(r"(?m)^#####\s+", raw_md_section_str)
    for block in blocks[1:]:
        _, _, body = block.partition("\n")
        # search for the tags block
        tag_block = re.search(r"```tags\s*\n(.*?)\n```", body, re.S)
        if not tag_block: continue
        # search for claim-type and tags in the tags block
        tag_text = tag_block.group(1)
        claim_match = re.search(r"(?m)^\[claim-type\]:\s*(.+)$", tag_text)
        tags_match = re.search(r"(?m)^\[tags\]:\s*(.+)$", tag_text)
        source_match = re.search(r"(?m)^\[source\]:\s*(.+)$", tag_text)
        ref_match = re.search(r"(?m)^\[ref\]:\s*(.+)$", tag_text)
        doi_match = re.search(r"(?m)^\[doi\]:\s*(.+)$", tag_text)

        taglist = []
        tags_text = tags_match.group(1)
        for tag in tags_text.split(","):
            tag = tag.strip()
            if tag: taglist.append(tag)
        # search for source in the tags block
        claim_kind = claim_match.group(1).strip()
        quote = (body[:tag_block.start()] + body[tag_block.end():]).strip()
        source = source_match.group(1).strip()
        quotations.append(
            {
                "claim-type": claim_kind,
                "tags": taglist,
                "quote": quote,
                "source": source,
                "ref": ref_match.group(1).strip() if ref_match else "",
                "doi": doi_match.group(1).strip() if doi_match else "",
            }
        )
    return quotations


# resolve quotations from certain section, and match them with claim-type and tags;

def resolve_section_quotations(raw_md_path: Path, sectionname: str) -> list[dict]:
    raw_md_str = raw_md_path.read_text(encoding="utf-8")
    raw_md_section_str = find_section(raw_md_str, sectionname)
    raw_md_section_quotations_list = split_section_quotations(raw_md_section_str)
    for raw_md_quotation in raw_md_section_quotations_list:
        raw_md_quotation["stem"] = raw_md_path.stem
        raw_md_quotation["source"] = raw_md_path.stem + ": " + raw_md_quotation["source"]
    return raw_md_section_quotations_list


def match_quotations(quotations_list: list[dict], claim_kind: str | None, tags_list: list[str]) -> list[dict]:
    raw_md_quotations_list_filtered = []
    for raw_md_quotation in quotations_list:
        if claim_kind is not None and raw_md_quotation["claim-type"] != claim_kind: continue
        matched = True
        for tag in tags_list:
            if tag not in raw_md_quotation["tags"]:
                matched = False
                break
        if matched:
            raw_md_quotations_list_filtered.append(raw_md_quotation)
    return raw_md_quotations_list_filtered


# Sort and format quotations for output


def sort_quotations(quotations_list: list[dict]) -> list[dict]:
    raw_md_quotations_list_sorted = sorted(quotations_list, key=lambda quote: quote["claim-type"])
    return raw_md_quotations_list_sorted


def format_quotations(quotations_list: list[dict]) -> list[str]:
    lines = []
    last_claim_kind = ""
    for quote in sort_quotations(quotations_list):
        if quote["claim-type"] != last_claim_kind:
            lines += [f"## {quote['claim-type']}", ""]
            last_claim_kind = quote["claim-type"]
        tags = ", ".join(quote["tags"])
        lines += ["```tags"]
        lines += [f"[tags]: {tags}", f"[source]: {quote['source']}"]
        if quote["ref"]: lines.append(f"[ref]: {quote['ref']}")
        if quote["doi"]: lines.append(f"[doi]: {quote['doi']}")
        lines += ["```", ""]
        lines += [quote["quote"], "---", ""]
    return lines


# Collect papers appearing in quotations, and sort them by year and author

def resolve_stem(stem: str) -> dict[str, str]:
    fields = ("author", "year", "journal", "volume", "number", "page")
    match = re.fullmatch(
        r"(?P<author>.+?)_Y\.(?P<year>[^_]*)_(?P<journal>.*?)_Vol\."
        r"(?P<volume>.*?)Nol\.(?P<number>.*?)P\.(?P<page>.*)",
        stem,
    )
    return match.groupdict() if match else dict.fromkeys(fields, "")


def get_paper_sort_key(stem: str) -> tuple[str, bool, str, str]:
    paper = resolve_stem(stem)
    return paper["author"], paper["year"] == "", paper["year"], stem


def format_stems(quotations_list: list[dict]) -> list[str]:
    papers = {}
    for quote in quotations_list:
        papers.setdefault(quote["stem"], []).append(quote)
    lines = ["## Papers", "", "| stem | year | quotes | claim-type | tags |", "|---|---:|---:|---|---|"]
    for stem in sorted(papers, key=get_paper_sort_key):
        year = resolve_stem(stem)["year"]
        quotes = papers[stem]
        claim_kinds = ", ".join(dict.fromkeys(quote["claim-type"] for quote in quotes))
        tags = ", ".join(dict.fromkeys(tag for quote in quotes for tag in quote["tags"]))
        lines.append(f"| {stem} | {year} | {len(quotes)} | {claim_kinds} | {tags} |")
    return lines + [""]


# Collect raw/*.md paths, and main function;

def get_raw_md_paths(raw_dir: str | Path = RAW_PATH) -> list[Path]:
    paths = Path(raw_dir).glob("*.md")
    return sorted(paths)


def main(tagnames: list[str]) -> list[Path]:
    sections = (
        ("Motivation", "Motivation"),
        ("Methods", "Methods"),
        ("Results", "Results"),
        ("Meanings", "Meanings"),
        ("Secondary Citations", "Secondary"),
    )
    for tagname in tagnames: check_tag(tagname)

    group_name = tagnames[0]
    output_paths = []
    for sectionname, suffix in sections:
        quotations = []
        for raw_md_path in get_raw_md_paths():
            quotations.extend(resolve_section_quotations(raw_md_path, sectionname))
        matched = [quote for quote in quotations if any(tag in quote["tags"] for tag in tagnames)]
        lines = [f"# {group_name}_{suffix}", ""] + format_stems(matched) + format_quotations(matched)
        output_path = TMP_PATH.joinpath(f"{group_name}_{suffix}.md")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(lines), encoding="utf-8")
        output_paths.append(output_path)
    return output_paths


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit("usage: python papers/scripts/search_a_tag.py TAG [TAG ...]")
    for output_path in main(sys.argv[1:]):
        print(output_path)
