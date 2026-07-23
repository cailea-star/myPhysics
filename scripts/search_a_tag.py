import os
from pathlib import Path
import re
import sys

ROOT_PATH = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def get_md_paths(raw_dir: str | Path = os.path.join(ROOT_PATH, "raw")) -> list[Path]:
    paths = Path(raw_dir).glob("*.md")
    return sorted(paths)


def search_section(fullmd: str, sectionname: str) -> str:
    sectionname = sectionname.strip()
    escaped_sectionname = re.escape(sectionname)
    header = re.search(rf"(?m)^(#+)\s+{escaped_sectionname}\s*$", fullmd)
    if not header: return ""
    body_start = header.end()
    header_marks = header.group(1)
    header_level = len(header_marks)
    section_tail = fullmd[body_start:]
    next_header = re.search(rf"(?m)^#{{1,{header_level}}}\s+", section_tail)
    if next_header: body_end = body_start + next_header.start()
    else: body_end = len(fullmd)
    section_body = fullmd[body_start:body_end]
    return section_body.strip()


def search_quotations(sectionmd: str) -> list[dict]:
    quotations = []

    # split the section into blocks by "##### " header, don't include the first block which is the section header
    blocks = re.split(r"(?m)^#####\s+", sectionmd)
    for block in blocks[1:]:
        _, _, body = block.partition("\n")
        # search for the tags block
        tag_block = re.search(r"```tags\s*\n(.*?)\n```", body, re.S)
        # search for claim_type and tags in the tags block
        tag_text = tag_block.group(1)
        claim_type_match = re.search(r"(?m)^\[claim_type\]:\s*(.+)$", tag_text)
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
        claim_type = claim_type_match.group(1).strip()
        quote = (body[:tag_block.start()] + body[tag_block.end():]).strip()
        source = source_match.group(1).strip()
        quotations.append(
            {
                "claim_type": claim_type,
                "tags": taglist,
                "quote": quote,
                "source": source,
                "ref": ref_match.group(1).strip() if ref_match else "",
                "doi": doi_match.group(1).strip() if doi_match else "",
            }
        )
    return quotations


def parse_quotations(md_path: Path, sectionlist: list[str]) -> list[dict]:
    fullmd = md_path.read_text(encoding="utf-8")
    quotations = []
    for sectionname in sectionlist:
        sectionmd = search_section(fullmd, sectionname)
        section_quotes = search_quotations(sectionmd)
        quotations.extend(section_quotes)
    for quote in quotations:
        quote["filename"] = md_path.stem
        quote["source"] = md_path.stem + ": " + quote["source"]
    return quotations


def filter_quotations(quotations: list[dict], claim_type: str | None, tags_list: list[str]) -> list[dict]:
    filtered = []
    for quote in quotations:
        if claim_type is not None and quote["claim_type"] != claim_type: continue
        matched = True
        for tag in tags_list:
            if tag not in quote["tags"]:
                matched = False
                break
        if matched:
            filtered.append(quote)
    return filtered


def sort_quotations(quotations: list[dict]) -> list[dict]:
    sorted_quotations = sorted(quotations, key=lambda quote: quote["claim_type"])
    return sorted_quotations


def summary_split_filename(filename: str) -> dict[str, str]:
    fields = ("author", "year", "journal", "volume", "number", "page")
    match = re.fullmatch(
        r"(?P<author>.+?)_Y\.(?P<year>[^_]*)_(?P<journal>.*?)_Vol\."
        r"(?P<volume>.*?)Nol\.(?P<number>.*?)P\.(?P<page>.*)",
        filename,
    )
    return match.groupdict() if match else dict.fromkeys(fields, "")


def summary_paper_sort_key(filename: str) -> tuple[bool, str, str]:
    year = summary_split_filename(filename)["year"]
    return year == "", year, filename


def summary_format(quotations: list[dict]) -> list[str]:
    papers = {}
    for quote in quotations:
        papers.setdefault(quote["filename"], []).append(quote)
    lines = ["## Papers", "", "| filename | year | quotes | claim-types | tags |", "|---|---:|---:|---|---|"]
    for filename in sorted(papers, key=summary_paper_sort_key):
        year = summary_split_filename(filename)["year"]
        quotes = papers[filename]
        claim_types = ", ".join(dict.fromkeys(quote["claim_type"] for quote in quotes))
        tags = ", ".join(dict.fromkeys(tag for quote in quotes for tag in quote["tags"]))
        lines.append(f"| {filename} | {year} | {len(quotes)} | {claim_types} | {tags} |")
    return lines + [""]


def format_quotations(quotations: list[dict]) -> list[str]:
    lines = []
    last_claim_type = ""
    for quote in sort_quotations(quotations):
        if quote["claim_type"] != last_claim_type:
            lines += [f"## {quote['claim_type']}", ""]
            last_claim_type = quote["claim_type"]
        tags = ", ".join(quote["tags"])
        lines += ["```tags"]
        lines += [f"[tags]: {tags}", f"[source]: {quote['source']}"]
        if quote["ref"]: lines.append(f"[ref]: {quote['ref']}")
        if quote["doi"]: lines.append(f"[doi]: {quote['doi']}")
        lines += ["```", ""]
        lines += [quote["quote"], "---", ""]
    return lines


def main(tagname: str, sectionname: str, outfile_path: str | Path) -> Path:
    quotations = []
    query_name = f"{tagname}_{sectionname}"
    for path in get_md_paths():
        path_quotes = parse_quotations(path, [sectionname])
        quotations.extend(path_quotes)
    target_quotes = filter_quotations(quotations, None, [tagname])
    lines = [f"# {query_name}", ""] + summary_format(target_quotes) + format_quotations(target_quotes)
    outfile_path = Path(outfile_path)
    outfile_path.parent.mkdir(parents=True, exist_ok=True)
    outfile_path.write_text("\n".join(lines), encoding="utf-8")
    return outfile_path


if __name__ == "__main__":
    tagname = sys.argv[1]
    motivation_path = os.path.join(ROOT_PATH, "tmp", f"{tagname}_Motivation.md")
    main(tagname, "Motivation", motivation_path)
    print(motivation_path)
    methods_path = os.path.join(ROOT_PATH, "tmp", f"{tagname}_Methods.md")
    main(tagname, "Methods", methods_path)
    print(methods_path)
    results_path = os.path.join(ROOT_PATH, "tmp", f"{tagname}_Results.md")
    main(tagname, "Results", results_path)
    print(results_path)
    meanings_path = os.path.join(ROOT_PATH, "tmp", f"{tagname}_Meanings.md")
    main(tagname, "Meanings", meanings_path)
    print(meanings_path)
    secondary_cpath = os.path.join(ROOT_PATH, "tmp", f"{tagname}_Secondary.md")
    main(tagname, "Secondary Citations", secondary_cpath)
    print(secondary_cpath)
