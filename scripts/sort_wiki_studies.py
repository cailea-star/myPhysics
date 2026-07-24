from pathlib import Path
import re
import sys


SECTION_HEADER = re.compile(r"(?m)^### Previous Studies\s*$")
NEXT_SECTION = re.compile(r"(?m)^###\s+")
STUDY_HEADER = re.compile(r"(?m)^#####\s+.*$")

def match_previous_studies(wiki_md_str: str) -> str:
    section = SECTION_HEADER.search(wiki_md_str)
    if not section:
        raise ValueError("Previous Studies section not found")
    next_section = NEXT_SECTION.search(wiki_md_str, section.end())
    section_end = next_section.start() if next_section else len(wiki_md_str)
    return wiki_md_str[section.start():section_end]


def split_studies(wiki_studies_str: str) -> list[str]:
    headers = list(STUDY_HEADER.finditer(wiki_studies_str))
    studies = []
    for index, header in enumerate(headers):
        study_end = headers[index + 1].start() if index + 1 < len(headers) else len(wiki_studies_str)
        studies.append(wiki_studies_str[header.start():study_end].rstrip())
    return studies


def resolve_stem(stem: str) -> dict[str, str]:
    fields = ("author", "year", "journal", "volume", "number", "page")
    match = re.fullmatch(
        r"(?P<author>.+?)_Y\.(?P<year>[^_]*)_(?P<journal>.*?)_Vol\."
        r"(?P<volume>.*?)Nol\.(?P<number>.*?)P\.(?P<page>.*)",
        stem,
    )
    return match.groupdict() if match else dict.fromkeys(fields, "")


def resolve_study(study_str: str) -> dict[str, str]:
    references = re.search(r"```references\s*\n(.*?)\n```", study_str, re.S)
    if not references: raise ValueError("Study has no references block")
    stems = re.findall(r"(?m)^\s*-\s+([^:\r\n]+):", references.group(1))
    if len(stems) != 1: raise ValueError("Study must reference exactly one Raw paper")
    stem = resolve_stem(stems[0].strip())
    return {"author": stem["author"], "year": stem["year"]}


def replace_study_title(study_str: str) -> str:
    study_data = resolve_study(study_str)
    author = study_data["author"].replace("_", " ")
    title = f"##### {author} ({study_data['year']})"
    return STUDY_HEADER.sub(title, study_str, count=1)


def sort_studies(studies_list: list[str]) -> list[str]:
    def sort_key(study: str) -> tuple[str, str]:
        study_data = resolve_study(study)
        return study_data["author"], study_data["year"]

    sorted_studies = sorted(studies_list, key=sort_key)
    for index, study in enumerate(sorted_studies):
        sorted_studies[index] = replace_study_title(study)
    return sorted_studies


def main(wiki_md_path: str | Path) -> Path:
    wiki_md_path = Path(wiki_md_path)
    wiki_md_str = wiki_md_path.read_text(encoding="utf-8")
    wiki_studies_str = match_previous_studies(wiki_md_str)
    studies_list = split_studies(wiki_studies_str)
    if not studies_list: return wiki_md_path

    prefix = wiki_studies_str[:STUDY_HEADER.search(wiki_studies_str).start()]
    suffix = wiki_studies_str[len(wiki_studies_str.rstrip()):]
    sortedwiki = prefix + "\n\n".join(sort_studies(studies_list)) + suffix
    wiki_md_path.write_text(wiki_md_str.replace(wiki_studies_str, sortedwiki, 1), encoding="utf-8")
    return wiki_md_path


if __name__ == "__main__":
    sorted_path = main(sys.argv[1])
    print(sorted_path)
