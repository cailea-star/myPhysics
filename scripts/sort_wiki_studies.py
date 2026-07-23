from pathlib import Path
import re
import sys

from search_a_tag import summary_split_filename


SECTION_HEADER = re.compile(r"(?m)^### Previous Studies\s*$")
NEXT_SECTION = re.compile(r"(?m)^###\s+")
STUDY_HEADER = re.compile(r"(?m)^#####\s+")


def match_previous_studies(fullwiki: str) -> str:
    section = SECTION_HEADER.search(fullwiki)
    if not section:
        raise ValueError("Previous Studies section not found")
    next_section = NEXT_SECTION.search(fullwiki, section.end())
    section_end = next_section.start() if next_section else len(fullwiki)
    return fullwiki[section.start():section_end]


def split_studies(studieswiki_str: str) -> list[str]:
    headers = list(STUDY_HEADER.finditer(studieswiki_str))
    studies = []
    for index, header in enumerate(headers):
        study_end = headers[index + 1].start() if index + 1 < len(headers) else len(studieswiki_str)
        studies.append(studieswiki_str[header.start():study_end].rstrip())
    return studies


def resolve_filename(study_str: str) -> dict[str, str]:
    references = re.search(r"```references\s*\n(.*?)\n```", study_str, re.S)
    if not references:
        raise ValueError("Study has no references block")
    stems = re.findall(r"(?m)^\s*-\s+([^:\r\n]+):", references.group(1))
    if len(stems) != 1:
        raise ValueError("Study must reference exactly one Raw paper")
    filename = summary_split_filename(stems[0].strip())
    return {"author": filename["author"], "year": filename["year"]}


def sort_studies(studies: list[str]) -> list[str]:
    def sort_key(study: str) -> tuple[bool, str, str]:
        filename = resolve_filename(study)
        return filename["year"] == "", filename["year"], filename["author"]

    return sorted(studies, key=sort_key)


def main(wikipath: str | Path) -> Path:
    wikipath = Path(wikipath)
    fullwiki = wikipath.read_text(encoding="utf-8")
    studieswiki = match_previous_studies(fullwiki)
    studies = split_studies(studieswiki)
    if not studies:
        return wikipath

    prefix = studieswiki[:STUDY_HEADER.search(studieswiki).start()]
    suffix = studieswiki[len(studieswiki.rstrip()):]
    sortedwiki = prefix + "\n\n".join(sort_studies(studies)) + suffix
    wikipath.write_text(fullwiki.replace(studieswiki, sortedwiki, 1), encoding="utf-8")
    return wikipath


if __name__ == "__main__":
    sorted_path = main(sys.argv[1])
    print(sorted_path)
