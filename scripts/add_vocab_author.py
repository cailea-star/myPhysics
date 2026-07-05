import json
import re
import sys
from pathlib import Path


def affiliation_names(affiliations):
    return [item.get("name") if isinstance(item, dict) else str(item) for item in affiliations]


def get_corresponding_authors(raw_data):
    return [
        {
            "given": corresponding.get("given") or "",
            "family": corresponding.get("family") or "",
            "affiliation": affiliation_names(corresponding.get("raw_affiliation_strings") or []),
            "ORCID": corresponding.get("ORCID") or corresponding.get("raw_orcid"),
        }
        for corresponding in raw_data.get("author-corresponding-openalex") or []
    ]


def has_author(authors, raw_author):
    name_key = lambda author: re.sub(r"\W+", "", " ".join(filter(None, [author["given"], author["family"]]))).lower()
    for author in authors:
        if (
            raw_author["ORCID"] == author["ORCID"]
            and name_key(raw_author) == name_key(author)
        ):
            return True
    return False


def main(json_path):
    raw_data = json.loads(Path(json_path).read_text(encoding="utf-8"))
    raw_corresponding_authors = get_corresponding_authors(raw_data)

    authors_path = Path(__file__).resolve().parents[1] / "vocab" / "authors.json"
    authors = json.loads(authors_path.read_text(encoding="utf-8"))

    changed = False
    for raw_author in raw_corresponding_authors:
        if not raw_author["given"] and not raw_author["family"]: 
            print(f"skipping author with no name: {raw_author['raw_author_name']}")
            continue
        if has_author(authors, raw_author):
            print(f"author unchanged: {raw_author['given']} {raw_author['family']}")
            continue
        authors.append(raw_author)
        changed = True
        print(f"author added: {raw_author['given']} {raw_author['family']}")

    if changed:
        authors_path.write_text(json.dumps(authors, ensure_ascii=False, indent=4) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
