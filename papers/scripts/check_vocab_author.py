import json
import re
import sys
from pathlib import Path

PAPERS_PATH = Path(__file__).resolve().parent.parent


def resolve_affiliation_names(raw_json_affiliations):
    return [item.get("name") if isinstance(item, dict) else str(item) for item in raw_json_affiliations]


def resolve_corresponding_authors(raw_json_data):
    return [
        {
            "given": corresponding.get("given") or "",
            "family": corresponding.get("family") or "",
            "affiliation": resolve_affiliation_names(corresponding.get("raw_affiliation_strings") or []),
            "ORCID": corresponding.get("ORCID") or corresponding.get("raw_orcid"),
        }
        for corresponding in raw_json_data.get("author-corresponding-openalex") or []
    ]


def is_recorded_author(vocab_json_authors, raw_json_author):
    name_key = lambda author: re.sub(r"\W+", "", " ".join(filter(None, [author["given"], author["family"]]))).lower()
    for author in vocab_json_authors:
        if (
            raw_json_author["ORCID"] == author["ORCID"]
            and name_key(raw_json_author) == name_key(author)
        ):
            return True
    return False


def main(raw_json_path):
    raw_json_data = json.loads(Path(raw_json_path).read_text(encoding="utf-8"))
    raw_json_corresponding_authors = resolve_corresponding_authors(raw_json_data)

    vocab_json_authors_path = PAPERS_PATH.joinpath("vocab", "authors.json")
    vocab_json_authors = json.loads(vocab_json_authors_path.read_text(encoding="utf-8"))

    changed = False
    for raw_author in raw_json_corresponding_authors:
        if not raw_author["given"] and not raw_author["family"]: 
            print("skipping author with no name")
            continue
        if is_recorded_author(vocab_json_authors, raw_author):
            print(f"author unchanged: {raw_author['given']} {raw_author['family']}")
            continue
        vocab_json_authors.append(raw_author)
        changed = True
        print(f"author added: {raw_author['given']} {raw_author['family']}")

    if changed:
        vocab_json_authors_path.write_text(json.dumps(vocab_json_authors, ensure_ascii=False, indent=4) + "\n", encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python papers/scripts/check_vocab_author.py RAW_JSON_PATH")
    main(sys.argv[1])
