import json
from pathlib import Path
import re
import sys

ROOT_PATH = Path(__file__).resolve().parent.parent
RAW_PATH = ROOT_PATH.joinpath("raw")


def get_raw_json_paths(raw_dir: str | Path = RAW_PATH) -> list[Path]:
    paths = Path(raw_dir).glob("*.json")
    return sorted(paths)


def find_author(given_name: str, family_name: str, raw_json_data: dict) -> list[str]:
    target = re.sub(r"\W+", "", f"{given_name} {family_name}").lower()
    for author in raw_json_data.get("author") or []:
        name = " ".join(filter(None, [author.get("given"), author.get("family")]))
        if re.sub(r"\W+", "", name).lower() == target:
            print(f"author: {raw_json_data.get('filename')}: {name}")



def find_corresponding(given_name: str, family_name: str, raw_json_data: dict) -> list[str]:
    target = re.sub(r"\W+", "", f"{given_name} {family_name}").lower()
    for author in raw_json_data.get("author-corresponding-openalex") or []:
        name = " ".join(filter(None, [author.get("given"), author.get("family")]))
        if re.sub(r"\W+", "", name).lower() == target:
            print(f"corresponding: {raw_json_data.get('filename')}: {name}")


def main(given_name: str, family_name: str) -> None:
    for raw_json_path in get_raw_json_paths():
        raw_json_data = json.load(raw_json_path.read_text(encoding="utf-8"))
        find_author(given_name, family_name, raw_json_data)
        find_corresponding(given_name, family_name, raw_json_data)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/search_a_author.py GIVEN FAMILY")
    main(sys.argv[1], sys.argv[2])
