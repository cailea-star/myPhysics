from pathlib import Path
import re
import sys



def find_quotation_lines(raw_md_lines: list[str]) -> list[int]:
    linenumberlist = []
    for linenumber, line in enumerate(raw_md_lines):
        if re.match(r"^##### quotation-\d+\s*$", line):
            linenumberlist.append(linenumber)
    return linenumberlist



def rewrite_quotation_number(raw_md_lines: list[str], linenumber_list: list[int]) -> str:
    for i, linenumber in enumerate(linenumber_list):
        raw_md_lines[linenumber] = f"##### quotation-{i+1:02d}\n"
    raw_md_str_new = "".join(raw_md_lines)
    return raw_md_str_new


def main(raw_md_path: str | Path) -> Path:
    raw_md_path = Path(raw_md_path)
    raw_md_lines_str_list = raw_md_path.read_text(encoding="utf-8").splitlines(keepends=True)
    linenumber_list = find_quotation_lines(raw_md_lines_str_list)
    raw_md_str_new = rewrite_quotation_number(raw_md_lines_str_list, linenumber_list)
    raw_md_path.write_text(raw_md_str_new, encoding="utf-8")
    return raw_md_path


if __name__ == "__main__":
    raw_md_path = main(sys.argv[1])
    print(raw_md_path)
