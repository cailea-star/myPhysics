from pathlib import Path
import sys


def record_quotation_linenumber(fullmd: str) -> list[int]:
    lines = fullmd.splitlines()
    linenumberlist = []
    for linenumber, line in enumerate(lines):
        if line.startswith("##### quotation-"):
            linenumberlist.append(linenumber)
    return linenumberlist


def rewrite_quotation_number(fullmd: str, linenumberlist: list[int]) -> str:
    lines = fullmd.splitlines()
    for i, linenumber in enumerate(linenumberlist):
        lines[linenumber] = f"##### quotation-{i+1:02d}"
    output = "\n".join(lines)
    if fullmd.endswith("\n"):
        output += "\n"
    return output


def main(mdpath: str | Path) -> Path:
    mdpath = Path(mdpath)
    fullmd = mdpath.read_text(encoding="utf-8")
    linenumberlist = record_quotation_linenumber(fullmd)
    output = rewrite_quotation_number(fullmd, linenumberlist)
    mdpath.write_text(output, encoding="utf-8")
    return mdpath


if __name__ == "__main__":
    sortmd_path = main(sys.argv[1])
    print(sortmd_path)
