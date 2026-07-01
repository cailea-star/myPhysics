import json
import sys
import urllib.parse
import urllib.request

def get_doi_metadata(doi):
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "doi-metadata-script/0.1 (mailto:your@email.com)"}
    )

    with urllib.request.urlopen(req, timeout=20) as r:
        msg = json.load(r)["message"]

    return msg

if __name__ == "__main__":
    doi = sys.argv[1]
    print(json.dumps(get_doi_metadata(doi), ensure_ascii=False, indent=2))