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
        msg = json.load(r)
    
    return msg["message"]

def get_doi_abstractdata(doi):
    msg = get_doi_metadata(doi)
    return {
        "source": msg["source"],
        "language": msg["language"],
        "title": msg["title"],
        "original-title": msg["original-title"],
        "short-title": msg["short-title"],
        "subtitle": msg["subtitle"],
        "author": msg["author"],
        "DOI": msg["DOI"],
        "URL": msg["URL"],
        "is-referenced-by-count": msg["is-referenced-by-count"],
        "reference-count": msg["reference-count"],
        "publisher": msg["publisher"],
        "container-title": msg["container-title"],
        "short-container-title": msg["short-container-title"],
        "type": msg["type"],
        "volume": msg["volume"],
        "page": msg["page"],
        "issued": msg["issued"],
        "article-number": msg["article-number"],
        "ISSN": msg["ISSN"],
        "funder": msg["funder"],
        "reference": msg["reference"],   
    }

if __name__ == "__main__":
    doi = sys.argv[1]
    
    print(json.dumps(get_doi_abstractdata(doi), ensure_ascii=False, indent=2))
