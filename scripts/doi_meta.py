import html
import json
import re
import sys
import urllib.parse
import urllib.request

HEADERS = {"User-Agent": "doi-metadata-script/0.1 (mailto:your@email.com)"}


def get_json(url):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)


def get_doi_metadata(doi):
    """ get metadata for a DOI from the CrossRef API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    return get_json(url)["message"]


def clean_abstract(text):
    if not text:
        return None
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text)
    text = re.sub(r"^\s*Abstract\s+", "", " ".join(text.split()), flags=re.I)
    return text or None


def get_abstract_openAlex(doi):
    """ get abstract for a DOI from the OpenAlex API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.openalex.org/works/https://doi.org/" + urllib.parse.quote(doi, safe="/")
    try:
        index = get_json(url).get("abstract_inverted_index") or {}
    except Exception:
        return None
    return clean_abstract(" ".join(word for _, word in sorted(
        (i, word) for word, positions in index.items() for i in positions
    )))


def get_abstract_inspire(doi):
    """ get abstract for a DOI from the INSPIRE API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    params = urllib.parse.urlencode({"q": f"doi:{doi}", "size": 1})
    try:
        abstracts = get_json("https://inspirehep.net/api/literature?" + params)["hits"]["hits"][0]["metadata"].get("abstracts") or []
    except Exception:
        return None
    return clean_abstract(abstracts[0].get("value") if abstracts else None)


def get_doi_abstractdata(doi):
    """ get metadata and abstract for a DOI.
    """
    msg = get_doi_metadata(doi)
    return {
        "source": msg["source"],
        "language": msg.get("language"),
        "title": msg.get("title"),
        "original-title": msg.get("original-title"),
        "short-title": msg.get("short-title"),
        "subtitle": msg.get("subtitle"),
        "author": msg.get("author"),
        "DOI": msg.get("DOI"),
        "URL": msg.get("URL"),
        "is-referenced-by-count": msg.get("is-referenced-by-count"),
        "reference-count": msg.get("reference-count"),
        "abstract": clean_abstract(msg.get("abstract")) or get_abstract_openAlex(doi) or get_abstract_inspire(doi),
        "publisher": msg.get("publisher"),
        "container-title": msg.get("container-title"),
        "short-container-title": msg.get("short-container-title"),
        "type": msg.get("type"),
        "volume": msg.get("volume"),
        "page": msg.get("page"),
        "issued": msg.get("issued"),
        "article-number": msg.get("article-number"),
        "ISSN": msg.get("ISSN"),
        "funder": msg.get("funder"),
        "reference": msg.get("reference"),   
    }


if __name__ == "__main__":
    doi = sys.argv[1]
    print(json.dumps(get_doi_abstractdata(doi), ensure_ascii=False, indent=2))
