import html
import json
import os
import re
import sys
import urllib.parse
import urllib.request

HEADERS = {"User-Agent": "doi-metadata-script/0.1 (mailto:your@email.com)"}


def get_json(url):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)


def metadata_from_CrossRef(doi):
    """ get metadata for a DOI from the CrossRef API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    return get_json(url)["message"]


def metadata_from_Inspire(doi):
    """ get metadata for a DOI from the INSPIRE API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    params = urllib.parse.urlencode({"q": f"doi:{doi}", "size": 1})
    try:
        return get_json("https://inspirehep.net/api/literature?" + params)["hits"]["hits"][0]["metadata"]
    except Exception:
        return {}


def metadata_from_openAlex(doi):
    """ get metadata for a DOI from the OpenAlex API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.openalex.org/works/https://doi.org/" + urllib.parse.quote(doi, safe="/")
    try:
        return get_json(url)
    except Exception:
        return {}


def clean_abstract(text):
    if not text:
        return None
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text)
    text = re.sub(r"^\s*Abstract\s+", "", " ".join(text.split()), flags=re.I)
    return text or None


def get_abstractdata(doi):
    """ get metadata and abstract for a DOI.
    """
    data_crossref = metadata_from_CrossRef(doi)
    data_inspire = metadata_from_Inspire(doi)
    data_openalex = metadata_from_openAlex(doi)
    abstract_crossref = clean_abstract(data_crossref.get("abstract"))
    abstract_inspire = clean_abstract(data_inspire.get("abstracts")[0].get("value") if data_inspire.get("abstracts") else None)
    abstract_openalex = clean_abstract(" ".join(word for _, word in sorted(
        (i, word) for word, positions in (data_openalex.get("abstract_inverted_index") or {}).items() for i in positions
    )))
    author_corresponding_openalex = [
        {
            "raw_author_name": authorship.get("raw_author_name") or (authorship.get("author") or {}).get("display_name"),
            "raw_affiliation_strings": authorship.get("raw_affiliation_strings"),
            "raw_orcid": authorship.get("raw_orcid") or (authorship.get("author") or {}).get("orcid"),
            "author_position": authorship.get("author_position"),
            "is_corresponding": authorship.get("is_corresponding"),
        }
        for authorship in data_openalex.get("authorships") or []
        if authorship.get("is_corresponding")
    ]
    return {
        "source": data_crossref["source"],
        "language": data_crossref.get("language"),
        "title": data_crossref.get("title"),
        "author": data_crossref.get("author"),
        "author-corresponding-openalex": author_corresponding_openalex,
        "DOI": data_crossref.get("DOI"),
        "URL": data_crossref.get("URL"),
        "is-referenced-by-count": data_crossref.get("is-referenced-by-count"),
        "reference-count": data_crossref.get("reference-count"),
        "abstract": abstract_crossref or abstract_inspire or abstract_openalex,
        "publisher": data_crossref.get("publisher"),
        "container-title": data_crossref.get("container-title"),
        "short-container-title": data_crossref.get("short-container-title"),
        "type": data_crossref.get("type"),
        "volume": data_crossref.get("volume"),
        "page": data_crossref.get("page"),
        "article-number": data_crossref.get("article-number"),
        "issued": data_crossref.get("issued"),
        "ISSN": data_crossref.get("ISSN"),
        "funder": data_crossref.get("funder"),
        "reference": data_crossref.get("reference"),   
    }


def generate_name(abstractdata):
    """ generate a file name for the DOI based on its metadata.
    """
    # Define a function to create a normalized key for name comparison
    name_key = lambda name: re.sub(r"\W+", "", name or "").lower()
    get_full_name = lambda author: author.get("name") or " ".join(filter(None, [author.get("given"), author.get("family")]))

    # get the corresponding author's name from OpenAlex, or fall back to the first author in CrossRef
    author_name_raw = ((abstractdata.get("author-corresponding-openalex") or [{}])[0]).get("raw_author_name")
    author_target = (abstractdata.get("author") or [{}])[0]
    for author_target in abstractdata.get("author") or []:
        if name_key(get_full_name(author_target)) == name_key(author_name_raw):
            author_target = author_target
            break
    author_name = get_full_name(author_target) or author_name_raw
    
    journal = abstractdata.get("short-container-title")[0]
    volume = abstractdata.get("volume")
    number = abstractdata.get("article-number") or ""
    page = abstractdata.get("page") or ""
    date_parts = (abstractdata.get("issued") or {}).get("date-parts") or [[]]
    year = (date_parts[0] or [""])[0]
    return f"{author_name.replace(' ', '_')}_Y.{year}_{journal.replace(' ', '')}_Vol.{volume}Nol.{number}P.{page}"

if __name__ == "__main__":
    doi = sys.argv[1]
    abstractdata = get_abstractdata(doi)
    abstractdata["filename"] = generate_name(abstractdata)

    raw_path = os.path.join(os.path.dirname(__file__), "..", "raw")
    os.makedirs(raw_path, exist_ok=True)
    json_path = os.path.join(raw_path, abstractdata["filename"] + ".json")

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(abstractdata, f, ensure_ascii=False, indent=2)
