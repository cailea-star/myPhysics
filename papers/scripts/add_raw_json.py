import html
import json
from pathlib import Path
import re
import sys
import urllib.parse
import urllib.request

HEADERS = {"User-Agent": "doi-metadata-script/0.1 (mailto:your@email.com)"}
PAPERS_PATH = Path(__file__).resolve().parent.parent
RAW_PATH = PAPERS_PATH.joinpath("raw")


def get_json(url):
    req = urllib.request.Request(url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)


def load_metadata_from_CrossRef(doi):
    """ get metadata for a DOI from the CrossRef API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="")
    return get_json(url)["message"]


def load_metadata_from_Inspire(doi):
    """ get metadata for a DOI from the INSPIRE API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    params = urllib.parse.urlencode({"q": f"doi:{doi}", "size": 1})
    try:
        return get_json("https://inspirehep.net/api/literature?" + params)["hits"]["hits"][0]["metadata"]
    except Exception:
        return {}


def load_metadata_from_openAlex(doi):
    """ get metadata for a DOI from the OpenAlex API.
    """
    doi = doi.strip().removeprefix("https://doi.org/")
    url = "https://api.openalex.org/works/https://doi.org/" + urllib.parse.quote(doi, safe="/")
    try:
        return get_json(url)
    except Exception:
        return {}


def clean_text(my_str):
    """Strip markup from metadata text.

    >>> clean_text("Neutron Radius of<mml:mi>P</mml:mi><mml:mn>208</mml:mn><mml:mi>b</mml:mi>")
    'Neutron Radius of 208Pb'
    """
    if not my_str:
        return None
    my_str = re.sub(r"<[^>]+>", " ", my_str)
    my_str = html.unescape(my_str)
    my_str = re.sub(r"\s+", " ", my_str).strip()
    my_str = re.sub(r"\b([A-Z])\s+(\d+)\s+([a-z])\b", r"\2\1\3", my_str)
    return my_str or None


def clean_abstract(abstract_str):
    abstract_str = clean_text(abstract_str)
    if not abstract_str:
        return None
    abstract_str = re.sub(r"^\s*Abstract\s+", "", abstract_str, flags=re.I)
    return abstract_str or None


def clean_titles(titles_str_list):
    cleaned_titles = []
    for title in titles_str_list or [""]:
        cleaned_titles.append(clean_text(title) or "")
    return cleaned_titles


def resolve_author_corresponding_from_openAlex(metadata_openAlex, metadata_CrossRef):
    name_key = lambda name: re.sub(r"\W+", "", name or "").lower()
    get_CrossRef_authorname = lambda author: " ".join(filter(None, [author.get("given"), author.get("family")]))
    author_corresponding_openalex = []
    for author_openAlex in metadata_openAlex.get("authorships") or []:
        # only include corresponding authors from OpenAlex
        if not author_openAlex.get("is_corresponding"): continue
        # find the corresponding author in CrossRef by matching names
        authorname_openAlex = author_openAlex.get("raw_author_name")
        author_CrossRef_target = (metadata_CrossRef.get("author") or [{}])[0]
        for author_CrossRef in metadata_CrossRef.get("author") or []:
            authorname_CrossRef = get_CrossRef_authorname(author_CrossRef)
            if name_key(authorname_CrossRef) == name_key(authorname_openAlex):
                author_CrossRef_target = author_CrossRef
                break
        # append the corresponding author information to the list
        author_corresponding_openalex.append({
            "raw_author_name": authorname_openAlex,
            "given": author_CrossRef_target.get("given") or "",
            "family": author_CrossRef_target.get("family") or "",
            "raw_orcid": author_openAlex.get("raw_orcid"),
            "ORCID": author_CrossRef_target.get("ORCID") or "",
            "raw_affiliation_strings": author_openAlex.get("raw_affiliation_strings"),
            "is_corresponding": author_openAlex.get("is_corresponding"),
        })
    return author_corresponding_openalex


def load_abstractdata(doi):
    """ get metadata and abstract for a DOI.
    """
    metadata_CrossRef   = load_metadata_from_CrossRef(doi)
    metadata_Inspire    = load_metadata_from_Inspire(doi)
    metadata_openAlex   = load_metadata_from_openAlex(doi)
    abstract_crossRef   = clean_abstract(metadata_CrossRef.get("abstract"))
    abstract_Inspire    = clean_abstract(metadata_Inspire.get("abstracts")[0].get("value") if metadata_Inspire.get("abstracts") else None)
    abstract_openAlex   = clean_abstract(" ".join(word for _, word in sorted(
        (i, word) for word, positions in (metadata_openAlex.get("abstract_inverted_index") or {}).items() for i in positions
    )))
    author_corresponding = resolve_author_corresponding_from_openAlex(metadata_openAlex, metadata_CrossRef)
    return {
        "source": metadata_CrossRef["source"],
        "language": metadata_CrossRef.get("language"),
        "title": clean_titles(metadata_CrossRef.get("title")),
        "author": metadata_CrossRef.get("author"),
        "author-corresponding-openalex": author_corresponding,
        "DOI": metadata_CrossRef.get("DOI"),
        "URL": metadata_CrossRef.get("URL"),
        "is-referenced-by-count": metadata_CrossRef.get("is-referenced-by-count"),
        "reference-count": metadata_CrossRef.get("reference-count"),
        "abstract": abstract_crossRef or abstract_Inspire or abstract_openAlex,
        "publisher": metadata_CrossRef.get("publisher"),
        "container-title": metadata_CrossRef.get("container-title"),
        "short-container-title": metadata_CrossRef.get("short-container-title"),
        "type": metadata_CrossRef.get("type"),
        "volume": metadata_CrossRef.get("volume"),
        "page": metadata_CrossRef.get("page"),
        "article-number": metadata_CrossRef.get("article-number"),
        "issued": metadata_CrossRef.get("issued"),
        "ISSN": metadata_CrossRef.get("ISSN"),
    }


def generate_stem(abstractdata):
    """ generate a file stem for the DOI based on its metadata.
    """
    author = (abstractdata.get("author-corresponding-openalex") or abstractdata.get("author") or [{}])[0]
    author_name = " ".join(filter(None, [author.get("given"), author.get("family")])) or author.get("raw_author_name") or ""
    
    journal = abstractdata.get("short-container-title")[0]
    journal = journal.replace(":", "_")
    volume = abstractdata.get("volume")
    number = abstractdata.get("article-number") or ""
    page = abstractdata.get("page") or ""
    date_parts = (abstractdata.get("issued") or {}).get("date-parts") or [[]]
    year = (date_parts[0] or [""])[0]
    return f"{author_name.replace(' ', '_')}_Y.{year}_{journal.replace(' ', '')}_Vol.{volume}Nol.{number}P.{page}"

if __name__ == "__main__":
    doi = sys.argv[1]
    abstractdata = load_abstractdata(doi)
    abstractdata["stem"] = generate_stem(abstractdata)

    RAW_PATH.mkdir(exist_ok=True)
    json_path = RAW_PATH.joinpath(f"{abstractdata['stem']}.json")
    json_path.write_text(json.dumps(abstractdata, ensure_ascii=False, indent=2), encoding="utf-8")
