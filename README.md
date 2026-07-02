# myWIKI

Evidence-first literature wiki for nuclear physics.

## File Structure

```text
raw/
  *.json   # DOI metadata from scripts
  *.pdf    # original paper
  *.tex    # optional source
  *.md     # tagged original quotes, not summaries

namelist/
  topics.md
  methods.md
  authors.md
  relations.md

scripts/
  doi_meta.py

wiki_docs/
  generated/  # optional query output/cache
```

`raw` is source of truth. `namelist` is project vocabulary. `wiki_docs` is rebuildable output.

## Operations

**Ingest**

DOI/PDF/TEX -> `raw`. Extract metadata, stable filename, abstract, authors, corresponding author, source path. Later: add tagged quotes.

**Query**

Keyword -> `namelist` terms -> search tagged raw evidence -> assemble answer or temporary wiki page.

**Lint**

Check missing pairs, broken filenames, duplicate topics/authors, stale tags, uncited claims.

## Principle

Do not maintain loose wiki prose as truth. Keep raw evidence + namelist. Generate wiki views when needed.
