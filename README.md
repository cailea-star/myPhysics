# myWIKI

Evidence-first literature wiki for nuclear physics.

## Structure

```text
raw/        # source evidence
  *.json    # DOI metadata
  *.pdf     # original paper, optional
  *.tex     # source text, optional
  *.md      # tagged original quotes, optional

taglist/   # controlled vocabulary
  entities.json
  conditions.json
  phenomena.json
  mechanisms.json
  methods.json
  systems.json
  quantities.json
  claims.json

scripts/    # ingest/query/lint tools
  add_raw_json.py

template/   # file templates
  paper.md
  taglist.json
```

`raw` is truth. `taglist` is vocabulary. Wiki pages are generated views, not maintained truth.

## Operations

**Ingest**

DOI/PDF/TEX -> `raw`: metadata, stable filename, abstract, corresponding author, source path, tagged quotes.

**Query**

Keyword -> canonical term in `taglist` -> tagged evidence in `raw` -> assembled answer/wiki view.

**Lint**

Check missing files, invalid tags, duplicate terms, stale metadata, uncited claims.

## Templates

Paper quote format lives in `template/paper.md`.

Taglist format lives in `template/taglist.json`.