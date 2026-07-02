# myWIKI

Evidence-first literature wiki for nuclear physics.

## Structure

```text
raw/        # source evidence
  *.json    # DOI metadata
  *.pdf     # original paper, optional
  *.tex     # source text, optional
  *.md      # tagged original quotes, optional

namelist/   # controlled vocabulary
  entities.json
  conditions.json
  phenomena.json
  mechanisms.json
  methods.json
  systems.json
  quantities.json
  claim_types.json

scripts/    # ingest/query/lint tools
  doi_meta.py
```

`raw` is truth. `namelist` is vocabulary. Wiki pages are generated views, not maintained truth.

## Operations

**Ingest**

DOI/PDF/TEX -> `raw`: metadata, stable filename, abstract, corresponding author, source path, tagged quotes.

**Query**

Keyword -> canonical term in `namelist` -> tagged evidence in `raw` -> assembled answer/wiki view.

**Lint**

Check missing files, invalid tags, duplicate terms, stale metadata, uncited claims.

## Tags

```text
entity:...
condition:...
phenomenon:...
mechanism:...
method:...
system:...
quantity:...
claim:...
```
