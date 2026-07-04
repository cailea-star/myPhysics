# myWIKI

Evidence-first literature wiki for nuclear physics.

## Structure

```text
raw/                 # source evidence
  *.json             # DOI metadata
  *.pdf              # original full text from network or user
  *.tex              # original source text from network or user
  *.md               # tagged quotations from the original full text

vocab/               # controlled vocabulary
  tags.json          # canonical tags
  types.json         # tag_types and claim_types
  authors.json       # author vocabulary

scripts/
  add_raw_json.py    # DOI -> raw/*.json
  add_raw_md.py      # raw/*.json -> raw/*.md
  add_raw_md.md      # raw md template
```

`raw` is truth. `vocab` names the vocabulary used to tag that truth.

## Ingest

Use the project `mywiki-ingest` skill for the ingest workflow.

```powershell
python scripts\add_raw_json.py [doi_number]
python scripts\add_raw_md.py raw\[json_filename].json
```

Put matching `*.pdf` or `*.tex` files in `raw/` with the same basename as the generated `*.json`.

## Raw Markdown

`raw/*.md` records quotations from the original `raw/*.pdf` or `raw/*.tex`.

Each quotation block uses:

````md
```tags
[claim_type]: ...
[tags]: ...
```

``` quote
quoted sentence from the original paper
```

source : section name.
````

`claim_type` must come from `vocab/types.json`.

Tags must come from `vocab/tags.json`.

Each quote needs a source section.

## Tags

Each tag has:

```json
{
  "tag": "tag_name",
  "definition": "definition",
  "types": ["quantity"],
  "aliases": ["seen alias"]
}
```

Keep noun-term compounds as precise tags, e.g. `alpha_decay_energy`.

Split adjective-like modifiers into property tags, e.g. `symmetry_energy` + `soft`.
