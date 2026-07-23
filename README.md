# myWIKI

Evidence-first literature wiki for nuclear physics.

```text
raw/                  # source evidence
  *.json              # DOI metadata
  *.pdf               # original full text from network or user
  *.tex               # original source text from network or user
  *.md                # tagged quotations from the original full text

vocab/                # controlled vocabulary
  tags.json           # canonical tags
  types.json          # tag_types and claim_types
  authors.json        # author vocabulary

scripts/
  add_raw_json.py           # DOI -> raw/*.json
  add_raw_md.py             # raw/*.json -> raw/*.md
  add_raw_md.md             # raw md template
  add_wiki_topic.md         # system and phenomenon wiki template
  add_wiki_method.md        # method wiki template
  check_vocab_author.py     # check and add author to vocab/authors.json
  search_a_author.py        # author -> matching raw metadata
  search_a_doi.py           # DOI -> recorded papers and citation links
  search_a_tag.py           # canonical tag -> tmp/* section quotation slices
  search_similar_tags.py    # query -> similar canonical tags
  sort_raw_quotations.py    # renumber raw/*.md quotation headings after edits
  sort_wiki_studies.py      # order Previous Studies by publication year

tmp/                  # generated temporary query outputs

.codex/skills/
  ingest/SKILL.md      # gated paper-ingest workflow

log.md                # completed paper log
log-PRL.md            # Physical Review Letters coverage log
```

`raw` is truth. `vocab` names the vocabulary. Use the project `ingest` skill for workflow rules.

Renumber a raw markdown file after adding or deleting quotations:

```powershell
python scripts\sort_raw_quotations.py raw\paper.md
```

Sort a wiki's `Previous Studies` after edits:

```powershell
python scripts\sort_wiki_studies.py wiki\topic.md
```
