# Papers

Evidence-first literature workspace for nuclear physics. `raw/` stores source evidence, `vocab/` defines the canonical vocabulary, and `wiki/` synthesizes verified quotations.

```text
papers/
├── raw/       # metadata, full text, and tagged quotations
├── vocab/     # tags, types, and authors
├── wiki/      # synthesized topic and method pages
├── scripts/   # ingestion, search, validation, and sorting tools
├── tmp/       # generated query results
├── log.md     # completed-paper log
└── README.md
```

## Scripts

Run from the repository root:

```powershell
python papers\scripts\add_raw_json.py DOI                                  # Fetch DOI metadata and create Raw JSON
python papers\scripts\add_raw_md.py papers\raw\paper.json                  # Create Raw Markdown from JSON metadata
python papers\scripts\add_wiki_md.py TAG                                   # Create a Wiki page from a canonical tag
python papers\scripts\check_vocab_author.py papers\raw\paper.json         # Check and add corresponding authors
python papers\scripts\rename_raw_tag.py OLD NEW                            # Replace a merged tag in Raw Markdown
python papers\scripts\search_a_author.py GIVEN FAMILY                      # Find papers by author name
python papers\scripts\search_a_paper.py DOI                                # Find recorded papers and citations by DOI
python papers\scripts\search_a_tag.py TAG [TAG ...]                         # Collect tagged quotations into tmp
python papers\scripts\search_similar_tags.py QUERY NUMBER                  # Rank canonical tags by similarity
python papers\scripts\sort_raw_quotations.py RAW_MD_PATH                   # Renumber Raw quotation headings
python papers\scripts\sort_wiki_studies.py WIKI_MD_PATH                    # Sort Previous Studies by year

# Templates
# papers\scripts\add_raw_md.md                                               Raw evidence and quotation template
# papers\scripts\add_wiki_method.md                                         Method Wiki template
# papers\scripts\add_wiki_topic.md                                          Topic and phenomenon Wiki template
```

Use [papers-ingest](../.codex/skills/papers-ingest/SKILL.md) to add papers and [papers-lint](../.codex/skills/papers-lint/SKILL.md) to audit Wiki coverage. The `papers-*-rules` skills define the canonical Raw, vocabulary, and Wiki constraints.
