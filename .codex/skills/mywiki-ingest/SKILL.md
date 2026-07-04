---
name: mywiki-ingest
description: Ingest one paper into myWIKI. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# myWIKI Ingest

## Rules

### Raw-Rules

`raw` is truth.

`raw/*.json` is metadata.

`raw/*.pdf` and `raw/*.tex` are original full text from network or user.

`raw/*.md` is quotation from original full text with tags.

### Vocab-Rules

Scope: Quotation tags are centered on physical quantities actually calculated, plotted, or compared in the paper; do not tag unrelated background mentions.

Types: Use tag types and claim types from `vocab/types.json`, tags from `vocab/tags.json`, and authors from `vocab/authors.json`.

Granularity: Keep noun-term compounds as precise tags, but split adjective-like modifiers into property tags; e.g. use `alpha_decay_energy`, but use `symmetry_energy` + `soft`.

Draft: If a needed tag is missing, draft `{tag, definition, types, aliases}` entries for `vocab/tags.json` one tag type at a time; include only aliases seen in the paper, metadata, or existing project vocabulary.

Write: Present one draft group, wait for user approval, write confirmed entries, then present the next draft group.

### Quotation-Rules

Source: Quote must come from `raw/*.pdf` or `raw/*.tex`.

Form: Each quote must be a complete sentence, not a phrase; it should support a clear claim and contain about 10-40 words.

Section: Discuss quotations under the template frames: Motivation, Methods, Results, Meanings. Each quote needs source section; each section must contain two to four quotations.

Definition: Each core tag used for a paper must have at least one quotation with `[claim_type]: definition` and `[tags]: tag_name`.

Write: Present and write quotation drafts one section at a time; wait for user review before writing each section to `raw/*.md`.

## Workflow

1. Confirm paper
   Check that the current directory is inside a git repository with `git rev-parse --is-inside-work-tree`.
   Check that the worktree is clean with `git status --short`; if it is not empty, stop and report the existing changes before ingest.
   User gives DOI, title, PDF, TEX, or JSON. Confirm exact paper before ingest.

2. Generate raw files
   ```powershell
   python scripts\add_raw_json.py [doi_number]
   python scripts\add_raw_md.py raw\[json_filename].json
   ```
   Use `git diff --name-only -- raw/*.json` to identify the new JSON filename; use `git status --short raw` if it is untracked.

3. Collect full text
   Get PDF/TEX from open source or user. Save into `raw/` with same basename as JSON.

4. Check tag & author
   Run `python scripts\add_vocab_author.py raw\[json_filename].json` as a routine check of the corresponding-author list, then read the PDF/TEX source text for corresponding-author information and report both the script terminal output and the source-text corresponding-author information to the user.
   Check rough paper-level keywords against `vocab/tags.json`; draft and confirm missing tags following [Vocab-Rules](#vocab-rules).

5. Discuss quotations
   Draft and write quotations under the template frames, following [Quotation-Rules](#quotation-rules).
