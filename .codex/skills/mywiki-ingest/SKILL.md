---
name: mywiki-ingest
description: Ingest one paper into myWIKI. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# myWIKI Ingest

## Rules

### Raw Rules

`raw` is truth.

`raw/*.json` is metadata.

`raw/*.pdf` and `raw/*.tex` are original full text from network or user.

`raw/*.md` is quotation from original full text with tags.

### Tag Rules

Tags must come from `taglist/*.json`.

For quotation tags, `quantities` must be physical quantities actually calculated, plotted, fitted, or numerically compared in the paper.

For quotation tags, `methods` must be theory, computational method, model, code, formula, or experimental method actually used to obtain results in the paper.

Before writing any new taglist entry, present its full `tag`, `description`, and `aliases` to the user for review.

Present new tag drafts one `taglist/*.json` category at a time; do not mix many tag categories in one review.

### Quote Rules

Quote must come from `raw/*.pdf` or `raw/*.tex`.

Each quote needs source section.

Each quote must be a complete sentence, not a phrase; it should support a clear claim and contain about 20-40 words.

Each section must contain at least two quotations.

Before adding new quotations in `raw/*.md`, present the draft and wait for user review.

## Workflow

1. Confirm paper
   User gives DOI, title, PDF, TEX, or JSON. Confirm exact paper before ingest.

2. Generate raw files
   ```powershell
   python scripts\add_raw_json.py [doi_number]
   python scripts\add_raw_md.py raw\[json_filename].json
   ```
   Use `git diff --name-only -- raw/*.json` to identify the new JSON filename; use `git status --short raw` if it is untracked.

3. Collect full text
   Get PDF/TEX from open source or user. Save into `raw/` with same basename as JSON.

4. Discuss quotations
   Discuss quotations under the template frames: Motivation, Methods, Results, Meanings. Section by section, choose quote, claim type, tags, source section.
   If a needed keyword is missing from `taglist`, draft a matching `taglist/*.json` entry using existing template. Present draft tag entries and quotations for user review before writing to `taglist/*.json` or `raw/*.md`.
