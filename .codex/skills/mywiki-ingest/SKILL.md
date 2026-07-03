---
name: mywiki-ingest
description: Ingest one paper into myWIKI. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# myWIKI Ingest

## Rules

`raw` is truth.

`raw/*.json` is metadata.

`raw/*.pdf` and `raw/*.tex` are original full text from network or user.

`raw/*.md` is quotation from original full text with tags.

Tags must come from `taglist/*.json`. Quote must come from `raw/*.pdf` or `raw/*.tex`. Each quote needs source section.

## Workflow

1. Confirm paper
   User gives DOI, title, PDF, TEX, or JSON. Confirm exact paper before ingest.

2. Generate raw files
   ```powershell
   python scripts\add_raw_json.py [doi_number]
   python scripts\add_raw_md.py raw\[json_filename].json
   ```

3. Collect full text
   Get PDF/TEX from open source or user. Save into `raw/` with same basename as JSON.

4. Discuss quotations
   Section by section, choose quote, claim type, tags, source section. Write into `raw/*.md`.
