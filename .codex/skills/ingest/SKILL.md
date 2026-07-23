---
name: ingest
description: Ingest one paper into myWIKI. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# Ingest

## Rules

Before Gate 1, read both [Vocab Rules](../vocab-rules/SKILL.md) and [Raw Rules](../raw-rules/SKILL.md) completely and apply all rules.

## Gated Workflow

Run gates strictly in order. At the start of each response, state the current gate, the last completed gate, and the next required user approval. Do not advance past a gate until its listed work and required review are complete.

### Gate 1 — Confirm Paper

- Check: Run `git rev-parse --is-inside-work-tree` and `git status --short`; dirty worktrees MUST STOP and report changes.
- Identify exactly one paper from the user’s DOI, title, PDF, TEX, or JSON.
- DOI: Run `python scripts\search_a_doi.py [doi_number]`; report `recorded:` reuse/re-ingest options and reference-only hits.
- Approval: PRINT the identified paper; STOP until explicit user approval.

### Gate 2 — Generate Raw Files / Collect Full Text

- Generate: Run `python scripts\add_raw_json.py [doi_number]`, then `python scripts\add_raw_md.py raw\[json_filename].json`.
- Identify: Run `git diff --name-only -- raw/*.json`; use `git status --short raw` for untracked JSON.
- Full text: Download ONLY from arXiv; save same-basename `raw/[json_basename].pdf` or `.tex`.
- Failure: If arXiv fails, NEVER use other sources; PRINT direct links and target path; STOP.
- Complete: JSON, Markdown, and same-basename PDF-or-TEX MUST exist before Gate 3.

### Gate 3 — Check Tag & Author

- Author: Apply [Vocab Author-Rules](../vocab-rules/SKILL.md#author-rules) to the current Raw JSON.
- Identify: Apply [Motivation coverage](../../../scripts/add_raw_md.md#motivation); PRINT exactly three core Abstract concepts.
- Fix: Keep all three for later Abstract coverage; NEVER substitute easier existing tags.
- Resolve: Process each concept one-by-one through [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules); completion requires verified canonical tags.

### Gate 4 — Discuss Quotations

- Draft: Apply [Raw Draft-Rules](../raw-rules/SKILL.md#draft-rules) until every quotation section completes.
- Advance: Gate 5 begins ONLY after all quotation sections pass Raw verification.

### Gate 5 — Confirm Wiki

- PRINT: Show the three fixed core TAGs from Abstract and whether each `wiki/TAG.md` exists.
- Approval: STOP and ask whether to create, update, or skip each Wiki.
- Route: Approved Wiki work MUST follow [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) before Gate 6; skipped work advances directly.

### Gate 6 — Summary & Recommend Next Paper(s)

- Log: Before recommendations, append Raw filename, DOI, title, and core tags to [log.md](../../../log.md).
- Summarize: Using completed Raw, give two-sentence Motivation, Methods, Results, and Meanings summaries; state core innovation.
- Recommend: PRINT at least three central Secondary Citations with DOI, tags, citation, and rationale.
- Approval: STOP until explicit user approval before selecting another paper.
