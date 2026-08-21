---
name: papers-ingest
description: Ingest one paper into myWIKI. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# Papers Ingest

## Rules

Before Gate 1, MUST read and follow both [Vocab Rules](../papers-vocab-rules/SKILL.md) and [Raw Rules](../papers-raw-rules/SKILL.md) completely.

## Gated Workflow

Run gates strictly in order. At the start of each response, state the current gate, the last completed gate, and the next required user approval. Do not advance past a gate until its listed work and required review are complete.

### Gate 1 — Confirm Paper

- Check: Run `git rev-parse --is-inside-work-tree` and `git status --short`; dirty worktrees MUST STOP and report changes.
- Identify exactly one paper from the user’s DOI, title, PDF, TEX, or JSON.
- DOI: Run `python papers\scripts\search_a_paper.py [doi_number]`; report `recorded:` reuse/re-ingest options and `md secondary:` hits.
- Approval: PRINT the identified paper; STOP until explicit user approval.

### Gate 2 — Generate Raw Files / Collect Full Text

- Generate: Run `python papers\scripts\add_raw_json.py [doi_number]`, then `python papers\scripts\add_raw_md.py papers\raw\[json_filename].json`.
- Identify: Run `git diff --name-only -- papers/raw/*.json`; use `git status --short papers/raw` for untracked JSON.
- Full text: Download ONLY from arXiv; save same-basename `papers/raw/[json_basename].pdf` or `.tex`.
- Failure: If arXiv fails, NEVER use other sources; PRINT direct links and target path; STOP.
- Complete: JSON, Markdown, and same-basename PDF-or-TEX MUST exist before Gate 3.

### Gate 3 — Check Tag & Author

- Author: MUST execute every required step in [Vocab Author-Rules](../papers-vocab-rules/SKILL.md#author-rules) for the current Raw JSON.
- Identify: Apply [Motivation coverage](../../../papers/scripts/add_raw_md.md#motivation); PRINT exactly three core Abstract concepts.
- Fix: Keep all three for later Abstract coverage; NEVER substitute easier existing tags.
- Resolve: MUST execute every required step in [Vocab Draft-Rules](../papers-vocab-rules/SKILL.md#draft-rules) for each concept one-by-one; completion requires verified canonical tags.

### Gate 4 — Discuss Quotations

- Draft: MUST execute every required step in [Raw Draft-Rules](../papers-raw-rules/SKILL.md#draft-rules) until every quotation section completes.
- Advance: Gate 5 begins ONLY after all quotation sections pass Raw verification.

### Gate 5 — Confirm Wiki

- PRINT: Show the three fixed core TAGs from Abstract and whether each `papers/wiki/TAG.md` exists.
- Approval: STOP and ask whether to create, update, or skip each Wiki.
- Route: Approved Wiki work MUST execute every required step in [Wiki Draft-Rules](../papers-wiki-rules/SKILL.md#draft-rules) before Gate 6; skipped work advances directly.

### Gate 6 — Summary & Recommend Next Paper(s)

- Log: Before recommendations, append Raw filename, DOI, title, and core tags to [papers/log.md](../../../papers/log.md).
- Summarize: Using completed Raw, give two-sentence Motivation, Methods, Results, and Meanings summaries; state core innovation.
- Recommend: PRINT at least three central Secondary Citations with DOI, tags, citation, and rationale.
- Approval: STOP until explicit user approval before selecting another paper.
