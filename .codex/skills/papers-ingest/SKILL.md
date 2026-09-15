---
name: papers-ingest
description: Ingest one paper into the literature workspace. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# Papers Ingest

## Rules

- Before Gate 1, MUST read and follow [Vocab Rules](../papers-vocab-rules/SKILL.md) completely.
- Before Gate 1, MUST read and follow [Raw Rules](../papers-raw-rules/SKILL.md) completely.

## Gated Workflow

- MUST run gates strictly in order.
- Response opening: MUST state current gate, last completed gate, and next required user approval.
- Advance ONLY after the current gate's listed work and required review are complete.

### Gate 1 — Confirm Paper

#### Identify Paper

- Check: MUST run `git rev-parse --is-inside-work-tree` and `git status --short`.
- Dirty worktrees: MUST STOP and report changes.
- Identify: MUST identify exactly one paper from the user's DOI, title, PDF, TEX, or JSON.
- DOI: MUST run `python papers\scripts\search_a_paper.py [doi_number]`.
- Report: MUST report `recorded:` reuse/re-ingest options and `md secondary:` hits.

#### Confirm Paper

- Approval: MUST PRINT the identified paper; STOP until explicit user approval.

### Gate 2 — Collect Full Text / Confirm Understanding

#### Collect Full Text

- Full text: Download ONLY from arXiv; MUST save PDF or TEX temporarily in `papers/tmp/`.
- Failure: If arXiv fails, NEVER use other sources; PRINT direct links and target path; STOP.

#### Determine Core Tags

- Core tags: The reader MUST determine exactly three core tags.
- Changes: NEVER change these core tags without the reader's approval.

#### Confirm Understanding

- Focus: MUST organize the understanding around the reader's three core tags.
- Summarize: Using full text, MUST give two-sentence Motivation, Methods, Results, and Meanings summaries.
- Assessment: MUST state the paper's core innovation and gaps, distinguishing explicit statements from inferred limitations.
- Approval: MUST discuss until the reader explicitly confirms shared understanding before entering Gate 3.

### Gate 3 — Generate Raw Files / Check Tags & Authors

#### Generate Raw Files

- JSON: MUST run `python papers\scripts\add_raw_json.py [doi_number]`.
- Markdown: MUST then run `python papers\scripts\add_raw_md.py papers\raw\[json_filename].json`.
- Identify: MUST run `git diff --name-only -- papers/raw/*.json`.
- Untracked JSON: MUST use `git status --short papers/raw`.
- Archive: MUST copy the collected full text to same-basename `papers/raw/[json_basename].pdf` or `.tex`.

#### Check Tags & Authors

- Author: MUST execute every required [Vocab Author-Rules](../papers-vocab-rules/SKILL.md#author-rules) step for the current Raw JSON.
- Resolve: For each concept one-by-one, MUST execute every required [Vocab Draft-Rules](../papers-vocab-rules/SKILL.md#draft-rules) step.

#### Complete Ingestion Setup

- Tags: Completion MUST require verified canonical tags.
- Files: JSON, Markdown, and same-basename PDF-or-TEX MUST exist before Gate 4.

### Gate 4 — Discuss Quotations

- Draft: MUST execute every required [Raw Draft-Rules](../papers-raw-rules/SKILL.md#draft-rules) step until every quotation section completes.
- Advance: Gate 5 begins ONLY after all quotation sections pass Raw verification.

### Gate 5 — Confirm Wiki

- PRINT: MUST show three fixed core TAGs from the reader and whether each `papers/wiki/TAG.md` exists.
- Approval: MUST STOP and ask whether to create, update, or skip each Wiki.
- Route: Approved Wiki work MUST execute every required [Wiki Draft-Rules](../papers-wiki-rules/SKILL.md#draft-rules) step before Gate 6.
- Skip: Skipped Wiki work advances directly to Gate 6.

### Gate 6 — Log & Recommend Next Paper(s)

- Log: Before recommendations, MUST append Raw filename, DOI, title, and core tags to [papers/log.md](../../../papers/log.md).
- Recommend: MUST PRINT at least three central Secondary Citations with DOI, tags, citation, and rationale.
- Approval: MUST STOP until explicit user approval before selecting another paper.
