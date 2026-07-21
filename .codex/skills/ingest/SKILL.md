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
   Check that the current directory is inside a git repository with `git rev-parse --is-inside-work-tree`.
   Check that the worktree is clean with `git status --short`; if it is not empty, stop and report the existing changes before ingest.
   Identify the exact paper from the user's DOI, title, PDF, TEX, or JSON.
   Known DOI: run `python scripts\search_a_doi.py [doi_number]` before Gate 2. On `recorded:`, report the match and reuse/re-ingest options; report reference-only hits as context.
   **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.

### Gate 2 — Generate Raw Files / Collect Full Text
   ```powershell
   python scripts\add_raw_json.py [doi_number]
   python scripts\add_raw_md.py raw\[json_filename].json
   ```
   Use `git diff --name-only -- raw/*.json` to identify the new JSON filename; use `git status --short raw` if it is untracked.
   Download PDF/TEX from arXiv only; save as `raw/[json_basename].pdf` or `.tex`.
   If arXiv fails, do not try other downloads. Give direct PDF link(s), the exact target path, and wait for the user file.

### Gate 3 — Check Tag & Author
   Run `python scripts\check_vocab_author.py raw\[json_filename].json` as a routine check of the corresponding-author list, then read the PDF/TEX source text for corresponding-author information and report both the script terminal output and the source-text corresponding-author information to the user.
   Select and print exactly three core Abstract concepts; keep them fixed for later abstract coverage and resolve each under [Vocab Resolve-Rules](../vocab-rules/SKILL.md#resolve-rules), applying Draft-Rules when needed—NEVER substitute easier existing tags.

### Gate 4 — Discuss Quotations
   Draft and discuss exactly one quotation section per response, in [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) order, following [Template-Rules](../raw-rules/SKILL.md#template-rules), [Draft-Rules](../raw-rules/SKILL.md#draft-rules), [Section-Rules](../raw-rules/SKILL.md#section-rules), and [Quotation-Rules](../raw-rules/SKILL.md#quotation-rules). Present one review batch. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Write only approved exact-mode quotations. Complete every quotation section before Gate 5.

### Gate 5 — Review Discuss Quotations
   Apply [Template-Rules](../raw-rules/SKILL.md#template-rules), [Draft-Rules](../raw-rules/SKILL.md#draft-rules), [Section-Rules](../raw-rules/SKILL.md#section-rules), and [Quotation-Rules](../raw-rules/SKILL.md#quotation-rules) to exactly one written quotation section per response in [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) order; present unresolved `fix` drafts in exact-mode batches of at most four. **🔴 CHECKPOINT · 🛑 STOP** — Write only explicitly approved drafts, then apply Re-review.

### Gate 6 — Summary & Recommend Next Paper(s)
   Log the completed paper before any recommendation: append one concise entry to [log.md](../../../log.md) with raw md filename, DOI, title, and core tags.
   Summarize from the completed `raw/*.md`: give exactly two sentences each for Motivation, Methods, Results, and Meanings, then state the paper's core innovation.
   Recommend one or more next papers from `### Secondary Citations` only: select cited references most central to the current paper's core tags and quotations, present each recommendation's citation information, DOI, matched tags, and why it is next. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.
