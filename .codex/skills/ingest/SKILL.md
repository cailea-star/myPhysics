---
name: ingest
description: Ingest one paper into myWIKI. Use when user asks to add a paper by DOI, title, PDF, TEX, or JSON; generate raw files; collect full text; discuss tagged quotations.
---

# Ingest

## Rules

### Raw-Rules

`raw` is truth.

`raw/*.json` is metadata.

`raw/*.pdf` and `raw/*.tex` are original full text from network or user.

`raw/*.md` is quotation from original full text with tags.

### Tag-Rules

Scope: Quotation tags are centered on physical quantities actually calculated, plotted, or compared in the paper; do not tag unrelated background mentions.

Types: Use tag types from [vocab/types.json](../../../vocab/types.json), tags from [vocab/tags.json](../../../vocab/tags.json), and authors from [vocab/authors.json](../../../vocab/authors.json).

Type Match: Every type assigned to a drafted tag MUST independently satisfy that tag type's `requirement` in [vocab/types.json](../../../vocab/types.json); otherwise remove the type or reject the draft.

Granularity: Keep noun-term compounds as precise tags, but split adjective-like modifiers into property tags; e.g. use `alpha_decay_energy`, but use `symmetry_energy` + `soft`.

Draft: If a needed tag is missing, draft `{tag, definition, types, aliases}` entries for [vocab/tags.json](../../../vocab/tags.json) one tag type at a time; include only aliases seen in the paper, metadata, or existing project vocabulary.

Similarity Check: Before drafting each new tag, run `python scripts\search_similar_tags.py "CANDIDATE" 3` and show all three results. Prefer reusing existing tags.

Write: Present one draft group with an explicit add, merge, or reuse proposal for each candidate. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Write only individually approved entries, then present the next draft group.

### Claim-Type-Rules

Types: Use only claim types from [vocab/types.json](../../../vocab/types.json). A quotation satisfies a claim type ONLY when it meets that type's `requirement`.

Match: Assign exactly one `[claim_type]` from the quotation's explicit primary claim. Keywords, tags, source section, and Coverage targets are not classification evidence.

Coverage: A required claim type is covered ONLY by a quotation that passes its requirement. If no valid quotation exists, report `gap`; NEVER relabel other evidence to satisfy Coverage.

Review: An unsupported or incorrect `[claim_type]` is `fix` and blocks section approval.

### Template-Rules

Template Integrity: Every quotation section MUST contain exactly one `claim-types`, one `coverage`, and one `quotation` declaration. Missing or duplicate declarations block drafting and review.

Template Claim Types: Before drafting or reviewing a section, read its `claim-types` declaration from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md). Give every declared requirement exactly one `pass`, `gap`, or `fix`. Missing evidence is `gap`; incorrectly classified evidence is `fix`. NEVER relabel evidence to satisfy the declaration.

Template Coverage: Before drafting or reviewing a section, read its `coverage` declaration from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md). Give every requirement exactly one `pass`, `gap`, or `fix`. Missing valid evidence is `gap`; violating an exclusion is `fix`. NEVER invent evidence.

Template Quotation: Read the section's `quotation` declaration from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md). Any violation is `fix` and blocks section approval.

Counting: `[and]` separates independently reviewed items; `[or]` joins alternatives within one item. A `for each` item expands to one verdict per identified target. `optional` and `none` contribute zero items, and their absence is NEVER `gap`.

### Section-Rules

Section: Discuss quotation sections from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) in template order. Present at most four candidates or fixes per response, repeating batches until the section passes, in [exact mode] (final raw/*.md block shape with exact quoted sentence(s), any required `math` block, `[claim_type]`, `[tags]`, and source section) or [summary mode] (source section, candidate `[claim_type]`, candidate `[tags]`, and one-sentence evidence summary before the checkpoint).

Secondary-Citations-Check:
- If a quotation contains an explicit external citation marker, such as `[55]` or `Ref. [45]`, it MUST be classified as a citation-marked quotation for this check.
- If the core information of a citation-marked quotation comes from the external cited work, it MUST be placed under `### Secondary Citations`.
- If the core information of a citation-marked quotation states the current paper's own problem setting, method use, calculated result, or conclusion judgment, it MUST remain eligible for `Motivation`, `Methods`, `Results`, or `Meanings`.

Review: Independently review exactly one written quotation section per response in [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) order. First print its `claim-types`, `coverage`, and `quotation` declarations. Give every quotation and every required template item exactly one `pass`, `gap`, or `fix` verdict. Report `Review Verdicts: <verdicts>/<quotation blocks + required template items>`; a count mismatch fails the section and blocks writing.

Review Pass: A section passes ONLY when every quotation and every declared claim-type, coverage, and quotation requirement passes. Any `fix` blocks advancement. A `gap` passes ONLY when no valid source evidence exists and the gap is explicitly recorded.

Re-review: After ANY write to `raw/*.md`, run `python scripts\sort_raw_md_quotations.py mdfile_path`, then re-review the ENTIRE current section. Do not advance until it passes.

### Quotation-Rules

Source: Quote must come from `raw/*.pdf` or `raw/*.tex`.

Form: Each quote must be one or more complete sentences, not a phrase; it should support a clear claim and contain about 10-40 words total.

Tag Co-occurrence: Each quotation must include at least two directly supported comma-separated tags in `[tags]: ...`. If one sentence supports only one tag, add a neighboring complete sentence from the same source section; if no second supported tag exists, report the gap instead of inventing a relation.

Math: Do not use standalone formulas as quote text. When a formula is important core `[claim_type]: definition` evidence, quote the complete explanatory sentence and add the formula in a following fenced `math` block.

Write: Present quotation drafts one section at a time. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Write each approved section to `raw/*.md` in [exact mode] matching [scripts/add_raw_md.md](../../../scripts/add_raw_md.md). Formatting is defective only if it violates that template or breaks parsing.

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
   Check rough paper-level keywords against [vocab/tags.json](../../../vocab/tags.json); draft and confirm missing tags following [Tag-Rules](#tag-rules).

### Gate 4 — Discuss Quotations
   Draft and discuss exactly one quotation section per response, in [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) order, following [Claim-Type-Rules](#claim-type-rules), [Template-Rules](#template-rules), [Section-Rules](#section-rules), and [Quotation-Rules](#quotation-rules). First print its three template declarations and requirement verdicts, then present one review batch. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Write only approved exact-mode quotations. Complete every quotation section before Gate 5.

### Gate 5 — Review Discuss Quotations
   Apply [Claim-Type-Rules](#claim-type-rules), [Template-Rules](#template-rules), [Section-Rules](#section-rules), and [Quotation-Rules](#quotation-rules) to exactly one written quotation section per response in [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) order; present unresolved `fix` drafts in exact-mode batches of at most four. **🔴 CHECKPOINT · 🛑 STOP** — Write only explicitly approved drafts, then apply Re-review.

### Gate 6 — Summary & Recommend Next Paper(s)
   Log the completed paper before any recommendation: append one concise entry to [log.md](../../../log.md) with raw md filename, DOI, title, and core tags.
   Summarize from the completed `raw/*.md`: give exactly two sentences each for Motivation, Methods, Results, and Meanings, then state the paper's core innovation.
   Recommend one or more next papers from `### Secondary Citations` only: select cited references most central to the current paper's core tags and quotations, present each recommendation's citation information, DOI, matched tags, and why it is next. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.
