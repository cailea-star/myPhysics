---
name: ingest
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

Similarity Check: Before drafting any new tag, check `vocab/tags.json` for semantic duplicates or related existing tags. For each candidate, include `similar-tags: ...`. Prefer reusing existing tags. The user must decide each candidate individually: add new, merge into existing, or use existing. Never write new tags before that decision.

Write: Present one draft group, wait for user approval, write confirmed entries, then present the next draft group.

### Section-Rules

Section: Discuss quotations under the template frames: Motivation, Methods, Results, Meanings. For each section, present three to four candidate quotations for review in [exact mode] (final raw/*.md block shape with exact quoted sentence(s), any required `math` block, `[claim_type]`, `[tags]`, and source section) or [summary mode] (source section, candidate `[claim_type]`, candidate `[tags]`, and one-sentence evidence summary before asking for approval).

Coverage: For each section's candidate set, satisfy section-level claim coverage from source evidence; if no valid source quote exists for a required claim, report the gap instead of inventing one.

Motivation-Coverage:
- Include at least one `[claim_type]: motivation`.
- Include `[claim_type]: definition` for core calculated quantity tags.
- Include `[claim_type]: comparison` or `[claim_type]: motivation` for core method tags.

Methods-Coverage:
- Include `[claim_type]: definition` for core method tags.
- Include `[claim_type]: comparison` for core method tags.
- Check whether the paper states an explicit `[claim_type]: assumption`.

Results-Coverage:
- Cover every figure and table: present >=1 valid complete-sentence quotation candidate for each item with its calculated or measured quantity tags, or report that item as a gap.
- Include `[claim_type]: definition` for core calculated or measured quantity tags.
- Include `[claim_type]: comparison` for core calculated or measured quantity tags.

Results-Quotation:
- Every quotation must include at least one directly supported calculation-method tag and at least one calculated or measured quantity tag.
- If either tag type is unsupported, reject the candidate.

Meanings-Coverage:
- Meanings may repeat Results, but ONLY for claims that state the paper's central significance.
- MUST cover high-level judgments from the abstract, summary, or conclusion when supported by valid source quotes.
- Reject minor results, figure-only details, and claims that do not affect the paper's main message.
- Include at least one `[claim_type]: comparison`.

Meanings-Quotation:
- Every quotation must include at least one directly supported method tag from `vocab/tags.json`.
- If no method tag is supported, reject the candidate.

Secondary-Citations-Coverage:
- Include secondary citations ONLY for external methods/formulas, parameter sets/models/data sources, or key physical conclusions/mechanism judgments that support core tags.
- Reject generic survey/background citations and citations that do not supply a reused method, formula, model/data input, or core interpretation.

Secondary-Citations-Quotation:
- Put a quote under `### Secondary Citations` ONLY when the quoted sentence depends on an external cited reference; use `secondary-tags`.
- Each quote MUST be a complete sentence with an explicit external citation marker, such as `[55]` or `Ref. [45]`; otherwise reject it.
- Do not use secondary quotes as primary evidence for Motivation, Methods, Results, or Meanings.
- Internal references to this paper's `Fig.`, `Table`, `Eq.`, or `section` are not secondary citations.

### Quotation-Rules

Source: Quote must come from `raw/*.pdf` or `raw/*.tex`.

Form: Each quote must be one or more complete sentences, not a phrase; it should support a clear claim and contain about 10-40 words total.

Tag Co-occurrence: Each quotation must include at least two directly supported comma-separated tags in `[tags]: ...`. If one sentence supports only one tag, add a neighboring complete sentence from the same source section; if no second supported tag exists, report the gap instead of inventing a relation.

Math: Do not use standalone formulas as quote text. When a formula is important core `[claim_type]: definition` evidence, quote the complete explanatory sentence and add the formula in a following fenced `math` block.

Write: Present quotation drafts one section at a time; wait for user review, then write each approved section to `raw/*.md` in [exact mode].

## Gated Workflow

Run gates strictly in order. At the start of each response, state the current gate, the last completed gate, and the next required user approval. Do not advance past a gate until its listed work and required review are complete.

### Gate 1 — Confirm Paper
   Check that the current directory is inside a git repository with `git rev-parse --is-inside-work-tree`.
   Check that the worktree is clean with `git status --short`; if it is not empty, stop and report the existing changes before ingest.
   User gives DOI, title, PDF, TEX, or JSON. Confirm exact paper before ingest.

### Gate 2 — Generate Raw Files / Collect Full Text
   ```powershell
   python scripts\add_raw_json.py [doi_number]
   python scripts\add_raw_md.py raw\[json_filename].json
   ```
   Use `git diff --name-only -- raw/*.json` to identify the new JSON filename; use `git status --short raw` if it is untracked.

   Get PDF/TEX from open source or user. Save into `raw/` with same basename as JSON.

### Gate 3 — Check Tag & Author
   Run `python scripts\add_vocab_author.py raw\[json_filename].json` as a routine check of the corresponding-author list, then read the PDF/TEX source text for corresponding-author information and report both the script terminal output and the source-text corresponding-author information to the user.
   Check rough paper-level keywords against `vocab/tags.json`; draft and confirm missing tags following [Vocab-Rules](#vocab-rules).

### Gate 4 — Discuss Quotations
   Draft and write quotations under the template frames, following [Section-Rules](#section-rules) and [Quotation-Rules](#quotation-rules).

### Gate 5 — Review Discuss Quotations
   Review the completed `raw/*.md` one section at a time in template-frame order: Motivation, Methods, Results, Meanings.
   For each section, first print that section's current Coverage rules and any section-specific Quotation rules from [Section-Rules](#section-rules).
   Review each section against [Section-Rules](#section-rules) and [Quotation-Rules](#quotation-rules); report pass/gap/fix items before moving to the next section.

### Gate 6 — Recommend Next Paper
   Read `### Secondary Citations` in the completed `raw/*.md`.
   Append a concise completed-paper entry to `log.md` before recommending next papers: raw md filename, DOI, title, and core tags.
   Recommend next-paper candidates from cited references that are most central to the current paper's tags and quotations.
   Present candidates with citation information and why they are next; wait for user approval before starting a new ingest.
