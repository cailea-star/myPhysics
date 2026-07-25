---
name: raw-rules
description: Apply canonical myWIKI raw evidence rules when creating, drafting, reviewing, sorting, or validating raw files, quotation sections, claim types, evidence links, templates, or tagged quotations in ingest or lint workflows.
---

# Raw Rules

## Input-Rules

### Raw-Inputs

- `raw/*.md` is evidence truth; `raw/*.json` is metadata.
- Same-basename `raw/*.pdf` or `raw/*.tex` is the primary source.
- PRINT: Missing, duplicate, or mismatched Raw inputs require STOP.

### Vocab-Inputs

- [vocab/tags.json](../../../vocab/tags.json) and [vocab/types.json](../../../vocab/types.json) are auxiliary inputs.
- [vocab/tags.json](../../../vocab/tags.json) defines canonical tags, aliases, definitions, and tag types.
- [vocab/types.json](../../../vocab/types.json) defines claim-type and tag-type requirements.
- Run `python scripts\search_similar_tags.py QUERY 3` for every proposed tag.
- PRINT: Report every unresolved proposed tag before input validation.

### Check-Inputs

- Scope: MUST validate inputs against [Raw-Inputs](#raw-inputs) and [Vocab-Inputs](#vocab-inputs) before Checkpoints.
- Checkpoint 1: MUST run `python scripts\search_a_paper.py DOI`; ONLY one matching JSON with same-basename PDF/TEX passes.
- Checkpoint 2: ONLY one canonical match passes; failure MUST STOP and execute [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules) completely.
- Checkpoint 3: Missing Markdown MUST enter create mode; existing Markdown passes ONLY with matching basename and `[source_json]`.
- PRINT: Any failure MUST STOP; report solutions; NEVER overwrite mismatched Markdown.


## Template-Rules

- Template: [scripts/add_raw_md.md](../../../scripts/add_raw_md.md); manual: `Motivation`, `Methods`, `Results`, `Meanings`, `Secondary Citations`; auxiliary: `Abstract`, `Gaps`.

### Section-Rules

- Each quotation section in the template MUST contain exactly one `claim-type`, `coverage`, and `quotation` HTML declaration.
- Collective section evidence MUST yield exactly one `pass`, `gap`, or `fix` per applicable target.
- A `gap` is valid ONLY after full-text review confirms no compliant evidence exists.
- Any `fix` MUST block section completion and require correction.
- A section completes ONLY when quotations pass and every target is `pass` or validated `gap`.

### Quotation-Rules

- Every quotation MUST contain `tags` and `quote`; `[claim-type]`, `[tags]`, and accurate `[source]` are required.
- Assign exactly one `[claim-type]` from explicit primary claims; satisfy its requirement and declared `quotation`.
- Quote wording MUST remain verbatim from source; standalone formulas are NEVER quotes; `math` MAY follow.
- Each quote MUST use 10–50 words from complete sentences supporting one explicit claim and every tag.
- Overlong evidence MUST join two explicitly-related same-paragraph sentences with `[...]`; NEVER omit meaning-changing text.

### Template-Check

- PRINT: Before every Draft and Re-review, MUST read [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) and print the current section’s `claim-type`, `coverage`, and `quotation` HTML declarations verbatim; NEVER source or infer them from `raw/*.md`.
- Scope: MUST process ONLY one section; review at most four candidate or written quotations per batch.
- Verdict: MUST expand every applicable target; assign exactly one `pass`, `gap`, or `fix`.
- Review: MUST apply [Section-Rules](#section-rules) and [Quotation-Rules](#quotation-rules); mark quotations `pass`/`fix`; record valid `gap`s once; fixes MUST STOP.
- Re-review: After every write, MUST rerun `PRINT`, `Scope`, `Verdict`, and `Review` on the persisted section.

## Draft-Rules

- Prepare: MUST execute [Check-Inputs](#check-inputs); any failure MUST STOP.
- Draft: Follow template order; present at most four candidates, then MUST execute [Template-Check](#template-check).
- Approval: PRINT exact candidate blocks and gaps; STOP until explicit user approval.
- Write: Write ONLY exact approved content; MUST then run `python scripts\sort_raw_quotations.py mdfile_path`.
- Verify: After EVERY write, MUST execute [Template-Check](#template-check); ONLY completed sections advance; incomplete MUST return to [Draft-Rules](#draft-rules).
