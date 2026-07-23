---
name: raw-rules
description: Apply canonical myWIKI raw evidence rules when creating, drafting, reviewing, sorting, or validating raw files, quotation sections, claim types, evidence links, templates, or tagged quotations in ingest or lint workflows.
---

# Raw Rules

## Input-Rules

### Raw-Inputs

- `raw` is truth; `raw/*.json` is metadata; matching `raw/*.pdf` or `raw/*.tex` is source.
- Before every draft, run `python scripts\search_a_doi.py DOI`.
- Exactly one `raw/*.json` MUST match the target DOI.
- At least one PDF or TEX MUST share the matched JSON basename.
- PRINT: Missing, duplicate, or mismatched Raw inputs require STOP.

### Vocab-Inputs

- [vocab/tags.json](../../../vocab/tags.json) and [vocab/types.json](../../../vocab/types.json) are auxiliary inputs.
- [vocab/tags.json](../../../vocab/tags.json) defines canonical tags, aliases, definitions, and tag types.
- [vocab/types.json](../../../vocab/types.json) defines claim-type and tag-type requirements.
- Run `python scripts\search_similar_tags.py QUERY 3` for every proposed tag.
- PRINT: Show required type requirements and canonical-tag results before approval.

### Check-Inputs

- Checkpoint 1: Source passes when exactly one DOI-matching JSON has same-basename PDF or TEX.
- Checkpoint 2: Vocabulary passes when every proposed tag resolves to one canonical tag.
- Checkpoint 3: Target mode is create if Markdown is absent; continue if basename and `[source_json]` match.
- PRINT: Any failure requires STOP; report solutions; NEVER overwrite mismatched Markdown.


## Template-Rules

- Template: [scripts/add_raw_md.md](../../../scripts/add_raw_md.md); manual: `Motivation`, `Methods`, `Results`, `Meanings`, `Secondary Citations`; auxiliary: `Abstract`, `Gaps`.

### Section-Rules

- Each quotation section MUST contain exactly one `claim-types`, `coverage`, and `quotation` HTML declaration.
- Section quotations collectively assign exactly one `pass`, `gap`, or `fix` to every applicable target.
- Missing valid full-text evidence MUST create exactly one matching `Gaps` entry.
- Any `fix` MUST block section completion and require correction.
- A section completes ONLY when quotations pass and every target is `pass` or validated `gap`.

### Quotation-Rules

- Every quotation MUST contain `tags` and `quote`; `[claim_type]`, `[tags]`, and accurate `[source]` are required.
- Assign exactly one `[claim_type]` from explicit primary claims; satisfy its requirement and declared `quotation`.
- Quote wording MUST remain verbatim from source; standalone formulas are NEVER quotes; `math` MAY follow.
- Each quote MUST use 10–50 words from complete sentences supporting one explicit claim and every tag.
- Overlong evidence MUST join two explicitly-related same-paragraph sentences with `[...]`; NEVER omit meaning-changing text.

### Template-Check

- PRINT: Show current section’s `claim-types`, `quotation`, and `coverage` declarations.
- Scope: Process ONLY one section; review at most four candidate or written quotations per batch.
- Verdict: Expand every applicable target; assign exactly one `pass`, `gap`, or `fix`.
- Review: Mark quotations `pass` or `fix`; record each valid `gap` once; any `fix` requires STOP.
- Re-review: After every write, rerun `PRINT`, `Scope`, `Verdict`, and `Review` on the persisted section.

## Draft-Rules

- Prepare: Apply all `Input-Rules`; any failure requires STOP.
- Draft: Follow template order; present at most four candidates for current section per batch.
- Approval: PRINT candidate quotations and gap entries, then STOP until explicit user approval.
- Write: Write ONLY approved quotations and gaps; run `python scripts\sort_raw_md_quotations.py mdfile_path`.
- Verify: Run `Template-Check` after every write; incomplete sections return to Draft, completed sections advance.
