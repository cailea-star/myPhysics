---
name: wiki-rules
description: Apply canonical myWIKI wiki synthesis rules when creating, drafting, reviewing, or validating wiki files, sections, templates, claims, formulas, references, tags, or evidence coverage in lint or other wiki workflows.
---

# Wiki Rules

## Input-Rules

### Quotation-Inputs

- `raw/*.md` contains quotations from exactly one specified paper.
- `tmp/*.md` contains read-only `TAG` quotation summaries generated ONLY from `raw/*.md`.
- Quotation-Inputs contain ONLY `Motivation`, `Methods`, `Results`, `Meanings`, and `Secondary Citations`; templates determine usage.
- Using `tmp/*.md` MUST run `python scripts\search_a_tag.py TAG` before input validation.
- PRINT: Report Wiki template requirements unsupported by current Quotation-Inputs.

### Vocab-Inputs

- [vocab/tags.json](../../../vocab/tags.json) and [vocab/types.json](../../../vocab/types.json) are auxiliary inputs.
- [vocab/tags.json](../../../vocab/tags.json) defines canonical tags, aliases, definitions, and tag types.
- [vocab/types.json](../../../vocab/types.json) defines claim-type and tag-type requirements.
- Run `python scripts\search_similar_tags.py QUERY 3` for every proposed tag.
- PRINT: Report every unresolved proposed tag before input validation.

### Check-Inputs

- Checkpoint 1: Specified Raw MUST run `python scripts\search_a_paper.py DOI`; matching JSON, Markdown, and PDF-or-TEX MUST exist.
- Checkpoint 2: One canonical match passes; failure MUST STOP and follow [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules).
- Checkpoint 3: Canonical `TAG` MUST run `python scripts\search_a_tag.py TAG`; generated `tmp/*.md` MUST contain matching quotations.
- Checkpoint 4: Raw quotations MUST cover Wiki needs; failure MUST STOP and follow [Raw Draft-Rules](../raw-rules/SKILL.md#draft-rules).
- PRINT: Any failure requires STOP; report the problem and recommended solution.

## Template-Rules

- Templates: [scripts/add_wiki_method.md](../../../scripts/add_wiki_method.md) for `method`; [scripts/add_wiki_topic.md](../../../scripts/add_wiki_topic.md) otherwise.
- Generate: ONLY after explicit creation approval, run `python scripts\add_wiki_md.py TAG TAG_TYPE`; NEVER create scaffolds manually.

### Section-Rules

- Each template `###` section MUST contain exactly one `claim-types`/`coverage` pair; NEVER copy declarations.
- MUST preserve template section order, heading levels, required blocks, and evidence-driven repeated `#####` instances.
- Claims and formulas MUST have direct verified Quotation-Inputs evidence, complete scope, and explained formula meaning.
- References MUST use adjacent existing Raw stems; prose MUST synthesize, preserve uncertainty, and separate conflicts.
- A section completes ONLY when items are `supported` or valid `not-applicable`; others MUST block completion.

### Paper-Rules

- MUST process exactly one specified Raw at a time.
- The specified Raw paper qualifies ONLY when its `Meanings` contains current `TAG`.
- The qualifying specified Raw MUST form exactly one Study using ONLY its quotations.
- After that Study passes, use its Quotation-Inputs to inspect other Wiki sections for evidence-supported improvements.
- Next Papers MUST run `python scripts\search_a_paper.py DOI` on each Secondary DOI in specified `raw/*.md`, rejecting `recorded:` papers.

### Template-Check

- PRINT: Show current section’s `claim-types` and `coverage` declarations.
- Scope: Review ONLY one entire section using all active Quotation-Inputs.
- Verdict: Assign exactly one `supported`, `weak`, `missing`, or `not-applicable` to every claim, formula, and target.
- Review: Apply [Section-Rules](#section-rules); `Previous Studies` MUST additionally pass [Paper-Rules](#paper-rules); `weak` or `missing` requires STOP.
- Sort: After every `Previous Studies` write, MUST run `python scripts\sort_wiki_studies.py wikipath` before Re-review.
- Re-review: After every write, rerun `PRINT`, `Scope`, `Verdict`, and `Review` on persisted section.

## Draft-Rules

- Prepare: Apply [Input-Rules](#input-rules); any failure MUST STOP before drafting.
- Draft: Process ONLY one section; MUST follow template order unless [Paper-Rules](#paper-rules) requires otherwise.
- MUST process the specified Raw through `Previous Studies` first; otherwise iterate `tmp/*_Meanings.md` papers one-by-one identically.
- Approve: PRINT at most four exact changes; STOP until explicit user approval.
- Write: Write ONLY exact approved content; NEVER alter unrelated Wiki content.
- Verify: Rerun [Template-Check](#template-check); `pass` advances, while `fix` requires another approved batch.
