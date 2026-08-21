---
name: papers-wiki-rules
description: Apply canonical wiki synthesis rules when creating, drafting, reviewing, or validating wiki files, sections, templates, claims, formulas, references, tags, or evidence coverage in papers-lint or other wiki workflows.
---

# Papers Wiki Rules

## Input-Rules

### Quotation-Inputs

- Raw mode MUST use exactly one specified paper’s `papers/raw/*.md`; NEVER treat `papers/tmp/*.md` as evidence.
- Tag mode MUST use one read-only grouped `papers/tmp/CURRENT_TAG_*.md` set for the user-selected ordered TAG set—current TAG first, then optional supporting TAGs—generated ONLY from `papers/raw/*.md`.
- Quotation-Inputs contain ONLY `Motivation`, `Methods`, `Results`, `Meanings`, and `Secondary Citations`; templates determine usage.
- Tag mode MUST run `python papers\scripts\search_a_tag.py CURRENT_TAG [SUPPORTING_TAG ...]` exactly once before validation; the current TAG MUST be first.
- PRINT: MUST report Wiki template requirements unsupported by selected-mode Quotation-Inputs.

### Vocab-Inputs

- [papers/vocab/tags.json](../../../papers/vocab/tags.json) and [papers/vocab/types.json](../../../papers/vocab/types.json) are auxiliary inputs.
- [papers/vocab/tags.json](../../../papers/vocab/tags.json) defines canonical tags, aliases, definitions, and tag types.
- [papers/vocab/types.json](../../../papers/vocab/types.json) defines claim-type and tag-type requirements.
- Run `python papers\scripts\search_similar_tags.py QUERY 10` for every proposed tag.
- PRINT: Report every unresolved proposed tag before input validation.

### Check-Inputs

- Checkpoint 1: Raw mode MUST run `python papers\scripts\search_a_paper.py DOI`; matching JSON, Markdown, and PDF/TEX MUST exist.
- Checkpoint 2: ONLY one canonical match passes; failure MUST STOP and execute [Vocab Draft-Rules](../papers-vocab-rules/SKILL.md#draft-rules) completely.
- Checkpoint 3: The grouped `papers/tmp/CURRENT_TAG_*.md` set MUST exist and contain matching quotations for every specified TAG; failure MUST STOP.
- Checkpoint 4: Selected-mode Quotation-Inputs MUST cover Wiki needs; failure MUST STOP and execute [Raw Draft-Rules](../papers-raw-rules/SKILL.md#draft-rules) completely.
- PRINT: MUST validate and report [Quotation-Inputs](#quotation-inputs) and [Vocab-Inputs](#vocab-inputs); failures MUST STOP with recommended solutions.

## Template-Rules

- Templates: [papers/scripts/add_wiki_method.md](../../../papers/scripts/add_wiki_method.md) for `method`; [papers/scripts/add_wiki_topic.md](../../../papers/scripts/add_wiki_topic.md) otherwise.
- Generate: ONLY approved current TAG may run `python papers\scripts\add_wiki_md.py TAG`; supporting TAGs NEVER generate Wikis; NEVER scaffold manually.

### Section-Rules

- Each template `###` section MUST contain exactly one `claim-type`/`coverage` pair; NEVER copy declarations.
- MUST preserve template section order, heading levels, required blocks, and evidence-driven repeated `#####` instances.
- Claims and formulas MUST have direct verified Quotation-Inputs evidence, complete scope, and explained formula meaning.
- References MUST use adjacent existing Raw stems; prose MUST synthesize, preserve uncertainty, and separate conflicts.
- A section completes ONLY when items are `supported` or valid `not-applicable`; others MUST block completion.

### Paper-Rules

- Raw mode MUST process exactly one specified Raw paper at a time.
- The specified Raw paper qualifies ONLY when its `Meanings` contains current `TAG`.
- The qualifying specified Raw MUST form exactly one Study using ONLY its quotations.
- After that Study passes, MUST inspect other Wiki sections using its [Quotation-Inputs](#quotation-inputs) for evidence-supported improvements.
- Next Papers MUST run `python papers\scripts\search_a_paper.py DOI` on each Secondary DOI in specified `papers/raw/*.md`, rejecting `recorded:` papers.

### Template-Check

- PRINT: Show current section’s `claim-type` and `coverage` declarations.
- Scope: MUST review ONLY one section—Tag mode: grouped `papers/tmp/CURRENT_TAG_*.md`; Raw mode: one specified `papers/raw/*.md`.
- Verdict: MUST assign exactly one `supported`, `weak`, `missing`, or `not-applicable` to every claim, formula, and target.
- Review: MUST apply [Section-Rules](#section-rules); Raw mode MUST apply [Paper-Rules](#paper-rules) to `Previous Studies`; `weak`/`missing` MUST STOP.
- Re-review: After EVERY write, MUST rerun `PRINT`, `Scope`, `Verdict`, and `Review` on persisted section.

## Draft-Rules

- Prepare: Tag mode MUST use ONLY the user-specified ordered TAG set with current TAG first; selected mode MUST execute [Check-Inputs](#check-inputs); approved missing Tag-mode Wiki MUST generate; failures MUST STOP.
- Draft: Process ONLY one Tag-mode template-order section or one Raw via [Paper-Rules](#paper-rules), starting with `Previous Studies`.
- Approve: MUST execute [Template-Check](#template-check); PRINT at most four exact changes; STOP until explicit user approval.
- Write: ONLY write approved content; `Previous Studies` MUST then run `python papers\scripts\sort_wiki_studies.py wikipath`; NEVER alter unrelated content.
- Verify: After EVERY write, MUST rerun [Template-Check](#template-check); ONLY `supported`/`not-applicable` advances; `weak`/`missing` requires another approved batch.
