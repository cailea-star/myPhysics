---
name: raw-rules
description: Apply canonical myWIKI raw evidence rules when creating, drafting, reviewing, sorting, or validating raw files, quotation sections, claim types, evidence links, templates, or tagged quotations in ingest or lint workflows.
---

# Raw Rules

## Rules

### Raw-Rules

#### Truth and File Roles

`raw` is truth.

`raw/*.json` is metadata.

`raw/*.pdf` and `raw/*.tex` are original full text from network or user.

Source Check: Before drafting or reviewing, run `python -c "from pathlib import Path; assert Path('raw/<basename>.pdf').is_file() or Path('raw/<basename>.tex').is_file()"`.

`raw/*.md` is quotation from original full text with tags.

`tmp/*.md` is a temporary rearrangement of raw quotation evidence generated from `raw/*.md` searches.

#### Global Prohibitions

NEVER edit `tmp/*.md`.

NEVER infer evidence or fill a template without direct source evidence.

### Template-Rules

#### Declaration and Classification

Declarations: Before drafting or review, verify that the corresponding template section contains exactly one `claim-types`, one `coverage`, and one `quotation` declaration.

Declaration Failure: Missing or duplicate declarations require STOP.

Claim Types: Before drafting or review, read the section's `claim-types` declaration from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md).

Vocabulary: Use ONLY claim types from [vocab/types.json](../../../vocab/types.json).

Validity: A quotation satisfies a claim type ONLY when it meets that type's `requirement`.

Definition Ban: Mention, use, effect, or implementation detail is NOT a definition.

Assignment: Assign every quotation exactly one `[claim_type]`.

Basis: Use ONLY the quotation's explicit primary claim.

Classification Ban: Keywords, tags, source section, and coverage targets are NEVER classification evidence.

Relabel Ban: NEVER relabel evidence to satisfy a declaration.

#### Coverage and Quotation

Coverage: Before drafting or review, read the section's `coverage` declaration from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md).

Coverage Verdict: Apply Verdict to every expanded coverage requirement.

Coverage Gap: Missing valid evidence is `gap`.

Coverage Fix: Any exclusion violation is `fix`.

Coverage Ban: NEVER invent evidence.

Quotation: Before drafting or review, read the section's `quotation` declaration from [scripts/add_raw_md.md](../../../scripts/add_raw_md.md).

Quotation Fix: Any quotation-declaration violation is `fix` and blocks approval.

#### Verdict and Expansion

Verdict: Give every declared requirement exactly one `pass`, `gap`, or `fix`.

Pass: Valid evidence MUST match the declared type and requirement.

Gap: No valid evidence is `gap`.

Fix: Any unsupported or incorrect `[claim_type]` is `fix` and blocks approval.

Grouping: Each outer parenthesized group is one independently reviewed item.

And: Plain `and` separates outer parenthesized items.

Or: Plain `or` joins alternatives within one parenthesized item.

Identification: An `identify` item passes ONLY when its target set is complete and correct; any missing, extra, or misclassified target is `fix`.

Expansion: `for each` creates one verdict per identified target.

Optional: `optional` and `none` create no items.

Optional Ban: Their absence is NEVER `gap`.

### Section-Rules

#### Review Contract

Scope: Independently review exactly one written quotation section per response.

Order: Follow [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) section order.

Review: Give every quotation and required template item exactly one `pass`, `gap`, or `fix` verdict.

Report: Print `Review Verdicts: <verdicts>/<quotation blocks + required template items>`.

Count Failure: A verdict-count mismatch fails the section and blocks advancement.

#### Gap-Rules

Eligibility: Before recording or retaining any `gap`, recheck its declared requirement against the full text and all written quotation sections. A `gap` is valid ONLY when the requirement applies, no valid evidence exists, and no necessary Missing Tag remains unresolved.

Record: Every `gap` MUST have exactly one matching `### Gaps` entry.

Form: Use the exact template; all fields required; `[reason]` is 1-10 words.

Identity: `[section]`, `[item]`, and `[target]` MUST identify one verdict.

Sync: Missing, duplicate, extra, malformed, or stale entries block approval.

Ban: NEVER record `pass` or `fix`.

Gap: A `gap` passes ONLY when all Gap-Rules pass.

#### Pass and Advancement

Pass: A section passes ONLY when every quotation and declared requirement passes.

Fix: Any `fix` blocks advancement.

Advance: Continue batches until the section passes; NEVER advance before it passes.

### Quotation-Rules

#### Core Evidence

Source: Every quote MUST come from `raw/*.pdf` or `raw/*.tex`.

Sentence: Every quote MUST contain one or more complete sentences; a phrase alone is invalid.

Claim: Every quote MUST support one clear claim.

Length: Every quote should contain about 10-50 words.

Tags: Every quotation MUST include at least one directly supported canonical tag in `[tags]: ...`.

Additional Tags: Include additional tags ONLY when the retained quotation directly supports each one.

Tag Ban: NEVER invent tag co-occurrence.

Formula Ban: NEVER use a standalone formula as quote text.

Definition Formula: When a formula is important core `[claim_type]: definition` evidence, quote its complete explanatory sentence.

Math Block: Add that formula in a following fenced `math` block.

Format: Formatting is defective ONLY when it violates the template or breaks parsing.

Violation: Any Quotation-Rule violation is `fix` and blocks approval.

#### Evidence-Link-Rules

Purpose: Link method evidence to a `[claim_type]: definition`, `[claim_type]: result`, or `[claim_type]: innovation` claim.

Claim Type: Preserve the core claim's true `[claim_type]`.

Method Coverage: A valid method sentence may satisfy `[tag_type]: method` requirements.

Scope: Both sentences MUST share one paragraph or figure/table discussion.

Scope Failure: Any cross-paragraph or cross-discussion link is invalid.

Relation: The method sentence MUST explicitly state that the `[tag_type]: method` is used to obtain, define, or realize the core claim.

Relation Failure: A method sentence without that direct relation is invalid.

Omission: Use `[...]` only for nonessential intervening text.

Integrity: NEVER omit conditions, negations, contrasts, or any other meaning-changing text.

Assignment: Assign `[claim_type]` from the core claim.

Method Tag: A retained `[tag_type]: method` tag MUST be directly supported by the retained method sentence.

Violation: Any invalid evidence link is `fix` and blocks approval.

#### Secondary-Citations-Check

Role: Secondary Citations are paper leads only.

Gap Lead: Prioritize an explicitly linked external citation when it supplies a core method definition, current-output quantity definition, or core mechanism missing under a validated gap; it remains a paper lead and NEVER satisfies that gap.

Evidence Ban: NEVER use their quoted claims or formulas as direct evidence.

Raw Ban: NEVER use them to backfill `raw/*.md`.

Wiki Ban: NEVER use them to update wiki.

Trigger: Any quotation containing an explicit external citation marker requires this check.

Classification: A triggered quotation MUST be classified as citation-marked for this check.

Examples: Citation markers include `[55]` and `Ref. [45]` and `\cite{...}`.

External Core: If a citation-marked quotation's core information comes from the cited work, it MUST be placed under `### Secondary Citations`.

Current Core: A citation-marked current-paper problem setting, method use, calculated result, or conclusion judgment MUST remain eligible for its main section.

Eligible Sections: Main sections are `Motivation`, `Methods`, `Results`, and `Meanings`.

Violation: Misplacing citation-marked evidence is `fix` and blocks approval.

### Draft-Rules

#### Preparation

Preflight: Complete all pre-write Draft-Rules for the current batch before ANY write.

Order: Follow [scripts/add_raw_md.md](../../../scripts/add_raw_md.md) section order; NEVER write before drafting and review are complete.

Section: Draft, present, and discuss exactly one quotation section per response.

Isolation: NEVER mix quotation sections.

Print Declarations: Print all three declarations.

Print Requirements: Print every relevant claim-type requirement.

Print Tag Types: Print every requirement for the proposed tags' tag types from [vocab/types.json](../../../vocab/types.json).

Introduction Map: For `definition` or `background` coverage, map every covered tag to its proposed claim type and exact qualifying clause.

Resolve: Before presenting or writing a candidate, resolve every proposed `[tags]` entry under [Vocab Rules](../vocab-rules/SKILL.md).

Missing: If the quotation directly supports a necessary unresolved concept, report it under `Missing Tags`.

Draft: Apply [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules) to every Missing Tag.

Block: Any unresolved Missing Tag blocks approval and writing.

Scope Ban: NEVER create tags for unrelated background mentions.

#### Batch Review and Approval

Batch: Present at most four candidates or fixes per response.

Exact Mode: In [exact mode], use the final `raw/*.md` block shape with exact quoted sentence(s), any required `math` block, `[claim_type]`, `[tags]`, and source section.

Summary Mode: In [summary mode], include source section, candidate `[claim_type]`, candidate `[tags]`, and one-sentence evidence summary.

Review: Before the checkpoint, apply Template-Rules to the current section; apply Quotation-Rules, any triggered Evidence-Link-Rules, and any triggered Secondary-Citations-Check to every candidate.

Pre-write Violation: Any pre-write Draft-Rule violation blocks approval and writing.

Approval: Present either mode, then **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; NEVER proceed without it.

#### Write and Re-review

Mode: Write approved quotations to `raw/*.md` in [exact mode] matching [scripts/add_raw_md.md](../../../scripts/add_raw_md.md).

Write: Write ONLY explicitly approved quotations.

Sort: After ANY write to `raw/*.md`, run `python scripts\sort_raw_md_quotations.py mdfile_path`.

Re-review: After sorting, re-review the ENTIRE current section under Section-Rules.

Post-write Violation: Any Re-review failure blocks advancement; any Re-review `fix` requires another fix batch.
