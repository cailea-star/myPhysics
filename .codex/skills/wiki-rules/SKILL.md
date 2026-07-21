---
name: wiki-rules
description: Apply canonical myWIKI wiki synthesis rules when creating, drafting, reviewing, or validating wiki files, sections, templates, claims, formulas, references, tags, or evidence coverage in lint or other wiki workflows.
---

# Wiki Rules

## Rules

### Wiki-Rules

#### Truth and File Roles

`raw/*.md` is the underlying evidence truth; generated primary `tmp/*.md` sections are read-only direct input for `wiki/*.md` synthesis. Wiki is NEVER original evidence.

A generated quotation may support Wiki ONLY when its source Raw quotation and entire containing section pass [Raw Section-Rules](../raw-rules/SKILL.md#section-rules); unreviewed, `gap`, or `fix` evidence MUST NOT support Wiki.

A wiki filename MUST equal its canonical tag resolved under [Vocab Rules](../vocab-rules/SKILL.md).

#### Evidence and Synthesis

Every factual claim and formula MUST be directly supported by verified primary Raw evidence. NEVER infer, generalize beyond evidence, fill a template from domain knowledge, or use `tmp/*_Secondary.md` claims or formulas as evidence.

Wiki prose MUST synthesize evidence; NEVER copy quotation text as wiki prose. Preserve every material condition, scope, uncertainty, negation, contrast, and conflicting conclusion.

Reference: Every factual claim and formula MUST have an adjacent `references` entry naming an existing `raw/*` file without path or suffix: `raw/PAPER.md` -> `PAPER`.

Support: Each cited Raw quotation MUST directly support the complete adjacent claim or formula, including its conditions and scope. A related topic, shared keyword, or paper-level relevance is insufficient.

Formula: A formula is supported ONLY when the cited Raw evidence contains the formula or its explicit equivalent and explains its meaning. NEVER reconstruct a formula from memory or unstated algebra.

Conflict: When verified sources disagree, preserve each conclusion, its conditions, and its references; NEVER collapse conflicts into a false consensus.

Tags: Resolve every frontmatter and fenced-block tag under [Vocab Rules](../vocab-rules/SKILL.md). Unresolved tags block writing.

### Template-Rules

#### Selection and Shape

Template Selection: A `[tag_type]: method` uses [scripts/add_wiki_method.md](../../../scripts/add_wiki_method.md); every other tag type uses [scripts/add_wiki_topic.md](../../../scripts/add_wiki_topic.md).

Template Integrity: Every `###` section in the selected template MUST contain exactly one `claim-types` and one `coverage` declaration. Declarations are not required for `#####` subsections. Missing or duplicate declarations block drafting and review.

Template Shape: Preserve the selected template's `###` section order, heading levels, and required block shapes; NEVER add, remove, reorder, or reshape its `###` sections. Add or remove repeated `#####` instances ONLY as verified evidence requires.

#### Declaration and Classification

Template Claim Types: Use only claim types from [vocab/types.json](../../../vocab/types.json). Before drafting or reviewing a wiki section, read its `claim-types` declaration. A required claim type is covered ONLY when sufficient verified Raw evidence has that `[claim_type]` and passes its `requirement`. Judge the quotation's explicit primary claim; keywords, headings, tags, and Coverage targets are not classification evidence. Give every required item exactly one verdict under Section-Rules; NEVER relabel evidence to satisfy the declaration.

Template Coverage: Before drafting or reviewing a wiki section, read its `coverage` declaration from the selected template. Give every required item exactly one verdict under Section-Rules. NEVER invent content or evidence.

#### Verdict and Expansion

Counting: `[and]` separates independently reviewed items; `[or]` joins alternatives within one item. A `for each` item expands to one verdict per identified target. `optional` and `none` contribute zero items, and their absence is NEVER `missing`.

Output: Template declarations are authoritative and MUST NOT be copied into `wiki/*.md`.

### Section-Rules

#### Review Contract

Review: Independently review exactly one current or proposed wiki section per response in selected-template order. List verdicts in source order and report `Wiki Verdicts: <verdicts>/<factual claims + formulas + required template items>`. A count mismatch fails the section and blocks writing.

#### Verdict Contract

Verdicts: Give every factual claim, formula, and required template item exactly one verdict:

- `supported`: sufficient verified Raw evidence directly supports it; do nothing.
- `weak`: content exists but its evidence is absent, indirect, insufficient, or mismatched; report it to the calling workflow.
- `missing`: required content is absent; report it to the calling workflow.
- `not-applicable`: a required template item does not apply; do not fill it.

#### Pass and Advancement

Review Pass: A section passes ONLY when counts match, every existing factual claim and formula is `supported`, and every required item is `supported` or `not-applicable`. Any `weak` or `missing` unresolved by the proposed change blocks approval, writing, and advancement.

Re-review: After ANY write to `wiki/*.md`, re-review the ENTIRE current section. Do not advance until it passes.

### Draft-Rules

#### Preparation

Before ANY write to `wiki/*.md`, complete these Draft-Rules for the affected section; NEVER write first and review afterward.

Before drafting or reviewing a section, verify that the selected template contains exactly one `claim-types` and one `coverage` declaration for it; otherwise stop. Print both declarations, then print every relevant claim-type requirement from [vocab/types.json](../../../vocab/types.json).

#### Batch Review and Approval

Section: Draft and discuss exactly one wiki section per response in selected-template order; NEVER mix sections. Write the smallest synthesis supported by evidence, separating conclusions when systems, assumptions, methods, parameter ranges, or uncertainties differ. Present at most four candidates or fixes per response, repeating batches until the section passes, in exact mode (final wiki block shape with exact prose, formulas, tags, and references) or summary mode (candidate claim, supporting Raw filename and quotation location, and one-sentence evidence summary before the checkpoint).

Evidence Search: Resolve exactly one canonical `TAG` under [Vocab Rules](../vocab-rules/SKILL.md) and run `python scripts\search_a_tag.py TAG`; if it fails or any expected primary output (`Motivation`, `Methods`, `Results`, or `Meanings`) is missing, STOP. Review every generated primary section as the read-only primary quotation reference before drafting or review. NEVER modify `raw/*.md` or generated sections, use `tmp/TAG_Secondary.md` as evidence, expand to another tag, sample sections, or stop early. If the evidence remains insufficient, report `weak` or `missing` to the calling workflow and block approval, writing, and advancement.

Validation: Before approval, validate every candidate under Wiki-Rules and Template-Rules, then apply Section-Rules to the resulting proposed section. Any violation or unresolved `weak` or `missing` blocks approval and writing.

Approval: Present the smallest section change set. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.

#### Write and Re-review

Apply only approved changes, show `git diff -- wiki`, then complete Section-Rules Re-review before the next section.
