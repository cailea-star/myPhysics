---
name: wiki-rules
description: Apply canonical myWIKI wiki synthesis rules when creating, drafting, reviewing, or validating wiki files, sections, templates, claims, formulas, references, tags, or evidence coverage in lint or other wiki workflows.
---

# Wiki Rules

## Rules

### Wiki-Rules

`wiki/*.md` is synthesis from verified `raw/*.md` quotation evidence; it is NEVER original evidence.

A wiki filename MUST equal its canonical tag resolved under [Vocab Rules](../vocab-rules/SKILL.md).

Every factual claim and formula MUST be directly supported by verified primary Raw evidence. NEVER infer, generalize beyond evidence, fill a template from domain knowledge, or use `tmp/*_Secondary.md` claims or formulas as evidence.

Wiki prose MUST synthesize evidence; NEVER copy quotation text as wiki prose. Preserve every material condition, scope, uncertainty, negation, contrast, and conflicting conclusion.

Reference: Every factual claim and formula MUST have an adjacent `references` entry naming an existing `raw/*` file without path or suffix: `raw/PAPER.md` -> `PAPER`.

Support: Each cited Raw quotation MUST directly support the complete adjacent claim or formula, including its conditions and scope. A related topic, shared keyword, or paper-level relevance is insufficient.

Formula: A formula is supported ONLY when the cited Raw evidence contains the formula or its explicit equivalent and explains its meaning. NEVER reconstruct a formula from memory or unstated algebra.

Conflict: When verified sources disagree, preserve each conclusion, its conditions, and its references; NEVER collapse conflicts into a false consensus.

Secondary Citations: Treat Secondary Citations only as paper leads. NEVER cite their claims or formulas in Wiki; use an already-ingested original paper's verified direct evidence, propose ingesting it, or report an explicit gap.

Tags: Resolve every frontmatter and fenced-block tag under [Vocab Rules](../vocab-rules/SKILL.md). Unresolved tags block writing.

### Template-Rules

Template Selection: A `[tag_type]: method` uses [scripts/add_wiki_method.md](../../../scripts/add_wiki_method.md); every other tag type uses [scripts/add_wiki_topic.md](../../../scripts/add_wiki_topic.md).

Template Integrity: Every `###` section in the selected template MUST contain exactly one `claim-types` and one `coverage` declaration. Declarations are not required for `#####` subsections. Missing or duplicate declarations block drafting and review.

Template Claim Types: Use only claim types from [vocab/types.json](../../../vocab/types.json). Before drafting or reviewing a wiki section, read its `claim-types` declaration. A required claim type is covered ONLY when sufficient verified Raw evidence has that `[claim_type]` and passes its `requirement`. Judge the quotation's explicit primary claim; keywords, headings, tags, and Coverage targets are not classification evidence. Give every required item exactly one verdict under Section-Rules; NEVER relabel evidence to satisfy the declaration.

Template Coverage: Before drafting or reviewing a wiki section, read its `coverage` declaration from the selected template. Give every required item exactly one verdict under Section-Rules. NEVER invent content or evidence.

Counting: `[and]` separates independently reviewed items; `[or]` joins alternatives within one item. A `for each` item expands to one verdict per identified target. `optional` and `none` contribute zero items, and their absence is NEVER `missing`.

Output: Template declarations are authoritative and MUST NOT be copied into `wiki/*.md`.

### Section-Rules

Verdicts: Give every factual claim, formula, and required template item exactly one verdict:

- `supported`: sufficient verified Raw evidence directly supports it; do nothing.
- `weak`: content exists but its evidence is absent, indirect, insufficient, or mismatched; support, qualify, or remove it.
- `missing`: required content is absent; find verified Raw evidence or record an explicit no-evidence gap.
- `not-applicable`: a required template item does not apply; do not fill it.

Review: Independently review exactly one written wiki section per response in selected-template order. List verdicts in source order and report `Wiki Verdicts: <verdicts>/<factual claims + formulas + required template items>`. A count mismatch fails the section and blocks writing.

Review Pass: A section passes ONLY when counts match, every existing factual claim and formula is `supported`, and every required item is `supported`, `not-applicable`, or an explicit no-evidence gap. Any unresolved `weak` or `missing` blocks advancement.

Re-review: After ANY write to `wiki/*.md`, re-review the ENTIRE current section. Do not advance until it passes.

### Draft-Rules

Before ANY write to `wiki/*.md`, complete these Draft-Rules for the affected section; NEVER write first and review afterward.

Before drafting or reviewing a section, verify that the selected template contains exactly one `claim-types` and one `coverage` declaration for it; otherwise stop. Print both declarations, then print every relevant claim-type requirement from [vocab/types.json](../../../vocab/types.json).

Section: Draft and discuss exactly one wiki section per response in selected-template order; NEVER mix sections. Write the smallest synthesis supported by evidence, separating conclusions when systems, assumptions, methods, parameter ranges, or uncertainties differ. Present at most four candidates or fixes per response, repeating batches until the section passes, in exact mode (final wiki block shape with exact prose, formulas, tags, and references) or summary mode (candidate claim, supporting Raw filename and quotation location, and one-sentence evidence summary before the checkpoint).

Evidence Search: For every `weak` or `missing`, search verified primary Raw evidence first. If none exists, inspect `tmp/TAG_Secondary.md` only for cited-paper leads: use an already-ingested original paper's verified direct evidence; otherwise propose ingesting it or record an explicit gap.

Approval: Present the smallest section change set. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Apply only approved changes, show `git diff -- wiki`, then complete Section-Rules Re-review before the next section.
