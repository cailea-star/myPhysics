---
name: wiki-rules
description: Apply canonical myWIKI wiki synthesis rules when creating, drafting, reviewing, or validating wiki files, sections, templates, claims, formulas, references, tags, or evidence coverage in lint or other wiki workflows.
---

# Wiki Rules

## Rules

### Wiki-Rules

#### Truth and File Roles

`raw/*.md` is the underlying evidence truth; generated primary `tmp/*.md` sections are read-only direct input for `wiki/*.md` synthesis. Wiki is NEVER original evidence.

A generated quotation may support Wiki ONLY when its source Raw quotation has a `pass` verdict under [Raw Section-Rules](../raw-rules/SKILL.md#section-rules); an unreviewed source quotation or one with a `gap` or `fix` verdict MUST NOT support Wiki.

A wiki filename MUST equal its canonical tag resolved under [Vocab Rules](../vocab-rules/SKILL.md).

#### Evidence and Synthesis

Every factual claim and formula MUST be directly supported by verified primary Raw evidence. NEVER infer, generalize beyond evidence, fill a template from domain knowledge, or use `tmp/*_Secondary.md` claims or formulas as evidence.

Wiki prose MUST synthesize evidence; NEVER copy quotation text as wiki prose. Preserve every material condition, scope, uncertainty, negation, contrast, and conflicting conclusion.

Reference: Every factual claim and formula MUST have an adjacent `references` entry naming an existing `raw/*` file without path or suffix: `raw/PAPER.md` -> `PAPER`.

Support: Each cited Raw quotation MUST directly support the complete adjacent claim or formula, including its conditions and scope. A related topic, shared keyword, or paper-level relevance is insufficient.

Formula: A formula is supported ONLY when the cited Raw evidence contains the formula or its explicit equivalent and explains its meaning. NEVER reconstruct a formula from memory or unstated algebra.

Conflict: When verified sources disagree, preserve each conclusion, its conditions, and its references; NEVER collapse conflicts into a false consensus.

Tags: Resolve every frontmatter and fenced-block tag under [Vocab Rules](../vocab-rules/SKILL.md). Unresolved tags block writing.

Next Papers: Run `python scripts\search_a_doi.py DOI` for every Secondary DOI; exclude any paper with a `recorded:` hit.

### Template-Rules

#### Selection and Shape

Template Selection: A `[tag_type]: method` uses [scripts/add_wiki_method.md](../../../scripts/add_wiki_method.md); every other tag type uses [scripts/add_wiki_topic.md](../../../scripts/add_wiki_topic.md).

Template Generation: After the calling workflow approves creation of an absent Wiki, run `python scripts\add_wiki_md.py TAG TAG_TYPE`; NEVER create the Wiki scaffold manually.

Template Integrity: Every `###` section in the selected template MUST contain exactly one `claim-types` and one `coverage` declaration. Declarations are not required for `#####` subsections. Missing or duplicate declarations block drafting and review.

Template Shape: Preserve the selected template's `###` section order, heading levels, and required block shapes; NEVER add, remove, reorder, or reshape its `###` sections. Add or remove repeated `#####` instances ONLY as verified evidence requires.

#### Declaration and Classification

Template Claim Types: Use only claim types from [vocab/types.json](../../../vocab/types.json). Before drafting or reviewing a wiki section, read its `claim-types` declaration. A required claim type is covered ONLY when sufficient verified Raw evidence has that `[claim_type]` and passes its `requirement`. Judge the quotation's explicit primary claim; keywords, headings, tags, and Coverage targets are not classification evidence. Give every required item exactly one verdict under Section-Rules; NEVER relabel evidence to satisfy the declaration.

Template Coverage: Before drafting or reviewing a wiki section, read its `coverage` declaration from the selected template. Give every required item exactly one verdict under Section-Rules. NEVER invent content or evidence.

#### Verdict and Expansion

Grouping: Each outer parenthesized group is one independently reviewed item.

And: Plain `and` separates outer parenthesized items.

Or: Plain `or` joins alternatives within one parenthesized item.

Identification: An `identify` item is `supported` ONLY when its target set is complete and correct; if it contains any extra or misclassified target, its verdict is `weak`; otherwise, if any target is missing, its verdict is `missing`.

Expansion: `for each` creates one verdict per identified target.

Equation Grouping:

- Group only parallel relations in one `math` block serving one definition or claim.
- Count derived, sequential, or independently argued relations as separate core equations.
- Review every displayed relation separately as a factual formula.

Optional: `optional` and `none` contribute zero items.

Optional Ban: Their absence is NEVER `missing`.

Empty Target Placeholder: For a complete empty target set without an explicit `not-applicable` group, preserve its fixed template placeholder verbatim; it is structural, not factual, and contributes zero review items. NEVER use it for an incomplete search or nonempty target.

Not-Applicable: An explicit `(if no TARGET qualifies, report not-applicable)` group produces one `not-applicable` verdict when the corresponding identified target set is complete and empty; when the set is nonempty, the group produces zero items.

Not-Applicable Ban: A required claim or coverage assignment for an identified target is NEVER `not-applicable`; absent valid content is `missing`.

Output: Template declarations are authoritative and MUST NOT be copied into `wiki/*.md`.

### Section-Rules

#### Review Contract

Review: Independently review exactly one current or proposed wiki section per response in selected-template order. List verdicts in source order and report `Wiki Verdicts: <verdicts>/<factual claims + formulas + required template items>`. A count mismatch fails the section and blocks writing.

#### Verdict and Pass

Verdicts: Give every factual claim, formula, and required template item exactly one verdict:

- `supported`: sufficient verified Raw evidence directly supports it; do nothing.
- `weak`: content exists but its evidence is absent, indirect, insufficient, or mismatched; report it to the calling workflow.
- `missing`: required content is absent; report it to the calling workflow.
- `not-applicable`: an explicitly authorized identified target set is complete and empty; do not fill it.

Review Pass: A section passes ONLY when counts match, every existing factual claim and formula is `supported`, and every required item is `supported` or `not-applicable`. Any `weak` or `missing` unresolved by the proposed change blocks approval, writing, and advancement.

Re-review: After ANY write to `wiki/*.md`, re-review the ENTIRE current section. Do not advance until it passes.

#### Previous Studies Rules

Target: Every unique Raw paper in generated Meanings MUST have exactly one Study.

Unit: Process exactly one paper at a time; after its required whole-section re-review, **🛑 STOP** before selecting another.

Gap: If validated Raw gaps prevent completion, ask whether to skip that paper.

Checkpoint: **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval.

Scope: Approval skips only the named paper for the current run.

No Approval: Omission is `missing` and blocks pass.

Isolation: Evaluate each Study independently. NEVER combine claims across papers.

### Draft-Rules

#### Preparation

Before ANY write to `wiki/*.md`, complete these Draft-Rules for the affected section; NEVER write first and review afterward.

Before drafting or reviewing a section, verify that the selected template contains exactly one `claim-types` and one `coverage` declaration for it; otherwise stop. Print both declarations, then print every relevant claim-type requirement from [vocab/types.json](../../../vocab/types.json).

Claim-Type Lead: First expand `claim-types` to establish the required claim categories. Then expand `coverage` to identify targets and map each target to those declared categories using ONLY verified Raw evidence that already satisfies the category requirement; coverage NEVER creates or changes a claim type.

Evidence Scope: Before drafting, the calling workflow MUST declare exactly one evidence scope: `all-primary` or `one-approved-primary-section`. For `one-approved-primary-section`, identify the exact generated primary section; this scope limits ONLY the evidence used to draft new or changed Wiki claims. Before approval validation and after ANY write, use `all-primary` to review the entire current Wiki section. A missing or ambiguous scope requires STOP.

Evidence Search: Resolve exactly one canonical `TAG` under [Vocab Rules](../vocab-rules/SKILL.md), then run `python scripts\search_a_tag.py TAG` to regenerate the primary quotation outputs. If it fails or any expected primary output (`Motivation`, `Methods`, `Results`, or `Meanings`) is missing, STOP. Review every generated primary section in the declared evidence scope as read-only quotation evidence. During Evidence Search, NEVER modify `raw/*.md` or generated sections, use `tmp/TAG_Secondary.md` as evidence, expand to another tag, sample within the declared scope, or stop early. If evidence remains insufficient, report each item under `Missing Quotations` and apply [Raw Draft-Rules](../raw-rules/SKILL.md#draft-rules). After any approved Raw change, regenerate the primary outputs; block Wiki approval, writing, and advancement until the quotation appears there and its source Raw quotation receives a `pass` verdict.

#### Batch Review and Approval

Section: Draft and discuss exactly one wiki section per response in selected-template order; NEVER mix sections. Write the smallest synthesis supported by evidence, separating conclusions when systems, assumptions, methods, parameter ranges, or uncertainties differ. Present at most four Wiki candidates or fixes per response after reviewing every generated primary section in the declared evidence scope, repeating batches until the current wiki section passes, in exact mode (final wiki block shape with exact prose, formulas, tags, and references) or summary mode (candidate claim, supporting Raw filename and quotation location, and one-sentence evidence summary before the checkpoint).

Validation: Before approval, validate every candidate under Wiki-Rules and Template-Rules, then apply Section-Rules to the resulting proposed section. Any violation or unresolved `weak` or `missing` blocks approval and writing.

Approval: Present the smallest section change set. Summary-mode approval authorizes ONLY a subsequent exact-mode proposal and NEVER authorizes writing. Write requires separate explicit approval of the final exact-mode change. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.

#### Write and Re-review

Write: After approval, immediately write ONLY the explicitly approved change to the current Wiki section in exact mode matching the selected template.

Re-review: Re-review the ENTIRE current Wiki section under Section-Rules, then show `git diff -- wiki`.

Post-write Violation: Any failed Re-review returns to a fix batch. Advance ONLY after the current Wiki section passes.
