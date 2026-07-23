---
name: wiki-rules
description: Apply canonical myWIKI wiki synthesis rules when creating, drafting, reviewing, or validating wiki files, sections, templates, claims, formulas, references, tags, or evidence coverage in lint or other wiki workflows.
---

# Wiki Rules

## Evidence Rules

### Truth and Inputs

- `raw/*.md` is the evidence truth.
- Generated `tmp/*.md` quotation files are read-only Wiki inputs.
- `wiki/*.md` is NEVER original evidence.
- If the user specifies a paper, use ONLY the quotations in that Raw file.
- Otherwise, run `python scripts\search_a_tag.py TAG` and use ONLY its generated quotation files.

### Claim Validity

- A quotation is eligible ONLY when its Raw quotation passes [Raw Section-Rules](../raw-rules/SKILL.md#section-rules).
- Unreviewed, `gap`, or `fix` quotations MUST NOT support Wiki.
- Every factual claim and formula MUST have direct verified primary Raw evidence.
- NEVER infer, reconstruct, generalize, or fill templates from domain knowledge.
- Wiki prose MUST synthesize evidence; NEVER copy quotations as prose.
- Preserve conditions, scope, uncertainty, negation, contrast, and conflicts.
- Every claim and formula MUST have an adjacent `references` block.
- Each reference MUST name an existing Raw file without path or suffix:
  `raw/PAPER.md` -> `PAPER`.
- Each cited quotation MUST support the complete claim, including conditions and scope.
- A related topic, shared keyword, or paper-level relevance is insufficient.
- A formula is valid ONLY when Raw evidence states it or an explicit equivalent
  and explains its meaning.
- Conflicting conclusions MUST remain separate with their conditions and references.
- NEVER collapse conflicting evidence into a false consensus.

### Candidate Papers

- With a specified paper, inspect ONLY its `### Secondary Citations`.
- Otherwise, inspect ONLY the generated Secondary quotation file.
- Deduplicate DOI values within the selected quotation input.
- Run `python scripts\search_a_doi.py DOI` for every unique Secondary DOI.
- Exclude every DOI with a `recorded:` hit.
- Remaining DOI records are Candidate Papers.
- Candidate Papers support discovery only; they NEVER support Wiki claims.
- Report `[doi]: unknown` as an unresolved Secondary reference.
- An unresolved Secondary reference is not a Candidate Paper until its DOI is verified.

## Rules

### Template-Rules

#### Selection and Shape

Wiki Filename: A wiki filename MUST equal its canonical tag resolved under [Vocab Rules](../vocab-rules/SKILL.md).

Template Selection: A `[tag_type]: method` uses [scripts/add_wiki_method.md](../../../scripts/add_wiki_method.md); every other tag type uses [scripts/add_wiki_topic.md](../../../scripts/add_wiki_topic.md).

Template Generation: After the calling workflow approves creation of an absent Wiki, run `python scripts\add_wiki_md.py TAG TAG_TYPE`; NEVER create the Wiki scaffold manually.

Template Integrity: Every `###` section in the selected template MUST contain exactly one `claim-types` and one `coverage` declaration. Declarations are not required for `#####` subsections. Missing or duplicate declarations block drafting and review.

Template Shape: Preserve the selected template's `###` section order, heading levels, and required block shapes; NEVER add, remove, reorder, or reshape its `###` sections. Add or remove repeated `#####` instances ONLY as verified evidence requires.

#### Declaration and Classification

Tags: Resolve every frontmatter and fenced-block tag under [Vocab Rules](../vocab-rules/SKILL.md). Unresolved tags block writing.

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

Gap: If validated Raw gaps prevent completion, ask whether to skip that paper.

Checkpoint: **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval.

Scope: Approval skips only the named paper for the current run.

No Approval: Omission is `missing` and blocks pass.

### Draft-Rules

#### Preparation

- Complete these Draft-Rules before ANY write to `wiki/*.md`.
- NEVER write first and review afterward.
- Apply Truth and Inputs before drafting.
- In specified-paper mode, select and declare exactly one Raw paper.
- Evaluate each Study independently.
- NEVER combine evidence or claims across papers.
- Verify the selected template before drafting or reviewing a section.
- It MUST contain exactly one `claim-types` and one `coverage` declaration.
- Missing or duplicate declarations require STOP.
- Print both declarations.
- Print every relevant claim-type requirement from [vocab/types.json](../../../vocab/types.json).
- Expand `claim-types` first to establish required categories.
- Expand `coverage` second to identify and map targets.
- Use ONLY verified Raw evidence satisfying each category requirement.
- Coverage NEVER creates or changes a claim type.
- In tag-search mode, declare exactly one evidence scope.
- Allowed scopes are `all-primary` or `one-approved-primary-section`.
- For `one-approved-primary-section`, identify the exact generated primary section.
- That scope limits ONLY evidence for new or changed Wiki claims.
- In specified-paper mode, the selected Raw paper is the quotation scope.
- Before approval validation, review the entire current Wiki section with the full active input.
- After ANY write, repeat that full-input review.
- A missing or ambiguous scope requires STOP.
- In tag-search mode, resolve exactly one canonical `TAG` under [Vocab Rules](../vocab-rules/SKILL.md).
- Run `python scripts\search_a_tag.py TAG`.
- A failed command requires STOP.
- Missing `Motivation`, `Methods`, `Results`, or `Meanings` output requires STOP.
- Review every generated quotation file under its Evidence Rules owner.
- In specified-paper mode, review every quotation in the selected Raw paper.
- During Evidence Search, NEVER modify `raw/*.md` or generated files.
- NEVER expand to another tag, sample within scope, or stop early.
- If evidence is insufficient, report every item under `Missing Quotations`.
- Apply [Raw Draft-Rules](../raw-rules/SKILL.md#draft-rules) to those items.
- After an approved Raw change, refresh the active quotation input.
- In tag-search mode, rerun `python scripts\search_a_tag.py TAG`.
- Block approval, writing, and advancement until the quotation appears in the active input.
- Its Raw quotation MUST also receive a `pass` verdict.

#### Review, Approval, and Write

- Draft and discuss exactly one Wiki section per response.
- In tag-search mode, follow selected-template order.
- NEVER mix sections.
- In specified-paper mode, run Gate 2 on `Previous Studies` first.
- This specified-paper step overrides template order.
- Screen the selected paper against every other Wiki section.
- Identify every affected section.
- Run Gate 2 for one affected section per response.
- Follow template order for affected sections.
- Write the smallest synthesis supported by evidence.
- Separate conclusions across systems, assumptions, methods, ranges, or uncertainties.
- Review every quotation in the active input under its Evidence Rules owner.
- Present at most four Wiki candidates or fixes per response.
- Repeat batches until the current Wiki section passes.
- Exact mode provides final block shape, prose, formulas, tags, and references.
- Summary mode provides the candidate claim and supporting Raw location.
- Summary mode also provides one evidence sentence before the checkpoint.
- Before approval, validate every candidate under Evidence Rules and Template-Rules.
- Then apply Section-Rules to the proposed section.
- Any violation, `weak`, or `missing` verdict blocks approval and writing.
- Present the smallest section change set.
- Summary-mode approval authorizes ONLY a later exact-mode proposal.
- Summary-mode approval NEVER authorizes writing.
- Writing requires separate explicit approval of the exact-mode change.
- **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.
- After approval, immediately write ONLY the explicitly approved change.
- Write in exact mode matching the selected template.

#### Verify, Pass, or Fix

- Re-review the ENTIRE current Wiki section under Section-Rules.
- Show `git diff -- wiki`.
- Any failed Re-review returns to a fix batch.
- Advance ONLY after the current Wiki section passes.
- Keep the selected paper active until its Study and every affected section pass.
- **🔴 CHECKPOINT · 🛑 STOP** — Await approval before selecting another paper.
