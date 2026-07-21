---
name: lint
description: Audit and update myWIKI from tagged quotations. Use when checking wiki coverage, backfilling missing evidence into raw, or reconciling wiki conclusions with vocab tags.
---

# Lint

## Comparison Rules

Before Gate 1, read [Vocab Rules](../vocab-rules/SKILL.md), [Raw Rules](../raw-rules/SKILL.md), and [Wiki Rules](../wiki-rules/SKILL.md) completely and apply them throughout every gate.

### Quotation-to-Wiki Rules

Give every generated primary quotation exactly one verdict:

- `covered`: Wiki states the conclusion and cites the paper.
- `miss-reference`: Wiki states the conclusion but does not cite the paper.
- `needs-inclusion`: Wiki omits an important supported conclusion.
- `conflict`: The quotation conflicts with Wiki.
- `irrelevant`: The quotation does not belong in this tag's Wiki.

Give every non-`irrelevant` quotation exactly one `[wiki-section]` destination using the selected Wiki template's exact `###` heading. Resolve it from the quotation's explicit primary claim against each Wiki section's `claim-types` and `coverage`; the generated primary section name is NEVER destination evidence. Use `[wiki-section]: none` only for `irrelevant`. Zero or multiple valid destinations is `routing-gap` and blocks advancement.

List one `{verdict, wiki-section}` pair per `quote` block in source order and report `Verdicts: <verdicts>/<quote blocks>`. A count mismatch fails the generated section and blocks advancement. Verdicts may share one Wiki change but MUST NOT be merged.

Queue every `miss-reference`, `needs-inclusion`, or `conflict` under its resolved `[wiki-section]`; `covered` and `irrelevant` require no Wiki change. Apply [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) only after all four generated primary sections pass verdict review and routing.

### Wiki-to-Quotation Rules

Apply [Wiki Section-Rules](../wiki-rules/SKILL.md#section-rules) to identify `weak` and `missing` items. Search all generated primary sections first.

If generated primary evidence is insufficient, inspect relevant already-ingested primary full text under [Raw Rules](../raw-rules/SKILL.md) and draft any valid direct evidence under [Raw Draft-Rules](../raw-rules/SKILL.md#draft-rules). Only after finding no valid direct evidence may you inspect `tmp/TAG_Secondary.md` for cited-paper leads. For an already-ingested original paper, apply Raw Draft-Rules; otherwise propose ingesting it. If no valid evidence or lead exists, report the evidence gap and STOP; a gap does not pass the Wiki section.

After ANY approved Raw change, complete [Raw Re-review](../raw-rules/SKILL.md#section-rules), then run:

```powershell
python scripts\search_a_tag.py TAG
```

The new quotation MUST appear in generated primary output and receive a Quotation-to-Wiki verdict before supporting Wiki; otherwise the section fails. Apply [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) only after this check passes.

## Gated Workflow

Run gates strictly in order. Start every response with the current gate, last completed gate, and next approval. NEVER advance from an unapproved or failed section or gate.

### Gate 1 — Confirm Tag & Search

1. Run; if either check fails or the worktree is dirty, **🛑 STOP** and report:
   ```powershell
   git rev-parse --is-inside-work-tree
   git status --short
   ```
2. Resolve exactly one canonical tag and its types under [Vocab Rules](../vocab-rules/SKILL.md), then confirm its Wiki path and required template under [Wiki Rules](../wiki-rules/SKILL.md); if no valid existing tag resolves, **🛑 STOP**.
3. If Wiki is absent, **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval to draft it from the required template; do not write, search, audit, or propose content changes before approval. After approval, run `python scripts\add_wiki_md.py TAG TAG_TYPE`; if it fails, **🛑 STOP**.
4. Complete [Wiki Evidence Search](../wiki-rules/SKILL.md#draft-rules); do not advance unless every expected primary output exists.

### Gate 2 — Quotations to Wiki

1. Apply Quotation-to-Wiki Rules, including destination routing, to exactly one generated primary section per response, in order: Motivation, Methods, Results, Meanings.
2. Report its complete verdict set. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next generated section; do not proceed.
3. After all four generated primary sections pass, process destination Wiki sections in template order. For exactly one Wiki section per response, collect every assigned actionable verdict and present the smallest exact-mode Wiki draft under [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules). **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not write or proceed.
4. After explicit approval, write that draft immediately, re-review the entire Wiki section, and show `git diff -- wiki`. A failed re-review returns to fixes; a pass advances to the next Wiki section. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next Wiki section; do not proceed.

### Gate 3 — Wiki to Quotations to Wiki

1. Apply [Wiki Section-Rules](../wiki-rules/SKILL.md#section-rules) to exactly one Wiki section per response in template order.
2. If the Wiki section passes, report pass. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next Wiki section; do not proceed.
3. Otherwise resolve every `weak` or `missing` under Wiki-to-Quotation Rules. Complete any required Raw Draft, Raw Re-review, regeneration, and Quotation-to-Wiki verdict before proposing a Wiki change.
4. Apply [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) to the resulting Wiki candidates or fixes. Do not advance until the current Wiki section passes.
5. Show `git diff -- raw wiki`. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next Wiki section; do not proceed.

### Gate 4 — Wiki to Tags

1. Begin only after every Wiki section passes and no actionable Quotation-to-Wiki verdict remains.
2. Reconcile vocab under [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules), requiring zero unsupported Wiki claims.
3. Present exactly one smallest vocab change. **🔴 CHECKPOINT · 🛑 STOP** — Await separate explicit approval; do not proceed.
4. Apply only the approved change, complete Vocab verification, and show `git diff -- raw wiki vocab`.
