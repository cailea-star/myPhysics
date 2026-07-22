---
name: lint
description: Audit and update myWIKI from tagged quotations. Use when checking wiki coverage, backfilling missing evidence into raw, or reconciling wiki conclusions with vocab tags.
---

# Lint

Before Gate 1, read [Vocab Rules](../vocab-rules/SKILL.md), [Raw Rules](../raw-rules/SKILL.md), and [Wiki Rules](../wiki-rules/SKILL.md) completely and apply them throughout every gate.

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
4. Declare Wiki evidence scope `all-primary`, then complete [Wiki Evidence Search](../wiki-rules/SKILL.md#draft-rules); do not advance unless every expected primary output exists.

### Gate 2 — Wiki to Quotations to Wiki

1. Declare Wiki evidence scope `all-primary`, then apply [Wiki Section-Rules](../wiki-rules/SKILL.md#section-rules) to exactly one Wiki section per response in template order.
2. If the Wiki section passes, report pass. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next Wiki section; do not proceed.
3. Otherwise resolve every `weak` or `missing` through [Wiki Evidence Search](../wiki-rules/SKILL.md#draft-rules).
   - Apply [Raw Draft-Rules](../raw-rules/SKILL.md#draft-rules) when new evidence is required.
   - Regenerate before Wiki drafting.
4. Apply [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) to the resulting Wiki candidates or fixes. Do not advance until the current Wiki section passes.
5. Show `git diff -- raw wiki`. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next Wiki section; do not proceed.

### Gate 3 — Process Previous Studies by Paper

1. Identify every unique Raw paper in generated `Meanings`, preserving source order.
2. Apply [Process-One-Paper-Rules](../wiki-rules/SKILL.md#process-one-paper-rules) to exactly one unreviewed paper.
3. Proceed to Gate 4 only after every identified paper passes.

### Gate 4 — Wiki to Tags

1. Begin only after `Previous Studies` passes Gate 3.
2. Reconcile vocab under [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules), requiring zero unsupported Wiki claims.
3. Present exactly one smallest vocab change. **🔴 CHECKPOINT · 🛑 STOP** — Await separate explicit approval; do not proceed.
4. Apply only the approved change, complete Vocab verification, and show `git diff -- raw wiki vocab`.
