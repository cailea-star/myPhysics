---
name: lint
description: Audit and update myWIKI from tagged quotations. Use when checking wiki coverage, backfilling missing evidence into raw, or reconciling wiki conclusions with vocab tags.
---

# Lint

Before Gate 1, read [Vocab Rules](../vocab-rules/SKILL.md), [Raw Rules](../raw-rules/SKILL.md), and [Wiki Rules](../wiki-rules/SKILL.md) completely and apply them throughout every gate.

## Gated Workflow

Run gates strictly in order. Start every response with the current gate, last completed gate, and next approval. NEVER advance from an unapproved or failed section or gate.

### Gate 1 — Confirm Tag & Prepare Inputs

- Check: Run `git rev-parse --is-inside-work-tree` and `git status --short`; failures or changes MUST STOP and report.
- Prepare: MUST execute every [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) Prepare step for one canonical `TAG`; failure MUST STOP.
- Mode: Existing `wiki/TAG.md` means update; absence MUST STOP until explicit creation approval under [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules).

### Gate 2 — Audit Wiki Sections

- Review: MUST process exactly one Wiki section through every [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) step in template order.
- Approval: After each section passes, PRINT `git diff -- raw wiki`; STOP until approval for next section.
- Advance: Gate 3 begins ONLY after every Wiki section passes.

### Gate 3 — Wiki to Tags

- Reconcile: After Gate 2 passes, MUST execute [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules) completely until every Wiki tag resolves canonically.
- Complete: ONLY after Vocab verification passes, PRINT `git diff -- raw wiki vocab`.
