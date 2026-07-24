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
- Prepare: MUST execute [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) Prepare, including `python scripts\search_a_tag.py TAG`; generated `tmp/*` quotations MUST pass validation.
- Mode: Existing `wiki/TAG.md` means update; absence MUST STOP until explicit creation approval under [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules).

### Gate 2 — Audit Wiki Sections

- Review: MUST execute every [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) section-mode step for exactly one section in template order.
- Approval: Each pass MUST PRINT `git diff -- raw wiki`; STOP until approval; ONLY all-section completion advances.

### Gate 3 — Scan Papers

- Scan: MUST queue every Raw indexed by `tmp/*_Meanings.md`; missing or duplicate mappings MUST STOP.
- Review: MUST execute every [Wiki Draft-Rules](../wiki-rules/SKILL.md#draft-rules) paper-mode step for exactly one queued Raw.
- Approval: Each pass MUST PRINT `git diff -- raw wiki`; STOP until approval; ONLY empty queue advances.

### Gate 4 — Wiki to Tags

- Reconcile: After Gate 3 passes, MUST execute [Vocab Draft-Rules](../vocab-rules/SKILL.md#draft-rules) completely until every Wiki tag resolves canonically.
- Complete: ONLY after Vocab verification passes, PRINT `git diff -- raw wiki vocab`.
