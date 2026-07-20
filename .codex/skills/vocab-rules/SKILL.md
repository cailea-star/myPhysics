---
name: vocab-rules
description: Apply canonical myWIKI vocabulary rules when selecting, drafting, reviewing, adding, merging, or reusing tags, types, aliases, or authors in vocab files, raw quotations, ingest, or lint workflows.
---

# Vocab Rules

## Rules

### Resolve-Rules

Scope: Quotation tags are centered on physical quantities actually calculated, plotted, or compared in the paper; do not tag unrelated background mentions.

Sources: Use tag types from [vocab/types.json](../../../vocab/types.json) and authors from [vocab/authors.json](../../../vocab/authors.json).

Lookup: For every candidate concept, run `python scripts\search_similar_tags.py "CANDIDATE" 3`, show all three results, and use their canonical tags, types, aliases, and definitions; NEVER load `vocab/tags.json` directly.

Type Match: Every type assigned to a drafted tag MUST independently satisfy that tag type's `requirement` in [vocab/types.json](../../../vocab/types.json); otherwise remove the type or reject the draft.

Granularity: Keep noun-term compounds as precise tags, but split adjective-like modifiers into `[tag_type]: property` tags; e.g. use `alpha_decay_energy`, but use `symmetry_energy` + `soft`.

### Draft-Rules

Proposal: Under Resolve-Rules, propose `reuse`, `add`, `merge`, or `rename` for each candidate. If no valid existing tag resolves, draft `{tag, definition, types, aliases}` for [vocab/tags.json](../../../vocab/tags.json) one tag type at a time; include only aliases seen in the paper, metadata, or existing project vocabulary.

Merge: Run `python scripts\search_similar_tags.py "TARGET" 3`, require the exact target in its output, and check its definition, aliases, and types against approved wiki and source evidence.

Write: Present one draft group with an explicit proposal for each candidate. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Apply only individually approved proposals: write approved additions; for an approved merge or rename, run `python scripts\rename_raw_tag.py OLD NEW` and update `tags.json` and affected wiki names and links in the same gate. After any required verification passes, present the next draft group.

Verify: After every approved add, merge, or rename, rerun `python scripts\search_a_tag.py TAG` for every affected tag; require zero unprocessed quotations, unresolved approved changes, or partial tag migrations. For merge or rename, also require zero OLD quotations. Re-audit the resulting tag and report evidence gaps without expanding scope.
