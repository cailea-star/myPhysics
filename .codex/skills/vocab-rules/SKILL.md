---
name: vocab-rules
description: Apply canonical myWIKI vocabulary rules when selecting, drafting, reviewing, adding, merging, or reusing tags, types, aliases, or authors in vocab files, raw quotations, ingest, or lint workflows.
---

# Vocab Rules

## Rules

### Tag-Rules

Scope: Quotation tags are centered on physical quantities actually calculated, plotted, or compared in the paper; do not tag unrelated background mentions.

Types: Use tag types from [vocab/types.json](../../../vocab/types.json), tags from [vocab/tags.json](../../../vocab/tags.json), and authors from [vocab/authors.json](../../../vocab/authors.json).

Type Match: Every type assigned to a drafted tag MUST independently satisfy that tag type's `requirement` in [vocab/types.json](../../../vocab/types.json); otherwise remove the type or reject the draft.

Granularity: Keep noun-term compounds as precise tags, but split adjective-like modifiers into `[tag_type]: property` tags; e.g. use `alpha_decay_energy`, but use `symmetry_energy` + `soft`.

### Draft-Rules

Draft: If a needed tag is missing, draft `{tag, definition, types, aliases}` entries for [vocab/tags.json](../../../vocab/tags.json) one tag type at a time; include only aliases seen in the paper, metadata, or existing project vocabulary.

Similarity Check: Before drafting each new tag, run `python scripts\search_similar_tags.py "CANDIDATE" 3` and show all three results. Prefer reusing existing tags.

Write: Present one draft group with an explicit add, merge, or reuse proposal for each candidate. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Write only individually approved entries, then present the next draft group.

