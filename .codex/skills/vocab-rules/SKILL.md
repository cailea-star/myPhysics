---
name: vocab-rules
description: Apply canonical myWIKI vocabulary rules when selecting, drafting, reviewing, adding, merging, or reusing tags, types, aliases, or authors in vocab files, raw quotations, ingest, or lint workflows.
---

# Vocab Rules

## Rules

### Tag-Rules

Scope: Quotation tags are centered on physical quantities actually calculated, plotted, or compared in the paper; do not tag unrelated background mentions.

Types: Use tag types from [vocab/types.json](../../../vocab/types.json) and authors from [vocab/authors.json](../../../vocab/authors.json). For every tag lookup, run `python scripts\search_similar_tags.py "QUERY" 3` and use its canonical tag, types, aliases, and definition output; NEVER load `vocab/tags.json` directly.

Type Match: Every type assigned to a drafted tag MUST independently satisfy that tag type's `requirement` in [vocab/types.json](../../../vocab/types.json); otherwise remove the type or reject the draft.

Granularity: Keep noun-term compounds as precise tags, but split adjective-like modifiers into `[tag_type]: property` tags; e.g. use `alpha_decay_energy`, but use `symmetry_energy` + `soft`.

### Verify-Rules

After every approved add, merge, or rename, rerun `python scripts\search_a_tag.py TAG` for every affected tag; require zero unprocessed quotations, unresolved approved changes, or partial tag migrations. Re-audit the resulting tag and report evidence gaps without expanding scope.


### Draft-Rules

Draft: If the script finds no valid existing tag, draft `{tag, definition, types, aliases}` entries for [vocab/tags.json](../../../vocab/tags.json) one tag type at a time; include only aliases seen in the paper, metadata, or existing project vocabulary.

Resolve: For each candidate concept, run `python scripts\search_similar_tags.py "CANDIDATE" 3`, show all three results, and propose `reuse`, `merge`, or `add`.

Write: Present one draft group with an explicit add, merge, or reuse proposal for each candidate. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed. Write only individually approved entries, verify every affected tag under [Verify-Rules](#verify-rules), then present the next draft group.

### Merge-Rules

Target: Run `python scripts\search_similar_tags.py "TARGET" 3`, require the exact target in its output, and check the returned definition, aliases, and types against the approved wiki and source evidence.

Apply: For an approved merge or rename, run `python scripts\rename_raw_tag.py OLD NEW`; update `tags.json` and affected wiki names and links in the same gate.

Verify: Verify OLD and NEW under [Verify-Rules](#verify-rules) and require zero OLD quotations.

