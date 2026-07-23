---
name: vocab-rules
description: Apply canonical myWIKI vocabulary rules when selecting, drafting, reviewing, adding, merging, or reusing tags, types, aliases, or authors in vocab files, raw quotations, ingest, or lint workflows.
---

# Vocab Rules

## Author-Rules

- Authors: [vocab/authors.json](../../../vocab/authors.json) stores corresponding authors’ `given`, `family`, `affiliation`, and `ORCID`.
- Check: For one `raw/*.json`, MUST run `python scripts\check_vocab_author.py raw\[json_filename].json`; report additions and unchanged authors.

## Tag-Rules

- Tags: [vocab/tags.json](../../../vocab/tags.json) defines canonical concepts, types, aliases, and definitions for Raw and Wiki classification.
- Types: [vocab/types.json](../../../vocab/types.json) defines tag types and their requirements.

### Input-Rules

- `raw/*.pdf` and `raw/*.tex` are Tag-Inputs; `raw/*.md` is NEVER Tag-Input.
- ONLY technical terms stated in Tag-Inputs and absent from canonical tags may become new-tag candidates.
- Definitions MUST be determined ONLY from Tag-Inputs.
- Aliases MUST come ONLY from actual Tag-Inputs usage; NEVER invent aliases.
- PRINT: If Tag-Inputs cannot determine a definition, MUST STOP and ask the user.

### Resolve-Rules

- For every candidate, MUST run `python scripts\search_similar_tags.py "CANDIDATE" 3` and show all three results.
- Resolve ONLY from returned canonical tags, types, aliases, and definitions; NEVER load [vocab/tags.json](../../../vocab/tags.json) directly.

### Draft-Rules

- Draft: PRINT tag-type requirement; read ONLY final 15 [vocab/tags.json](../../../vocab/tags.json) lines; present exact proposals.
- Validate: Every tag MUST satisfy its tag-type requirement; assign `reuse`, `add`, `merge`, or `rename`.
- Approve: STOP until explicit user approval.
- Write: Apply ONLY approved changes across [vocab/tags.json](../../../vocab/tags.json), `raw/*.md`, and Wiki files/links; merge/rename MUST run `python scripts\rename_raw_tag.py OLD NEW`.
- Verify: Run `python scripts\search_a_tag.py TAG`; OLD MUST vanish from [vocab/tags.json](../../../vocab/tags.json), `raw/*.md`, and Wiki files/links; otherwise redraft.
