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
- Resolve: MUST run `python scripts\search_similar_tags.py "CANDIDATE" 3` for each candidate and resolve ONLY from its three printed tags, types, aliases, and definitions; NEVER read [vocab/tags.json](../../../vocab/tags.json) directly.

### Input-Rules

- Existing: Every existing tag MUST run `python scripts\search_a_tag.py TAG` as its modification reference.
- Unknown: `raw/*.pdf` and `raw/*.tex` are Tag-Inputs; `raw/*.md` is NEVER Tag-Input.
- Candidate: ONLY technical terms stated in Tag-Inputs and absent from canonical tags may become new-tag candidates.
- Definition: Definitions MUST be determined ONLY from Tag-Inputs; unresolved definitions MUST STOP and ask the user.
- Aliases MUST come ONLY from actual Tag-Inputs usage; NEVER invent aliases.

### Edit-Rules

- Scope: Handle ONLY `add`, `merge`, and `rename`; `reuse` MUST NOT enter Edit-Rules.
- PRINT: MUST read [vocab/types.json](../../../vocab/types.json) and ONLY final 15 [vocab/tags.json](../../../vocab/tags.json) lines; show every affected tag-type requirement.
- Input: MUST execute every required step in [Input-Rules](#input-rules).
- Write: STOP until explicit approval; apply ONLY exact approved changes; merge/rename MUST run `python scripts\rename_raw_tag.py OLD NEW`.
- Verify: After EVERY write, NEW MUST exist and OLD MUST vanish everywhere; failures MUST STOP and restart Edit-Rules from Input; ONLY passing verification may RETURN.

### Draft-Rules

- Resolve: MUST execute the [Tag-Rules](#tag-rules) Resolve action for every candidate; unresolved candidates MUST STOP.
- Label: PRINT exact proposals; assign every candidate exactly one `reuse`, `add`, `merge`, or `rename`.
- Route: ONLY when every label is `reuse`, PRINT mappings and RETURN; otherwise pass ONLY `add`, `merge`, and `rename` proposals to [Edit-Rules](#edit-rules).
