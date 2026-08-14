---
name: papers-vocab-rules
description: Apply canonical myWIKI vocabulary rules when selecting, drafting, reviewing, adding, updating, merging, or reusing tags, types, aliases, or authors in vocab files, raw quotations, papers-ingest, or papers-lint workflows.
---

# Papers Vocab Rules

## Author-Rules

- Authors: [papers/vocab/authors.json](../../../papers/vocab/authors.json) stores corresponding authors’ `given`, `family`, `affiliation`, and `ORCID`.
- Check: For one `papers/raw/*.json`, MUST run `python papers\scripts\check_vocab_author.py papers\raw\[json_filename].json`; report additions and unchanged authors.

## Tag-Rules

- Tags: [papers/vocab/tags.json](../../../papers/vocab/tags.json) defines canonical concepts, types, aliases, and definitions for Raw and Wiki classification.
- Types: [papers/vocab/types.json](../../../papers/vocab/types.json) defines tag types and their requirements.
- Resolve: MUST run `python papers\scripts\search_similar_tags.py "CANDIDATE" 10`; use ONLY top-three records; NEVER read [papers/vocab/tags.json](../../../papers/vocab/tags.json).
- Alias-Resolve: MUST run `python papers\scripts\search_similar_tags.py "ALIAS" 10` for EVERY alias; MUST match ONLY selected-tag.

### Input-Rules

- Existing: MUST run `python papers\scripts\search_a_tag.py TAG` and review EVERY output for each existing tag.
- Unknown: `papers/raw/*.pdf` and `papers/raw/*.tex` are Tag-Inputs; `papers/raw/*.md` is NEVER Tag-Input.
- Candidate: ONLY Tag-Input technical terms absent from canonical tags qualify for `add`.
- Definition: MUST fit candidate/type; changes use ONLY Tag-Input evidence/location; failures MUST STOP; ask user.
- Aliases: ONLY verbatim Tag-Input synonyms pass; missing persisted aliases MUST STOP; ask `update` deletion.

### Edit-Rules

- Scope: ONLY `update`, `merge`, `add` enter; `reuse` MUST NEVER enter; `update` MUST preserve `tag`.
- PRINT: MUST read [papers/vocab/types.json](../../../papers/vocab/types.json) and ONLY final 15 [papers/vocab/tags.json](../../../papers/vocab/tags.json) lines; show EVERY affected requirement.
- Write: STOP until approval; apply ONLY approved-content; merge MUST run `python papers\scripts\rename_raw_tag.py OLD NEW`.
- Verify: After EVERY write, persisted `tag`/`definition`/`types`/`aliases` MUST match approval; mismatches MUST STOP; restart Input.
- Identity: MUST enforce add(absent→exists), update(exists→exists), merge(OLD+NEW→NEW; OLD absent everywhere); failures MUST STOP; restart Input.

### Draft-Rules

- Resolve: MUST execute the [Tag-Rules](#tag-rules) Resolve action for every candidate; unresolved candidates MUST STOP.
- Evidence: EVERY proposal, including `reuse`, MUST execute [Input-Rules](#input-rules).
- Review: PRINT exact `definition`/`types`/`aliases` with evidence/locations/similarities and one `pass|fix`; `fix` MUST STOP before approval.
- Label: For EVERY proposal, PRINT exact `tag`/`definition`/`types`/`aliases`; assign EXACTLY one: unchanged=`reuse`, same-tag-change=`update`, duplicates=`merge`, absent=`add`.
- Route: PRINT/RETURN mappings ONLY when all=`reuse`; otherwise send ONLY `update|merge|add` to [Edit-Rules](#edit-rules).
