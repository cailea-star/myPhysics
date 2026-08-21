---
name: skill-refactor
description: Refactor rule-heavy SKILL.md files without losing semantics or constraints. Use when simplifying structure, execution flow, ownership, or duplication in a skill.
---

# Skill Refactor

## Refactoring Constraints

- Separate mechanical moves from semantic edits; NEVER shorten commands, paths, formats, or required semantics.
- Format: Checkpoints, commands, paths, fields, and templates MUST NOT change.
- Modal force MUST remain; NEVER change `MUST`, `NEVER`, or `ONLY` without approval.
- Align Rules with Gates; maximum three subsections, five rules each, fifteen words each.
- Report gaps; NEVER write semantic changes before approval; verify approved scope afterward.

## Workflow

1. Confirm the minimal workflow with the user; select exactly one current gate.
2. Present the current gate's existing rules; propose at most three changes using the required format.
3. Obtain the user's decision; semantic changes MUST receive explicit approval.
4. Apply approved deprecations, merges, or deletions; complete post-write review.

### Required Comparison Format

````text
Workflow: ...

Current Gate: ...

Current Section：## Parent → ### Target

Current Section Content：
```text
...
```

Original1：...
Modified1：...
Reason1：...

Original2：...
Modified2：...
Reason2：...
````
