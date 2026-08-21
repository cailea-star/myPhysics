---
name: skill-refactor
description: Refactor rule-heavy SKILL.md files without losing semantics or constraints, strengthening modal force where possible. Use for structural, workflow, ownership, duplication, or normative-language changes.
---

# Skill Refactor

## Refactoring Constraints

- Separate mechanical moves from semantic edits; NEVER shorten commands, paths, formats, or required semantics.
- Format: Checkpoints, commands, paths, fields, and templates MUST NOT change.
- Strengthen modal force: requirements use `MUST`, limits `ONLY`, prohibitions `NEVER`; semantic changes MUST receive approval.
- Align Rules with Gates; maximum three subsections, five rules each, fifteen words each.
- Report gaps; NEVER write semantic changes before approval; verify approved scope afterward.

## Workflow

1. Confirm the minimal workflow with the user; select exactly one current gate.
2. NEVER add entries by default; additions MUST show full subsection; propose ≤3 changes.
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
