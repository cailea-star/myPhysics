---
name: skill-refactor
description: Refactor rule-heavy SKILL.md files without losing semantics or constraints. Use when simplifying structure, execution flow, ownership, or duplication in a skill.
---

# Skill Refactor

## Workflow

1. Preserve key meaning, constraint strength, fixed formats, and scope.
2. Assign each rule one owner: truth, schema, validity, review, or execution.
3. Make the workflow explicit: prepare → validate → approve → write → verify → pass or fix.
4. Deduplicate only after every path applies the surviving owner.
5. Restructure last; keep few main sections and at most three subsections each.

## Checks

- Separate mechanical moves from semantic edits.
- Preserve `MUST`, `NEVER`, `ONLY`, checkpoints, commands, paths, fields, and templates.
- Check non-heading content, links, callers, and pre-write/post-write gates.

## Output

Before proposing or applying changes, print the resulting one-way workflow once, then one exact minimal comparison per independent change:

```text
Workflow: ...

Original1：...
Modified1：...
Reason1：...

Original2：...
Modified2：...
Reason2：...
```
