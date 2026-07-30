---
name: skill-refactor
description: Refactor rule-heavy SKILL.md files without losing semantics or constraints. Use when simplifying structure, execution flow, ownership, or duplication in a skill.
---

# Skill Refactor

## Workflow

1. Confirm the minimal gate workflow with the user.
2. Present one gate's core cautions at a time; use at most five items.
3. Map current rules to each gate; show the user differences and redundancies against core cautions.
4. Use gate-led decisions to deprecate, merge, or delete rules.

## Checks

- Separate mechanical moves from semantic edits; NEVER shorten commands, paths, formats, or required semantics.
- Format: Checkpoints, commands, paths, fields, and templates MUST NOT change.
- Constraints: Preserve and prefer `MUST`, `NEVER`, and `ONLY`.
- Align static Rules with Gates; allow three subsections, five rules each, fifteen words per rule.
- Strengthen constraints; close gaps; require inputs → validation → STOP → approval → write boundaries → post-write review.

## Output

Before proposing or applying changes, select exactly one current section, print the one-way workflow once, then show one exact minimal comparison per independent change:

```text
Workflow: ...

Current Section：## Parent → ### Target

Original1：...
Modified1：...
Reason1：...

Original2：...
Modified2：...
Reason2：...
```
