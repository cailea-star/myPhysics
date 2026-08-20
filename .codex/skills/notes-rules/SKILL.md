---
name: notes-rules
description: Apply canonical myWIKI notes rules when creating, drafting, reviewing, or updating personal notes and sections under notes/, including English prose, heading levels, mathematical formatting, symbol definitions, canonical links, and approval-gated writing.
---

# Notes Rules

## Writing Rules

### Text Rules

- Write note prose in English.
- Use `#` for the document title and `###` for topic sections.
- Present drafts as rendered Markdown; NEVER wrap note drafts in fenced code blocks.
- Make equations primary and prose secondary; use prose only to define symbols, state assumptions, and connect equations.
- Prefer short phrase lists for parallel properties; avoid long prose paragraphs.
- Maintain general theory in exactly one canonical note. Link to that note instead of repeating its derivation elsewhere.

### Concision Rules

- Remove unnecessary prose.
- Merge only consecutive formulas forming one short calculation; preserve distinct steps and useful itemized structures.
- Remove explanations outside the main line.
- Refine one section at a time; present its draft before writing.

### Formula Rules

- Use `$...$` for inline mathematics. For display mathematics, place the opening and closing `$$` delimiters on separate lines.
- Keep the formula between the delimiters on a single line, regardless of length. Before using `\begin{aligned}` or any multiline formula, ask the user and wait for explicit approval.
- Use `\boxed{...}` to emphasize core formulas; do not box intermediate derivation steps.
- Give each symbol exactly one meaning. Define every symbol before its first use.

## Section Workflow

Skip approved stages and resume from the relevant stage.

1. **Reference Confirmation:** Confirm and read all reference materials before planning.
2. **Section Structure Confirmation:** Confirm the document `#` title and the `###` sections in their order.
3. **Subsection Outline Confirmation:** Confirm each `###` section's scope and formula sequence.
4. **Draft Confirmation:** Revise the complete subsection draft until explicitly approved.
5. **Write:** Reread the target file, write only approved content, and verify formatting.
