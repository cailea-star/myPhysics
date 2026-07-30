---
name: notes-rules
description: Apply canonical myWIKI notes rules when creating, drafting, reviewing, or updating personal notes and sections under notes/, including English prose, heading levels, mathematical formatting, symbol definitions, canonical links, and approval-gated writing.
---

# Notes Rules

## Writing Rules

### Text Rules

- Write note prose in English.
- Use `###` for topic sections and `#####` for internal blocks.
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

Work on exactly one note section at a time.

1. Ask whether reference documents exist. If provided, read them before drafting.
2. Require the user to determine the section title before drafting. If no title has been determined, ask for it and stop.
3. Present the complete proposed section directly in the conversation without a fenced code block. Do not write it to a file.
4. Wait for explicit user approval. Treat comments, corrections, and layout changes as requests for a revised draft, not as approval to write.
5. After approval, reread the latest target file, write only the approved section, and preserve unrelated user changes.
