---
name: notes-rules
description: Apply canonical myWIKI notes rules when creating, drafting, reviewing, or updating personal notes and sections under notes/, including English prose, heading levels, mathematical formatting, symbol definitions, canonical links, and approval-gated writing.
---

# Notes Rules

## Writing Rules

- Write note prose in English.
- Use `###` for topic sections and `#####` for internal blocks.
- Use `$...$` for inline mathematics and `$$...$$` for display mathematics. Never use `\begin{aligned}`.
- Give each symbol exactly one meaning. Define every symbol before its first use.
- Maintain general theory in exactly one canonical note. Link to that note instead of repeating its derivation elsewhere.

## Section Workflow

Work on exactly one note section at a time.

1. Require the user to determine the section title before drafting. If no title has been determined, ask for it and stop.
2. Present the complete proposed section directly in the conversation without a fenced code block. Do not write it to a file.
3. Wait for explicit user approval. Treat comments, corrections, and layout changes as requests for a revised draft, not as approval to write.
4. After approval, reread the latest target file, write only the approved section, and preserve unrelated user changes.
