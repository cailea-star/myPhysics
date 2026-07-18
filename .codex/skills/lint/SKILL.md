---
name: lint
description: Audit and update myWIKI from tagged quotations. Use when checking wiki coverage, backfilling missing evidence into raw, or reconciling wiki conclusions with vocab tags.
---

# Lint

## Rules

### Truth Rules

- `raw/*.pdf` and `raw/*.tex` are source truth; `raw/*.md` is quotation evidence; `wiki/*.md` is synthesis; `vocab/tags.json` is canonical schema; never edit generated `tmp/*.md`.
- A wiki filename MUST equal its canonical tag. A `method` uses `scripts/add_wiki_method.md`; every other type uses `scripts/add_wiki_topic.md`.
- Wiki MUST synthesize evidence; never copy a quotation into wiki prose.
- Every factual claim and formula MUST have an adjacent reference. Reference MUST be an existing `raw/*` filename without path or suffix: `raw/PAPER.md` -> `PAPER`.
- Every new `raw/*.md` quotation MUST follow the Section-Rules, Quotation-Rules, and applicable section-specific rules in `.codex/skills/ingest/SKILL.md` and receive approval before writing.
- Never infer evidence, fill a template without evidence, or modify `raw`, `wiki`, or `vocab` before approval.

### Quotation Verdict Rules

Give every quotation exactly one verdict:

- `covered`: wiki states the conclusion and cites the paper; do nothing.
- `miss-reference`: conclusion exists but lacks the paper; add only its reference.
- `miss-conclusion`: an important conclusion is absent; add the smallest supported entry and reference.
- `conflict`: quotation conflicts with wiki; preserve both conclusions, conditions, and references.
- `irrelevant`: quotation does not belong in this tag wiki; do not add it.

Every quotation MUST receive exactly one verdict; otherwise the section fails. Propose the smallest change for every non-`covered` verdict; never apply it before approval.

### Wiki Verdict Rules

Give every factual claim, formula, and required template item exactly one verdict:

- `supported`: evidence is sufficient; do nothing.
- `weak`: content exists but evidence is insufficient; find supporting quotations.
- `missing`: required content is absent; find candidate papers and quotations.
- `not-applicable`: the template item does not apply; do not fill it.

For every `weak` or `missing` verdict, identify candidate sources or report an evidence gap. After writing an approved quotation, run:

```powershell
python scripts\sort_raw_md_quotations.py raw\PAPER.md
python scripts\search_a_tag.py TAG
```

The new quotation MUST appear in generated `tmp` and receive a Quotation Verdict before entering wiki; otherwise the section fails.

### Vocab Rules

- Reconcile vocab only after wiki approval. Definition defines the tag; aliases MUST occur in source, metadata, or existing vocab.
- Before every add, merge, or rename, show three `search_similar_tags.py` results and obtain separate approval.
- For an approved merge or rename, run `python scripts\rename_raw_tag.py OLD NEW`, update `tags.json` and affected wiki names/links in the same gate, then require zero OLD quotations and re-audit NEW.

## Gated Workflow

Run gates strictly in order. Start every response with current gate, last completed gate, and next approval. Never advance from an unapproved section or gate. After any raw change, rerun `search_a_tag` and re-audit the affected quotation section before continuing.

### Gate 1 — Confirm Tag & Search

1. Run:
   ```powershell
   git rev-parse --is-inside-work-tree
   git status --short
   ```
2. If the worktree is not clean, stop and report existing changes.
3. Confirm the canonical tag, types, wiki path, and required template under Truth Rules.
4. If the tag is absent from `vocab/tags.json`, stop.
5. If wiki is absent, **🔴 CHECKPOINT · 🛑 STOP**: do not search, audit, or propose changes; obtain approval to create it from the required template.
6. Run `python scripts\search_a_tag.py TAG`; if it fails or any expected `tmp` section is absent, stop.

### Gate 2 — Quotations to Wiki

1. Audit exactly one generated section per response, in order: Motivation, Methods, Results, Meanings, Secondary Citations.
2. Apply Quotation Verdict Rules to every quotation.
3. Present the section's smallest wiki change set and stop for approval.
4. Apply only approved references and entries; show `git diff -- wiki`.
5. Advance only after approval of the updated section.

### Gate 3 — Wiki to Raw to Wiki

1. Audit exactly one wiki section per response, in template order.
2. Apply Wiki Verdict Rules to every claim, formula, and required item.
3. For each `weak` or `missing`, find exact source evidence under Truth Rules; report a gap if none exists.
4. Review every proposed raw quotation under the referenced ingest rules with `pass`, `gap`, or `fix`; only `pass` may enter the change set.
5. Present the section's raw-and-wiki change set and stop for approval.
6. After approval, write raw, renumber, rerun `search_a_tag`, and give every new quotation a Quotation Verdict.
7. Only verified evidence may enter wiki; show `git diff -- raw wiki`.
8. Advance only after approval of the updated section.

### Gate 4 — Wiki to Tags

1. Check the target tag's definition, aliases, and types against the approved wiki and source evidence.
2. Before every add, merge, or rename, run:
   ```powershell
   python scripts\search_similar_tags.py "CANDIDATE" 3
   ```
3. Show all three results and the smallest vocab change; obtain separate approval for each change.
4. Apply approved changes under Vocab Rules.
5. Rerun `search_a_tag` for all affected tags and require zero unprocessed quotations, unresolved approved changes, unsupported wiki claims, or partial tag migrations.
6. Show:
   ```powershell
   git diff -- raw wiki vocab
   ```
7. Report evidence gaps without expanding scope.
