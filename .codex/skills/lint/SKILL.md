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
- Never infer evidence, fill a template without evidence, or modify `raw`, `wiki`, or `vocab` before approval.

### Raw Backfill Rules

- Quote exact complete sentence(s) from `raw/*.pdf` or `raw/*.tex`, normally 10–40 words; include `[claim_type]`, `[tags]`, and `[source]`.
- Use at least two directly supported canonical tags. Results MUST include a method and quantity tag; Meanings MUST include a method tag.
- Put externally attributed claims in Secondary Citations. Pair any core formula with its explanatory sentence and a fenced `math` block.
- Show the final exact-mode block and obtain approval before writing `raw/*.md`.

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

For every `weak` or `missing` verdict, identify candidate sources or report an evidence gap. New quotations MUST follow Raw Backfill Rules. After writing, run:

```powershell
python scripts\sort_raw_md_quotations.py raw\PAPER.md
python scripts\search_a_tag.py TAG
```

The new quotation MUST appear in generated `tmp`; then reopen its Gate 2 section. Otherwise the section fails.

### Vocab Rules

- Reconcile vocab only after wiki approval. Definition defines the tag; aliases MUST occur in source, metadata, or existing vocab.
- Before every add, merge, or rename, show three `search_similar_tags.py` results and obtain separate approval.
- For an approved merge or rename, run `python scripts\rename_raw_tag.py OLD NEW`, update `tags.json` and affected wiki names/links in the same gate, then require zero OLD quotations and re-audit NEW.

## Gated Workflow

Run gates strictly in order. Start every response with current gate, last completed gate, and next approval. Never advance from an unapproved section or gate. If raw changes after search, rerun `search_a_tag` and restart the affected Gate 2 section.

### Gate 1 — Confirm Tag

1. Run:
   ```powershell
   git rev-parse --is-inside-work-tree
   git status --short
   ```
2. If the worktree is not clean, stop and report existing changes.
3. Confirm the canonical tag, types, wiki path, and required template under Truth Rules.
4. If the tag is absent from `vocab/tags.json`, stop.
5. If wiki is absent, obtain approval and create it from the required template before Gate 2.

### Gate 2 — Audit Quotations

1. Run:
   ```powershell
   python scripts\search_a_tag.py TAG
   ```
2. If the command fails or any expected `tmp` section is absent, stop; never audit partial output.
3. Audit exactly one generated section per response, in order: Motivation, Methods, Results, Meanings, Secondary Citations.
4. Apply Quotation Verdict Rules; stop after each section and advance only after approval.

### Gate 3 — Audit Wiki Gaps

1. Audit exactly one wiki section per response, in template order.
2. Apply Wiki Verdict Rules; stop after each section and advance only after approval.

### Gate 4 — Update Wiki

1. Apply only changes approved in Gates 2 and 3.
2. Add no unsupported claim and no change outside the approved set.
3. Show:
   ```powershell
   git diff -- wiki raw
   ```
4. Stop and obtain wiki-update approval.

### Gate 5 — Reconcile Vocab

1. Check only the target tag's definition, aliases, and types under Vocab Rules.
2. Run before every add, merge, or rename:
   ```powershell
   python scripts\search_similar_tags.py "CANDIDATE" 3
   ```
3. Obtain separate approval, then apply the complete approved vocab change under Vocab Rules.

### Gate 6 — Verify

1. Rerun `python scripts\search_a_tag.py TAG`.
2. Verify every added quotation, reference, wiki entry, and vocab change; references MUST satisfy Truth Rules.
3. Require zero unprocessed quotations, unresolved approved changes, unsupported wiki claims, or partial tag migrations.
4. Show:
   ```powershell
   git diff -- raw wiki vocab
   ```
5. Report evidence gaps without expanding scope.
