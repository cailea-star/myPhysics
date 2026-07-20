---
name: lint
description: Audit and update myWIKI from tagged quotations. Use when checking wiki coverage, backfilling missing evidence into raw, or reconciling wiki conclusions with vocab tags.
---

# Lint

## Rules

### Truth Rules

- `raw/*.pdf` and `raw/*.tex` are source truth; `raw/*.md` is quotation evidence; `wiki/*.md` is synthesis; `vocab/tags.json` is canonical schema; never edit generated `tmp/*.md`.
- Every new `raw/*.md` quotation MUST follow ingest [Claim-Type-Rules](../ingest/SKILL.md#claim-type-rules), [Section-Rules](../ingest/SKILL.md#section-rules), [Quotation-Rules](../ingest/SKILL.md#quotation-rules), and applicable section-specific rules, and receive approval before writing.
- Never infer evidence, fill a template without evidence, or modify `raw`, `wiki`, or `vocab` before approval.

### Quotation Verdict Rules

Give every quotation exactly one verdict:

- `covered`: wiki states the conclusion and cites the paper; do nothing.
- `miss-reference`: conclusion exists but lacks the paper; add only its reference.
- `miss-conclusion`: an important conclusion is absent; add the smallest supported entry and reference.
- `conflict`: quotation conflicts with wiki; preserve both conclusions, conditions, and references.
- `irrelevant`: quotation does not belong in this tag wiki; do not add it.

List one verdict per `quote` block in source order and report `Verdicts: <verdicts>/<quote blocks>`. If counts differ, the section fails and no change set may be proposed; verdicts may share a change but MUST NOT be merged. Propose the smallest change for every non-`covered` verdict; never apply it before approval.

### Wiki-Rules

- A wiki filename MUST equal its canonical tag. A `method` uses `scripts/add_wiki_method.md`; every other type uses `scripts/add_wiki_topic.md`.
- Wiki MUST synthesize evidence; never copy a quotation into wiki prose.
- Every factual claim and formula MUST have an adjacent reference. Reference MUST be an existing `raw/*` filename without path or suffix: `raw/PAPER.md` -> `PAPER`.
- Secondary Citations are source leads, not wiki evidence; never use their quoted claims or formulas to update wiki or backfill raw.

Give every factual claim, formula, and required template item exactly one verdict:

- `supported`: evidence is sufficient; do nothing.
- `weak`: content exists but evidence is insufficient; find supporting quotations.
- `missing`: required content is absent; find candidate papers and quotations.
- `not-applicable`: the template item does not apply; do not fill it.

List one verdict per factual claim, formula, and required template item and report `Wiki Verdicts: <verdicts>/<claims + formulas + required items>`. A section passes ONLY when counts match and every `weak` or `missing` is resolved or recorded as an explicit no-evidence gap; otherwise no change set may be proposed.

For each `weak` or `missing`, search primary raw evidence first. If none exists, inspect `tmp/TAG_Secondary.md` only for cited-paper leads: use an already-ingested original paper's direct evidence; otherwise propose ingesting it or report an explicit gap.

Review every proposed raw quotation under ingest [Claim-Type-Rules](../ingest/SKILL.md#claim-type-rules), [Section-Rules](../ingest/SKILL.md#section-rules), and [Quotation-Rules](../ingest/SKILL.md#quotation-rules) with `pass`, `gap`, or `fix`; only `pass` may enter the change set.

After any approved raw change, run:

```powershell
python scripts\sort_raw_md_quotations.py raw\PAPER.md
python scripts\search_a_tag.py TAG
```

Re-audit the affected quotation section. A new quotation MUST appear in generated `tmp` and receive a Quotation Verdict before entering wiki; otherwise the section fails. Only verified evidence may enter wiki.

### Tag-Rules

- Every vocab proposal MUST follow ingest [Tag-Rules](../ingest/SKILL.md#tag-rules).
- Reconcile vocab only after wiki approval; check the target definition, aliases, and types against the approved wiki and source evidence.
- For an approved merge or rename, run `python scripts\rename_raw_tag.py OLD NEW`, update `tags.json` and affected wiki names/links in the same gate, then require zero OLD quotations and re-audit NEW.
- After every approved vocab change, rerun `search_a_tag` for all affected tags and require zero unprocessed quotations, unresolved approved changes, unsupported wiki claims, or partial tag migrations. Report evidence gaps without expanding scope.

## Gated Workflow

Run gates strictly in order. Start every response with current gate, last completed gate, and next approval. Never advance from an unapproved section or gate.

### Gate 1 — Confirm Tag & Search

1. Run; if either check fails or the worktree is dirty, **🛑 STOP** and report:
   ```powershell
   git rev-parse --is-inside-work-tree
   git status --short
   ```
2. Confirm the canonical tag, types, wiki path, and required template under [Truth Rules](#truth-rules) and [Wiki-Rules](#wiki-rules); if the tag is absent from `vocab/tags.json`, **🛑 STOP**.
3. If wiki is absent, do not search, audit, or propose changes. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval to create it from the required template; do not proceed.
4. Run `python scripts\search_a_tag.py TAG`; if it fails or any expected `tmp` section is absent, **🛑 STOP**.

### Gate 2 — Quotations to Wiki

1. Apply [Quotation Verdict Rules](#quotation-verdict-rules) to exactly one generated section per response, in order: Motivation, Methods, Results, Meanings.
2. Present the smallest wiki change set. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.
3. Apply only approved references and entries; show `git diff -- wiki`.
4. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next section; do not proceed.

### Gate 3 — Wiki to Raw to Wiki

1. Apply [Wiki-Rules](#wiki-rules) to exactly one wiki section per response, in template order.
2. Present the raw-and-wiki change set. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval; do not proceed.
3. Apply only approved changes under [Wiki-Rules](#wiki-rules); show `git diff -- raw wiki`.
4. **🔴 CHECKPOINT · 🛑 STOP** — Await explicit approval before the next wiki section; do not proceed.

### Gate 4 — Wiki to Tags

1. Reconcile vocab under [Tag-Rules](#tag-rules) after wiki approval.
2. Present exactly one smallest vocab change. **🔴 CHECKPOINT · 🛑 STOP** — Await separate explicit approval; do not proceed.
3. Apply only the approved change, verify under [Tag-Rules](#tag-rules), and show `git diff -- raw wiki vocab`.
