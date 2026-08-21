---
name: academic-writing
description: Structure academic and technical writing through approval-gated goals, keyword Wiki maintenance, section and paragraph functions, keyword chains, and verified Raw quotations. Use for introductions, research proposals, literature-grounded outlines, or evidence-backed academic writing.
---

# Academic Writing

## Gated Workflow

### Gate 1 — Keywords and Writing Goal

Gate 1 MUST present `title`, `method`, `topics`, and `goal`; NEVER write before approval.

- `title` MUST identify the article subject; `goal` MUST define the writing objective.
- `method` MUST identify the research method; `topics` MUST contain ONLY 2–5 core keywords.
- After approval, `papers-lint` MUST process every `topics` entry; Gate 2 MUST NEVER start until all pass.
- Gate 1 MUST fill ONLY article-level fields; section entries MUST remain blank until Gate 2.

### Gate 2 — Section Chain and Section Goals

Gate 2 MUST present every section's `keyword` and `section-goal`; NEVER write before approval.

- Each `keyword` MUST come ONLY from `topics`; sections MUST NEVER repeat keywords.
- Each `section-goal` MUST define one `keyword` logic and exit pressure; together MUST chain toward `goal`.
- Concepts MUST follow ONLY “introduction → definition → emphasis”; methods MUST follow ONLY “introduction → definition → advantage.”
- Gate 2 MUST fill ONLY section fields; paragraph fields MUST remain blank until Gate 3.

### Gate 3 — Paragraph Chain and Paragraph Functions

Gate 3 MUST present every paragraph's title, `keyword`, and `function`; NEVER write before approval.

- Each paragraph MUST use ONLY one `keyword`; `function` MUST follow `[claim-type] — ...`.
- Each paragraph MUST use ONLY one `claim-type`; `function` MUST state how it advances `section-goal`.
- Paragraphs MUST chain toward section `keyword` and `section-goal`; otherwise MUST return to Gate 2.
- Gate 3 MUST fill ONLY paragraph fields; `chain` MUST remain blank until Gate 4.

### Gate 4 — Paragraph Keyword Chains

Gate 4 MUST present every paragraph's `chain`; NEVER write before approval.

- Each `chain` MUST center on its paragraph's `keyword` and `claim-type`.
- Every relation MUST follow “A → B (claim-type)”; its type MUST exist in `types.json`.
- Missing concepts, premature definitions, overstrong or duplicated relations, or misplaced emphasis MUST return to Gate 3.
- Gate 4 MUST fill ONLY `chain`; `Raw` MUST remain blank until Gate 5.

### Gate 5 — Raw Evidence

Gate 5 MUST retrieve `Raw` for every relation via relevant Wikis; NEVER write before approval.

- Before retrieval, MUST reread the current paragraph's `keyword`, `function`, and `chain`.
- Use Wikis ONLY to locate papers; every claim MUST trace to its `Raw` quote.
- Every relation MUST contain at least two strong `Raw` entries from different core references.
- NEVER write insufficient evidence; adjust relations in Gate 4 or supplement Wikis in Gate 1.

## Template

````md
# [title]

> method: ...
> topics: [keyword1], [keyword2], ...
> goal: ...

### [section1]

> keyword: [keyword1]
> section-goal: ...

##### [paragraph1-1]

> keyword: [keyword1]
> function: [claim-type] — ...

chain:
- [keyword1] → [keywordtmp1] (claim-type)

```tags
[source]: [raw-stem]: [source].
[tags]: [keyword1], [keywordtmp1]
```

```quote
...
```

```tags
[source]: [raw-stem]: [source].
[tags]: [keyword1], [keywordtmp1]
```

```quote
...
```

```tags
[source]: [raw-stem]: [source].
[tags]: [keyword1], [keywordtmp1]
```

```quote
...
```

- [keywordtmp1] → [keywordtmp2] (claim-type)

```tags
[source]: [raw-stem]: [source].
[tags]: [keywordtmp1], [keywordtmp2]
```

```quote
...
```

```tags
[source]: [raw-stem]: [source].
[tags]: [keywordtmp1], [keywordtmp2]
```

```quote
...
```

...

##### [paragraph1-2]

...
````
