---
[source_json]: [this_source_json]
[doi]: [this_doi]
[first_author]: [this_first_name], [this_family_name]
[affiliations]: [this_affiliations]
---

# [title_name]

### Abstract

[abstract]

<!-- claim-types: ([claim_type]: motivation is required) and ([claim_type]: definition or [claim_type]: background is required for core-tag coverage; prefer definition) -->
<!-- coverage: (identify the three most central tags from the Abstract, excluding [tag_type]: method and current-work output [tag_type]: quantity) and (for each identified tag, assign >=1 valid quotation to [claim_type]: definition or [claim_type]: background; otherwise report a gap) -->
<!-- quotation: (none) -->
### Motivation

##### quotation-01

```tags
[claim_type]: ...
[tags]: ...
[source]: section 1.
```

```quote
This is a quotation from the original paper, which is tagged with the above tags.
```

##### quotation-02

...

<!-- claim-types: ([claim_type]: definition or [claim_type]: background is required for core-method coverage; prefer definition) and ([claim_type]: method is required) and ([claim_type]: comparison is optional) and ([claim_type]: assumption is required if the paper states an explicit assumption) -->
<!-- coverage: (identify the core [tag_type]: method tags) and (for each identified method, assign >=1 valid quotation to [claim_type]: definition or [claim_type]: background; otherwise report a gap) and (for each identified method, require >=1 valid method-use quotation under its true [claim_type]; otherwise report a gap) -->
<!-- quotation: (each quotation requires >=1 directly supported [tag_type]: method tag; reject otherwise) -->
### Methods

##### quotation-xx

...

<!-- claim-types: ([claim_type]: result is required) and ([claim_type]: definition or [claim_type]: background is required if output-quantity coverage has at least one target; prefer definition) and ([claim_type]: comparison is optional) -->
<!-- coverage: (identify current-work output [tag_type]: quantity tags not defined earlier) and (for each identified quantity, assign >=1 valid quotation to [claim_type]: definition or [claim_type]: background; otherwise report a gap) and (for each figure or table, require >=1 valid complete-sentence quotation using its current-work output quantity tags; otherwise report a gap) -->
<!-- quotation: (each quotation requires >=1 directly supported current-work output [tag_type]: quantity tag; reject otherwise) and (each quotation not assigned to [claim_type]: definition or [claim_type]: background requires >=1 directly supported [tag_type]: method tag; reject otherwise) -->
### Results

...

<!-- claim-types: ([claim_type]: comparison is required) and ([claim_type]: innovation is required) -->
<!-- coverage: (identify supported high-level judgments from the abstract, summary, or conclusion) and (for each identified judgment, require >=1 valid quotation) and (repeat Results only for central-significance claims) and (reject minor results, figure-only details, and claims outside the main message) -->
<!-- quotation: (each quotation requires >=1 directly supported [tag_type]: method tag; reject otherwise) -->
### Meanings

...


<!-- claim-types: (none: no section-specific [claim_type] is required) and (each quotation must match its explicit primary claim) -->
<!-- coverage: (identify eligible external citations linked to validated gaps in core method definitions, current-output quantity definitions, or core mechanisms, or to reused external models, formulas, parameter sets, or data sources) and (for each identified citation, require >=1 valid quotation; otherwise report a gap) and (reject generic background citations and citations without a gap or reused-dependency link) -->
<!-- quotation: (include only quotations dependent on an explicit external citation) and (each quotation requires directly supported tags) and (internal Fig., Table, Eq., or section references do not qualify) -->
### Secondary Citations

##### quotation-xx

```tags
[claim_type]: ...
[tags]: ...
[source]: ...
[ref]: [ab] reference index infromation from the original paper
[doi]: doi number of the reference or unknown.
```

```quote
This is a quotation with secondary citations from the original paper[ab], which is tagged with the above tags.
```

<!-- gaps: record only validated gap verdicts after full-text review; use one bullet per independently counted template item; limit [reason] to 10 words; delete the placeholder when there are no gaps; do not record fix verdicts -->
### Gaps

- `[section]: [section name] | [item]: [claim-types|coverage|quotation] | [target]: [exact requirement or expanded target] | [reason]: [why no valid evidence exists, up to 10 words]`
