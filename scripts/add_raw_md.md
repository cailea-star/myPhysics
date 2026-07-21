---
[source_json]: [this_source_json]
[doi]: [this_doi]
[first_author]: [this_first_name], [this_family_name]
[affiliations]: [this_affiliations]
---

# [title_name]

### Abstract

[abstract]

<!-- claim-types: require [claim_type]: motivation [and] require [claim_type]: definition if this section has at least one definition coverage target -->
<!-- coverage: for each of the three most central tags identified from the abstract that is neither a core [tag_type]: method tag nor a current-work output [tag_type]: quantity tag, require >=1 valid definition quotation -->
<!-- quotation: no section-specific requirements -->
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

<!-- claim-types: require [claim_type]: definition [and] require [claim_type]: method [and] [claim_type]: comparison is optional [and] require [claim_type]: assumption if the paper states an explicit assumption -->
<!-- coverage: for each core [tag_type]: method tag, require >=1 valid definition quotation [and] require >=1 valid method-use quotation -->
<!-- quotation: every quotation requires at least one directly supported [tag_type]: method tag from vocab/tags.json; reject otherwise -->
### Methods

##### quotation-xx

...

<!-- claim-types: require [claim_type]: result [and] [claim_type]: definition is optional [and] [claim_type]: comparison is optional -->
<!-- coverage: for each [tag_type]: quantity tag calculated or measured as a current-work output and not defined earlier, require >=1 valid definition quotation [and] for each figure or table, require >=1 valid complete-sentence quotation using its calculated or measured [tag_type]: quantity tags; otherwise report a gap -->
<!-- quotation: every quotation requires at least one directly supported [tag_type]: method tag and one calculated or measured [tag_type]: quantity tag; reject otherwise -->
### Results

...

<!-- claim-types: require [claim_type]: comparison [and] require [claim_type]: innovation -->
<!-- coverage: may repeat Results only for central-significance claims [and] cover supported high-level judgments from the abstract, summary, or conclusion [and] reject minor results, figure-only details, and claims outside the main message -->
<!-- quotation: every quotation requires at least one directly supported [tag_type]: method tag from vocab/tags.json; reject otherwise -->
### Meanings

...


<!-- claim-types: no section-specific type required [and] every quotation must still match its explicit primary claim -->
<!-- coverage: prioritize citations explicitly linked to validated gaps in core method definitions, current-output quantity definitions, or core mechanisms [and] include actually reused external models, formulas, parameter sets, or data sources [and] reject generic background citations and citations with no gap or reused-dependency link -->
<!-- quotation: include only quotations dependent on an external cited reference [and] use tags [and] internal Fig., Table, Eq., or section references are not secondary citations -->
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
