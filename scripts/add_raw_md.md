---
[source_json]: [this_source_json]
[doi]: [this_doi]
[first_author]: [this_first_name], [this_family_name]
[affiliations]: [this_affiliations]
---

# [title_name]

### Abstract

[abstract]

<!-- claim-types: require definition [and] require motivation -->
<!-- coverage: for each of the three most central tags identified from the abstract, require >=1 valid definition quotation -->
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

<!-- claim-types: require definition [and] require comparison [and] require assumption if the paper states an explicit assumption -->
<!-- coverage: for each core method tag, require >=1 valid definition quotation [and] for each core method tag, require >=1 valid comparison quotation -->
<!-- quotation: no section-specific requirements -->
### Methods

##### quotation-xx

...

<!-- claim-types: require result [and] definition is optional [and] comparison is optional -->
<!-- coverage: for each calculated or measured quantity tag not defined earlier, require >=1 valid definition quotation [and] for each figure or table, require >=1 valid complete-sentence quotation using its calculated or measured quantity tags; otherwise report a gap -->
<!-- quotation: every quotation requires at least one directly supported calculation-method tag and one calculated or measured quantity tag; reject otherwise -->
### Results

...

<!-- claim-types: require comparison [and] require innovation -->
<!-- coverage: may repeat Results only for central-significance claims [and] cover supported high-level judgments from the abstract, summary, or conclusion [and] reject minor results, figure-only details, and claims outside the main message -->
<!-- quotation: every quotation requires at least one directly supported method tag from vocab/tags.json; reject otherwise -->
### Meanings

...


<!-- claim-types: no section-specific type required [and] every quotation must still match its explicit primary claim -->
<!-- coverage: include only external methods/formulas, parameter sets/models/data sources, or key physical conclusions/mechanism judgments supporting core tags [and] reject generic survey/background citations and citations that supply no reused method, formula, model/data input, or core interpretation -->
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
