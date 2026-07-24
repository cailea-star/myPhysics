---
[tags]: [this_tags_name]
[tag_types]: [this_tag_types]
---

# [this_tag_name]

<!-- claim-types: ([claim_type]: definition or [claim_type]: background) -->
<!-- coverage: (from Motivation and Results, select the smallest sufficient root-topic quotation set with [claim_type]: definition or [claim_type]: background; prefer [claim_type]: definition) and (origin with [claim_type]: background is optional) -->
### Definition & Origin

```references
  - source1: Abstract.
  - ...
```

[What the topic is, where it appears, and how the concept originated.]

<!-- claim-types: ([claim_type]: definition or [claim_type]: background or [claim_type]: comparison) -->
<!-- coverage: (from Motivation and Results, select the smallest sufficient root-topic quotation set with [claim_type]: definition or [claim_type]: background or [claim_type]: comparison that describes composition or mechanism details) and (exclude quotations used in Definition & Origin) -->
### Composition & Mechanism

##### Composition & Mechanism 1

```references
  - source_xx: ...
  - ...
```

- **description**: This is a physical picture.

<!-- Optional: add a directly supported math block for this item. -->

##### Composition & Mechanism 2

...

<!-- claim-types: ([claim_type]: definition or [claim_type]: background or [claim_type]: result) -->
<!-- coverage: (from Motivation and Results, select the smallest sufficient quotation set for root-topic-related derived [tag_type]: quantity or [tag_type]: property tags with [claim_type]: definition or [claim_type]: background or [claim_type]: result) and (exclude quotations used in Definition & Origin or Composition & Mechanism) -->
### Quantities & Properties

```tags
[tags]: ...
```

##### Quantity or Property 1

```references
  - source_xx: ...
  - ...
```

```math
\textrm{This is quantity-expression-1.}
```
- **Definition**: variable description.
- **description**: This is a physical picture of quantity-expression-1.


##### Quantity or Property 2

...

<!-- claim-types: ([claim_type]: method) and ([claim_type]: result) -->
<!-- coverage: (from Methods, select root-topic-related core [tag_type]: method quotations with [claim_type]: method) and (from Results, select their output quotations with [claim_type]: result) -->
### Related Methods

##### Method 1

```tags
[tags]: ...
```

```references
  - source_xx: ...
  - ...
```

- **description**: How this core method obtains or analyzes the root topic.
- **Outputs**: Outputs provided by Results quotations and its precision(must have, at least claim it is none).
- **link**: [method](wiki_path)

<!-- Optional: remove the link when no method Wiki exists; add a directly supported math block when available. -->

##### Method 2

...

<!-- claim-types: ([claim_type]: motivation is required for each Study) and ([claim_type]: comparison is required for each Study) and ([claim_type]: innovation or [claim_type]: result is required for each Study) and ([claim_type]: result is required for each Study) -->
<!-- coverage: (each Study describes exactly one raw paper) and (for each Study, its Motivation provides [claim_type]: motivation) and (for each Study, its Meanings provides [claim_type]: comparison) and (for each Study, its Meanings provides [claim_type]: innovation or [claim_type]: result) and (for each Study, its Meanings provides [claim_type]: result) -->
### Previous Studies

##### Study 1

```tags
[tags]: ...
```

```references
  - source_xx: ...
```

- **Motivation**: This is the motivation of the study 1.
- **Innovation or Result**: This is the innovation or result of the study 1.
- **Quantities**: This is the quantities of the study 1.
- **Precision**: exact precision or none.
- **Meaning**: This is the meaning of the study 1.

##### Study 2

...

<!-- claim-types: (none: Candidate Papers do not support Wiki claims) -->
<!-- coverage: (from Secondary, identify every unique Candidate Paper) and (for each Candidate Paper, include tags, title, DOI, journal, year, reason) and (if no Candidate Paper qualifies, report not-applicable) -->
### Next Papers

##### Candidate Paper 1

```tags
[tags]: ...
[doi]: ...
[title]: ...
```

- **Citation**: journal, year.
- **Reason**: why this paper should be ingested next.

##### Candidate Paper 2

...
