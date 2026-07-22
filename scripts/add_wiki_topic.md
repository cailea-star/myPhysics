---
[tags]: [this_tags_set]
[tag_types]: [this_tag_types]
---

# [tag_name]

<!-- claim-types: ([claim_type]: definition or [claim_type]: background is required for the root topic; prefer definition) and ([claim_type]: background or [claim_type]: result is required for each applicability condition under its true claim type) -->
<!-- coverage: (identify the root canonical tag and every explicit applicability condition) and (assign >=1 supported Wiki claim defining the root topic to [claim_type]: definition or [claim_type]: background; otherwise report missing) and (for each identified applicability condition, assign it to [claim_type]: background or [claim_type]: result according to its explicit primary claim; otherwise report missing) and (origin is optional and requires direct primary evidence) -->
### Definition & Origin

```references
  - source1: Abstract.
  - ...
```

[What the topic is, where it appears, and how the concept originated.]

<!-- claim-types: ([claim_type]: definition or [claim_type]: background or [claim_type]: result is required for each core component, state, condition, or mechanism under its true claim type) and ([claim_type]: definition or [claim_type]: background is required for each core equation; prefer definition) -->
<!-- coverage: (identify every core component or state, every existence condition, every core or competing mechanism, and every core equation) and (for each identified component or state, assign it to [claim_type]: definition, [claim_type]: background, or [claim_type]: result; otherwise report missing) and (for each identified condition, assign it to [claim_type]: background or [claim_type]: result; otherwise report missing) and (for each identified mechanism, assign it to [claim_type]: definition, [claim_type]: background, or [claim_type]: result; otherwise report missing) and (for each identified core equation, assign it to [claim_type]: definition or [claim_type]: background with its physical meaning; otherwise report missing) and (if no competing mechanism qualifies, report not-applicable) -->
### Composition & Mechanism

```tags
[tags]: topic, composition, mechanism.
```

##### Composition

```references
  - source_xx: ...
  - ...
```

- **Composition**: Components, participants, or states.
- **Conditions**: Conditions under which the system or phenomenon exists.

##### Mechanism

```references
  - source_xx: ...
  - ...
```

- **Competing Mechanisms**: Alternative mechanisms or none.

---

```math
\textrm{This is core-equation-1.}
```

- **description**: This is a physical picture of core-equation-1.

---

```math
\textrm{This is core-equation-2.}
```

- **description**: This is a physical picture of core-equation-2.

...

<!-- claim-types: ([claim_type]: definition or [claim_type]: background is required for the meaning and variables of each core [tag_type]: quantity or [tag_type]: property; prefer definition) and ([claim_type]: definition or [claim_type]: background is required for the physical interpretation of each core [tag_type]: quantity or [tag_type]: property) and ([claim_type]: background or [claim_type]: result is required for the applicable conditions of each core [tag_type]: quantity or [tag_type]: property under their true claim type) and ([claim_type]: definition or [claim_type]: background is required for each associated equation; prefer definition) -->
<!-- coverage: (identify every core [tag_type]: quantity or [tag_type]: property and every associated equation) and (for each identified quantity or property, assign its meaning and variables to [claim_type]: definition or [claim_type]: background; otherwise report missing) and (for each identified quantity or property, assign its physical interpretation to [claim_type]: definition or [claim_type]: background; otherwise report missing) and (for each identified quantity or property, assign its applicable conditions to [claim_type]: background or [claim_type]: result; otherwise report missing) and (for each identified equation, assign it to [claim_type]: definition or [claim_type]: background with its meaning; otherwise report missing) and (if no core quantity or property qualifies, report not-applicable) -->
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

<!-- claim-types: ([claim_type]: method is required for each related [tag_type]: method) and ([claim_type]: result or [claim_type]: comparison is required for the precision of each related [tag_type]: method under the true claim type) and ([claim_type]: background or [claim_type]: result is required for the applicability of each related [tag_type]: method under the true claim type) and ([claim_type]: definition or [claim_type]: background is required for each core equation; prefer definition) -->
<!-- coverage: (identify every related [tag_type]: method centrally used to obtain or analyze the root topic and every associated core equation) and (for each identified method, assign its use, inputs, and outputs to [claim_type]: method; otherwise report missing) and (for each identified method, assign its precision to [claim_type]: result or [claim_type]: comparison; otherwise report missing) and (for each identified method, assign its applicability to [claim_type]: background or [claim_type]: result; otherwise report missing) and (for each identified core equation, assign it to [claim_type]: definition or [claim_type]: background with its meaning; otherwise report missing) and (if no related method qualifies, report not-applicable) -->
### Related Methods

##### Method 1

```tags
[tags]: ...
```

```references
  - source_xx: ...
  - ...
```

- **Inputs**:
- **Outputs**:
- **Precision**:
- **Applicability**:


```math
\textrm{ This is core-equation of the method1.}
```

- **description**: This is a physical picture of equation-1.

- **link**: [wiki_file_name](wiki_file_path)

##### Method 2

...

<!-- claim-types: ([claim_type]: motivation or [claim_type]: background is required for each study) and ([claim_type]: result or [claim_type]: comparison is required for each study) and ([claim_type]: innovation is required for each study) -->
<!-- coverage: (identify every primary study that centrally investigates, measures, computes, or compares the root topic) and (for each identified study, identify its core method and quantity or property tags) and (for each identified study, assign its purpose to [claim_type]: motivation or [claim_type]: background; otherwise report missing) and (for each identified study, assign its quantities, precision, and main meaning to [claim_type]: result or [claim_type]: comparison; otherwise report missing) and (for each identified study, assign its explicitly supported innovation to [claim_type]: innovation; otherwise report missing) and (if no primary study qualifies, report not-applicable) -->
### Previous Studies

##### Study Type 1

```tags
[tags]: ...
```

```references
  - source_xx: ...
  - ...
```

- **Motivation**: This is the motivation of the study type 1.
- **Innovation**: This is the innovation of the study type 1.
- **Quantities**: This is the quantities of the study type 1.
- **Precision**: exact precision or none.
- **Meaning**: This is the meaning of the study type 1.

##### Study Type 2

...

<!-- claim-types: ([claim_type]: motivation or [claim_type]: background or [claim_type]: comparison or [claim_type]: result is required for each limitation or open question under its true claim type) -->
<!-- coverage: (identify every explicit current-precision limitation, unmeasured or unexplained quantity, competing explanation, methodological gap, and open question) and (for each identified limitation or open question, retain its supported quantity, condition, phenomenon, or method tags) and (for each identified limitation or open question, assign it to [claim_type]: motivation, [claim_type]: background, [claim_type]: comparison, or [claim_type]: result according to its explicit primary claim; otherwise report missing) and (if no limitation or open question qualifies, report not-applicable) -->
### Limitations & Open Questions

```references
  - source_xx: ...
  - ...
```

- **Current Precision**: Current experimental, numerical, or theoretical precision; use `not quantified` when unavailable.
- **Unmeasured or Unexplained Quantities**: Quantities or behavior not yet measured, calculated, or explained.
- **Competing Explanations**: Unresolved interpretations or physical mechanisms.
- **Methodological Gaps**: Missing methods or limitations shared by current approaches.
- **Open Questions**: Questions not resolved by current studies.
