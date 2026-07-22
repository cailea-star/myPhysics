---
[tags]: [this_tags_set]
[tag_types]: [this_tag_types]
---

# [tag_name]

<!-- claim-types: ([claim_type]: definition or [claim_type]: background is required for the root method; prefer definition) and ([claim_type]: motivation or [claim_type]: background is required for the problem it solves) -->
<!-- coverage: (identify the root [tag_type]: method tag) and (assign >=1 supported Wiki claim defining the method to [claim_type]: definition or [claim_type]: background; otherwise report missing) and (assign >=1 supported Wiki claim explaining the problem it solves to [claim_type]: motivation or [claim_type]: background; otherwise report missing) and (origin is optional and requires direct primary evidence) -->
### Definition & Origin

```references
  - source1: Abstract.
  - ...
```

[What the method is, what problem it solves, and where it originated.]

<!-- claim-types: ([claim_type]: method is required) and ([claim_type]: assumption is required for each explicit assumption, approximation, or truncation) and ([claim_type]: definition or [claim_type]: background is required for each core equation; prefer definition) -->
<!-- coverage: (identify the operational framework, every parameter-dependence target, every explicit assumption, approximation, or truncation, and every core equation of the root method) and (assign >=1 supported Wiki claim describing the framework to [claim_type]: method; otherwise report missing) and (for each identified parameter-dependence target, assign it to [claim_type]: method or [claim_type]: assumption under its true claim type; otherwise report missing) and (for each identified assumption, approximation, or truncation, assign it to [claim_type]: assumption; otherwise report missing) and (for each identified core equation, assign it to [claim_type]: definition or [claim_type]: background with its meaning; otherwise report missing) -->
### Framework

```tags
[tags]: method, conditions, approximations, and assumptions
```

```references
  - source_xx: ...
  - ...
```

##### Assumptions and Approximations

- **Assumptions**:
- **Truncations or approximations**:
- **Parameter dependence**:

##### Core Equations

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

<!-- claim-types: ([claim_type]: definition or [claim_type]: background is required for each core input or output tag; prefer definition) and ([claim_type]: method is required for each input-output relation) -->
<!-- coverage: (identify every canonical tag explicitly used as an input, parameter, initial condition, or output of the root method and every input-output relation) and (for each identified input or output tag, assign it to [claim_type]: definition or [claim_type]: background; otherwise report missing) and (for each identified input-output relation, assign it to [claim_type]: method; otherwise report missing) -->
### Inputs & Outputs

```tags
[tags]: ...
```

##### Inputs

```references
  - source_xx: ...
  - ...
```

core inputs, parameters, and initial conditions.

##### Outputs

```references
  - source_xx: ...
  - ...
```

basic outputs, derived quantities, and final results.


<!-- claim-types: ([claim_type]: motivation is required for each correction problem) and ([claim_type]: method is required for each correction method) and ([claim_type]: result or [claim_type]: innovation is required for each correction effect; use innovation only when explicit) and ([claim_type]: definition or [claim_type]: background is required for each correction equation; prefer definition) -->
<!-- coverage: (identify every [tag_type]: method explicitly used to correct the root method and every associated correction equation) and (for each identified correction method, assign its problem to [claim_type]: motivation; otherwise report missing) and (for each identified correction method, assign its procedure to [claim_type]: method; otherwise report missing) and (for each identified correction method, assign its effect to [claim_type]: result or [claim_type]: innovation; otherwise report missing) and (for each identified correction equation, assign it to [claim_type]: definition or [claim_type]: background with its meaning; otherwise report missing) and (if no correction method qualifies, report not-applicable) -->
### Correction Algorithm

##### Correction Algorithm 1

```tags
[tags]: ...
```

```references
  - source_xx: ...
  - ...
```

- **Problem**: Limitations or deficiencies of the main framework.
- **Solution**: How the correction algorithm addresses the problem.
- **Effect**: How the correction algorithm improves the main framework.

---

```math
\textrm{ This is correction-equation-1 of the solution1.}
```

- **description**: This is a physical picture of correction-equation-1.

---

```math
\textrm{This is correction-equation-2 of the solution1.}
```

- **description**: This is a physical picture of correction-equation-2.

---

...


##### Correction Algorithm 2

...

<!-- claim-types: ([claim_type]: motivation or [claim_type]: background is required for each study) and ([claim_type]: result or [claim_type]: comparison is required for each study) and ([claim_type]: innovation is required for each study) -->
<!-- coverage: (identify every primary study that centrally uses, tests, or compares the root method) and (for each identified study, identify its core method and output quantity or property tags) and (for each identified study, assign its purpose to [claim_type]: motivation or [claim_type]: background; otherwise report missing) and (for each identified study, assign its main finding to [claim_type]: result or [claim_type]: comparison; otherwise report missing) and (for each identified study, assign its explicitly supported innovation to [claim_type]: innovation; otherwise report missing) and (if no primary study qualifies, report not-applicable) -->
### Previous Studies

#### Study Type 1

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

#### Study Type 2

...

<!-- claim-types: ([claim_type]: motivation or [claim_type]: background or [claim_type]: comparison or [claim_type]: result is required for each limitation under its true claim type) -->
<!-- coverage: (identify every explicit precision, missing-physics, formal, computational, or applicability limitation of the root method) and (for each identified limitation, retain its supported quantity, condition, or method tags) and (for each identified limitation, assign it to [claim_type]: motivation, [claim_type]: background, [claim_type]: comparison, or [claim_type]: result according to its explicit primary claim; otherwise report missing) and (if no limitation qualifies, report not-applicable) -->
### Limitations

```references
  - source_xx: ...
  - ...
```

- **Current Precision**: Current numerical precision, experimental error, or theoretical uncertainty; use `not quantified` when unavailable.
- **Missing Quantities or Physics**: Quantities, phenomena, or problems the method cannot yet describe.
- **Formal Limitations**: Limitations caused by assumptions, truncations, or parameter dependence.
- **Computational Limitations**: Computational cost, numerical stability, or scaling limits.
- **Applicability Limits**: Systems, conditions, or regimes where the method becomes unreliable.
