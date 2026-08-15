---
name: codes-rules
description: Establishes theory-document, variable-naming, comment, API-design, implementation-approval, and testing rules for C++ numerical tools under codes/. Use when adding, migrating, modifying, or reviewing code under codes/.
---

# Codes Rules

## Variable Naming Rules

- `name`: the meaning of the variable.
- `T`: the element type, such as `{int: I, double: F, complex<double>: C, template: T}`.
- `dD`: the variable rank, such as `{1D, 2D, 3D}`.
- `_indexes`: the index variables for each dimension, such as `{_x, _x_y, _x_y_z}`.
- Scalars use `name_T`, such as `{hbar_F, x_F, i_I}`.
- Vectors use `name_T1D_index`, such as `{x_F1D_x, n_I1D_n}`.
- Multidimensional arrays use `name_TdD_indexes`, such as `{H_F2D_x_n, psi_C2D_x_n}`.
- Index suffixes MUST match the physical meaning and ordering of the dimensions.

## Formal Theory Document Rules

- A formal theory document MUST share the code file's basename, such as `integration_gauss.hpp` and `integration_gauss.md`.
- The document MUST follow [notes-rules](../notes-rules/SKILL.md).
- Keep only the necessary formal theory and make equations primary.
- NEVER repeat implementation details from the code.

## Comment Rules

- All comments MUST be written in English.
- Keep comments concise and use a single line whenever possible.

### File Header Comments

- Every C++ file MUST contain `@file`, `@author`, `@date`, and `@brief`.
- Every `.hpp` file MUST place `#pragma once` after the file header comment.

```cpp
/**
 * @file    filename.hpp
 * @author  author
 * @date    YYYY-MM-DD
 * @brief   File purpose
 */

#pragma once
```

### Function Comments

- Every function declaration MUST be preceded by `@brief`, `@math`, and `@output`.
- `@math` MUST be a concise equation using Unicode Greek letters where applicable.
- All explanatory prose beyond `@brief` MUST use `@note`.
- If a function has no separate declaration, place the comment immediately before its definition.

```cpp
/**
 * @brief  Function purpose
 * @math   y = αx + β
 * @output Outputs or side effects
 */
```

## Workflow

1. Confirm and read the reference code and documentation completely; if either is missing or the scope is unclear, STOP and ask the user.
2. Draft the same-basename formal theory document under `notes-rules`; write each section ONLY after explicit approval.
3. Propose the API from the approved formal theory; write ONLY function declarations after explicit approval.
4. Present one function's computation outline at a time; implement ONLY that function after explicit approval.
5. Run the corresponding test after each implementation; after the module is complete, run `codes\run.bat` and the full CTest suite, and NEVER claim completion if tests fail.
