---
name: codes-rules
description: Establishes theory-document, variable-naming, comment, API-design, implementation-approval, and testing rules for C++ numerical tools under codes/. Use when adding, migrating, modifying, or reviewing code under codes/.
---

# Codes Rules

## Variable Naming Rules

- `name`: the actual mathematical or physical meaning of the variable; it MUST match the corresponding formula symbol, such as {$N_x$ → `Nx_I`, $x_i$ → `x_F1D_x`, $H_n(x_i)$ → `H_F2D_n_x`}.
- `T`: the element type, such as `{int: I, double: F, complex<double>: C, template: T}`.
- `dD`: the variable rank, such as `{1D, 2D, 3D}`.
- `_indexes`: the index variables for each dimension, such as `{_x, _x_y, _x_y_z}`.
- Scalars use `name_T`, such as `{hbar_F, x_F, i_I}`.
- Vectors use `name_T1D_index`, such as `{x_F1D_x, n_I1D_n}`.
- Multidimensional arrays use `name_TdD_indexes`, such as `{H_F2D_x_n, psi_C2D_x_n}`.
- Index suffixes MUST match the physical meaning and ordering of the dimensions.
- Before every code write, MUST list all added or renamed variables in a change table; write ONLY after explicit approval.

## Variable Construction Rules

- Numerical lists, vectors, matrices, and arrays SHOULD use Eigen types whenever practical.
- Every `Eigen::Tensor` MUST explicitly use `Eigen::ColMajor`.

## Code Formatting Rules

- Every declaration and every function definition signature MUST be written on a single line.
- Every simple statement MUST occupy exactly one line: NEVER split one statement across multiple lines and NEVER place multiple statements on the same line.
- Single-statement `if` or `for` blocks MAY remain inline: `if (...) {...}`, `for (...) {...}`.
- Developer-facing numerical preconditions MUST use single-line `assert`; recoverable runtime failures MUST use explicit handling.

## Comment Rules

- All comments MUST be written in English.
- Keep comments concise; trailing comments SHOULD use formulas, such as `Scale: x → x / b`.
- Standalone comment lines MAY separate code blocks when necessary.

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

- `@brief`, `@math`, and `@output` MUST immediately precede each declaration or declarationless definition.
- `@brief` MUST name any numerical algorithm that determines behavior or accuracy.
- `@math` MUST be a concise equation using Unicode Greek letters where applicable.
- All explanatory prose beyond `@brief` MUST use `@note`.
- `@note` MUST state any required algorithm dependency, assumption, or limitation.

```cpp
/**
 * @brief  Integrate sampled data using Simpson quadrature.
 * @math   I ≈ Q_S[f]
 * @output Approximate integral.
 * @note   Requires an odd, panelwise-uniform grid.
 */
```

## Testing Rules

- Every test MUST contain, in order: function usage, labeled input/output, acceptance asserts.
- Every test MUST print the labeled reference result first—an analytic solution or copied reference value—then print the computed result.
- Vectors and tables MUST print first five, ellipsis, last five; print all when $N \leq 10$.
- Matrices MUST print four $5\times5$ corner blocks joined by ellipses; print all when each dimension is at most 10.
- Tests MUST remain linear; NEVER add fixtures or wrappers; MUST report files exceeding 100 lines.
- Prints NEVER determine success; use one-to-five single-line asserts on key values with explicit tolerances.

## Workflow

1. Confirm and read the reference code and documentation completely; if either is missing or the scope is unclear, STOP and ask the user.
2. Draft same-basename theory under `notes-rules`; keep equations primary; NEVER repeat implementation; write ONLY after approval.
3. Propose the API from the approved formal theory; write ONLY function declarations after explicit approval.
4. Present one function's computation outline at a time; implement ONLY that function after explicit approval.
5. Test each implementation; finally audit code/theory, run `codes\run.bat` and full CTest; NEVER complete if failing.
