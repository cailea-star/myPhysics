---
name: codes-rules
description: Establishes variable-naming, construction, formatting, comment, API-design, and implementation-approval rules for C++ numerical tools under codes/. Use when adding, migrating, modifying, or reviewing functional code under codes/.
---

# Codes Rules

## Variable Naming Rules

- `name` MUST encode the mathematical or physical meaning and match its formula symbol, such as {$N_x$ → `Nx_I`, $x_i$ → `x_F1D_x`, $H_n(x_i)$ → `H_F2D_n_x`}.
- `T` MUST encode the element type using `{int: I, double: F, complex<double>: C, template: T}`.
- `dD` MUST encode rank when present: scalars use `name_T`, such as `{hbar_F, x_F, i_I}`; vectors use `name_T1D_index`, such as `{x_F1D_x, n_I1D_n}`; multidimensional arrays use `name_TdD_indexes`, such as `{H_F2D_x_n, psi_C2D_x_n}`.
- `_indexes` MUST list each dimension's physical index in array order, such as `{_x, _x_y, _x_y_z}`.
- Function names MUST use `snake_case`; mathematical symbols MAY preserve canonical capitalization, such as `{calc_energy, build_tables, read_Gz}`.
- Before every code write, MUST list all added or renamed variables—or `None`—and write ONLY after explicit approval.

## Variable Construction Rules

- Numerical containers SHOULD use Eigen when practical; every `Eigen::Tensor` MUST explicitly use `Eigen::ColMajor`.
- Matrix products MUST use `.noalias()` except when initializing a new variable; aliased products MUST use a temporary.
- Large numerical buffers MUST be allocated outside loops and reused, such as {Vector, Matrix, Tensor}.
- Numerical computation blocks MUST NOT use runtime conditional branches; exceptions require explicit approval.

## Code Formatting Rules

- Within each class access section, declare data members first; define constructors inline before member functions.
- Constructor parameters matching members SHOULD append `_`; initializer lists MUST follow declaration order on one separate line.
- Declarations, function signatures, and simple statements MUST each occupy one line and NEVER share lines.
- Single-statement `if` or `for` blocks MAY remain inline but MUST retain braces: `if (...) {...}`, `for (...) {...}`.
- Floating-point-to-integer conversion MUST use `std::floor` or `std::ceil`.
- Developer-facing numerical preconditions MUST use single-line `assert`; recoverable runtime failures MUST use explicit handling.

## Comment Rules

- All comments MUST be concise English.
- Comments SHOULD prefer formulas or transformations over prose, such as `x → x / b`.
- Comments MUST add information and NEVER merely restate symbols, such as `{x_i}`.
- Standalone comments MUST begin and delimit every logical computation block.
- End-of-line comments MUST state variable meaning, units, normalization, or domain.

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

## Workflow

1. MUST validate user-approved [Reference Code] and theory dependencies; reproduce ONLY the smallest independently verifiable module.
2. MUST draft same-basename theory via `$notes-rules`; prioritize equations; NEVER duplicate implementation; write ONLY after approval.
3. MUST derive APIs ONLY from approved theory; write ONLY approved declarations and inline constructor definitions.
4. Present one function's computation outline; implement ONLY that function after explicit approval.
5. MUST get approval before `$codes-debug` creates test files. Post-write, MUST audit code/theory, run `codes\run.bat` and full CTest; NEVER finish failing.
