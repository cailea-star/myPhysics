---
name: codes-rules
description: Establishes theory-document, variable-naming, comment, API-design, implementation-approval, and testing rules for C++ numerical tools under codes/. Use when adding, migrating, modifying, or reviewing code under codes/.
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

- Numerical lists, vectors, matrices, and arrays SHOULD use Eigen types whenever practical.
- Every `Eigen::Tensor` MUST explicitly use `Eigen::ColMajor`.
- Matrix products MUST use `.noalias()` except when initializing a new variable; aliased products MUST use a temporary.
- Large numerical buffers MUST be allocated outside loops and reused, such as {Vector, Matrix, Tensor}.

## Code Formatting Rules

- Within each class access section, declare data members first; define constructors inline before member functions.
- Constructor parameters matching members SHOULD append `_`; initializer lists MUST follow declaration order on one separate line.
- Declarations, function signatures, and simple statements MUST each occupy one line and NEVER share lines.
- Single-statement `if` or `for` blocks MAY remain inline but MUST retain braces: `if (...) {...}`, `for (...) {...}`.
- Developer-facing numerical preconditions MUST use single-line `assert`; recoverable runtime failures MUST use explicit handling.

## Comment Rules

- All comments MUST be written in English and be as concise as possible.
- Comments SHOULD use formulas or transformations instead of prose whenever possible, such as `x → x / b`.
- Comments MUST add information; NEVER merely restate symbols, such as `{x_i}`.
- Every logical computation block MUST begin with one standalone comment line.

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

### Structure

- Tests MUST remain linear; NEVER use fixtures or wrappers; report files exceeding 120 lines.
- Tests MUST sequentially stream: standard inputs → function usage → labeled prints → acceptance asserts.
- Prints MUST follow: inputs → reference outputs → computed outputs → errors.
- Reference outputs MUST be analytic solutions or copied reference values.
- Use one-to-five single-line asserts with explicit tolerances; prints NEVER determine success.

### Printing

- Vectors MUST print horizontally: first five, ellipsis, last five.
- Tables MUST print labeled columns vertically: first five rows, standalone ellipsis, last five rows.
- Matrices MUST print four $5\times5$ corner blocks joined by ellipses.
- Print all entries when every relevant dimension is at most 10.

```text
[Input] x = 0, 1, 2, 3, 4, ..., 96, 97, 98, 99, 100
[Reference] I = 3.141593
[Computed] I = 3.141592
[Computed] f(x):
x    f
0    0.000000
1    0.099833
2    0.198669
3    0.295520
4    0.389418
...  ...
96  -0.174327
97  -0.271761
98  -0.366479
99  -0.457536
100 -0.544021
[Error] |I - I_ref| = 1.000000e-6
```

## Workflow

1. Confirm and read the reference code, theory, target files, and direct callers completely; if anything is missing or the scope is unclear, STOP and ask the user.
2. Draft same-basename theory under `notes-rules`; keep equations primary; NEVER repeat implementation; write ONLY after approval.
3. Propose the API from the approved formal theory; write ONLY approved API declarations, while constructors MUST include their approved inline definitions.
4. Present one function's computation outline at a time; implement ONLY that function after explicit approval.
5. After each write, review the diff and run focused tests; finally audit code/theory, run `codes\run.bat` and full CTest; NEVER complete if failing.
