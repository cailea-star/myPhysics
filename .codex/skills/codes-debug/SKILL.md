---
name: codes-debug
description: Debug C++ scientific and numerical code under codes/ against reference implementations. Use for theory approval, computational data-flow tracing, standard debug-printing utilities, upstream-to-downstream numerical matching, and validation after refactoring or migration.
---

# Codes Debug

## Print Rules

- Prints MUST follow: inputs → reference outputs → computed outputs → errors.
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

## Test Rules

- Tests MUST remain linear; NEVER use fixtures or wrappers; report files exceeding 120 lines.
- Tests MUST sequentially stream: standard inputs → function usage → labeled prints → acceptance asserts.
- Reference outputs MUST be analytic solutions or copied reference values.
- Use one-to-five single-line asserts with explicit tolerances; prints NEVER determine success.

## Debug Workflow

1. MUST obtain user approval that current-code and [Reference Code] theories are correct and equivalent.
2. MUST map the [Reference Code] [Computational Data Flow] upstream-to-downstream, e.g., grid -> potential -> wave function.
3. MUST build `tools_debug.hpp/tools_debug.cpp` upstream-to-downstream using the Print Rules above.
4. MUST verify upstream values against [Reference Code]; proceed downstream ONLY after they match.
5. MUST get approval before `$codes-debug` creates test files. Post-write, MUST audit code/theory, run `codes\run.bat` and full CTest; NEVER finish failing.
