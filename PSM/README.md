# Original PSM examples

`adgj` is the CPC archive for Y. Sun and K. Hara, Computer Physics
Communications 104 (1997) 245. It contains four independent Fortran programs.
The two parity labels refer to neutron and proton numbers, in that order.

| Package | Neutrons | Protons | Example N, Z | Final spin |
| --- | --- | --- | --- | --- |
| EE | even | even | 96, 70 | 30 |
| EO | even | odd | 96, 69 | 59/2 |
| OE | odd | even | 95, 70 | 59/2 |
| OO | odd | odd | 95, 69 | 25 |

Each package contains `XX_DATA` (input), `XX_PARA` (compile-time dimensions),
`PSM_XX.f` (complete source), and `XX_OUT1`, `XX_OUT2`, `XX_OUT3`
(published reference outputs). There are 24 files in total.

From the repository root, extract the archive with:

```powershell
python PSM/split_adgj.py
```

The script writes `PSM/split/{EE,EO,OE,OO}/`, preserving original bytes,
line endings, and card columns. It checks section boundaries and inventory,
and refuses to replace an extracted file whose contents have changed.

With Python and GNU Fortran on PATH, compile and run all four examples:

```powershell
python PSM/run_examples.py
```

To run only the even-even example:

```powershell
python PSM/run_examples.py EE
```

The runner uses `-O2 -std=legacy -fallow-argument-mismatch`. Lowercase `.f`
selects fixed-form Fortran without preprocessing. Each program runs in
`PSM/build/XX/`, containing its executable, copied input, `run.log`, and new
outputs. Repeated runs replace these generated results. The original outputs
in `split/` remain available for comparison. `build/` is ignored by Git.
The runner checks process exit; it does not certify
numerical agreement by itself.

Verified on 2026-09-16 with GNU Fortran 13.2.0 on Windows: all four examples
completed through their final requested spin without source changes.
All extracted lines retain 80 columns. All 12 output files have the same
line counts as the reference. Ignoring trailing whitespace, changed line
counts for OUT1/OUT2/OUT3 are EE: 4/0/1, EO: 2/1/0, OE: 2/2/0,
OO: 3/0/0. The observed changes affect the last printed numeric digit.
OUT3 agrees exactly at printed precision for EO/OE/OO; EE has one value
145.7824 instead of 145.7823 (absolute difference 0.0001).
