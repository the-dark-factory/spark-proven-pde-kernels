# SPARK-proven 1-D PDE kernels: linear advection and inviscid Burgers

Lanyon AI publishes numerical schemes alongside Lean proofs of their properties. The Lean proofs are
about real numbers, but the C code runs in IEEE-754 binary64. This repository ports two of those 1-D
schemes to SPARK (Ada) and proves properties of the **floating-point code itself** with GNATprove. For
each Lean theorem, the question is which of three forms survives in floats:

- **exact**: the same statement holds bit for bit;
- **bounded**: it holds up to a proved rounding error ε;
- **not provable as stated**: it is false in floats, and we explain why.

The answers, theorem by theorem, are in [RESULTS.md](RESULTS.md). They include the limits.

- `src/advection_1d/`: linear advection f_t + a f_x = 0. Flux, Lax-Friedrichs wave split,
  fluctuations, minmod reconstruction, one upwind cell update over a periodic grid, an exact
  integer variant, and the ghost (proof-only) rounding lemmas behind the two bounded theorems.
- `src/burgers_1d/`: inviscid Burgers f_t + (f²/2)_x = 0, in the same shape. Here the wave split
  divides by the state itself, and the sonic point (both states 0) is defined rather than NaN.

Every unit is a package specification whose subprograms are expression functions, so there are no
`.adb` bodies. Each file is byte-identical to the version our admission process re-proved and
recorded; the digests are in `SHA256SUMS`, and the ledger numbers are in RESULTS.md.

## How to prove it

You need **gnatprove FSF 15.0**, with the alt-ergo, cvc5 and z3 provers it ships with. One way to get
it is `alr get gnatprove`, then put its `bin/` on your PATH. You don't need an Alire project for this
tree.

```sh
gnatprove --version            # expect: FSF 15.0
shasum -a 256 -c SHA256SUMS    # optional: check the sources are the recorded ones
gnatprove -P pde_kernels.gpr --level=2 --report=all -j4
```

`make prove` runs the same command, and `make checksums` runs the check.

What to expect:
- `gnatprove.out` (under `obj/gnatprove/`) reports **0 unproved checks**.
- It also lists four Big_Real conversion subprograms of the GNAT runtime as **skipped**, once for
  each of the two rounding units (8 lines). That is expected: RESULTS.md, "Honest limits", explains
  why. Nothing else is skipped.
- The run takes a few minutes on four cores. A heavily loaded machine can push one of the costlier
  goals over its time limit; if that happens, re-run on a quieter machine.

**Compiling into an executable.** The proofs model one rounding per floating-point operation. If you
build these units into a program, keep `-ffp-contract=off` (it's in `pde_kernels.gpr`). Without it
GCC may fuse `a*b + c` into an FMA, and the binary no longer computes what was proved. RESULTS.md
gives the measured divergence.

## Credits

Schemes from Lanyon AI's open repositories: AdvectionDiffusion (no licence file — we ported the
published scheme; no Lanyon code is included) and BurgersEquation (MIT, attribution below).

- **AdvectionDiffusion** (github.com/lanyonai/AdvectionDiffusion, `linear_advection_1d`). The
  repository has no licence file, so we ported the mathematics of the published scheme. None of its
  C, Lean or Racket text is in this repository. Its Lean theorems are cited by name only.
- **BurgersEquation** (github.com/lanyonai/BurgersEquation, `inviscid_burgers_1d`, at commit
  dd9b081). MIT licence. The Burgers units follow the reference implementation's scheme and its
  order of floating-point operations. None of its source files are included. Its notice:

```
MIT License

Copyright (c) 2026 Lanyon AI Inc

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

Related public discussion: lanyonai/AdvectionDiffusion#4 and lanyonai/BurgersEquation#2.

## Licence

The SPARK code and documents in this repository are Copyright (C) 2026 The Dark Factory Ltd, licensed under the
GNU Affero General Public License v3.0 or later (see `LICENSE`). The MIT notice above covers the
BurgersEquation scheme these units follow.
