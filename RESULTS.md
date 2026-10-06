# Results: which real-number theorems survive into proven floating-point code

Two of Lanyon AI's 1-D schemes, ported to SPARK and proved with gnatprove FSF 15.0 at level 2. For
every Lean theorem, this file records the SPARK contract on the executable code and the form that
survives in IEEE binary64: **exact**, **bounded** (with the proved ε), or **not provable as stated**.

Notation:
- u = 2⁻⁵³, the unit roundoff of binary64;
- η = 2⁻¹⁰⁷⁴, the smallest subnormal;
- R(x) is the exact real value of the double x (a ghost `Big_Real`).

The Lean theorems are cited by **name** and paraphrased. No Lanyon text is reproduced.

## How these results were produced, and what "admitted" means

1. Each unit was specified in prose. A local language model emitted the SPARK spec from that prose.
2. gnatprove decided at `--level=2`, with the switches in `pde_kernels.gpr`. A round passes only with
   0 unproved checks, no assumed bodies (with one narrow, stated exemption: see "Honest limits"), and
   an emission identical to what the prose specified.
3. A **separate machine** then re-proved the unit (step limit 100,000, timeout 120 s) and appended a
   signed receipt to an append-only ledger. The ledger number is the entry for that unit.

Every `.ads` here is byte-identical to the file the ledger records. The sha256 of each file equals
the ledger's `source_sha256`. `SHA256SUMS` lists them, and the table below pairs each with its ledger
entry.

### The units, with ledger entries

**Advection** (`src/advection_1d/`)

| Package | Carries (rows below) | Ledger |
|---|---|---|
| Lanyon_Adv1d_Types_Pkg | A1-A4, A17 | 1003 |
| Lanyon_Adv1d_Waves_Pkg | A5, A17 | 1004 |
| Lanyon_Adv1d_Fluct_Pkg | A7, A8, A17 | 1005 |
| Lanyon_Adv1d_Update_Pkg | A19, A21 | 1006 |
| Lanyon_Adv1d_Exact_Pkg | A20 | 1007 |
| Lanyon_Adv1d_Minmod_Pkg | A10-A13 | 1008 |
| Lanyon_Adv1d_Rounding_Pkg (ghost: per-operation rounding lemmas) | basis of A6, A9 | 1009 |
| Lanyon_Adv1d_Jump_Bound_Pkg | **A6** | 1010 |
| Lanyon_Adv1d_Flux_Bound_Pkg (facts) | A9, part 1 of 3 | 1018 |
| Lanyon_Adv1d_Flux_Pieces_Pkg (pieces) | A9, part 2 of 3 | 1019 |
| Lanyon_Adv1d_Flux_Final_Pkg (the theorem) | **A9** | 1020 |

**Burgers** (`src/burgers_1d/`)

| Package | Carries | Ledger |
|---|---|---|
| Lanyon_Burg1d_Types_Pkg | B1, B2, B4; the sonic-point characterisation | 1011 |
| Lanyon_Burg1d_Exact_Pkg | level 3, exact variant | 1012 |
| Lanyon_Burg1d_Waves_Pkg | B5; the sonic point defined | 1014 |
| Lanyon_Burg1d_Fluct_Pkg | B7, B8; Sonic_Quiet | 1015 |
| Lanyon_Burg1d_Update_Pkg | interface shared, constant fixed point | 1016 |
| Lanyon_Rb_Core_Pkg (ghost: per-operation rounding lemmas) | basis of B6, B9 | 1021 |
| Lanyon_Rb_Jump_Pkg (ghost) | basis of B6 | 1022 |
| Lanyon_Rb_Flux_Pkg (ghost) | basis of B9 | 1023 |
| Lanyon_Burg1d_Quot_Pkg (ghost) | level 1 with no floor on the state | 1024 |
| Lanyon_Burg1d_Jump_Bound_Pkg | **B6** | 1025 |
| Lanyon_Burg1d_Flux_Facts_Pkg | B9, part 1 of 3 | 1026 |
| Lanyon_Burg1d_Flux_Pieces_Pkg | B9, part 2 of 3 | 1027 |
| Lanyon_Burg1d_Flux_Final_Pkg (the theorem) | **B9** | 1028 |

**Earlier ledger entries.** Ledger 1013 and 1017 are earlier admissions of `Lanyon_Burg1d_Quot_Pkg` and
`Lanyon_Burg1d_Jump_Bound_Pkg`. Those versions `with` the advection rounding unit. The Burgers rounding
lemmas were later split into their own three units (1021-1023), so that every proof unit stays within
a size limit, and the two packages were re-admitted on them as 1024 and 1025. The texts differ only in
their `with` lines. The two versions share package names and can't sit in one project, so this tree
ships the current ones.

## Proof run of this tree

Run with exactly the command in README.md, on a 4-core build machine (CPU only):

```
$ gnatprove --version          # FSF 15.0 (Why3 1.7.1+git, Alt-Ergo 2.6.0)
$ gnatprove -P pde_kernels.gpr --level=2 --report=all -j4
                               Total   Flow        Provers                                  Justified  Unproved
Run-time Checks                1063      .    1063 (CVC5 10%, Trivial 85%, Z3 4%, altergo 1%)     .          .
Functional Contracts            691      .     691 (CVC5 43%, Trivial 1%, Z3 53%, altergo 3%)     .          .
Termination                     302    298       4 (CVC5 25%, Z3 75%)                             .          .
Total                          2056    298 (14%) 1758 (86%)                                       .          .
max steps used for successful proof: 50802
wall clock 183 s (user 878 s), 4 jobs, CPU only, 2026-10-06

Skipped: exactly 8 lines, the four Big_Real conversions once for each of the two
Float_Conversions instantiations (Lanyon_Adv1d_Rounding_Pkg.LF and Lanyon_Rb_Core_Pkg.LF).
No other subprogram is skipped. 0 justified, 0 unproved.
```

## Advection 1-D (Lanyon `linear_advection_1d`)

Declared ranges: the state |f| ≤ 10⁶ and the speed |a| ≤ 10³. Validity requires |a| ≥ A_Min = 10⁻⁶;
see A17 for why it isn't a ≠ 0.

### Part 1: the Lean theorems on the executable

| # | Lean theorem | Says, over the reals | SPARK contract (package . function) | Form in floats |
|---|---|---|---|---|
| A1 | xHyperbolicity | the flux Jacobian has a real eigenvalue | Types . X_Eigen returns a finite Speed; **plus** Types . Flux_Is_Linear: flux = eigenvalue × state, exactly | **exact**. The added theorem ties the eigenvalue to the flux, so it is not vacuous (cf. AdvectionDiffusion#2). A mutant eigenvalue of 42 is refused by the prover |
| A2 | xWaveStability | \|wave speed\| ≥ \|eigenvalue\| | Types . Wave_Stability | **exact** |
| A3 | xDiffusiveFluxConsistency | zero gradient gives zero diffusive flux | Types . Diffusive_Consistency | **exact** |
| A4 | (speed definitions) | the two speeds are −\|a\| and \|a\| | Types . Speeds_Opposite | **exact** |
| A5 | xWaveConsistency | equal states give two zero waves | Waves . Wave_Consistency | **exact** (x − x = 0 and 0 / \|a\| = 0 in IEEE). It describes the binary only for a build without FMA contraction (A22) |
| A6 | xWaveJumpCondition | jump = wave1 + wave2, for **all** a | Jump_Bound . Wave_Jump_Bounded: \|R(w1) + R(w2) − (R(f_R) − R(f_L))\| ≤ 6u(\|f_L\| + \|f_R\|) + 6η | **bounded**: 6u proved, 2.00u the largest measured over 2,000,000 cases. **Not as stated**: at a = 0 the float waves are NaN. Over the reals it holds at a = 0 only because Mathlib defines x / 0 = 0 |
| A7 | xLeftFluctuationsConsistent | equal states give a zero left fluctuation | Fluct . Left_Consistent | **exact** |
| A8 | xRightFluctuationsConsistent | as A7, right | Fluct . Right_Consistent | **exact** |
| A9 | xFluxConservative | flux jump = left + right fluctuation (a ≠ 0) | Flux_Final . Flux_Conservative_Bounded: \|R(𝓛) + R(𝓡) − R(a)(R(f_R) − R(f_L))\| ≤ 10u\|a\|(\|f_L\| + \|f_R\|) + 10⁴η | **bounded**: 10u proved, 5.71u the largest measured |
| A10 | xLeftReconstructionConsistent | equal states reconstruct f | Minmod . Left_Consistent | **exact** |
| A11 | xRightReconstructionConsistent | as A10 | Minmod . Right_Consistent | **exact** |
| A12 | xLeftReconstructionLinearityPreservation | inputs f − D, f, f + D reconstruct f − D/2 | Minmod . Left_Linear, **restated**: if the two computed slopes are equal, the result is f − slope/2 | **not provable as stated**: in floats, f − D and f + D are themselves rounded, so the hypothesis has no float counterpart. **Exact** in the restated form |
| A13 | xRightReconstructionLinearityPreservation | as A12, + D/2 | Minmod . Right_Linear, restated as A12 | as A12 |
| A14 | xReconstructionSymmetric | left recon of (L, U, R) = right recon of (R, U, L) | (not in this tree) | **exact** in floats: it rests on fl(−x) = −fl(x) for a difference, and it was proved. **On hold, not included**: see "Honest limits" |

Beyond the Lean:

| # | Property | Form |
|---|---|---|
| A15 | reversing the speed swaps the two waves | **exact**, bit for bit. Proved, but on hold with A14 |
| A16 | upwinding: for a > 0 the whole jump goes to wave 2 (wave 1 = 0) | **false in floats**: wave 1 is exactly 0 in only 73.75 % of a > 0 cases (measured). AdvectionDiffusion#5 is the 1-ulp case. Not a contract |

### Part 2: absence of runtime errors, conservation, and the build

| # | Obligation | SPARK contract | Result |
|---|---|---|---|
| A17 | **Level 1**: no overflow, no division by zero, no NaN, for every input in the declared ranges | every function, every unit | **achieved**, with validity \|a\| ≥ 10⁻⁶. With a ≠ 0 instead, subnormal a is admitted, and the quotient's range becomes a relational float fact. It proves only at level 4 (cvc5, about 7 s per check), not at level 2. a = 0, −0, subnormal a and 5·10⁻⁷ are refused by name |
| A18 | **Level 3, floats**: the cell sum after a step equals the sum before (periodic) | not a contract | **not provable as stated: false in floats.** Summation and update rounding break it in 197 of 200 random grids. A bound would need error accumulation over N cells |
| A19 | a face's outgoing flux = its neighbour's incoming flux | Update . Interface_Shared | **exact**: the same value, computed once per face |
| A20 | **Level 3, exact variant**: integer cells (\|q\| ≤ 2⁴⁰), Courant number C/2¹⁶; the sum after a step **equals** the sum before, for every grid and Courant number | Exact . Conserved, by telescoping induction (Exact . Telescope) | **exact** |
| A21 | a constant grid is a fixed point of the float step | Update . Constant_Fixed | **exact** |
| A22 | the proof describes the executable | build with -ffp-contract=off | **only for a contraction-off build.** GNAT's default build fuses 54 multiply-add sites and then diverges from the proved one-rounding-per-operation model: 92 of 200,006 random pairs, and 160 of 200 grid steps |
| A23 | the maximum principle (CFL) | not attempted | relational, like A6. One step is proved safe into a wider cell type only |

## Inviscid Burgers 1-D (Lanyon `inviscid_burgers_1d`)

Why it is harder:
- The flux f²/2 is nonlinear.
- The wave split divides by amax = max(|f_L|, |f_R|), the state itself, so there is no parameter to
  put a floor on.
- The NaN case is the **sonic point**, where both states are exactly 0. It occurs on the standard
  sin(2πx) test, at interface 99|100.

Declared range |f| ≤ 10⁶. **There is no floor on the state.** At the sonic point the waves are
defined as (0, 0), their exact real limit, since |q| ≤ |f_R − f_L|.

| # | Lean theorem | SPARK contract | Form in floats |
|---|---|---|---|
| B1 | xHyperbolicity | Types . Flux_Is_Jacobian: flux = (0.5 × eigenvalue) × state, exactly | **exact**. It ties the eigenvalue f to the flux f²/2 |
| B2 | xWaveStability | Types . Wave_Stability | **exact** |
| B3 | xDiffusiveFluxConsistency | not ported: the diffusive flux is identically 0 | **exact**, trivially |
| B4 | (speed definitions) ±amax | Types . Speeds_Opposite | **exact** |
| B5 | xWaveConsistency | Waves . Wave_Consistency | **exact**, including the sonic state |
| B6 | xWaveJumpCondition (no hypothesis) | Jump_Bound . Wave_Jump_Bounded: \|R(w1) + R(w2) − (R(f_R) − R(f_L))\| ≤ 6u(\|f_L\| + \|f_R\|) + 6η, for amax ≥ 2⁻⁵¹¹; **plus** Wave_Jump_Sonic: exact at the sonic point | **bounded** for amax ≥ 2⁻⁵¹¹. **Exact** at the sonic point, where the reference gives NaN. **Open** for 0 < amax < 2⁻⁵¹¹ (see "Honest limits") |
| B7 | xLeftFluctuationsConsistent | Fluct . Left_Consistent | **exact** |
| B8 | xRightFluctuationsConsistent | Fluct . Right_Consistent | **exact** |
| B9 | xFluxConservative (amax ≠ 0) | Flux_Final . Flux_Conservative_Bounded: \|R(𝓛) + R(𝓡) − (R(f_R²/2) − R(f_L²/2))\| ≤ 10u·amax·(\|f_L\| + \|f_R\|) + 10⁷η, for amax ≥ 2⁻⁵¹¹; **plus** Fluct . Sonic_Quiet: both fluctuations are exactly 0 at the sonic point | **bounded** for amax ≥ 2⁻⁵¹¹; **open** below |
| B10-B14 | the five minmod theorems | Lanyon's Burgers minmod is the same formula as the advection one, so the advection Minmod package covers it | as A10-A14 |

The brief's three levels, for Burgers:

| Level | Obligation | Result |
|---|---|---|
| 1 | no overflow, no division by zero, no NaN, for every input in range | **achieved with no floor on the state.** The quotient's float range is proved from a ghost real fact, \|N\| ≤ 2·amax² + 8η (Quot . Quot_Safe), which alt-ergo turns into the division's overflow check. Consequence: the Burgers cores are not `Pure`, because they `with` a ghost unit that uses Big_Reals |
| 2 | a face's outgoing flux = its neighbour's incoming flux; a constant grid is a fixed point | Update . Interface_Shared and Update . Constant_Fixed: **exact** |
| 3, floats | the cell sum is conserved | not exact, as for advection |
| 3, exact variant | integer cells (\|q\| ≤ 2²⁰), Murman-Roe upwind face, scaled flux C·(q²/2)/2¹⁶ | Exact . Conserved: **exact** conservation for every grid and every Courant number |

## Differential test: proven Ada against Lanyon's C

Method:
- The same random inputs go to Lanyon's C and to a build of these SPARK units, and the outputs are
  compared bit for bit.
- Both sides are compiled with contraction off.
- The seed is 20261005.

The harness is **not included** here, because it compiles Lanyon's C into the test driver. To
reproduce it, fetch their sources from their repositories and compare your own driver of these units
against them.

| Scheme | Cases | DIVERGENT | Notes |
|---|---|---|---|
| Advection | 252,213 non-refused | **0** | Bit-identical. The planted a = 0 cases are refused here and give NaN waves in the reference |
| Burgers | 252,219 | **0** | 12,354 pairs and 172 grids are **sonic-defined**: the reference gives NaN, this code gives 0. Every one of those pairs has both inputs exactly ±0. The minmod comparison is 50,002 identical. The exact variant conserves the sum in all 2,000 cases |

Related measurements:
- **FMA contraction** (AdvectionDiffusion#6, consistency on equal states). Over 60,000 equal-state
  pairs, the reference built with default contraction gives nonzero waves in all 60,000 (at bit level;
  its 1e-8 tolerance hides most of them). Contraction off gives 0 of 60,000, on both sides.
- **AdvectionDiffusion#5** (a 1-ulp jump). This code reproduces the reference bit for bit:
  w1 = w2 = 2.22·10⁻¹⁶ on both sides. See the next section.

## What the proofs rule out, against the upstream issues

| Issue | Ruled out? |
|---|---|
| AdvectionDiffusion#2: hyperbolicity theorem not tied to the flux | **Yes, for these 1-D units.** Flux_Is_Linear (advection) and Flux_Is_Jacobian (Burgers) are exact |
| AdvectionDiffusion#4: a = 0 gives NaN waves; the validity check accepts it | **Yes, provably.** Every result is finite and in a constant-bounded range for every valid input. Advection refuses \|a\| < 10⁻⁶ by name. Burgers defines the sonic point instead, so no in-range input is refused |
| AdvectionDiffusion#5: a 1-ulp jump makes the flux difference exactly 0 and splits the jump evenly | **No.** That is a finite but wrong answer, not a runtime error. The property that would exclude it (upwinding, A16) is false in floats. What does hold: the bounds A6/A9 and the exact theorems |
| AdvectionDiffusion#6: FMA contraction breaks exact consistency | **By build flag and test, not by proof.** gnatprove's model can't see compiler flags. Build with -ffp-contract=off (A22) |
| AdvectionDiffusion#7: 2-D DG ordering | not applicable (2-D) |

## Honest limits

1. **Whole-domain float conservation is false** (A18). Only the integer variants conserve exactly
   (A20 and Burgers level 3). We prove no bound for the float sum.
2. **Burgers' bounded theorems need amax ≥ 2⁻⁵¹¹** (B6, B9).
   - Below that threshold both squares are subnormal. The per-operation lemma |err| ≤ u|x| + η then
     gives only η/amax, which can be up to 1: an absolute bound, not a relative one.
   - The fact that would close the gap, |fl(x) − x| ≤ |x| under round-to-nearest, is true. But
     cvc5, z3 and alt-ergo don't discharge it at level 2 or at level 4.
   - So the band 0 < amax < 2⁻⁵¹¹ is **open**. The sonic point amax = 0 itself is proved exact.
3. **The symmetry unit is on hold, so A14 and A15 are not in this tree.**
   - Both rest on fl(−x) = −fl(x) for a sum or difference. Only z3 proves that, in 2-4 s per goal,
     against a 5 s limit, so the unit proves reliably only on a lightly loaded machine.
   - It was proved but never admitted, so it isn't shipped.
4. **The proved constants are loose.**
   - Advection: 6u and 10u proved; 2.00u and 5.71u the largest measured. That is about 3× and 1.75×.
   - Burgers' additive terms are also large in η units (10⁷η). They are still astronomically small in
     absolute terms.
5. **Four runtime conversions are listed as "skipped" (the narrow exemption).**
   - The ghost rounding lemmas express exact values through
     `Ada.Numerics.Big_Numbers.Big_Reals.Float_Conversions (Long_Float)`. That instantiation pulls in
     four GNAT runtime subprograms whose bodies are not in SPARK: `To_Big_Real`, `From_Big_Real`, and
     `Conv.To_Big_Integer` / `Conv.From_Big_Integer`. gnatprove reports them as skipped.
   - **Why this is not an assumption.** gnatprove does not assume a contract for them. It interprets
     the float-to-real conversion natively, as the exact value of the double. `To_Big_Real` has no
     postcondition, yet the per-operation lemmas |R(x op y) − (R(x) op R(y))| ≤ u|R(x) op R(y)| + η
     prove directly from gnatprove's IEEE model.
   - **Why exempt them.** Our pipeline refuses any unit with a skipped (assumed) body. It exempts
     exactly these four, and only for gnatprove FSF 15.0 with the runtime specs at pinned hashes.
     Canaries confirm that `pragma Assume`, any other body-less subprogram, a user wrapper of the
     conversion, and a changed toolchain are all still refused.
   - **Measured alternative.** SPARKlib's Big_Reals does not avoid the issue: it gives 4 skips
     (regular) or 3 (light). The conversion has no SPARK body in any library.
   - This tree has two instantiations (one per rounding unit), so a run lists 8 skipped lines: these
     four, twice. No others.
6. **The proofs are about the spec-level model, under the declared ranges** (|f| ≤ 10⁶, |a| ≤ 10³).
   Inputs outside them are not covered. A real driver must refuse them, as the differential harness's
   driver does.
7. **One step, not many.** A step is proved safe into a wider cell type. Chaining steps would need
   the maximum principle (A23), which is not proved.
