with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
with Lanyon_Rb_Flux_Pkg; use Lanyon_Rb_Flux_Pkg;
with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Quot_Pkg; use Lanyon_Burg1d_Quot_Pkg;
with Lanyon_Burg1d_Waves_Pkg; use Lanyon_Burg1d_Waves_Pkg;
with Lanyon_Burg1d_Jump_Bound_Pkg; use Lanyon_Burg1d_Jump_Bound_Pkg;
package Lanyon_Burg1d_Flux_Facts_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081.
   --  THE LEAN xFluxConservative IN FLOATS for Burgers, part 1 of 3: the facts as exact reals.

   --  the local copy of the jump unit's regime floor
   Floor_B : constant := 2.0 ** (-511);

   --  the jump unit's s1 / s2 facts, re-exported so that the later units need not with it
   function S1_Fx (UL, UR : State) return Boolean is
     (S1_Fact (UL, UR))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => S1_Fx'Result
                  and then abs (R (Jump (UL, UR) - Quotient (UL, UR)) - (R (Jump (UL, UR)) - R (Quotient (UL, UR))))
                           <= Un * abs (R (Jump (UL, UR)) - R (Quotient (UL, UR))) + Eta;

   function S2_Fx (UL, UR : State) return Boolean is
     (S2_Fact (UL, UR))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => S2_Fx'Result
                  and then abs (R (Jump (UL, UR) + Quotient (UL, UR)) - (R (Jump (UL, UR)) + R (Quotient (UL, UR))))
                           <= Un * abs (R (Jump (UL, UR)) + R (Quotient (UL, UR))) + Eta;



   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. THE LEAN xFluxConservative IN FLOATS for Burgers: the computed fluctuations sum to the exact flux jump
   --  f_R^2/2 - f_L^2/2 within 10u A (|f_L| + |f_R|) + 1e7 eta, A = max (|f_L|, |f_R|) >= 2**-511.
   --  shorthands for the computed values, as reals
   function A_R (UL, UR : State) return Valid_Big_Real is (R (Amax (UL, UR)));

   function M (UL, UR : State) return Valid_Big_Real is (abs R (UR.F) + abs R (UL.F));

   function K (UL, UR : State) return Valid_Big_Real is
     ((1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M (UL, UR) + 5000000.0 * Eta);

   function W1 (UL, UR : State) return Wave_Value is (X_Wave_Family (UL, UR).Wave1);

   function W2 (UL, UR : State) return Wave_Value is (X_Wave_Family (UL, UR).Wave2);

   function Rd (UL, UR : State) return Valid_Big_Real is (R (Jump (UL, UR)));

   function Rq (UL, UR : State) return Valid_Big_Real is (R (Quotient (UL, UR)))
     with Pre => Amax (UL, UR) > 0.0 and then Quot_Safe (UL.F, UR.F);

   function Rs1 (UL, UR : State) return Valid_Big_Real is (R (Jump (UL, UR) - Quotient (UL, UR)))
     with Pre => Amax (UL, UR) > 0.0 and then Quot_Safe (UL.F, UR.F);

   function Rs2 (UL, UR : State) return Valid_Big_Real is (R (Jump (UL, UR) + Quotient (UL, UR)))
     with Pre => Amax (UL, UR) > 0.0 and then Quot_Safe (UL.F, UR.F);

   function Rn (UL, UR : State) return Valid_Big_Real is (R (Num (UL.F, UR.F)));

   function Rh (F : Field) return Valid_Big_Real is (R (0.5 * (F * F)));

   function Rp (F : Field) return Valid_Big_Real is (R (F * F));

   --  the exact flux jump of the Lean: f_R^2/2 - f_L^2/2
   function Jx (UL, UR : State) return Valid_Big_Real is
     (0.5 * (R (UR.F) * R (UR.F)) - 0.5 * (R (UL.F) * R (UL.F)));

   --  the waves, as reals
   function W1_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Half_Exact and then T1 (UL, UR) and then Mul_Err (0.5, Jump (UL, UR) - Quotient (UL, UR)))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => W1_Fact'Result
                  and then abs (R (W1 (UL, UR)) - 0.5 * Rs1 (UL, UR)) <= Un * abs (0.5 * Rs1 (UL, UR)) + Eta;

   function W2_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Half_Exact and then T1b (UL, UR) and then Mul_Err (0.5, Jump (UL, UR) + Quotient (UL, UR)))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => W2_Fact'Result
                  and then abs (R (W2 (UL, UR)) - 0.5 * Rs2 (UL, UR)) <= Un * abs (0.5 * Rs2 (UL, UR)) + Eta;

   --  one square, as a real
   function Sq_Fact (F : Field) return Boolean is
     (U_Is_Un and then Mul_Err (F, F))
     with Post => Sq_Fact'Result and then abs (Rp (F) - R (F) * R (F)) <= Un * abs (R (F) * R (F)) + Eta;

   --  one half, as a real
   function Half_Fact (F : Field) return Boolean is
     (U_Is_Un and then Half_Exact and then Mul_Err (0.5, F * F))
     with Post => Half_Fact'Result and then abs (Rh (F) - 0.5 * Rp (F)) <= Un * abs (0.5 * Rp (F)) + Eta;

   --  the difference, as a real
   function Diff_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Sub_Err (0.5 * (UR.F * UR.F), 0.5 * (UL.F * UL.F)))
     with Post => Diff_Fact'Result
                  and then abs (Rn (UL, UR) - (Rh (UR.F) - Rh (UL.F))) <= Un * abs (Rh (UR.F) - Rh (UL.F)) + Eta;

   --  the numerator chain, as reals: squares, halves, difference
   function Num_Chain (UL, UR : State) return Boolean is
     (Sq_Fact (UR.F) and then Sq_Fact (UL.F) and then Half_Fact (UR.F) and then Half_Fact (UL.F)
      and then Diff_Fact (UL, UR))
     with Post => Num_Chain'Result
                  and then abs (Rp (UR.F) - R (UR.F) * R (UR.F)) <= Un * abs (R (UR.F) * R (UR.F)) + Eta
                  and then abs (Rp (UL.F) - R (UL.F) * R (UL.F)) <= Un * abs (R (UL.F) * R (UL.F)) + Eta
                  and then abs (Rh (UR.F) - 0.5 * Rp (UR.F)) <= Un * abs (0.5 * Rp (UR.F)) + Eta
                  and then abs (Rh (UL.F) - 0.5 * Rp (UL.F)) <= Un * abs (0.5 * Rp (UL.F)) + Eta
                  and then abs (Rn (UL, UR) - (Rh (UR.F) - Rh (UL.F))) <= Un * abs (Rh (UR.F) - Rh (UL.F)) + Eta;

   --  pure real: the computed numerator is the exact flux jump, to rounding
   function N_Near_B (N, H1, H2, P1, P2, FR, FL, E : Valid_Big_Real) return Boolean is
     (abs (N - (0.5 * (FR * FR) - 0.5 * (FL * FL)))
        <= Un * abs (H1 - H2) + 0.5 * Un * (abs P1 + abs P2) + 0.5 * Un * (FR * FR + FL * FL) + 4.0 * E)
     with Pre  => E >= 0.0
                  and then abs (P1 - FR * FR) <= Un * abs (FR * FR) + E
                  and then abs (P2 - FL * FL) <= Un * abs (FL * FL) + E
                  and then abs (H1 - 0.5 * P1) <= Un * abs (0.5 * P1) + E
                  and then abs (H2 - 0.5 * P2) <= Un * abs (0.5 * P2) + E
                  and then abs (N - (H1 - H2)) <= Un * abs (H1 - H2) + E,
          Post => N_Near_B'Result;

   --  pure real: the magnitudes of the numerator chain against A^2
   function Chain_Mags (H1, H2, P1, P2, FR, FL, A, E : Valid_Big_Real) return Boolean is
     (abs (H1 - H2) <= 0.51 * (A * A) + 2.0 * E and then abs P1 + abs P2 <= 2.01 * (A * A) + 2.0 * E)
     with Pre  => E >= 0.0 and then E <= 0.000000000000000000001 and then H1 >= 0.0 and then H2 >= 0.0
                  and then FR * FR <= A * A and then FL * FL <= A * A
                  and then abs (P1 - FR * FR) <= Un * abs (FR * FR) + E
                  and then abs (P2 - FL * FL) <= Un * abs (FL * FL) + E
                  and then abs (H1 - 0.5 * P1) <= Un * abs (0.5 * P1) + E
                  and then abs (H2 - 0.5 * P2) <= Un * abs (0.5 * P2) + E,
          Post => Chain_Mags'Result;

   --  the magnitudes, against M and K
   function Mag_Fact (UL, UR : State) return Boolean is
     (Quot_Safe (UL.F, UR.F)
      and then D_Fact (UL, UR) and then S1_Fact (UL, UR) and then S2_Fact (UL, UR) and then Q_Fact (UL, UR)
      and then W1_Fact (UL, UR) and then W2_Fact (UL, UR) and then Eta_Small
      and then abs (R (UR.F) - R (UL.F)) <= M (UL, UR)
      and then abs Rd (UL, UR) <= K (UL, UR)
      and then Mags (R (W1 (UL, UR)), R (W2 (UL, UR)), Rs1 (UL, UR), Rs2 (UL, UR),
                     Rd (UL, UR), Rq (UL, UR), M (UL, UR), Eta, K (UL, UR)))
     with Pre  => Amax (UL, UR) >= Floor_N,
          Post => Mag_Fact'Result
                  and then abs R (W1 (UL, UR)) + abs R (W2 (UL, UR)) <= 1.001 * K (UL, UR) + 4.0 * Eta
                  and then abs (Rd (UL, UR) - Rq (UL, UR)) + abs (Rd (UL, UR) + Rq (UL, UR)) <= 2.0 * K (UL, UR)
                  and then abs Rs1 (UL, UR) + abs Rs2 (UL, UR) <= 2.001 * K (UL, UR) + 2.0 * Eta;

   --  the numerator chain against A^2, and against the exact flux jump
   function Chain_Fact (UL, UR : State) return Boolean is
     (Num_Chain (UL, UR) and then Eta_Small and then Max_Real (UL.F, UR.F)
      and then Sq_Le (R (UR.F), A_R (UL, UR)) and then Sq_Le (R (UL.F), A_R (UL, UR))
      and then R_Nonneg (0.5 * (UR.F * UR.F)) and then R_Nonneg (0.5 * (UL.F * UL.F))
      and then Chain_Mags (Rh (UR.F), Rh (UL.F), Rp (UR.F), Rp (UL.F), R (UR.F), R (UL.F), A_R (UL, UR), Eta)
      and then N_Near_B (Rn (UL, UR), Rh (UR.F), Rh (UL.F), Rp (UR.F), Rp (UL.F), R (UR.F), R (UL.F), Eta)
      and then Num_Real_B (Rn (UL, UR), Rh (UR.F), Rh (UL.F), Rp (UR.F), Rp (UL.F), R (UR.F), R (UL.F),
                           A_R (UL, UR), Eta))
     with Post => Chain_Fact'Result
                  and then R (UR.F) * R (UR.F) <= A_R (UL, UR) * A_R (UL, UR)
                  and then R (UL.F) * R (UL.F) <= A_R (UL, UR) * A_R (UL, UR)
                  and then abs (Rh (UR.F) - Rh (UL.F)) <= 0.51 * (A_R (UL, UR) * A_R (UL, UR)) + 2.0 * Eta
                  and then abs Rp (UR.F) + abs Rp (UL.F) <= 2.01 * (A_R (UL, UR) * A_R (UL, UR)) + 2.0 * Eta
                  and then abs Rn (UL, UR) <= 0.6 * (A_R (UL, UR) * A_R (UL, UR)) + 4.0 * Eta
                  and then abs (Rn (UL, UR) - Jx (UL, UR))
                           <= Un * abs (Rh (UR.F) - Rh (UL.F)) + 0.5 * Un * (abs Rp (UR.F) + abs Rp (UL.F))
                              + 0.5 * Un * (R (UR.F) * R (UR.F) + R (UL.F) * R (UL.F)) + 4.0 * Eta;

   --  the divisor: 0 < A <= 1e6, A <= M
   function A_Fact (UL, UR : State) return Boolean is
     (Eta_Le (Amax (UL, UR)) and then Max_Le_Sum (UL, UR))
     with Pre  => Amax (UL, UR) >= Floor_N,
          Post => A_Fact'Result and then A_R (UL, UR) > 0.0 and then A_R (UL, UR) <= 1000000.0
                  and then A_R (UL, UR) <= M (UL, UR);

   --  NEW STONE (the extra lemma for 3/3): the quotient error, as reals, with its premises stated once
   function QD_Fx (UL, UR : State) return Boolean is
     (Quot_Safe (UL.F, UR.F) and then Num_Fact (UL.F, UR.F) and then Eta_Le (Amax (UL, UR))
      and then U_Is_Un and then Div_Any (Num (UL.F, UR.F), Amax (UL, UR)))
     with Pre  => Amax (UL, UR) >= Floor_B,
          Post => QD_Fx'Result
                  and then abs (Rq (UL, UR) - Rn (UL, UR) / A_R (UL, UR)) <= Un * abs (Rn (UL, UR) / A_R (UL, UR)) + Eta;


end Lanyon_Burg1d_Flux_Facts_Pkg;
