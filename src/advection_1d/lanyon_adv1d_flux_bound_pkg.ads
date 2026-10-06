with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Adv1d_Rounding_Pkg; use Lanyon_Adv1d_Rounding_Pkg;
with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
with Lanyon_Adv1d_Waves_Pkg; use Lanyon_Adv1d_Waves_Pkg;
with Lanyon_Adv1d_Jump_Bound_Pkg; use Lanyon_Adv1d_Jump_Bound_Pkg;
package Lanyon_Adv1d_Flux_Bound_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  THE LEAN xFluxConservative IN FLOATS, part 1 of 3: the facts about the waves, the products a f, the flux
   --  difference and the quotient, as exact reals; the accessors; the magnitude facts.

   --  shorthands for the computed floats
   function W1 (UL, UR : State; P : Parameters) return Long_Float is (X_Wave_Family (UL, UR, P).Wave1)
     with Pre => Parameters_Valid (P);
   function W2 (UL, UR : State; P : Parameters) return Long_Float is (X_Wave_Family (UL, UR, P).Wave2)
     with Pre => Parameters_Valid (P);
   function A_R (P : Parameters) return Valid_Big_Real is (R (Abs_A (P)))
     with Pre => Parameters_Valid (P);

   --  |a| as a real lies in [9e-7, 1000.1]
   function A_Range (P : Parameters) return Boolean is
     (Abs_Real (P) and then R (Abs_A (P)) <= R (1.0E3) and then R (1.0E3) = 1000.0)
     with Pre  => Parameters_Valid (P),
          Post => A_Range'Result and then A_R (P) >= 0.0000009 and then A_R (P) <= 1000.1
                  and then A_R (P) = abs R (P.A);

   --  the waves, as reals
   function W1_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Half_Exact and then T1 (UL, UR, P) and then Mul_Err (0.5, Jump (UL, UR) - Quotient (UL, UR, P)))
     with Pre  => Parameters_Valid (P),
          Post => W1_Fact'Result
                  and then abs (R (W1 (UL, UR, P)) - 0.5 * R (Jump (UL, UR) - Quotient (UL, UR, P)))
                           <= Un * abs (0.5 * R (Jump (UL, UR) - Quotient (UL, UR, P))) + Eta;

   function W2_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Half_Exact and then T1b (UL, UR, P) and then Mul_Err (0.5, Jump (UL, UR) + Quotient (UL, UR, P)))
     with Pre  => Parameters_Valid (P),
          Post => W2_Fact'Result
                  and then abs (R (W2 (UL, UR, P)) - 0.5 * R (Jump (UL, UR) + Quotient (UL, UR, P)))
                           <= Un * abs (0.5 * R (Jump (UL, UR) + Quotient (UL, UR, P))) + Eta;

   --  the two products, the flux difference and the quotient, as reals
   function P1_Fact (UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Mul_Err (P.A, UR.F))
     with Post => P1_Fact'Result
                  and then abs (R (P.A * UR.F) - R (P.A) * R (UR.F)) <= Un * abs (R (P.A) * R (UR.F)) + Eta;

   function N_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Flux_Diff_Real (UL, UR, P) and then Sub_Err (P.A * UR.F, P.A * UL.F))
     with Post => N_Fact'Result
                  and then abs (R (Flux_Difference (UL, UR, P)) - (R (P.A * UR.F) - R (P.A * UL.F)))
                           <= Un * abs (R (P.A * UR.F) - R (P.A * UL.F)) + Eta;

   function Q_Div_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Div_Err (Flux_Difference (UL, UR, P), Abs_A (P)))
     with Pre  => Parameters_Valid (P),
          Post => Q_Div_Fact'Result
                  and then abs (R (Quotient (UL, UR, P)) - R (Flux_Difference (UL, UR, P)) / A_R (P))
                           <= Un * abs (R (Flux_Difference (UL, UR, P)) / A_R (P)) + Eta;


   --  the jump unit's facts, re-exported so that the later units need not with it
   function D_Fx (UL, UR : State) return Boolean is
     (D_Fact (UL, UR))
     with Post => D_Fx'Result
                  and then abs (R (Jump (UL, UR)) - (R (UR.F) - R (UL.F))) <= Un * abs (R (UR.F) - R (UL.F)) + Eta;

   function S1_Fx (UL, UR : State; P : Parameters) return Boolean is
     (S1_Fact (UL, UR, P))
     with Pre  => Parameters_Valid (P),
          Post => S1_Fx'Result
                  and then abs (R (Jump (UL, UR) - Quotient (UL, UR, P)) - (R (Jump (UL, UR)) - R (Quotient (UL, UR, P))))
                           <= Un * abs (R (Jump (UL, UR)) - R (Quotient (UL, UR, P))) + Eta;

   function S2_Fx (UL, UR : State; P : Parameters) return Boolean is
     (S2_Fact (UL, UR, P))
     with Pre  => Parameters_Valid (P),
          Post => S2_Fx'Result
                  and then abs (R (Jump (UL, UR) + Quotient (UL, UR, P)) - (R (Jump (UL, UR)) + R (Quotient (UL, UR, P))))
                           <= Un * abs (R (Jump (UL, UR)) + R (Quotient (UL, UR, P))) + Eta;

   function M (UL, UR : State) return Valid_Big_Real is (abs R (UR.F) + abs R (UL.F));
   function K (UL, UR : State) return Valid_Big_Real is
     ((1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M (UL, UR) + 5000000.0 * Eta);
   function Rd (UL, UR : State) return Valid_Big_Real is (R (Jump (UL, UR)));
   function Rq (UL, UR : State; P : Parameters) return Valid_Big_Real is (R (Quotient (UL, UR, P)))
     with Pre => Parameters_Valid (P);
   function Rs1 (UL, UR : State; P : Parameters) return Valid_Big_Real is (R (Jump (UL, UR) - Quotient (UL, UR, P)))
     with Pre => Parameters_Valid (P);
   function Rs2 (UL, UR : State; P : Parameters) return Valid_Big_Real is (R (Jump (UL, UR) + Quotient (UL, UR, P)))
     with Pre => Parameters_Valid (P);
   function Rw1 (UL, UR : State; P : Parameters) return Valid_Big_Real is (R (W1 (UL, UR, P)))
     with Pre => Parameters_Valid (P);
   function Rw2 (UL, UR : State; P : Parameters) return Valid_Big_Real is (R (W2 (UL, UR, P)))
     with Pre => Parameters_Valid (P);
   function Rn (UL, UR : State; P : Parameters) return Valid_Big_Real is (R (Flux_Difference (UL, UR, P)));
   function Rp1 (UR : State; P : Parameters) return Valid_Big_Real is (R (P.A * UR.F));


   --  the magnitudes in terms of M and K
   function Mag_Fact (UL, UR : State; P : Parameters) return Boolean is
     (D_Fact (UL, UR) and then S1_Fact (UL, UR, P) and then S2_Fact (UL, UR, P) and then Q_Fact (UL, UR, P)
      and then W1_Fact (UL, UR, P) and then W2_Fact (UL, UR, P) and then Eta_Small
      and then abs (R (UR.F) - R (UL.F)) <= M (UL, UR)
      and then abs Rd (UL, UR) <= K (UL, UR)
      and then Mags (Rw1 (UL, UR, P), Rw2 (UL, UR, P), Rs1 (UL, UR, P), Rs2 (UL, UR, P),
                     Rd (UL, UR), Rq (UL, UR, P), M (UL, UR), Eta, K (UL, UR)))
     with Pre  => Parameters_Valid (P),
          Post => Mag_Fact'Result
                  and then abs Rw1 (UL, UR, P) + abs Rw2 (UL, UR, P) <= 1.001 * K (UL, UR) + 4.0 * Eta
                  and then abs (Rd (UL, UR) - Rq (UL, UR, P)) + abs (Rd (UL, UR) + Rq (UL, UR, P)) <= 2.0 * K (UL, UR)
                  and then abs Rs1 (UL, UR, P) + abs Rs2 (UL, UR, P) <= 2.001 * K (UL, UR) + 2.0 * Eta;

   --  the flux-side magnitudes: |N| and |P1 - P2|
   function Num_Fact (UL, UR : State; P : Parameters) return Boolean is
     (Quotient_Bound (UL, UR, P) and then A_Range (P) and then P1_Fact (UR, P) and then P1_Fact (UL, P)
      and then Prod_Mag (Rp1 (UR, P), R (P.A), R (UR.F), Eta) and then Prod_Mag (Rp1 (UL, P), R (P.A), R (UL.F), Eta)
      and then U_Is_Un and then N_Fact (UL, UR, P)
      and then Num_Mag (Rn (UL, UR, P), Rp1 (UR, P), Rp1 (UL, P), A_R (P), R (UR.F), R (UL.F), Eta)
      and then Diff_Mag (Rp1 (UR, P), Rp1 (UL, P), R (P.A), A_R (P), R (UR.F), R (UL.F), Eta))
     with Pre  => Parameters_Valid (P),
          Post => Num_Fact'Result
                  and then abs Rn (UL, UR, P)
                           <= (1.0 + Un) * (1.0 + Un) * (A_R (P) * abs R (UR.F) + A_R (P) * abs R (UL.F)) + 4.0 * Eta
                  and then abs (Rp1 (UR, P) - Rp1 (UL, P))
                           <= (1.0 + Un) * (A_R (P) * abs R (UR.F) + A_R (P) * abs R (UL.F)) + 2.0 * Eta;

end Lanyon_Adv1d_Flux_Bound_Pkg;
