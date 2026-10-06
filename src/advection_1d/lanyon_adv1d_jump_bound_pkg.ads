with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Adv1d_Rounding_Pkg; use Lanyon_Adv1d_Rounding_Pkg;
with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
with Lanyon_Adv1d_Waves_Pkg; use Lanyon_Adv1d_Waves_Pkg;
package Lanyon_Adv1d_Jump_Bound_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  THE LEAN xWaveJumpCondition IN FLOATS: the computed waves sum to the exact jump within
   --  6u (|f_L| + |f_R|) + 6 eta, proved from per-operation rounding lemmas.

   function T0 (UL, UR : State; P : Parameters) return Boolean is
     (0.5 * (Jump (UL, UR) - Quotient (UL, UR, P)) in Wave_Value
      and then X_Wave_Family (UL, UR, P).Wave1 = 0.5 * (Jump (UL, UR) - Quotient (UL, UR, P)))
     with Pre => Parameters_Valid (P), Post => T0'Result;
   function T1 (UL, UR : State; P : Parameters) return Boolean is
     (T0 (UL, UR, P) and then R_Cong (X_Wave_Family (UL, UR, P).Wave1, 0.5 * (Jump (UL, UR) - Quotient (UL, UR, P)))
      and then R (X_Wave_Family (UL, UR, P).Wave1) = R (0.5 * (Jump (UL, UR) - Quotient (UL, UR, P))))
     with Pre => Parameters_Valid (P), Post => T1'Result;
   function T2 (X : Big) return Boolean is
     (R (0.5) * R (X) = 0.5 * R (X))
     with Post => T2'Result;
   function T3 (X : Big) return Boolean is
     (Mul_Err (0.5, X) and then abs (R (0.5 * X) - 0.5 * R (X)) <= U * abs (0.5 * R (X)) + Eta)
     with Post => T3'Result;
   function T0b (UL, UR : State; P : Parameters) return Boolean is
     (0.5 * (Jump (UL, UR) + Quotient (UL, UR, P)) in Wave_Value
      and then X_Wave_Family (UL, UR, P).Wave2 = 0.5 * (Jump (UL, UR) + Quotient (UL, UR, P)))
     with Pre => Parameters_Valid (P), Post => T0b'Result;
   function T1b (UL, UR : State; P : Parameters) return Boolean is
     (T0b (UL, UR, P) and then R_Cong (X_Wave_Family (UL, UR, P).Wave2, 0.5 * (Jump (UL, UR) + Quotient (UL, UR, P)))
      and then R (X_Wave_Family (UL, UR, P).Wave2) = R (0.5 * (Jump (UL, UR) + Quotient (UL, UR, P))))
     with Pre => Parameters_Valid (P), Post => T1b'Result;


   function Jump_In_Reals (UL, UR : State; P : Parameters) return Boolean is
     (Half_Exact and then T1 (UL, UR, P) and then T1b (UL, UR, P)
      and then Sub_Err (UR.F, UL.F)
      and then Sub_Err (Jump (UL, UR), Quotient (UL, UR, P))
      and then Add_Err (Jump (UL, UR), Quotient (UL, UR, P))
      and then Mul_Err (0.5, (Jump (UL, UR) - Quotient (UL, UR, P)))
      and then Mul_Err (0.5, (Jump (UL, UR) + Quotient (UL, UR, P)))
      and then Comb (R (UR.F) - R (UL.F), R (Jump (UL, UR)), R (Quotient (UL, UR, P)),
                     R ((Jump (UL, UR) - Quotient (UL, UR, P))), R ((Jump (UL, UR) + Quotient (UL, UR, P))),
                     R (X_Wave_Family (UL, UR, P).Wave1), R (X_Wave_Family (UL, UR, P).Wave2), U, Eta))
     with Pre  => Parameters_Valid (P),
          Post => Jump_In_Reals'Result;


   function Abs_Real (P : Parameters) return Boolean is
     (R (Abs_A (P)) = abs R (P.A) and then R (Abs_A (P)) >= 0.0000009)
     with Pre => Parameters_Valid (P), Post => Abs_Real'Result;

   function Prod_Mem (A : Speed; F : Field) return Boolean is
     (A * F in Flux_Value)
     with Post => Prod_Mem'Result;

   function Diff_Mem (X : Flux_Value; Y : Flux_Value) return Boolean is
     (X - Y in Num_Value)
     with Post => Diff_Mem'Result;

   function FD_Mem (UL, UR : State; P : Parameters) return Boolean is
     (Prod_Mem (P.A, UR.F) and then Prod_Mem (P.A, UL.F)
      and then Diff_Mem (P.A * UR.F, P.A * UL.F)
      and then P.A * UR.F - P.A * UL.F in Num_Value)
     with Post => FD_Mem'Result;

   function FD_Eq (UL, UR : State; P : Parameters) return Boolean is
     (FD_Mem (UL, UR, P) and then Flux_Difference (UL, UR, P) = P.A * UR.F - P.A * UL.F)
     with Post => FD_Eq'Result;

   function Flux_Diff_Real (UL, UR : State; P : Parameters) return Boolean is
     (FD_Eq (UL, UR, P)
      and then R_Cong (Flux_Difference (UL, UR, P), P.A * UR.F - P.A * UL.F)
      and then R (Flux_Difference (UL, UR, P)) = R (P.A * UR.F - P.A * UL.F))
     with Post => Flux_Diff_Real'Result;

   function Jump_Cong (UL, UR : State) return Boolean is
     (UR.F - UL.F in Jump_Value and then Jump (UL, UR) = UR.F - UL.F
      and then R_Cong (Jump (UL, UR), UR.F - UL.F)
      and then R (Jump (UL, UR)) = R (UR.F - UL.F))
     with Post => Jump_Cong'Result;

   function Quotient_Bound (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Eta_Small and then Abs_Real (P) and then Flux_Diff_Real (UL, UR, P)
      and then Mul_Err (P.A, UR.F) and then Mul_Err (P.A, UL.F)
      and then Prod_Mag (R (P.A * UR.F), R (P.A), R (UR.F), Eta)
      and then Prod_Mag (R (P.A * UL.F), R (P.A), R (UL.F), Eta)
      and then Sub_Err (P.A * UR.F, P.A * UL.F)
      and then Num_Mag (R (Flux_Difference (UL, UR, P)), R (P.A * UR.F), R (P.A * UL.F),
                        R (Abs_A (P)), R (UR.F), R (UL.F), Eta)
      and then Div_Err (Flux_Difference (UL, UR, P), Abs_A (P))
      and then Div_Mag (R (Quotient (UL, UR, P)), R (Flux_Difference (UL, UR, P)), R (Abs_A (P)),
                        abs R (UR.F) + abs R (UL.F), Eta))
     with Pre => Parameters_Valid (P), Post => Quotient_Bound'Result;


   --  the jump, as a real: |d - (f_R - f_L)| <= u |f_R - f_L| + eta
   function D_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Jump_Cong (UL, UR) and then Sub_Err (UR.F, UL.F))
     with Post => D_Fact'Result
                  and then abs (R (Jump (UL, UR)) - (R (UR.F) - R (UL.F))) <= Un * abs (R (UR.F) - R (UL.F)) + Eta;

   --  s1 = d - q, as a real
   function S1_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Sub_Err (Jump (UL, UR), Quotient (UL, UR, P)))
     with Pre  => Parameters_Valid (P),
          Post => S1_Fact'Result
                  and then abs (R (Jump (UL, UR) - Quotient (UL, UR, P)) - (R (Jump (UL, UR)) - R (Quotient (UL, UR, P))))
                           <= Un * abs (R (Jump (UL, UR)) - R (Quotient (UL, UR, P))) + Eta;

   --  s2 = d + q, as a real
   function S2_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Add_Err (Jump (UL, UR), Quotient (UL, UR, P)))
     with Pre  => Parameters_Valid (P),
          Post => S2_Fact'Result
                  and then abs (R (Jump (UL, UR) + Quotient (UL, UR, P)) - (R (Jump (UL, UR)) + R (Quotient (UL, UR, P))))
                           <= Un * abs (R (Jump (UL, UR)) + R (Quotient (UL, UR, P))) + Eta;

   --  |q| <= (1+u)^3 (|f_R| + |f_L|) + 5e6 eta
   function Q_Fact (UL, UR : State; P : Parameters) return Boolean is
     (Quotient_Bound (UL, UR, P))
     with Pre  => Parameters_Valid (P),
          Post => Q_Fact'Result
                  and then abs R (Quotient (UL, UR, P))
                           <= (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * (abs R (UR.F) + abs R (UL.F)) + 5000000.0 * Eta;

   --  the sum of the waves, against the jump, in terms of d, q, s1, s2
   function W_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Jump_In_Reals (UL, UR, P))
     with Pre  => Parameters_Valid (P),
          Post => W_Fact'Result
                  and then abs (R (X_Wave_Family (UL, UR, P).Wave1) + R (X_Wave_Family (UL, UR, P).Wave2)
                                - (R (UR.F) - R (UL.F)))
                           <= Un * abs (R (UR.F) - R (UL.F))
                              + Un * (abs (R (Jump (UL, UR)) - R (Quotient (UL, UR, P)))
                                      + abs (R (Jump (UL, UR)) + R (Quotient (UL, UR, P))))
                              + Un * (abs R (Jump (UL, UR) - Quotient (UL, UR, P))
                                      + abs R (Jump (UL, UR) + Quotient (UL, UR, P)))
                              + 5.0 * Eta;

   --  THEOREM, the Lean xWaveJumpCondition in floats: the two computed waves sum to the exact jump within
   --  6u (|f_L| + |f_R|) + 6 eta, where u = 2**-53 and eta = 2**-1074
   function Wave_Jump_Bounded (UL, UR : State; P : Parameters) return Boolean is
     (D_Fact (UL, UR) and then S1_Fact (UL, UR, P) and then S2_Fact (UL, UR, P)
      and then Q_Fact (UL, UR, P) and then W_Fact (UL, UR, P) and then Eta_Small
      and then Jump_Final (R (X_Wave_Family (UL, UR, P).Wave1) + R (X_Wave_Family (UL, UR, P).Wave2),
                           R (UR.F) - R (UL.F), R (Jump (UL, UR)), R (Quotient (UL, UR, P)),
                           R (Jump (UL, UR) - Quotient (UL, UR, P)), R (Jump (UL, UR) + Quotient (UL, UR, P)),
                           abs R (UL.F) + abs R (UR.F), Eta))
     with Pre  => Parameters_Valid (P),
          Post => Wave_Jump_Bounded'Result
                  and then abs (R (X_Wave_Family (UL, UR, P).Wave1) + R (X_Wave_Family (UL, UR, P).Wave2)
                                - (R (UR.F) - R (UL.F)))
                           <= 6.0 * Un * (abs R (UL.F) + abs R (UR.F)) + 6.0 * Eta;

end Lanyon_Adv1d_Jump_Bound_Pkg;
