with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
with Lanyon_Rb_Jump_Pkg; use Lanyon_Rb_Jump_Pkg;
with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Quot_Pkg; use Lanyon_Burg1d_Quot_Pkg;
with Lanyon_Burg1d_Waves_Pkg; use Lanyon_Burg1d_Waves_Pkg;
package Lanyon_Burg1d_Jump_Bound_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. THE LEAN xWaveJumpCondition IN FLOATS for Burgers: the computed waves sum to the exact jump within
   --  6u (|f_L| + |f_R|) + 6 eta whenever max (|f_L|, |f_R|) >= 2**-511, and exactly at the sonic point. Below
   --  2**-511 (both states' squares subnormal) the per-operation lemmas give only an absolute bound: the provers lack
   --  the round-to-nearest fact |fl (x) - x| <= |x| (refused at levels 2 and 4).

   Floor_N : constant := 2.0 ** (-511);

   --  the waves away from the sonic point, as the expressions they are
   function T0 (UL, UR : State) return Boolean is
     (0.5 * (Jump (UL, UR) - Quotient (UL, UR)) in Wave_Value
      and then X_Wave_Family (UL, UR).Wave1 = 0.5 * (Jump (UL, UR) - Quotient (UL, UR)))
     with Pre => Amax (UL, UR) > 0.0, Post => T0'Result;

   function T1 (UL, UR : State) return Boolean is
     (T0 (UL, UR) and then R_Cong (X_Wave_Family (UL, UR).Wave1, 0.5 * (Jump (UL, UR) - Quotient (UL, UR)))
      and then R (X_Wave_Family (UL, UR).Wave1) = R (0.5 * (Jump (UL, UR) - Quotient (UL, UR))))
     with Pre => Amax (UL, UR) > 0.0, Post => T1'Result;

   function T0b (UL, UR : State) return Boolean is
     (0.5 * (Jump (UL, UR) + Quotient (UL, UR)) in Wave_Value
      and then X_Wave_Family (UL, UR).Wave2 = 0.5 * (Jump (UL, UR) + Quotient (UL, UR)))
     with Pre => Amax (UL, UR) > 0.0, Post => T0b'Result;

   function T1b (UL, UR : State) return Boolean is
     (T0b (UL, UR) and then R_Cong (X_Wave_Family (UL, UR).Wave2, 0.5 * (Jump (UL, UR) + Quotient (UL, UR)))
      and then R (X_Wave_Family (UL, UR).Wave2) = R (0.5 * (Jump (UL, UR) + Quotient (UL, UR))))
     with Pre => Amax (UL, UR) > 0.0, Post => T1b'Result;

   function Jump_In_Reals (UL, UR : State) return Boolean is
     (Half_Exact and then T1 (UL, UR) and then T1b (UL, UR)
      and then Sub_Err (UR.F, UL.F)
      and then Sub_Err (Jump (UL, UR), Quotient (UL, UR))
      and then Add_Err (Jump (UL, UR), Quotient (UL, UR))
      and then Mul_Err (0.5, (Jump (UL, UR) - Quotient (UL, UR)))
      and then Mul_Err (0.5, (Jump (UL, UR) + Quotient (UL, UR)))
      and then Comb (R (UR.F) - R (UL.F), R (Jump (UL, UR)), R (Quotient (UL, UR)),
                     R ((Jump (UL, UR) - Quotient (UL, UR))), R ((Jump (UL, UR) + Quotient (UL, UR))),
                     R (X_Wave_Family (UL, UR).Wave1), R (X_Wave_Family (UL, UR).Wave2), U, Eta))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => Jump_In_Reals'Result;

   function Jump_Cong (UL, UR : State) return Boolean is
     (UR.F - UL.F in Jump_Value and then Jump (UL, UR) = UR.F - UL.F
      and then R_Cong (Jump (UL, UR), UR.F - UL.F)
      and then R (Jump (UL, UR)) = R (UR.F - UL.F))
     with Post => Jump_Cong'Result;

   --  a nonnegative float has a nonnegative real value
   function R_Nonneg (X : Long_Float) return Boolean is (if X >= 0.0 then R (X) >= 0.0)
     with Post => R_Nonneg'Result;

   --  away from the deep-subnormal regime, eta is negligible against A^2
   function Eta_Vs_A (A : Any_Pos) return Boolean is
     (Eta <= 0.0000000000000003 * (R (A) * R (A)))
     with Pre => A >= Floor_N, Post => Eta_Vs_A'Result;

   --  pure real: the numerator is at most 0.6 A^2 + 4 E (the two halves are nonnegative)
   function Num_Real_B (N, H1, H2, P1, P2, FR, FL, A, E : Valid_Big_Real) return Boolean is
     (abs N <= 0.6 * (A * A) + 4.0 * E)
     with Pre  => E >= 0.0 and then E <= 0.000000000000000000001 and then H1 >= 0.0 and then H2 >= 0.0
                  and then FR * FR <= A * A and then FL * FL <= A * A
                  and then abs (P1 - FR * FR) <= Un * abs (FR * FR) + E
                  and then abs (P2 - FL * FL) <= Un * abs (FL * FL) + E
                  and then abs (H1 - 0.5 * P1) <= Un * abs (0.5 * P1) + E
                  and then abs (H2 - 0.5 * P2) <= Un * abs (0.5 * P2) + E
                  and then abs (N - (H1 - H2)) <= Un * abs (H1 - H2) + E,
          Post => Num_Real_B'Result;

   --  pure real: |N / A| <= 0.6 A + 4 E / A, and E / A <= 3e-16 A
   function NA_B (N, A, E : Valid_Big_Real) return Boolean is
     (abs (N / A) <= 0.6 * A + 4.0 * (E / A) and then E / A <= 0.0000000000000003 * A)
     with Pre  => A > 0.0 and then E >= 0.0 and then abs N <= 0.6 * (A * A) + 4.0 * E
                  and then E <= 0.0000000000000003 * (A * A),
          Post => NA_B'Result;

   --  pure real: |q| <= (1+u)^3 M + 5e6 E, from |N / A| <= 0.61 A and A <= M
   function Q_Rel_B (Q, N, A, M, E : Valid_Big_Real) return Boolean is
     (abs Q <= (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M + 5000000.0 * E)
     with Pre  => A > 0.0 and then E >= 0.0 and then A <= M
                  and then abs (N / A) <= 0.61 * A
                  and then abs (Q - N / A) <= Un * abs (N / A) + E,
          Post => Q_Rel_B'Result;

   --  the divisor is one of |f_L|, |f_R|, so its real value is at most |f_L| + |f_R|
   function Max_Le_Sum (UL, UR : State) return Boolean is
     (R_Abs (UL.F) and then R_Abs (UR.F)
      and then (if abs UL.F >= abs UR.F
                then Amax (UL, UR) = abs UL.F and then R_Cong (Amax (UL, UR), abs UL.F)
                else Amax (UL, UR) = abs UR.F and then R_Cong (Amax (UL, UR), abs UR.F)))
     with Post => Max_Le_Sum'Result and then R (Amax (UL, UR)) <= abs R (UR.F) + abs R (UL.F);

   function Q_Fact (UL, UR : State) return Boolean is
     (Quot_Safe (UL.F, UR.F) and then Eta_Le (Amax (UL, UR))
      and then U_Is_Un and then Half_Exact and then Eta_Small and then Max_Real (UL.F, UR.F)
      and then Eta_Vs_A (Amax (UL, UR)) and then Max_Le_Sum (UL, UR)
      and then Sq_Le (R (UR.F), R (Amax (UL, UR))) and then Sq_Le (R (UL.F), R (Amax (UL, UR)))
      and then Mul_Err (UR.F, UR.F) and then Mul_Err (UL.F, UL.F)
      and then Mul_Err (0.5, UR.F * UR.F) and then Mul_Err (0.5, UL.F * UL.F)
      and then Sub_Err (0.5 * (UR.F * UR.F), 0.5 * (UL.F * UL.F))
      and then R_Nonneg (0.5 * (UR.F * UR.F)) and then R_Nonneg (0.5 * (UL.F * UL.F))
      and then Num_Real_B (R (Num (UL.F, UR.F)), R (0.5 * (UR.F * UR.F)), R (0.5 * (UL.F * UL.F)),
                           R (UR.F * UR.F), R (UL.F * UL.F), R (UR.F), R (UL.F), R (Amax (UL, UR)), Eta)
      and then NA_B (R (Num (UL.F, UR.F)), R (Amax (UL, UR)), Eta)
      and then Div_Any (Num (UL.F, UR.F), Amax (UL, UR))
      and then Q_Rel_B (R (Quotient (UL, UR)), R (Num (UL.F, UR.F)), R (Amax (UL, UR)),
                        abs R (UR.F) + abs R (UL.F), Eta))
     with Pre  => Amax (UL, UR) >= Floor_N,
          Post => Q_Fact'Result
                  and then abs R (Quotient (UL, UR))
                           <= (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * (abs R (UR.F) + abs R (UL.F)) + 5000000.0 * Eta;

   function D_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Jump_Cong (UL, UR) and then Sub_Err (UR.F, UL.F))
     with Post => D_Fact'Result
                  and then abs (R (Jump (UL, UR)) - (R (UR.F) - R (UL.F))) <= Un * abs (R (UR.F) - R (UL.F)) + Eta;

   function S1_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Sub_Err (Jump (UL, UR), Quotient (UL, UR)))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => S1_Fact'Result
                  and then abs (R (Jump (UL, UR) - Quotient (UL, UR)) - (R (Jump (UL, UR)) - R (Quotient (UL, UR))))
                           <= Un * abs (R (Jump (UL, UR)) - R (Quotient (UL, UR))) + Eta;

   function S2_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Add_Err (Jump (UL, UR), Quotient (UL, UR)))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => S2_Fact'Result
                  and then abs (R (Jump (UL, UR) + Quotient (UL, UR)) - (R (Jump (UL, UR)) + R (Quotient (UL, UR))))
                           <= Un * abs (R (Jump (UL, UR)) + R (Quotient (UL, UR))) + Eta;

   function W_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Jump_In_Reals (UL, UR))
     with Pre  => Amax (UL, UR) > 0.0,
          Post => W_Fact'Result
                  and then abs (R (X_Wave_Family (UL, UR).Wave1) + R (X_Wave_Family (UL, UR).Wave2)
                                - (R (UR.F) - R (UL.F)))
                           <= Un * abs (R (UR.F) - R (UL.F))
                              + Un * (abs (R (Jump (UL, UR)) - R (Quotient (UL, UR)))
                                      + abs (R (Jump (UL, UR)) + R (Quotient (UL, UR))))
                              + Un * (abs R (Jump (UL, UR) - Quotient (UL, UR))
                                      + abs R (Jump (UL, UR) + Quotient (UL, UR)))
                              + 5.0 * Eta;

   --  THEOREM, the Lean xWaveJumpCondition in floats for Burgers, away from the deep-subnormal regime
   function Wave_Jump_Bounded (UL, UR : State) return Boolean is
     (Quot_Safe (UL.F, UR.F) and then D_Fact (UL, UR) and then S1_Fact (UL, UR) and then S2_Fact (UL, UR)
      and then Q_Fact (UL, UR) and then W_Fact (UL, UR) and then Eta_Small
      and then Jump_Final (R (X_Wave_Family (UL, UR).Wave1) + R (X_Wave_Family (UL, UR).Wave2),
                           R (UR.F) - R (UL.F), R (Jump (UL, UR)), R (Quotient (UL, UR)),
                           R (Jump (UL, UR) - Quotient (UL, UR)), R (Jump (UL, UR) + Quotient (UL, UR)),
                           abs R (UL.F) + abs R (UR.F), Eta))
     with Pre  => Amax (UL, UR) >= Floor_N,
          Post => Wave_Jump_Bounded'Result
                  and then abs (R (X_Wave_Family (UL, UR).Wave1) + R (X_Wave_Family (UL, UR).Wave2)
                                - (R (UR.F) - R (UL.F)))
                           <= 6.0 * Un * (abs R (UL.F) + abs R (UR.F)) + 6.0 * Eta;

   --  THEOREM, at the sonic point the jump condition holds exactly
   function Wave_Jump_Sonic (UL, UR : State) return Boolean is
     (Sonic_Defined (UL, UR)
      and then X_Wave_Family (UL, UR).Wave1 + X_Wave_Family (UL, UR).Wave2 = Jump (UL, UR))
     with Pre  => UL.F = 0.0 and then UR.F = 0.0,
          Post => Wave_Jump_Sonic'Result;

end Lanyon_Burg1d_Jump_Bound_Pkg;
