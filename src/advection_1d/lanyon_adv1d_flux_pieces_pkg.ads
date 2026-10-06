with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Adv1d_Rounding_Pkg; use Lanyon_Adv1d_Rounding_Pkg;
with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
with Lanyon_Adv1d_Waves_Pkg; use Lanyon_Adv1d_Waves_Pkg;
with Lanyon_Adv1d_Fluct_Pkg; use Lanyon_Adv1d_Fluct_Pkg;
with Lanyon_Adv1d_Flux_Bound_Pkg; use Lanyon_Adv1d_Flux_Bound_Pkg;
package Lanyon_Adv1d_Flux_Pieces_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  THE LEAN xFluxConservative IN FLOATS, part 2 of 3: the fluctuations as exact reals, and the four pieces.
   --  the speed parts: min (-|a|, 0) = -|a|, min (|a|, 0) = 0, max (-|a|, 0) = 0, max (|a|, 0) = |a|
   function Parts (P : Parameters) return Boolean is
     (A_Range (P)
      and then Neg_Part (X_Speed_Family (P).Speed1) = -Abs_A (P)
      and then Neg_Part (X_Speed_Family (P).Speed2) = 0.0
      and then Pos_Part (X_Speed_Family (P).Speed1) = 0.0
      and then Pos_Part (X_Speed_Family (P).Speed2) = Abs_A (P)
      and then R_Cong (Neg_Part (X_Speed_Family (P).Speed1), -Abs_A (P))
      and then R_Cong (Neg_Part (X_Speed_Family (P).Speed2), 0.0)
      and then R_Cong (Pos_Part (X_Speed_Family (P).Speed1), 0.0)
      and then R_Cong (Pos_Part (X_Speed_Family (P).Speed2), Abs_A (P))
      and then R (-Abs_A (P)) = 0.0 - R (Abs_A (P)) and then R (0.0) = 0.0)
     with Pre  => Parameters_Valid (P),
          Post => Parts'Result
                  and then R (Neg_Part (X_Speed_Family (P).Speed1)) = 0.0 - A_R (P)
                  and then R (Neg_Part (X_Speed_Family (P).Speed2)) = 0.0
                  and then R (Pos_Part (X_Speed_Family (P).Speed1)) = 0.0
                  and then R (Pos_Part (X_Speed_Family (P).Speed2)) = A_R (P);


   --  a Term as the float product it is
   function Term_Cong (W : Wave_Value; S : Speed) return Boolean is
     (W * S in Term_Value and then Term (W, S) = W * S and then R_Cong (Term (W, S), W * S)
      and then R (Term (W, S)) = R (W * S))
     with Post => Term_Cong'Result;

   --  a fluctuation as the float sum it is
   function Sum_Cong (F : Fluct_Value; X : Term_Value; Y : Term_Value) return Boolean is
     (if F = X + Y then R_Cong (F, X + Y) and then R (F) = R (X + Y))
     with Post => Sum_Cong'Result;

   function NP1 (P : Parameters) return Speed is (Neg_Part (X_Speed_Family (P).Speed1));
   function NP2 (P : Parameters) return Speed is (Neg_Part (X_Speed_Family (P).Speed2));
   function PP1 (P : Parameters) return Speed is (Pos_Part (X_Speed_Family (P).Speed1));
   function PP2 (P : Parameters) return Speed is (Pos_Part (X_Speed_Family (P).Speed2));

   function TL1 (UL, UR : State; P : Parameters) return Long_Float is (Term (W1 (UL, UR, P), NP1 (P)))
     with Pre => Parameters_Valid (P);
   function TL2 (UL, UR : State; P : Parameters) return Long_Float is (Term (W2 (UL, UR, P), NP2 (P)))
     with Pre => Parameters_Valid (P);
   function TR1 (UL, UR : State; P : Parameters) return Long_Float is (Term (W1 (UL, UR, P), PP1 (P)))
     with Pre => Parameters_Valid (P);
   function TR2 (UL, UR : State; P : Parameters) return Long_Float is (Term (W2 (UL, UR, P), PP2 (P)))
     with Pre => Parameters_Valid (P);

   function NP_Real (P : Parameters) return Boolean is
     (Parts (P))
     with Pre  => Parameters_Valid (P),
          Post => NP_Real'Result and then R (NP1 (P)) = 0.0 - A_R (P) and then R (NP2 (P)) = 0.0
                  and then R (PP1 (P)) = 0.0 and then R (PP2 (P)) = A_R (P);

   function T1_Real (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then NP_Real (P) and then Term_Cong (W1 (UL, UR, P), NP1 (P)) and then Mul_Err (W1 (UL, UR, P), NP1 (P)))
     with Pre  => Parameters_Valid (P),
          Post => T1_Real'Result
                  and then R (TL1 (UL, UR, P)) = R (W1 (UL, UR, P) * NP1 (P))
                  and then abs (R (W1 (UL, UR, P) * NP1 (P)) - R (W1 (UL, UR, P)) * R (NP1 (P)))
                           <= Un * abs (R (W1 (UL, UR, P)) * R (NP1 (P))) + Eta;

   function T4_Real (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then NP_Real (P) and then Term_Cong (W2 (UL, UR, P), PP2 (P)) and then Mul_Err (W2 (UL, UR, P), PP2 (P)))
     with Pre  => Parameters_Valid (P),
          Post => T4_Real'Result
                  and then R (TR2 (UL, UR, P)) = R (W2 (UL, UR, P) * PP2 (P))
                  and then abs (R (W2 (UL, UR, P) * PP2 (P)) - R (W2 (UL, UR, P)) * R (PP2 (P)))
                           <= Un * abs (R (W2 (UL, UR, P)) * R (PP2 (P))) + Eta;

   --  the left fluctuation, as reals: T1 ~ -A w1, |T2| <= eta, Lf ~ T1 + T2
   function LF_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Parts (P) and then A_Range (P) and then NP_Real (P) and then T1_Real (UL, UR, P)
      and then Term_Cong (W1 (UL, UR, P), NP1 (P)) and then Term_Cong (W2 (UL, UR, P), NP2 (P))
      and then Mul_Err (W1 (UL, UR, P), NP1 (P)) and then Mul_Err (W2 (UL, UR, P), NP2 (P))
      and then Sum_Cong (X_Left_Fluctuation (UL, UR, P).Fluctuation_F, TL1 (UL, UR, P), TL2 (UL, UR, P))
      and then X_Left_Fluctuation (UL, UR, P).Fluctuation_F = TL1 (UL, UR, P) + TL2 (UL, UR, P)
      and then Add_Err (TL1 (UL, UR, P), TL2 (UL, UR, P))
      and then Subst (R (TL1 (UL, UR, P)), R (W1 (UL, UR, P) * NP1 (P)), R (W1 (UL, UR, P)), R (NP1 (P)), 0.0 - A_R (P), Eta))
     with Pre  => Parameters_Valid (P),
          Post => LF_Fact'Result
                  and then abs (R (TL1 (UL, UR, P)) - R (W1 (UL, UR, P)) * (0.0 - A_R (P)))
                           <= Un * abs (R (W1 (UL, UR, P)) * (0.0 - A_R (P))) + Eta
                  and then abs R (TL2 (UL, UR, P)) <= Eta
                  and then abs (R (X_Left_Fluctuation (UL, UR, P).Fluctuation_F) - (R (TL1 (UL, UR, P)) + R (TL2 (UL, UR, P))))
                           <= Un * abs (R (TL1 (UL, UR, P)) + R (TL2 (UL, UR, P))) + Eta;

   --  the right fluctuation, as reals: |T3| <= eta, T4 ~ A w2, Rf ~ T3 + T4
   function RF_Fact (UL, UR : State; P : Parameters) return Boolean is
     (U_Is_Un and then Parts (P) and then A_Range (P) and then NP_Real (P) and then T4_Real (UL, UR, P)
      and then Term_Cong (W1 (UL, UR, P), PP1 (P)) and then Term_Cong (W2 (UL, UR, P), PP2 (P))
      and then Mul_Err (W1 (UL, UR, P), PP1 (P)) and then Mul_Err (W2 (UL, UR, P), PP2 (P))
      and then Sum_Cong (X_Right_Fluctuation (UL, UR, P).Fluctuation_F, TR1 (UL, UR, P), TR2 (UL, UR, P))
      and then X_Right_Fluctuation (UL, UR, P).Fluctuation_F = TR1 (UL, UR, P) + TR2 (UL, UR, P)
      and then Add_Err (TR1 (UL, UR, P), TR2 (UL, UR, P))
      and then Subst (R (TR2 (UL, UR, P)), R (W2 (UL, UR, P) * PP2 (P)), R (W2 (UL, UR, P)), R (PP2 (P)), A_R (P), Eta))
     with Pre  => Parameters_Valid (P),
          Post => RF_Fact'Result
                  and then abs R (TR1 (UL, UR, P)) <= Eta
                  and then abs (R (TR2 (UL, UR, P)) - R (W2 (UL, UR, P)) * A_R (P))
                           <= Un * abs (R (W2 (UL, UR, P)) * A_R (P)) + Eta
                  and then abs (R (X_Right_Fluctuation (UL, UR, P).Fluctuation_F) - (R (TR1 (UL, UR, P)) + R (TR2 (UL, UR, P))))
                           <= Un * abs (R (TR1 (UL, UR, P)) + R (TR2 (UL, UR, P))) + Eta;

   function Lfx (UL, UR : State; P : Parameters) return Valid_Big_Real is
     (R (X_Left_Fluctuation (UL, UR, P).Fluctuation_F)) with Pre => Parameters_Valid (P);
   function Rfx (UL, UR : State; P : Parameters) return Valid_Big_Real is
     (R (X_Right_Fluctuation (UL, UR, P).Fluctuation_F)) with Pre => Parameters_Valid (P);

   --  the four pieces
   function Piece_Fact (UL, UR : State; P : Parameters) return Boolean is
     (LF_Fact (UL, UR, P) and then RF_Fact (UL, UR, P) and then A_Range (P)
      and then D_Fx (UL, UR) and then S1_Fx (UL, UR, P) and then S2_Fx (UL, UR, P)
      and then W1_Fact (UL, UR, P) and then W2_Fact (UL, UR, P)
      and then P1_Fact (UR, P) and then P1_Fact (UL, P) and then N_Fact (UL, UR, P) and then Q_Div_Fact (UL, UR, P)
      and then Left_Near (Lfx (UL, UR, P), R (TL1 (UL, UR, P)), R (TL2 (UL, UR, P)), Rw1 (UL, UR, P), A_R (P), Eta)
      and then Right_Near (Rfx (UL, UR, P), R (TR1 (UL, UR, P)), R (TR2 (UL, UR, P)), Rw2 (UL, UR, P), A_R (P), Eta)
      and then Waves_Differ (Rw1 (UL, UR, P), Rw2 (UL, UR, P), Rs1 (UL, UR, P), Rs2 (UL, UR, P),
                             Rd (UL, UR), Rq (UL, UR, P), Eta)
      and then Aq_Step (Rq (UL, UR, P), Rn (UL, UR, P), A_R (P), Eta)
      and then Abs_Prod (R (P.A), A_R (P), R (UR.F)) and then Abs_Prod (R (P.A), A_R (P), R (UL.F))
      and then N_Near (Rn (UL, UR, P), Rp1 (UR, P), Rp1 (UL, P), R (P.A), A_R (P), R (UR.F), R (UL.F), Eta)
      and then Aq_Near (Rq (UL, UR, P), Rn (UL, UR, P), Rp1 (UR, P), Rp1 (UL, P), R (P.A), A_R (P),
                        R (UR.F), R (UL.F), Eta))
     with Pre  => Parameters_Valid (P),
          Post => Piece_Fact'Result
                  and then abs (Lfx (UL, UR, P) + A_R (P) * Rw1 (UL, UR, P)) <= 2.1 * Un * (A_R (P) * abs Rw1 (UL, UR, P)) + 4.0 * Eta
                  and then abs (Rfx (UL, UR, P) - A_R (P) * Rw2 (UL, UR, P)) <= 2.1 * Un * (A_R (P) * abs Rw2 (UL, UR, P)) + 4.0 * Eta
                  and then abs (Rw2 (UL, UR, P) - Rw1 (UL, UR, P) - Rq (UL, UR, P))
                           <= Un * (abs (Rd (UL, UR) - Rq (UL, UR, P)) + abs (Rd (UL, UR) + Rq (UL, UR, P)))
                              + Un * (abs Rs1 (UL, UR, P) + abs Rs2 (UL, UR, P)) + 4.0 * Eta
                  and then abs (A_R (P) * Rq (UL, UR, P) - R (P.A) * (R (UR.F) - R (UL.F)))
                           <= Un * abs Rn (UL, UR, P) + Un * abs (Rp1 (UR, P) - Rp1 (UL, P))
                              + Un * (A_R (P) * abs R (UR.F) + A_R (P) * abs R (UL.F)) + A_R (P) * Eta + 3.0 * Eta;

end Lanyon_Adv1d_Flux_Pieces_Pkg;
