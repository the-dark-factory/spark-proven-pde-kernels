with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
with Lanyon_Rb_Flux_Pkg; use Lanyon_Rb_Flux_Pkg;
with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Waves_Pkg; use Lanyon_Burg1d_Waves_Pkg;
with Lanyon_Burg1d_Fluct_Pkg; use Lanyon_Burg1d_Fluct_Pkg;
with Lanyon_Burg1d_Flux_Facts_Pkg; use Lanyon_Burg1d_Flux_Facts_Pkg;
package Lanyon_Burg1d_Flux_Pieces_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081.
   --  THE LEAN xFluxConservative IN FLOATS for Burgers, part 2 of 3: the fluctuations and the four pieces.

   function Lfx (UL, UR : State) return Valid_Big_Real is (R (X_Left_Fluctuation (UL, UR).Fluctuation_F));

   function Rfx (UL, UR : State) return Valid_Big_Real is (R (X_Right_Fluctuation (UL, UR).Fluctuation_F));

   --  the speed parts: min (-A, 0) = -A, min (A, 0) = 0, max (-A, 0) = 0, max (A, 0) = A
   function NP1 (UL, UR : State) return Field is (Neg_Part (X_Speed_Family (UL, UR).Speed1));

   function NP2 (UL, UR : State) return Field is (Neg_Part (X_Speed_Family (UL, UR).Speed2));

   function PP1 (UL, UR : State) return Field is (Pos_Part (X_Speed_Family (UL, UR).Speed1));

   function PP2 (UL, UR : State) return Field is (Pos_Part (X_Speed_Family (UL, UR).Speed2));

   function Parts (UL, UR : State) return Boolean is
     (NP1 (UL, UR) = -Amax (UL, UR) and then NP2 (UL, UR) = 0.0
      and then PP1 (UL, UR) = 0.0 and then PP2 (UL, UR) = Amax (UL, UR)
      and then R_Cong (NP1 (UL, UR), -Amax (UL, UR)) and then R_Cong (NP2 (UL, UR), 0.0)
      and then R_Cong (PP1 (UL, UR), 0.0) and then R_Cong (PP2 (UL, UR), Amax (UL, UR))
      and then R (-Amax (UL, UR)) = 0.0 - R (Amax (UL, UR)) and then R (0.0) = 0.0)
     with Post => Parts'Result
                  and then R (NP1 (UL, UR)) = 0.0 - A_R (UL, UR) and then R (NP2 (UL, UR)) = 0.0
                  and then R (PP1 (UL, UR)) = 0.0 and then R (PP2 (UL, UR)) = A_R (UL, UR);

   --  a Term as the float product it is; a fluctuation as the float sum it is
   function Term_Cong (W : Wave_Value; S : Field) return Boolean is
     (W * S in Term_Value and then Term (W, S) = W * S and then R_Cong (Term (W, S), W * S)
      and then R (Term (W, S)) = R (W * S))
     with Post => Term_Cong'Result;

   function Sum_Cong (F : Fluct_Value; X : Term_Value; Y : Term_Value) return Boolean is
     (if F = X + Y then R_Cong (F, X + Y) and then R (F) = R (X + Y))
     with Post => Sum_Cong'Result;

   function TL1 (UL, UR : State) return Long_Float is (Term (W1 (UL, UR), NP1 (UL, UR)));

   function TL2 (UL, UR : State) return Long_Float is (Term (W2 (UL, UR), NP2 (UL, UR)));

   function TR1 (UL, UR : State) return Long_Float is (Term (W1 (UL, UR), PP1 (UL, UR)));

   function TR2 (UL, UR : State) return Long_Float is (Term (W2 (UL, UR), PP2 (UL, UR)));

   --  a wave times a speed part, as a real
   function WS_Fact (W : Wave_Value; S : Field) return Boolean is
     (U_Is_Un and then Mul_Err (W, S))
     with Post => WS_Fact'Result and then abs (R (W * S) - R (W) * R (S)) <= Un * abs (R (W) * R (S)) + Eta;

   function T1_Real (UL, UR : State) return Boolean is
     (U_Is_Un and then Parts (UL, UR) and then Term_Cong (W1 (UL, UR), NP1 (UL, UR))
      and then WS_Fact (W1 (UL, UR), NP1 (UL, UR)))
     with Post => T1_Real'Result
                  and then R (TL1 (UL, UR)) = R (W1 (UL, UR) * NP1 (UL, UR))
                  and then abs (R (W1 (UL, UR) * NP1 (UL, UR)) - R (W1 (UL, UR)) * R (NP1 (UL, UR)))
                           <= Un * abs (R (W1 (UL, UR)) * R (NP1 (UL, UR))) + Eta;

   function T4_Real (UL, UR : State) return Boolean is
     (U_Is_Un and then Parts (UL, UR) and then Term_Cong (W2 (UL, UR), PP2 (UL, UR))
      and then WS_Fact (W2 (UL, UR), PP2 (UL, UR)))
     with Post => T4_Real'Result
                  and then R (TR2 (UL, UR)) = R (W2 (UL, UR) * PP2 (UL, UR))
                  and then abs (R (W2 (UL, UR) * PP2 (UL, UR)) - R (W2 (UL, UR)) * R (PP2 (UL, UR)))
                           <= Un * abs (R (W2 (UL, UR)) * R (PP2 (UL, UR))) + Eta;

   --  the left fluctuation: T1 ~ -A w1, |T2| <= eta, Lf ~ T1 + T2
   function LF_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Parts (UL, UR) and then T1_Real (UL, UR)
      and then Term_Cong (W2 (UL, UR), NP2 (UL, UR)) and then Mul_Err (W2 (UL, UR), NP2 (UL, UR))
      and then Sum_Cong (X_Left_Fluctuation (UL, UR).Fluctuation_F, TL1 (UL, UR), TL2 (UL, UR))
      and then X_Left_Fluctuation (UL, UR).Fluctuation_F = TL1 (UL, UR) + TL2 (UL, UR)
      and then Add_Err (TL1 (UL, UR), TL2 (UL, UR))
      and then Subst (R (TL1 (UL, UR)), R (W1 (UL, UR) * NP1 (UL, UR)), R (W1 (UL, UR)), R (NP1 (UL, UR)),
                      0.0 - A_R (UL, UR), Eta))
     with Post => LF_Fact'Result
                  and then abs (R (TL1 (UL, UR)) - R (W1 (UL, UR)) * (0.0 - A_R (UL, UR)))
                           <= Un * abs (R (W1 (UL, UR)) * (0.0 - A_R (UL, UR))) + Eta
                  and then abs R (TL2 (UL, UR)) <= Eta
                  and then abs (Lfx (UL, UR) - (R (TL1 (UL, UR)) + R (TL2 (UL, UR))))
                           <= Un * abs (R (TL1 (UL, UR)) + R (TL2 (UL, UR))) + Eta;

   --  the right fluctuation: |T3| <= eta, T4 ~ A w2, Rf ~ T3 + T4
   function RF_Fact (UL, UR : State) return Boolean is
     (U_Is_Un and then Parts (UL, UR) and then T4_Real (UL, UR)
      and then Term_Cong (W1 (UL, UR), PP1 (UL, UR)) and then Mul_Err (W1 (UL, UR), PP1 (UL, UR))
      and then Sum_Cong (X_Right_Fluctuation (UL, UR).Fluctuation_F, TR1 (UL, UR), TR2 (UL, UR))
      and then X_Right_Fluctuation (UL, UR).Fluctuation_F = TR1 (UL, UR) + TR2 (UL, UR)
      and then Add_Err (TR1 (UL, UR), TR2 (UL, UR))
      and then Subst (R (TR2 (UL, UR)), R (W2 (UL, UR) * PP2 (UL, UR)), R (W2 (UL, UR)), R (PP2 (UL, UR)),
                      A_R (UL, UR), Eta))
     with Post => RF_Fact'Result
                  and then abs R (TR1 (UL, UR)) <= Eta
                  and then abs (R (TR2 (UL, UR)) - R (W2 (UL, UR)) * A_R (UL, UR))
                           <= Un * abs (R (W2 (UL, UR)) * A_R (UL, UR)) + Eta
                  and then abs (Rfx (UL, UR) - (R (TR1 (UL, UR)) + R (TR2 (UL, UR))))
                           <= Un * abs (R (TR1 (UL, UR)) + R (TR2 (UL, UR))) + Eta;

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. The assembly of the Burgers flux theorem from the pieces in Lanyon_Burg1d_Flux_Bound_Pkg.
   --  pure real: multiplying the quotient error by any A > 0
   function Aq_Step_B (Q, N, A, E : Valid_Big_Real) return Boolean is
     (abs (A * Q - N) <= Un * abs N + A * E)
     with Pre  => A > 0.0 and then E >= 0.0 and then abs (Q - N / A) <= Un * abs (N / A) + E,
          Post => Aq_Step_B'Result;

   --  the four pieces
   function Piece_Fact (UL, UR : State) return Boolean is
     (A_Fact (UL, UR)
      and then LF_Fact (UL, UR) and then RF_Fact (UL, UR)
      and then S1_Fx (UL, UR) and then S2_Fx (UL, UR) and then W1_Fact (UL, UR) and then W2_Fact (UL, UR)
      and then Chain_Fact (UL, UR)
      and then QD_Fx (UL, UR)
      and then Left_Near (Lfx (UL, UR), R (TL1 (UL, UR)), R (TL2 (UL, UR)), R (W1 (UL, UR)), A_R (UL, UR), Eta)
      and then Right_Near (Rfx (UL, UR), R (TR1 (UL, UR)), R (TR2 (UL, UR)), R (W2 (UL, UR)), A_R (UL, UR), Eta)
      and then Waves_Differ (R (W1 (UL, UR)), R (W2 (UL, UR)), Rs1 (UL, UR), Rs2 (UL, UR),
                             Rd (UL, UR), Rq (UL, UR), Eta)
      and then Aq_Step_B (Rq (UL, UR), Rn (UL, UR), A_R (UL, UR), Eta))
     with Pre  => Amax (UL, UR) >= Floor_B,
          Post => Piece_Fact'Result
                  and then abs (Lfx (UL, UR) + A_R (UL, UR) * R (W1 (UL, UR)))
                           <= 2.1 * Un * (A_R (UL, UR) * abs R (W1 (UL, UR))) + 4.0 * Eta
                  and then abs (Rfx (UL, UR) - A_R (UL, UR) * R (W2 (UL, UR)))
                           <= 2.1 * Un * (A_R (UL, UR) * abs R (W2 (UL, UR))) + 4.0 * Eta
                  and then abs (R (W2 (UL, UR)) - R (W1 (UL, UR)) - Rq (UL, UR))
                           <= Un * (abs (Rd (UL, UR) - Rq (UL, UR)) + abs (Rd (UL, UR) + Rq (UL, UR)))
                              + Un * (abs Rs1 (UL, UR) + abs Rs2 (UL, UR)) + 4.0 * Eta
                  and then abs (A_R (UL, UR) * Rq (UL, UR) - Rn (UL, UR)) <= Un * abs Rn (UL, UR) + A_R (UL, UR) * Eta;

end Lanyon_Burg1d_Flux_Pieces_Pkg;
