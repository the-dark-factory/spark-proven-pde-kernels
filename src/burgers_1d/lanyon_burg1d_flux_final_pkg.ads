with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
with Lanyon_Rb_Flux_Pkg; use Lanyon_Rb_Flux_Pkg;
with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Flux_Facts_Pkg; use Lanyon_Burg1d_Flux_Facts_Pkg;
with Lanyon_Burg1d_Flux_Pieces_Pkg; use Lanyon_Burg1d_Flux_Pieces_Pkg;
package Lanyon_Burg1d_Flux_Final_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081.
   --  THE LEAN xFluxConservative IN FLOATS for Burgers, part 3 of 3: the theorem.

   --  pure real: A q is the exact flux jump, to rounding
   function Aq_Near_B (AQ, N, J, Bn, Ba : Valid_Big_Real) return Boolean is
     (abs (AQ - J) <= Ba + Bn)
     with Pre  => abs (AQ - N) <= Ba and then abs (N - J) <= Bn,
          Post => Aq_Near_B'Result;

   --  pure real: the final collapse for Burgers
   function Flux_Final_B (Err, A, W1, W2, S1, S2, D, Q, N, H1, H2, P1, P2, FR, FL, M, E, K : Valid_Big_Real)
     return Boolean is
     (Err <= 10.0 * Un * (A * M) + 10000000.0 * E)
     with Pre  => E >= 0.0 and then E <= 0.000000000000000000001 and then M >= 0.0
                  and then A > 0.0 and then A <= 1000000.0 and then A <= M
                  and then FR * FR <= A * A and then FL * FL <= A * A
                  and then K = (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M + 5000000.0 * E
                  and then abs W1 + abs W2 <= 1.001 * K + 4.0 * E
                  and then abs (D - Q) + abs (D + Q) <= 2.0 * K
                  and then abs S1 + abs S2 <= 2.001 * K + 2.0 * E
                  and then abs N <= 0.6 * (A * A) + 4.0 * E
                  and then abs (H1 - H2) <= 0.51 * (A * A) + 2.0 * E
                  and then abs P1 + abs P2 <= 2.01 * (A * A) + 2.0 * E
                  and then Err <= (2.1 * Un * (A * abs W1) + 4.0 * E) + (2.1 * Un * (A * abs W2) + 4.0 * E)
                                  + A * (Un * (abs (D - Q) + abs (D + Q)) + Un * (abs S1 + abs S2) + 4.0 * E)
                                  + ((Un * abs N + A * E)
                                     + (Un * abs (H1 - H2) + 0.5 * Un * (abs P1 + abs P2)
                                        + 0.5 * Un * (FR * FR + FL * FL) + 4.0 * E)),
          Post => Flux_Final_B'Result;

   --  THEOREM, the Lean xFluxConservative in floats for Burgers: the computed fluctuations sum to the exact flux
   --  jump f_R^2/2 - f_L^2/2 within 10u A (|f_L| + |f_R|) + 1e7 eta, A = max (|f_L|, |f_R|) >= 2**-511
   function Flux_Conservative_Bounded (UL, UR : State) return Boolean is
     (Piece_Fact (UL, UR) and then Mag_Fact (UL, UR) and then Chain_Fact (UL, UR)
      and then A_Fact (UL, UR) and then Eta_Small
      and then Aq_Near_B (A_R (UL, UR) * Rq (UL, UR), Rn (UL, UR), Jx (UL, UR),
                          Un * abs (Rh (UR.F) - Rh (UL.F)) + 0.5 * Un * (abs Rp (UR.F) + abs Rp (UL.F))
                          + 0.5 * Un * (R (UR.F) * R (UR.F) + R (UL.F) * R (UL.F)) + 4.0 * Eta,
                          Un * abs Rn (UL, UR) + A_R (UL, UR) * Eta)
      and then Scale (A_R (UL, UR), R (W2 (UL, UR)) - R (W1 (UL, UR)) - Rq (UL, UR),
                      Un * (abs (Rd (UL, UR) - Rq (UL, UR)) + abs (Rd (UL, UR) + Rq (UL, UR)))
                      + Un * (abs Rs1 (UL, UR) + abs Rs2 (UL, UR)) + 4.0 * Eta)
      and then Flux_Sum (Lfx (UL, UR), Rfx (UL, UR), R (W1 (UL, UR)), R (W2 (UL, UR)), Rq (UL, UR), A_R (UL, UR),
                         Jx (UL, UR),
                         2.1 * Un * (A_R (UL, UR) * abs R (W1 (UL, UR))) + 4.0 * Eta,
                         2.1 * Un * (A_R (UL, UR) * abs R (W2 (UL, UR))) + 4.0 * Eta,
                         Un * (abs (Rd (UL, UR) - Rq (UL, UR)) + abs (Rd (UL, UR) + Rq (UL, UR)))
                         + Un * (abs Rs1 (UL, UR) + abs Rs2 (UL, UR)) + 4.0 * Eta,
                         (Un * abs Rn (UL, UR) + A_R (UL, UR) * Eta)
                         + (Un * abs (Rh (UR.F) - Rh (UL.F)) + 0.5 * Un * (abs Rp (UR.F) + abs Rp (UL.F))
                            + 0.5 * Un * (R (UR.F) * R (UR.F) + R (UL.F) * R (UL.F)) + 4.0 * Eta))
      and then Flux_Final_B (abs (Lfx (UL, UR) + Rfx (UL, UR) - Jx (UL, UR)), A_R (UL, UR),
                             R (W1 (UL, UR)), R (W2 (UL, UR)), Rs1 (UL, UR), Rs2 (UL, UR), Rd (UL, UR), Rq (UL, UR),
                             Rn (UL, UR), Rh (UR.F), Rh (UL.F), Rp (UR.F), Rp (UL.F), R (UR.F), R (UL.F),
                             M (UL, UR), Eta, K (UL, UR)))
     with Pre  => Amax (UL, UR) >= Floor_B,
          Post => Flux_Conservative_Bounded'Result
                  and then abs (Lfx (UL, UR) + Rfx (UL, UR) - Jx (UL, UR))
                           <= 10.0 * Un * (A_R (UL, UR) * M (UL, UR)) + 10000000.0 * Eta;

end Lanyon_Burg1d_Flux_Final_Pkg;
