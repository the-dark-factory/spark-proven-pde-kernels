with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Adv1d_Rounding_Pkg; use Lanyon_Adv1d_Rounding_Pkg;
with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
with Lanyon_Adv1d_Flux_Bound_Pkg; use Lanyon_Adv1d_Flux_Bound_Pkg;
with Lanyon_Adv1d_Flux_Pieces_Pkg; use Lanyon_Adv1d_Flux_Pieces_Pkg;
package Lanyon_Adv1d_Flux_Final_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  THE LEAN xFluxConservative IN FLOATS, part 3 of 3: the computed fluctuations sum to the exact flux jump
   --  a (f_R - f_L) within 10u |a| (|f_L| + |f_R|) + 1e4 eta.

   --  THEOREM, the Lean xFluxConservative in floats: the computed left and right fluctuations sum to the exact flux
   --  jump a (f_R - f_L) within 10u |a| (|f_L| + |f_R|) + 1e4 eta, where u = 2**-53 and eta = 2**-1074
   function Flux_Conservative_Bounded (UL, UR : State; P : Parameters) return Boolean is
     (Piece_Fact (UL, UR, P) and then Mag_Fact (UL, UR, P) and then Num_Fact (UL, UR, P) and then A_Range (P)
      and then Eta_Small
      and then Scale (A_R (P), Rw2 (UL, UR, P) - Rw1 (UL, UR, P) - Rq (UL, UR, P),
                      Un * (abs (Rd (UL, UR) - Rq (UL, UR, P)) + abs (Rd (UL, UR) + Rq (UL, UR, P)))
                      + Un * (abs Rs1 (UL, UR, P) + abs Rs2 (UL, UR, P)) + 4.0 * Eta)
      and then Flux_Sum (Lfx (UL, UR, P), Rfx (UL, UR, P), Rw1 (UL, UR, P), Rw2 (UL, UR, P), Rq (UL, UR, P), A_R (P),
                         R (P.A) * (R (UR.F) - R (UL.F)),
                         2.1 * Un * (A_R (P) * abs Rw1 (UL, UR, P)) + 4.0 * Eta,
                         2.1 * Un * (A_R (P) * abs Rw2 (UL, UR, P)) + 4.0 * Eta,
                         Un * (abs (Rd (UL, UR) - Rq (UL, UR, P)) + abs (Rd (UL, UR) + Rq (UL, UR, P)))
                         + Un * (abs Rs1 (UL, UR, P) + abs Rs2 (UL, UR, P)) + 4.0 * Eta,
                         Un * abs Rn (UL, UR, P) + Un * abs (Rp1 (UR, P) - Rp1 (UL, P))
                         + Un * (A_R (P) * abs R (UR.F) + A_R (P) * abs R (UL.F)) + A_R (P) * Eta + 3.0 * Eta)
      and then Flux_Final (abs (Lfx (UL, UR, P) + Rfx (UL, UR, P) - R (P.A) * (R (UR.F) - R (UL.F))),
                           A_R (P), Rw1 (UL, UR, P), Rw2 (UL, UR, P), Rs1 (UL, UR, P), Rs2 (UL, UR, P),
                           Rd (UL, UR), Rq (UL, UR, P), Rn (UL, UR, P), Rp1 (UR, P), Rp1 (UL, P),
                           R (UR.F), R (UL.F), M (UL, UR), Eta, K (UL, UR)))
     with Pre  => Parameters_Valid (P),
          Post => Flux_Conservative_Bounded'Result
                  and then abs (Lfx (UL, UR, P) + Rfx (UL, UR, P) - R (P.A) * (R (UR.F) - R (UL.F)))
                           <= 10.0 * Un * (A_R (P) * M (UL, UR)) + 10000.0 * Eta;
end Lanyon_Adv1d_Flux_Final_Pkg;
