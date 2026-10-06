with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
package Lanyon_Rb_Jump_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The pure real lemmas of the bounded jump condition.

   --  the jump theorem collapsed to the Lean shape
   function Jump_Final (W, J, D, Q, S1, S2, M, E : Valid_Big_Real) return Boolean is
     (abs (W - J) <= 6.0 * Un * M + 6.0 * E)
     with Pre  => E >= 0.0 and then M >= 0.0 and then E <= 0.000000000000000000001
                  and then abs J <= M
                  and then abs (D - J) <= Un * abs J + E
                  and then abs (S1 - (D - Q)) <= Un * abs (D - Q) + E
                  and then abs (S2 - (D + Q)) <= Un * abs (D + Q) + E
                  and then abs Q <= (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M + 5000000.0 * E
                  and then abs (W - J) <= Un * abs J + Un * (abs (D - Q) + abs (D + Q)) + Un * (abs S1 + abs S2) + 5.0 * E,
          Post => Jump_Final'Result;

   --  the pure real step: given the five per-operation error facts, the jump is reproduced within the sum
   function Comb (J, D, Q, S1, S2, W1, W2, Uv, Ev : Valid_Big_Real) return Boolean is
     (abs (W1 + W2 - J) <= Uv * abs J + Uv * (abs (D - Q) + abs (D + Q)) + Uv * (abs S1 + abs S2) + 5.0 * Ev)
     with Pre  => Uv >= 0.0 and then Ev >= 0.0
                  and then abs (D - J) <= Uv * abs J + Ev
                  and then abs (S1 - (D - Q)) <= Uv * abs (D - Q) + Ev
                  and then abs (S2 - (D + Q)) <= Uv * abs (D + Q) + Ev
                  and then abs (W1 - 0.5 * S1) <= Uv * abs (0.5 * S1) + Ev
                  and then abs (W2 - 0.5 * S2) <= Uv * abs (0.5 * S2) + Ev,
          Post => Comb'Result;

end Lanyon_Rb_Jump_Pkg;
