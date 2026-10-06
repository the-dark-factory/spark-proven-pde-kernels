with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
package Lanyon_Rb_Flux_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The pure real lemmas of the bounded flux conservation.

   --  (a) the left fluctuation is -A w1, to rounding: Lf = fl (T1 + T2), T1 = fl (w1 * (-A)), T2 = fl (w2 * 0)
   function Left_Near (Lf, T1, T2, W1, A, E : Valid_Big_Real) return Boolean is
     (abs (Lf + A * W1) <= 2.1 * Un * (A * abs W1) + 4.0 * E)
     with Pre  => A >= 0.0 and then E >= 0.0
                  and then abs (T1 - W1 * (0.0 - A)) <= Un * abs (W1 * (0.0 - A)) + E
                  and then abs T2 <= E
                  and then abs (Lf - (T1 + T2)) <= Un * abs (T1 + T2) + E,
          Post => Left_Near'Result;

   --  (b) the right fluctuation is A w2, to rounding: Rf = fl (T3 + T4), T3 = fl (w1 * 0), T4 = fl (w2 * A)
   function Right_Near (Rf, T3, T4, W2, A, E : Valid_Big_Real) return Boolean is
     (abs (Rf - A * W2) <= 2.1 * Un * (A * abs W2) + 4.0 * E)
     with Pre  => A >= 0.0 and then E >= 0.0
                  and then abs T3 <= E
                  and then abs (T4 - W2 * A) <= Un * abs (W2 * A) + E
                  and then abs (Rf - (T3 + T4)) <= Un * abs (T3 + T4) + E,
          Post => Right_Near'Result;

   --  (c) the two waves differ by q, to rounding
   function Waves_Differ (W1, W2, S1, S2, D, Q, E : Valid_Big_Real) return Boolean is
     (abs (W2 - W1 - Q) <= Un * (abs (D - Q) + abs (D + Q)) + Un * (abs S1 + abs S2) + 4.0 * E)
     with Pre  => E >= 0.0
                  and then abs (S1 - (D - Q)) <= Un * abs (D - Q) + E
                  and then abs (S2 - (D + Q)) <= Un * abs (D + Q) + E
                  and then abs (W1 - 0.5 * S1) <= Un * abs (0.5 * S1) + E
                  and then abs (W2 - 0.5 * S2) <= Un * abs (0.5 * S2) + E,
          Post => Waves_Differ'Result;

   --  scaling a bound by A >= 0
   function Scale (A, X, B : Valid_Big_Real) return Boolean is
     (A * abs X <= A * B and then abs (A * X) = A * abs X)
     with Pre  => A >= 0.0 and then abs X <= B,
          Post => Scale'Result;

   --  (e) the four pieces assembled: Lf + Rf - J = (Lf + A w1) + (Rf - A w2) + A (w2 - w1 - q) + (A q - J)
   function Flux_Sum (Lf, Rf, W1, W2, Q, A, J, Ba, Bb, Bc, Bd : Valid_Big_Real) return Boolean is
     (abs (Lf + Rf - J) <= Ba + Bb + A * Bc + Bd)
     with Pre  => A >= 0.0
                  and then abs (Lf + A * W1) <= Ba
                  and then abs (Rf - A * W2) <= Bb
                  and then A * abs (W2 - W1 - Q) <= A * Bc
                  and then abs (A * Q - J) <= Bd,
          Post => Flux_Sum'Result;

   --  (f) the magnitudes, in terms of M = |fL| + |fR| and K = (1+u)^3 M + 5e6 E
   function Mags (W1, W2, S1, S2, D, Q, M, E, K : Valid_Big_Real) return Boolean is
     (abs W1 + abs W2 <= 1.001 * K + 4.0 * E
      and then abs (D - Q) + abs (D + Q) <= 2.0 * K
      and then abs S1 + abs S2 <= 2.001 * K + 2.0 * E)
     with Pre  => E >= 0.0 and then M >= 0.0
                  and then K = (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M + 5000000.0 * E
                  and then abs D <= K and then abs Q <= K
                  and then abs (S1 - (D - Q)) <= Un * abs (D - Q) + E
                  and then abs (S2 - (D + Q)) <= Un * abs (D + Q) + E
                  and then abs (W1 - 0.5 * S1) <= Un * abs (0.5 * S1) + E
                  and then abs (W2 - 0.5 * S2) <= Un * abs (0.5 * S2) + E,
          Post => Mags'Result;

   --  substituting X = Y under the product error bound
   function Subst (T, Pr, W, X, Y, E : Valid_Big_Real) return Boolean is
     (abs (T - W * Y) <= Un * abs (W * Y) + E)
     with Pre  => X = Y and then T = Pr and then abs (Pr - W * X) <= Un * abs (W * X) + E,
          Post => Subst'Result;

end Lanyon_Rb_Flux_Pkg;
