with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
package Lanyon_Adv1d_Rounding_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  Per-operation rounding lemmas for IEEE binary64 (round to nearest) as exact real facts, and the pure
   --  real lemmas the bounded theorems are assembled from. u = 2**-53, eta = 2**-1074.

   package LF is new Float_Conversions (Long_Float);

   subtype Big is Long_Float range -1.0E19 .. 1.0E19;
   subtype Pos is Long_Float range 1.0E-6 .. 1.0E3;

   function R (X : Long_Float) return Valid_Big_Real is (LF.To_Big_Real (X));
   function U return Valid_Big_Real is (R (2.0 ** (-53)));
   function Eta return Valid_Big_Real is (R (2.0 ** (-1074)));

   function Sub_Err (X : Big; Y : Big) return Boolean is
     (abs (R (X - Y) - (R (X) - R (Y))) <= U * abs (R (X) - R (Y)) + Eta)
     with Post => Sub_Err'Result;

   function Add_Err (X : Big; Y : Big) return Boolean is
     (abs (R (X + Y) - (R (X) + R (Y))) <= U * abs (R (X) + R (Y)) + Eta)
     with Post => Add_Err'Result;

   function Mul_Err (X : Big; Y : Big) return Boolean is
     (abs (R (X * Y) - R (X) * R (Y)) <= U * abs (R (X) * R (Y)) + Eta)
     with Post => Mul_Err'Result;

   function Div_Err (X : Big; Y : Pos) return Boolean is
     (abs (R (X / Y) - R (X) / R (Y)) <= U * abs (R (X) / R (Y)) + Eta)
     with Post => Div_Err'Result;


   --  2**-53, exactly
   Un : constant Big_Real := 0.00000000000000011102230246251565404236316680908203125;

   --  a product is within (1+u) of its magnitude
   function Prod_Mag (P, Av, F, E : Valid_Big_Real) return Boolean is
     (abs P <= (1.0 + Un) * (abs Av * abs F) + E)
     with Pre  => E >= 0.0 and then abs (P - Av * F) <= Un * abs (Av * F) + E,
          Post => Prod_Mag'Result;

   --  the flux difference: |N| <= (1+u)^2 A (|fR| + |fL|) + 4E
   function Num_Mag (N, P1, P2, A, FR, FL, E : Valid_Big_Real) return Boolean is
     (abs N <= (1.0 + Un) * (1.0 + Un) * (A * abs FR + A * abs FL) + 4.0 * E)
     with Pre  => A >= 0.0 and then E >= 0.0
                  and then abs P1 <= (1.0 + Un) * (A * abs FR) + E
                  and then abs P2 <= (1.0 + Un) * (A * abs FL) + E
                  and then abs (N - (P1 - P2)) <= Un * abs (P1 - P2) + E,
          Post => Num_Mag'Result;

   --  dividing by A >= 1e-6
   function Div_Mag (Q, N, A, M, E : Valid_Big_Real) return Boolean is
     (abs Q <= (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M + 5000000.0 * E)
     with Pre  => A >= 0.0000009 and then E >= 0.0 and then M >= 0.0
                  and then abs N <= (1.0 + Un) * (1.0 + Un) * (A * M) + 4.0 * E
                  and then abs (Q - N / A) <= Un * abs (N / A) + E,
          Post => Div_Mag'Result;

   --  |x - y| + |x + y| = 2 max (|x|, |y|) <= 2 K
   function Two_Max (X, Y, K : Valid_Big_Real) return Boolean is
     (abs (X - Y) + abs (X + Y) <= 2.0 * K)
     with Pre  => abs X <= K and then abs Y <= K,
          Post => Two_Max'Result;

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


   --  float equality is not substitution (0.0 = -0.0), but equal floats have equal real values
   function R_Cong (X : Long_Float; Y : Long_Float) return Boolean is
     (if X = Y then R (X) = R (Y))
     with Post => R_Cong'Result;

   function Half_Exact return Boolean is (R (0.5) = 0.5) with Post => Half_Exact'Result;

   function U_Is_Un return Boolean is (U = Un) with Post => U_Is_Un'Result;

   function Eta_Small return Boolean is (Eta >= 0.0 and then Eta <= 0.000000000000000000001)
     with Post => Eta_Small'Result;


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

   --  (d1) multiplying the quotient error by A > 0
   function Aq_Step (Q, N, A, E : Valid_Big_Real) return Boolean is
     (abs (A * Q - N) <= Un * abs N + A * E)
     with Pre  => A >= 0.0000009 and then E >= 0.0
                  and then abs (Q - N / A) <= Un * abs (N / A) + E,
          Post => Aq_Step'Result;

   --  |Av * F| = A |F| when A = |Av|
   function Abs_Prod (Av, A, F : Valid_Big_Real) return Boolean is
     (abs (Av * F) = A * abs F)
     with Pre  => A = abs Av,
          Post => Abs_Prod'Result;

   --  (d2) the rounded flux difference is the exact flux jump, to rounding
   function N_Near (N, P1, P2, Av, A, FR, FL, E : Valid_Big_Real) return Boolean is
     (abs (N - Av * (FR - FL)) <= Un * abs (P1 - P2) + Un * (A * abs FR + A * abs FL) + 3.0 * E)
     with Pre  => E >= 0.0 and then A = abs Av
                  and then abs (Av * FR) = A * abs FR and then abs (Av * FL) = A * abs FL
                  and then abs (P1 - Av * FR) <= Un * abs (Av * FR) + E
                  and then abs (P2 - Av * FL) <= Un * abs (Av * FL) + E
                  and then abs (N - (P1 - P2)) <= Un * abs (P1 - P2) + E,
          Post => N_Near'Result;

   --  (d) A q is the exact flux jump, to rounding
   function Aq_Near (Q, N, P1, P2, Av, A, FR, FL, E : Valid_Big_Real) return Boolean is
     (abs (A * Q - Av * (FR - FL)) <= Un * abs N + Un * abs (P1 - P2) + Un * (A * abs FR + A * abs FL) + A * E + 3.0 * E)
     with Pre  => E >= 0.0
                  and then abs (A * Q - N) <= Un * abs N + A * E
                  and then abs (N - Av * (FR - FL)) <= Un * abs (P1 - P2) + Un * (A * abs FR + A * abs FL) + 3.0 * E,
          Post => Aq_Near'Result;

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

   --  (g) the flux theorem collapsed to the Lean shape
   function Flux_Final (Err, A, W1, W2, S1, S2, D, Q, N, P1, P2, FR, FL, M, E, K : Valid_Big_Real) return Boolean is
     (Err <= 10.0 * Un * (A * M) + 10000.0 * E)
     with Pre  => E >= 0.0 and then E <= 0.000000000000000000001 and then M >= 0.0
                  and then A >= 0.0 and then A <= 1000.1
                  and then M = abs FR + abs FL
                  and then K = (1.0 + Un) * (1.0 + Un) * (1.0 + Un) * M + 5000000.0 * E
                  and then abs W1 + abs W2 <= 1.001 * K + 4.0 * E
                  and then abs (D - Q) + abs (D + Q) <= 2.0 * K
                  and then abs S1 + abs S2 <= 2.001 * K + 2.0 * E
                  and then abs N <= (1.0 + Un) * (1.0 + Un) * (A * abs FR + A * abs FL) + 4.0 * E
                  and then abs (P1 - P2) <= (1.0 + Un) * (A * abs FR + A * abs FL) + 2.0 * E
                  and then Err <= (2.1 * Un * (A * abs W1) + 4.0 * E) + (2.1 * Un * (A * abs W2) + 4.0 * E)
                                  + A * (Un * (abs (D - Q) + abs (D + Q)) + Un * (abs S1 + abs S2) + 4.0 * E)
                                  + (Un * abs N + Un * abs (P1 - P2) + Un * (A * abs FR + A * abs FL) + A * E + 3.0 * E),
          Post => Flux_Final'Result;

   --  substituting X = Y under the product error bound
   function Subst (T, Pr, W, X, Y, E : Valid_Big_Real) return Boolean is
     (abs (T - W * Y) <= Un * abs (W * Y) + E)
     with Pre  => X = Y and then T = Pr and then abs (Pr - W * X) <= Un * abs (W * X) + E,
          Post => Subst'Result;

   --  |P1 - P2| from the two product magnitudes, A = |Av|
   function Diff_Mag (P1, P2, Av, A, FR, FL, E : Valid_Big_Real) return Boolean is
     (abs (P1 - P2) <= (1.0 + Un) * (A * abs FR + A * abs FL) + 2.0 * E)
     with Pre  => A = abs Av
                  and then abs P1 <= (1.0 + Un) * (abs Av * abs FR) + E
                  and then abs P2 <= (1.0 + Un) * (abs Av * abs FL) + E,
          Post => Diff_Mag'Result;

end Lanyon_Adv1d_Rounding_Pkg;
