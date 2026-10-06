with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
with Lanyon_Rb_Core_Pkg; use Lanyon_Rb_Core_Pkg;
with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
package Lanyon_Burg1d_Quot_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. The level-1 bound on the wave split's quotient q = (f_R^2/2 - f_L^2/2) / max (|f_L|, |f_R|), for EVERY
   --  pair of states with a nonzero divisor, subnormal included: no floor on the state. Relational, so it is carried
   --  by exact real (big real) facts and brought back to a float range at the end.

   subtype Any_Pos is Long_Float range 2.0 ** (-1074) .. Field_Bound;
   subtype Sq_Value is Long_Float range 0.0 .. 1.0E12 + 1.0E6;

   --  dividing by any positive double, given |x| <= 2 y^2 + 8 eta (which keeps x / y finite)
   function Div_Any (X : Num_Value; Y : Any_Pos) return Boolean is
     (abs (R (X / Y) - R (X) / R (Y)) <= Un * abs (R (X) / R (Y)) + Eta)
     with Pre  => abs R (X) <= 2.0 * (R (Y) * R (Y)) + 8.0 * Eta,
          Post => Div_Any'Result;

   function Eta_Le (Y : Any_Pos) return Boolean is (Eta > 0.0 and then Eta <= R (Y) and then R (Y) <= 1000000.0)
     with Post => Eta_Le'Result;

   function R_Abs (X : Field) return Boolean is (R (abs X) = abs R (X))
     with Post => R_Abs'Result;

   function R_Mono (X, Y : Field) return Boolean is (if X >= Y then R (X) >= R (Y))
     with Post => R_Mono'Result;

   --  the real value of the divisor bounds both states
   function Max_Real (FL, FR : Field) return Boolean is
     (R_Abs (FL) and then R_Abs (FR)
      and then R_Mono (Long_Float'Max (abs FL, abs FR), abs FL) and then R_Mono (Long_Float'Max (abs FL, abs FR), abs FR))
     with Post => Max_Real'Result
                  and then R (Long_Float'Max (abs FL, abs FR)) >= abs R (FL)
                  and then R (Long_Float'Max (abs FL, abs FR)) >= abs R (FR);

   function Sq_Le (F, A : Valid_Big_Real) return Boolean is (F * F <= A * A)
     with Pre => abs F <= A, Post => Sq_Le'Result;

   --  pure real: the numerator is at most 2 A^2 + 8 eta
   function Num_Real (N, H1, H2, P1, P2, FR, FL, A, E : Valid_Big_Real) return Boolean is
     (abs N <= 2.0 * (A * A) + 8.0 * E)
     with Pre  => E >= 0.0 and then E <= 0.000000000000000000001
                  and then FR * FR <= A * A and then FL * FL <= A * A
                  and then abs (P1 - FR * FR) <= Un * abs (FR * FR) + E
                  and then abs (P2 - FL * FL) <= Un * abs (FL * FL) + E
                  and then abs (H1 - 0.5 * P1) <= Un * abs (0.5 * P1) + E
                  and then abs (H2 - 0.5 * P2) <= Un * abs (0.5 * P2) + E
                  and then abs (N - (H1 - H2)) <= Un * abs (H1 - H2) + E,
          Post => Num_Real'Result;

   --  pure real: E / A <= 1 when 0 < E <= A
   function E_Over_A (E, A : Valid_Big_Real) return Boolean is (E / A <= 1.0)
     with Pre => E > 0.0 and then E <= A, Post => E_Over_A'Result;

   --  pure real: |N / A| <= 2 A + 8 E / A
   function N_Over_A (N, A, E : Valid_Big_Real) return Boolean is
     (abs (N / A) <= 2.0 * A + 8.0 * (E / A))
     with Pre  => A > 0.0 and then E >= 0.0 and then abs N <= 2.0 * (A * A) + 8.0 * E,
          Post => N_Over_A'Result;

   --  pure real: the quotient is at most 2 A (1 + u) + 9 (eta <= A <= 1e6)
   function Quot_Real (Q, N, A, E : Valid_Big_Real) return Boolean is
     (E_Over_A (E, A) and then N_Over_A (N, A, E) and then abs (N / A) <= 2000008.0 and then abs Q <= 2000032.0)
     with Pre  => E > 0.0 and then E <= A and then A <= 1000000.0 and then E <= 0.000000000000000000001
                  and then abs N <= 2.0 * (A * A) + 8.0 * E
                  and then abs (Q - N / A) <= Un * abs (N / A) + E,
          Post => Quot_Real'Result;

   function R_Back (X : Long_Float) return Boolean is
     (if abs R (X) <= 2000032.0 then X >= -Quot_Bound and then X <= Quot_Bound)
     with Post => R_Back'Result;

   --  the numerator's real-value facts
   function Num_Fact (FL, FR : Field) return Boolean is
     (U_Is_Un and then Half_Exact and then Eta_Small and then Max_Real (FL, FR)
      and then Sq_Le (R (FR), R (Long_Float'Max (abs FL, abs FR)))
      and then Sq_Le (R (FL), R (Long_Float'Max (abs FL, abs FR)))
      and then Mul_Err (FR, FR) and then Mul_Err (FL, FL)
      and then Mul_Err (0.5, FR * FR) and then Mul_Err (0.5, FL * FL)
      and then Sub_Err (0.5 * (FR * FR), 0.5 * (FL * FL))
      and then Num_Real (R (Num (FL, FR)), R (0.5 * (FR * FR)), R (0.5 * (FL * FL)), R (FR * FR), R (FL * FL),
                         R (FR), R (FL), R (Long_Float'Max (abs FL, abs FR)), Eta))
     with Post => Num_Fact'Result
                  and then abs R (Num (FL, FR)) <= 2.0 * (R (Long_Float'Max (abs FL, abs FR)) * R (Long_Float'Max (abs FL, abs FR))) + 8.0 * Eta;

   --  THEOREM (level 1, no floor): for every pair of states with a nonzero divisor, q lies in the bounded range
   function Quot_Safe (FL, FR : Field) return Boolean is
     (Num_Fact (FL, FR) and then Eta_Le (Long_Float'Max (abs FL, abs FR))
      and then Div_Any (Num (FL, FR), Long_Float'Max (abs FL, abs FR))
      and then Quot_Real (R (Num (FL, FR) / Long_Float'Max (abs FL, abs FR)), R (Num (FL, FR)),
                          R (Long_Float'Max (abs FL, abs FR)), Eta)
      and then R_Back (Num (FL, FR) / Long_Float'Max (abs FL, abs FR)))
     with Pre  => Long_Float'Max (abs FL, abs FR) > 0.0,
          Post => Quot_Safe'Result
                  and then Num (FL, FR) / Long_Float'Max (abs FL, abs FR) >= -Quot_Bound
                  and then Num (FL, FR) / Long_Float'Max (abs FL, abs FR) <= Quot_Bound;

end Lanyon_Burg1d_Quot_Pkg;
