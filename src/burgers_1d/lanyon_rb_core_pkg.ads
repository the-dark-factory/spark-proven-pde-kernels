with Ada.Numerics.Big_Numbers.Big_Reals; use Ada.Numerics.Big_Numbers.Big_Reals;
package Lanyon_Rb_Core_Pkg with SPARK_Mode, Ghost is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  Per-operation rounding lemmas for IEEE binary64 as exact real facts. u = 2**-53, eta = 2**-1074.

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

   --  2**-53, exactly
   Un : constant Big_Real := 0.00000000000000011102230246251565404236316680908203125;

   --  float equality is not substitution (0.0 = -0.0), but equal floats have equal real values
   function R_Cong (X : Long_Float; Y : Long_Float) return Boolean is
     (if X = Y then R (X) = R (Y))
     with Post => R_Cong'Result;

   function Half_Exact return Boolean is (R (0.5) = 0.5) with Post => Half_Exact'Result;

   function U_Is_Un return Boolean is (U = Un) with Post => U_Is_Un'Result;

   function Eta_Small return Boolean is (Eta >= 0.0 and then Eta <= 0.000000000000000000001)
     with Post => Eta_Small'Result;

end Lanyon_Rb_Core_Pkg;
