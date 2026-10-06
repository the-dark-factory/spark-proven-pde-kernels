with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Quot_Pkg; use Lanyon_Burg1d_Quot_Pkg;
package Lanyon_Burg1d_Waves_Pkg with SPARK_Mode is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. The Lax-Friedrichs wave family, in the reference's operation order, for EVERY pair of states. At the
   --  sonic point (both states exactly zero, so the divisor is zero) the reference divides 0 by 0 and returns NaN
   --  (reproduced on jasisz's sin(2 pi x) case); here the waves are DEFINED as zero there, the exact real limit
   --  (|q| <= |f_R - f_L|), seat ruling R-BU-1.

   Jump_Bound : constant := 2.0E6;
   Wave_Bound : constant := 2.0E6 + 64.0;

   subtype Jump_Value is Long_Float range -Jump_Bound .. Jump_Bound;
   subtype Quot_Value is Long_Float range -Quot_Bound .. Quot_Bound;
   subtype Wave_Value is Long_Float range -Wave_Bound .. Wave_Bound;

   type Wave_Family is record
      Wave1 : Wave_Value;
      Wave2 : Wave_Value;
   end record;

   --  the jump: f_R - f_L
   function Jump (UL : State; UR : State) return Jump_Value is
     (UR.F - UL.F);

   --  the quotient (f_R^2/2 - f_L^2/2) / max (|f_L|, |f_R|), bounded for every nonzero divisor (no floor)
   function Quotient (UL : State; UR : State) return Quot_Value is
     (Num (UL.F, UR.F) / Amax (UL, UR))
     with Pre => Amax (UL, UR) > 0.0 and then Quot_Safe (UL.F, UR.F);

   --  the Lax-Friedrichs waves: (jump - q) / 2 and (jump + q) / 2; both zero at the sonic point
   function X_Wave_Family (UL : State; UR : State) return Wave_Family is
     (if Amax (UL, UR) = 0.0 then (Wave1 => 0.0, Wave2 => 0.0)
      else (Wave1 => 0.5 * (Jump (UL, UR) - Quotient (UL, UR)),
            Wave2 => 0.5 * (Jump (UL, UR) + Quotient (UL, UR))));

   --  THEOREM (Lean xWaveConsistency): equal states give two zero waves, exactly
   function Wave_Consistency (U : State) return Boolean is
     (X_Wave_Family (U, U).Wave1 = 0.0 and then X_Wave_Family (U, U).Wave2 = 0.0)
     with Post => Wave_Consistency'Result;

   --  THEOREM (the sonic point): two zero states give zero waves and a zero jump: the jump condition holds exactly
   function Sonic_Defined (UL : State; UR : State) return Boolean is
     (if UL.F = 0.0 and then UR.F = 0.0
      then X_Wave_Family (UL, UR).Wave1 = 0.0 and then X_Wave_Family (UL, UR).Wave2 = 0.0
           and then Jump (UL, UR) = 0.0)
     with Post => Sonic_Defined'Result;

end Lanyon_Burg1d_Waves_Pkg;
