with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
package Lanyon_Adv1d_Waves_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The Lax-Friedrichs wave family, in the reference's operation order, with every intermediate bounded.

   Jump_Bound : constant := 2.0E6;
   Num_Bound  : constant := 2.0E9;
   Quot_Bound : constant := 2.0E15 + 1.0;
   Wave_Bound : constant := 1.0E15 + 2.0E6;

   subtype Jump_Value is Long_Float range -Jump_Bound .. Jump_Bound;
   subtype Num_Value is Long_Float range -Num_Bound .. Num_Bound;
   subtype Quot_Value is Long_Float range -Quot_Bound .. Quot_Bound;
   subtype Wave_Value is Long_Float range -Wave_Bound .. Wave_Bound;

   type Wave_Family is record
      Wave1 : Wave_Value;
      Wave2 : Wave_Value;
   end record;

   --  the jump: f_R - f_L
   function Jump (UL : State; UR : State) return Jump_Value is
     (UR.F - UL.F);

   --  (a * f_R) - (a * f_L)
   function Flux_Difference (UL : State; UR : State; P : Parameters) return Num_Value is
     (P.A * UR.F - P.A * UL.F);

   --  ((a * f_R) - (a * f_L)) / |a|
   function Quotient (UL : State; UR : State; P : Parameters) return Quot_Value is
     (Flux_Difference (UL, UR, P) / Abs_A (P))
     with Pre => Parameters_Valid (P);

   --  the Lax-Friedrichs waves: (jump - q) / 2 and (jump + q) / 2
   function X_Wave_Family (UL : State; UR : State; P : Parameters) return Wave_Family is
     ((Wave1 => 0.5 * (Jump (UL, UR) - Quotient (UL, UR, P)),
       Wave2 => 0.5 * (Jump (UL, UR) + Quotient (UL, UR, P))))
     with Pre => Parameters_Valid (P);

   --  THEOREM [Lean xWaveConsistency]: equal states give two zero waves, exactly
   function Wave_Consistency (U : State; P : Parameters) return Boolean is
     (X_Wave_Family (U, U, P).Wave1 = 0.0 and then X_Wave_Family (U, U, P).Wave2 = 0.0)
     with Pre  => Parameters_Valid (P),
          Post => Wave_Consistency'Result;

end Lanyon_Adv1d_Waves_Pkg;
