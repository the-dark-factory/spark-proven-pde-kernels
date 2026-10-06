with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
with Lanyon_Adv1d_Waves_Pkg; use Lanyon_Adv1d_Waves_Pkg;
package Lanyon_Adv1d_Fluct_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The left and right fluctuations.

   Term_Bound  : constant := 1.0E18 + 1.0E10;
   Fluct_Bound : constant := 2.0E18 + 1.0E10;

   subtype Term_Value is Long_Float range -Term_Bound .. Term_Bound;
   subtype Fluct_Value is Long_Float range -Fluct_Bound .. Fluct_Bound;

   type Fluctuation is record
      Fluctuation_F : Fluct_Value;
   end record;

   --  min (speed, 0)
   function Neg_Part (S : Speed) return Speed is
     (Long_Float'Min (S, 0.0));

   --  max (speed, 0)
   function Pos_Part (S : Speed) return Speed is
     (Long_Float'Max (S, 0.0));

   --  one wave times one speed part
   function Term (W : Wave_Value; S : Speed) return Term_Value is
     (W * S);

   --  the left fluctuation: the sum of each wave times min (its speed, 0)
   function X_Left_Fluctuation (UL : State; UR : State; P : Parameters) return Fluctuation is
     ((Fluctuation_F =>
         Term (X_Wave_Family (UL, UR, P).Wave1, Neg_Part (X_Speed_Family (P).Speed1))
         + Term (X_Wave_Family (UL, UR, P).Wave2, Neg_Part (X_Speed_Family (P).Speed2))))
     with Pre => Parameters_Valid (P);

   --  the right fluctuation: the sum of each wave times max (its speed, 0)
   function X_Right_Fluctuation (UL : State; UR : State; P : Parameters) return Fluctuation is
     ((Fluctuation_F =>
         Term (X_Wave_Family (UL, UR, P).Wave1, Pos_Part (X_Speed_Family (P).Speed1))
         + Term (X_Wave_Family (UL, UR, P).Wave2, Pos_Part (X_Speed_Family (P).Speed2))))
     with Pre => Parameters_Valid (P);

   --  THEOREM [Lean xLeftFluctuationsConsistent]: equal states give a zero left fluctuation, exactly
   function Left_Consistent (U : State; P : Parameters) return Boolean is
     (Wave_Consistency (U, P) and then X_Left_Fluctuation (U, U, P).Fluctuation_F = 0.0)
     with Pre  => Parameters_Valid (P),
          Post => Left_Consistent'Result;

   --  THEOREM [Lean xRightFluctuationsConsistent]: equal states give a zero right fluctuation, exactly
   function Right_Consistent (U : State; P : Parameters) return Boolean is
     (Wave_Consistency (U, P) and then X_Right_Fluctuation (U, U, P).Fluctuation_F = 0.0)
     with Pre  => Parameters_Valid (P),
          Post => Right_Consistent'Result;

end Lanyon_Adv1d_Fluct_Pkg;
