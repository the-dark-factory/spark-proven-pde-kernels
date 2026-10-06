with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Waves_Pkg; use Lanyon_Burg1d_Waves_Pkg;
package Lanyon_Burg1d_Fluct_Pkg with SPARK_Mode is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. The left and right fluctuations.

   Term_Bound  : constant := 2.1E12;
   Fluct_Bound : constant := 4.2E12;
   subtype Term_Value is Long_Float range -Term_Bound .. Term_Bound;
   subtype Fluct_Value is Long_Float range -Fluct_Bound .. Fluct_Bound;

   type Fluctuation is record
      Fluctuation_F : Fluct_Value;
   end record;

   --  min (speed, 0) and max (speed, 0)
   function Neg_Part (S : Field) return Field is (Long_Float'Min (S, 0.0));
   function Pos_Part (S : Field) return Field is (Long_Float'Max (S, 0.0));

   --  one wave times one speed part
   function Term (W : Wave_Value; S : Field) return Term_Value is (W * S);

   --  the left fluctuation: the sum of each wave times min (its speed, 0)
   function X_Left_Fluctuation (UL : State; UR : State) return Fluctuation is
     ((Fluctuation_F =>
         Term (X_Wave_Family (UL, UR).Wave1, Neg_Part (X_Speed_Family (UL, UR).Speed1))
         + Term (X_Wave_Family (UL, UR).Wave2, Neg_Part (X_Speed_Family (UL, UR).Speed2))));

   --  the right fluctuation: the sum of each wave times max (its speed, 0)
   function X_Right_Fluctuation (UL : State; UR : State) return Fluctuation is
     ((Fluctuation_F =>
         Term (X_Wave_Family (UL, UR).Wave1, Pos_Part (X_Speed_Family (UL, UR).Speed1))
         + Term (X_Wave_Family (UL, UR).Wave2, Pos_Part (X_Speed_Family (UL, UR).Speed2))));

   --  THEOREM (Lean xLeftFluctuationsConsistent): equal states give a zero left fluctuation, exactly
   function Left_Consistent (U : State) return Boolean is
     (Wave_Consistency (U) and then X_Left_Fluctuation (U, U).Fluctuation_F = 0.0)
     with Post => Left_Consistent'Result;

   --  THEOREM (Lean xRightFluctuationsConsistent): equal states give a zero right fluctuation, exactly
   function Right_Consistent (U : State) return Boolean is
     (Wave_Consistency (U) and then X_Right_Fluctuation (U, U).Fluctuation_F = 0.0)
     with Post => Right_Consistent'Result;

   --  THEOREM (the sonic point): two zero states give zero fluctuations, exactly
   function Sonic_Quiet (UL : State; UR : State) return Boolean is
     (if UL.F = 0.0 and then UR.F = 0.0
      then Sonic_Defined (UL, UR)
           and then X_Left_Fluctuation (UL, UR).Fluctuation_F = 0.0
           and then X_Right_Fluctuation (UL, UR).Fluctuation_F = 0.0)
     with Post => Sonic_Quiet'Result;

end Lanyon_Burg1d_Fluct_Pkg;
