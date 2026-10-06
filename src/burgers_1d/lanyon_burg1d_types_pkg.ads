package Lanyon_Burg1d_Types_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081: the state, flux f*f/2, wave speed f, and the Lax-Friedrichs speed family of f_t + (f^2/2)_x = 0.

   Field_Bound : constant := 1.0E6;
   subtype Field is Long_Float range -Field_Bound .. Field_Bound;
   subtype Mag is Long_Float range 0.0 .. Field_Bound;
   subtype Flux_Value is Long_Float range 0.0 .. 5.0E11 + 1.0E6;

   --  the wave split's numerator and quotient ranges (the quotient bound is proved in Lanyon_Burg1d_Quot_Pkg)
   Quot_Bound : constant := 2.0E6 + 32.0;
   subtype Num_Value is Long_Float range -5.0E11 - 1.0E6 .. 5.0E11 + 1.0E6;

   type State is record
      F : Field;
   end record;

   type Flux is record
      Flux_F : Flux_Value;
   end record;

   type Speed_Family is record
      Speed1 : Field;
      Speed2 : Field;
   end record;

   --  the flux, in the reference's order: (0.5 * f) * f
   function X_Flux (U : State) return Flux is
     ((Flux_F => 0.5 * U.F * U.F))
     with Post => (if U.F = 0.0 then X_Flux'Result.Flux_F = 0.0);

   --  the wave split's numerator, in the reference's order: 0.5 * (f_R * f_R) - 0.5 * (f_L * f_L)
   function Num (FL : Field; FR : Field) return Num_Value is
     (0.5 * (FR * FR) - 0.5 * (FL * FL));

   --  the wave speed and the flux Jacobian's eigenvalue: both f
   function X_Wavespeed (U : State) return Field is (U.F);
   function X_Eigen (U : State) return Field is (U.F);

   --  the Lax-Friedrichs bound: max (|f_L|, |f_R|), the divisor of the wave split
   function Amax (UL : State; UR : State) return Mag is
     (Long_Float'Max (abs UL.F, abs UR.F));

   function X_Speed_Family (UL : State; UR : State) return Speed_Family is
     ((Speed1 => -Amax (UL, UR), Speed2 => Amax (UL, UR)));

   --  THEOREM (Lean xWaveStability): |mu1| >= |lambda1|
   function Wave_Stability (U : State) return Boolean is
     (abs X_Wavespeed (U) >= abs X_Eigen (U))
     with Post => Wave_Stability'Result;

   --  THEOREM (ties the eigenvalue to the flux; upstream #2 in the sister repository): flux = (0.5 * eigenvalue) * state
   function Flux_Is_Jacobian (U : State) return Boolean is
     (X_Flux (U).Flux_F = 0.5 * X_Eigen (U) * U.F)
     with Post => Flux_Is_Jacobian'Result;

   --  THEOREM (the speeds): opposite signs, equal magnitude
   function Speeds_Opposite (UL : State; UR : State) return Boolean is
     (X_Speed_Family (UL, UR).Speed1 <= 0.0 and then X_Speed_Family (UL, UR).Speed2 >= 0.0
      and then X_Speed_Family (UL, UR).Speed1 = -X_Speed_Family (UL, UR).Speed2)
     with Post => Speeds_Opposite'Result;

   --  THE SONIC POINT: the divisor is zero exactly when both states are exactly zero
   function Sonic_Exact (UL : State; UR : State) return Boolean is
     ((Amax (UL, UR) = 0.0) = (UL.F = 0.0 and then UR.F = 0.0))
     with Post => Sonic_Exact'Result;

end Lanyon_Burg1d_Types_Pkg;
