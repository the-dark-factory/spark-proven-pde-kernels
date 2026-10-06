package Lanyon_Adv1d_Types_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The state, parameters, flux, wave speed and speed family of 1-D linear advection, f_t + a f_x = 0.

   --  the declared bounds: every state value and every speed lies within them
   Field_Bound : constant := 1.0E6;
   Speed_Bound : constant := 1.0E3;
   Flux_Bound  : constant := 1.0E9;

   --  the smallest valid speed magnitude (seat ruling R-LA-3; Lanyon's validity check accepts any a)
   A_Min : constant := 1.0E-6;

   subtype Field is Long_Float range -Field_Bound .. Field_Bound;
   subtype Speed is Long_Float range -Speed_Bound .. Speed_Bound;
   subtype Abs_Speed is Long_Float range A_Min .. Speed_Bound;
   subtype Flux_Value is Long_Float range -Flux_Bound .. Flux_Bound;

   type State is record
      F : Field;
   end record;

   type Parameters is record
      A : Speed;
   end record;

   type Flux is record
      Flux_F : Flux_Value;
   end record;

   type Wave_Speed is record
      Mu1 : Speed;
   end record;

   type Eigen is record
      Lambda1 : Speed;
   end record;

   type Speed_Family is record
      Speed1 : Speed;
      Speed2 : Speed;
   end record;

   --  Lanyon's validity check accepts every a; here a speed below A_Min is invalid (R-LA-3)
   function Parameters_Valid (P : Parameters) return Boolean is
     (abs P.A >= A_Min);

   --  the flux: a * f
   function X_Flux (P : Parameters; U : State) return Flux is
     ((Flux_F => P.A * U.F))
     with Post => (if U.F = 0.0 or else P.A = 0.0 then X_Flux'Result.Flux_F = 0.0)
                  and then (if P.A >= 0.0 and then U.F >= 0.0 then X_Flux'Result.Flux_F >= 0.0)
                  and then (if P.A >= 0.0 and then U.F <= 0.0 then X_Flux'Result.Flux_F <= 0.0)
                  and then (if P.A <= 0.0 and then U.F >= 0.0 then X_Flux'Result.Flux_F <= 0.0)
                  and then (if P.A <= 0.0 and then U.F <= 0.0 then X_Flux'Result.Flux_F >= 0.0);

   --  the wave speed: a
   function X_Wavespeed (P : Parameters) return Wave_Speed is
     ((Mu1 => P.A));

   --  the flux Jacobian's eigenvalue: a
   function X_Eigen (P : Parameters) return Eigen is
     ((Lambda1 => P.A));

   --  the diffusive flux: zero (pure advection has no diffusive part)
   function X_Diffusive_Flux return Long_Float is
     (0.0);

   --  the Lax-Friedrichs speeds: -|a| and |a|
   function X_Speed_Family (P : Parameters) return Speed_Family is
     ((Speed1 => -(abs P.A), Speed2 => abs P.A));

   --  |a| of a valid parameter, as a positive bounded value
   function Abs_A (P : Parameters) return Abs_Speed is
     (abs P.A)
     with Pre => Parameters_Valid (P);

   --  THEOREM [Lean xWaveStability]: |mu1| >= |lambda1|
   function Wave_Stability (P : Parameters) return Boolean is
     (abs X_Wavespeed (P).Mu1 >= abs X_Eigen (P).Lambda1)
     with Post => Wave_Stability'Result;

   --  THEOREM (ties the eigenvalue to the flux; upstream #2 notes that Lean xHyperbolicity does not):
   --  the flux is the eigenvalue times the state, exactly
   function Flux_Is_Linear (P : Parameters; U : State) return Boolean is
     (X_Flux (P, U).Flux_F = X_Eigen (P).Lambda1 * U.F)
     with Post => Flux_Is_Linear'Result;

   --  THEOREM [Lean xDiffusiveFluxConsistency]: the diffusive flux of a zero gradient is zero
   function Diffusive_Consistency return Boolean is
     (X_Diffusive_Flux = 0.0)
     with Post => Diffusive_Consistency'Result;

   --  THEOREM (the speeds): the two speeds have opposite signs and equal magnitude
   function Speeds_Opposite (P : Parameters) return Boolean is
     (X_Speed_Family (P).Speed1 <= 0.0
      and then X_Speed_Family (P).Speed2 >= 0.0
      and then X_Speed_Family (P).Speed1 = -X_Speed_Family (P).Speed2)
     with Post => Speeds_Opposite'Result;

end Lanyon_Adv1d_Types_Pkg;
