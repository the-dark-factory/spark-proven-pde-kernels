with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
package Lanyon_Adv1d_Minmod_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The minmod reconstruction, and the Lean reconstruction theorems in the forms that hold in floats.

   Slope_Bound : constant := 2.0E6;
   Recon_Bound : constant := 2.0E6;
   subtype Slope is Long_Float range -Slope_Bound .. Slope_Bound;
   subtype Recon_Value is Long_Float range -Recon_Bound .. Recon_Bound;

   --  f - f_L
   function Back_Slope (UL : State; U : State) return Slope is
     (U.F - UL.F);

   --  f_R - f
   function Fore_Slope (U : State; UR : State) return Slope is
     (UR.F - U.F);

   --  max(0, min(d1, d2)) + min(0, max(d1, d2))
   function Minmod (D1 : Slope; D2 : Slope) return Slope is
     (Long_Float'Max (0.0, Long_Float'Min (D1, D2)) + Long_Float'Min (0.0, Long_Float'Max (D1, D2)));

   --  the left reconstruction: f - minmod / 2
   function X_Left_Reconstruction (UL : State; U : State; UR : State) return Recon_Value is
     (U.F - 0.5 * Minmod (Back_Slope (UL, U), Fore_Slope (U, UR)));

   --  the right reconstruction: f + minmod / 2
   function X_Right_Reconstruction (UL : State; U : State; UR : State) return Recon_Value is
     (U.F + 0.5 * Minmod (Back_Slope (UL, U), Fore_Slope (U, UR)));

   --  THEOREM [Lean xLeftReconstructionConsistent], exact
   function Left_Consistent (U : State) return Boolean is
     (X_Left_Reconstruction (U, U, U) = U.F)
     with Post => Left_Consistent'Result;

   --  THEOREM [Lean xRightReconstructionConsistent], exact
   function Right_Consistent (U : State) return Boolean is
     (X_Right_Reconstruction (U, U, U) = U.F)
     with Post => Right_Consistent'Result;

   --  minmod of two equal slopes is that slope
   function Minmod_Equal (D1 : Slope; D2 : Slope) return Boolean is
     (if D1 = D2 then Minmod (D1, D2) = D1)
     with Post => Minmod_Equal'Result;

   --  THEOREM [Lean xLeftReconstructionLinearityPreservation], RESTATED: when the two computed
   --  slopes are equal the reconstruction is f - 0.5 * slope, exactly (the Lean f - D and f + D are rounded)
   function Left_Linear (UL : State; U : State; UR : State) return Boolean is
     (if Back_Slope (UL, U) = Fore_Slope (U, UR)
      then Minmod_Equal (Back_Slope (UL, U), Fore_Slope (U, UR))
           and then X_Left_Reconstruction (UL, U, UR) = U.F - 0.5 * Back_Slope (UL, U))
     with Post => Left_Linear'Result;

   --  THEOREM [Lean xRightReconstructionLinearityPreservation], RESTATED as Left_Linear
   function Right_Linear (UL : State; U : State; UR : State) return Boolean is
     (if Back_Slope (UL, U) = Fore_Slope (U, UR)
      then Minmod_Equal (Back_Slope (UL, U), Fore_Slope (U, UR))
           and then X_Right_Reconstruction (UL, U, UR) = U.F + 0.5 * Back_Slope (UL, U))
     with Post => Right_Linear'Result;

end Lanyon_Adv1d_Minmod_Pkg;
