with Lanyon_Adv1d_Types_Pkg; use Lanyon_Adv1d_Types_Pkg;
with Lanyon_Adv1d_Fluct_Pkg; use Lanyon_Adv1d_Fluct_Pkg;
package Lanyon_Adv1d_Update_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  The first-order upwind cell update in flux form on a periodic grid. Lanyon's code has no update;
   --  the interface flux is f(Q_L) + A-dQ (wave-propagation form), built from Lanyon's fluctuation.

   N_Cells     : constant := 256;
   Iface_Bound : constant := 3.0E18;
   Wide_Bound  : constant := 7.0E18;

   subtype Ratio is Long_Float range 0.0 .. 1.0;
   subtype Iface_Value is Long_Float range -Iface_Bound .. Iface_Bound;
   subtype Wide_Value is Long_Float range -Wide_Bound .. Wide_Bound;
   subtype Cell_Index is Positive range 1 .. N_Cells;
   type Cells is array (Cell_Index) of Field;
   type Wide_Cells is array (Cell_Index) of Wide_Value;

   --  periodic neighbours
   function Left_Of (I : Cell_Index) return Cell_Index is
     (if I = 1 then N_Cells else I - 1);

   function Right_Of (I : Cell_Index) return Cell_Index is
     (if I = N_Cells then 1 else I + 1);

   --  cell I as a state
   function Cell (Q : Cells; I : Cell_Index) return State is
     ((F => Q (I)));

   --  the interface flux between QL and QR: f(QL) + A-dQ, the one value both cells use
   function Interface_Flux (QL : State; QR : State; P : Parameters) return Iface_Value is
     (X_Flux (P, QL).Flux_F + X_Left_Fluctuation (QL, QR, P).Fluctuation_F)
     with Pre => Parameters_Valid (P);

   --  the flux leaving cell I through its right face
   function Out_Flux (Q : Cells; P : Parameters; I : Cell_Index) return Iface_Value is
     (Interface_Flux (Cell (Q, I), Cell (Q, Right_Of (I)), P))
     with Pre => Parameters_Valid (P);

   --  the flux entering cell I through its left face
   function In_Flux (Q : Cells; P : Parameters; I : Cell_Index) return Iface_Value is
     (Interface_Flux (Cell (Q, Left_Of (I)), Cell (Q, I), P))
     with Pre => Parameters_Valid (P);

   --  Q_i - nu * (F_{i+1/2} - F_{i-1/2}), nu = dt / dx
   function Update_Cell (Q : Cells; P : Parameters; Nu : Ratio; I : Cell_Index) return Wide_Value is
     (Q (I) - Nu * (Out_Flux (Q, P, I) - In_Flux (Q, P, I)))
     with Pre => Parameters_Valid (P);

   --  one step of every cell, periodic: cell 1's left neighbour is cell N_Cells
   function Step (Q : Cells; P : Parameters; Nu : Ratio) return Wide_Cells is
     ((for I in Cell_Index => Update_Cell (Q, P, Nu, I)))
     with Pre => Parameters_Valid (P);

   --  THEOREM: the value cell I sends through its right face is the value cell I+1 receives
   function Interface_Shared (Q : Cells; P : Parameters; I : Cell_Index) return Boolean is
     (Out_Flux (Q, P, I) = In_Flux (Q, P, Right_Of (I)))
     with Pre  => Parameters_Valid (P),
          Post => Interface_Shared'Result;

   --  the grid whose every cell holds U
   function Uniform (U : State) return Cells is
     ((for I in Cell_Index => U.F));

   --  THEOREM: one cell of a constant state is unchanged by the update, exactly
   function Cell_Fixed (U : State; P : Parameters; Nu : Ratio; I : Cell_Index) return Boolean is
     (Left_Consistent (U, P)
      and then Interface_Flux (U, U, P) = X_Flux (P, U).Flux_F
      and then Update_Cell (Uniform (U), P, Nu, I) = U.F)
     with Pre  => Parameters_Valid (P),
          Post => Cell_Fixed'Result;

   --  THEOREM: a constant state is a fixed point of the step, exactly
   function Constant_Fixed (U : State; P : Parameters; Nu : Ratio) return Boolean is
     (for all I in Cell_Index =>
        Cell_Fixed (U, P, Nu, I) and then Step (Uniform (U), P, Nu) (I) = U.F)
     with Pre  => Parameters_Valid (P),
          Post => Constant_Fixed'Result;

end Lanyon_Adv1d_Update_Pkg;
