with Lanyon_Burg1d_Types_Pkg; use Lanyon_Burg1d_Types_Pkg;
with Lanyon_Burg1d_Fluct_Pkg; use Lanyon_Burg1d_Fluct_Pkg;
package Lanyon_Burg1d_Update_Pkg with SPARK_Mode is

   --  After Lanyon AI's BurgersEquation (MIT licence, Copyright (c) 2026 Lanyon AI Inc), inviscid_burgers_1d at
   --  dd9b081. The first-order cell update in flux form on a periodic grid (the reference has no update): the
   --  interface flux f(Q_L) + A-dQ, computed once per face.

   N_Cells     : constant := 256;
   Iface_Bound : constant := 5.0E12;
   Wide_Bound  : constant := 1.1E13;

   subtype Ratio is Long_Float range 0.0 .. 1.0;
   subtype Iface_Value is Long_Float range -Iface_Bound .. Iface_Bound;
   subtype Wide_Value is Long_Float range -Wide_Bound .. Wide_Bound;
   subtype Cell_Index is Positive range 1 .. N_Cells;
   type Cells is array (Cell_Index) of Field;
   type Wide_Cells is array (Cell_Index) of Wide_Value;

   --  periodic neighbours
   function Left_Of (I : Cell_Index) return Cell_Index is (if I = 1 then N_Cells else I - 1);
   function Right_Of (I : Cell_Index) return Cell_Index is (if I = N_Cells then 1 else I + 1);

   --  cell I as a state
   function Cell (Q : Cells; I : Cell_Index) return State is ((F => Q (I)));

   --  the interface flux between QL and QR: f(QL) + A-dQ, the one value both cells use
   function Interface_Flux (QL : State; QR : State) return Iface_Value is
     (X_Flux (QL).Flux_F + X_Left_Fluctuation (QL, QR).Fluctuation_F);

   --  the flux leaving cell I through its right face, and entering it through its left face
   function Out_Flux (Q : Cells; I : Cell_Index) return Iface_Value is
     (Interface_Flux (Cell (Q, I), Cell (Q, Right_Of (I))));
   function In_Flux (Q : Cells; I : Cell_Index) return Iface_Value is
     (Interface_Flux (Cell (Q, Left_Of (I)), Cell (Q, I)));

   --  Q_i - nu * (F_{i+1/2} - F_{i-1/2}), nu = dt / dx
   function Update_Cell (Q : Cells; Nu : Ratio; I : Cell_Index) return Wide_Value is
     (Q (I) - Nu * (Out_Flux (Q, I) - In_Flux (Q, I)));

   --  one step of every cell, periodic
   function Step (Q : Cells; Nu : Ratio) return Wide_Cells is
     ((for I in Cell_Index => Update_Cell (Q, Nu, I)));

   --  THEOREM: the value cell I sends through its right face is the value cell I+1 receives
   function Interface_Shared (Q : Cells; I : Cell_Index) return Boolean is
     (Out_Flux (Q, I) = In_Flux (Q, Right_Of (I)))
     with Post => Interface_Shared'Result;

   --  the grid whose every cell holds U
   function Uniform (U : State) return Cells is ((for I in Cell_Index => U.F));

   --  every cell of the uniform grid is U, and so are its neighbours
   function Cells_Are_U (U : State; I : Cell_Index) return Boolean is
     (Cell (Uniform (U), I) = U and then Cell (Uniform (U), Left_Of (I)) = U and then Cell (Uniform (U), Right_Of (I)) = U)
     with Post => Cells_Are_U'Result;

   --  both face fluxes of a cell of the uniform grid are the flux of U
   function Faces_Are_F (U : State; I : Cell_Index) return Boolean is
     (Cells_Are_U (U, I) and then Left_Consistent (U)
      and then Interface_Flux (U, U) = X_Flux (U).Flux_F
      and then Out_Flux (Uniform (U), I) = Interface_Flux (U, U)
      and then In_Flux (Uniform (U), I) = Interface_Flux (U, U)
      and then Out_Flux (Uniform (U), I) - In_Flux (Uniform (U), I) = 0.0)
     with Post => Faces_Are_F'Result;

   --  THEOREM: one cell of a constant state is unchanged by the update, exactly
   function Cell_Fixed (U : State; Nu : Ratio; I : Cell_Index) return Boolean is
     (Faces_Are_F (U, I)
      and then Interface_Flux (U, U) = X_Flux (U).Flux_F
      and then Update_Cell (Uniform (U), Nu, I) = U.F)
     with Post => Cell_Fixed'Result;

   --  THEOREM: a constant state is a fixed point of the step, exactly
   function Constant_Fixed (U : State; Nu : Ratio) return Boolean is
     (for all I in Cell_Index => Cell_Fixed (U, Nu, I) and then Step (Uniform (U), Nu) (I) = U.F)
     with Post => Constant_Fixed'Result;

end Lanyon_Burg1d_Update_Pkg;
