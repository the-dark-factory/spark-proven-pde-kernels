package Lanyon_Adv1d_Exact_Pkg with SPARK_Mode, Pure is

   --  After Lanyon AI's AdvectionDiffusion research (github.com/lanyonai/AdvectionDiffusion, which carries no
   --  licence): the SCHEME is ported from its mathematics; none of its text is copied.
   --  the first-order upwind scheme on a scaled integer
   --  grid. One interface flux per face, subtracted from one cell and added to the next, so the cell sum
   --  telescopes and is conserved EXACTLY on a periodic grid (brief level 3, the exact variant).

   --  the Courant number is C / Den, C in -Den .. Den; cell values are scaled integers
   Den        : constant := 2 ** 16;
   Cell_Bound : constant := 2 ** 40;
   Wide_Bound : constant := 3 * Cell_Bound;
   N_Cells    : constant := 256;

   subtype Cell_Value is Long_Long_Integer range -Cell_Bound .. Cell_Bound;
   subtype Wide_Value is Long_Long_Integer range -Wide_Bound .. Wide_Bound;
   subtype Courant is Long_Long_Integer range -Den .. Den;
   subtype Cell_Index is Positive range 1 .. N_Cells;
   subtype Cell_Count is Natural range 0 .. N_Cells;
   type Cells is array (Cell_Index) of Cell_Value;
   type Wide_Cells is array (Cell_Index) of Wide_Value;

   --  periodic neighbours
   function Left_Of (I : Cell_Index) return Cell_Index is
     (if I = 1 then N_Cells else I - 1);

   function Right_Of (I : Cell_Index) return Cell_Index is
     (if I = N_Cells then 1 else I + 1);

   --  floor (X / Den), for either sign of X
   function Floor_Div (X : Long_Long_Integer) return Cell_Value is
     (if X >= 0 then X / Den else -((Den - 1 - X) / Den))
     with Pre => X >= -(Den * Cell_Bound) and then X <= Den * Cell_Bound;

   --  the flux through the right face of cell I, taken from the upwind cell
   function Face_Flux (Q : Cells; C : Courant; I : Cell_Index) return Cell_Value is
     (if C >= 0 then Floor_Div (C * Q (I)) else Floor_Div (C * Q (Right_Of (I))));

   --  Q_i - F_{i+1/2} + F_{i-1/2}
   function Step_Cell (Q : Cells; C : Courant; I : Cell_Index) return Wide_Value is
     (Q (I) - Face_Flux (Q, C, I) + Face_Flux (Q, C, Left_Of (I)));

   --  one step of every cell
   function Step (Q : Cells; C : Courant) return Wide_Cells is
     ((for I in Cell_Index => Step_Cell (Q, C, I)));

   --  the grid in the wide type, unchanged
   function Widen (Q : Cells) return Wide_Cells is
     ((for I in Cell_Index => Q (I)));

   --  the sum of the first K cells
   function Sum (W : Wide_Cells; K : Cell_Count) return Long_Long_Integer is
     (if K = 0 then 0 else Sum (W, K - 1) + W (K))
     with Post               => abs Sum'Result <= Long_Long_Integer (K) * Wide_Bound,
          Subprogram_Variant => (Decreases => K);

   --  the prefix sums after the step telescope: the first K cells lose F_K and gain F_N (cell 1's inflow)
   function Telescope (Q : Cells; C : Courant; K : Cell_Index) return Boolean is
     (if K = 1 then True else Telescope (Q, C, K - 1))
     with Post               => Telescope'Result
                                  and then Sum (Step (Q, C), K)
                                    = Sum (Widen (Q), K) - Face_Flux (Q, C, K)
                                      + Face_Flux (Q, C, N_Cells),
          Subprogram_Variant => (Decreases => K);

   --  THEOREM (brief level 3, exact variant): the step conserves the cell sum exactly
   function Conserved (Q : Cells; C : Courant) return Boolean is
     (Telescope (Q, C, N_Cells)
      and then Sum (Step (Q, C), N_Cells) = Sum (Widen (Q), N_Cells))
     with Post => Conserved'Result;

end Lanyon_Adv1d_Exact_Pkg;
