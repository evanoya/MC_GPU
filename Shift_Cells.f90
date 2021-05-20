!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
Subroutine Shift_Cells(ivec,disp)

   Use set_precision
   Use Cell
   Use Configuration
   Implicit None
   Integer     :: ivec(0:ndim-1)
   Real(wp)    :: disp(0:ndim-1)



   call shiftcellsgpu(Ndim,Nat_Cell_Max, Ntot_Cell, Nl_Cell, Sphere_Cell, List_Cell, Nat_Cell, W_RU, ivec, disp)


End Subroutine Shift_Cells
