
Subroutine Volume_Move

   Use set_precision
   Use rundata
   Use configuration
   Use potential
   Use utils
   Use properties
   Use Cell       
   Implicit None

   Integer                         :: I, J, istate
   Integer, Dimension(0:ndim-1)    :: Nl_aux, Nl_Cell_old      
   Integer                         :: Indx_Cell, indx_atm_aux, indx_atm_aux2, Nat_Cell_ijk
   Integer                         :: jmol, jmol2
   Real(wp)                            :: En_old, Vol_old, Vol_new, Arg, Prob, rnd
   Real(wp)                            :: wv, rangeP_min
   Real(wp), Dimension(0:ndim*ndim-1)  :: h_old
   Real(wp), Dimension(0:ndim-1)       :: Box_old
   Real(wp), Dimension(0:ndim-1)       :: W_RU_old, W_RU_aux, W_aux, W_old
   Logical :: overlap, update_cell
 
   En_old= En_tot
   Vol_old = box(0)*box(1)*box(2)
   Box_old =Box
   W_RU_old = W_RU
   W_old = W
   Nl_Cell_old = Nl_Cell
   Call Random_number(rnd)
   IF(iscale.eq.0) Then
     box(0)=box(0)+(rnd-0.5d00)*vmax(0)
     box(1)=box(0)
     box(2)=box(0)
   ElseIF(iscale.eq.1) Then
     box(0)=box(0)+(rnd-0.5d00)*vmax(0)
     box(1)=box(0)
     Call Random_number(rnd)
     box(2)=box(2)+(rnd-0.5d00)*vmax(2)
   ElseIF(iscale.eq.2) Then
     box(0)=box(0)+(rnd-0.5d00)*vmax(0)
     Call Random_number(rnd)
     box(1)=box(1)+(rnd-0.5d00)*vmax(1)
     Call Random_number(rnd)
     box(2)=box(2)+(rnd-0.5d00)*vmax(2)
   EndIf
   H_old= H

   H(0)=H(0)*Box(0)/Box_old(0)
   H(3)=H(3)*Box(0)/Box_old(0)
   H(6)=H(6)*Box(0)/Box_old(0)

   H(1)=H(1)*Box(1)/Box_old(1)
   H(4)=H(4)*Box(1)/Box_old(1)
   H(7)=H(7)*Box(1)/Box_old(1)

   H(2)=H(2)*Box(2)/Box_old(2)
   H(5)=H(5)*Box(2)/Box_old(2)
   H(8)=H(8)*Box(2)/Box_old(2)

   Vol_new = box(0)*box(1)*box(2)
   rangeP_min=minval(rangeP)
   wv = 1.001d0*rangeP_min
   w_aux(:) = (/ wv, wv, wv /)
   W_RU_aux = W_aux /box 
   Nl_aux(:) = 1.d0 / W_RU_aux (:)
   
   If (mod(Nl_aux(0),2)  == 1)  Nl_aux(0) = Nl_aux(0) - 1
   If(Nl_aux(0) > 64) Nl_aux(0) = 64*((Nl_aux(0)/64))  ! round up to nearest multiple of 64

   If (mod(Nl_aux(1),2)  == 1)  Nl_aux(1) = Nl_aux(1) - 1
   If(Nl_aux(1) > 64) Nl_aux(1) = 64*((Nl_aux(1)/64))  ! round up to nearest multiple of 64

   If (mod(Nl_aux(2),2)  == 1)  Nl_aux(2) = Nl_aux(2) - 1
   If(Nl_aux(2) > 64) Nl_aux(2) = 64*((Nl_aux(2)/64))  ! round up to nearest multiple of 64

   if(Nl_aux(0)==2 .or. NL_aux(1) == 2 .or. Nl_aux(2)==2) Then
       print*, 'your system is too small, sorry!'
       stop
   endif

   update_cell=.false.
   istate=1
   If((Nl_aux(0).ne.Nl_Cell(0)).or.(Nl_aux(1).ne.Nl_Cell(1)).or.(Nl_aux(2).ne.Nl_Cell(2))) Then
       update_cell=.true.
       W_RU =W_RU_old
       Call Convert_CB_RU  ! aqui uso Nl_viejo y W_RU viejo
       Call Update_checkerboard(Nl_aux)
       istate=3
   Else
       W(:)=W_RU(:)*Box(:)
   EndIf 
   Call CB_system_energy_gpu(istate)

   Arg= -beta*((En_tot-En_old)+pres*(Vol_new-Vol_old)-dble(Natoms)*log(Vol_new/Vol_old)/beta)
   Prob = min (1.d0, exp(Arg))
   Call Random_number(rnd)
   If(rnd < Prob ) Then ! accept
        nvol_accept=nvol_accept+1 
        if(update_cell) then
             call energyinitialize(Ndim, Nat_cell_max, Ntot_Cell, &
                     &   h, Sphere_Cell, List_Cell, &
                     &   Map_Cell, Nl_Cell, Nat_Cell, DR_Cell,W_RU)
         endif
   Else   ! reject move
       If(update_cell) Then  
         Call Convert_CB_RU  
         h = h_old
         Box = Box_old
         Nl_Cell = Nl_Cell_old
         Call Update_checkerboard(Nl_Cell)
         call energyinitialize(Ndim, Nat_cell_max, Ntot_Cell, &
                     &   h, Sphere_Cell, List_Cell, &
                     &   Map_Cell, Nl_Cell, Nat_Cell, DR_Cell, W_RU)
       Else
         h = h_old
         Box = Box_old
         W_RU = W_RU_old
         W = W_old
         Nl_Cell = Nl_Cell_old
       EndIf
       En_tot=En_old
   Endif

End Subroutine Volume_Move
