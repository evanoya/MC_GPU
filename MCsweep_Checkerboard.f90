Subroutine MCsweep()

     Use set_precision
     Use potential
     Use Cell
     Use rundata
     Use utils
     Use configuration
     Use properties     
     Integer    :: I, J
     Integer    :: Indx, Nat_aux
     Integer, Dimension(0:Ndim-1)        :: Offset

     Do I=1,Nset_Cell

        Call Random_number(rnd)
        Indx = int(rnd*8.0)+1
        if(Indx.gt.8) Indx=8
        if(Indx.lt.1) Indx=1
        Offset(0)=Indx_Set(Indx,1)
        Offset(1)=Indx_Set(Indx,2)
        Offset(2)=Indx_Set(Indx,3)
        
        Call subsweepgpu(bool_tor, seed, ntrans_ac, nrot_ac, Ndim, Nmove, &
          & Nat_Cell_Max, Ntot_Cell, Nsubset_Cell,  &
          & Nrb_sites_max, Ntor_Max, Npart_types, & 
          & Nrb_sites, patch, ref_tor_vec, sigma_jon, sigma_tor_jon, sigma_LJ, rangeP, &
          & VL0, xop, &
          & Ntor_angle, tor_angle, Vpot_Matrix, & 
          & h, beta, hmax, omax, Sphere_Cell, List_Cell, &
          & Map_Cell, Nl_Cell, Nat_Cell, DR_Cell, W_RU, Offset, En_tot)

     EndDo
     Do I=0, Ntot_Cell-1
         If( Nat_Cell(i) > 0 ) Then
             no_mc_real = no_mc_real +1 
         EndIf
     EndDo

End Subroutine MCsweep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
Subroutine rng_shuffle(N,vec)

     Use set_precision
     Implicit None 
     Integer                 :: I, J
     Integer, Intent(IN)     :: N
     Integer, Intent(INOUT) :: vec(N,3)
     Integer                 :: aux(3)
     Real(wp)                :: rndp    
     

     Do I=N,2,-1

          Call Random_number(rndp)
          J=Int(i*rndp)+1
          aux(:)=vec(J,:)
          vec(J,:)=vec(I,:)
          vec(I,:)=aux(:)
     EndDo

End Subroutine rng_shuffle
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
Subroutine shuffle_indx(N,indx,vec1,vec2)

     Use set_precision
     Use Configuration, Only  : ndim
     Implicit None 
     Integer                 :: I, J, iaux1
     Integer, Intent(IN)     :: N
     Integer, Intent(INOUT)  :: indx(1:N)
     Real(wp)                :: rndp, vec1(1:N,1:ndim), Vec2(1:N,1:4), aux    
     

     Do I=N,2,-1

          Call Random_number(rndp)
          J=Int(i*rndp)+1
!
          iaux1=indx(J)
          indx(J)=indx(I)
          indx(I)=iaux1

          aux=vec1(J,1)
          vec1(J,1)=vec1(I,1)
          vec1(I,1)=aux

          aux=vec1(J,2)
          vec1(J,2)=vec1(I,2)
          vec1(I,2)=aux

          aux=vec1(J,3)
          vec1(J,3)=vec1(I,3)
          vec1(I,3)=aux

          aux=vec2(J,1)
          vec2(J,1)=vec2(I,1)
          vec2(I,1)=aux

          aux=vec2(J,2)
          vec2(J,2)=vec2(I,2)
          vec2(I,2)=aux

          aux=vec2(J,3)
          vec2(J,3)=vec2(I,3)
          vec2(I,3)=aux

          aux=vec2(J,4)
          vec2(J,4)=vec2(I,4)
          vec2(I,4)=aux
!
     EndDo

End Subroutine shuffle_indx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Convert_CB_RU

   Use set_precision
   Use configuration
   Use cell            
   Implicit None
   Integer      :: J0, J1, J2, J3, Nat_JCell
   Integer      :: indx, Indice, Indice2, ICell
   Integer      :: Indx2, IndxQ, Ncount, Indice3, indx3
   Real(wp)     :: Xcell, Ycell, Zcell

   Ncount = 0
   Do J3 = 0, Nl_Cell(2)-1
     Do J2 = 0, Nl_Cell(1)-1
       Nx: Do J1 = 0, Nl_Cell(0)-1

           Indice = ICell (J1, J2, J3)
           Nat_JCell = Nat_Cell( Indice )
           If(Nat_JCell == 0) CYCLE Nx
           Indice2 = Indice*Nat_cell_max*(Ndim+4)
           Indice3 = Indice*Nat_cell_max*2
           Xcell = J1 * w_RU(0) + w_RU(0)/2.d0 !+ Offset_CB(0)
           Ycell = J2 * w_RU(1) + w_RU(1)/2.d0 !+ Offset_CB(1)
           Zcell = J3 * w_RU(2) + w_RU(2)/2.d0 !+ Offset_CB(2)


           Do J0=0,Nat_JCell-1
              indx= Indice2 + J0*(Ndim+4)
              indx3= Indice3 + J0*2
              Indx2 = List_Cell(indx3)* Ndim
              R(indx2) = Sphere_Cell(indx) + Xcell
              R(indx2+1) = Sphere_Cell(indx+1) + Ycell
              R(indx2+2) = Sphere_Cell(indx+2) + Zcell
              IndxQ = List_Cell(indx3)* 4
              Q(IndxQ ) = Sphere_Cell(indx+3)
              Q(IndxQ+1 ) = Sphere_Cell(indx+4)
              Q(IndxQ+2 ) = Sphere_Cell(indx+5)
              Q(IndxQ+3 ) = Sphere_Cell(indx+6)
              Ncount = Ncount + 1
           EndDo

       EndDo Nx
     EndDo
   EndDo

End Subroutine Convert_CB_RU

