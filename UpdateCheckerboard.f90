Subroutine Update_checkerboard(Nl_aux)

    Use set_precision
    Use configuration 
    Use potential , Only     : rangeP, Npart_types
    Use rundata  , Only     : hmax   
    Use Cell      
!
    Implicit none
!
    Integer       ::   NCell, JCell, KCell, LCell, I, J, K, L, M, N, NxNy
    Integer       ::   NN, Nx, Indx_Nat, counter, indx_x, indx_y, indx_z
    Integer       ::   Ix, Iy, Iz, counter2, Indice, Indice2
    Integer       ::   Imap, Imap2, Indice1, Indice_Q
    Integer       ::   Ix_aux,IY_aux,IZ_aux, Indx_Cell
    Integer, Dimension(ndim)                ::  Indx
    Integer                                 ::  ICell, nat
    Integer                                 ::  Indx_x0, Indx_y0, Indx_z0
    Integer                                 ::  Nl_aux(0:Ndim-1)
    Real(wp)      ::   XCell, YCell, ZCell, W_RU_min
    Real(wp)      ::   CosAB, SinAB, XI, YI, ZI, wv, rangeP_min

    
    ! divide the simulation box into sets of square cells
    ! data that needs to be provided:
    !   w:  width of the cell
    !   box:  dimensions of the unit cell

    Nl_Cell(:) = Nl_aux(:)
    W_RU(:) =  1.d0 / dble(Nl_Cell (:))
    W(:) = w_RU(:) * Box(:)
    disp_max(:)=0.4*W_RU(:)

    W_RU_min=minval(W_RU)
    if(hmax > W_RU_min/2.0) then
          hmax= W_RU_min/2.0
     endif

    Ntot_Cell=Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)
    if(Ntot_Cell.gt.NCell_max) Then
        print*, '*****************************'
        print*, '*****************************'
        print*, ''
        print*, ' You exceeded the maximum number of Cells',NCell_max
        print*, ' Please, try  to start a simulation with a volume'
        print*, ' closer to the equilibrium value'
        print*, ''
        print*, '*****************************'
        print*, '*****************************'
        stop
    Endif


    If(ndim == 2) Then ! Initially I choose Cells sequentially
      Nset_cell=4
      Nsubset_cell=Ntot_Cell/Nset_cell
      Allocate(CB_Set(0:Nset_cell-1) )
      Allocate(Offset_CB(0:1))
      Do I=1,Nset_cell
         CB_Set(I)=I
      Enddo
    ElseIf(ndim == 3) Then
      Nset_cell=8
      Nsubset_cell=Ntot_Cell/Nset_cell

      Deallocate(List_CB)
      Allocate(List_CB(0:Ntot_Cell-1))
      counter2=0
      Do I=1,Nset_cell   ! defino un List_CB 
         CB_Set(I)=I
         Indx_z=Indx_Set(I,3)
         DO Iz=0,Nl_Cell(2)-1,2
            Indx_z = Indx_z +Iz
            Indx_y=Indx_Set(I,2)
            DO Iy=0,Nl_Cell(1)-1,2
               Indx_y = Indx_y + Iy
               Indx_x=Indx_Set(I,1)
               DO Ix=0,Nl_Cell(0)-1,2
                  Indx_x = Indx_x + Ix
                  List_CB(counter2)= ICell( Indx_x, Indx_y, Indx_z )
                  counter2=counter2+1
               EndDo
            EndDo
         EndDo
      Enddo

    Endif
     
    !  make a list of the atoms in each cell list_cell(i=1,Ntot_cell)
    Nat_cell_max=int(sqrt(2.d0)*W(0)*W(1)*W(2))+8 ! at close packing

    ! defino el mapa de celdas vecinas

    Deallocate( Map_Cell,DR_Cell)
    Allocate( Map_Cell(0:Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26-1))
    Allocate( DR_Cell(0:Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26*3-1))
    Do IZ = 0, Nl_Cell(2)-1
        Do IY = 0, Nl_Cell(1)-1
           Do Ix = 0, Nl_Cell(0)-1
              Indx_Cell=Icell(Ix,Iy,Iz)
              call Get_Cell_Indices( Indx_Cell , IX_aux, IY_aux, IZ_aux)
              Imap = ICell(Ix, Iy, Iz) * 26
              Imap2 = Imap*3

              Indice1= ICell ( Ix-1, Iy-1, Iz-1)
              Map_Cell ( Imap )   = Indice1
              DR_Cell( Imap2 )   = - W_RU(0)
              DR_Cell( Imap2+1 ) = - W_RU(1)
              DR_Cell( Imap2+2 ) = - W_RU(2)

              Indice1= ICell ( Ix  , Iy-1, Iz-1)
              Map_Cell ( Imap+1 ) = Indice1
              DR_Cell( Imap2+3 ) = 0.d0
              DR_Cell( Imap2+4 ) = -W_RU(1)
              DR_Cell( Imap2+5 ) = -W_RU(2)

              Indice1= ICell ( Ix+1, Iy-1, Iz-1)
              Map_Cell ( Imap+2 ) = Indice1
              DR_Cell( Imap2+6 ) = + W_RU(0)
              DR_Cell( Imap2+7 ) = - W_RU(1)
              DR_Cell( Imap2+8 ) = - W_RU(2)

              Indice1= ICell ( Ix-1, Iy  , Iz-1)
              Map_Cell ( Imap+3 ) = Indice1
              DR_Cell( Imap2+9 ) = - W_RU(0)
              DR_Cell( Imap2+10 ) = 0.00
              DR_Cell( Imap2+11 ) = - W_RU(2)

              Indice1= ICell ( Ix  , Iy  , Iz-1)
              Map_Cell ( Imap+4 ) = Indice1
              DR_Cell( Imap2+12 ) = 0.d0
              DR_Cell( Imap2+13 ) = 0.d0
              DR_Cell( Imap2+14 ) = - W_RU(2)

              Indice1= ICell ( Ix+1, Iy  , Iz-1)
              Map_Cell ( Imap+5 ) = Indice1
              DR_Cell( Imap2+15 ) = + W_RU(0)
              DR_Cell( Imap2+16 ) = 0.d0      
              DR_Cell( Imap2+17 ) = - W_RU(2)

              Indice1= ICell ( Ix-1, Iy+1, Iz-1)
              Map_Cell ( Imap+6 ) = Indice1
              DR_Cell( Imap2+18 ) = - W_RU(0)
              DR_Cell( Imap2+19 ) = + W_RU(1)
              DR_Cell( Imap2+20 ) = - W_RU(2)

              Indice1= ICell ( Ix  , Iy+1, Iz-1)
              Map_Cell ( Imap+7 ) = Indice1
              DR_Cell( Imap2+21 ) = 0.d0      
              DR_Cell( Imap2+22 ) = + W_RU(1)
              DR_Cell( Imap2+23 ) = - W_RU(2)

              Indice1= ICell ( Ix+1, Iy+1, Iz-1)
              Map_Cell ( Imap+8 ) = Indice1
              DR_Cell( Imap2+24 ) = + W_RU(0)
              DR_Cell( Imap2+25 ) = + W_RU(1)
              DR_Cell( Imap2+26 ) = - W_RU(2)

              Indice1= ICell ( Ix-1, Iy-1, Iz  )
              Map_Cell ( Imap+9 ) = Indice1
              DR_Cell( Imap2+27 ) = - W_RU(0)
              DR_Cell( Imap2+28 ) = - W_RU(1)
              DR_Cell( Imap2+29 ) = 0.d0     

              Indice1= ICell ( Ix  , Iy-1, Iz  )
              Map_Cell ( Imap+10) = Indice1
              DR_Cell( Imap2+30 ) = 0.d0      
              DR_Cell( Imap2+31 ) = - W_RU(1)
              DR_Cell( Imap2+32 ) = 0.d0      

              Indice1= ICell ( Ix+1, Iy-1, Iz  )
              Map_Cell ( Imap+11) = Indice1
              DR_Cell( Imap2+33 ) = + W_RU(0)
              DR_Cell( Imap2+34 ) = - W_RU(1)
              DR_Cell( Imap2+35 ) = 0.d0      

              Indice1= ICell ( Ix-1, Iy  , Iz  )
              Map_Cell ( Imap+12) = Indice1
              DR_Cell( Imap2+36 ) = - W_RU(0)
              DR_Cell( Imap2+37 ) = 0.d0      
              DR_Cell( Imap2+38 ) = 0.d0       

              Indice1= ICell ( Ix+1, Iy  , Iz  )
              Map_Cell ( Imap+13) = Indice1
              DR_Cell( Imap2+39 ) = + W_RU(0)
              DR_Cell( Imap2+40 ) = 0.d0      
              DR_Cell( Imap2+41 ) = 0.d0     

              Indice1= ICell ( Ix-1, Iy+1, Iz  )
              Map_Cell ( Imap+14) = Indice1
              DR_Cell( Imap2+42 ) = - W_RU(0)
              DR_Cell( Imap2+43 ) = + W_RU(1)
              DR_Cell( Imap2+44 ) = 0.d0     

              Indice1= ICell ( Ix  , Iy+1, Iz  )
              Map_Cell ( Imap+15) = Indice1
              DR_Cell( Imap2+45 ) = 0.d0      
              DR_Cell( Imap2+46 ) = + W_RU(1)
              DR_Cell( Imap2+47 ) = 0.d0      

              Indice1= ICell ( Ix+1, Iy+1, Iz  )
              Map_Cell ( Imap+16) = Indice1
              DR_Cell( Imap2+48 ) = + W_RU(0)
              DR_Cell( Imap2+49 ) = + W_RU(1)
              DR_Cell( Imap2+50 ) = 0.d0      

              Indice1= ICell ( Ix-1, Iy-1, Iz+1)
              Map_Cell ( Imap+17) = Indice1
              DR_Cell( Imap2+51 ) = - W_RU(0)
              DR_Cell( Imap2+52 ) = - W_RU(1)
              DR_Cell( Imap2+53 ) = + W_RU(2)

              Indice1= ICell ( Ix  , Iy-1, Iz+1)
              Map_Cell ( Imap+18) = Indice1
              DR_Cell( Imap2+54 ) = 0.d0     
              DR_Cell( Imap2+55 ) = - W_RU(1)
              DR_Cell( Imap2+56 ) = + W_RU(2)

              Indice1= ICell ( Ix+1, Iy-1, Iz+1)
              Map_Cell ( Imap+19) = Indice1
              DR_Cell( Imap2+57 ) = + W_RU(0)
              DR_Cell( Imap2+58 ) = - W_RU(1)
              DR_Cell( Imap2+59 ) = + W_RU(2)

              Indice1= ICell ( Ix-1, Iy  , Iz+1)
              Map_Cell ( Imap+20) = Indice1
              DR_Cell( Imap2+60 ) = - W_RU(0)
              DR_Cell( Imap2+61 ) = 0.d0      
              DR_Cell( Imap2+62 ) = + W_RU(2)

              Indice1= ICell ( Ix  , Iy  , Iz+1)
              Map_Cell ( Imap+21) = Indice1
              DR_Cell( Imap2+63 ) = 0.d0       
              DR_Cell( Imap2+64 ) = 0.d0       
              DR_Cell( Imap2+65 ) = + W_RU(2)

              Indice1= ICell ( Ix+1, Iy  , Iz+1)
              Map_Cell ( Imap+22) = Indice1
              DR_Cell( Imap2+66 ) = + W_RU(0)
              DR_Cell( Imap2+67 ) = 0.d0     
              DR_Cell( Imap2+68 ) = + W_RU(2)

              Indice1= ICell ( Ix-1, Iy+1, Iz+1)
              Map_Cell ( Imap+23) = Indice1
              DR_Cell( Imap2+69 ) = - W_RU(0)
              DR_Cell( Imap2+70 ) = + W_RU(1)
              DR_Cell( Imap2+71 ) = + W_RU(2)

              Indice1= ICell ( Ix  , Iy+1, Iz+1)
              Map_Cell ( Imap+24) = Indice1
              DR_Cell( Imap2+72 ) = 0.d0      
              DR_Cell( Imap2+73 ) = + W_RU(1)
              DR_Cell( Imap2+74 ) = + W_RU(2)

              Indice1= ICell ( Ix+1, Iy+1, Iz+1)
              Map_Cell ( Imap+25) = Indice1
              DR_Cell( Imap2+75 ) = + W_RU(0)
              DR_Cell( Imap2+76 ) = + W_RU(1)
              DR_Cell( Imap2+77 ) = + W_RU(2)
           EndDo
        EndDo
    EndDo


    Deallocate (List_Cell, Nat_Cell)
    Deallocate (Sphere_Cell)
    Allocate( List_Cell (0: Nl_cell(0)*Nl_Cell(1)*Nl_Cell(2)*Nat_cell_max*2-1) )

    Allocate( Nat_Cell (0:(Nl_cell(0)*Nl_Cell(1)*Nl_Cell(2))-1 ) )
    Allocate( Sphere_Cell (0:(Nl_cell(0)*Nl_Cell(1)*Nl_Cell(2)*Nat_Cell_max*(Ndim+4))-1 ) )

    Nat_Cell = 0
    List_Cell = 0
    Sphere_Cell = 0.d0

    Offset_CB = (/ 0., 0., 0. /)

    Do I=0,Natoms-1
       
       Indice=I*Ndim
       XI = R(Indice) 
       YI = R(Indice+1) 
       ZI = R(Indice+2) 

       If( XI < 0 ) XI = XI + INT(XI) + 1.d0 
       If( XI > 1 ) XI = XI - INT(XI) 

       If( YI < 0 ) YI = YI + INT(YI) + 1.d0 
       If( YI > 1 ) YI = YI - INT(YI) 

       If( ZI < 0 ) ZI = ZI + INT(ZI) + 1.d0 
       If( ZI > 1 ) ZI = ZI - INT(ZI) 

       Indx(1) = int( XI / w_RU(0) ) 
       Indx(2) = int( YI / w_RU(1) ) 
       Indx(3) = int( ZI / w_RU(2) ) 
 
       XCell = Indx(1)  * w_RU(0) + w_RU(0)*0.5d0!+ Offset_CB(0)
       YCell = Indx(2)  * w_RU(1) + w_RU(1)*0.5d0! Offset_CB(1) 
       ZCell = Indx(3)  * w_RU(2) + w_RU(2)*0.5d0 !Offset_CB(2) 

       Indice= ICell ( Indx(1), Indx(2), Indx(3))
       Indx_Nat = Nat_Cell( Indice ) + 1


       If(Indx_Nat > Nat_cell_max) Then
           print*, '*************************************************'
           print*, '*************************************************'
           print*, ''
           print*, '     ERROR in Initialize_neighbor'
           print*, '    Indx_Nat larger than Nat_cell_Max', Indx_Nat,Nat_Cell_Max
           print*, '    Total number of atoms', Natoms
           print*, '    Range of interactions', RangeP
           print*, '    Cell dimensions    ', W(:)
           print*, ''
           print*, '*************************************************'
           print*, '*************************************************'
           print*, '*************************************************'
           stop
       Endif

       Nat_Cell( Indice ) = Indx_Nat

       Indice2 = Indice*Nat_cell_max*(Ndim+4)+(Indx_Nat-1)*(Ndim+4)

       Indice_Q = I*4
       Sphere_Cell(Indice2) =  XI - XCell 
       Sphere_Cell(Indice2 + 1 ) = YI - YCell 
       Sphere_Cell(Indice2 + 2 ) = ZI - ZCell 
       Sphere_Cell(Indice2 + 3 ) = Q(Indice_Q)
       Sphere_Cell(Indice2 + 4 ) = Q(Indice_Q +1 )
       Sphere_Cell(Indice2 + 5 ) = Q(Indice_Q +2 )
       Sphere_Cell(Indice2 + 6 ) = Q(Indice_Q +3 )

       Indice2 = Indice*Nat_cell_max*2+(Indx_Nat-1)*2

       List_Cell(Indice2     ) = i   ! guardo el indice del atomo 
       List_Cell(Indice2 + 1 ) = itype(i)   ! y su tipo

    EndDo


End Subroutine Update_checkerboard
