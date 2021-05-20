Subroutine Initialize_checkerboard

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
    Real(wp)      ::   AB_P(ndim)
    Real(wp)      ::   XCell, YCell, ZCell 
    Real(wp)      ::   CosAB, SinAB, AB_P_N, XI, YI, ZI, wv, rangeP_min

    
    rangeP_min=minval(rangeP)
    print*, 'rangeP_min',rangeP_min,minval(rangeP)

    wv = 1.001d0*rangeP_min
    w(:) = (/ wv, wv, wv /)

    print*, 'Entering Initialize Checkerboard'
    print*, 'w=',w

    Call Get_edges_angles
   
    If(abs(bangles(0))>0.001) Then
        print*, 'non cubic box not implemented yet'
        stop
    EndIf
    If(abs(bangles(1))>0.001) Then
        print*, 'non cubic box not implemented yet'
        stop
    EndIf
    If(abs(bangles(2))>0.001) Then
        print*, 'non cubic box not implemented yet'
        stop
    EndIf

    W_RU(:) = w(:) / Box(:) 

    print*, 'W_RU' ,W(:),W_RU(:)

    Nl_Cell(0) = 1. / W_RU(0)
    If (mod(Nl_Cell(0),2)  == 1)  Nl_Cell(0) = Nl_Cell(0) - 1
    If(Nl_Cell(0) > 64) Nl_cell(0) = 64*((Nl_cell(0)/64))  ! round up to nearest multiple of 64

    Nl_Cell(1) = 1. / W_RU(1)
    If (mod(Nl_Cell(1),2)  == 1)  Nl_Cell(1) = Nl_Cell(1) - 1
    If(Nl_Cell(1) > 64) Nl_cell(1) = 64*((Nl_cell(1)/64))  ! round up to nearest multiple of 64

    Nl_Cell(2) = 1. / W_RU(2)
    If (mod(Nl_Cell(2),2)  == 1)  Nl_Cell(2) = Nl_Cell(2) - 1
    If(Nl_Cell(2) > 64) Nl_cell(2) = 64*((Nl_cell(2)/64))  ! round up to nearest multiple of 64
    print*, 'Nl_Cell',NL_Cell
    if(Nl_Cell(0)==2 .or. NL_Cell(1) == 2 .or. Nl_Cell(2)==2) Then
        print*, 'your system is too small, sorry!'
        stop
    endif

    ! Correct the values of W and W_RU

    W_RU(:) =  1.d0 / dble(Nl_Cell (:))
    W(:) = w_RU(:) * Box(:) 
    disp_max(:)=0.4*W_RU(:)
    print*, 'W_RU 2' ,W(:),W_RU(:),W_RU(:)/2.d0
 
    if(hmax > W_RU(0)/2.0) then
          hmax= W_RU(0)/2.0
          print*, 'adjusting hmax=',hmax
     endif

    Ntot_Cell=Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)


    If(ndim == 2) Then ! Initially I choose Cells sequentially
      Nset_cell=4
      Nsubset_cell=Ntot_Cell/Nset_cell
      Allocate(CB_Set(0:Nset_cell-1) )
      Allocate(Offset_CB(0:1))
      Do I=1,Nset_cell
         CB_Set(I)=I
      Enddo
    ElseIf(ndim == 3) Then
      print*, 'defining the CB sets'
      ! I am assuming that:
      ! CB_Set  = 1--->  x impar  ,  y impar , z impar
      ! CB_Set  = 2--->  x impar  ,  y impar , z par
      ! CB_Set  = 3--->  x impar  ,  y par   , z impar
      ! CB_Set  = 4--->  x impar  ,  y par   , z par
      ! CB_Set  = 5--->  x par    ,  y impar , z impar
      ! CB_Set  = 6--->  x par    ,  y impar , z par
      ! CB_Set  = 7--->  x par    ,  y par   , z impar
      ! CB_Set  = 8--->  x par    ,  y par   , z par
      Nset_cell=8
      Nsubset_cell=Ntot_Cell/Nset_cell
      print*, 'Ntot_Cell=',Ntot_Cell
      print*, 'Nsubset_Cell=',Nsubset_Cell
      Allocate(CB_Set(Nset_cell))
      Allocate(Indx_Set(Nset_cell,3))
      Allocate(Offset_CB(0:2))
      Indx_Set(1,:)= (/ 0, 0, 0 /)
      Indx_Set(2,:)= (/ 0, 0, 1 /)
      Indx_Set(3,:)= (/ 0, 1, 0 /)
      Indx_Set(4,:)= (/ 0, 1, 1 /)
      Indx_Set(5,:)= (/ 1, 0, 0 /)
      Indx_Set(6,:)= (/ 1, 0, 1 /)
      Indx_Set(7,:)= (/ 1, 1, 0 /)
      Indx_Set(8,:)= (/ 1, 1, 1 /)


      !print*, 'List_CB'
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
     
    Nat_cell_max=int(sqrt(2.d0)*W(0)*W(1)*W(2))+8 ! at close packing
    print*, 'Nat_cell_max=',Nat_cell_max

    ! defino el mapa de celdas vecinas

    print*, 'Before allocate Map_Cell'
    Allocate( Map_Cell(0:Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26-1))
    print*, 'Before allocate DR_Cell',Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26
    Allocate( DR_Cell(0:Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26*3-1))
    print*, 'Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26*3-1',Nl_Cell(0)*Nl_Cell(1)*Nl_Cell(2)*26*3-1
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

   print*, 'here'
    print*, 'Before allocate List_Cell',Nl_cell(0)*Nl_Cell(1)*Nl_Cell(2)*Nat_cell_max
    Allocate( List_Cell (0: Nl_cell(0)*Nl_Cell(1)*Nl_Cell(2)*Nat_cell_max*2-1) )

    Allocate( Nat_Cell (0:(Nl_cell(0)*Nl_Cell(1)*Nl_Cell(2))-1 ) )
    print*, 'Before allocate Sphere_Cell'
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
 
       XCell = Indx(1)  * w_RU(0) + w_RU(0)/2.d0 + Offset_CB(0) 
       YCell = Indx(2)  * w_RU(1) + w_RU(1)/2.d0 + Offset_CB(1) 
       ZCell = Indx(3)  * w_RU(2) + w_RU(2)/2.d0 + Offset_CB(2) 

       Indice= ICell ( Indx(1), Indx(2), Indx(3))
       Indx_Nat = Nat_Cell( Indice ) + 1

       If(Indx_Nat > Nat_cell_max) Then
           print*, '*************************************************'
           print*, '*************************************************'
           print*, ''
           print*, '     ERROR in Initialize_neighbor'
           print*, '    Indx_Nat larger than Nat_cell_Max', Indx_Nat,Nat_Cell_Max
           print*, '    Total number of atoms', Natoms
           print*, '    Range of interactions', RangeP,W(:)
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
    print*, 'checkerboard initialized'


End Subroutine Initialize_checkerboard

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Function Cross (a,b)
!   
!   Use configuration, Only  : ndim 
!   Implicit None
!   Real(wp), Intent (in)    :: A(ndim), B(ndim)
!   Real(wp), Intent (out)   :: Cross(ndim)
!
!   If(ndim == 2) Then
!
!        print*, 'ERROR'
!        print*, 'Function Cross not defined for dimension',ndim
!        Stop
!        
!   ElseIf(ndim == 3) Then
!
!       Cross(1) = A(2)*B(3) - A(3)*B(2)
!       Cross(2) = A(3)*B(1) - A(1)*B(3)
!       Cross(3) = A(1)*B(2) - A(2)*B(1)
!
!   Else
!        print*, 'ERROR'
!        print*, 'Function Cross not defined for dimension',ndim
!        Stop
!   Endif
!
!End Function Cross 
    Subroutine Get_Cell_Indices( ICell , IX, IY, IZ)
    Use Cell      
    Implicit none
    Integer ICell, IX, IY, IZ

    Iz = int(Icell /Nl_Cell(0)/Nl_Cell(1))
    Iy = int((Icell - Nl_Cell(0)*Nl_Cell(1)*Iz) /Nl_Cell(0))
    If(Iz.gt.0.or.Iy.gt.0) then
          Ix =  mod( Icell, Nl_Cell(0)*Iy+Nl_cell(0)*Nl_Cell(1)*Iz) 
    else
          Ix=ICell
    Endif

    return
    End Subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer function icell(Ix, Iy, Iz)
    Use Cell, only : Nl_Cell
    Implicit none
    Integer ix,iy,iz

      ICell  =     MOD (IX+Nl_Cell(0), Nl_Cell(0)) + &
                 & MOD (IY+Nl_Cell(1), Nl_Cell(1)) * Nl_Cell(0) + &
                 & MOD (IZ+Nl_Cell(2), Nl_Cell(2)) * Nl_Cell(0)*Nl_Cell(1)
    Return
    End function 
