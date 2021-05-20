
Subroutine Write_Input_File(istep)

   Use set_precision
   Use rundata
   Use configuration, Only : ndim
   Use potential
   Use utils
   Implicit None
   Integer :: istep, istep_ini_new
   Integer :: J, Indx, Indx_tor, Indx_tor2, I, Indx2, K, Indx1, M

   Open(21,file='input-restart.d')

   write(21,'(a20)') base

   write(21,*)  "Step_ini Step_fin Neq  Nmove  Nsave Nsave2 Nrestart"
   istep_ini_new=istep_fin-istep
   write(21,'(7I15)') istep,istep_fin,0, Nmove, Nsave, Nsave2, Nrestart

   write(21,*)  "Hmax Max_or vmax adjust_displ  temp_initial temp_final npt iscale pres seed"
   temp0=1./beta
   write(21,'(5F8.4,x,L6,x,2F12.5,x,L6,i3,x,F12.5,x,i10)') hmax,omax,vmax,displ_update,temp0,temp1,npt,iscale,pres,seed
   write(21,*) "(not implemented)"
   write(21,*) "XXXXXXXXXXXXXX"
   write(21,*) "imovie"
   write(21,'(L2)') imovie  

   write(21,*) "Model  (KF or LJG)"
   write(21,*) "LJG"
   write(21,*) "Parameters LJG model     "
   write(21,*) "cosmax ctorsiona_max deltaa  bool_torsion"
   write(21,'(3F8.4,3x,I2)') sigma_jon_aux, dsqrt(sigma_tor_jon/2.d0),rangepp, Bool_tor
   write(21,*) "(not implemented)"
   write(21,*) ""
   write(21,*) " Number  of particle types"
   Write(21,*) Npart_types

   Write(21,*) "Number of patches in each particle"
   Write(21,*) (Nrb_sites(J),J=0,Npart_types-1)
   Do J = 0, Npart_types-1

      Indx = J * Nrb_sites_max*Ndim
      Indx_tor = J * Nrb_sites_max
      Indx_tor2 = J * Nrb_sites_max*Ntor_max

      If( Nrb_sites(J) >= 2 ) Then ! particle J has at least two patches
         Do I = 0, Nrb_sites(J)-1

            Indx2 = Indx + I*Ndim

            Write(21,'(3F12.8)',advance='no') patch( Indx2 ),patch( Indx2+1 ),patch( Indx2+2 )

            If (Ntor_angle(Indx_tor+i) == 0 ) Then 

                Write(21,'(3F12.8,3x,I3)')      &
                   &  ref_tor_vec( Indx2 ), ref_tor_vec( Indx2+1 ), ref_tor_vec( Indx2+2 ), &
                   &   Ntor_angle( Indx_tor + I )
            Else
                Write(21,'(3F12.8,3x,I3)',advance='no')      &
                   &  ref_tor_vec( Indx2 ), ref_tor_vec( Indx2+1 ), ref_tor_vec( Indx2+2 ), &
                   &   Ntor_angle( Indx_tor + I )

                Do K = 0, Ntor_angle( Indx_tor+I )-1
                   Write(21,'(F12.6)',advance='no') tor_angle( Indx_tor2 +I*Ntor_max + K )
                EndDo
                Write(21,*)
            EndIf

         End Do
      ElseIf( Nrb_sites(J) == 1) Then
            Write(21,'(3F12.8)')  ( patch( Indx+K ), K = 0, Ndim -1 )
      Endif

   EndDo

   Write(21,*) "-vpot matrix --"

   Do I = 0, Npart_types -1

     Indx1= I * Npart_types * Nrb_sites_max *Nrb_sites_max

     Do J = 0, Nrb_sites(I) -1

       Indx2 = Indx1 +  J * (Nrb_sites_max*Npart_types) ! me situo en la fila correspondiente al parche j

         Do K = 0, Npart_types -1
            Do M = 0, Nrb_sites(K) -1

                 write(21,'(F5.2)',advance='no') Vpot_Matrix( Indx2+k*Nrb_sites_max+M )
            EndDo
         EndDo
         write(21,*) 

     EndDo
   EndDo

   ! Write the sigmas
   If(sigma_jon_aux.gt.0.d0) then
      ! Do nothing 
   Else

      Write(21,*) 'Sigma_ang' 
      Do I= 0, Npart_types-1
         Indx= I *Nrb_sites_max
         Do J = 0, Nrb_sites(i)-1
            write(21,'(F6.3)',advance='no') sigma_jon(indx+j)
         EndDo
         write(21,*)
      EndDo

      Write(21,*) 'Sigma_LJ'
      Do I=0,Npart_types-1
           Indx = I*Npart_types
           Do J = 0, Npart_types-1
              Write(21,'(F6.3)',advance='no') sigma_LJ(indx+j)
           EndDo
           Write(21,*) 
      EndDo


   EndIf

   Close(21)

End Subroutine Write_Input_File

Subroutine Read_Input_Data
 
   Use set_precision
   Use rundata
   Use configuration, Only : ndim
   Use potential 
   Use utils
   Implicit None
   Integer :: i, j, k, l, m, Indx, Indx1, Indx2, Indx3, IndxT1, IndxT2
   Integer :: IndxT3, itypei,ipatchi,itypej,ipatchj,num_tor, ios, Indx4
   Integer :: Indx_tor, Indx_tor2, i_fila, i_columna
   Real(wp) :: tor_aux(0:Ntor_max-1), norm, px, py, pz
   Real(wp) :: sigma_jon_ij, xopp
   real(wp), Allocatable :: Matrix(:,:)
   character*3  :: model

   Open(20,file='input.d')
   read(20,*) base 
   read(20,*) 
   print*, 'Frequency for averages', Nsave, Nsave2
   read(20,*) istep_ini,istep_fin,Neq, Nmove, Nsave, Nsave2, Nrestart
   print*, 'Initial and final MC cycles',istep_ini,istep_fin
   print*, 'Cycles for equilibration', Neq
   print*, 'Frequency for averages', Nsave, Nsave2, Nrestart
   read(20,*)
   read(20,*) hmax,omax,vmax(0),vmax(1),vmax(2),displ_update,temp0,temp1,npt,iscale,pres,seed
   print*, 'maximum translational displacement',hmax
   print*, 'maximum orientational change',omax
   print*, 'maximum displacement adjusted along the simulation?',displ_update
   print*, 'seed to rng',seed
   print*, 'temperature',temp0,temp1
   If(NpT) Then
       print*, 'Performing NpT simulation at pressure',pres
       print*, 'Maximum volume displacement',vmax
   Endif
   read(20,*)
   read(20,*)   ! leave this for furture extension to free energy calculations
   read(20,*)
   read(20,*)  imovie !  logical variable to decide whether to store the movie
   print*, 'imovie=', imovie !  logical variable to decide whether to store the movie
   read(20,*)
   read(20,*)  model
   read(20,*)
   print*, 'Modelo ',model
   if( model == 'K-F') then

      read(20,*)
      read(20,*) cosmax_KF, ctor_max_KF, delta_KF, Bool_tor

      print*, ''
      print*, 'Model parameters'
      print*, 'Opening angle',cosmax_KF,acos(cosmax_KF)*360./dospi
      print*, 'Interaction range',delta_KF
      rangeP=1.d0+delta_KF
      print*, 'range_KF=',rangeP

      print*, 'Torsional term =',Bool_tor
      If(Bool_tor.eq.1) print*, 'Torsional angle=',ctor_max_KF

   elseif( model == 'LJG') then

      read(20,*)
      read(20,*) sigma_jon_aux, sigma_tor_jon, rangepp, Bool_tor
      sigma_tor_jon = 2.d0*sigma_tor_jon*sigma_tor_jon;

   else
         
      print*, 'Model not yet implemented'

   endif

   read(20,*)  !leave this for future extension to AVB moves 
   read(20,*)
   read(20,*)
   read(20,*) Npart_types
   print*, 'Number of particle types',Npart_types

   Allocate( Nrb_sites(0:Npart_types)) 

   read(20,*)
   read(20,*) (Nrb_sites(J),J=0,Npart_types-1)

   Do i = 0, Npart_types-1
       print*, 'Particle type',i,'number of patches',Nrb_sites(i)
       If (Nrb_sites(i) > Nrb_sites_max) Then
            print*, ''
            print*, '***************************************'
            print*, '***************************************'
            print*, ''
            print*, 'ERROR in particle type',i
            print*, 'Number of patches',Nrb_sites(i)
            print*, 'exceeds the maximum number of patches',Nrb_sites_max
            print*, 'increase Nrb_sites_max in Definitions'
            print*, ''
            print*, '***************************************'
            print*, '***************************************'
            print*, ''
            stop
       Endif
   EndDo

   Allocate ( patch(0:Nrb_sites_max*Npart_types*ndim-1) )
   Allocate ( ref_tor_vec(0:Nrb_sites_max*Npart_types*ndim-1) )

   Allocate( Ntor_angle(0:Npart_types*Nrb_sites_max-1)) 
   Allocate( tor_angle(0:Npart_types*Nrb_sites_max*Ntor_max-1)) 

   Ntor_angle = 0
   tor_angle = 0.d0

   patch=0.d0

   Do J = 0, Npart_types-1

      Indx = J * Nrb_sites_max*Ndim
      Indx_tor = J * Nrb_sites_max
      Indx_tor2 = J * Nrb_sites_max*Ntor_max
      print*, '********************************'
      print*, '********************************'
      print*, 'particle j (number of patches)=',j, Nrb_sites(J)

      If( Nrb_sites(J) >= 2 ) Then ! particle J has at least two patches
         print*, 'particle type j 2',j,Nrb_sites(j)
         Do I = 0, Nrb_sites(J)-1

            Indx2 = Indx + I*Ndim

            Read(20,'(3F12.8)',advance='no')  px,py,pz
            norm= sqrt(px*px+py*py+pz*pz)
            patch( Indx2 ) =Px/norm
            patch( Indx2+1 ) =Py/norm
            patch( Indx2+2 ) =Pz/norm


            Read(20,'(3F12.8,3x,I3)',advance='no',IOSTAT=ios)      &
                   &  px, py, pz,    &
                   &   Ntor_angle( Indx_tor + I )
            norm= sqrt(px*px+py*py+pz*pz)
            ref_tor_vec( Indx2 ) =Px/norm
            ref_tor_vec( Indx2+1 ) =Py/norm
            ref_tor_vec( Indx2+2 ) =Pz/norm

            print*, 'i, Ntor=',i,Ntor_angle( Indx_tor + I ), ios
            If(ios.ge.1) Then
                Ntor_angle(Indx_tor+i) = 0
            Else

                Do K = 0, Ntor_angle( Indx_tor+I )-1 
                   Read(20,'(F12.6)',advance='no') tor_angle( Indx_tor2 +I*Ntor_max + K )
                EndDo
                Read(20,*) 
            EndIf

            print*, ''
            print*, ''
            print*, 'particle patch',j,i,(patch(Indx2+k),k=0,Ndim-1)
            print*, 'ref_vector',(ref_tor_vec(Indx2+k),k=0,Ndim-1)
            print*, 'Ntor=',Ntor_angle( Indx_tor + I )
            Do K = 0, Ntor_angle( Indx_tor+I )-1 
               Write(*,'(F12.6)') tor_angle( Indx_tor2 + I*Ntor_max + K )
            EndDo
            print*, ''

         End Do
      ElseIf( Nrb_sites(J) == 1) Then
            print*, 'particle type j 1',j,Nrb_sites(j)
            Read(20,'(3F12.8)')  ( patch( Indx+K ), K = 0, Ndim -1 )
            print*, 'patch=',  ( patch( Indx+K ), K = 0, Ndim -1 )
      Endif

   EndDo
   print*, 'after the double loop'

   If ( Maxval(Ntor_angle) > Ntor_max) Then
        print*, ''
        print*, '***************************************'
        print*, '***************************************'
        print*, ''
        print*, 'Error: maximum number of torsional angles allowed:',Ntor_max
        print*, 'increase in Definitions Ntor_max to ',Maxval(Ntor_angle)
        print*, ''
        print*, '***************************************'
        print*, '***************************************'
        print*, ''
        stop
   EndIf


   ! read in the interaction matrix

   Allocate ( Vpot_Matrix( 0 : Npart_types*Npart_types*Nrb_sites_max*Nrb_sites_max -1 ) )
   Vpot_Matrix = 0.

   Read(20,*) 

   print*, ''
   print*, ''
   print*, 'Vpot_Matrix(type1,patch1,type2,patch2)'
   Do I = 0, Npart_types -1

     Indx1= I * Npart_types * Nrb_sites_max *Nrb_sites_max ! me pongo en el indice de fila donde empieza molec i
     !print*, 'Particle I=',I+1

     Do J = 0, Nrb_sites(I) -1 

       Indx2 = Indx1 +  J * (Nrb_sites_max*Npart_types) ! me situo en la fila correspondiente al parche j
     !print*, 'Patch=',J+1

       !Do K = 0, Npart_types -1 

       !    Indx3 = Indx2 + K * Nrb_sites_max ! me muevo dentro de la fila para ver ia con particula de tipo K
    ! print*, 'Indx3=',Indx3, Nrb_sites(k)

           Read(20,*) ((Vpot_Matrix( Indx2+k*Nrb_sites_max+M ), M=0, Nrb_sites(K) -1), k=0, Npart_types-1)
           write(*,*) ((Vpot_Matrix( Indx2+k*Nrb_sites_max+M ), M=0, Nrb_sites(K) -1), k=0, Npart_types-1)

       !Do K = 0, Npart_types -1 
       !    Indx3 = Indx2 + K * Nrb_sites_max ! me muevo dentro de la fila para ver ia con particula de tipo K
       !    Do M = 0, Nrb_sites(K)-1 
       !       print*, ''
       !       write(*,'("Vpot_Matrix(",i3,",",i3,",",i3,",",i3,")=",f10.3)'),i,j,k,m,Vpot_Matrix( Indx3+M )
       !    EndDo
!
       !EndDo

     EndDo
   EndDo
   print*, ''
   print*, ''

   Allocate ( Matrix( Npart_types*Nrb_sites_max, Npart_types*Nrb_sites_max))


   print*, 'matrix dimensions',Npart_types*Nrb_sites_max

   ! aqui estaria bien meter un chequeo de que Vpot_matrix es simetrica
   Do I = 0, Npart_types*Npart_types*Nrb_sites_max*Nrb_sites_max -1
        
      i_fila=int( I / (Npart_types*Nrb_sites_max) ) + 1
      i_columna = mod(i,Npart_types*Nrb_sites_max) + 1
      matrix( i_fila, i_columna)=Vpot_Matrix( I) 

   EndDo

   Do I= 1, Npart_types*Nrb_sites_max
      Do J= 1, Npart_types*Nrb_sites_max
         If( matrix(i, j) .ne. matrix(j,i) ) Then
             print*, 'error, vpot_matrix not symmetric'
             print*, 'i,j',i,j
             print*, 'matrix(i, j)',matrix(i, j)
             print*, 'matrix(j, i ',matrix(j, i)
         EndIf
      EndDo
   EndDo

   ! Read the sigmas
   Allocate ( sigma_jon(0:Npart_types*Nrb_sites_max-1) )
   If(sigma_jon_aux.gt.0.d0) then

      Do I = 0, Npart_types-1
         Indx = I*Nrb_sites_max
         Do J = 0, Nrb_sites(I)-1
             sigma_jon(Indx+J) = sigma_jon_aux
         EndDo
      EndDo

   Else

      Read(20,*) 
      Do I= 0, Npart_types-1
         Indx = I*Nrb_sites_max
         print*, 'i',i,Nrb_sites(i)
         Do J = 0, Nrb_sites(i)-1
            Read(20,'(F6.3)',advance='no') sigma_jon_ij
            sigma_jon(indx+j) = sigma_jon_ij 
            print*, 'sigma_LJ i,j',i,j,indx+j,sigma_jon(indx+j)
         EndDo
         Read(20,*) 
      EndDo

      Allocate ( sigma_LJ(0:Npart_types*Npart_types-1) )
      Allocate ( rangep(0:Npart_types*Npart_types-1) )
      Allocate ( VL0(0:Npart_types*Npart_types-1) )
      Allocate ( xop(0:Npart_types*Npart_types-1) )

      Read(20,*) 
      Do I=0,Npart_types-1

           Indx = I*Npart_types
           Do J = 0, Npart_types-1
              Read(20,'(F6.3)',advance='no') sigma_jon_ij
              sigma_LJ(indx+j) = sigma_jon_ij 
              rangep(indx+j) = rangepp * sigma_jon_ij
              VL0(indx+j)=(sigma_jon_ij/rangep(indx+j))**12.d0-(sigma_jon_ij/rangep(indx+j))**6.0
              xopp = (1.d0+sqrt(1.d0+4.d00*VL0(indx+j)))/2.d0
              xop(indx+j) = sigma_jon_ij/(xopp**(1.d0/6.d0))
              print*, 'sigma_LJ',sigma_LJ(indx+j),rangep(indx+j),VL0(indx+j)
           EndDo
           Read(20,*) 
      EndDo

   EndIf

   close(20)

End Subroutine Read_input_data

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Read_init_conf

   Use set_precision
   Use configuration
   Use rundata, Only : base
   Implicit None
   Integer :: i, Indx, J
   Character(len=3) :: el

   print*, 'reading initial configuration'
   Open(21,file='coords-'//base)

   read(21,*) natoms
   Do I=1,ndim
      Indx=(I-1)*Ndim
      read(21,*) (h(Indx+J),J=0,Ndim-1)
   End Do

   Allocate(r(0:natoms*ndim-1))
   Allocate(q(0:natoms*4-1))

   Allocate(itype(0:natoms-1))

   print*, 'natoms',natoms

   Do I=0,natoms-1
        Indx=I*Ndim
        read(21,*) el,(r(Indx+J),J=0,Ndim-1)
        if(el.eq.'C1') then  ! change this when introducing LA
            itype(i)=0
        elseif(el.eq.'C2') then
            itype(i)=1
        elseif(el.eq.'C3') then
            itype(i)=2
        elseif(el.eq.'C4') then
            itype(i)=3
        elseif(el.eq.'C5') then
            itype(i)=4
        elseif(el.eq.'C6') then
            itype(i)=5
        elseif(el.eq.'C7') then
            itype(i)=6
        elseif(el.eq.'C8') then
            itype(i)=7
        elseif(el.eq.'C9') then
            itype(i)=8
        elseif(el.eq.'C10') then
            itype(i)=9  
        elseif(el.eq.'C11') then
            itype(i)=10 
        elseif(el.eq.'C12') then
            itype(i)=11 
        elseif(el.eq.'C13') then
            itype(i)=12 
        elseif(el.eq.'C14') then
            itype(i)=13 
        elseif(el.eq.'C15') then
            itype(i)=14 
        elseif(el.eq.'C16') then
            itype(i)=15 
        elseif(el.eq.'C17') then
            itype(i)=16 
        elseif(el.eq.'C18') then
            itype(i)=17 
        elseif(el.eq.'C19') then
            itype(i)=18
        elseif(el.eq.'C20') then
            itype(i)=19 
        elseif(el.eq.'C21') then
            itype(i)=20 
        elseif(el.eq.'C22') then
            itype(i)=21 
        elseif(el.eq.'C23') then
            itype(i)=22 
        elseif(el.eq.'C24') then
            itype(i)=23 
        elseif(el.eq.'C25') then
            itype(i)=24 
        elseif(el.eq.'C26') then
            itype(i)=25 
        elseif(el.eq.'C27') then
            itype(i)=26 
        elseif(el.eq.'C28') then
            itype(i)=27 
        elseif(el.eq.'C29') then
            itype(i)=28
        elseif(el.eq.'C30') then
            itype(i)=29 
        elseif(el.eq.'C31') then
            itype(i)=30 
        elseif(el.eq.'C32') then
            itype(i)=31 
        elseif(el.eq.'C33') then
            itype(i)=32 
        elseif(el.eq.'C34') then
            itype(i)=33 
        elseif(el.eq.'C35') then
            itype(i)=34 
        elseif(el.eq.'C36') then
            itype(i)=35 
        elseif(el.eq.'C37') then
            itype(i)=36 
        elseif(el.eq.'C38') then
            itype(i)=37 
        elseif(el.eq.'C39') then
            itype(i)=38
        elseif(el.eq.'C40') then
            itype(i)=39 
        else
            print*, 'this type of particle is not implemented yet'
            write(*,'(a3)') el
            stop
        endif

   End Do

   Do I=0,natoms-1
        Indx=I*4
        read(21,*) (q(Indx+J),J=0,3)
   End Do
   
   Close (21)
   print*, 'exit Read_init_conf'

End Subroutine Read_init_conf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Get_edges_angles

   Use set_precision
   Use configuration, Only : ndim, h, box, bangles
   Implicit None
   Integer   :: l, Indx

   Box=0.
   Bangles=0.

   do l=1,3

        Indx=(L-1)*Ndim
        Box(0) = Box(0)+h(Indx)**2
        Box(1) = Box(1)+h(Indx+1)**2
        Box(2) = Box(2)+h(Indx+2)**2

        Bangles(0) = Bangles(0)+h(Indx)*h(Indx+1)
        Bangles(1) = Bangles(1)+h(Indx)*h(Indx+2)
        Bangles(2) = Bangles(2)+h(Indx+1)*h(Indx+2)

   enddo

   Box(:) = sqrt(Box(:))
   print*, 'Box', Box(:)

   Bangles(0) = Bangles(0)/Box(0)/Box(1)
   Bangles(1) = Bangles(1)/Box(0)/Box(2)
   Bangles(2) = Bangles(2)/Box(1)/Box(2)
   print*, 'Angles', acos(Bangles(:))*180./acos(-1.0)

End Subroutine Get_edges_angles

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Write_final_conf2

   Use set_precision
   Use configuration
   Use potential        
   Use rundata, Only : base
   Implicit None
   Integer :: i, Indx, J, itypei, Indxi, Indxp, Num_tot
   Real(wp)   :: q0, q1, q2, q3,  R11, R12, R13, R21, R22, R23
   Real(wp)   :: R31, R32, R33, px_ref, py_ref, pz_ref, px, py, pz
   Character(len=3) :: el
   Character(len=5) :: el2

   write(674,*) natoms
   Do I=1,ndim
      Indx=(I-1)*Ndim
      write(674,*) (h(Indx+J),J=0,Ndim-1)
      !write(*,*) (h(Indx+J),J=0,Ndim-1)
   End Do

   Num_tot=0
   Do I=0,natoms-1
        Indx=I*Ndim
        if(itype(i).eq.0) then  ! change this when introducing LA
            el='C1'
        elseif(itype(i).eq.1) then  
            el='C2'
        elseif(itype(i).eq.2) then  
            el='C3'
        elseif(itype(i).eq.3) then  
            el='C4'
        elseif(itype(i).eq.4) then  
            el='C5'
        elseif(itype(i).eq.5) then  
            el='C6'
        elseif(itype(i).eq.6) then  
            el='C7'
        elseif(itype(i).eq.7) then  
            el='C8'
        elseif(itype(i).eq.8) then  
            el='C9'
        elseif(itype(i).eq.9) then  
            el='C10'
        elseif(itype(i).eq.10) then  
            el='C11'
        elseif(itype(i).eq.11) then  
            el='C12'
        elseif(itype(i).eq.12) then  
            el='C13'
        elseif(itype(i).eq.13) then  
            el='C14'
        elseif(itype(i).eq.14) then  
            el='C15'
        elseif(itype(i).eq.15) then  
            el='C16'
        elseif(itype(i).eq.16) then  
            el='C17'
        elseif(itype(i).eq.17) then  
            el='C18'
        elseif(itype(i).eq.18) then  
            el='C19'
        elseif(itype(i).eq.19) then  
            el='C20'
        elseif(itype(i).eq.20) then  
            el='C21'
        elseif(itype(i).eq.21) then  
            el='C22'
        elseif(itype(i).eq.22) then  
            el='C23'
        elseif(itype(i).eq.23) then  
            el='C24'
        elseif(itype(i).eq.24) then  
            el='C25'
        elseif(itype(i).eq.25) then  
            el='C26'
        elseif(itype(i).eq.26) then  
            el='C27'
        elseif(itype(i).eq.27) then  
            el='C28'
        elseif(itype(i).eq.28) then  
            el='C29'
        elseif(itype(i).eq.29) then  
            el='C30'
        elseif(itype(i).eq.30) then  
            el='C31'
        elseif(itype(i).eq.31) then  
            el='C32'
        elseif(itype(i).eq.32) then  
            el='C33'
        elseif(itype(i).eq.33) then  
            el='C34'
        elseif(itype(i).eq.34) then  
            el='C35'
        elseif(itype(i).eq.35) then  
            el='C36'
        elseif(itype(i).eq.36) then  
            el='C37'
        elseif(itype(i).eq.37) then  
            el='C38'
        elseif(itype(i).eq.38) then  
            el='C39'
        elseif(itype(i).eq.39) then  
            el='C40'
        else
            print*, 'this type of particle is not implemented yet'
            stop
        endif
        write(674,'(A3,5x,3F20.8)') el,(r(Indx+J),J=0,Ndim-1)
        itypei= itype(i)
        Num_tot = Num_tot +1
        Do J=0,Nrb_sites(itypei)-1
            Num_tot = Num_tot +1
        EndDo
   EndDo
   Do I=0,natoms-1
        Indx=I*4
        write(674,'(4F20.8)') (q(Indx+J),J=0,3)
   End Do
   
End Subroutine Write_final_conf2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Write_final_conf

   Use set_precision
   Use configuration
   Use potential        
   Use rundata, Only : base
   Implicit None
   Integer :: i, Indx, J, itypei, Indxi, Indxp, Num_tot
   Real(wp)   :: q0, q1, q2, q3,  R11, R12, R13, R21, R22, R23
   Real(wp)   :: R31, R32, R33, px_ref, py_ref, pz_ref, px, py, pz
   Character(len=3) :: el
   Character(len=5) :: el2

   Open(21,file='coords.out-'//base)
   write(21,*) natoms
   Do I=1,ndim
      Indx=(I-1)*Ndim
      write(21,*) (h(Indx+J),J=0,Ndim-1)
      !write(*,*) (h(Indx+J),J=0,Ndim-1)
   End Do

   Num_tot=0
   Do I=0,natoms-1
        Indx=I*Ndim
        if(itype(i).eq.0) then  ! change this when introducing LA
            el='C1'
        elseif(itype(i).eq.1) then  
            el='C2'
        elseif(itype(i).eq.2) then  
            el='C3'
        elseif(itype(i).eq.3) then  
            el='C4'
        elseif(itype(i).eq.4) then  
            el='C5'
        elseif(itype(i).eq.5) then  
            el='C6'
        elseif(itype(i).eq.6) then  
            el='C7'
        elseif(itype(i).eq.7) then  
            el='C8'
        elseif(itype(i).eq.8) then  
            el='C9'
        elseif(itype(i).eq.9) then  
            el='C10'
        elseif(itype(i).eq.10) then  
            el='C11'
        elseif(itype(i).eq.11) then  
            el='C12'
        elseif(itype(i).eq.12) then  
            el='C13'
        elseif(itype(i).eq.13) then  
            el='C14'
        elseif(itype(i).eq.14) then  
            el='C15'
        elseif(itype(i).eq.15) then  
            el='C16'
        elseif(itype(i).eq.16) then  
            el='C17'
        elseif(itype(i).eq.17) then  
            el='C18'
        elseif(itype(i).eq.18) then  
            el='C19'
        elseif(itype(i).eq.19) then  
            el='C20'
        elseif(itype(i).eq.20) then  
            el='C21'
        elseif(itype(i).eq.21) then  
            el='C22'
        elseif(itype(i).eq.22) then  
            el='C23'
        elseif(itype(i).eq.23) then  
            el='C24'
        elseif(itype(i).eq.24) then  
            el='C25'
        elseif(itype(i).eq.25) then  
            el='C26'
        elseif(itype(i).eq.26) then  
            el='C27'
        elseif(itype(i).eq.27) then  
            el='C28'
        elseif(itype(i).eq.28) then  
            el='C29'
        elseif(itype(i).eq.29) then  
            el='C30'
        elseif(itype(i).eq.30) then  
            el='C31'
        elseif(itype(i).eq.31) then  
            el='C32'
        elseif(itype(i).eq.32) then  
            el='C33'
        elseif(itype(i).eq.33) then  
            el='C34'
        elseif(itype(i).eq.34) then  
            el='C35'
        elseif(itype(i).eq.35) then  
            el='C36'
        elseif(itype(i).eq.36) then  
            el='C37'
        elseif(itype(i).eq.37) then  
            el='C38'
        elseif(itype(i).eq.38) then  
            el='C39'
        elseif(itype(i).eq.39) then  
            el='C40'
        else
            print*, 'this type of particle is not implemented yet'
            stop
        endif
        write(21,'(A3,5x,3F20.8)') el,(r(Indx+J),J=0,Ndim-1)
        itypei= itype(i)
        Num_tot = Num_tot +1
        Do J=0,Nrb_sites(itypei)-1
            Num_tot = Num_tot +1
        EndDo
   EndDo
   Do I=0,natoms-1
        Indx=I*4
        write(21,'(4F20.8)') (q(Indx+J),J=0,3)
        !write(*,'(4F20.8)') i,(q(Indx+J),J=0,3)
   End Do
   
   Close (21)

End Subroutine Write_final_conf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine System_energy(overlap)

   Use set_precision
   Use configuration
   Use properties, Only : en_tot
   Use rundata, Only: istep_fin
   Implicit None
   Integer :: i, j
   Real(wp) :: en, ener_part(0:natoms-1), en_tot1
   Logical :: overlap

   overlap=.false.
   en_tot=0.e0
   ener_part=0
   Do I=0,Natoms-2
      Do J=i+1,Natoms-1
          call ener(I,J,en,overlap)
          ener_part(i)=ener_part(i)+en
          ener_part(j)=ener_part(j)+en
          if(overlap) Then
               en_tot=9.e+33
               Return
          EndIf
          en_tot= en_tot + en
      End Do
   End Do
   en_tot1=0.d0
   Do i=0,Natoms-1
       !write(61,*) i,ener_part(i),itype(i)
       en_tot1=en_tot1+ener_part(i)
   EndDo
   en_tot1=en_tot1/dble(Natoms)/2.d0
   print*, 'En_tot1 = ',En_tot1
   If(istep_fin.eq.0) then
        print*, ''
        print*, ''
        call print_ener
        print*, ''
        print*, ''
        print*, 'Printing configuration in pdb format and exiting'
        print*, ''
        print*, ''
        Call Convert_Real_Units
        Call Write_pdb(ener_part)
        stop
   Endif

End Subroutine System_energy

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Write_pdb(ener_part)        

   Use set_precision
   Use configuration
   Use potential        
   Use rundata, Only : base
   Implicit None
   Integer :: i, Indx, J, itypei, Indxi, Indxp, Num_tot
   Real(wp)   :: q0, q1, q2, q3,  R11, R12, R13, R21, R22, R23
   Real(wp)   :: R31, R32, R33, px_ref, py_ref, pz_ref, px, py, pz
   Real(wp)   :: angle,alpha
   Real(wp)   :: ener_part(0:natoms-1)
   Character(len=3) :: el
   Character(len=5) :: el2

   Open(21,file='coords.pdb')
   write(21,*) "MODEL    "
   angle=90.0
   write(21,'("CRYST1",6f9.3)') h(0),h(4),h(8),angle,angle,angle

   alpha=1.0
   Do I=0,natoms-1
        Indx=I*Ndim
        if(itype(i).eq.0) then  ! change this when introducing LA
            el='C1'
        elseif(itype(i).eq.1) then  
            el='C2'
        elseif(itype(i).eq.2) then  
            el='C3'
        elseif(itype(i).eq.3) then  
            el='C4'
        elseif(itype(i).eq.4) then  
            el='C5'
        elseif(itype(i).eq.5) then  
            el='C6'
        elseif(itype(i).eq.6) then  
            el='C7'
        elseif(itype(i).eq.7) then  
            el='C8'
        elseif(itype(i).eq.8) then  
            el='C9'
        elseif(itype(i).eq.9) then  
            el='C10'
        elseif(itype(i).eq.10) then  
            el='C11'
        elseif(itype(i).eq.11) then  
            el='C12'
        elseif(itype(i).eq.12) then  
            el='C13'
        elseif(itype(i).eq.13) then  
            el='C14'
        elseif(itype(i).eq.14) then  
            el='C15'
        elseif(itype(i).eq.15) then  
            el='C16'
        elseif(itype(i).eq.16) then  
            el='C17'
        elseif(itype(i).eq.17) then  
            el='C18'
        elseif(itype(i).eq.18) then  
            el='C19'
        elseif(itype(i).eq.19) then  
            el='C20'
        elseif(itype(i).eq.20) then  
            el='C21'
        elseif(itype(i).eq.21) then  
            el='C22'
        elseif(itype(i).eq.22) then  
            el='C23'
        elseif(itype(i).eq.23) then  
            el='C24'
        elseif(itype(i).eq.24) then  
            el='C25'
        elseif(itype(i).eq.25) then  
            el='C26'
        elseif(itype(i).eq.26) then  
            el='C27'
        elseif(itype(i).eq.27) then  
            el='C28'
        elseif(itype(i).eq.28) then  
            el='C29'
        elseif(itype(i).eq.29) then  
            el='C30'
        elseif(itype(i).eq.30) then  
            el='C31'
        elseif(itype(i).eq.31) then  
            el='C32'
        elseif(itype(i).eq.32) then  
            el='C33'
        elseif(itype(i).eq.33) then  
            el='C34'
        elseif(itype(i).eq.34) then  
            el='C35'
        elseif(itype(i).eq.35) then  
            el='C36'
        elseif(itype(i).eq.36) then  
            el='C37'
        elseif(itype(i).eq.37) then  
            el='C38'
        elseif(itype(i).eq.38) then  
            el='C39'
        elseif(itype(i).eq.39) then  
            el='C40'
        else
            print*, 'this type of particle is not implemented yet'
            stop
        endif
        write(21,'("ATOM",i7,x,a2,3x,"MOL",10x,3f8.3,2f8.4)') i+1,el,(r(Indx+J),J=0,Ndim-1),ener_part(i),ener_part(i)
        !itypei= itype(i)
       ! Num_tot = Num_tot +1
        !Do J=0,Nrb_sites(itypei)-1
        !    Num_tot = Num_tot +1
        !EndDo
   EndDo
   write(21,*) "ENDMDL    "
   
   Close (21)

End Subroutine Write_pdb        

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine print_ener

   Use set_precision
   Use properties, Only : en_tot
   Use configuration, Only : natoms
   Implicit None

   Print*, ''
   Print*, ''
   Print*, 'Total energy of the initial configuration',en_tot
   Print*, 'Total energy per particle',en_tot/natoms
   Print*, ''
   Print*, ''

End Subroutine Print_ener

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Initialize_variables

   Use set_precision
   Use properties
   Use configuration, Only : natoms
   Use rundata
   Implicit None

   ntrans = 0 
   nrot = 0 
   ntrans_ac = 0 
   nrot_ac = 0 
   nvol = 0
   nvol_accept = 0

   rho_av=0.0
   en_av=0.0
   lx_av=0.0
   ly_av=0.0
   lz_av=0.0
   open(40,file='run-data-'//base,access='append')

End Subroutine Initialize_variables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Accumulate_averages

   Use set_precision
   Use properties 
   Use rundata, Only: temp
   Use configuration
   Implicit None

   if(ndim.eq.2) then
     vol=h(0)*h(3)-h(1)*h(2)
     print*, 'so far the code is just ready for 3D'
     stop
   else
     vol=h(0)*(h(4)*h(8)-h(5)*h(7))-h(3)*(h(1)*h(8) &
        & -h(2)*h(7))+h(6)*(h(1)*h(5)-h(2)*h(4))

   endif
   rho=natoms/vol
   rho_av=rho_av+natoms/vol
   en_av=en_av+en_tot
   lx_av=lx_av+box(0)
   ly_av=ly_av+box(1)
   lz_av=lz_av+box(2)
   write(40,'(7f18.5)') temp,rho,vol,box,en_tot/natoms

End Subroutine Accumulate_averages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Calculate_averages

   Use set_precision
   Use properties
   Use rundata
   Use configuration, Only : natoms
   Implicit None
   Integer   :: Npoints
   Real(wp)  :: prob_trans, prob_rot

   Npoints=(istep_fin-istep_ini-Neq)/Nsave

   rho_av=rho_av/Npoints
   en_av=en_av/Npoints
   lx_av=lx_av/Npoints
   ly_av=ly_av/Npoints
   lz_av=lz_av/Npoints

   print*, 'Average density = ',rho_av
   print*, 'Average energy = ',en_av/natoms,en_av
   print*, 'Average edge lx= ',lx_av
   print*, 'Average edge ly= ',ly_av
   print*, 'Average edge lz= ',lz_av

End Subroutine Calculate_averages

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Write_final_properties

   Use set_precision
   Use properties, Only : en_tot
   Use configuration, Only : natoms 
   Implicit None
   Logical :: overlap

   print*, ''
   print*, ''
   print*, 'Final energy',en_tot,en_tot/dble(natoms)
   Call system_energy(overlap)
   If (overlap) Then
      print*, 'Overlap in final configuration'
   EndIf
   print*, 'New calculation of final energy',en_tot,en_tot/natoms
   print*, ''
   print*, ''

End Subroutine Write_final_properties

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!Subroutine Random_orientation(qN)
!
!   Use set_precision
!   Real(wp), Dimension(4) :: qN
!   Real(wp), Dimension(2) :: dice
!   Real(wp)               :: ransq12, ransq34, ranh, xran1, xran2, xran3, xran4
!
!   ransq12=2.d00
!   do while (ransq12 >= 1.d00)
!        Call Random_number(dice)
!        xran1=1.d00-2.d00*dice(1)
!        xran2=1.d00-2.d00*dice(2)
!        ransq12=xran1*xran1+xran2*xran2
!   enddo
!   ransq34=2.d00
!   do while (ransq34 >= 1.d00)
!        Call Random_number(dice)
!        xran3=1.d00-2.d00*dice(1)
!        xran4=1.d00-2.d00*dice(2)
!        ransq34=xran3*xran3+xran4*xran4
!   enddo
!   ranh=sqrt((1.d00-ransq12)/ransq34)
!
!   qN(1)=xran1
!   qN(2)=xran2
!   qN(3)=xran3*ranh
!   qN(4)=xran4*ranh
!
!End Subroutine Random_orientation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Get_Inverse_matrix

   Use set_precision
   Use configuration
   Real(wp) :: detH

   If (ndim == 2) Then

       print*, 'Two dimensional case not implemented yet'
       stop

   Elseif (ndim == 3) Then  

       detH=h(0)*(h(4)*h(8)-h(5)*h(7))-h(3)*(h(1)*h(8) &
     &  -h(2)*h(7))+h(6)*(h(1)*h(5)-h(2)*h(4))

       hinv(0)=(h(4)*h(8)-h(5)*h(7))/detH
       hinv(3)=-(h(3)*h(8)-h(6)*h(5))/detH
       hinv(6)=(h(3)*h(7)-h(6)*h(4))/detH
       hinv(1)=-(h(1)*h(8)-h(7)*h(2))/detH
       hinv(4)=(h(0)*h(8)-h(6)*h(2))/detH
       hinv(7)=-(h(0)*h(7)-h(6)*h(1))/detH
       hinv(2)=(h(1)*h(5)-h(4)*h(2))/detH
       hinv(5)=-(h(0)*h(5)-h(3)*h(2))/detH
       hinv(8)=(h(0)*h(4)-h(3)*h(1))/detH
       print*, 'h'
       print*, h(0),h(1),h(2)
       print*, h(3),h(4),h(5)
       print*, h(6),h(7),h(8)
       print*, 'hinv'
       print*, hinv(0),hinv(1),hinv(2)
       print*, hinv(3),hinv(4),hinv(5)
       print*, hinv(6),hinv(7),hinv(8)

   Else

       print*, 'This dimension is not implemented yet', ndim
       stop

   Endif 

End Subroutine Get_Inverse_matrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Convert_Box_Units

   Use set_precision
   Use configuration
   Integer :: i,Indx
   Real(wp), Dimension(ndim) :: x1, x1p  

   !open(45,file='coords-reduced.xyz')
   !write(45,*) Natoms
   !write(45,*) 
   Do i=0,Natoms-1
       Indx=I*Ndim
       x1(1) = r(Indx)
       x1(2) = r(Indx+1)
       x1(3) = r(Indx+2)
       x1p(1) = hinv(0)*x1(1)+hinv(1)*x1(2)+hinv(2)*x1(3)
       x1p(2) = hinv(3)*x1(1)+hinv(4)*x1(2)+hinv(5)*x1(3)
       x1p(3) = hinv(6)*x1(1)+hinv(7)*x1(2)+hinv(8)*x1(3)
       r(Indx) = x1p(1)
       r(Indx+1) = x1p(2)
       r(Indx+2) = x1p(3)
       !write(45,'("C",3x,3f16.6)') x1p(1),x1p(2),x1p(3)
   EndDo

End Subroutine Convert_Box_Units

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Convert_Real_Units

   Use set_precision
   Use configuration
   Integer :: i, Indx
   Real(wp), Dimension(ndim) :: x1, x1p  

   Do i=0,Natoms-1
       Indx = I*Ndim
       x1(1) = r(Indx)
       x1(2) = r(Indx+1)
       x1(3) = r(Indx+2)
       x1(:) = x1(:) - int( x1(:) )
       x1p(1) = h(0)*x1(1)+h(1)*x1(2)+h(2)*x1(3)
       x1p(2) = h(3)*x1(1)+h(4)*x1(2)+h(5)*x1(3)
       x1p(3) = h(6)*x1(1)+h(7)*x1(2)+h(8)*x1(3)
       r(Indx) = x1p(1)
       r(Indx+1) = x1p(2)
       r(Indx+2) = x1p(3)
   EndDo

End Subroutine Convert_Real_Units


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subroutine Write_movie(istep)

   Use set_precision
   Use configuration
   Use potential
   Use rundata, Only : base, beta
   Implicit None
   Integer :: i, Indx, J, itypei, Indxi, Indxp, Num_tot, istep
   Real(wp)   :: q0, q1, q2, q3,  R11, R12, R13, R21, R22, R23
   Real(wp)   :: R31, R32, R33, px_ref, py_ref, pz_ref, px, py, pz
   Character(len=3) :: el
   Character(len=5) :: el2

   Num_tot=0
   Do I=0,natoms-1
        Indx=I*Ndim
        itypei= itype(i)
        Num_tot = Num_tot +1
        Do J=0,Nrb_sites(itypei)-1
            Num_tot = Num_tot +1
        EndDo
   EndDo

   !print*, 'Num_tot',Num_tot

   write(133,*) Num_tot
   write(133,*) h(0),istep,1.0/beta

   Do I=0,natoms-1
        Indx=I*Ndim
        if(itype(i).eq.0) then  ! change this when introducing LA
            el='C1'
        elseif(itype(i).eq.1) then
            el='C2'
        elseif(itype(i).eq.2) then
            el='C3'
        elseif(itype(i).eq.3) then
            el='C4'
        elseif(itype(i).eq.4) then
            el='C5'
        elseif(itype(i).eq.5) then
            el='C6'
        elseif(itype(i).eq.6) then
            el='C7'
        elseif(itype(i).eq.7) then
            el='C8'
        elseif(itype(i).eq.8) then
            el='C9'
        elseif(itype(i).eq.9) then
            el='C10'
        elseif(itype(i).eq.10) then
            el='C11'
        elseif(itype(i).eq.11) then
            el='C12'
        elseif(itype(i).eq.12) then
            el='C13'
        elseif(itype(i).eq.13) then
            el='C14'
        elseif(itype(i).eq.14) then
            el='C15'
        elseif(itype(i).eq.15) then
            el='C16'
        elseif(itype(i).eq.16) then
            el='C17'
        elseif(itype(i).eq.17) then
            el='C18'
        elseif(itype(i).eq.18) then
            el='C19'
        elseif(itype(i).eq.19) then
            el='C20'
        elseif(itype(i).eq.20) then
            el='C21'
        elseif(itype(i).eq.21) then
            el='C22'
        elseif(itype(i).eq.22) then
            el='C23'
        elseif(itype(i).eq.23) then
            el='C24'
        elseif(itype(i).eq.24) then
            el='C25'
        elseif(itype(i).eq.25) then
            el='C26'
        elseif(itype(i).eq.26) then
            el='C27'
        elseif(itype(i).eq.27) then
            el='C28'
        elseif(itype(i).eq.28) then
            el='C29'
        elseif(itype(i).eq.29) then
            el='C30'
        elseif(itype(i).eq.30) then
            el='C31'
        elseif(itype(i).eq.31) then
            el='C32'
        elseif(itype(i).eq.32) then
            el='C33'
        elseif(itype(i).eq.33) then
            el='C34'
        elseif(itype(i).eq.34) then
            el='C35'
        elseif(itype(i).eq.35) then
            el='C36'
        elseif(itype(i).eq.36) then
            el='C37'
        elseif(itype(i).eq.37) then
            el='C38'
        elseif(itype(i).eq.38) then
            el='C39'
        elseif(itype(i).eq.39) then
            el='C40'
        else
            print*, 'this type of particle is not implemented yet'
            stop
        endif
        write(133,'(A3,5x,3F20.8,5x,i6)') el,(r(Indx+J),J=0,Ndim-1)

        itypei= itype(i)

        Indxi=I*4

        q0=q(Indxi)
        q1=q(Indxi+1)
        q2=q(Indxi+2)
        q3=q(Indxi+3)
        R11=q0*q0+q1*q1-q2*q2-q3*q3
        R12=2.0*(q1*q2-q0*q3)
        R13=2.0*(q1*q3+q0*q2)
        R21=2.0*(q1*q2+q0*q3)
        R22=q0*q0-q1*q1+q2*q2-q3*q3
        R23=2.0*(q2*q3-q0*q1)
        R31=2.0*(q1*q3-q0*q2)
        R32=2.0*(q2*q3+q0*q1)
        R33=q0*q0-q1*q1-q2*q2+q3*q3

        Do J=0,Nrb_sites(itypei)-1

            Indxp = Nrb_sites_max*Ndim*itypei + J*Ndim

            px_ref= PATCH( Indxp )
            py_ref= PATCH( Indxp + 1 )
            pz_ref= PATCH( Indxp + 2 )
            px = r( Indx )    + 0.35 * ( R11*px_ref + R12*py_ref + R13*pz_ref )
            py = r( Indx +1 ) + 0.35 * ( R21*px_ref + R22*py_ref + R23*pz_ref )
            pz = r( Indx +2 ) + 0.35 * ( R31*px_ref + R32*py_ref + R33*pz_ref )

            if((itypei+1.lt.10).and.(j+1.lt.10)) then
                 write( el2, '(a1,I1,I1)') "P",itypei+1,J+1
             elseif((itypei+1.ge.10).and.(j+1.lt.10)) then
                 write( el2, '(a1,I2,I1)') "P",itypei+1,J+1
             elseif((itypei+1.lt.10).and.(j+1.ge.10)) then
                 write( el2, '(a1,I1,I2)') "P",itypei+1,J+1
             elseif((itypei+1.ge.10).and.(j+1.ge.10)) then
                 write( el2, '(a1,I2,I2)') "P",itypei+1,J+1
             endif
             el2= trim(el2)
             !print*, 'el2=',el2
             write(133,'(a5,4x,3F20.8)') el2,px,py,pz
        EndDo
   End Do

End Subroutine Write_movie


                                                                                         
