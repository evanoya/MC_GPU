
Program GPU_MC

 Use set_precision

 Use configuration
 Use properties

 Use rundata
 Use cell       

 Implicit None
 Integer :: i, j, ll, Iend, Ieq, vec_shift(6,0:ndim-1), vec(0:ndim-1)

 Integer :: ngpus,gpuid,Njumps, Nsteps_DeltaT, Natoms_aux

 Real(wp) :: rnd, disp(0:ndim-1), start, finish, omaxo, hmaxo, prob_vol, vmaxo(0:ndim-1)
 Real(wp) :: rnd_vec(ndim), prob_trans, prob_rot, En_tot_old, boxmin, W_RU_min
 Logical :: overlap

 ! Check if there are GPU's available

 call gpucount(ngpus)
 write(*,*) '*********************************************************'
 write(*,*) '*********************************************************'
 write(*,*) 
 write(*,*) 
 write(*,"('*** No. of GPU devices available:',i2)") ngpus
 write(*,*) 
 write(*,*) 
 write(*,*) '*********************************************************'
 write(*,*) '*********************************************************'

 If ( ngpus == 0) Then
     print*, 'ERROR: There is not any GPU device available'
     stop
 ElseIf ( ngpus >= 2) Then
     read(*,*) gpuid
     call choosegpu(gpuid)
 Endif
 ! Read the input file

 Call Read_input_data
  
 print*, 'input file read'
 print*, ''

 Call Read_init_conf
 
 print*, ''
 print*, 'initial configuration read'

 Call Get_Inverse_matrix
 Call Convert_Box_Units
 
 print*, ''
 print*, 'Change to Box Units'

 Call Initialize_variables

 print*, ''
 print*, 'variables initialized'

 call cpu_time(start)
 Call system_energy(overlap) 
                             
 call cpu_time(finish)
 print*, ''
 print*, ''
 write(*,'("Time initial energy=",e13.6,"   seconds.")') finish-start
 print*, ''
 print*, ''
 If (overlap) Then

    print*, 'Overlap in initial configuration'
    stop

 EndIf

 Call print_ener


 Call Initialize_checkerboard

 NCell_max=Ntot_Cell*10

 print*, ''
 print*, 'checkerboard initialized'
! 
 call cpu_time(start)

 En_tot =0
 print*, 'about to enter CB_system_energy'
 Call CB_system_energy_gpu(init)
 print*, 'Initial energy evaluated with GPU',En_tot
 call cpu_time(finish)
 print*, ''
 print*, ''
 write(*,'("Time initial energy (checkerboard) =",e13.6,"   seconds.")') finish-start
 print*, ''
 print*, ''

 If (overlap) Then

    print*, 'Overlap in initial configuration'
    stop

 EndIf
 Iend=istep_ini+istep_fin
 IEq=istep_ini+Neq
 vec_shift(1,:) = (/  1,  0,  0 /)
 vec_shift(2,:) = (/ -1,  0,  0 /)
 vec_shift(3,:) = (/  0,  1,  0 /)
 vec_shift(4,:) = (/  0, -1,  0 /)
 vec_shift(5,:) = (/  0,  0,  1 /)
 vec_shift(6,:) = (/  0,  0, -1 /)

 init=1

 print*, ''
 print*, 'Starting the main simulation loop'
 print*, ''

 open(133,file='movie.xyz',access='append')
 open(674,file='movie-quat.dat',access='append')

 no_mc_real=0
 ntrans_ac=0
 nrot_ac=0

 Njumps=(istep_fin-istep_ini)/Nsave
 DeltaT=((Temp1-Temp0)/dble(Njumps-1))
 Temp = Temp0
 beta=1.d0/Temp
 print*, ''
 print*, 'Number of temperature jumps',Njumps
 print*, 'i-istep_ini',istep_fin-istep_ini
 print*, 'DeltaT',DeltaT
 print*, ''
 print*, ''
 print*, 'Step,En_tot, prob_trans,prob_rot,ntrans,hmax*box(0),omax,prob_vol,vmax,No_cell'
 print*, ''

 call cpu_time(start)
 Do I=1+istep_ini,istep_ini+istep_fin

      Call MCsweep()

      Call Random_number(rnd)
      disp(:) = rnd*disp_max(:)
      Call Random_number(rnd)
      vec(:)  = vec_shift ( int(rnd*size(vec_shift,1))+1, :)
!
      Call Shift_Cells(vec,disp)
      If(NpT) then
             Call Volume_move
             nvol=nvol+1
      EndIf

      If( (i > Ieq) .AND. (mod(i,Nsave)==0) ) Then
           Call Accumulate_averages
           Ntrans=no_mc_real*Nmove/2  
           Nrot=no_mc_real*Nmove/2
           prob_trans=dble(ntrans_ac)/dble(ntrans)
           prob_rot=dble(nrot_ac)/dble(nrot)
           prob_vol=dble(nvol_accept)/dble(nvol)

           If ( displ_update) then

               hmaxo=hmax
               hmax=hmax*dabs(prob_trans/0.4d0)
               if(hmax/hmaxo > 1.5d0) hmax=hmaxo*1.5
               if(hmax/hmaxo < 0.5d0) hmax=hmaxo*0.5
               W_RU_min=minval(W_RU)
               if(hmax > W_RU_min/2.d0 ) hmax =W_RU_min/2.d0
               omaxo=omax
               omax=omax*dabs(prob_rot/0.4d00)
               if(omax/omaxo.gt.1.5d00) omax=omaxo*1.5d00
               if(omax/omaxo.lt.0.5d00) omax=omaxo*0.5d00
               if(omax.gt.0.5d0) omax=0.5d0

               vmaxo(0:2)=vmax(0:2)
               boxmin=minval(box)
               if(prob_vol.eq.0) then
                    vmax(0:2)=vmaxo(0:2)*0.5
               else
                    if(vmax(0).ne.0) then
                      vmax(0)=vmax(0)*dabs(prob_vol/0.4d0)
                      vmax(1)=vmax(0)*box(1)/boxmin
                      vmax(2)=vmax(0)*box(2)/boxmin
                    else
                      vmax(1)=vmax(1)*dabs(prob_vol/0.4d0)
                      vmax(2)=vmax(2)*dabs(prob_vol/0.4d0)
                    endif
                    if(vmax(0)/vmaxo(0) > 1.5d0) vmax(0)=vmaxo(0)*1.5
                    if(vmax(1)/vmaxo(1) > 1.5d0) vmax(1)=vmaxo(1)*1.5
                    if(vmax(2)/vmaxo(2) > 1.5d0) vmax(2)=vmaxo(2)*1.5
                    if(vmax(0)/vmaxo(0) < 0.5d0) vmax(0)=vmaxo(0)*0.5
                    if(vmax(1)/vmaxo(1) < 0.5d0) vmax(1)=vmaxo(1)*0.5
                    if(vmax(2)/vmaxo(2) < 0.5d0) vmax(2)=vmaxo(2)*0.5
               endif

           EndIf

           write(*,'(i10,x,"En_tot",f18.8,x,2f12.8,x,i12,x,6f12.8,x,i10)') i,En_tot, &
        &   prob_trans,prob_rot,ntrans,hmax*box(0),omax,prob_vol,vmax,Ntot_Cell
           no_mc_real=0
           ntrans_ac=0
           nrot_ac=0
           nvol=0
           nvol_accept=0

           Temp=Temp+DeltaT
           beta=1.d0/Temp

      EndIf
      If( (i > Ieq) .AND. (mod(i,Nsave2)==0) ) Then
           Call Convert_CB_RU
           Call Convert_Real_Units
           Call Write_movie(i)
           Call Write_final_conf2
      EndIf
      If(mod(i,Nrestart)==0) Then
           Call Write_final_conf
           Call Write_Input_File(i)
           En_tot_old=En_tot
           Call CB_system_energy_gpu(1)
           print*, ''
           print*, 'Check point: En_tot=',En_tot_old,En_tot
           print*, ''
      EndIf

 End Do
 print*, ''
 print*, 'Final energy after MC moves',En_tot
 print*, ''
 call cpu_time(finish)
 close(133)

 En_tot =0.d0
 print*, ''
 print*, ''
 write(*,'("Time spent in the main loop =",e13.6,"   seconds.")') finish-start
 print*, ''
 print*, ''
  
 print*, 'call to CB_system_energy_gpu'
 overlap=.false.
 Call CB_system_energy_gpu(2)
 print*, 'Final energy evaluated with GPU',En_tot

 print*, 'before convert_CB_ru'
 Call Convert_CB_RU
 print*, 'after convert_CB_ru'
 
 Call Calculate_averages

 Call Write_final_properties

 print*, 'before convert_ru'
 Call Convert_Real_Units
 print*, 'after convert_ru'

 Call Write_final_conf
 Call limpio
 
 print*, 'cleaning finished'
 !print*, 'Final energy evaluated with CB',En_tot
 !print*, 'after write final conf'
 stop
End Program GPU_MC
