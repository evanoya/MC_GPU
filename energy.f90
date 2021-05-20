Subroutine ener(i1,i2,en,overlap)
 
   Use set_precision
   Use configuration
   Use potential
   Implicit None
   Integer :: i1,i2, i, k, Indxi, Indxj
   Integer :: itypei, itypej, j3, j4
   Integer :: I1MIN, I2MIN, INDX, INDX2, INDX3, INDX4
   Integer :: INDXT1, INDXT2, INDXT3
   Real(wp) :: en
   Real(wp) :: X1, Y1, Z1, X2, Y2, Z2, DX, DY, DZ
   Real(wp) :: Q10, Q11, Q12, Q13, Q20, Q21, Q22, Q23
   Real(wp) :: DXP, DYP, DZP, RSQ, R12
   Real(wp) :: PX, PY, PZ, CSM1, CSM2, CSM1M, CSM2M, CSOMG1, CSOMG2, OMG1, OMG2 
   Real(wp) :: Q102, Q112, Q122, Q132, Q10Q11, Q10Q12, Q10Q13, Q11Q12, Q11Q13, Q12Q13
   Real(wp) :: R111, R112, R113, R121, R122, R123, R131, R132, R133
   Real(wp) :: Q202, Q212, Q222, Q232, Q20Q21, Q20Q22, Q20Q23, Q21Q22, Q21Q23, Q22Q23
   Real(wp) :: R211, R212, R213, R221, R222, R223, R231, R232, R233
   Real(wp) :: VTOR, factor_energy, PXP, PYP, PZP
   Real(wp) :: PX1_ref(0:Nrb_sites_max-1), PX2_ref(0:Nrb_sites_max-1)
   Real(wp) :: PY1_ref(0:Nrb_sites_max-1), PY2_ref(0:Nrb_sites_max-1)
   Real(wp) :: PZ1_ref(0:Nrb_sites_max-1), PZ2_ref(0:Nrb_sites_max-1)
   Real(wp) :: PX1(0:Nrb_sites_max-1), PX2(0:Nrb_sites_max-1)
   Real(wp) :: PY1(0:Nrb_sites_max-1), PY2(0:Nrb_sites_max-1)
   Real(wp) :: PZ1(0:Nrb_sites_max-1), PZ2(0:Nrb_sites_max-1)
   Real(wp) :: cos1(0:Nrb_sites_max-1), cos2(0:Nrb_sites_max-1)
   Real(wp) :: v1x, v1y, v1z, v2x, v2y, v2z, mod_v1, mod_v2, Utor, Utor_max
   Real(wp) :: v12x, v12y, v12z, a , R126, R1212, VRAD, VANG, Vangtor, FX, FXI
   Real(wp) :: v12b2, cangle, torsion, ctor_max, cos_tor, Vangtor_max 
   Real(wp) :: Vang_max, Vtor_max, dtor, pi, dospi, sigma_ij
   Real(wp) :: sigma_LJ_aux, xop_aux, vl0_aux, rangep_aux
   Logical :: overlap

   pi = acos(-1.d0)
   dospi = 2.d0*pi
   a = 1000.d0 

   Indxi=I1*Ndim
   X1=r(Indxi)
   Y1=r(Indxi+1)
   Z1=r(Indxi+2)

   Indxi=I1*4
   Q10=q(Indxi)
   Q11=q(Indxi+1)
   Q12=q(Indxi+2)
   Q13=q(Indxi+3)

   itypei = itype(I1)

   Indxj = I2*Ndim
   X2=r(Indxj)
   Y2=r(Indxj+1)
   Z2=r(Indxj+2)

   Indxj = I2*4
   Q20=q(Indxj)
   Q21=q(Indxj+1)
   Q22=q(Indxj+2)
   Q23=q(Indxj+3)

   itypej = itype(I2)

   DX=X2-X1
   DY=Y2-Y1
   DZ=Z2-Z1

   DX=DX-NINT(DX)
   DY=DY-NINT(DY)
   DZ=DZ-NINT(DZ)

   DXP=H(0)*DX+H(1)*DY+H(2)*DZ
   DYP=H(3)*DX+H(4)*DY+H(5)*DZ
   DZP=H(6)*DX+H(7)*DY+H(8)*DZ

   RSQ=DXP**2.0+DYP**2.0+DZP**2.
   R12=SQRT(RSQ)

   EN=0.0
   VTOR=0.d0
 
   Indx = itypei*Npart_types+itypej
   sigma_LJ_aux = sigma_LJ(Indx)
   xop_aux = xop(Indx)
   vl0_aux = vl0(Indx)
   rangep_aux = rangep(Indx)

   IF(R12.lt.0.5) then

        print*, 'dist 0',r12,i1,i2
        EN=9.e+33

   ELSEIF(R12.LT.0.9800000*sigma_LJ_aux) THEN

           R126=(RSQ/sigma_LJ_aux/sigma_LJ_aux)**3.d0
           R1212=R126*R126
           EN=-4.d0*((r1212 - r126)/(r1212*r126) + VL0_aux)

   ELSEIF(R12.LT.RANGEP_AUX) THEN

       ! rotation matrix for i1
        Q102=Q10*Q10
        Q112=Q11*Q11
        Q122=Q12*Q12
        Q132=Q13*Q13
        Q10Q11=Q10*Q11
        Q10Q12=Q10*Q12
        Q10Q13=Q10*Q13
        Q11Q12=Q11*Q12
        Q11Q13=Q11*Q13
        Q12Q13=Q12*Q13

        R111=Q102+Q112-Q122-Q132
        R112=2.d00*(Q11Q12-Q10Q13)
        R113=2.d00*(Q11Q13+Q10Q12)
        R121=2.d00*(Q11Q12+Q10Q13)
        R122=Q102-Q112+Q122-Q132
        R123=2.d00*(Q12Q13-Q10Q11)
        R131=2.d00*(Q11Q13-Q10Q12)
        R132=2.d00*(Q12Q13+Q10Q11)
        R133=Q102-Q112-Q122+Q132

       ! rotation matrix for i2

        Q202=Q20*Q20
        Q212=Q21*Q21
        Q222=Q22*Q22
        Q232=Q23*Q23
        Q20Q21=Q20*Q21
        Q20Q22=Q20*Q22
        Q20Q23=Q20*Q23
        Q21Q22=Q21*Q22
        Q21Q23=Q21*Q23
        Q22Q23=Q22*Q23

        R211=Q202+Q212-Q222-Q232
        R212=2.d00*(Q21Q22-Q20Q23)
        R213=2.d00*(Q21Q23+Q20Q22)
        R221=2.d00*(Q21Q22+Q20Q23)
        R222=Q202-Q212+Q222-Q232
        R223=2.d00*(Q22Q23-Q20Q21)
        R231=2.d00*(Q21Q23-Q20Q22)
        R232=2.d00*(Q22Q23+Q20Q21)
        R233=Q202-Q212-Q222+Q232

        CSM1M=-100.
        CSM2M=-100.
        I1MIN = -1
        I2MIN = -1


        DO I=0 , NRB_SITES(itypei)-1

            Indx = Nrb_sites_max*Ndim*itypei + I*Ndim

            PX=PATCH(Indx)
            PY=PATCH(Indx+1)
            PZ=PATCH(Indx+2)

            PXP = R111*PX+R112*PY+R113*PZ
            PYP = R121*PX+R122*PY+R123*PZ
            PZP = R131*PX+R132*PY+R133*PZ

            PX1(I) = PXP
            PY1(I) = PYP
            PZ1(I) = PZP

            CSM1 = (PXP*DXP+ PYP*DYP + PZP*DZP) / R12

            COS1(I) = CSM1

            If (CSM1>CSM1M) Then
               CSM1M=CSM1
               I1MIN = I
            End If
 
            PX=ref_tor_vec(Indx)
            PY=ref_tor_vec(Indx+1)
            PZ=ref_tor_vec(Indx+2)

            PX1_ref(I) = R111*PX+R112*PY+R113*PZ
            PY1_ref(I) = R121*PX+R122*PY+R123*PZ
            PZ1_ref(I) = R131*PX+R132*PY+R133*PZ

        EndDo 
       

        DO I=0 , NRB_SITES(itypej)-1

            Indx = Nrb_sites_max*Ndim*itypej + I*Ndim

            PX=PATCH(Indx)
            PY=PATCH(Indx+1)
            PZ=PATCH(Indx+2)

            PXP = R211*PX+R212*PY+R213*PZ
            PYP = R221*PX+R222*PY+R223*PZ
            PZP = R231*PX+R232*PY+R233*PZ

            PX2(I) = PXP
            PY2(I) = PYP
            PZ2(I) = PZP

            CSM2 = -(PXP*DXP+ PYP*DYP + PZP*DZP) / R12

            COS2(I) = CSM2

            If (CSM2>CSM2M) Then
               CSM2M=CSM2
               I2MIN = I
            End If

            PX=ref_tor_vec(Indx)
            PY=ref_tor_vec(Indx+1)
            PZ=ref_tor_vec(Indx+2)

            PX2_ref(I) = R211*PX+R212*PY+R213*PZ
            PY2_ref(I) = R221*PX+R222*PY+R223*PZ
            PZ2_ref(I) = R231*PX+R232*PY+R233*PZ

        EndDo 

        Vangtor_max = -1.d00
        Vang_max = -1.d00
        Utor_max = -1.d00
      
        Do J3 = 0, Nrb_sites(itypei) -1
 
           v1x=pz1_ref(J3)*dyp-py1_ref(J3)*dzp
           v1y=px1_ref(J3)*dzp-pz1_ref(J3)*dxp
           v1z=py1_ref(J3)*dxp-px1_ref(J3)*dyp
           mod_v1=sqrt(v1x*v1x+v1y*v1y+v1z*v1z)
           CSOMG1 = cos1 (J3)

           If(abs(CSOMG1) > 1.0) Then
                If( abs(CSOMG1) - 1.0 < 1.d-4) Then
                   CSOMG1 = sign(1.0,CSOMG1)
                Else
                   print*, 'error in csomg1',CSOMG1
                   print*, 'itypei',itypei
                   print*, 'j3=',j3
                   print*, 'px_ref',px1(J3),py1(J3),pz1(J3)
                   stop
                Endif
           Endif
           OMG1 = acos(CSOMG1)

           Do J4 = 0, Nrb_sites(itypej) -1

              Indx2 = Itypei * Nrb_sites_max + J3 ! elijo fila
              Indx3 = Indx2 * Nrb_sites_max * Npart_types  ! pongo el indice al inicio de la fila
              Indx4 = Indx3 + itypej * Nrb_sites_max  ! me muevo en la fila hasta el atomo de tipo Itypej
              factor_energy = Vpot_matrix( Indx4+ J4)

              if ( factor_energy > 1.d-07) Then

                        
                  CSOMG2 = cos2(J4)

                  If(abs(CSOMG2) > 1.0) Then
                       If( abs(CSOMG2) - 1.0 < 1.d-4) Then
                          CSOMG2 = sign(1.0,CSOMG2)
                       Else       
                          print*, 'error in csomg2',CSOMG2
                          stop       
                       Endif      
                  Endif  

                  OMG2 = acos(CSOMG2)

                  IndxT1 = Nrb_sites_max*itypei + j3
                  IndxT2 = Nrb_sites_max*itypej + j4

                  sigma_ij=2.d0*sigma_jon(IndxT1)*sigma_jon(IndxT2)

                  VANG = dexp(-(OMG1**2.d0+OMG2**2.d0)/sigma_ij)

                  if (bool_tor .eq. 1 .and. (Nrb_sites(itypei).gt.1) .and. &
               & (Nrb_sites(itypej).gt.1) .and. Ntor_angle( IndxT1) > 0 .and.  &
               &         Ntor_angle( IndxT2 ) > 0 ) then


                     v2x=pz2_ref(J4)*dyp-py2_ref(J4)*dzp
                     v2y=px2_ref(J4)*dzp-pz2_ref(J4)*dxp
                     v2z=py2_ref(J4)*dxp-px2_ref(J4)*dyp

                     mod_v2=sqrt(v2x*v2x+v2y*v2y+v2z*v2z)

                     If(mod_v1.ne.0.0.and.mod_v2.ne.0.0) then

                         cangle=(v1x*v2x+v1y*v2y+v1z*v2z)/mod_v1/mod_v2

                         If(abs(cangle) > 1.0) Then
                              If( abs(cangle) - 1.0 < 1.d-4) Then
                                 cangle = sign(1.0,cangle)
                              Else
                                 print*, 'error in cangle',cangle
                                 stop
                              Endif
                         Endif

                         v12x=v1y*v2z-v1z*v2y
                         v12y=v1z*v2x-v1x*v2z
                         v12z=v1x*v2y-v1y*v2x
                         v12b2=(v12x*dxp+v12y*dyp+v12z*dzp)/r12/mod_v1/mod_v2

                         If(abs(v12b2) > 1.0) Then
                             If( abs(v12b2) - 1.0 < 1.d-4) Then
                                   v12b2 = sign(1.0,v12b2)
                             Else
                                   print*, 'error in v12b2',v12b2
                                   stop
                             Endif
                         Endif

                         torsion=atan2(v12b2,cangle)

                         Utor_max = -100.d0

                         IndxT1 = Nrb_sites_max*itypei + j3

                         Do K= 0, Ntor_angle( IndxT1 )-1

                               dtor = torsion-tor_angle(IndxT1*Ntor_max + K)
                               if(dtor.gt.pi) then
                                  dtor = dtor - dospi
                               elseif(dtor.lt.-pi) then
                                  dtor = dtor + dospi
                               elseif(abs(dtor).gt.dospi) then
                                  print*, 'error in fdtor'
                                  stop
                               endif

                               Utor = dexp (-dtor**2.d0/sigma_tor_jon)
                               If (Utor  > Utor_max) Utor_max = Utor 

                         EndDo

                         Vangtor = factor_energy*Vang*Utor_max
                     Else
  
                     EndIf ! If(mod_v1.ne.0.0.and.mod_v2.ne.0.0) then

                  Else ! if (bool_tor .eq. 0) : No hay torsion

                         Vangtor = factor_energy*Vang

                  EndIf
                         
                  If( Vangtor  > Vangtor_max) Then
                      Vangtor_max = Vangtor
                      Vang_max = Vang
                      Vtor_max = Utor_max
                  EndIf


              EndIf ! if ( factor_energy > 1.d-08) Then
           EndDo ! Do J4 = 0, Nrb_sites(itypej) -1
        EndDo ! Do J3 = 0, Nrb_sites(itypei) -1

        if(Vangtor_max.eq.-1.d00) then
           En= 0.0 
        else 
           VANGTOR=Vangtor_max
           FX=(1 + Tanh(a*(-xop_aux + r12)))/2.d00
           FXI=(1 - Tanh(a*(-xop_aux + r12)))/2.d00
           R126=(RSQ/sigma_LJ_aux/sigma_LJ_aux)**3.d0
           R1212=R126*R126
           VRAD=-4.d0*((r1212 - r126)/(r1212*r126) + VL0_aux)
           En = VRAD * (FXI + VANGTOR * FX) 
        endif

   END IF  ! R12 < rangep

End Subroutine ener
