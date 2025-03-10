	program BOOP
        implicit none
        integer natoms,i,j,k,nmax,nvmax,npoints,npmax,indx
        parameter(nmax=500000,nvmax=200,npmax=10000)
        integer neigh(nmax),nneigh(nmax,nvmax),ibin,nevents
        integer itype(nmax)
        integer irotate, iaxis,nframes,kk,ibinx,ibiny,nbmax
        integer nevents2,nevents2_II,nevents2_OO,nevents2_IO
        integer Npart_types_max,Nrb_sites_max, Ntor_max, Ntmax
        parameter(Npart_types_max=5, Nrb_sites_max =15, Ntor_max = 6, Ntmax=100)
        integer Nrb_sites(Npart_types_max),itipo
        integer iref_vec(Npart_types_max,Nrb_sites_max)
        integer Ntor_angle(Npart_types_max,Nrb_sites_max),Npart_types

        double precision xn_av,xnII_av,xnOO_av,xnIO_av
        double precision patchx(nmax,Nrb_sites_max)
        double precision patchy(nmax,Nrb_sites_max)
        double precision patchz(nmax,Nrb_sites_max)
        double precision x(nmax),y(nmax),z(nmax),lx,ly,lz,dum
        double precision xi,yi,zi,dx,dy,dz,dr,px,py,rcut,pnorm
        double precision dxn(nmax,nvmax),dyn(nmax,nvmax),en
        double precision dzn(nmax,nvmax),drn(nmax,nvmax)
        double precision delta_gr,delta_gri,gr(npmax)
        double precision xpi,vol,const,ralto,rbajo,cideal
        double precision distmax,angle,xri,yri,zri,matriz(3,3)
        double precision cangle,sangle,histo_BOOD_av(1000,1000)
        double precision histo_BOOD(1000,1000),delta,xp,yp,xdA
        double precision BOOD_Lambert(1000,1000),BOOD_Lambert_av(1000,1000)
        double precision BOOD_Lambert_II(1000,1000),BOOD_Lambert_II_av(1000,1000)
        double precision BOOD_Lambert_OO(1000,1000),BOOD_Lambert_OO_av(1000,1000)
        double precision BOOD_Lambert_IO(1000,1000),BOOD_Lambert_IO_av(1000,1000)
        double precision vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
        double precision sigma_jon(Npart_types_max,Nrb_sites_max)
        double precision tor_angle(Npart_types_max,Nrb_sites_max,Ntor_max)
        double precision sigma_tor_jon,en_tot
        double precision sigma_LJ(Npart_types_max,Npart_types_max)
        double precision Vl0(Npart_types_max,Npart_types_max)
        double precision rangeP(Npart_types_max,Npart_types_max)
        double precision xop(Npart_types_max,Npart_types_max)

        double precision xp2,yp2
        character*2 el
        logical overlap,bool_tor

       common /nsites/ Nrb_sites, iref_vec, Ntor_angle
        common /vpot_par/ sigma_tor_jon
        common /vpot_par2/ sigma_LJ,VL0,xop,rangeP

       common /scalar_nsites/ Npart_types, bool_tor
       common /array_vpot_par/ vpot_matrix, sigma_jon
       common /array_tor/ tor_angle
        common /coordspacth/ patchx, patchy, patchz
        common /coords/ x, y, z
        common /coordstype/ itype


        print*, '¿Quieres rotar la configuracion (0-no, 1-si)?'
        read(*,*) irotate
        print*, 'Elige el eje (x-1, y-2, z-3)'
        read(*,*) iaxis
        print*, 'Cual es el angulo de rotacion (en grados)?'
        read(*,*) angle
        if(irotate.eq.1) then
          xpi=acos(-1.d0)
          cangle=dcos(angle*xpi/180.)
          sangle=dsin(angle*xpi/180.)
          if(iaxis.eq.1) then
              matriz(1,1)=1.00
              matriz(1,2)=0.00
              matriz(1,3)=0.00
              matriz(2,1)=0.00
              matriz(2,2)=cangle
              matriz(2,3)=-sangle
              matriz(3,1)=0.00
              matriz(3,2)=sangle
              matriz(3,3)=cangle
          elseif(iaxis.eq.2) then
              matriz(1,1)=cangle
              matriz(1,2)=0.00
              matriz(1,3)=sangle
              matriz(2,1)=0.00
              matriz(2,2)=1.d0  
              matriz(2,3)=0.d0   
              matriz(3,1)=-sangle 
              matriz(3,2)=0.d0  
              matriz(3,3)=cangle
          elseif(iaxis.eq.3) then
              matriz(1,1)=cangle
              matriz(1,2)=-sangle 
              matriz(1,3)=0.00
              matriz(2,1)=sangle
              matriz(2,2)=cangle
              matriz(2,3)=0.d0   
              matriz(3,1)=0.d0
              matriz(3,2)=0.d0   
              matriz(3,3)=1.d00 
          else
              print*, 'eje de rotacion no definido'
              stop
          endif
        endif

        delta=0.01
        nbmax=int(4.d0/0.01)+1
        do i=1,nbmax
           do j=1,nbmax
               histo_BOOD_av(i,j)=0.d0
               BOOD_Lambert_av(i,j)=0.d0
               BOOD_Lambert_II_av(i,j)=0.d0
               BOOD_Lambert_OO_av(i,j)=0.d0
               BOOD_Lambert_IO_av(i,j)=0.d0
           enddo
        enddo
        print*, 'enter rcut'
        read(*,*) rcut

        Call Read_Input_Data

        xn_av=0.d0
        xnII_av=0.d0
        xnOO_av=0.d0
        xnIO_av=0.d0
        
        open(20,file='coords.dat')
        open(44,file='BondType.dat')
        !open(21,file='coords-rotated.xyz')
        lx=2000.
        print*, 'enter number of frames'
        read(*,*) nframes
        do kk=1,nframes
           print*, 'kk=',kk
           read(20,*) natoms
           read(20,*)  
           !write(21,*) natoms
           !write(21,*)  lx

        indx=0
        do i=1,natoms
            read(20,*) el,xi,yi,zi
            !xi=xi/10.d0
            !yi=yi/10.d0
            !zi=zi/10.d0
            if (irotate.eq.1) then
              xri=matriz(1,1)*xi+matriz(1,2)*yi+matriz(1,3)*zi
              yri=matriz(2,1)*xi+matriz(2,2)*yi+matriz(2,3)*zi
              zri=matriz(3,1)*xi+matriz(3,2)*yi+matriz(3,3)*zi
              xi=xri
              yi=yri
              zi=zri
            endif
            if(el.eq.'C') then
               indx=indx+1
               itype(indx)=1
               x(indx)=xi
               y(indx)=yi
               z(indx)=zi
               j=0
            !elseif(el.eq.'T') then
            !   indx=indx+1
            !   x(indx)=xi
            !   y(indx)=yi
            !   z(indx)=zi
            !   itype(indx)=2
            !   j=0
            !elseif(el.eq.'U') then
            !   indx=indx+1
            !   x(indx)=xi
            !   y(indx)=yi
            !   z(indx)=zi
            !   itype(indx)=3
            !   j=0
            !elseif(el.eq.'V') then
            !   indx=indx+1
            !   x(indx)=xi
            !   y(indx)=yi
            !   z(indx)=zi
            !   itype(indx)=4
            !   j=0
            !elseif(el.eq.'X') then
            !   indx=indx+1
            !   x(indx)=xi
            !   y(indx)=yi
            !   z(indx)=zi
            !   itype(indx)=5
            !   j=0
            else
               j=j+1
               patchx(indx,j)=(xi-x(indx))/0.35d0
               patchy(indx,j)=(yi-y(indx))/0.35d0
               patchz(indx,j)=(zi-z(indx))/0.35d0
               pnorm=dsqrt(patchx(indx,j)**2.+patchy(indx,j)**2.+patchz(indx,j)**2.)
               if(dabs(pnorm-1.d0).gt.1.d-03) then
                 print*, 'norm deviates too much from unity',pnorm
                 stop
               endif
               patchx(indx,j)=patchx(indx,j)/pnorm
               patchy(indx,j)=patchy(indx,j)/pnorm
               patchz(indx,j)=patchz(indx,j)/pnorm
            endif
            !write(21,*) el,xi,yi,zi
            !write(21,*) el,xi*10.,yi*10.,zi*10.
        enddo
        print*, 'natoms=',indx
        natoms=indx
        !go to 83

        npoints=1000
        do i=1,npoints
              gr(i)=0.d0
        enddo
        delta_gr=0.1
        delta_gri=1.d0/delta_gr
c
        do i=1,natoms
           neigh(i)=0
        enddo
        ly=lx
        lz=lx
c
        do i=1,nbmax
          do j=1,nbmax
            histo_BOOD(i,j)=0.d0                           
            BOOD_Lambert(i,j)=0.d0                           
            BOOD_Lambert_II(i,j)=0.d0                           
            BOOD_Lambert_OO(i,j)=0.d0                           
            BOOD_Lambert_IO(i,j)=0.d0                           
          enddo
        enddo
        nevents=0
        nevents2=0
        nevents2_II=0
        nevents2_OO=0
        nevents2_IO=0
        distmax=0.d0
        en_tot=0.d0
        do i=1,natoms-1
             xi=x(i)
             yi=y(i)
             zi=z(i)
             dr=xi*xi+yi*yi+zi*zi
             if(dr.gt.distmax) distmax=dr
             if(dr.lt.20.d0**2.d0) then
             
             do j=i+1,natoms
                if(i.ne.j) then
                   dx=x(j)-xi
                   dy=y(j)-yi
                   dz=z(j)-zi
                   dx=dx-lx*dnint(dx/lx)
                   dy=dy-ly*dnint(dy/ly)
                   dz=dz-lz*dnint(dz/lz)
                   dr=dsqrt(dx*dx+dy*dy+dz*dz)
                   if(dr.lt.lx/2.d0) then
                     ibin=int(dr*delta_gri)+1
                     gr(ibin)=gr(ibin)+2.d0
                   endif

                   if(dr.lt.rcut) then

                      call ener(i,j,en,itipo,overlap)

                      en_tot=en_tot+en

                      if(en.lt.-0.20) then
                         neigh(i)=neigh(i)+1
                         neigh(j)=neigh(j)+1
                         nneigh(i,neigh(i))=j
                         nneigh(j,neigh(j))=i
                         dxn(i,neigh(i))=dx
                         dyn(i,neigh(i))=dy
                         dzn(i,neigh(i))=dz
                         drn(i,neigh(i))=dr
                         dxn(j,neigh(j))=-dx
                         dyn(j,neigh(j))=-dy
                         dzn(j,neigh(j))=-dz
                         drn(j,neigh(j))=dr
                         dx=dx/dr
                         dy=dy/dr
                         dz=dz/dr
                         px=dx/(1.d0-dz)
                         py=dy/(1.d0-dz)
                         if(px*px+py*py< =1.d0) then
                           ! paso a esfericas
                            ibinx=int((1.d0+px)/delta)+1
                            ibiny=int((1.d0+py)/delta)+1
                            histo_BOOD(ibinx,ibiny)=histo_BOOD(ibinx,ibiny)+1 
                            nevents=nevents+1
                         endif
                         ! Lambert projection
                         px=dsqrt(2.d0/(1.d0-dz))*dx
                         py=dsqrt(2.d0/(1.d0-dz))*dy
                         ibinx=int((2.d0+px)/delta)+1
                         ibiny=int((2.d0+py)/delta)+1
                         BOOD_Lambert(ibinx,ibiny)=BOOD_Lambert(ibinx,ibiny)+1
                         nevents2=nevents2+1
                         if(itipo.eq.1) then
                            BOOD_Lambert_II(ibinx,ibiny)=BOOD_Lambert_II(ibinx,ibiny)+1
                            nevents2_II=nevents2_II+1
                         elseif(itipo.eq.2) then
                            BOOD_Lambert_OO(ibinx,ibiny)=BOOD_Lambert_OO(ibinx,ibiny)+1
                            nevents2_OO=nevents2_OO+1
                         elseif(itipo.eq.3) then
                            BOOD_Lambert_IO(ibinx,ibiny)=BOOD_Lambert_IO(ibinx,ibiny)+1
                            nevents2_IO=nevents2_IO+1
                         else
                             print*, 'unknown bond type',itipo
                             stop
                         endif

                         dx=-dx
                         dy=-dy
                         dz=-dz
                         px=dx/(1.d0-dz)
                         py=dy/(1.d0-dz)
                         if(px*px+py*py< 1.d0) then
                            ibinx=int((1.d0+px)/delta)+1
                            ibiny=int((1.d0+py)/delta)+1
                            histo_BOOD(ibinx,ibiny)=histo_BOOD(ibinx,ibiny)+1 
                            nevents=nevents+1
                         endif
                         ! Lambert projection
                         px=dsqrt(2.d0/(1.d0-dz))*dx
                         py=dsqrt(2.d0/(1.d0-dz))*dy
                         ibinx=int((2.d0+px)/delta)+1
                         ibiny=int((2.d0+py)/delta)+1
                         BOOD_Lambert(ibinx,ibiny)=BOOD_Lambert(ibinx,ibiny)+1
                         nevents2=nevents2+1
                         if(itipo.eq.1) then
                            BOOD_Lambert_II(ibinx,ibiny)=BOOD_Lambert_II(ibinx,ibiny)+1
                            nevents2_II=nevents2_II+1
                         elseif(itipo.eq.2) then
                            BOOD_Lambert_OO(ibinx,ibiny)=BOOD_Lambert_OO(ibinx,ibiny)+1
                            nevents2_OO=nevents2_OO+1
                         elseif(itipo.eq.3) then
                            BOOD_Lambert_IO(ibinx,ibiny)=BOOD_Lambert_IO(ibinx,ibiny)+1
                            nevents2_IO=nevents2_IO+1
                         else
                             print*, 'unknown bond type',itipo
                             stop
                         endif

                      endif !(en.lt.-0.20) then
                   endif
                endif
             enddo
             endif
        enddo
        do i=1,natoms
            print*, 'neigh',i,neigh(i)
        enddo

             print*, '***************************************'
             print*, '***************************************'
             print*, '***************************************'
             print*, 'ener',en_tot,en_tot/dble(natoms)
             print*, 'nevents2',nevents2,nevents2_II+nevents2_OO+nevents2_IO
             print*, 'nevents2_II (intra icos)',nevents2_II
             print*, 'nevents2_OO (inter icos)',nevents2_OO
             print*, 'nevents2_IO',nevents2_IO

             xn_av=(xn_av*dble(kk-1)+dble(nevents2))/dble(kk)
             xnII_av=(xnII_av*dble(kk-1)+dble(nevents2_II))/dble(kk)
             xnOO_av=(xnOO_av*dble(kk-1)+dble(nevents2_OO))/dble(kk)
             xnIO_av=(xnIO_av*dble(kk-1)+dble(nevents2_IO))/dble(kk)
             write(44,'(3f10.4,2x,4f12.1)') xnII_av/xn_av,xnOO_av/xn_av,
     +            xnIO_av/xn_av,xn_av,xnII_av,xnOO_av,xnIO_av

             ! normalizo
             do k=1,nbmax
                do j=1,nbmax
                   histo_BOOD(k,j)=histo_BOOD(k,j)/dble(nevents)
                   BOOD_Lambert(k,j)=BOOD_Lambert(k,j)/dble(nevents2)
                   BOOD_Lambert_II(k,j)=BOOD_Lambert_II(k,j)/dble(nevents2_II)
                   BOOD_Lambert_OO(k,j)=BOOD_Lambert_OO(k,j)/dble(nevents2_OO)
                   BOOD_Lambert_IO(k,j)=BOOD_Lambert_IO(k,j)/dble(nevents2_IO)

                   histo_BOOD_av(k,j)=histo_BOOD_av(k,j)+histo_BOOD(k,j)
                   BOOD_Lambert_av(k,j)=BOOD_Lambert_av(k,j)+BOOD_Lambert(k,j)
                   BOOD_Lambert_II_av(k,j)=BOOD_Lambert_II_av(k,j)+BOOD_Lambert_II(k,j)
                   BOOD_Lambert_OO_av(k,j)=BOOD_Lambert_OO_av(k,j)+BOOD_Lambert_OO(k,j)
                   BOOD_Lambert_IO_av(k,j)=BOOD_Lambert_IO_av(k,j)+BOOD_Lambert_IO(k,j)
                enddo
             enddo

 83          continue
        enddo
        print*, 'main loop finished'
c
        open(44,file='BOOD-av.dat')
        open(45,file='BOOD-Lambert-av.dat')
        do k=1,nbmax
           xp=-1.d0+dble(k-1)*delta
           xp2=-2.d0+dble(k-1)*delta
           do j=1,nbmax
              yp=-1.d0+dble(j-1)*delta
              yp2=-2.d0+dble(j-1)*delta
              xdA=(4.d0*delta**2.d0)/(1.d0+xp**2.d0+yp**2.0)**2.d0
              histo_BOOD_av(k,j)=histo_BOOD_av(k,j)/dble(nframes)/xdA
              BOOD_Lambert_av(k,j)=BOOD_Lambert_av(k,j)/dble(nframes)
              BOOD_Lambert_II_av(k,j)=BOOD_Lambert_II_av(k,j)/dble(nframes)
              BOOD_Lambert_OO_av(k,j)=BOOD_Lambert_OO_av(k,j)/dble(nframes)
              BOOD_Lambert_IO_av(k,j)=BOOD_Lambert_IO_av(k,j)/dble(nframes)
              if(xp2*xp2+yp2*yp2<=4.d0) then
              write(45,'(6f16.8)') xp2,yp2,BOOD_Lambert_av(k,j),BOOD_Lambert_II_av(k,j),
     +              BOOD_Lambert_OO_av(k,j),BOOD_Lambert_IO_av(k,j)
              endif
              if(xp*xp+yp*yp<=1.d0) then
                 write(44,*) xp,yp,histo_BOOD_av(k,j)
              endif
           enddo
           write(44,*) 
           write(45,*) 
        enddo
        close(44)
        close(45)
        close(21)
c
        distmax=dsqrt(distmax)
        print*, 'distmax=',distmax
        xpi=acos(-1.d0)
        vol=(4./3.)*xpi*distmax**(3.d0/2.d0)
        vol=lx*ly*lz
        const=4.00*xpi*dble(natoms)/3.d00/vol
        open(20,file='PDF.dat')
        npoints=int((lx/2.d0)*delta_gri)
        do i=1,npoints
             rbajo=dfloat(i-1)*delta_gr
             ralto=rbajo+delta_gr
             dr=rbajo+delta_gr/2.d00
             cideal=const*(ralto*ralto*ralto-rbajo*rbajo*rbajo)
             gr(i)=gr(i)/dble(natoms)/cideal
             write(20,*) dr,gr(i)
        enddo

        open(34,file='BOOD-sphere.dat')
        open(35,file='BOOD-half-sphere.dat')
        open(33,file='BOOD-disk.dat')
        do i=1,natoms
             xi=x(i)
             yi=y(i)
             zi=z(i)
             dr=dsqrt(xi*xi+yi*yi+zi*zi)
             if(dr.lt.distmax-10.d0) then
                do j=1,neigh(i)
                      dr=drn(i,j)
                      dx=dxn(i,j)/dr
                      dy=dyn(i,j)/dr
                      dz=dzn(i,j)/dr
                      px=dx/(1.d0-dz)
                      py=dy/(1.d0-dz)
                      write(33,*) px,py
                      write(34,*) dx,dy,dz
                      if(dz.le.0.d0) write(35,*) dx,dy,dz
                    
                enddo
             endif
        enddo

        stop
        end
c**********************************************************************
       subroutine ener(i,j,en,itipo,overlap)
       implicit none
        integer Npart_types_max,Nrb_sites_max, Ntor_max, nmax
        parameter(Npart_types_max=5, Nrb_sites_max =15, Ntor_max = 6)
        parameter(nmax=500000)
       integer Npart_types,bool_tor, i, j, itipo
       integer Nrb_sites(Npart_types_max)
       integer iref_vec(Npart_types_max,Nrb_sites_max)
       integer Ntor_angle(Npart_types_max,Nrb_sites_max)
       integer itypei, itypej, I1MIN, I2MIN
       integer i1,j3, J4, Indx_patch, k,j3_max,j4_max
       integer itype(nmax)
       double precision Vangtor_max, Vang_max, Utor_max, v1x, v1y, v1z
       double precision mod_v1, CSOMG1, OMG1, factor_energy, Vtor_max
       double precision CSOMG2, OMG2, VANG, sigma_ij, v2x, v2y, v2z
       double precision mod_v2, cangle, v12x, v12y, v12z, v12b2, torsion
       double precision dtor, Utor, Vangtor, FX, FXI, VRAD
       double precision pi, dospi, a, en
       double precision sigma_tor_jon
       double precision sigma_LJ(Npart_types_max,Npart_types_max)
       double precision Vl0(Npart_types_max,Npart_types_max)
       double precision sigma_LJ_aux,xop_aux,vl0_aux,rangep_aux
       double precision rangeP(Npart_types_max,Npart_types_max)
       double precision xop(Npart_types_max,Npart_types_max)

       double precision vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
       double precision sigma_jon(Npart_types_max,Nrb_sites_max)
       double precision X1, Y1, Z1, X2, Y2, Z2, DX, DY, DZ
       double precision DXP, DYP, DZP, RSQ, R12, R126, R1212
       double precision CSM1M, CSM2M, CSM1, CSM2, px, py, pz
       double precision cos1(Nrb_sites_max), cos2(Nrb_sites_max)
       double precision x(nmax),y(nmax),z(nmax), h(3,3)
       double precision patchx(nmax,Nrb_sites_max)
       double precision patchy(nmax,Nrb_sites_max)
       double precision patchz(nmax,Nrb_sites_max)
       double precision PX1_ref(Nrb_sites_max),PY1_ref(Nrb_sites_max)
       double precision PZ1_ref(Nrb_sites_max),PX2_ref(Nrb_sites_max)
       double precision PY2_ref(Nrb_sites_max),PZ2_ref(Nrb_sites_max)
       double precision tor_angle(Npart_types_max,Nrb_sites_max,Ntor_max)
       logical overlap

        ! parameters for calculation of vpot
        common /scalar_nsites/ Npart_types, bool_tor
        common /vpot_par/ sigma_tor_jon
        common /vpot_par2/ sigma_LJ,VL0,xop,rangeP

        common /array_vpot_par/ vpot_matrix, sigma_jon
        common /nsites/ Nrb_sites, iref_vec, Ntor_angle
       common /array_tor/ tor_angle

        ! coordinates
        common /coords/ x, y, z
        common /coordspacth/ patchx, patchy, patchz
        common /coordstype/ itype

       pi = acos(-1.d0)
       dospi = 2.d0*pi
       a = 1000.d0

       X1=x(i)
       Y1=y(i)
       Z1=z(i)

       itypei = itype(I)

       X2=x(j)
       Y2=y(j)
       Z2=z(j)

       itypej = itype(j)

       DX=X2-X1
       DY=Y2-Y1
       DZ=Z2-Z1

       RSQ=DX**2.0+DY**2.0+DZ**2.
       R12=SQRT(RSQ)

       sigma_LJ_aux = sigma_LJ(itypei,itypej)
       xop_aux = xop(itypei,itypej)
       vl0_aux = vl0(itypei,itypej)
       rangep_aux = rangep(itypei,itypej)

       if(r12.lt.0.5) then
            overlap=.true.
       elseif(r12.lt.0.98*sigma_LJ_aux) then
            R126=(RSQ/sigma_LJ_aux/sigma_LJ_aux)**3.d0
            R1212=R126*R126
            EN=-4.d0*((r1212 - r126)/(r1212*r126) + VL0_aux)
       elseif(r12.lt.rangep_aux) then
            CSM1M=-100.
            CSM2M=-100.
            I1MIN = -1
            I2MIN = -1

            do i1= 1,nrb_sites(itypei)
               px=patchx(i,i1)
               py=patchy(i,i1)
               pz=patchz(i,i1)
               CSM1 = (PX*DX+ PY*DY + PZ*DZ) / R12
               cos1(i1)=csm1
               If (CSM1>CSM1M) Then
                  CSM1M=CSM1
                  I1MIN = I1
               End If
               If( Nrb_sites(itypei) >= 2 .and. Ntor_angle(itypei,i1) >= 1) Then ! particle J has at least two patches
                 Indx_patch=iref_vec(itypei,i1)
                 PX1_ref(i1)=patchx(i,Indx_patch)
                 PY1_ref(i1)=patchy(i,Indx_patch)
                 PZ1_ref(i1)=patchz(i,Indx_patch)
               Endif
            enddo

            do i1= 1,nrb_sites(itypej)
               px=patchx(j,i1)
               py=patchy(j,i1)
               pz=patchz(j,i1)
               CSM2 = -(PX*DX+ PY*DY + PZ*DZ) / R12
               cos2(i1)=csm2

               If (CSM2>CSM2M) Then
                  CSM2M=CSM2
                  I2MIN = I1
               End If

               If( Nrb_sites(itypej) >= 2. and. Ntor_angle(itypej,i1) >= 1 ) Then ! particle J has at least two patches
                 Indx_patch=iref_vec(itypej,i1)
                 PX2_ref(i1)=patchx(j,Indx_patch) ! mirar como va lo del ref_tor_vec (read_input)
                 PY2_ref(i1)=patchy(j,Indx_patch)
                 PZ2_ref(i1)=patchz(j,Indx_patch)
               Endif
            enddo
            Vangtor_max = -1.d00
            Vang_max = -1.d00
            Utor_max = -1.d00
            itipo=1
            if(i1min.eq.4.and.i2min.eq.4) itipo=2

            do j3=1,nrb_sites(itypei)
               v1x=pz1_ref(J3)*dyp-py1_ref(J3)*dzp
               v1y=px1_ref(J3)*dzp-pz1_ref(J3)*dxp
               v1z=py1_ref(J3)*dxp-px1_ref(J3)*dyp
               mod_v1=sqrt(v1x*v1x+v1y*v1y+v1z*v1z)
               CSOMG1 = cos1 (J3)

               If(abs(CSOMG1) > 1.0) Then
                    If( abs(CSOMG1) - 1.0 < 1.d-4) Then
                       CSOMG1 = sign(1.d0,CSOMG1)
                    Else
                       print*, 'error in csomg1',CSOMG1
                       print*, 'itypei',itypei
                       print*, 'j3=',j3
                       print*, 'px_ref',px1_ref(J3),py1_ref(J3),pz1_ref(J3)
                       stop
                    Endif
               Endif
               OMG1 = acos(CSOMG1)

               Do J4 = 1, Nrb_sites(itypej)

                    factor_energy = Vpot_matrix( itypei, j3, itypej, j4)

                    if ( factor_energy > 1.d-07) Then

                        CSOMG2 = cos2(J4)
                        If(abs(CSOMG2) > 1.0) Then
                             If( abs(CSOMG2) - 1.0 < 1.d-4) Then
                                CSOMG2 = sign(1.d0,CSOMG2)
                             Else
                                print*, 'error in csomg2',CSOMG2
                                stop
                             Endif
                        Endif

                        OMG2 = acos(CSOMG2)
                        sigma_ij=2.d0*sigma_jon(itypei, j3)*sigma_jon(itypej, j4)
                        VANG = dexp(-(OMG1**2.d0+OMG2**2.d0)/sigma_ij)
                        

                        if (bool_tor.eq.1 .and. (Nrb_sites(itypei).gt.1) .and. 
     +          (Nrb_sites(itypej).gt.1) .and. Ntor_angle( itypei, j3) > 0 .and.  
     +                  Ntor_angle( itypej, j4 ) > 0 ) then

                              v2x=pz2_ref(J4)*dyp-py2_ref(J4)*dzp
                              v2y=px2_ref(J4)*dzp-pz2_ref(J4)*dxp
                              v2z=py2_ref(J4)*dxp-px2_ref(J4)*dyp
                              mod_v2=sqrt(v2x*v2x+v2y*v2y+v2z*v2z)

                              If(mod_v1.ne.0.0.and.mod_v2.ne.0.0) then

                                 cangle=(v1x*v2x+v1y*v2y+v1z*v2z)/mod_v1/mod_v2

                                 If(abs(cangle) > 1.0) Then
                                      If( abs(cangle) - 1.0 < 1.d-4) Then
                                         cangle = sign(1.d0,cangle)
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
                                           v12b2 = sign(1.d0,v12b2)
                                     Else
                                           print*, 'error in v12b2',v12b2
                                           stop
                                     Endif
                                 Endif

                                 torsion=atan2(v12b2,cangle)

                                 Utor_max = -100.d0

                                 Do K= 1, Ntor_angle( itypei, j3 )
                                    dtor = torsion-tor_angle(itypei, j3, k)
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
                              EndIf ! If(mod_v1.ne.0.0.and.mod_v2.ne.0.0) then
                        Else ! if (bool_tor .eq. 0) : No hay torsion
                              Vangtor = factor_energy*Vang
                        EndIf

                        If( Vangtor  > Vangtor_max) Then
                            Vangtor_max = Vangtor
                            Vang_max = Vang
                            Vtor_max = Utor_max
                            j3_max=j3
                            j4_max=j4
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

               if(j4_max.le.5.and.j3_max.le.5)  then
                   itipo=1
               elseif(j4_max.ge.6.and.j3_max.ge.6) then
                   itipo=3
               else
                   itipo=2
               endif
            endif
       else
            en=0.d0
       endif

       return
       end

c**********************************************************************
       subroutine Read_Input_Data !Npart_types_max,Nrb_sites_max, Ntor_max)
       implicit none
!      faltan por asignar los valores de las variables:
!                 Npart_types_max,Nrb_sites_max, Ntor_max
        integer Npart_types_max,Nrb_sites_max, Ntor_max
        parameter(Npart_types_max=5, Nrb_sites_max =15, Ntor_max = 6)
       integer istep_ini,istep_fin,Neq, Nmove, Nsave, Nsave2, Nrestart
       integer Nrb_sites(Npart_types_max),i,j,k,bool_tor
       integer iref_vec(Npart_types_max,Nrb_sites_max)
       integer Ntor_angle(Npart_types_max,Nrb_sites_max)
       integer Indx, Indx_tor, Indx_tor2, Indx2, ncount, indx_patch
       integer i1, i2, Indx3, ios, Npart_types
       double precision hmax,omax,temp0,temp1, seed, xopp
       double precision sigma_jon_aux, sigma_tor_jon 
       double precision px, py, pz, px_ref, py_ref, pz_ref
       double precision patch_vec(3*Nrb_sites_max*Npart_types_max)
       double precision ref_tor_vec(3*Nrb_sites_max*Npart_types_max)
       double precision vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
       double precision sigma_jon(Npart_types_max,Nrb_sites_max)
       double precision tor_angle(Npart_types_max,Nrb_sites_max,Ntor_max)
       double precision Vl0(Npart_types_max,Npart_types_max)
       double precision rangeP(Npart_types_max,Npart_types_max)
       double precision xop(Npart_types_max,Npart_types_max)
       double precision sigma_LJ(Npart_types_max,Npart_types_max)
       double precision rangepp,sigma_jon_ij

       logical imovie,displ_update, perfect_match
!      parametros que voy a necesitar pasar para el calculo de energia
!      integer arrays: Nrb_sites(Npart_types_max), iref_vec(Nrb_sites_max*Npart_types_max)
!                      Ntor_angle(Nrb_sites_max*Npart_types_max)
!      double scalar: VL0,xop,sigma_tor_jon,rangeP
!      integer scalar: Npart_types, bool_tor
!      double array: vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
!                    sigma_jon(Nrb_sites_max*Npart_types_max)
       common /nsites/ Nrb_sites, iref_vec, Ntor_angle
       common /vpot_par/ sigma_tor_jon
       common /vpot_par2/ sigma_LJ,VL0,xop,rangeP
       common /scalar_nsites/ Npart_types, bool_tor
       common /array_vpot_par/ vpot_matrix, sigma_jon
       common /array_tor/ tor_angle

       Character*20 base
       character*3 model

       Open(20,file='input.d')
       read(20,*) base
       read(20,*)
       read(20,*) istep_ini,istep_fin,Neq, Nmove, Nsave, Nsave2, Nrestart
       read(20,*)
       read(20,*) 
       read(20,*)
       read(20,*)   ! leave this for furture extension to free energy calculations
       read(20,*)
       read(20,*)  imovie !  logical variable to decide whether to store the movie
       read(20,*)
       read(20,*)  model
       read(20,*)
       print*, 'Modelo ',model
       if( model == 'K-F') then
             print*, 'Model not yet implemented'
       elseif( model == 'LJG') then
             read(20,*)
             read(20,*) sigma_jon_aux, sigma_tor_jon, rangepp, Bool_tor
             sigma_tor_jon = 2.d0*sigma_tor_jon*sigma_tor_jon
             Bool_tor=0
       else
             print*, 'Model not yet implemented'
       endif
       read(20,*)  !leave this for future extension to AVB moves 
       read(20,*)
       read(20,*)
       read(20,*) Npart_types
       print*, 'Number of particle types',Npart_types

       read(20,*)
       read(20,*) (Nrb_sites(j),j=1,Npart_types)

       Do i = 1, Npart_types
          print*, 'Particle type',i,'number of patches',Nrb_sites(i)
          If (Nrb_sites(i) > Nrb_sites_max) Then
            print*, 'ERROR in particle type',i
            print*, 'Number of patches',Nrb_sites(i)
            print*, 'exceeds the maximum number of patches',Nrb_sites_max
            stop
          EndIf
       Enddo

       Ntor_angle=0
       tor_angle=0.d0

       Do J= 1, Npart_types

          Indx = (J-1)*Nrb_sites_max*3
          Indx_tor = (J-1)*Nrb_sites_max
          Indx_tor2 = (J-1)*Nrb_sites_max*Ntor_max

          If( Nrb_sites(J) >= 2 ) Then ! particle J has at least two patches
              Do I= 1, Nrb_sites(J)
                 Indx2 = Indx + (i-1)*3
                 Read(20,'(3F12.8)',advance='no')  px,py,pz
                 patch_vec( Indx2+1 ) =Px
                 patch_vec( Indx2+2 ) =Py
                 patch_vec( Indx2+3 ) =Pz
                 print*, 'px,py,pz=',px,py,pz
                 Read(20,'(3F12.8,3x,I3)',advance='no',IOSTAT=ios)
     +               px, py, pz, Ntor_angle( J, I )
                 ref_tor_vec( Indx2+1 ) =Px
                 ref_tor_vec( Indx2+2 ) =Py
                 ref_tor_vec( Indx2+3 ) =Pz
                 print*, 'i, Ntor=',i,Ntor_angle( J, I ), ios
                 if(Ntor_angle( J, I ).gt.Ntor_max) then
                       print*, ''
                       print*, 'Error: :',Ntor_angle( J, I )
                       print*, 'increase Ntor_max to ',Ntor_max
                       print*, ''
                       stop
                 endif
                 if(ios.eq.1) then
                      Ntor_angle(J,i) = 0
                 else
                      Do k= 1, Ntor_angle( J,I )
                         Read(20,'(F12.6)',advance='no') tor_angle( J, I, K )
                      EndDo
                      Read(20,*)
                 endif
                 print*, ''
              EndDo
              ! una vez que he leido todos los parches, veo cual es el de referencia
              Do i=1, nrb_sites(j)
                  If(Ntor_angle(j,i).gt.0) then
                     Indx2= Indx + (i-1)*3
                     px_ref=ref_tor_vec( Indx2+1 )
                     py_ref=ref_tor_vec( Indx2+2 )
                     pz_ref=ref_tor_vec( Indx2+3 )
                     ncount=0
                     do k=1, nrb_sites(j)
                        Indx3= Indx + (k-1)*3
                        perfect_match=.true.
                        if(dabs(px_ref-patch_vec(Indx3+1)).gt.1.d-03) then
                            perfect_match=.false.
                        endif
                        if(dabs(py_ref-patch_vec(Indx3+2)).gt.1.d-03) then
                            perfect_match=.false.
                        endif
                        if(dabs(pz_ref-patch_vec(Indx3+3)).gt.1.d-03) then
                            perfect_match=.false.
                        endif
                        if(perfect_match) then
                            indx_patch=k
                            ncount=ncount+1
                        endif
                     enddo
                     if(ncount.eq.0) then
                        print*, 'ERROR in reading data'
                        print*, 'patch used as reference for torsion'
                        print*, 'does not match any of the ang. patches'
                        !stop
                        if(i.gt.1) then
                           iref_vec(j,i)=i-1
                        else
                           iref_vec(j,i)=nrb_sites(j)
                        endif
                     elseif(ncount.eq.1) then
                        iref_vec(j,i)=indx_patch
                     else
                        print*, 'ERROR in reading data'
                        print*, 'patch used as reference for torsion'
                        print*, 'matches more than one ang patch',ncount
                        stop
                     endif
                  Endif
              Enddo
          Elseif ( Nrb_sites(J) == 1) Then
            print*, 'particle type j 1',j,Nrb_sites(j)
            Read(20,'(3F12.8)')  ( patch_vec( Indx+K ), K = 1, 3 )
            print*, 'patch=',  ( patch_vec( Indx+K ), K = 1, 3)
          EndIf
       EndDo

       ! read in the interaction matrix
       Read(20,*)

       do i=1,npart_types_max
          do i1 = 1, nrb_sites_max
             do j=1,npart_types_max
                do i2=1, nrb_sites_max
                  vpot_matrix(i,i1,j,i2)=0.d0
                enddo
             enddo
          enddo
       enddo

       print*, ''
       print*, ''
       print*, 'Vpot_Matrix(type1,patch1,type2,patch2)'

       do i=1,npart_types
          do i1 = 1, nrb_sites(i)
             read(20,*) ((vpot_matrix(i,i1,j,i2),i2=1,nrb_sites(j)),j=1,npart_types)
             print*, ((vpot_matrix(i,i1,j,i2),i2=1,nrb_sites(j)),j=1,npart_types)
          enddo
       enddo

       ! asigna valores al array sigma
       If(sigma_jon_aux.gt.0.d0) then
!
          do i=1, Npart_types
              Indx = (i-1)*Nrb_sites_max
              do j=1, Nrb_sites(i)
                 sigma_jon(I,j) = sigma_jon_aux
              Enddo
          Enddo
!
       Else

          Read(20,*)
          Do I= 1, Npart_types
             Do J = 1, Nrb_sites(i)
                Read(20,'(F6.3)',advance='no') sigma_jon_ij
                sigma_jon(i,j) = sigma_jon_ij
                print*, 'sigma_LJ i,j',i,j,sigma_jon(i,j)
             EndDo
             Read(20,*)
          EndDo

          Read(20,*)
          Do I=1,Npart_types
               Do J = 1, Npart_types
                  Read(20,'(F6.3)',advance='no') sigma_jon_ij
                  sigma_LJ(i,j) = sigma_jon_ij
                  rangep(i,j) = rangepp * sigma_jon_ij
                  print*, 'i,j',i,j,rangepp,sigma_jon_ij, rangep(i,j)
                  VL0(i,j)=(sigma_jon_ij/rangep(i,j))**12.d0-(sigma_jon_ij/rangep(i,j))**6.0
                  xopp = (1.d0+sqrt(1.d0+4.d00*VL0(i,j)))/2.d0
                  xop(i,j) = sigma_jon_ij/(xopp**(1.d0/6.d0))
                  print*, 'sigma_LJ',sigma_LJ(i,j),rangep(i,j),VL0(i,j)
               EndDo
               Read(20,*)
          EndDo

       Endif

       print*, ''
       print*, 'Input file read'
       print*, ''
       close(20)

       return
       end

