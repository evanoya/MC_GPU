
        program lifting2
        implicit none
        integer ni,nip,nim
        integer i,natoms,i3,nmax,k,j,kk,nbonds,nassign,nvmax
        integer j3,j6,i6,indxj,i1,nzero,n_proj,nun
        parameter(nmax=100000,nvmax=20)
        integer iassign(nmax),neigh(nmax),nneigh(nmax,nvmax)
        integer ncount,unassigned,ibin,ncounter(5000),nbin
        integer nassigned, nmisalign, nisolated,imin
        integer num_tot,current_layer,layer(nmax)
        integer ll, nconflict, conflict(nmax,10)
        integer nlayer,nalayer, head(nmax)
	integer niter,iloop, nuna

        double precision x(3*nmax),bond(200,4),b6D(200,7)
        double precision xn(nmax,nvmax),yn(nmax,nvmax),zn(nmax,nvmax)
        double precision dn(nmax,nvmax),x6D(nmax,7)
        double precision xi,yi,zi,dx,dy,dz,dr,drn,proj
        double precision xp,yp,zp,xpp,ypp,zpp
        double precision a6D,rcut,ctheta,norma,norma2,delta2
        double precision d6D(6),x6D1,x6D2,x6D3,x6D4,x6D5,x6D6
        double precision bond_length,tau,dperp(5000),dens(5000)
        double precision dmin,dr2,xtest(6),delta,r_out,r_in,vol
	double precision vol2,rperp,anglex,angley,anglez,xpi
        double precision ax(3,3),ay(3,3),az(3,3)
        logical disagree,ana
        character*2 el

        read(*,*)
        read(*,*) anglex
        read(*,*)
        read(*,*) angley
        read(*,*)
        read(*,*) anglez

        xpi=acos(-1.d0)

        ! rotate y
        angley= angley*xpi/180.d0
        ay(1,1)=dcos(angley)
        ay(1,2)=0.d0
        ay(1,3)=dsin(angley)
        ay(2,1)=0.d0
        ay(2,2)=1.d0
        ay(2,3)=0.d0
        ay(3,1)=-dsin(angley)
        ay(3,2)=0.d0
        ay(3,3)=dcos(angley)

        ! rotate x
        anglex=anglex*xpi/180.d0
        ax(1,1)=1.d0
        ax(1,2)=0.d0
        ax(1,3)=0.d0
        ax(2,1)=0.d0
        ax(2,2)=dcos(anglex)
        ax(2,3)=-dsin(anglex)
        ax(3,1)=0.d0
        ax(3,2)=dsin(anglex)
        ax(3,3)=dcos(anglex)

        ! rotate z
        anglez=anglez*xpi/180.d0
        az(1,1)=dcos(anglez)
        az(1,2)=-dsin(anglez)
        az(1,3)=0.d0
        az(2,1)=dsin(anglez)
        az(2,2)=dcos(anglez)
        az(2,3)=0.d0
        az(3,1)=0.d0
        az(3,2)=0.d0
        az(3,3)=1.d0


        tau=(1.d0+dsqrt(5.d0))/2.d0
        open(20,file='coords.xyz')
        read(20,*) natoms
        read(20,*)
        open(30,file='coords-rotated.xyz')
        write(30,*) natoms
        write(30,*)
        do i=1,natoms
             i3=3*(i-1)
             read(20,*) el,xi,yi,zi
!
              xp=xi*ay(1,1)+yi*ay(1,2)+zi*ay(1,3)
              yp=xi*ay(2,1)+yi*ay(2,2)+zi*ay(2,3)
              zp=xi*ay(3,1)+yi*ay(3,2)+zi*ay(3,3)
!
              xpp=xp*ax(1,1)+yp*ax(1,2)+zp*ax(1,3)
              ypp=xp*ax(2,1)+yp*ax(2,2)+zp*ax(2,3)
              zpp=xp*ax(3,1)+yp*ax(3,2)+zp*ax(3,3)
!
              x(i3+1)=xpp*az(1,1)+ypp*az(1,2)+zpp*az(1,3)
              x(i3+2)=xpp*az(2,1)+ypp*az(2,2)+zpp*az(2,3)
              x(i3+3)=xpp*az(3,1)+ypp*az(3,2)+zpp*az(3,3)
              write(30,'("C1",3f16.8)') x(i3+1),x(i3+2),x(i3+3)

        enddo
        close(20)
        close(30)
        open(20,file='bonds.xyz')
        read(20,*) nbonds
        do i=1,nbonds
         read(20,*) k,(bond(i,j),j=1,4),
     +        (b6D(i,j),j=1,7),a6D
         write(*,'(i3,4f10.5,6f5.1,2f10.5)') k,(bond(k,j),j=1,4),
     +        (b6D(k,j),j=1,7),a6D
        enddo
        close(20)

        rcut=1.3d0
        ! construir lista de vecinos
        do i=1,natoms
            neigh(i)=0
            do j=1,nvmax
               nneigh(i,j)=0
            enddo
        enddo
        do i=1,natoms
            i3=3*(i-1)
            xi=x(i3+1)
            yi=x(i3+2)
            zi=x(i3+3)
            do j=1,natoms
             if(i.ne.j) then
                 j3=3*(j-1)
                 dx=x(j3+1)-xi
                 dy=x(j3+2)-yi
                 dz=x(j3+3)-zi
                 dr=dsqrt(dx*dx+dy*dy+dz*dz)
                 if(dr.le.rcut) then
                    neigh(i)=neigh(i)+1
                    nneigh(i,neigh(i))=j
                    xn(i,neigh(i))=dx
                    yn(i,neigh(i))=dy
                    zn(i,neigh(i))=dz
                    dn(i,neigh(i))=dr
           endif
                 endif
            enddo
        enddo
        !do i=1,natoms
        !    print*, 'i=',i,neigh(i)
        !enddo
        print*, 'lista de vecinos construida'
        ! me centro en un atomo
        dmin=100000.
        do i=1,natoms
            i3=3*(i-1)
            xi=x(i3+1)
            yi=x(i3+2)
            zi=x(i3+3)
            dr2=xi*xi+yi*yi+zi*zi
            if(dr2.lt.dmin.and.neigh(i).gt.0) then
               dmin=dr2
               imin=i
           endif
        enddo
        i1=imin 
        print*, 'dmin',dmin,imin
        print*, 'origin at atom',i1
        print*, 'xi',x(3*(i1-1)+1),x(3*(i1-1)+2),x(3*(i1-1)+3)

        do i=1,natoms
           do j=1,10
              conflict(i,j)=0
           enddo
        enddo

	niter=2

	do iloop=1,niter

	print*, '55555555555555555555555555555555555'
	print*, ''
	print*, ''
        print*, 'doing iteration number',iloop
	print*, ''
	print*, ''
	print*, '55555555555555555555555555555555555'

        do i=1,natoms
           iassign(i)=-1
        enddo
        iassign(i1)=1
        i6=6*(i1-1)
        i3=3*(i1-1)
        print*, 'x',x(i3+1),x(i3+2),x(i3+3)
	do i=1,natoms
           x6D(i,1)=0.d0
           x6D(i,2)=0.d0
           x6D(i,3)=0.d0
           x6D(i,4)=0.d0
           x6D(i,5)=0.d0
           x6D(i,6)=0.d0
        enddo
        nassign=1
        !bond_length=1.0000*(1.d0+2.d0*tau)/dsqrt(4.d0+2.0*tau)
        !bond_length=0.714399136159273*(1.d0+2.d0*tau)/dsqrt(4.d0+2.0*tau)
        bond_length=2.d0**(1.d0/6.d0)
        do i=1,natoms
           layer(i)=-20
           head(i)=-1  
        enddo
        layer(i1)=1
        head(i1)=0
        current_layer=1
        nzero=0

        do kk=1,200
            print*, '**************'
            print*, 'kk=',kk
            print*, 'current_layer=',current_layer
            print*,  'nassign=',nassign
            print*, 'isolated atoms',nzero
            print*, 'number of atoms in the previous layer', nlayer
            print*, 'number assignments in the previous layer', nalayer
            print*, 'total number of atoms', natoms
            nlayer=0
            nalayer=0
            ni=0
            nip=0
            nim=0

            do i=1,natoms
              if(neigh(i).eq.0) then
               if(iassign(i).eq.-1) then
                iassign(i)=0
                nzero=nzero+1
                x6D(i,1)=0.d0
                x6D(i,2)=0.d0
                x6D(i,3)=0.d0
                x6D(i,4)=0.d0
                x6D(i,5)=0.d0
                x6D(i,6)=0.d0
               endif
              else
               if(iassign(i).eq.1.and.layer(i).eq.current_layer) then
                 ana=.false.
                 if(iloop.eq.1) then
                   ana=.true.
                 else
                   if(conflict(i,iloop-1).eq.0) ana=.true.
                 endif
                 if(ana) then
                  do j=1,neigh(i)
                      ! chequeo si es paralelo a algun bond
                      dx=xn(i,j)
                      dy=yn(i,j)
                      dz=zn(i,j)
                      dr=dn(i,j)
                      indxj=nneigh(i,j)
                      ncount=0
                      do k=1,nbonds
                         proj=dx*bond(k,1)+dy*bond(k,2)+dz*bond(k,3)
                         ctheta=proj/dr/bond(k,4)
                         proj=proj/bond(k,4)
                         n_proj=nint(proj/bond_length)
                         if(dabs(ctheta).ge.0.97) then
                            if(abs(n_proj).ne.1) print*, 'n_proj',n_proj
                            ncount=ncount+1
                            xtest(1)=x6D(i,1)+b6D(k,1)*n_proj
                            xtest(2)=x6D(i,2)+b6D(k,2)*n_proj
                            xtest(3)=x6D(i,3)+b6D(k,3)*n_proj
                            xtest(4)=x6D(i,4)+b6D(k,4)*n_proj
                            xtest(5)=x6D(i,5)+b6D(k,5)*n_proj
                            xtest(6)=x6D(i,6)+b6D(k,6)*n_proj
                            if(iassign(indxj).eq.1) then ! check if coincides with previous assignment
                               disagree=.false.
                               do ll=1,6
                                 if(xtest(ll).ne.x6D(indxj,ll)) disagree=.true.
                               enddo
                               if(disagree) then 
                                 conflict(indxj,iloop)=1
                               endif
                            else ! do the assignment
                               nassign=nassign+1
                               nalayer=nalayer+1
                               layer(indxj)=current_layer+1
                               nlayer=nlayer+1
                               iassign(indxj)=1
                               if(layer(indxj).eq.current_layer) then
                                  ni=ni+1
                               elseif(layer(indxj).eq.current_layer+1) then
                                  nip=nip+1
                               elseif(layer(indxj).eq.current_layer-1) then
                                  nim=nim+1
                               else
                                 !print*, 'layer(i)',layer(i)
                                 !print*, 'current_later',current_layer
                                 !print*, 'layer(indxj)',layer(indxj)
                                 !print*, 'something weird is going on'
                               endif
                         
                               x6D(indxj,1)=xtest(1)
                               x6D(indxj,2)=xtest(2)
                               x6D(indxj,3)=xtest(3)
                               x6D(indxj,4)=xtest(4)
                               x6D(indxj,5)=xtest(5)
                               x6D(indxj,6)=xtest(6)
                               !print*, 'new assignment',layer(i),layer(indxj)
                               head(indxj)=i
                            endif
                         
                         endif

                      enddo
                      if(ncount.gt.1) then
                        print*, 'ncount=',ncount
                        stop
                      elseif(ncount.eq.0) then
                        if(iassign(indxj).eq.-1) iassign(indxj)=0
                      endif
                  enddo
                 endif ! if(iloop.eq.1.or.(iloop.gt.1.and.conflict(i,iloop-1).eq.0)) then
               endif ! if iassign(i).eq.1.and.layer(i).eq.current_layer
              endif ! if nneigh(i) > 0
            enddo
            !print*, 'after analyzing current layer I found'
            !print*, 'new assigments with the current layer',ni
            !print*, 'new assigments with the next layer',nip
            !print*, 'new assigments with the previous layer',nim

            if(ni.eq.0.and.nip.eq.0.and.nim.eq.0) then
               print*, 'finished after ',kk,'iterations'
               print*, 'assigned',nassign,'out of',natoms
               open(40,file='6D.xyz')
                  do i=1,natoms
                     i3=3*(i-1) 
                     i6=6*(i-1) 
                     write(40,'(3f10.5,6f6.1)') x(i3+1),x(i3+2),x(i3+3),x6D(i,1),
     +             x6D(i,2),x6D(i,3),x6D(i,4),x6D(i,5),x6D(i,6)
                  enddo
               close(40) 
               call write_lifted(natoms,x,x6D,iassign,conflict,neigh,iloop)
               go to 100
            endif

            if(nassign.eq.natoms) then
               print*, 'finished after ',kk,'iterations'
               print*, 'assigned',nassign,'out of',natoms
               open(40,file='6D.xyz')
                  do i=1,natoms
                     i3=3*(i-1) 
                     i6=6*(i-1) 
                     write(40,'(3f10.5,6f6.1)') x(i3+1),x(i3+2),x(i3+3),x6D(i,1),
     +             x6D(i,2),x6D(i,3),x6D(i,4),x6D(i,5),x6D(i,6)
                  enddo
               close(40) 
               go to 100
            endif

            current_layer=current_layer+1

        enddo
 100    continue

	enddo



        nun=0
        nuna=0
        do i=1,natoms
            if(iassign(i).eq.-1) then
               nun=nun+1
            elseif(iassign(i).eq.1) then
               nuna=nuna+1
            endif
        enddo
        print*, 'Particles that I have never tried assigning',nun
        print*, 'Assigned particles',nuna

        delta=1.1d0
        delta2=0.02
        nbin=5000

        do i=1,nbin
           dperp(i)=0.d0
           ncounter(i)=0
           dens(i)=0.d0
        enddo
        num_tot=0
        do i=1,natoms-1
            x6D1=x6D(i,1)
            x6D2=x6D(i,2)
            x6D3=x6D(i,3)
            x6D4=x6D(i,4)
            x6D5=x6D(i,5)
            x6D6=x6D(i,6)
            if(iassign(i).eq.1) then
            do j=i+1,natoms
                d6D(1)=x6D(j,1)-x6D1
                d6D(2)=x6D(j,2)-x6D2
                d6D(3)=x6D(j,3)-x6D3
                d6D(4)=x6D(j,4)-x6D4
                d6D(5)=x6D(j,5)-x6D5
                d6D(6)=x6D(j,6)-x6D6
                if(iassign(j).eq.1) then
                  call project(d6D,norma)
                  call project_perp(d6D,norma2)
                  
                  ibin=int(norma/delta)+1
                  dperp(ibin)=dperp(ibin)+norma2
                  ncounter(ibin)=ncounter(ibin)+1
                  ibin=int(norma2/delta2)+1
                  dens(ibin)=dens(ibin)+1.d0
	          num_tot=num_tot+1
                endif
            enddo
            endif
        enddo
        print*, 'num_tot=',num_tot
        print*, 'natoms=',natoms
        open(80,file='histo_dperp_pairs.dat')
        do i=1,nbin
           rperp=dble(i-1)*delta2+delta2/2.0
           r_out=dble(i)*delta2
           r_in=dble(i-1)*delta2
           if(i.eq.1) r_in=0.d0
           vol=(4./3.)*acos(-1.0)*(r_out**3.0-r_in**3.0)
           r_in=dble(i-1)*delta2+delta2*0.5d0
           vol2=4.*acos(-1.0)*r_in**2.0*delta2
           write(80,'(4e16.8)') rperp,dens(i),dens(i)/vol,dens(i)/vol/num_tot
        enddo
        close(80)

        do i=1,nbin
            dens(i)=0.d0
        enddo
        open(60,file='histo_dpar_dperp_pairs.dat')
        do i=1,nbin
           if(ncounter(i).gt.0) then
           dperp(i)=dperp(i)/dble(ncounter(i))
           ibin=int(dperp(i)/delta2)+1
           dens(i)=dens(i)+dble(ncounter(i))
           write(60,*) dble(i-1)*delta+delta/2.d0,dperp(i),ncounter(i)
           endif
        enddo
        open(80,file='histo_dperp_pairs2.dat')
        do i=1,nbin
           rperp=dble(i-1)*delta2+delta2/2.0
           r_out=dble(i)*delta2
           r_in=dble(i-1)*delta2
           if(i.eq.1) r_in=0.d0
           vol=(4./3.)*acos(-1.0)*(r_out**3.0-r_in**3.0)
           r_in=dble(i-1)*delta2+delta2*0.5d0
           vol2=4.*acos(-1.0)*r_in**2.0*delta2
           write(80,'(4e16.8)') rperp,dens(i),dens(i)/vol,dens(i)/vol2
        enddo
        close(80)
        print*, 'delta2=',delta2
        stop
        end
c
c*************************************************************************
c
        subroutine project(b6D,norm)
        implicit none
        integer i,j
        double precision qpar(3,6),b6d(6),norm
        double precision sum,tau,a

        tau=(1.d0+dsqrt(5.d0))/2.d0

        qpar(1,1)=0.d0
        qpar(2,1)=tau
        qpar(3,1)=1.d0

        qpar(1,2)=tau
        qpar(2,2)=1.d0
        qpar(3,2)=0.d0

        qpar(1,3)=1.d0
        qpar(2,3)=0.d0
        qpar(3,3)=tau

        qpar(1,4)=0.d0
        qpar(2,4)=-tau
        qpar(3,4)=1.d0

        qpar(1,5)=-tau
        qpar(2,5)=1.d0
        qpar(3,5)=0.d0

        qpar(1,6)=1.d0
        qpar(2,6)=0.d0
        qpar(3,6)=-tau

        norm=0.d0
        do i=1,3
           sum=0.d0
           do j=1,6
              sum=sum+qpar(i,j)*b6d(j)
           enddo
           norm=norm+sum*sum
        enddo
        a=(2.d0**(1.d0/6.0))*dsqrt(2.d0*(2.d0+tau))/2.d0
        norm=a*dsqrt(norm)/dsqrt(4.d0+2.d0*tau)

        return
        end

c
c*************************************************************************
c
        subroutine project_perp(b6D,norm)
        implicit none
        integer i,j
        double precision qperp(3,6),b6d(6),norm
        double precision sum,tau,a,xp,xq

        tau=(1.d0+dsqrt(5.d0))/2.d0
        xp=2.d0
        xq=3.d0
        xp=1.d0
        xq=tau

        qperp(1,1)=0.d0
        qperp(2,1)=-xp !-1.d0
        qperp(3,1)=xq ! tau

        qperp(1,2)=-xp ! -1.0
        qperp(2,2)=xq ! tau 
        qperp(3,2)=0.d0

        qperp(1,3)=xq ! tau 
        qperp(2,3)=0.d0
        qperp(3,3)=-xp ! -1.d0

        qperp(1,4)=0.d0
        qperp(2,4)=xp ! 1.d0
        qperp(3,4)=xq ! tau 

        qperp(1,5)=xp ! 1.d0
        qperp(2,5)=xq ! tau 
        qperp(3,5)=0.d0

        qperp(1,6)=xq! tau 
        qperp(2,6)=0.d0
        qperp(3,6)=xp ! 1.d0 

        norm=0.d0
        do i=1,3
           sum=0.d0
           do j=1,6
              sum=sum+qperp(i,j)*b6d(j)
           enddo
           norm=norm+sum*sum
        enddo
        a=(2.d0**(1.d0/6.0))*dsqrt(2.d0*(2.d0+tau))/2.d0
        norm=a*dsqrt(norm)/dsqrt(4.d0+2.d0*tau)

        return
        end

c************************************************************************
        subroutine write_lifted(natoms,x,x6D,iassign,conflict,neigh,niter)
	implicit none
	integer natoms,niter,nmax,nvmax
        parameter(nmax=100000,nvmax=20)
        integer iassign(nmax),conflict(nmax,10)
        integer unassigned,nassigned, nisolated, nconflict, nbin
        integer nmisalign,neigh(nmax)
        integer i,i3
        double precision x(3*nmax),x6D(nmax,7)
	double precision delta,delta2,d6D(6)
	double precision norma,norma2
	character*25 filename
        
        write(filename,'("coords-lifted-",I1,".xyz")')niter 
        open(80,file=filename)
        write(80,*) natoms
        write(80,*) 
        unassigned=0
        nassigned=0
        nmisalign=0
        nisolated=0
        nconflict=0
        nbin=200
        delta=1.1d0
        delta2=0.002
        do i=1,natoms
            i3=3*(i-1)       
           if(iassign(i).eq.1) then
            if(conflict(i,niter).eq.1) then
              write(80,'("C",3x,3f16.4)') x(i3+1),x(i3+2),x(i3+3)
              nconflict=nconflict+1
            else
              write(80,'("A",3x,3f16.4)') x(i3+1),x(i3+2),x(i3+3)
            endif
            d6D(1)=x6D(i,1)
            d6D(2)=x6D(i,2)
            d6D(3)=x6D(i,3)
            d6D(4)=x6D(i,4)
            d6D(5)=x6D(i,5)
            d6D(6)=x6D(i,6)
            call project(d6D,norma)
            call project_perp(d6D,norma2)
            !ibin=int(norma2/delta2)+1
            !dens(ibin)=dens(ibin)+1.d0
            nassigned=nassigned+1
           elseif(iassign(i).eq.0) then
             if(neigh(i).eq.0) then
                write(80,'("U",3x,3f16.4)') x(i3+1),x(i3+2),x(i3+3)
                nisolated=nisolated+1
             else
                write(80,'("M",3x,3f16.4)') x(i3+1),x(i3+2),x(i3+3)
                nmisalign=nmisalign+1
             endif
           else
            write(80,'("X",3x,3f16.4)') x(i3+1),x(i3+2),x(i3+3)
            unassigned=unassigned+1
           endif
        enddo
	close(80)
        print*, '**************************************'
        print*, '**************************************'
        print*, '**************************************'
	print*, ''
	print*, 'Results for iteraction number', niter
	print*, ''
        print*, 'Unassigned',unassigned
        print*, 'nassigned',nassigned
        print*, 'nmisalign',nmisalign
        print*, 'nisolated',nisolated
        prinT*, 'total unassigned',unassigned+nmisalign+nisolated
        print*, 'nconflict',nconflict
        print*, 'natoms',natoms
        print*, 'check',unassigned+nassigned+nmisalign+nisolated
	print*, ''
	print*, ''
        print*, '**************************************'
        print*, '**************************************'
        print*, '**************************************'
        return
	end
c**************************************************************************
