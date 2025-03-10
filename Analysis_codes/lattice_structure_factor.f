	program lattice_structure_factor	
	implicit none
        integer i,j,k,kk,nx,ny,nz,indx,natoms,indxq,nsteps
        integer natoms_uc,natoms1,natoms2,nsnapshots,ll,indxj
        integer nucx,nucy,nucz,ngr,ngrid,n_ox,nsteps_block,nblock
        integer i1,i2
        double precision pi,h(3,3),x(500000),qx,qy,qz,dq,qmax
        double precision dqx,dqy,dqz,qr,cosq1,sinq1,q,sq(1000,1000)
        double precision ux,uy,uz,v_unit,b(12000),sq_av(1000,1000)
        double precision sq_block(1000,1000) 
        double precision sq2_av(1000,1000)
        double precision cosq2,sinq2
        double precision gr(10000,7),box(3),vol,dist,ralto,rbajo
        double precision cideal,xpi,rmax,deltagr,deltagri,rho,theta
        double precision const_ar,const_ox,cideal_ar,cideal_ox
        double precision const_tot,cideal_tot,delta_anglei,rcut
        double precision sigma2,sigma2_tot,sigma2_ar,sigma2_mixed
        double precision delta_angle,x_ox(200000),h_q
        character*2 label(200000)
        	
        natoms1=0
        h(1,1)=1000.       
        h(2,2)=100.       
        h(3,3)=100.       
        pi=acos(-1.d0)
        deltagr=0.10d0
        deltagri=1.d0/deltagr

        dqx=0.1! 2.d0*pi/ux  ! red reciproca (ojo de la caja de simulacion)
        dqy=0.1! 2.d0*pi/uy
        dqz=0.1! 2.d0*pi/uz
        print*, 'red reciproca:',dqx,dqy,dqz
c
        qmax=20.0d0
        nx=int(qmax/dqx)
        ny=int(qmax/dqy)
        nz=int(qmax/dqz)
        print*, 'nx,ny,nz',nx,ny,nz
        dq=dqx 
c
        do i=1,2*nx+1
          do j=1,2*ny+1
            sq_av(i,j)=0.d00
            sq_block(i,j)=0.d0
            sq2_av(i,j)=0.d0
          enddo
        enddo

        open(33,file='coords.dat')
        open(89,file='Sq-block.dat')
        read(*,*) nsnapshots,nblock
        print*, 'nsapshots',nsnapshots
        nsteps_block=nsnapshots/nblock
        print*, 'nsteps_block=',nsteps_block
        print*, 'nblock=',nblock
c
        do ll=1,nsnapshots
c
            read(33,*) natoms2
            read(33,*) 
            print*, 'll=',ll
            print*, 'natoms2=',natoms2
!
            do i=1,natoms2
                 indx=3*(i-1)
                 read(33,*) label(i),x(indx+1),x(indx+2),x(indx+3)
            enddo
!
            natoms=natoms2
            print*, 'total number of atoms',natoms

            do i=1,2*nx+1
             do j=1,2*ny+1
               sq(i,j)=0.d00
             enddo
            enddo

            i1=0
            do i=-nx,nx
              qx=dble(i)*dqx
              i1=i1+1
              if(mod(i,50).eq.0) print*, 'i=',i
              i2=0
              do j=-ny,ny
                qy=dble(j)*dqy
                i2=i2+1
                   h_q=0
                     q=dsqrt(qx*qx+qy*qy)
                     if(q.ge.dqx) then
                       cosq2=0.d00  ! for argon
                       sinq2=0.d00
                       natoms2=0
                       do kk=1,natoms
                          indx=3*(kk-1)
                          qr=x(indx+1)*qx+x(indx+2)*qy!+x(indx+3)*qz
                          cosq2=cosq2+dcos(qr)
                          sinq2=sinq2+dsin(qr)
                          natoms2=natoms2+1
                       enddo
                       sq(i1,i2)=(cosq2**2.+sinq2**2.)/dble(natoms)
                     endif
              enddo
            enddo
            print*, 'nx,ny',nx,ny
            print*, 'i1,i2',i1,i2
            ! take average
            do i=1,i1     
              do j=1,i2     
                sq_block(i,j)=sq_block(i,j)+sq(i,j)
              enddo
            enddo
            print*, 'll,nsteps_block',ll,nsteps_block
            if(mod(ll,nsteps_block).eq.0) then
 
                print*, 'calculating block averages'
                do i=1,i1
                   qx=dqx*dble(-nx+1+i)+dqx/2.d00
                   do j=1,i2
                      qy=dqy*dble(-ny+1+i)+dqy/2.d00
                      sq_block(i,j)=sq_block(i,j)/dble(nsteps_block)
                      write(89,*) qx,qy,sq_block(i,j)
                      sq_av(i,j)=sq_av(i,j)+sq_block(i,j)
                      sq2_av(i,j)=sq2_av(i,j)+(sq_block(i,j)**2.)
                      sq_block(i,j)=0.d0
                    enddo
                    write(89,*) 
                enddo
             endif
        enddo
        print*, 'here'

        open(83,file='Sq-av.dat')
        i1=0
        do i=-nx,nx    
            qx=dqx*dble(i)
            i1=i1+1
            i2=0
            do j=-ny,ny   
               qy=dqy*dble(j)
               i2=i2+1
               sq_av(i1,i2)=sq_av(i1,i2)/dble(nblock)
               sq2_av(i1,i2)=sq2_av(i1,i2)/dble(nblock)
               sigma2=(sq2_av(i1,i2)-sq_av(i1,i2)*sq_av(i1,i2))/dble(nblock)
               write(83,'(3f16.8)') qx,qy,sq_av(i1,i2)
            enddo
            write(83,*)
        enddo
	stop
	end
c**********************************************************************
c
        subroutine rdf(natoms1,natoms2,label,x,box,deltagri,
     +               delta_anglei,dmaxgr,rcut,gr,g_theta)
        implicit none
        integer natoms,i,j,k,indx,ibin,bmax,natoms1,natoms2
        integer jp,kp,indxj,indxk,neigh_max
        integer nneigh(10000),neigh(10000,5000)
        double precision x(50000),box(3),gr(10000,7)
        double precision xi,yi,zi,dx,dy,dz,dr,r12,rcut2
        double precision deltagri,dmaxgr,g_theta(180)
        double precision sigma,eps,sigma_wf,eps_wf,rcut,vl0
        double precision rij,rik,ctheta,theta,delta_anglei
        double precision dxj,dyj,dzj,dxk,dyk,dzk,x_ox(30000)
        common /model/ sigma,eps,sigma_wf,eps_wf,rcut2,vl0
        character*2 label(20000)

        neigh_max=20
        do i=1,natoms2
             nneigh(i)=0
        enddo
        do j=1,neigh_max
             do i=1,natoms2
                   neigh(i,j)=0
             enddo
        enddo
        
        bmax=0
        natoms=natoms1+natoms2
        do i=1,natoms-1
            indx=3*(i-1)
            xi=x(indx+1)
            yi=x(indx+2)
            zi=x(indx+3)
            do j=i+1,natoms
                indxj=3*(j-1)
                dx=xi-x(indxj+1)
                dy=yi-x(indxj+2)
                dz=zi-x(indxj+3)
c
                dx=dx-box(1)*dnint(dx/box(1))
                dy=dy-box(2)*dnint(dy/box(2))
                dz=dz-box(3)*dnint(dz/box(3))
c
               ! if(label(i).eq.'N2'.and.label(j).eq.'N2') then 
               !   ibin=idint(dx*deltagri)+1
               !   gr(ibin,5)=gr(ibin,5)+2.d00
               !   ibin=idint(dy*deltagri)+1
               !   gr(ibin,6)=gr(ibin,6)+2.d00
               !   ibin=idint(dz*deltagri)+1
               !   gr(ibin,7)=gr(ibin,7)+2.d00
               ! endif
c
                dr=dsqrt(dx*dx+dy*dy+dz*dz)
                if(dr.lt.dmaxgr) then
                  ibin=idint(dr*deltagri)+1
                  if(ibin.gt.bmax) bmax=ibin
                  if(label(i).eq.'LA'.and.label(j).eq.'LA') then 
                       gr(ibin,1)=gr(ibin,1)+2.d00
                       gr(ibin,4)=gr(ibin,4)+2.d00
                       if(dr.lt.rcut) then
                           nneigh(i-natoms1)=nneigh(i-natoms1)+1
                           nneigh(j-natoms1)=nneigh(j-natoms1)+1
                           neigh(i-natoms1,nneigh(i-natoms1))=j
                           neigh(j-natoms1,nneigh(j-natoms1))=i
                      endif
                  elseif(label(i).eq.' O'.and.label(j).eq.'LA') then
                       gr(ibin,2)=gr(ibin,2)+1.d00
                       gr(ibin,3)=gr(ibin,3)+1.d00
                       gr(ibin,4)=gr(ibin,4)+2.d00
                  !elseif(label(i).eq.' O'.and.label(j).eq.' O') then
                  !     gr(ibin,4)=gr(ibin,4)+2.d00
                  endif
                endif
c
c
            enddo
        enddo
c
        ! calculate the bond-angle distribution fucntion
        !do i=1,natoms2
        !     print*, 'nne',i,nneigh(i)
        !enddo
        !stop
        do i=natoms1+1,natoms
            indx=3*(i-1)
            xi=x(indx+1)
            yi=x(indx+2)
            zi=x(indx+3)
            do j=1,nneigh(i-natoms1)
               jp=neigh(i-natoms1,j)
               indxj=3*(jp-1)
               dxj=x(indxj+1)-xi
               dyj=x(indxj+2)-yi
               dzj=x(indxj+3)-zi
               dxj=dxj-box(1)*dnint(dxj/box(1))
               dyj=dyj-box(2)*dnint(dyj/box(2))
               dzj=dzj-box(3)*dnint(dzj/box(3))
               rij=dsqrt(dxj*dxj+dyj*dyj+dzj*dzj)
               do k=1,nneigh(i-natoms1)
                   if(k.ne.j) then
                     kp=neigh(i-natoms1,k)
                     indxk=3*(kp-1)
                     dxk=x(indxk+1)-xi
                     dyk=x(indxk+2)-yi
                     dzk=x(indxk+3)-zi
                     dxk=dxk-box(1)*dnint(dxk/box(1))
                     dyk=dyk-box(2)*dnint(dyk/box(2))
                     dzk=dzk-box(3)*dnint(dzk/box(3))
                     rik=dsqrt(dxk*dxk+dyk*dyk+dzk*dzk)
                     ctheta=(dxj*dxk+dyj*dyk+dzj*dzk)/rij/rik
                     if(ctheta.gt.1.d00) ctheta=1.d0
                     if(ctheta.lt.-1.d00) ctheta=-1.d0
                     theta=acos(ctheta)
                     ibin=idint(theta*delta_anglei)+1
                     g_theta(ibin)=g_theta(ibin)+1
                   endif
               enddo
            enddo
        enddo
c
        return
        end
c
c**********************************************************************
c

