        program OP_MC
        implicit none
        integer i,j,k,l,m,indxj,nmax,nnmax,nmol,nconf,ios,ngrmax,ngrmax2,kl,nnmax2
        parameter(nmax=120000,nnmax=50,nnmax2=500,ngrmax=50000,ngrmax2=401)
        integer nneigh(nmax),neigh(nmax,nnmax),ibin,nsize,nc3
        integer nneigh2(nmax),neigh2(nmax,nnmax2)
        integer nneigh_av(nmax),n_vH,num_points,npoints(6)
        integer nlarge,Nsizemax,nshort,nsol,nClust_min
        integer nClust,num_clusters,nthres,timestep1,nClust_NN
        integer inClust(120000),clustHead(120000),clustSize(120000)
        integer inClustNN(120000),clustSize_NN(120000)
        integer clustNext(120000),clustSize_order(120000)
        integer nfcc, nhcp, indx,natoms,idum,icore(120000),nc1
        integer indx_x,indx_y,indx_z,ngr,npar,iclust
        integer nmol_h,nmol_s,nmol_Z,nmol_un,angt(nnmax)
        integer n60,n90,nun,nmol_c1,nmol_few_nei,nmol_many_nei
        integer itype(nmax),icenter,itypei,nt_vanHove,ncounter
        ! parameters for calculation of vpot
        integer Npart_types_max,Nrb_sites_max, Ntor_max, Ntmax
        parameter(Npart_types_max=5, Nrb_sites_max = 10, Ntor_max = 3, Ntmax=410)
        integer Npart_types,bool_tor,nc(Npart_types_max)
        integer Nrb_sites(Npart_types_max),num_neigh
        integer iref_vec(Npart_types_max,Nrb_sites_max)
        integer Ntor_angle(Npart_types_max,Nrb_sites_max)
        integer n_c1,n_c2,n_c3,n_c4,n_c5,n_tot,ncont,imin
        integer nc1_write,n_bound,n_isolated,jaux,npc
        integer ncount_par(nmax),isol(nmax),npcNN,npcc,npcs

        ! end of: parameters for calculation of vpot
        integer nbin_x,nbin_y,nbin_z,nbin,ntime,nveces
        double precision gr_av(ngrmax),dmin
        double precision xn_c1_av,xn_c2_av,xn_c3_av,xn_c4_av,xn_c5_av,n_tot_av
        double precision ax(3,3),ay(3,3),anglex,angley,anglez,xp,yp,zp
        double precision x_min,x_max,y_min,y_max,z_min,z_max,xpp,ypp,zpp
        double precision x_medium,y_medium,z_medium,xo,yo,zo,xside
        double precision x(nmax),y(nmax),z(nmax),x0,y0,z0,az(3,3)
        double precision xav(nmax),yav(nmax),zav(nmax),xn_av
        double precision xinit(nmax),yinit(nmax),zinit(nmax)
        double precision x_store(nmax,Ntmax),y_store(nmax,Ntmax),
     +                     z_store(nmax,NTmax),delta_gr_vH,xcheck_norm(6)
        double precision rmsd_i(nmax),rmsd_type(Npart_types_max)
        double precision patchx(nmax,Nrb_sites_max)
        double precision patchy(nmax,Nrb_sites_max)
        double precision patchz(nmax,Nrb_sites_max)
        double precision xn(nmax),yn(nmax),zn(nmax),xvec_prom
        double precision vn(3),v1(3),v2(3),v3(3),pos_n(3)
        double precision xi,yi,zi,dx,dy,dz,box(3),dr,rmax
        double precision dr2_n(nmax,nnmax2),ener_part(nmax)
        double precision dx_n_proj(nnmax,3),dr_n_proj(nnmax)

        double precision pl(0:12,-12:12)
        double precision hs_re(0:12,-12:12),hs_im(0:12,-12:12)
        double precision q_av(nmax,0:12)
        double precision q(nmax,0:12),factor(0:12)
        double precision q_re(nmax,0:12,-12:12),q_im(nmax,0:12,-12:12)
        double precision q_re_av(0:12,-12:12),q_im_av(0:12,-12:12)
        double precision dij,di,dj,dot(0:12)

        double precision vy(3),cthetap,norm,norm1,normj
        double precision rcut,proj_norm(nnmax),theta_harm
        double precision theta(nnmax),ctheta,stheta,phi,xpi,xpi4
        double precision thetar(nnmax),G_vanHove(ngrmax2,ngrmax2,6)
        double precision G_vanHove_aux(ngrmax2,ngrmax2,6)
        double precision G_vanHove_scalar(ngrmax2)
        double precision G_vanHove_scalar_aux(ngrmax2)
        double precision timestep2,dum,delta_gr2
        double precision deltar,delta_gr,gr(ngrmax),q6_limit,q4_limit
        double precision gx(ngrmax),gy(ngrmax),gz(ngrmax),rhoR(ngrmax)
        double precision enerR(ngrmax),enerR_av(ngrmax)
        double precision vol,rho,const,const2,rbajo,ralto,dist,cideal
        double precision h(3,3),h_inv(3,3),dxp,dyp,dzp,rhoR_av(ngrmax)
        double precision x1,y1,z1,x2,y2,z2,deth,tinit,tend,dist2,dist3
        double precision alpha,beta,gama,xcdm,ycdm,zcdm,pnorm
        double precision px1,py1,pz1,en,en_tot,rmsd,rmsd_av
        double precision px,py,pz,r_esfera_central
        double precision histo_vec(Npart_types_max,0:nnmax)
        double precision xix_sum,xiy_sum,xiz_sum,zetax_sum,zetay_sum,zetaz_sum
        double precision thetaxi,thetayi,thetazi,thetabarx,thetabary,thetabarz
        double precision xcom,ycom,zcom
        ! parameters for calculation of vpot
        double precision sigma_tor_jon, rangep, VL0, xop
        double precision vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
        double precision sigma_jon(Npart_types_max,Nrb_sites_max)
        double precision tor_angle(Npart_types_max,Nrb_sites_max,Ntor_max)
        ! end of: parameters for calculation of vpot
        double precision histo_Clust(nmax),histo_Clust_av(nmax)
        character*30 filename,filename2,fname
        character*4 el
        character*4 el2
        character*3 elem(120000)
        character*1 env(120000)
        logical iprint,success, overlap
        common /constants/ xpi,xpi4
        common /limits/ q6_limit
        common /clust1/ inClust,clustHead,clustSize,clustNext
        common /clust2/ clustSize_order

        ! parameters for calculation of vpot
        common /scalar_nsites/ Npart_types, bool_tor
        common /vpot_par/ VL0, xop, sigma_tor_jon, rangeP
        common /array_vpot_par/ vpot_matrix, sigma_jon
        common /nsites/ Nrb_sites, iref_vec, Ntor_angle
       common /array_tor/ tor_angle


        ! coordinates
        common /box1/ h
        common /coords/ x, y, z
        common /coordspacth/ patchx, patchy, patchz
        common /coordstype/ itype

        xpi=acos(-1.d0)
        xpi4=4.d0*xpi
        const2=(4./3.)*acos(-1.d0)
        print*, 'Enter the rcut, esfera central'
        read(*,*) 
        read(*,*) rcut,r_esfera_central
        print*, 'Enter the minimum size for a cluster'
        read(*,*) 
        read(*,*) nClust_min
        print*, 'Print configurations (T or F)?'
        read(*,*) 
        read(*,*) iprint
        print*, 'Print initial and final time of saved configurations'
        read(*,*) 
        read(*,*) tinit,tend
        print*, 'time for evalutation of vanHove function (number of frames)'
        read(*,*) 
        read(*,*) nt_vanHove

        alpha=90.
        beta=90.
        gama=90.
        ! rotate y
        angley=0.0! -3.5*xpi/180.d0
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
        anglex=0.0! -10.*xpi/180.d0
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
        anglez=0.0! 9.0*xpi/180.d0
        az(1,1)=dcos(anglez)
        az(1,2)=-dsin(anglez)
        az(1,3)=0.d0         
        az(2,1)=dsin(anglez)
        az(2,2)=dcos(anglez)
        az(2,3)=0.d0           
        az(3,1)=0.d0          
        az(3,2)=0.d0           
        az(3,3)=1.d0          

        print*, 'before input'
        Call Read_Input_Data !(Npart_types_max,Nrb_sites_max, Ntor_max)
        print*, 'after input'

        do i=1,nmax
           nneigh_av(i)=0
        enddo

        open(27,file='conf-core.xyz')
        open(20,file='movie.xyz')
        open(28,file='conf-core-patch.xyz')
        open(29,file='movie-centred.xyz')
        open(36,file='Cluster-size.dat')
        open(35,file='Local-environment.dat')
        open(47,file='conc-vs-time.dat')
        open(57,file='Ener-vs-time.dat')
        ios=0
        nconf=0
        do i=1,nmax
           rmsd_i(i)=0.d0
           xav(i)=0.d0
           yav(i)=0.d0
           zav(i)=0.d0
        enddo
        do i=1,Npart_types_max
           rmsd_type(i)=0.d0
          do j=1,nnmax
              histo_vec(i,j)=0.d0
          enddo
        enddo

        do i=1,ngrmax
            gr_av(i)=0.d0
            rhoR(i)=0.d0
            rhoR_av(i)=0.d0
            enerR_av(i)=0.d0
        enddo
        do i=1,ngrmax2
            G_vanHove_scalar(i)=0.d0
        enddo
        do i=1,ngrmax2
           do j=1,ngrmax2
              do  k=1,6 
                G_vanHove(i,j,k)=0.d0
              enddo
           enddo
        enddo

        xn_c1_av=0
        xn_c2_av=0
        xn_c3_av=0
        xn_c4_av=0
        xn_c5_av=0
        n_tot_av=0

        call init_clusters(nClust)
        rmsd_av=0.d0
        deltar=1.0
        delta_gr=0.01
        delta_gr2=2.00
        n_vH=101
        delta_gr_vH=4.0/dble(n_vH)
        do i=1,nmax   
            histo_Clust_av(i)=0.d0
        enddo
c
        do while(ios.eq.0) 
c
           print*, 'starting',natoms
           read(20,*,IOSTAT=ios) natoms              
           print*, 'natoms=',natoms,ios
           if(ios.eq.0) then
           read(20,*,IOSTAT=ios) box(1),idum
           indx=0
           do i=1,natoms
            read(20,*,IOSTAT=ios) el,x1,y1,z1

            if(el.eq.'C1') then
               indx=indx+1
               x(indx)=x1
               y(indx)=y1
               z(indx)=z1
               elem(indx)=el
               itype(indx)=1
               j=0
            elseif(el.eq.'C2') then
               indx=indx+1
               x(indx)=x1
               y(indx)=y1
               z(indx)=z1
               elem(indx)=el
               itype(indx)=2
               j=0
            elseif(el.eq.'C3') then
               indx=indx+1
               x(indx)=x1
               y(indx)=y1
               z(indx)=z1
               elem(indx)=el
               itype(indx)=3
               j=0
            elseif(el.eq.'C4') then
               indx=indx+1
               x(indx)=x1
               y(indx)=y1
               z(indx)=z1
               elem(indx)=el
               itype(indx)=4
               j=0
            elseif(el.eq.'C5') then
               indx=indx+1
               x(indx)=x1
               y(indx)=y1
               z(indx)=z1
               elem(indx)=el
               itype(indx)=5
               j=0
            else
               j=j+1
               patchx(indx,j)=(x1-x(indx))/0.35d0
               patchy(indx,j)=(y1-y(indx))/0.35d0
               patchz(indx,j)=(z1-z(indx))/0.35d0
               pnorm=dsqrt(patchx(indx,j)**2.+patchy(indx,j)**2.+patchz(indx,j)**2.)
               if(dabs(pnorm-1.d0).gt.1.d-03) then
                 print*, 'norm deviates too much from unity',pnorm
                 stop
               endif
               patchx(indx,j)=patchx(indx,j)/pnorm
               patchy(indx,j)=patchy(indx,j)/pnorm
               patchz(indx,j)=patchz(indx,j)/pnorm
            endif
           enddo
           nmol=indx
           h(1,1)=box(1)
           h(2,2)=box(1)
           h(3,3)=box(1)
           h(1,2)=0.d0
           h(1,3)=0.d0
           h(2,1)=0.d0
           h(2,3)=0.d0
           h(3,1)=0.d0
           h(3,2)=0.d0
           print*, 'nmol=',nmol
           natoms=nmol

           nconf=nconf+1

           print*, 'ios',ios
           print*, 'nconf',nconf
           print*, 'box=',h(1,1),h(2,2),h(3,3),h(2,1),
     +   h(3,1),h(1,2),h(3,2),h(1,3),h(2,3)
            print*, '*************',timestep1
c
           call calculo_matriz_inversa(h,h_inv)
           call cambio_unidades_caja(nmol,nmax,x,y,z,h_inv)


           do i=1,nmax
                nneigh(i)=0
                nneigh2(i)=0
                ener_part(i)=0.d0
                do j=1,nnmax
                   neigh(i,j)=0
                   neigh2(i,j)=0
                enddo
           enddo
           print*, 'finding neighbours'

           en_tot=0.d0
           do i=1,nmol
               xi=x(i)
               yi=y(i)
               zi=z(i)
               do j=1,nmol
                  if( i.ne. j) then
                    dx=x(j)-xi
                    dy=y(j)-yi
                    dz=z(j)-zi
                    dxp=dx-dnint(dx)
                    dyp=dy-dnint(dy)
                    dzp=dz-dnint(dz)
                    dx=h(1,1)*dxp+h(1,2)*dyp+h(1,3)*dzp
                    dy=h(2,1)*dxp+h(2,2)*dyp+h(2,3)*dzp
                    dz=h(3,1)*dxp+h(3,2)*dyp+h(3,3)*dzp
                    dr=dsqrt(dx*dx+dy*dy+dz*dz)
                    if(dr.lt.2.5d0) then
                         call ener(i,j,en,overlap)
                         en_tot=en_tot+en
                         ener_part(i)=ener_part(i)+en
                         if(dr.lt.rcut.and.en.lt.-0.30d0) then
                           nneigh(i)=nneigh(i)+1
                           neigh(i,nneigh(i))=j
                         endif
                    endif
                    if(dr.lt.2.00d0) then
                       nneigh2(i)=nneigh2(i)+1
                       neigh2(i,nneigh2(i))=j
                       dr2_n(i,nneigh2(i))=dr
                    endif
                  endif
               enddo
           enddo
           print*, 'neigbours found'
           print*, 'Total energy=',en_tot,en_tot/dble(nmol)
           write(57,*) nconf,en_tot/dble(nmol)/2.d0
              
           deth=h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))-h(2,1)*(h(1,2)*h(3,3)
     *    -h(1,3)*h(3,2))+h(3,1)*(h(1,2)*h(2,3)-h(1,3)*h(2,2))
           vol=deth
           rho=dble(nmol)/vol
           const=4.00*xpi*dble(nmol)/3.d00/vol

!#########################################################################3#######
           print*, 'before cluster analysis'
           !do i=1,natoms
           !    histo_Clust(i)=0.d0
           !enddo

           call cluster_analysis(neigh,nneigh,nmax,nnmax,nmol,nclust_min,nClust)
           print*, 'after cluster analysis'
           nClust_NN=nClust
           do i=1,nmol
             inClustNN(i)=inClust(i)
           enddo
           do i=1,nClust_NN
                clustSize_NN(i)=clustSize_order(i)
           enddo

           write(*,'(I10,3x,6I6)') nconf,(clustSize_order(k),k=1,5),Nsol
           ! distribucion de clusters
           !n_bound=0
           !do i=1,nClust
           !   nSize=clustSize_order(i)
           !   n_bound=n_bound+nSize
           !   histo_Clust(nSize)=histo_Clust(nSize)+nSize
           !enddo
           !n_isolated=natoms-n_bound
           !histo_Clust(1)=n_isolated

           !do i=1,natoms
           !   histo_Clust(i)=histo_Clust(i)/dble(natoms)
           !   histo_Clust_av(i)=histo_Clust_av(i)+histo_Clust(i)
           !enddo

           ! method to try to center the largest cluster
           print*, 'centering, nconf',nconf

           !if(nconf.eq.1) then
               xix_sum=0.d0
               xiy_sum=0.d0
               xiz_sum=0.d0
               zetax_sum=0.d0
               zetay_sum=0.d0
               zetaz_sum=0.d0
               ncont=0
               do i=1,nmol
                  if(inClust(i).eq.1) then
                    thetaxi= x(i)*2.*xpi
                    thetayi= y(i)*2.*xpi
                    thetazi= z(i)*2.*xpi
                    xix_sum=xix_sum+dcos(thetaxi)
                    xiy_sum=xiy_sum+dcos(thetayi)
                    xiz_sum=xiz_sum+dcos(thetazi)
                    zetax_sum=zetax_sum+dsin(thetaxi)
                    zetay_sum=zetay_sum+dsin(thetayi)
                    zetaz_sum=zetaz_sum+dsin(thetazi)
                    ncont=ncont+1
                 endif
              enddo
              print*, 'ncont=',ncont
              xix_sum=xix_sum/dble(ncont)
              xiy_sum=xiy_sum/dble(ncont)
              xiz_sum=xiz_sum/dble(ncont)
              zetax_sum=zetax_sum/dble(ncont)
              zetay_sum=zetay_sum/dble(ncont)
              zetaz_sum=zetaz_sum/dble(ncont)
              thetabarx=datan2(-zetax_sum,-xix_sum)+xpi
              thetabary=datan2(-zetay_sum,-xiy_sum)+xpi
              thetabarz=datan2(-zetaz_sum,-xiz_sum)+xpi
              xcom=thetabarx/2./xpi
              ycom=thetabary/2./xpi
              zcom=thetabarz/2./xpi
              print*, 'com',xcom,ycom,zcom

               do i=1,nmol
                 if(inClust(i).eq.1) then
                     xi=x(i)-xcom
                     yi=y(i)-ycom
                     zi=z(i)-zcom
                     xi=xi-dnint(xi)
                     yi=yi-dnint(yi)
                     zi=zi-dnint(zi)
                     dx=h(1,1)*xi+h(1,2)*yi+h(1,3)*zi
                     dy=h(2,1)*xi+h(2,2)*yi+h(2,3)*zi
                     dz=h(3,1)*xi+h(3,2)*yi+h(3,3)*zi
                     dr=dsqrt(dx*dx+dy*dy+dz*dz)
                     if(dr.gt.rmax) rmax=dr
                  endif
              enddo

              print*, 'radius of the cluster',rmax
              nc1_write=0
              dmin=1000.
              do i=1,nmol
                  icore(i)=0
                  if(inClust(i).eq.1) then
                     xi=x(i)-xcom
                     yi=y(i)-ycom
                     zi=z(i)-zcom
                     xi=xi-dnint(xi)
                     yi=yi-dnint(yi)
                     zi=zi-dnint(zi)
                     dx=h(1,1)*xi+h(1,2)*yi+h(1,3)*zi
                     dy=h(2,1)*xi+h(2,2)*yi+h(2,3)*zi
                     dz=h(3,1)*xi+h(3,2)*yi+h(3,3)*zi
                     dr=dsqrt(dx*dx+dy*dy+dz*dz)
                     if(dr.lt.dmin) then
                          dmin=dr
                          imin=i
                     endif
                     if(dr.lt.0.65*rmax) then  ! 0.52=0.80*0.65
                        icore(i)=1
                        nc1_write=nc1_write+1
                     endif
                  endif
              enddo
              print*, 'nc1=',nc1_write
              print*, 'particle closer to com',imin
          !endif
          !imin=9953

          xcom=x(imin)
          ycom=y(imin)
          zcom=z(imin)
          print*, 'com',xcom,ycom,zcom
          do i=1,nmol
                x(i)=x(i)-xcom
                y(i)=y(i)-ycom
                z(i)=z(i)-zcom
                dx=x(i)-dnint(x(i))
                dy=y(i)-dnint(y(i))
                dz=z(i)-dnint(z(i))
                x1=h(1,1)*dx+h(1,2)*dy+h(1,3)*dz
                y1=h(2,1)*dx+h(2,2)*dy+h(2,3)*dz
                z1=h(3,1)*dx+h(3,2)*dy+h(3,3)*dz
                xp=x1*ay(1,1)+y1*ay(1,2)+z1*ay(1,3)  
                yp=x1*ay(2,1)+y1*ay(2,2)+z1*ay(2,3)  
                zp=x1*ay(3,1)+y1*ay(3,2)+z1*ay(3,3)  
                xpp=xp*ax(1,1)+yp*ax(1,2)+zp*ax(1,3)
                ypp=xp*ax(2,1)+yp*ax(2,2)+zp*ax(2,3)
                zpp=xp*ax(3,1)+yp*ax(3,2)+zp*ax(3,3)
                x(i)=xpp*az(1,1)+ypp*az(1,2)+zpp*az(1,3)
                y(i)=xpp*az(2,1)+ypp*az(2,2)+zpp*az(2,3)
                z(i)=xpp*az(3,1)+ypp*az(3,2)+zpp*az(3,3)
                itypei=itype(i)
                do j=1,nrb_sites(itypei)
                      px=x1+patchx(i,j)*0.35d0
                      py=y1+patchy(i,j)*0.35d0
                      pz=z1+patchz(i,j)*0.35d0
                      xp=px*ay(1,1)+py*ay(1,2)+pz*ay(1,3)  
                      yp=px*ay(2,1)+py*ay(2,2)+pz*ay(2,3)  
                      zp=px*ay(3,1)+py*ay(3,2)+pz*ay(3,3)  
                      xpp=xp*ax(1,1)+yp*ax(1,2)+zp*ax(1,3)
                      ypp=xp*ax(2,1)+yp*ax(2,2)+zp*ax(2,3)
                      zpp=xp*ax(3,1)+yp*ax(3,2)+zp*ax(3,3)
                      xp=xpp*az(1,1)+ypp*az(1,2)+zpp*az(1,3)
                      yp=xpp*az(2,1)+ypp*az(2,2)+zpp*az(2,3)
                      zp=xpp*az(3,1)+ypp*az(3,2)+zpp*az(3,3)
                      patchx(i,j)=(xp-x(i))/0.35d0
                      patchy(i,j)=(yp-y(i))/0.35d0
                      patchz(i,j)=(zp-z(i))/0.35d0
                      pnorm=dsqrt(patchx(i,j)**2.+patchy(i,j)**2.+patchz(i,j)**2.)
                      if(dabs(pnorm-1.d0).gt.1.d-03) then
                        print*, 'norm deviates too much from unity',pnorm
                        stop
                      endif
                      patchx(i,j)=patchx(i,j)/pnorm
                      patchy(i,j)=patchy(i,j)/pnorm
                      patchz(i,j)=patchz(i,j)/pnorm
                enddo
          enddo

               
           !una vez que esta centrado, calculo Rho(r)

           ngr=int(0.5*box(1)/delta_gr2)
           print*, 'delta_gr,ngr=',delta_gr2,ngr
           do i=1,nmol
               dist=dsqrt(x(i)**2.+y(i)**2.+z(i)**2.)
               ibin=int(dist/delta_gr2)+1
               rhoR(ibin)=rhoR(ibin)+1.d0
               enerR(ibin)=enerR(ibin)+ener_part(i)
           enddo
           if(nconf.lt.10) then
                 write(fname,'("rhoR",I1,".dat")')nconf    
           elseif(nconf.lt.100) then
                 write(fname,'("rhoR",I2,".dat")')nconf    
           elseif(nconf.lt.1000) then
                 write(fname,'("rhoR",I3,".dat")')nconf    
           endif
           open(48,file=fname)
           ngr=int(0.5*box(1)/delta_gr2)
           do i=1,ngr
             rbajo=dfloat(i-1)*delta_gr2
             ralto=rbajo+delta_gr2
             dist=rbajo+delta_gr2/2.d00
             vol=const2*(ralto**3-rbajo**3.d0)
             enerR(i)=enerR(i)/rhoR(i)/2.d0
             rhoR(i)=rhoR(i)/vol
             enerR_av(i)=enerR_av(i)+enerR(i)
             rhoR_av(i)=rhoR_av(i)+rhoR(i)
             if(dist.gt.0.50.and.rhoR_av(i).gt.1.d-07) write(48,*) dist,rhoR(i),enerR(i)
             rhoR(i)=0.d0
             enerR(i)=0.d0
           enddo
           close(48)
!#########################################################################

           ! calculate the order parameters
           do i=1,nmol
             isol(i)=0
             do l=0,12
             do m=-12,12
                  q_re(i,l,m)=0.d0
                  q_im(i,l,m)=0.d0
                enddo
             enddo
           enddo

           call cambio_unidades_caja(nmol,nmax,x,y,z,h_inv)

           do i=1,nmol            
               xi=x(i)
               yi=y(i)
               zi=z(i)
               do j=1,nneigh2(i)
                  indxj=neigh2(i,j)
                  dx=x(indxj)-xi
                  dy=y(indxj)-yi
                  dz=z(indxj)-zi
                  dxp=dx-dnint(dx)
                  dyp=dy-dnint(dy)
                  dzp=dz-dnint(dz)
                  dx=h(1,1)*dxp+h(1,2)*dyp+h(1,3)*dzp
                  dy=h(2,1)*dxp+h(2,2)*dyp+h(2,3)*dzp
                  dz=h(3,1)*dxp+h(3,2)*dyp+h(3,3)*dzp
                  dr=dsqrt(dx*dx+dy*dy+dz*dz)
                  dx=dx/dr
                  dy=dy/dr
                  dz=dz/dr
                  call spherical(dx,dy,dz,theta_harm,phi)
                  ctheta=cos(theta_harm)
                  
                  call legendre(pl,ctheta)
                  call SphericalHarmonics(pl,phi,hs_re,hs_im)
                  do l=0,12
                     do m=-l,l
                        q_re(i,l,m)=q_re(i,l,m)+hs_re(l,m)
                        q_im(i,l,m)=q_im(i,l,m)+hs_im(l,m)
                     enddo
                  enddo

                  call spherical(-dx,-dy,-dz,theta_harm,phi)
                  ctheta=cos(theta_harm)
                  call legendre(pl,ctheta)
                  call SphericalHarmonics(pl,phi,hs_re,hs_im)
                  do l=0,12
                     do m=-l,l
                        q_re(indxj,l,m)=q_re(indxj,l,m)+hs_re(l,m)
                        q_im(indxj,l,m)=q_im(indxj,l,m)+hs_im(l,m)
                      enddo
                  enddo

               enddo ! j=1,neigh(i)
           enddo ! i=1,nmol

           do l=0,12
             factor(l)=xpi4/(2.d0*dble(l)+1.d0)
           enddo

           do i=1,nmol
             do l=0,12
               do m=-l,l
               if (nneigh2(i).gt.2) then
                 q_re(i,l,m)=q_re(i,l,m)/dble(nneigh2(i))
                 q_im(i,l,m)=q_im(i,l,m)/dble(nneigh2(i))
               else
                 q_re(i,l,m)=0.d0
                 q_im(i,l,m)=0.d0
               endif
               enddo
             enddo
           enddo

           do i=1,nmol
               do l=0,12
                 q(i,l)=0.d0
                 do m=-l,l
                  q(i,l)=q(i,l)+q_re(i,l,m)**2.d0+q_im(i,l,m)**2.d0
                 enddo
                 q(i,l)=dsqrt(factor(l)*q(i,l))
               enddo
           enddo

           ! qv Lechner-Dellago

           open(73,file='OP-Dellago.dat')
           do i=1,nmol
              xi=x(i)
              yi=y(i)
              zi=z(i)
              dr=box(1)*dsqrt(xi*xi+yi*yi+zi*zi)
              do l=0,12
                 do m=-l,l
                    q_re_av(l,m)=q_re(i,l,m)
                    q_im_av(l,m)=q_im(i,l,m)
                    do j=1,nneigh2(i)
                       indxj=neigh2(i,j)
                       q_re_av(l,m)=q_re_av(l,m)+q_re(indxj,l,m)
                       q_im_av(l,m)=q_im_av(l,m)+q_im(indxj,l,m)
                    enddo
                    q_re_av(l,m)=q_re_av(l,m)!/dble(nneigh(i)+1)
                    q_im_av(l,m)=q_im_av(l,m)!/dble(nneigh(i)+1)
                 enddo !m=-l,l

                 q_av(i,l)=0.d0
                 do m=-l,l
                     q_av(i,l)=q_av(i,l)+q_re_av(l,m)**2.d0+
     +              q_im_av(l,m)**2.d0
                 enddo
                 if(q_av(i,l).gt.1.d-08) then
                    q_av(i,l)=dsqrt(factor(l)*q_av(i,l))
                 endif
              enddo !l=0,12
              if(dr.lt.r_esfera_central) write(73,'(12f16.8)') (q_av(i,l),l=0,12)
           enddo
           close(73)

           ! paso a unidades reales
           do i=1,nmol
                dx=x(i)-dnint(x(i))
                dy=y(i)-dnint(y(i))
                dz=z(i)-dnint(z(i))
                x(i)=h(1,1)*dx+h(1,2)*dy+h(1,3)*dz
                y(i)=h(2,1)*dx+h(2,2)*dy+h(2,3)*dz
                z(i)=h(3,1)*dx+h(3,2)*dy+h(3,3)*dz
           enddo

           ! calculate the dot products
           open(73,file='OP.dat')
           do i=1,nmol
              xi=x(i)
              yi=y(i)
              zi=z(i)
              dr=dsqrt(xi*xi+yi*yi+zi*zi)
              ncount_par(i)=0
              if(dr.lt.r_esfera_central) then
              !write(23,'("C",3x,3f16.8)') xi,yi,zi
              if(nneigh2(i).gt.2) then
              do j=1,nneigh2(i)
                 indxj=neigh2(i,j)
                 do l=0,12
                    dij=0.d0
                    di=0.d0
                    dj=0.d0
                    do m=-l,l
                      dij=dij+q_re(i,l,m)*q_re(indxj,l,m)+
     +                   q_im(i,l,m)*q_im(indxj,l,m)
                      di=di+q_re(i,l,m)*q_re(i,l,m)+
     +                   q_im(i,l,m)*q_im(i,l,m)
                      dj=dj+q_re(indxj,l,m)*q_re(indxj,l,m)+
     +                   q_im(indxj,l,m)*q_im(indxj,l,m)
                    enddo !m=-l,l
                    if(di.gt.0.000000001.and.dj.gt.0.0000000001) then
                      dot(l)=dij/dsqrt(di)/dsqrt(dj)
                    endif
                 enddo
                 if(dot(12).gt.0.40) ncount_par(i)=ncount_par(i)+1
!                 write(73,'(13f16.4,3x,3i6)') 
!     +          (dot(l),l=0,12),nneigh2(i),nneigh(i)
              enddo !j=1,nneigh(i)
!              write(74,'(4i7)') i,ncount_par(i),nneigh2(i),nneigh(i)
              if(ncount_par(i).ge.5) isol(i)=1
              endif
              endif
           enddo ! i=1,nmol
           print*, 'q(0:12) dot products calculated'
           close(73)
           

!#########################################################################3#######
           print*, 'before cluster analysis'
           do i=1,natoms
               histo_Clust(i)=0.d0
           enddo

           call cluster_analysis_dotq(isol,neigh,nneigh,nmax,nnmax,nmol,nclust_min,nClust)
           print*, 'after cluster analysis'

           write(*,'(I10,3x,6I6)') nconf,(clustSize_order(k),k=1,5),Nsol
           ! distribucion de clusters
           n_bound=0
           do i=1,nClust
              nSize=clustSize_order(i)
              n_bound=n_bound+nSize
              histo_Clust(nSize)=histo_Clust(nSize)+nSize
           enddo
           n_isolated=natoms-n_bound
           histo_Clust(1)=n_isolated

           do i=1,natoms
              histo_Clust(i)=histo_Clust(i)/dble(natoms)
              histo_Clust_av(i)=histo_Clust_av(i)+histo_Clust(i)
           enddo

           ! van Hove
           print*, '888888888888888888888888888888888888888'
           print*, '888888888888888888888888888888888888888'
           print*, 'nconf',nconf,nt_vanHove
           print*, '888888888888888888888888888888888888888'
           print*, '888888888888888888888888888888888888888'
           if(nconf.le.nt_vanHove) then   ! just store the coordinates
                   print*, 'here',nmol
              do i=1,nmol
                 if(icore(i).eq.1) then
                     x_store(i,nconf)=x(i)
                     y_store(i,nconf)=y(i)
                     z_store(i,nconf)=z(i)
                 endif
              enddo
              print*, 'exiting loop'
           else  ! evaluate histo and update coordinates
              nveces=(nconf-1)/nt_vanHove
              ntime=nconf-nt_vanHove*nveces
              print*, '**************************'
              print*, '**************************'
              print*, 'finally evaluatimng van Hove'
              print*, 'nveces=',nveces
              print*, 'ntime=',ntime 
              num_points=0
              do k=1,6
                 npoints(k)=0
              enddo
              do i=1,ngrmax2
                 G_vanHove_scalar_aux(i)=0.d0
                 do j=1,ngrmax2
                    do k=1,6
                       G_vanHove_aux(i,j,k)=0.d0
                    enddo
                 enddo
              enddo
              do i=1,nmol

                 if(icore(i).eq.1) then
                     dx=x(i)-x_store(i,ntime)
                     dy=y(i)-y_store(i,ntime)
                     dz=z(i)-z_store(i,ntime)
                     if(dx.gt.-2.d0.and.dy.gt.-2.0.and.dz.gt.-2.0
     +               .and.dx.lt.2.d0.and.dy.lt.2.0.and.dz.lt.2.0) then
                       nbin_x=int((dx+2.0)/delta_gr_vH)+1
                       nbin_y=int((dy+2.0)/delta_gr_vH)+1
                       nbin_z=int((dz+2.0)/delta_gr_vH)+1
                       num_points=num_points+1
                       !print*, 'aqui',dx,dy,dz,itype(i)
                       npoints(itype(i))=npoints(itype(i))+1
                       npoints(6)=npoints(6)+1
                       G_vanHove_aux(nbin_x,nbin_y,itype(i))=G_vanHove_aux(nbin_x,nbin_y,itype(i))+1.d0
                       G_vanHove_aux(nbin_x,nbin_y,6)=G_vanHove_aux(nbin_x,nbin_y,6)+1.d0
                       dr=dsqrt(dx*dx+dy*dy+dz*dz)
                       nbin=int(dr/delta_gr_vH)+1
                       G_vanHove_scalar_aux(nbin)=G_vanHove_scalar_aux(nbin)+1.d0
                       x_store(i,ntime)=x(i)
                       y_store(i,ntime)=y(i)
                       z_store(i,ntime)=z(i)
                     else
                         !print*, 'i',i,dx,dy,dz
                     endif
                 endif
              enddo
              print*, '******************'
              print*, '******************'
              print*, 'num_points=',num_points
              print*, '******************'
              print*, '******************'
              do i=1,ngrmax2
                G_vanHove_scalar(i)=G_vanHove_scalar(i)+G_vanHove_scalar_aux(i)/dble(num_points)      
                do j=1, ngrmax2
                    do k=1,6
                      G_vanHove(i,j,k)=G_vanHove(i,j,k)+G_vanHove_aux(i,j,k)/dble(npoints(k))
                    enddo
                enddo
              enddo
           endif

c*************************************************************************************

           if(nconf.lt.10) then
                 write(filename,'("OP-conf",I1,".xyz")')nconf    
           elseif(nconf.lt.100) then
                 write(filename,'("OP-conf",I2,".xyz")')nconf    
           elseif(nconf.lt.1000) then
                 write(filename,'("OP-conf",I3,".xyz")')nconf    
           endif

           if(iprint) then
           if(nconf.gt.tinit.and.nconf.le.tend) then
           open(26,file=filename)
           npc=0
           npcc=0
           npcs=0
           npcNN=0
           do i=1,nmol
              if(inClust(i).eq.1.or.inClustNN(i).eq.1) then
                  npc=npc+1
                  if(inClust(i).eq.1.and.inClustNN(i).eq.1) then
                      npcs=npcs+1
                  else
                      npcc=npcc+1
                  endif
              endif
              if(inClustNN(i).eq.1) then
                  npcNN=npcNN+1
              endif
           enddo
           write(26,'(4i8)') npc, npcNN,npcs,npcc   ! npc same, npc conflict
           write(27,*) nc1_write
           write(28,*) nc1_write+nc1_write*nrb_sites(1)
           write(26,*) box(1)
           write(27,181) box(1),0.,0.,0.,box(1),0.,0.,0.,box(1),dble(nconf)
           write(28,181) box(1),0.,0.,0.,box(1),0.,0.,0.,box(1),dble(nconf)
           write(29,*) nmol
           write(29,*)
 181       format('Lattice="',9f5.2,'" Properties=species:S:1:pos:R:3 Time=',f7.0)
           nsol=0
           nc1=0
           nc3=0
           n_c1=0
           n_c2=0
           n_c3=0
           n_c4=0
           n_c5=0
           n_tot=0

           do i=1,ngrmax
            gr(i)=0.d0
           enddo
           do i=1,5
              nc(i)=0
           enddo

           do i=1,nmol
                  if(inClust(i).eq.1) then
                     nsol=nsol+1
                     if(elem(i).eq.'C1') then 
                        n_c1=n_c1+1
                     endif
                     n_tot=n_tot+1
                     dx=x(i)!-box(1)*dnint(x(i)/box(1))
                     dy=y(i)!-box(1)*dnint(y(i)/box(1))
                     dz=z(i)!-box(1)*dnint(z(i)/box(1))
                  endif
                  itypei=itype(i)
                  if(elem(i).eq.'C1') el='S'
                  if(elem(i).eq.'C2') el='T'
                  if(elem(i).eq.'C3') el='U'
                  if(elem(i).eq.'C4') el='V'
                  if(elem(i).eq.'C5') el='X'
                  if(icore(i).eq.1) then
                     dx=x(i)!-box(1)*dnint(x(i)/box(1))
                     dy=y(i)!-box(1)*dnint(y(i)/box(1))
                     dz=z(i)!-box(1)*dnint(z(i)/box(1))
                       num_neigh=nneigh(i)
                       histo_vec(itypei,num_neigh)=histo_vec(itypei,num_neigh)+1
                       nc(itypei)=nc(itypei)+1
                       write(27,'(a3,3x,5f16.8)') el,dx,dy,dz
                       itypei=itype(i)
                       write(28,'(a3,3x,5f16.8)') elem(i),dx,dy,dz
                     do kl=1,nrb_sites(itypei)
                      px=dx+patchx(i,kl)*0.35d0
                      py=dy+patchy(i,kl)*0.35d0
                      pz=dz+patchz(i,kl)*0.35d0
                      if(kl.le.9) then
                        write(el,'("P",i1,i1)') itypei,kl 
                      write(28,'(a3,3x,5f16.8)') el,px,py,pz
                      else
                        write(el2,'("P",i1,i2)') itypei,kl 
                      write(28,'(a4,3x,5f16.8)') el2,px,py,pz
                      endif
                     enddo
                  endif
                     dx=x(i)!-box(1)*dnint(x(i)/box(1))
                     dy=y(i)!-box(1)*dnint(y(i)/box(1))
                     dz=z(i)!-box(1)*dnint(z(i)/box(1))
                     iclust=0
                     if(inClustNN(i).gt.1) then
                        iclust=clustSize_NN(inClustNN(i))
                     endif
                     if(inClust(i).eq.1.and.inClustNN(i).eq.1) then
                           el='SS'
                     elseif(inClustNN(i).eq.1) then
                           el='I'
                     elseif(inClust(i).eq.1) then
                           el='S'
                           print*, 'aqui=',el,inClustNN(i),inClust(i),iclust
                     elseif(inClustNN(i).gt.1.and.iclust.lt.9) then
                           el='P'
                     elseif(inClustNN(i).gt.1.and.iclust.gt.9) then
                           el='G'
                     else 
                           el='F'
                     endif
                     write(29,'(a3,3x,5f16.8)') el,dx,dy,dz
                  if(inClust(i).eq.1.or.inClustNN(i).eq.1) then 
                     write(26,'(a3,3x,5f16.8)') el,dx,dy,dz
                     nc1=nc1+1
                     xav(nc1)=xav(nc1)+x(i)
                     yav(nc1)=yav(nc1)+y(i)
                     zav(nc1)=zav(nc1)+z(i)
                     dr=dsqrt(dx*dx+dy*dy+dz*dz)
                     if(dr.lt.17.d0) then
                        nc3=nc3+1
                        do j=1,nneigh2(i)  ! cambio aqui
                            dist=dr2_n(i,j)
                            if(dist.lt.5.d0) then
                              ibin=int(dist/delta_gr)+1
                              gr(ibin)=gr(ibin)+1.d0
                            endif
                        enddo
                     endif
                  if(elem(i).eq.'C1') el='S'
                  if(elem(i).eq.'C2') el='T'
                  if(elem(i).eq.'C3') el='U'
                  if(elem(i).eq.'C4') el='V'
                  if(elem(i).eq.'C5') el='X'
                  dx=x(i)!-box(1)*dnint(x(i)/box(1))
                  dy=y(i)!-box(1)*dnint(y(i)/box(1))
                  dz=z(i)!-box(1)*dnint(z(i)/box(1))
                  endif
           enddo
           print*, 'number of particles in core',(nc(i),i=1,5)

           rho=dble(nc3)/((4./3.)*xpi*(17.d0)**3.)
           do i=1,ngrmax
                gr_av(i)=gr_av(i)+gr(i)/rho/dble(nc3)
           enddo
           if(mod(nconf,100).eq.0) then
               print*, 'eoeoeoeoeoeoeo',nconf
               write(filename2,'("histo_vec1-",I3,".dat")')nconf      
               open(70,file=filename2)
               xn_av=0.d0
               do j=1,Nrb_sites_max
                 write(70,*) j,histo_vec(1,j)/dble(nconf)/dble(nc(1))
                 xn_av=xn_av+dble(j)*histo_vec(1,j)/dble(nconf)/dble(nc(1))
               enddo
               write(70,*) 10.,xn_av
               close(70)
           endif


           xn_c1_av=xn_c1_av+dble(n_c1)/dble(n_tot)
           xn_c2_av=xn_c2_av+dble(n_c2)/dble(n_tot)
           xn_c3_av=xn_c3_av+dble(n_c3)/dble(n_tot)
           xn_c4_av=xn_c4_av+dble(n_c4)/dble(n_tot)
           xn_c5_av=xn_c5_av+dble(n_c5)/dble(n_tot)
           write(47,'(7i10)') nconf,n_c1,n_c2,n_c3,n_c4,n_c5,n_tot
           write(*,*) 'nconf,n_c1,n_c2,n_c3,n_c4,n_c5,n_tot'
           write(*,'(7i10)') nconf,n_c1,n_c2,n_c3,n_c4,n_c5,n_tot
           write(38,'(6i6)') nconf,xn_c1_av/dble(nconf),xn_c2_av/dble(nconf),
     +            xn_c3_av/dble(nconf),xn_c4_av/dble(nconf),xn_c5_av/dble(nconf)

           write(37,'(i5,3x,6f10.4)') nconf,dble(n_c1)/dble(n_tot),
     +                 dble(n_c2)/dble(n_tot),dble(n_c3)/dble(n_tot),
     +                 dble(n_c4)/dble(n_tot),dble(n_c5)/dble(n_tot)
           close(26)
           endif
           endif

           write(36,'(I10,3x,7I7)') nconf,npcNN,(clustSize_order(k),k=1,5),Nsol
        endif
        enddo
 105    format(a2,3f10.4)


        do i=1,5
            write(filename2,'("histo_vec",I1,".dat")')i        
            open(70,file=filename2)
            xn_av=0.d0
            do j=1,Nrb_sites_max
                 write(70,*) j,histo_vec(i,j)/dble(nconf)/dble(nc(i))
                 xn_av=xn_av+dble(j)*histo_vec(i,j)/dble(nconf)/dble(nc(i))
            enddo
            write(70,*) 10.,xn_av

            close(70)
        enddo
        ntime=nconf-nt_vanHove
        print*, ''
        print*, 'ntime=',ntime
        print*, 'nconf=',nconf 
        print*, 'nt_vanHove=',nt_vanHove 
        print*, ''
        open(81,file='vanHove_scalar_autocor.dat')
        open(82,file='vanHove_autorcor-C1.dat')
        open(83,file='vanHove_autorcor-C2.dat')
        open(84,file='vanHove_autorcor-C3.dat')
        open(85,file='vanHove_autorcor-C4.dat')
        open(86,file='vanHove_autorcor-C5.dat')
        open(87,file='vanHove_autorcor-tot.dat')
        do k=1,6
           xcheck_norm(k)=0.d0
        enddo
        do i=1,n_vH    
             dist=dble(i-1)*delta_gr_vH+delta_gr_VH*0.50d0
             G_vanHove_scalar(i)=G_vanHove_scalar(i)/dble(ntime)
             write(81,*) dist,G_vanHove_scalar(i)
             dist=-2.d0+dble(i-1)*delta_gr_vH+delta_gr_VH*0.50d0
             do j=1,n_vH    
                dist2=-2.d0+dble(j-1)*delta_gr_vH+delta_gr_VH*0.50d0
                do k=1,6
                if(G_vanHove(i,j,k).gt.1.d-16) then
                  G_vanHove(i,j,k)=G_vanHove(i,j,k)/dble(ntime)
                  xcheck_norm(k)=xcheck_norm(k)+G_vanHove(i,j,k)
                endif
                write(81+k,'(4f16.8)') dist,dist2,G_vanHove(i,j,k)
                enddo
             enddo
             write(82,*) 
             write(83,*) 
             write(84,*) 
             write(85,*) 
             write(86,*) 
             write(87,*) 
        enddo
        print*, 'chequeo la norma de la funcion de autocor. de VanHove'
        do k=1,5
           print*, 'xcheck_norm=',xcheck_norm(k)
        enddo
        close(81)
        close(82)


        open(30,file='PDF.dat')
        open(32,file='rhoR_av.dat')
        ngr=int(17.0/delta_gr)
        print*, 'ngr=',ngr
        do i=1,ngr
             rbajo=dfloat(i-1)*delta_gr
             ralto=rbajo+delta_gr
             dist=rbajo+delta_gr/2.d00
             cideal=const2*(ralto*ralto*ralto-rbajo*rbajo*rbajo)
             gr_av(i)=gr_av(i)/cideal/dble(nconf)
             write(30,*) dist,gr_av(i)
        enddo
        ngr=int(0.5*box(1)/delta_gr2)
        do i=1,ngr
             rbajo=dfloat(i-1)*delta_gr2
             ralto=rbajo+delta_gr2
             dist=rbajo+delta_gr2/2.d00
             cideal=const*(ralto*ralto*ralto-rbajo*rbajo*rbajo)
             vol=const2*(ralto**3-rbajo**3.d0)
             enerR_av(i)=enerR_av(i)/dble(nconf)
             rhoR_av(i)=rhoR_av(i)/dble(nconf)
             vol=const2*dist**3
             if(dist.gt.0.50) write(32,*) dist,rhoR_av(i),enerR_av(i)
        enddo
        close(30)
        close(27)
        close(28)
         open(74,file='concentration.dat')
         write(74,*) 4,xn_c1_av,xn_c1_av/dble(nconf)
         write(74,*) 4,xn_c2_av,xn_c2_av/dble(nconf)
         write(74,*) 5,xn_c3_av,xn_c3_av/dble(nconf)
         write(74,*) 6,xn_c4_av,xn_c4_av/dble(nconf)
         write(74,*) 7,xn_c5_av,xn_c5_av/dble(nconf)
         close (74)
         open(74,file='histo_NClust.dat')
         do i=1,natoms
              if(histo_Clust_av(i).gt.0.d0) then
                 histo_Clust_av(i)=histo_Clust_av(i)/dble(nconf)
                 write(74,*) i,histo_Clust_av(i)
              endif
         enddo
         close(74)

	stop   
        end
c
c************************************************o 
c
        subroutine spherical(dx,dy,dz,theta,phi)
        implicit none
        double precision dx,dy,dz,theta,phi,xpi,xpi4
        common /constants/ xpi,xpi4
        
        if(dz.gt.0.d0) then
           theta=atan(dsqrt(dx*dx+dy*dy)/dz)
        elseif(dz.eq.0.d0) then
           theta=xpi/2.d0
        else
           theta=xpi+atan(dsqrt(dx*dx+dy*dy)/dz)
        endif
        if(dx.gt.0.d0) then
           if(dy.gt.0.d0) then
              phi=atan(dy/dx)
           else
              phi=2.d0*xpi+atan(dy/dx)
           endif
        elseif(dx.eq.0.d0) then
           if(dy.gt.0.d0) then
              phi=xpi/2.d0
           else
              phi=-xpi/2.d0
           endif
        else
           phi=xpi+atan(dy/dx)
        endif
         
        return
        end
c
c*************************************************************************************
c
        subroutine legendre(pl,x)
        implicit none
        integer m
        double precision x,dx
        double precision pl(0:12,-12:12)

        dx=1.d0-x*x

        ! l=0
        pl(0,0)=1.d0

        ! l=1
        pl(1,0)=x
        pl(1,1)=-dsqrt(dx)
        pl(1,-1)=-pl(1,1)/2.d0

        ! l=2
        pl(2,0)=(3.d0*x*x-1.d0)/2.d0
        pl(2,1)=-3.d0*x*dsqrt(dx)
        pl(2,-1)=-pl(2,1)/6.d0
        pl(2,2)=2.d0*dx
        pl(2,-2)=pl(2,2)/24.d0

        ! l=3
        pl(3,0)=(5.d0*x*x*x-3.d0*x)/2.d0
        pl(3,1)=-(3.d0/2.d0)*(5.d0*x*x-1.d0)*dsqrt(dx)
        pl(3,-1)=-pl(3,1)/12.d0
        pl(3,2)=15.d0*x*dx
        pl(3,-2)=pl(3,2)/120.d0
        pl(3,3)=-15.d0*dx**(3./2.)
        pl(3,-3)=-pl(3,3)/720.

        ! l=4
        pl(4,0)=(35.d0*x**4.d0-30.d0*x**2.+3.d0)/8.d0
        pl(4,1)=-(5.d0/2.d0)*(7.d0*x**3.-3.d0*x)*dsqrt(dx)
        pl(4,2)=(15.d0/2.0)*(7.d0*x**2.d0-1.d0)*dx
        pl(4,3)=-105.d0*x*dx**(3.d0/2.d0)
        pl(4,4)=105.d0*dx**2.d0
        pl(4,-1)=-pl(4,1)/20.d0
        pl(4,-2)=pl(4,2)/360.d0
        pl(4,-3)=-pl(4,3)/5040.d0
        pl(4,-4)=pl(4,4)/40320.d0

        ! l=5
        pl(5,0)=(63.d0*x**5.-70.d0*x**3.+15.d0*x)/8.d0
        pl(5,1)=-(15./8.)*dsqrt(dx)*(21.*x**4.-14.*x**2.+1.)
        pl(5,-1)=-pl(5,1)/30.d0
        pl(5,2)=(105.d0/2.d0)*dx*(3.d0*x**3.d0-x)
        pl(5,-2)=pl(5,2)/840.d0
        pl(5,3)=-(105.d0/2.d0)*dx**(3.d0/2.d0)*(9.d0*x**2-1.d0)
        pl(5,-3)=-pl(5,3)/20160.d0
        pl(5,4)=945.d0*dx**2.d0*x
        pl(5,-4)=pl(5,4)/362880.d0
        pl(5,5)=-945.d0*dx**(5.d0/2.d0)
        pl(5,-5)=-pl(5,5)/3628800.d0
!
        ! l=6
        pl(6,0)=(231.d0*x**6.d0-315.d0*x**4.d0+105.d0*x**2.d0-5.0)/16.0
        pl(6,1)=-(21./8.)*dsqrt(dx)*(33.*x**5.-30.*x**3.+5.d0*x)
        pl(6,-1)=-pl(6,1)/42.d0
        pl(6,2)=(105.d0/8.d0)*dx*(33.d0*x**4.d0-18.d0*x**2.d0+1.d0)
        pl(6,-2)=pl(6,2)/1680.d0
        pl(6,3)=-(315.d0/2.d0)*dx**(3./2.)*(11.d0*x**3.d0-3.d0*x)
        pl(6,-3)=-pl(6,3)/60480.d0
        pl(6,4)=(945.d0/2.d0)*dx**2.d0*(11.d0*x**2.d0-1.d0)
        pl(6,-4)=pl(6,4)/1814400.d0
        pl(6,5)=-10395.d0*x*dx**(5.d0/2.d0)
        pl(6,-5)=-pl(6,5)/39916800.
        pl(6,6)=10395.*dx**3.d0
        pl(6,-6)=pl(6,6)/479001600.d0
       ! l=7
        pl(7,0)=(-35.d0*x + 315.d0*x**3. - 693.d0*x**5. +
     +          429.d0*x**7.)/16.d0
        pl(7,1)=(-7.d0*dsqrt(dx)*(-5.d0+135.d0*x**2.d0-495.d0*x**4.d0+429.d0*x**6.d0))/16.d0
        pl(7,-1)=(dsqrt(dx)*(-5.d0+135.d0*x**2.d0-495.d0*x**4.d0+429.d0*x**6.d0))/128.d0
        pl(7,2)=(63.d0*dx*(15.d0*x-110.d0*x**3.d0+143.d0*x**5.d0))/8.d0
        pl(7,-2)=-(-dx*(15.d0*x-110.d0*x**3.d0+143.d0*x**5.d0))/384.d0
        pl(7,3)=(-315.d0*dx**1.5d0*(3.d0-66.d0*x**2.d0+143.d0*x**4.d0))/8.d0
        pl(7,-3)=(dx**1.5*(3.d0-66.d0*x**2.0+143.d0*x**4.d0))/3840.d0
        pl(7,4)=(3465.d0*(-dx)**2.d0*(-3.d0*x+13.d0*x**3.d0))/2.d0
        pl(7,-4)=((-dx)**2.d0*(-3.d0*x+13.d0*x**3.d0))/3480.d0
        pl(7,5)=(-10395.d0*dx**2.5d0*(-1.d0+13.d0*x**2.d0))/2.d0
        pl(7,-5)=(dx**2.5*(-1.d0+13.d0*x**2.d0))/46080.d0
        pl(7,6)=-135135.d0*x*(-dx)**3.d0
        pl(7,-6)=-(x*(-dx)**3.d0)/46080.d0
        pl(7,7)=-135135.d0*dx**3.5d0
        pl(7,-7)=dx**3.5/645120.d0

        ! l=8
        pl(8,0)=(35.d0-1260.d0*x**2.d0+6930.d0*x**4.d0-12012.d0*x**6.d0+
     +          6435.d0*x**8.d0)/128.d0
        pl(8,1)=(-9.d0*dsqrt(dx)*(-35.d0*x+385.d0*x**3.d0-1001.d0*x**5.d0+715.d0*x**7.d0))/16.d0
        pl(8,-1)=(dsqrt(dx)*(-35.d0*x+385.d0*x**3.d0-1001.d0*x**5.d0+715.d0*x**7.d0))/128.d0
        pl(8,2)=(-315.d0*(-dx)*(-1.d0+33.d0*x**2.d0-143.d0*x**4.d0+143.d0*x**6.d0))/16.d0
        pl(8,-2)=(dx*(-1.d0+33.d0*x**2.d0-143.d0*x**4.d0+143.d0*x**6.d0))/256.d0
        pl(8,3)=(-3465.d0*dx**1.5d0*(3.d0*x-16.d0*x**3.d0+39.d0*x**5.d0))/8.d0
        pl(8,-3)=(dx**1.5d0*(3.d0*x-26.d0*x**3.d0+39.d0*x**5.d0))/768.d0
        pl(8,4)=(10395.d0*(-dx)**2.d0*(1.d0-26.d0*x**2.d0+65.d0*x**4.d0))/8.d0
        pl(8,-4)=(dx**2.d0*(1.d0-26.d0*x**2.d0+65.d0*x**4.d0))/15360.d0
        pl(8,5)=(-135135.d0*dx**2.5d0*(-x+5.d0*x**3.d0))/2.d0
        pl(8,-5)=(dx**2.5d0*(-x+5.d0*x**3.d0))/15360.d0
        pl(8,6)=(-135135.d0*(-dx)**3.d0*(-1.d0+15.d0*x**2.d0))/2.d0
        pl(8,-6)=-((-dx)**3.d0*(-1.d0+15.d0*x**2.0))/645120.d0
        pl(8,7)=-2027025.d0*x*dx**3.d5
        pl(8,-7)=(x*dx**3.5d0)/645120.d0
        pl(8,8)=2027025.d0*(-dx)**4.d0
        pl(8,-8)=((-dx)**4.d0)/1.032192d7
        ! l=9
        pl(9,0)=(315.d0*x-4620.d0*x**3.d0+18018.d0*x**5.d0-25740.d0*x**7.d0+12155.d0*x**9.d0)/128.d0
        pl(9,1)=(-45.d0*dsqrt(dx)*(7.d0-308.d0*x**2.d0+2002.d0*x**4.d0-4004.d0*x**6.d0+2431.d0*x**8.d0))/128.d0
        pl(9,-1)=(dsqrt(dx)*(7.d0-308.d0*x**2.d0+2002.d0*x**4.d0-4004.d0*x**6.d0+2431.d0*x**8.d0))/256.d0
        pl(9,2)=(-495.d0*(-dx)*(-7.d0*x+91.d0*x**3.d0-273.d0*x**5.d0+221.d0*x**7.d0))/16.d0
        pl(9,-2)=-((-dx)*(-7.d0*x+91.d0*x**3.d0-273.d0*x**5.d0+221.d0*x**7.d0))/256.d0
        pl(9,3)=(-3465.d0*dx**1.5d0*(-1.d0+39.d0*x**2.d0-195.d0*x**4.d0+221.d0*x**6.d0))/16.d0
        pl(9,-3)=(dx**1.5d0*(-1.d0+39.d0*x**2.0-195.d0*x**4.d0+221.d0*x**6.d0))/3072.d0
        pl(9,4)=(135135.d0*(-dx)**2.d0*(x-10.d0*x**3.d0+17.d0*x**5.d0))/8.d0
        pl(9,-4)=((-dx)**2.d0*(x-10.d0*x**3.d0+17.d0*x**5.d0))/3072.d0
        pl(9,5)=(-135135.d0*dx**2.5d0*(1.d0-30.d0*x**2.d0+85.d0*x**4.d0))/8.d0
        pl(9,-5)=(dx**2.5d0*(1.d0-30.d0*x**2.d0+85.d0*x**4.d0))/215040.d0
        pl(9,6)=(-675675.d0*(-dx)**3.d0*(-3.d0*x+17.d0*x**3.d0))/2.d0
        pl(9,-6)=-((-dx)**3.d0*(-3.d0*x+17.d0*x**3.d0))/645120.d0
        pl(9,7)=(-2027025.d0*dx**3.5*(-1.d0+17.d0*x**2.d0))/2.d0
        pl(9,-7)=(dx**3.5d0*(-1.d0+17.d0*x**2.d0))/1.032192d7
        pl(9,8)=34459425.d0*x*(-dx)**4.d0
        pl(9,-8)=(x*(-dx)**4.d0)/1.032192d7
        pl(9,9)=-34459425.d0*dx**4.5
        pl(9,-9)=dx**4.5d0/1.8579456d8
!l=10
        pl(10,0)=(-63.d0+3465.d0*x**2.d0-30030.d0*x**4.d0+90090.d0*x**6.d0-109395.d0*x**8.d0+46189.d0*x**10.d0)/256.d0
        pl(10,1)=(-55.d0*dsqrt(dx)*(63.d0*x-1092.d0*x**3.d0+4914.d0*x**5.d0-7956.d0*x**7.d0+4199.d0*x**9.0))/128.d0
        pl(10,-1)=(dsqrt(dx)*(63.d0*x-1092.d0*x**3.d0+4914.d0*x**5.d0-7956.d0*x**7.d0+4199.d0*x**9.d0))/256.d0
        pl(10,2)=(-495.d0*(-dx)*(7.d0-364.d0*x**2.d0+2730.d0*x**4.d0-6188.d0*x**6.d0+4199.d0*x**8.d0))/128.d0
        pl(10,-2)=-((-dx)*(7.d0-364.d0*x**2.d0+2730.d0*x**4.d0-6188.0*x**6.d0+4199.d0*x**8.d0))/3072.d0
        pl(10,3)=(-6435.d0*dx**1.5d0*(-7.d0*x+105.d0*x**3.d0-357.d0*x**5.d0+323.d0*x**7.d0))/16.d0
        pl(10,-3)=(dx**1.5d0*(-7.d0*x+105.d0*x**3.d0-357.d0*x**5.d0+323.d0*x**7.d0))/3072.d0
        pl(10,4)=(45045.d0*(-dx)**2.d0*(-1.d0+45.d0*x**2.0-255.d0*x**4.d0+323.d0*x**6.d0))/16.d0
        pl(10,-4)=((-dx)**2.d0*(-1.d0+45.d0*x**2.d0-255.d0*x**4.d0+323.d0*x**6.d0))/43008.d0
        pl(10,5)=(-135135.d0*dx**2.5*(15.d0*x-170.d0*x**3.d0+323.d0*x**5.d0))/8.d0
        pl(10,-5)=(dx**2.5*(15.d0*x-170.d0*x**3.d0+323.d0*x**5.d0))/645120.d0
        pl(10,6)=(-675675.d0*(-dx)**3.d0*(3.d0-102.d0*x**2.d0+323.d0*x**4.d0))/8.d0
        pl(10,-6)=-((-dx)**3.d0*(3.d0-102.d0*x**2.0+323.d0*x**4.d0))/1.032192d7
        pl(10,7)=(-11486475.d0*dx**3.5d0*(-3.d0*x+19.d0*x**3.d0))/2.d0
        pl(10,-7)=(dx**3.5d0*(-3.d0*x+19.d0*x**3.d0))/1.032192d7
        pl(10,8)=(34459425.d0*(-dx)**4.d0*(-1.d0+19.d0*x**2.d0))/2.d0
        pl(10,-8)=((-dx)**4.d0*(-1.d0+19.d0*x**2.d0))/1.8579456e8
        pl(10,9)=-654729075.d0*x*dx**4.5d0
        pl(10,-9)=(x*dx**4.5d0)/1.8579456d8
        pl(10,10)=-654729075.d0*(-dx)**5.d0
        pl(10,-10)=-((-dx)**5.d0)/3.7158912d9

        !l=11
        pl(11,0)=(-693.d0*x+15015.d0*x**3.d0-90090.d0*x**5.d0+218790.d0*x**7.d0-230945.d0*x**9.d0+88179.d0*x**11.d0)/256.d0
        pl(11,1)=(-33.d0*dsqrt(dx)*(-21.d0+1365.d0*x**2.d0-13650.d0*x**4.d0+46410.d0*x**6.d0-62985.d0*x**8.d0+29393.d0*x**10.d0))/256.d0
        pl(11,-1)=(dsqrt(dx)*(-21.d0+1365.d0*x**2.d0-13650.d0*x**4.d0+46410.d0*x**6.d0-62985.d0*x**8.d0+29393.d0*x**10.d0))/1024.d0
        pl(11,2)=(-2145.d0*(-dx)*(21.d0*x-420.d0*x**3.d0+2142.d0*x**5.d0-3876.d0*x**7.d0+2261.d0*x**9.d0))/128.d0
        pl(11,-2)=-((-dx)*(21.d0*x-420.d0*x**3.d0+2142.d0*x**5.d0-3876.d0*x**7.d0+2261.d0*x**9.d0))/1024.d0
        pl(11,3)=(-45045.d0*dx**1.5d0*(1.d0-60.d0*x**2.d0+510.d0*x**4.d0-1292.d0*x**6.d0+969.d0*x**8.d0))/128.d0
        pl(11,-3)=(dx**1.5d0*(1.d0-60.d0*x**2.d0+510.d0*x**4.d0-1292.d0*x**6.d0+969.d0*x**8.d0))/6144.d0
        pl(11,4)=(135135.d0*(-dx)**2.d0*(-5.d0*x+85.d0*x**3.d0-323.d0*x**5.d0+323.d0*x**7.d0))/16.d0
        pl(11,-4)=((-dx)**2.d0*(-5.d0*x+85.d0*x**3.d0-323.d0*x**5.d0+323.d0*x**7.d0))/30720.d0
        pl(11,5)=(-135135.d0*dx**2.5d0*(-5.d0+255.d0*x**2.d0-1615.d0*x**4.d0+2261.d0*x**6.d0))/16.d0
        pl(11,-5)=(dx**2.5*(-5.d0+255.d0*x**2.d0-1615.d0*x**4.d0+2261.d0*x**6.d0))/3.44064d6
        pl(11,6)=(-2297295.d0*(-dx)**3.d0*(15.d0*x-190.d0*x**3.d0+399.d0*x**5.d0))/8.d0
        pl(11,-6)=-((-dx)**3.d0*(15.d0*x-190.d0*x**3.d0+399.d0*x**5.d0))/1.032192d7
        pl(11,7)=(-34459425.d0*dx**3.5*(1.d0-38.d0*x**2.d0+133.d0*x**4.d0))/8.d0
        pl(11,-7)=(dx**3.5d0*(1.d0-38.d0*x**2.d0+133.d0*x**4.d0))/6.193152d7
        pl(11,8)=(654729075.d0*(-dx)**4.d0*(-x+7.d0*x**3.d0))/2.d0
        pl(11,-8)=((-dx)**4.d0*(-x+7.d0*x**3.d0))/6.193152d7
        pl(11,9)=(-654729075.d0*dx**4.5*(-1.d0+21.d0*x**2.d0))/2.d0
        pl(11,-9)=(dx**4.5d0*(-1.d0+21.d0*x**2.d0))/3.7158912d9
        pl(11,10)=-13749310575.d0*x*(-dx)**5.d0
        pl(11,-10)=-(x*(-dx)**5.d0)/3.7158912d9
        pl(11,11)=-13749310575.d0*dx**5.5d0
        pl(11,-11)=dx**5.5d0/8.17496064d10


        ! l12
        pl(12,0)= (231. - 18018.*x**2. + 225225.*x**4. - 1021020.*x**6. +
     -    2078505.*x**8. - 1939938.*x**10. + 676039.*x**12.)/1024.
        pl(12,1)=    (-39.*Sqrt(1. - x**2.)*
     -    (-231.*x + 5775.*x**3. - 39270.*x**5. + 106590.*x**7. -
     -      124355.*x**9. + 52003.*x**11.))/256.
        pl(12,-1)=   (Sqrt(1. - x**2.)*(-231.*x + 5775.*x**3. - 39270.*x**5. +
     -      106590.*x**7. - 124355.*x**9. + 52003.*x**11.))/1024.
        pl(12,2)=     (-3003.*(-1. + x**2.)*
     -    (-3. + 225.*x**2. - 2550.*x**4. + 9690.*x**6. - 14535.*x**8. +
     -      7429.*x**10.))/256.
        pl(12,-2)=   -((-1. + x**2.)*(-3. + 225.*x**2. - 2550.*x**4. + 9690.*x**6. -
     -       14535.*x**8. + 7429.*x**10.))/2048.
        pl(12,3)=        (15015.*Sqrt(1. - x**2.)*(-1. + x**2.)*
     -    (45.*x - 1020.*x**3. + 5814.*x**5. - 11628.*x**7. +
     -      7429.*x**9.))/128.
        pl(12,-3)= -(Sqrt(1. - x**2.)*(-1. + x**2.)*
     -     (45.*x - 1020.*x**3. + 5814.*x**5. - 11628.*x**7. +
     -       7429.*x**9.))/30720.
        pl(12,4)=        (135135.*(-1. + x**2.)**2.*
     -    (5. - 340.*x**2. + 3230.*x**4. - 9044.*x**6. + 7429.*x**8.))/
     -  128.
        pl(12,-4)=   ((-1. + x**2.)**2.*(5. - 340.*x**2. + 3230.*x**4. - 9044.*x**6. +
     -      7429.*x**8.))/491520.
        pl(12,5)=        (-2297295.*Sqrt(1. - x**2.)*(-1. + x**2.)**2.*
     -    (-5.*x + 95.*x**3. - 399.*x**5. + 437.*x**7.))/16.
        pl(12,-5)=        (Sqrt(1. - x**2.)*(-1. + x**2.)**2.*
     -    (-5.*x + 95.*x**3. - 399.*x**5. + 437.*x**7.))/491520.
        pl(12,6)=        (-2297295.*(-1. + x**2.)**3.*
     -    (-5. + 285.*x**2. - 1995.*x**4. + 3059.*x**6.))/16.
        pl(12,-6)= -((-1. + x**2.)**3.*(-5. + 285.*x**2. - 1995.*x**4. +
     -       3059.*x**6.))/6.193152e7
        pl(12,7)=        (130945815.*Sqrt(1. - x**2.)*(-1. + x**2.)**3.*
     -    (5.*x - 70.*x**3. + 161.*x**5.))/8.
        pl(12,-7)=   -(Sqrt(1. - x**2.)*(-1. + x**2.)**3.*
     -     (5.*x - 70.*x**3. + 161.*x**5.))/6.193152e7
        pl(12,8)=(654729075.*(-1. + x**2.)**4.*(1. - 42.*x**2. + 161.*x**4.))/8.
        pl(12,-8)=((-1. + x**2.)**4.*(1. - 42.*x**2. + 161.*x**4.))/1.2386304e9
        pl(12,9)=(-4583103525.*Sqrt(1. - x**2.)*(-1. + x**2.)**4.*(-3.*x + 23.*x**3.))/2.
        pl(12,-9)=(Sqrt(1. - x**2.)*(-1. + x**2.)**4.*(-3.*x + 23.*x**3.))/3.7158912e9
        pl(12,10)=(-13749310575.*(-1. + x**2.)**5.*(-1. + 23.*x**2.))/2.
        pl(12,-10)=-((-1. + x**2.)**5.*(-1. + 23.*x**2.))/8.17496064e10
        pl(12,11)=-316234143225.*x*(1. - x**2.)**5.5
        pl(12,-11)=(x*(1. - x**2.)**5.5)/8.17496064e10
        pl(12,12)=316234143225.*(-1 + x**2.)**6.
        pl(12,-12)=(-1. + x**2.)**6./1.9619905536e12

        return
        end
c
c*****************************************************************************************
c
        subroutine SphericalHarmonics(pl,phi,hs_re,hs_im)
        implicit none
        integer l,m
        double precision xpi,xpi4
        double precision a,b,c,d,e
        double precision hs_re(0:12,-12:12)
        double precision hs_im(0:12,-12:12)
        double precision pl(0:12,-12:12),phi
        double precision factorial

        common /constants/ xpi,xpi4

        do l=0,12
           a=(2.d0*dble(l)+1.d0)
           do m=-l,l
               b=(-1.d0)**dble(m)
               c=factorial(l-m)
               d=factorial(l+m)
               e=dble(m)*phi
               hs_re(l,m)=b*dsqrt(a*c/d/xpi4)*pl(l,m)*dcos(e)
               hs_im(l,m)=b*dsqrt(a*c/d/xpi4)*pl(l,m)*dsin(e)
               !write(*,'("l,m",2i5,3f10.3)') l,m,hs_re(l,m),hs_im(l,m),phi
           enddo
        enddo
        return
        end

c
c**************************************************************************
c
        function factorial(n)
        implicit none
        integer n,i
        double precision factorial

        factorial=1.d0
        do i=1,n
           factorial=factorial*dble(i)
        enddo

        return 
        end
c
c**************************************************************************
c
        subroutine init_clusters(nClust)
        implicit none
        integer i,nClust,nmax,nthres
        integer inClust(120000),clustHead(120000),clustSize(120000)
        integer clustNext(120000),clustSize_order(120000)
        common /clust1/ inClust,clustHead,clustSize,clustNext
        common /clust2/ clustSize_order


        do i=1,120000 !/* initialize at the beginning
            inClust(i)=-1
            clustHead(i)=0
            clustNext(i)=0
            clustSize(i)=0
            clustSize_order(i)=0
        enddo
        nClust=0
        nthres=0

        return
        end
c
c**************************************************************************
c
        subroutine cluster_analysis(neigh,nneigh,nmax,nvmax,nmol,nClust_min,nClust)
        implicit none
        integer i,j,indxj,nmol,nmax,nvmax,nthres
        integer inClust(120000),clustHead(120000),clustSize(120000)
        integer clustNext(120000),neigh(nmax,nvmax),nneigh(nmax)
        integer clustSize_order(120000),clustIndex_order(nmax)
        integer nc,nClust,nc1,nc2,m,mp,Nliq,nclust_min
        integer nSize,nSizeMax,cBig,cSmall,nClustSize(120000)
        integer clustHead_order(120000)
        common /clust1/ inClust,clustHead,clustSize,clustNext
        common /clust2/ clustSize_order
c
        do i=1,nmax
            clustSize_order(i)=0
        enddo
        call init_clusters(nClust)

        do i=1,nmol 
c
              do j=1,nneigh(i)
                indxj=neigh(i,j)
           
                   if(inClust(i).eq.-1) then ! i has not yet been assigned
                     if(inClust(indxj).eq.-1) then ! indxj has not yet been assigned
                       nClust=nClust+1
                       inClust(i)=nClust
                       inClust(indxj)=nClust
                       clustSize(nClust)=2
                       clustHead(nClust)=i
                       clustNext(i)=indxj
                       clustNext(indxj)=0
                     else !indxj already belongs to a cluster
                       nc=inClust(indxj)
                       clustSize(nc)=clustSize(nc)+1
                       inClust(i)=nc
                       clustNext(i)=clustHead(nc)
                       clustHead(nc)=i
                     endif

                   else   ! i already belongs to a cluster

                     if(inClust(indxj).eq.-1) then ! indxj has not yet been assigned
                       nc=inClust(i)
                       clustSize(nc)=clustSize(nc)+1
                       inClust(indxj)=nc
                       clustNext(indxj)=clustHead(nc)
                       clustHead(nc)=indxj  
                     else  ! both already belong to a cluster
                              ! careful here, we might be following one
                              ! or the two clusters
                       nc1=inClust(i)
                       nc2=inClust(indxj)
                       if(nc1.ne.nc2) then
                          if(nc1.le.nthres.and.nc2.gt.nthres) then
                                cBig=nc1
                          elseif(nc1.gt.nthres.and.nc2.le.nthres) then
                                cBig=nc2
                          elseif(nc1.le.nthres.and.nc2.le.nthres) then  ! fo now, I will compress all the clusters
                             if(clustSize(nc1).ge.clustSize(nc2)) then
                                cBig=nc1
                             else
                                cBig=nc2
                             endif
                          else
                             if(clustSize(nc1).ge.clustSize(nc2)) then
                                cBig=nc1
                             else
                                cBig=nc2
                             endif
                          endif
                          cSmall=nc1+nc2-cBig
 
                          m=clustHead(cSmall)
                          do while(m.gt.0) 
                             inClust(m)=cBig
                             mp=m
                             m=clustNext(m)
                          enddo
                          clustNext(mp)=clustHead(cBig)
                          clustHead(cBig)=clustHead(cSmall)
                          clustSize(cBig)=clustSize(cBig)+
     +                      clustSize(cSmall)
                          clustSize(cSmall)=0
                     endif
                    
                   endif
                   endif

              enddo
        enddo 
        !print*, 'assignment done'
        !print*, 'now compressing the clusters'
        !print*, 'nClust',nClust
c
        ! compressing the clusters
c
        nc=0
        do i=1,nClust
           !print*, 'clustSize(i)',clustSize(i)
           if(clustSize(i).gt.0) then
               nc=nc+1
               clustSize(nc)=clustSize(i)
               clustHead(nc)=clustHead(i)
               m=clustHead(nc)
               do while(m.gt.0) 
                  inClust(m)=nc
                  m=clustNext(m)
               enddo
           endif
        enddo
        nClust=nc
        print*, 'nClust=',nClust

        ! analysis of cluster size
        cBig=1
        do i=1,nClust
           nSize=clustSize(i)
           if(nSize.gt.clustSize(cBig)) cBig=i
           nClustSize(nSize)=nClustSize(nSize)+1
        enddo
        nSizeMax=clustSize(cBig)
        print*, 'Larger Cluster ',NSizeMax

        Nliq=0
        do i=1,nmol
          if(inClust(i).eq.0) then
             Nliq=Nliq+1
          endif
        enddo
        print*, 'Nliq=',Nliq

        do i=1,nClust

           CBig=1
           do j=1,nClust!+1-i
               nSize=clustSize(j)
               if(nSize.gt.clustSize(cBig)) cBig=j
           enddo
           clustSize_order(i)=clustSize(cBig)
           clustHead_order(i)=clustHead(cBig)
           m=clustHead(cBig)  ! y tambien cambio la etiqueta en los atomos
           do while(m.gt.0)
               inClust(m)=i
               m=clustNext(m)
           enddo
           
           clustSize(cBig)=0

        enddo

        !do i=1,nClust
        !    print*, 'i, size=',i,clustSize_order(i),clustIndex_order(i)
        !enddo
c
!       just keep assigned those atoms belonging to the first largest five clusters
c
        if(clustSize_order(5).ge.nClust_min) then
            nthres=5
        elseif(clustSize_order(4).ge.nClust_min) then
            nthres=4
        elseif(clustSize_order(3).ge.nClust_min) then
            nthres=3
        elseif(clustSize_order(2).ge.nClust_min) then
            nthres=2
        elseif(clustSize_order(1).ge.nClust_min) then
            nthres=1
        else
            nthres=0
        endif
        print*, 'following',nthres,'clusters',nclust_min

!        if(nthres.ge.1) then
!            do i=1,nmol
!               if(inClust(i).gt.nthres) then
!                 inClust(i)=-1  ! desassign those belonging to smaller clusters
!               elseif(inClust(i).eq.0) then
!                 inClust(i)=-1  ! or liquid
!               endif
!            enddo
!            do i=1,nthres
!               clustHead(i)=clustHead_order(i)
!               clustSize(i)=clustSize_order(i)
!            enddo
!            do i=nthres+1,nmax   
!               clustHead(i)=0
!               clustSize(i)=0
!            enddo
!            nClust=nthres
!        else
!            do i=1,nmax   
!               clustHead(i)=0
!               clustSize(i)=0
!            enddo
!            nClust=0
!            do i=1,nmol
!                inClust(i)=-1
!            enddo
!        endif

        return
        end


c**************************************************************************
c
        subroutine cluster_analysis_dotq(isol,neigh,nneigh,nmax,nvmax,nmol,nClust_min,nClust)
        implicit none
        integer i,j,indxj,nmol,nmax,nvmax,nthres
        integer inClust(120000),clustHead(120000),clustSize(120000)
        integer clustNext(120000),neigh(nmax,nvmax),nneigh(nmax)
        integer clustSize_order(120000),clustIndex_order(nmax)
        integer nc,nClust,nc1,nc2,m,mp,Nliq,nclust_min
        integer nSize,nSizeMax,cBig,cSmall,nClustSize(120000)
        integer clustHead_order(120000),isol(nmax)
        common /clust1/ inClust,clustHead,clustSize,clustNext
        common /clust2/ clustSize_order
c
        do i=1,nmax
            clustSize_order(i)=0
        enddo
        call init_clusters(nClust)

        do i=1,nmol 
c
          if(isol(i).eq.0) then ! es líquido
                 
                 inClust(i)=0

          else  ! es sólido

              do j=1,nneigh(i)
                indxj=neigh(i,j)

                if(isol(indxj).eq.0)  then

                    inClust(indxj)=0

                else
           
                   if(inClust(i).eq.-1) then ! i has not yet been assigned
                     if(inClust(indxj).eq.-1) then ! indxj has not yet been assigned
                       nClust=nClust+1
                       inClust(i)=nClust
                       inClust(indxj)=nClust
                       clustSize(nClust)=2
                       clustHead(nClust)=i
                       clustNext(i)=indxj
                       clustNext(indxj)=0
                     else !indxj already belongs to a cluster
                       nc=inClust(indxj)
                       clustSize(nc)=clustSize(nc)+1
                       inClust(i)=nc
                       clustNext(i)=clustHead(nc)
                       clustHead(nc)=i
                     endif

                   else   ! i already belongs to a cluster

                     if(inClust(indxj).eq.-1) then ! indxj has not yet been assigned
                       nc=inClust(i)
                       clustSize(nc)=clustSize(nc)+1
                       inClust(indxj)=nc
                       clustNext(indxj)=clustHead(nc)
                       clustHead(nc)=indxj  
                     else  ! both already belong to a cluster
                              ! careful here, we might be following one
                              ! or the two clusters
                       nc1=inClust(i)
                       nc2=inClust(indxj)
                       if(nc1.ne.nc2) then
                          if(nc1.le.nthres.and.nc2.gt.nthres) then
                                cBig=nc1
                          elseif(nc1.gt.nthres.and.nc2.le.nthres) then
                                cBig=nc2
                          elseif(nc1.le.nthres.and.nc2.le.nthres) then  ! fo now, I will compress all the clusters
                             if(clustSize(nc1).ge.clustSize(nc2)) then
                                cBig=nc1
                             else
                                cBig=nc2
                             endif
                          else
                             if(clustSize(nc1).ge.clustSize(nc2)) then
                                cBig=nc1
                             else
                                cBig=nc2
                             endif
                          endif
                          cSmall=nc1+nc2-cBig
 
                          m=clustHead(cSmall)
                          do while(m.gt.0) 
                             inClust(m)=cBig
                             mp=m
                             m=clustNext(m)
                          enddo
                          clustNext(mp)=clustHead(cBig)
                          clustHead(cBig)=clustHead(cSmall)
                          clustSize(cBig)=clustSize(cBig)+
     +                      clustSize(cSmall)
                          clustSize(cSmall)=0
                     endif
                    
                   endif
                   endif

                   endif !if(isol(indxj).eq.0)  then
              enddo
            endif !if(isol(i).eq.1) then ! es sólido
        enddo 
        !print*, 'assignment done'
        !print*, 'now compressing the clusters'
        !print*, 'nClust',nClust
c
        ! compressing the clusters
c
        nc=0
        do i=1,nClust
           !print*, 'clustSize(i)',clustSize(i)
           if(clustSize(i).gt.0) then
               nc=nc+1
               clustSize(nc)=clustSize(i)
               clustHead(nc)=clustHead(i)
               m=clustHead(nc)
               do while(m.gt.0) 
                  inClust(m)=nc
                  m=clustNext(m)
               enddo
           endif
        enddo
        nClust=nc
        print*, 'nClust=',nClust

        ! analysis of cluster size
        cBig=1
        do i=1,nClust
           nSize=clustSize(i)
           if(nSize.gt.clustSize(cBig)) cBig=i
           nClustSize(nSize)=nClustSize(nSize)+1
        enddo
        nSizeMax=clustSize(cBig)
        print*, 'Larger Cluster ',NSizeMax

        Nliq=0
        do i=1,nmol
          if(inClust(i).eq.0) then
             Nliq=Nliq+1
          endif
        enddo
        print*, 'Nliq=',Nliq

        do i=1,nClust

           CBig=1
           do j=1,nClust!+1-i
               nSize=clustSize(j)
               if(nSize.gt.clustSize(cBig)) cBig=j
           enddo
           clustSize_order(i)=clustSize(cBig)
           clustHead_order(i)=clustHead(cBig)
           m=clustHead(cBig)  ! y tambien cambio la etiqueta en los atomos
           do while(m.gt.0)
               inClust(m)=i
               m=clustNext(m)
           enddo
           
           clustSize(cBig)=0

        enddo

        !do i=1,nClust
        !    print*, 'i, size=',i,clustSize_order(i),clustIndex_order(i)
        !enddo
c
!       just keep assigned those atoms belonging to the first largest five clusters
c
        if(clustSize_order(5).ge.nClust_min) then
            nthres=5
        elseif(clustSize_order(4).ge.nClust_min) then
            nthres=4
        elseif(clustSize_order(3).ge.nClust_min) then
            nthres=3
        elseif(clustSize_order(2).ge.nClust_min) then
            nthres=2
        elseif(clustSize_order(1).ge.nClust_min) then
            nthres=1
        else
            nthres=0
        endif
        print*, 'following',nthres,'clusters',nclust_min

!        if(nthres.ge.1) then
!            do i=1,nmol
!               if(inClust(i).gt.nthres) then
!                 inClust(i)=-1  ! desassign those belonging to smaller clusters
!               elseif(inClust(i).eq.0) then
!                 inClust(i)=-1  ! or liquid
!               endif
!            enddo
!            do i=1,nthres
!               clustHead(i)=clustHead_order(i)
!               clustSize(i)=clustSize_order(i)
!            enddo
!            do i=nthres+1,nmax   
!               clustHead(i)=0
!               clustSize(i)=0
!            enddo
!            nClust=nthres
!        else
!            do i=1,nmax   
!               clustHead(i)=0
!               clustSize(i)=0
!            enddo
!            nClust=0
!            do i=1,nmol
!                inClust(i)=-1
!            enddo
!        endif

        return
        end


c*************************************************************************   
        subroutine calculo_matriz_inversa(h,hinv)
        implicit none
        double precision h(3,3),hinv(3,3),deth

        detH=h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))-h(2,1)*(h(1,2)*h(3,3)
     *  -h(1,3)*h(3,2))+h(3,1)*(h(1,2)*h(2,3)-h(1,3)*h(2,2))
ccj
c
        hinv(1,1)=(h(2,2)*h(3,3)-h(2,3)*h(3,2))/detH
        hinv(2,1)=-(h(2,1)*h(3,3)-h(3,1)*h(2,3))/detH
        hinv(3,1)=(h(2,1)*h(3,2)-h(3,1)*h(2,2))/detH
        hinv(1,2)=-(h(1,2)*h(3,3)-h(3,2)*h(1,3))/detH
        hinv(2,2)=(h(1,1)*h(3,3)-h(3,1)*h(1,3))/detH
        hinv(3,2)=-(h(1,1)*h(3,2)-h(3,1)*h(1,2))/detH
        hinv(1,3)=(h(1,2)*h(2,3)-h(2,2)*h(1,3))/detH
        hinv(2,3)=-(h(1,1)*h(2,3)-h(2,1)*h(1,3))/detH
        hinv(3,3)=(h(1,1)*h(2,2)-h(2,1)*h(1,2))/detH

        return
        end
c****************************************************************************
        subroutine cambio_unidades_caja(natoms,nmax,x,y,z,hinv)
        implicit none
        integer i,natoms,indx,nmax
        double precision x1,y1,z1,x1p,y1p,z1p
        double precision hinv(3,3),x(nmax),y(nmax),z(nmax)

        do i=1,natoms
           x1=x(i)
           y1=y(i)
           z1=z(i)
           x1p=hinv(1,1)*x1+hinv(1,2)*y1+hinv(1,3)*z1
           y1p=hinv(2,1)*x1+hinv(2,2)*y1+hinv(2,3)*z1
           z1p=hinv(3,1)*x1+hinv(3,2)*y1+hinv(3,3)*z1
           x(i)=x1p
           y(i)=y1p
           z(i)=z1p
        enddo

        return
        end

c*******************************************************************
       subroutine cross_product(v1,v2,v3)
       implicit none
       double precision v1(3),v2(3),v3(3)

       v3(1)= v1(2)*v2(3)-v2(2)*v1(3)
       v3(2)= v1(3)*v2(1)-v2(3)*v1(1)
       v3(3)= v1(1)*v2(2)-v2(1)*v1(2)

       return
       end
c*********************************************************************
       subroutine order(n,nnmax,v1)
       implicit none
       integer i,j,n,indx_min,nnmax
       double precision vmin,v1(nnmax),v1_order(nnmax)

       do j=1,n
          vmin=+1000.d0
          do i=1,n
             if(v1(i).lt.vmin) then
                vmin=v1(i)
                indx_min=i
             endif
          enddo
          v1_order(j)=vmin
          v1(indx_min)=+2000.d0
       enddo

       do i=1,n
          v1(i)=v1_order(i)
       enddo
       return
       end
c**********************************************************************
       subroutine ener(i,j,en,overlap)
       implicit none
        integer Npart_types_max,Nrb_sites_max, Ntor_max, nmax
        parameter(Npart_types_max=5, Nrb_sites_max = 10, Ntor_max = 3)
        parameter(nmax=120000)
       integer Npart_types,bool_tor, i, j
       integer Nrb_sites(Npart_types_max)
       integer iref_vec(Npart_types_max,Nrb_sites_max)
       integer Ntor_angle(Npart_types_max,Nrb_sites_max)
       integer itypei, itypej, I1MIN, I2MIN
       integer i1,j3, J4, Indx_patch, k
       integer itype(nmax)
       double precision Vangtor_max, Vang_max, Utor_max, v1x, v1y, v1z
       double precision mod_v1, CSOMG1, OMG1, factor_energy, Vtor_max
       double precision CSOMG2, OMG2, VANG, sigma_ij, v2x, v2y, v2z
       double precision mod_v2, cangle, v12x, v12y, v12z, v12b2, torsion
       double precision dtor, Utor, Vangtor, FX, FXI, VRAD
       double precision pi, dospi, a, en
       double precision sigma_tor_jon, rangep, VL0, xop
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
        common /vpot_par/ VL0, xop, sigma_tor_jon, rangeP
        common /array_vpot_par/ vpot_matrix, sigma_jon
        common /nsites/ Nrb_sites, iref_vec, Ntor_angle
       common /array_tor/ tor_angle

        ! coordinates
        common /box1/ h
        common /coords/ x, y, z
        common /coordspacth/ patchx, patchy, patchz
        common /coordstype/ itype

       pi = acos(-1.d0)
       dospi = 2.d0*pi
       a = 1000.d0
       bool_tor=0

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

       DX=DX-NINT(DX)
       DY=DY-NINT(DY)
       DZ=DZ-NINT(DZ)

       DXP=H(1,1)*DX+H(1,2)*DY+H(1,3)*DZ
       DYP=H(2,1)*DX+H(2,2)*DY+H(2,3)*DZ
       DZP=H(3,1)*DX+H(3,2)*DY+H(3,3)*DZ

       RSQ=DXP**2.0+DYP**2.0+DZP**2.
       R12=SQRT(RSQ)

       if(r12.lt.0.8) then
            overlap=.true.
       elseif(r12.lt.0.98) then
            R126=RSQ**3.d0
            R1212=R126*R126
            EN=-4.d0*((r1212 - r126)/(r1212*r126) + VL0)
       elseif(r12.lt.rangep) then
            CSM1M=-100.
            CSM2M=-100.
            I1MIN = -1
            I2MIN = -1

            do i1= 1,nrb_sites(itypei)
               px=patchx(i,i1)
               py=patchy(i,i1)
               pz=patchz(i,i1)
               CSM1 = (PX*DXP+ PY*DYP + PZ*DZP) / R12
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
               CSM2 = -(PX*DXP+ PY*DYP + PZ*DZP) / R12
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

            do j3=1,nrb_sites(itypei)
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
                                CSOMG2 = sign(1.0,CSOMG2)
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
                        EndIf

                    EndIf ! if ( factor_energy > 1.d-08) Then
               EndDo ! Do J4 = 0, Nrb_sites(itypej) -1
            EndDo ! Do J3 = 0, Nrb_sites(itypei) -1

            if(Vangtor_max.eq.-1.d00) then
               En= 0.0
            else
               VANGTOR=Vangtor_max
               FX=(1 + Tanh(a*(-xop + r12)))/2.d00
               FXI=(1 - Tanh(a*(-xop + r12)))/2.d00
               R126=RSQ**3.d0
               R1212=R126*R126
               VRAD=-4.d0*((r1212 - r126)/(r1212*r126) + VL0)
               En = VRAD * (FXI + VANGTOR * FX)
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
        parameter(Npart_types_max=5, Nrb_sites_max = 10, Ntor_max = 3)
       integer istep_ini,istep_fin,Neq, Nmove, Nsave, Nsave2, Nrestart
       integer Nrb_sites(Npart_types_max),i,j,k,bool_tor
       integer iref_vec(Npart_types_max,Nrb_sites_max)
       integer Ntor_angle(Npart_types_max,Nrb_sites_max)
       integer Indx, Indx_tor, Indx_tor2, Indx2, ncount, indx_patch
       integer i1, i2, Indx3, ios, Npart_types
       double precision hmax,omax,temp0,temp1, seed
       double precision sigma_jon_aux, sigma_tor_jon, rangep, VL0, xop
       double precision px, py, pz, px_ref, py_ref, pz_ref
       double precision patch_vec(3*Nrb_sites_max*Npart_types_max)
       double precision ref_tor_vec(3*Nrb_sites_max*Npart_types_max)
       double precision vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
       double precision sigma_jon(Npart_types_max,Nrb_sites_max)
       double precision tor_angle(Npart_types_max,Nrb_sites_max,Ntor_max)

       logical imovie,displ_update, perfect_match
!      parametros que voy a necesitar pasar para el calculo de energia
!      integer arrays: Nrb_sites(Npart_types_max), iref_vec(Nrb_sites_max*Npart_types_max)
!                      Ntor_angle(Nrb_sites_max*Npart_types_max)
!      double scalar: VL0,xop,sigma_tor_jon,rangeP
!      integer scalar: Npart_types, bool_tor
!      double array: vpot_matrix(Npart_types_max,Nrb_sites_max,Npart_types_max,Nrb_sites_max)
!                    sigma_jon(Nrb_sites_max*Npart_types_max)
       common /nsites/ Nrb_sites, iref_vec, Ntor_angle
       common /vpot_par/ VL0, xop, sigma_tor_jon, rangeP
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
       read(20,*) hmax,omax,displ_update,temp0,temp1, seed
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
       elseif( model == 'Jon'.or.model == 'LJG') then
             read(20,*)
             read(20,*) sigma_jon_aux, sigma_tor_jon, rangep, Bool_tor
             VL0=(1.d0/(rangep**12.d0))-(1.d0/(rangep**6.0))
             xop=(1.d0+sqrt(1.d0+4.d00*VL0))/2.d0
             xop=xop**(-1.d0/6.d0)
             sigma_tor_jon = 2.d0*sigma_tor_jon*sigma_tor_jon;
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
                       ! stop
                     elseif(ncount.eq.1) then
                        iref_vec(j,i)=indx_patch
                     else
                        print*, 'ERROR in reading data'
                        print*, 'patch used as reference for torsion'
                        print*, 'matches more than one ang patch',ncount
                       ! stop
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
       do i=1, Npart_types
           Indx = (i-1)*Nrb_sites_max
           do j=1, Nrb_sites(i)
              sigma_jon(I,j) = sigma_jon_aux
           Enddo
       Enddo
       print*, ''
       print*, 'Input file read'
       print*, ''
       close(20)

       return
       end
