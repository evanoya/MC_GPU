MODULE set_precision
! ..
! .. Intrinsic Functions ..
      INTRINSIC KIND
! .. Parameters ..
! Define the standard precisions
      INTEGER, PARAMETER :: skind = KIND(0.0E0)
      INTEGER, PARAMETER :: dkind = KIND(0.0D0)
! Set the precision for the whole package
      INTEGER, PARAMETER :: wp = skind
! To change the default package precision to single precision change
! the parameter assignment to wp above to
!     INTEGER, PARAMETER :: wp = skind
! and recompile the complete package.
END MODULE  set_precision

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE  configuration

   use set_precision
   implicit none

   Integer, Parameter :: ndim=3
   Integer            :: ntypes, natoms

   Integer, Allocatable                :: itype(:)
   Real(wp), Allocatable               :: r(:), q(:)
   Real (wp), Dimension(0:ndim*ndim-1) :: h, hinv
   Real(wp), Dimension(0:ndim-1)       :: Box, Bangles

END MODULE configuration

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE rundata
   use set_precision
   implicit none
   Integer :: istep_ini,istep_fin
   Integer :: ntrans, ntrans_ac, nrot, nrot_ac, no_mc_real
   Integer :: nvol_accept, iscale, nvol
   Integer :: length, Nmove, Neq, Nsave, Nsave2, Nrestart
   Integer :: seed, init=0
   Real(wp) :: temp, temp0, temp1, rate, beta, omax, hmax, DeltaT
   Real(wp) :: vmax(0:2), pres
   Logical  :: imovie, displ_update, npt
   Character*20 base


END MODULE rundata

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE potential
   Use set_precision
   Implicit None

   Integer :: Npart_types, Nrb_sites_max = 15, Ntor_max = 6
   Integer, Allocatable :: Nrb_sites(:), Ntor_angle(:)
   Integer    ::   Bool_tor

   Real(wp) :: cosmax_KF, delta_KF, ctor_max_KF, rangepp
   Real(wp) :: sigma_jon_aux,sigma_tor_jon, rcut_jon
   real(wp), Allocatable :: patch(:), ref_tor_vec(:)
   real(wp), Allocatable :: Vpot_Matrix(:), tor_angle(:)
   real(wp), Allocatable :: sigma_jon(:), sigma_LJ(:)
   real(wp), Allocatable :: xop(:), VL0(:), rangeP(:)
   !character(len=2), Allocatable :: el(:)

END MODULE potential

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE properties
   use set_precision
   Implicit None
   Real (wp) :: en_tot, rho, vol
   Real (wp) :: en_av, rho_av, lx_av, ly_av, lz_av

END MODULE properties

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE utils
   use set_precision
   Implicit None
   real(wp), parameter :: dospi=2.d0*acos(-1.d0)

END MODULE utils

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE Cell 

   use set_precision
   use configuration , only :  ndim
   Implicit None
   Integer                               :: Ntot_Cell, Nat_Cell_Max, NCell_max  ! total number of cells AND maximum number of atoms per cell
   Integer, Dimension(0:ndim-1)          :: Nl_Cell  ! stores the number of cells along x,y, and z
   Integer                               :: Nset_Cell, Nsubset_Cell  ! number of checkerboar sets (as given by dimension) AND number of cell in each set
   Integer, Allocatable, Dimension(:)        :: Nat_Cell  ! stores the number of atoms in each cell
   Integer, Allocatable, Dimension(:)        :: List_Cell ! for each cell stores the indices of the atoms belonging to it
   Integer, Allocatable, Dimension(:)        :: List_CB   ! for each checkerboar set stores the indices of cells belonging to it
   Integer, Allocatable, Dimension(:)        :: Map_Cell  ! for each cell stores the index of the neighbouring cells
   Integer, Allocatable, Dimension(:,:)      :: Indx_Set
   real(wp),Allocatable, Dimension(:)        :: CB_Set   ! indices of the checkerboard sets
   real(wp),Allocatable, Dimension(:)        :: Offset_CB   ! origin of the Cell
   real(wp),Allocatable, Dimension(:)        :: Sphere_Cell
   real(wp),Allocatable, Dimension(:)        :: DR_Cell     
   real(wp)                                  :: W(0:ndim-1)
   real(dkind)                               :: W_RU(0:ndim-1) 
   real(wp)                                  :: disp_max(0:ndim-1)

END MODULE Cell 
