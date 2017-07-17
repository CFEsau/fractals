!******************************************************************************!
!
! Modules for analysis code
!
!******************************************************************************!
!
! This module is for variables needed to read in the kira output file.
   MODULE sl_input_module
!
     INTEGER :: numstars         ! Number of stars
     INTEGER :: nsnaps           ! Number of snapshots
     INTEGER :: nmax             ! Max number of 'particles' in snapshot (guess)
     INTEGER :: npart            ! Actual number of 'particles' in snapshot
     INTEGER, DIMENSION(:), ALLOCATABLE :: multN ! Num stars in each particle
     CHARACTER*40 :: descriptor  ! Name of value being read in
     CHARACTER*4 :: equals       ! The '='
     CHARACTER*20 :: partname    ! Particle name (e.g. 'root', or '8')
     
! top_info:
     INTEGER :: Ntop     ! The number of stars N in each snap from top info
     DOUBLE PRECISION :: system_time ! System_time for snapshot from top info
     DOUBLE PRECISION :: topm                ! Total system mass
     DOUBLE PRECISION :: topr1,topr2,topr3   ! Total system lengths
     DOUBLE PRECISION :: topv1,topv2,topv3   ! Total system velocities
     DOUBLE PRECISION :: topt                ! Total system time
     DOUBLE PRECISION :: top_energy   ! Total_energy value for the snapshot
     
! multiple_info:
     INTEGER, DIMENSION(:), ALLOCATABLE :: nMult ! Number of multiples in snap
     
! star info:
     INTEGER, DIMENSION(:), ALLOCATABLE :: nstars ! Number of stars in snap
!
     END MODULE sl_input_module
!
!******************************************************************************!
!
! Another module for general use which has the arrays properties of
! the individual stars per snapshot, in non-N-body units, hopefully making it
! a bit easier to use by breaking the program up a bit more.
     MODULE parameters_module
!====================
! Stellar parameters
!====================
! munit,runit,vunit,tunit = conversion factors between N-Body units and
! solar masses, pc, km/s and Myr
       DOUBLE PRECISION :: munit,runit,vunit,tunit
! tmax = time in Myr of each particle in each snapshot
! mmax = mass in solar masses of each particle in each snapshot
! rmax = position in pc of each particle in each snapshot
! vmax = velocity in km/s of each particle in each snapshot
       INTEGER, DIMENSION(:,:), ALLOCATABLE :: idmax  ! IDs of particles
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tmax,mmax
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: rmax,vmax
! As above but for each star in each snapshot (after multiples analysed)
       INTEGER, DIMENSION(:,:), ALLOCATABLE :: idstar
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tstar,mstar
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: rstar,vstar

!====================
! Cluster parameters
!====================
       DOUBLE PRECISION :: totalmass    ! Total mass of the distribution
! com_cluster = the centre of mass of the cluster (x, y, z positions)
! ri_com = distance between each star and cluster centre of mass
       double precision, dimension(:,:),allocatable :: com_cluster
       double precision, dimension(:,:,:),allocatable :: ri_com
! r_halfmass = the half mass radius of the distribution
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r_halfmass
       logical, dimension(:,:), allocatable :: incluster  ! Star escaped if F
       CHARACTER*20 :: limittype    ! Condition in which stars are in cluster.
!           Currently either field of view (=FoV) or half-mass radius (=rhalf)
! Fov_lim = the limit of field of view in pc
! rfac = the factor by which half-mass radius is multiplied
! to establish the boundary of stars that are still in  the cluster
       integer :: FoV_lim, rfac
! save the final value of half-mass radius calculated from 
! whole cluster for use in the rfac*rhalf cluster calculation
       double precision, dimension(4) :: rhalf_all
       character*2 :: thisproj ! Projection of cluster (xy/yz/xz/3D)
       integer :: projnum      ! Integer representing projection type

!===============
! Calculations
!===============
! Total kinetic & potential energy, and total overall energy:
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: kinetic_energy
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: potential_energy
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: total_energy

!============== Lambda ===============
! Lengths of 'object' edges in MST for CDFs
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: edgelengths
! Number of median values to include in l_Nmed (e.g. 2 or 3 for even/odd nedge)
       INTEGER :: Nmed
       INTEGER :: nCDF       ! Number of CDF plots of random MSTs
! lambda = measure of mass segregation & errors
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda, l_up, l_low
! mean MST length of nloop random MSTs, object MST length
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: l_avranmst, l_objmst
       double precision, dimension(:,:), allocatable :: l_allmsts
! lambda_bar uses mean length of MST (basically same as lambda)
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_bar, l_up_bar, l_low_bar
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lbar_avranmst, lbar_objmst
       double precision, dimension(:,:), allocatable :: lbar_allmsts
! lambda_rms uses root mean square length of MST (generalised mean with p=2)
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_rms, l_up_rms, l_low_rms
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lrms_avranmst, lrms_objmst
       double precision, dimension(:,:), allocatable :: lrms_allmsts
! lambda_smr uses square mean root length of MST (generalised mean with p=1/2)
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_smr, l_up_smr, l_low_smr
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lsmr_avranmst, lsmr_objmst
       double precision, dimension(:,:), allocatable :: lsmr_allmsts
! lambda_har uses harmonic mean (generalised mean with p=-1)
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_har, l_up_har, l_low_har
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lhar_avranmst, lhar_objmst
       double precision, dimension(:,:), allocatable :: lhar_allmsts
! lambda_tilde uses median length of MST
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_til, l_up_til, l_low_til
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ltil_avranmst, ltil_objmst
       double precision, dimension(:,:), allocatable :: ltil_allmsts
! lambda_Nmed uses the mean of the N median lengths of MST (Nmed)
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_Nmed, l_up_Nmed, l_low_Nmed
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lNmed_avranmst, lNmed_objmst
       double precision, dimension(:,:), allocatable :: lNmed_allmsts
! lambda_star uses median length of MST and adds this to the actual length of the MST
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_star, l_up_star, l_low_star
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lstar_avranmst, lstar_objmst
       double precision, dimension(:,:), allocatable :: lstar_allmsts
! lambda_gam uses the geometric mean
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_gam, l_up_gam, l_low_gam
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lgam_avranmst, lgam_objmst
       double precision, dimension(:,:), allocatable :: lgam_allmsts
! lamda_ln: sum the exponents and take the log (sort of inverse of geometric mean)
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda_ln, l_up_ln, l_low_ln
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lln_avranmst, lln_objmst
       double precision, dimension(:,:), allocatable :: lln_allmsts
       
! state which types of lambda you want to find, to save doing all every time
       LOGICAL :: findlam, findlambar, findlamrms, findlamsmr, findlamhar
       LOGICAL :: findlamtil, findlamNmed, findlamstar, findgam, findlamln
       
       integer :: nmst          ! number of stars in the minimum spanning tree
       INTEGER :: nloop         ! number of random MSTs calculated in loop
       double precision, dimension(:,:), allocatable :: obj_mass
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: time

!===============
! Outputs
!===============
!
! unit# = keep track of file units by setting as variables at top of file
       integer :: fileunit, sepunit1, sepunit2
       integer :: objmunit, cdfobjunit, cdfranunit
       CHARACTER*150 :: outarg  ! Destination directory (e.g. 'outputs')
       CHARACTER*150 :: newpath ! Output path with cluster type appended
!                                 (e.g. all, FoV)

     END MODULE parameters_module
!
!******************************************************************************!
!
     MODULE constants_module
       REAL :: pi,twopi
       DOUBLE PRECISION :: au=1.5d11                   ! m
       DOUBLE PRECISION :: rsun=6.955d8                ! m
       DOUBLE PRECISION :: msun=2.d30                  ! kg
       DOUBLE PRECISION :: pc=3.086d16                 ! m
       DOUBLE PRECISION :: yr=365.25*24.*60.*60.       ! s
       DOUBLE PRECISION :: day=24*60*60                ! s
       DOUBLE PRECISION :: Myr=365.25*24.*60.*60.*1.d6 ! s
       DOUBLE PRECISION :: G=6.673d-11                 ! m^3 kg^-1 s^-2
       DOUBLE PRECISION :: mjup=1.8986d27              ! kg
     END MODULE constants_module
