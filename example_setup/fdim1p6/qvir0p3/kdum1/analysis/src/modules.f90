!******************************************************************************!
!
! Modules for analysis code
!
!******************************************************************************!
!
! This module is for variables needed to read in the kira output file.
   MODULE sl_input_module
!
! snapnum = number of snapshots. Starts at zero and is increased incrementally.
       INTEGER :: snapnum
! desciptor = The first part of a line being read in, generally descibing the data
! on the line, e.g. '(Particle' or 'name'
! equals = the '=' sign found on most lines.
! partname = the particle name, e.g. 'root' for topinfo, maybe '(6,100006)' for sysinfo, stars have an ID not a name
     CHARACTER*40 :: descriptor  ! name of value being read in
     CHARACTER*4 :: equals       ! the '='
     CHARACTER*20 :: partname    ! particle name (e.g. 'root', or '8')
! top_info:
! nTop = the number of stars N given in the topinfo
     INTEGER :: nTop
! system_time = the system_time value for the snapshot
! top_energy = the total_energy value for the snapshot
! mass_scale = the mass_scale value for the snapshot
! size_scale = the size_scale value for the snapshot
! time_scale = the time_scale value for the snapshot
     REAL :: system_time, top_energy, mass_scale, size_scale, time_scale
! multiple_info:
! nMult = the number of multiples in the snapshot
     INTEGER, DIMENSION(:), ALLOCATABLE :: nMult
! mult_nstars = number of stars in each  multiple
! mult_ids = ids of stars in each multiple
     INTEGER, DIMENSION (:,:), ALLOCATABLE :: mult_nstars
     INTEGER, DIMENSION (:,:,:), ALLOCATABLE :: mult_ids
! These find the IDs of stars in a binary. What about a multiple?
! parop, parcl = Opening and closing parenthesis
! cstar1, cstar2 = ids of star1 and star2 in character format
! cstar1_len, cstar2_len = length of strings cstar1 and cstar2 in character format (used for formatting purposes)
     CHARACTER*1 :: parop, parcl
     CHARACTER*10, DIMENSION(5) :: cstar, cstar_len
! mult_t = time of multiple
! mult_m = mass of multiple
! mult_r = center of mass of multiple
! mult_v = center of velocity of multiple
! mult_childof = whether the multiple is a child of another multiple (given by ID of first star)
     REAL, DIMENSION(:,:), ALLOCATABLE :: mult_t,mult_m
     REAL, DIMENSION(:,:,:), ALLOCATABLE :: mult_r,mult_v
! star info:
! nstars = the total number of stars in the snapshot
! ids = the ids of the individual stars
     INTEGER, DIMENSION(:), ALLOCATABLE :: nstars
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: ids
! star_t = time of star
! star_m = mass of star
! star_r = position of star
! star_v = velocity of star
     REAL, DIMENSION(:,:), ALLOCATABLE :: star_t,star_m
     REAL, DIMENSION(:,:,:), ALLOCATABLE :: star_r,star_v
!
     END MODULE sl_input_module
!
!******************************************************************************!
!
! I'm going to create another module for general use which has the arrays properties of
! the individual stars per snapshot, in non-N-body units, hopefully making it a bit
! easier to use by breaking the program up a bit more. Plus we want double precision.
     MODULE parameters_module
!
! munit,runit,vunit,tunit = conversion factors between N-Body units and
! solar masses,pc,km/s and Myr
       DOUBLE PRECISION :: munit,runit,vunit,tunit
! t = time in Myr of each star in each snapshot
! m = mass in solar masses of each star in each snapshot
! r = position in pc of each star in each snapshot
! v = velocity in km/s of each star in each snapshot
       DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: t,m
       DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: r,v
! totalmass = the total mass of the distribution
! r_com = distance between each star and cluster centre of mass
       DOUBLE PRECISION :: totalmass
       double precision, dimension(:,:,:),allocatable :: r_com
! These next variables are for finding the kinetic, potential and total energy:
! kinetic_energy = the total kinetic energy
! potential_energy = the total gravitational potential energy
! total_energy = the total energy
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: kinetic_energy,potential_energy,total_energy
! r_halfmass = the half mass radius of the distribution
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: r_halfmass
! lambda = measure of mass segregation & errors
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda, l_up, l_low
! Time
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: time
!
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
