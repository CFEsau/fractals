!******************************************************************************!
!
! Modules for initial conditions code
!
!******************************************************************************!
!
! These are variables that are input into initials.f90 from an input file.
!
   MODULE input_module
!
     CHARACTER*25 :: infilename !the name of the input file
     INTEGER :: nsys            !number of systems (not individual stars)
     CHARACTER :: distrib       !mass distribution (f=fractal, p=plummer)
     REAL :: fdim               !the fractal dimension
     REAL :: tend               !the simulation length in Myr
     REAL :: tout               !the time interval between snapshots in Myr
     CHARACTER*5 :: filestem    !the general filestem for one batch of runs
     REAL :: runit  !conversion factor for position between pc & N-body units
!
   END MODULE input_module
!
!******************************************************************************!
!
! These are all variables that characterise the properties of the distribution
! within the simulation, the number of stars and the individual positions and
! velocities of those stars.
! These are all properties that are calculated within initials.f90.
!
   MODULE properties_module
!
     INTEGER :: nsim       !the number of simulations
     INTEGER :: nstars     !the total number of stars
     INTEGER :: nmax       !a guess at the max possible number of stars
!                          !(used to allocate arrays)
     INTEGER :: numsingle  !number of stars in single systems
     REAL, DIMENSION(:, :), ALLOCATABLE :: rsys  !the system positions
     REAL, DIMENSION(:, :), ALLOCATABLE :: vsys  !the system velocities
     REAL, DIMENSION(:), ALLOCATABLE :: msys     !the system masses
     INTEGER, DIMENSION(:), ALLOCATABLE :: sinfo !number of stars in a system
     REAL, DIMENSION(:, :), ALLOCATABLE :: r     !stellar positions
!                                                !(for single-star systems this
!                                                ! will be the same as rsys)
     REAL, DIMENSION(:, :), ALLOCATABLE :: v     !stellar velocities
     REAL, DIMENSION(:), ALLOCATABLE :: m        !stellar masses
! com, cov = The centres of mass and velocity of the distribution
     DOUBLE PRECISION :: com(1:3),cov(1:3) !the centres of mass & velocity
!                                          !of the distributions
     REAL :: qvir           !virial ratio of the distribution
! Binaries:
     INTEGER :: numbinary   !the total number of binary systems
     REAL :: fbinary        !the fraction of systems that are in binaries
     CHARACTER*5  :: pairing   !the method by which binaries are paired*
!                              !(e.g. both masses randomly selected from IMF or
                               !mass of 2ndary is fraction of mass of primary)
!* TODO: need to vary length of string: 'ratio' (*5) or 'IMF' (*3)
     REAL :: fmult                 !the multiple fraction parameter.
!                                  !Used to see whether a system is a binary
     REAL :: rsep                  !the binary separation (real)
     REAL :: P                     !the orbital period of a binary system
     REAL :: ecc                   !the orbital eccentricity of a binary system
     DOUBLE PRECISION :: asep      !the binary separation (double precision)
     DOUBLE PRECISION :: mtotsys   !total mass of all systems (double precision)
     DOUBLE PRECISION :: mtotstars !total mass of all stars (double precision)
     DOUBLE PRECISION :: ekin  !total kinetic energy of all systems (dble)
     DOUBLE PRECISION :: epoti !total potential energy of a system (dble)
     DOUBLE PRECISION :: epot  !total potential energy of all systems (dble)
     REAL :: rkin              !total kinetic energy of all systems (real)
     REAL :: rpot              !total potential energy of all systems (real)
     REAL :: rij               !the distance between two systems
     REAL :: cenc              !close encounted distance in au
! Let's have some variables for the Maschberger IMF in here
! (from Maschberger 2011)
     REAL :: alpha,beta       !values given in the paper
     REAL :: mu,mUpper,mLower !scale factor, and upper & lower mass limits
!
   END MODULE properties_module
!
!******************************************************************************!
!
! These are variables that are used to convert variables into
! proper outputs for kira and to create the output files.
!
   MODULE output_module
!
     REAL :: ntend, ntout !simulation length & snapshot interval (N-body units)
     CHARACTER*10 :: ntendchar, ntoutchar !ntend & ntout converted into strings
     REAL :: ncenc                  !close encounter distance in N-body units
     CHARACTER*10 :: ncencchar      !ncenc converted into a string
     REAL :: munit  !conversion factor for position between pc & N-body units
     REAL :: tunit  !conversion factor for time between Myr & N-body units
     REAL :: vunit  !conversion factor for velocity between km/s & N-body units
     CHARACTER*240 :: kiracomm !the kira command line for the runstart script
! File outputs:
! restartfilestem/restartfilename = 
! scriptname = 
     CHARACTER*2 :: id            !the run ID number
     CHARACTER*10 :: outfilestem  !the output file that contains the
     CHARACTER*14 :: outfilename  !       overall properties of the simulation
     CHARACTER*9 :: icfilestem    !the output file that contains the stellar
     CHARACTER*12 :: icfilename   !  positions & velocities for input into KIRA
     CHARACTER*10 :: runfilestem  !the file that will contain
     CHARACTER*13 :: runfilename  !                           KIRA output      
     CHARACTER*10 :: restartfilestem 
     CHARACTER*13 :: restartfilename
     CHARACTER*11 :: scriptname
!
   END MODULE output_module
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
     DOUBLE PRECISION :: Myr=365.25*24.*60.*60.*1.d6 ! s
     DOUBLE PRECISION :: G=6.673d-11                 ! m^3 kg^-1 s^-2
     DOUBLE PRECISION :: mjup=1.8986d27              ! kg
   END MODULE constants_module
