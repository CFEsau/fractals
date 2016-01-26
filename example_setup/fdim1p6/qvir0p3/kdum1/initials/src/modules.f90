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
! infilename = the name of the input file
     CHARACTER*25 :: infilename
! nsys = the number of systems (Not individual stars)
     INTEGER :: nsys
! distrib = the mass distribution (f=fractal, p=plummer)
! fdim = the fractal dimension
     CHARACTER :: distrib
     REAL :: fdim
! tend = The simulation length in Myr
! tout = The time interval between snapshots in Myr
     REAL :: tend, tout
! filestem = the general filestem for one batch of runs
     CHARACTER*6 :: filestem
! runit = The conversion factor for position between pc and N-body units
     REAL :: runit
!
   END MODULE input_module
!
!******************************************************************************!
!
! These are all variables that characterise the properties of the distribution within the
! simulation, the number of stars and the individual positions and velocities of those stars.
! These are all properties that are calculated within initials.f90.
!
   MODULE properties_module
!
! nsim = the number of simulations
     INTEGER :: nsim
! nstars = the total number of stars
! numsingle = the number of stars in single systems
! nmax = a guess as to the maximum possible number of stars (used to allocate arrays)
! rsys = The system positions
! vsys = the system velocities
! msys = the system masses
! sinfo = the number of stars in a system
     INTEGER :: nstars, nmax, numsingle
     REAL, DIMENSION(:, :), ALLOCATABLE :: rsys,vsys
     REAL, DIMENSION(:), ALLOCATABLE :: msys
     INTEGER, DIMENSION(:), ALLOCATABLE :: sinfo   
! r = the stellar positions (for single-star systems, this will be the same as rsys)
! v = the stellar velocities
! m = the stellar masses
     REAL, DIMENSION(:, :), ALLOCATABLE :: r,v
     REAL, DIMENSION(:), ALLOCATABLE :: m
! com, cov = The centres of mass and velocity of the distribution
     DOUBLE PRECISION :: com(1:3),cov(1:3)
! qvir = the virial ratio of the distribution.
     REAL :: qvir
! Binaries:
! numbinary = the total number of binary systems
     INTEGER :: numbinary
! fbinary = the fraction of systems that are binaries
     REAL :: fbinary
! pairing = the method by which binaries are paired
! e.g. both masses randomly selected from IMF,
! or mass of secondry is a given fraction of mass of primary
     CHARACTER :: pairing
! fmult = the multiple fraction parameter, used so see if a system is a binary
! asep = the binary separation (double precision)
! rsep = the binary separation (real)
! P = the orbital period of a binary system
! ecc = the orbital eccentricity of a binary system
! mtotsys = the total mass of all systems (double precision)
! mtotstars =  the total mass of all systems (double precision)
     REAL :: fmult, rsep, P, ecc
     DOUBLE PRECISION :: asep, mtotsys, mtotstars
! ekin = total kinetic energy of all systems (double precision)
! epoti = potential energy of a system (double precision)
! epot = total potential energy of all systems (double precision)
! rkin = total kinetic energy of all systems (real)
! rpot = total potential energy of all systems (real)
! rij = the distance between two systems
     DOUBLE PRECISION :: ekin, epoti, epot
     REAL :: rkin, rpot, rij
! cenc = close encounter distance in au
     REAL :: cenc
! let's have some variables for the Maschberger IMF in here (from Maschberger 2011)
! alpha,beta = values given in the paper
! mu = scale factor
! mUpper,mLower = upper and lower mass limits
     REAL :: alpha,beta,mu,mLower,mUpper
!
   END MODULE properties_module
!
!******************************************************************************!
!
! These are variables that are used to convert variables into proper outputs for 
! kira and to create the output files.
!
   MODULE output_module
!
! ntend = The simulation length in N-body units
! ntout The snapshot interval in N-body units
     REAL :: ntend, ntout
! ntendchar = ntend converted into a string
! ntoutchar = ntout converted into a string
     CHARACTER*10 :: ntendchar, ntoutchar
! ncenc = close encounter distance in N-body units
! ncencchar = ncenc converted into a string
     REAL :: ncenc
     CHARACTER*10 :: ncencchar
! runit = The conversion factor for position between pc and N-body units
! munit = The conversion factor for mass between Msun and N-body units
! tunit = The conversion factor for time between Myr and N-body units
! vunit = The conversion factor for velocity between km/s and N-body units
     REAL :: munit, tunit, vunit
! kiracomm = the kira command line for the runstart script
     CHARACTER*240 :: kiracomm
! File outputs:
! id = the run ID number
! outfilestem/outfilename = The output file that contains the overall properties
! of the simulation
! icfilestem/icfilename = The output file that contains the stellar positions and
! velocities, for input into KIRA.
! runfilestem/runfilename = The file that will contain KIRA output.
! restartfilestem/restartfilename = 
! scriptname = 
     CHARACTER*2 :: id
     CHARACTER*11 :: outfilestem
     CHARACTER*15 :: outfilename
     CHARACTER*10 :: icfilestem
     CHARACTER*13 :: icfilename
     CHARACTER*11 :: runfilestem
     CHARACTER*14 :: runfilename
     CHARACTER*11 :: restartfilestem
     CHARACTER*14 :: restartfilename
     CHARACTER*12 :: scriptname
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
