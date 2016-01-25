!******************************************************************************!
!
! 08/12/14 Dan Griffiths
! 
! reduce.f90
!
!******************************************************************************!
PROGRAM reduce

! Declare modules used
  USE sl_input_module
  USE parameters_module
  USE constants_module
  IMPLICIT NONE
! Declare variables:
! iterators
  INTEGER :: i,j,k,l
! Name and path of runfile
  CHARACTER*150 :: inarg, outarg
  CHARACTER*4 :: ofilen
  CHARACTER*8 :: outfile
! Get the name and path of the runfile from the command line
  CALL GETARG(1,inarg)
  PRINT *, 'Input file: ',TRIM(inarg)

  CALL GETARG(2,outarg)
  PRINT *, 'Output directory: ',TRIM(outarg)
!
! read_sl reads in the sl file, converting the relevant data into arrays.
! These include the overall properties of each individual snapshot and, for each snapshot,
! the properties (mass, position, velocity) of the individual stars.
!
  CALL read_sl_out(inarg)
!
! Next step is to convert these masses, positions and velocities from N-body units to usable units
! (Solar masses, au and km/s?). Can do this here, don't need a separate subroutine.
! So let's get runit, munit, and tunit (NOT the same as runit, munit and tunit in initials.f90):
  munit=mass_scale
  runit=size_scale
  tunit=time_scale
!
! Now let's just multiply the masses, positions, velocities and times by these values
! to get m in solar masses, r in pc, v in km/s and t in Myr
  star_t=star_t/tunit
  star_m=star_m/munit
  star_r=star_r*(2.255e-8/runit)  
  star_v=star_v*(2.255e-8/runit)*tunit
  ! Convert from pc/Myr into km/s
  star_v=star_v*pc*1.e-9/yr
! These data are Reals at the moment. I want Double Precision.
! First allocate arrays
  ALLOCATE(t(1:snapnum,1:nstars(1)))
  ALLOCATE(m(1:snapnum,1:nstars(1)))
  ALLOCATE(r(1:snapnum,1:nstars(1),1:3))
  ALLOCATE(v(1:snapnum,1:nstars(1),1:3))
! Now fill the arrays
  t=DBLE(star_t)
  m=DBLE(star_m)
  r=DBLE(star_r)
  v=DBLE(star_v)
!
!******************************************************************************!
!
! Find total kinetic, gravitational potential and total energy.
  ALLOCATE(kinetic_energy(1:snapnum))
  ALLOCATE(potential_energy(1:snapnum))
  ALLOCATE(total_energy(1:snapnum))
! Loop over all snapshots
  DO i=1, snapnum
     CALL find_energy(i,nstars(i))
!!$     PRINT *, kinetic_energy(i),potential_energy(i),total_energy(i)
  END DO
!
!******************************************************************************!
!
! Find the Half mass radius.
  ALLOCATE(r_halfmass(1:snapnum))
  r_halfmass=0.
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_halfmass(i,nstars(i))
!!$     PRINT *, i, r_halfmass(i)
  END DO
!
!******************************************************************************!
! Make a directory for this simulation data
  CALL SYSTEM('mkdir -p '//TRIM(outarg))
!
! Write out data from snapshots. One file for each snapshot.
!
  CALL SYSTEM('mkdir -p '//TRIM(outarg)//'/snapshots')
  DO i=1,snapnum
     IF (i<10) THEN
        WRITE(ofilen,'(i1)')i
        outfile='snap'//'000'//ofilen
     ELSE IF (i<100) THEN
        WRITE(ofilen,'(i2)')i
        outfile='snap'//'00'//ofilen
     ELSE IF (i<1000) THEN
        WRITE(ofilen,'(i3)')i
        outfile='snap'//'0'//ofilen
     ELSE
        WRITE(ofilen,'(i4)')i
        outfile='snap'//ofilen
     END IF
     OPEN(4,file=TRIM(outarg)//'/snapshots/'//outfile,status='new')
     DO j=1,nstars(i)
        WRITE(4,*) j, ids(i,j),t(i,j),m(i,j),r(i,j,1:3),v(i,j,1:3)
     END DO
     CLOSE(4)
  END DO
!
! Write out the half mass radius and energy data
  OPEN(4,file=TRIM(outarg)//'/macro',status='new')
  DO i=1,snapnum
     WRITE(4,*) i,kinetic_energy(i),potential_energy(i),total_energy(i),r_halfmass(i)
  END DO
  CLOSE(4)
!
!  
!******************************************************************************!
!
! Deallocate global arrays
  CALL DEALLOCATE
END PROGRAM reduce

SUBROUTINE find_energy(snapshoti,ni)
  ! Find total kinetic, gravitational potential and total energy.
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri, vi
! ke = the total kinetic energy
! rij = the magnitude of the separation between two stars
! epoti = the gravitational potential energy between two stars
! epot = the total gravitational potential energy
! etot = the total energy
  DOUBLE PRECISION :: ekin,rij,epoti,epot,etot
  INTEGER :: i,j,k
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(vi(1:ni,1:3))
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  vi(1:ni,1:3)=v(snapshoti,1:ni,1:3)
! First convert to SI units
  mi=mi*msun
  ri=ri*pc
  vi=vi*1000
  ekin=0.
! Loop over all stars
  DO i=1,ni
! Find kinetic energy
     ekin=ekin+(0.5*mi(i)*(vi(i,1)**2+vi(i,2)**2+vi(i,3)**2))
  END DO
  epot=0.
  rij=0.
! Loop over all stars
  DO i=1,ni-1
     epoti=0.d0
     DO j=i+1,ni
! Find separation between two stars
        rij=(ri(i,1)-ri(j,1))**2 + (ri(i,2)-ri(j,2))**2 + & 
             &              (ri(i,3)-ri(j,3))**2
! Find gravitational potential energy between two stars
        epoti=epoti + DBLE(mi(j)/SQRT(rij))
     END DO
! Find total gravitational potential energy
     epot=epot + G*DBLE(mi(i))*epoti
  END DO
  etot=ekin-epot
! add energies to array
  kinetic_energy(snapshoti)=ekin
  potential_energy(snapshoti)=-epot
  total_energy(snapshoti)=etot
! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(vi)
!!$  Print *, ekin, -epot, etot
END SUBROUTINE find_energy

SUBROUTINE find_halfmass(snapshoti,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
! com_dist = the centre of mass of the distribution
! totalmass = the total mass of the distribution
! rmag = array with the magnitude of the separation of all stars from the distribution com
! rlist = a list of all stars in the distribution
  DOUBLE PRECISION, DIMENSION(3) :: com_dist
  DOUBLE PRECISION :: totalmass, massi, massj
  REAL, DIMENSION(:), ALLOCATABLE :: rmag
  INTEGER, DIMENSION(:), ALLOCATABLE :: rlist
  INTEGER :: i,j,k
! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(rmag(1:ni))
  ALLOCATE(rlist(1:ni))
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
! First find the com of the cluster.
! initialise variables
  totalmass=0.
  com_dist=0.
  rlist=0
  rmag=0.
! Loop over all stars
  DO i=1,ni
     totalmass=totalmass+mi(i)
     com_dist(1:3)=com_dist(1:3)+(mi(i)*ri(i,1:3))
  END DO
  com_dist(1:3)=com_dist(1:3)/totalmass
!
! Next find the distance between each star and the com of the cluster
  DO i=1,ni
     rlist(i)=i
     rmag(i)=REAL(SQRT((ri(i,1)-com_dist(1))**2+(ri(i,2)-com_dist(2))**2+(ri(i,3)-com_dist(3))**2))
  END DO
!
! Sort the stars by distance from the com
  CALL heapsort(ni,rmag,rlist)
! The half mass radius is then the radius at which half of the stelalr mass in inside and half is outside.
! So basically keep adding up the masses of the stars, starting from the centre moving out, until
! half of the mass is inside the radius. There'll be a radius ri-1 where M < totalM/2 and ri > totalM/2.
! Use linear interpolation of these two values to find the halfmass radius.
  i=0
  j=0
  massi=0.
  massj=0.
  DO WHILE(massi<(0.5*totalmass))
     massj=massi
     j=i
     i=i+1
     massi=massi+mi(rlist(i))
  END DO
  ! Linear interpolation to find the halfmass radius
  r_halfmass(snapshoti)=rmag(rlist(j))+((rmag(rlist(i))-rmag(rlist(j)))*(((0.5*totalmass)-massj)/(massi-massj)))
!
  ! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(rmag)
  DEALLOCATE(rlist)
END SUBROUTINE find_halfmass

SUBROUTINE DEALLOCATE
  USE sl_input_module
  USE parameters_module
    
  DEALLOCATE(nMult)
  DEALLOCATE(nstars)
! Multiple information
  DEALLOCATE(mult_nstars)
  DEALLOCATE(mult_ids)
  DEALLOCATE(mult_t)
  DEALLOCATE(mult_m)
  DEALLOCATE(mult_r)
  DEALLOCATE(mult_v)
  
! Stellar information
  DEALLOCATE(ids)
  DEALLOCATE(star_t)
  DEALLOCATE(star_m)
  DEALLOCATE(star_r)
  DEALLOCATE(star_v)

  DEALLOCATE(t)
  DEALLOCATE(m)
  DEALLOCATE(r)
  DEALLOCATE(v)

  DEALLOCATE(kinetic_energy)
  DEALLOCATE(potential_energy)
  DEALLOCATE(total_energy)

  DEALLOCATE(r_halfmass)
  
END SUBROUTINE DEALLOCATE

