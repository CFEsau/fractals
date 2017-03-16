!******************************************************************************!
!
! 08/12/14 Dan Griffiths
! 
! reduce.f90
!
!******************************************************************************!
PROGRAM reduce

!
!TODO: parallelise!
!
! Declare modules used
  USE sl_input_module
  USE parameters_module
  USE constants_module
  IMPLICIT NONE
  
! Declare variables:
! iterators
  INTEGER :: i,j,k,l
! ignore certain stars?
  LOGICAL :: ignore
! Name and path of runfile
  CHARACTER*150 :: inarg
  CHARACTER*4 :: ofilen
  CHARACTER*8 :: outfile
  LOGICAL :: writesnap !write data to snapshots? T/F

  writesnap=.FALSE.
  
  
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
! Make a directory for this simulation data
  CALL SYSTEM('mkdir -p '//TRIM(outarg))

!
!******************************************************************************!
!
! Write out data from snapshots. One file for each snapshot.
!
  IF (writesnap) THEN
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
           WRITE(4,104) j, ids(i,j),t(i,j),m(i,j),r(i,j,1:3),v(i,j,1:3)
        END DO
104     FORMAT (2(2X,I4),2X,F8.4,2X,F7.3,6(2X,F8.3))
        CLOSE(4)
     END DO
  END IF
!
!
!******************************************************************************!
!
! Allocate arrays for calculations.

! logical array: is the star in the cluster
  ALLOCATE(incluster(1:snapnum,1:nstars(1)))

! x, y, z positions of centre of mass:
  ALLOCATE(com_cluster(1:snapnum,1:3))

! x, y, z distances of each star from com:
  ALLOCATE(ri_com(1:snapnum,1:nstars(1),1:3))

! halfmass radius of cluster:
  ALLOCATE(r_halfmass(1:snapnum))

! Energy of each snapshot
  ALLOCATE(kinetic_energy(1:snapnum))
  ALLOCATE(potential_energy(1:snapnum))
  ALLOCATE(total_energy(1:snapnum))


! Define the number of stars in the MST:
  nmst=10

! Find energy, c of m, rhalf, mass segregation for all stars in cluster:
!  CALL reduce_cluster(nstars(1))

! Find energy, c of m, rhalf, mass segregation for stars
! within a field of view of FoV_lim pc:

  FoV_lim = 5
  CALL reduce_FoV(nstars(1))

! Find energy, c of m, rhalf, mass segregation for stars
! within rfac*r_halfmass(snapnum,1:4) of all stars:

!  rfac = 2
!  CALL reduce_rhalf(nstars(1))

!  rfac = 3
!  CALL reduce_rhalf(nstars(1))
!
!  
!*************************!
!
! Deallocate global arrays
  CALL DEALLOCATE
  
END PROGRAM reduce


!******************************************************************************!

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

  DEALLOCATE(incluster)

  DEALLOCATE(com_cluster)
  DEALLOCATE(ri_com)

  DEALLOCATE(r_halfmass)

  DEALLOCATE(kinetic_energy)
  DEALLOCATE(potential_energy)
  DEALLOCATE(total_energy)

END SUBROUTINE DEALLOCATE


!CHARACTER(len=20) FUNCTION str(k)
! Converts an integer to a string
!  INTEGER :: k
!  WRITE(str,*) k
!END FUNCTION str
