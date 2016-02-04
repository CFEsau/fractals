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
  ! Find distance between each star and cluster centre of mass.
  !r_com(:,1:3) gives x, y, z from com, and col 4 gives magnitude.
  ALLOCATE(r_com(1:snapnum,1:nstars(1),1:4))
  r_com=0.
! Loop over all snapshots
! Calculate com in each case and populate the array
  DO i=1, snapnum
     CALL c_of_m(i,nstars(i))
  END DO
  
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
!
! Find the degree of mass segregation (lambda).
  ALLOCATE(lambda(1:snapnum))
  ALLOCATE(l_low(1:snapnum))
  ALLOCATE(l_up(1:snapnum))
  lambda=0.
  l_up=0.
  l_low=0.
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_lambda(i,nstars(i))
     if ((i==1).or.(i==900)) then
        !PRINT *, 'Lambda:', i, lambda(i)
     end if
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
     WRITE(4,*) i,kinetic_energy(i),potential_energy(i),total_energy(i),r_halfmass(i),lambda(i),l_low(i),l_up(i)
  END DO
  CLOSE(4)
!
!  
!******************************************************************************!
!
! Deallocate global arrays
  CALL DEALLOCATE
END PROGRAM reduce


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

  DEALLOCATE(r_com)
  
  DEALLOCATE(kinetic_energy)
  DEALLOCATE(potential_energy)
  DEALLOCATE(total_energy)

  DEALLOCATE(r_halfmass)
  
  DEALLOCATE(lambda)
  DEALLOCATE(l_up)
  DEALLOCATE(l_low)
  
END SUBROUTINE DEALLOCATE

