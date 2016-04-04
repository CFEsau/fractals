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
  logical :: ignore
! Name and path of runfile
  CHARACTER*150 :: inarg
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
! Make a directory for this simulation data
  CALL SYSTEM('mkdir -p '//TRIM(outarg))
 
 
!
!******************************************************************************!
!
! Find distance between each star and cluster centre of mass.
! r_com(:,1:3) gives x, y, z from com.
! 4:6 gives 2D distance magnitude xy, yz, xz
! column 7 gives 3D distance magnitude
  ALLOCATE(r_com(1:snapnum,1:nstars(1),1:7))
  r_com=0.

! Loop over all snapshots
! Calculate com in each case and populate the array
  write(6,*)"       Calculating centre of mass..."
  DO i=1, snapnum
     CALL c_of_m(i,nstars(i))
  END DO

!
!******************************************************************************!
!
  ! Find the Half mass radius.

! 1:3 gives halfmass radius in xy, yz, and xz
! column 4 gives 3D halfmass radius
  ALLOCATE(r_halfmass(1:snapnum,1:4))
  r_halfmass=0.
  write(6,*)"       Calculating half-mass radius..."
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_halfmass(i,nstars(i))
!!$     PRINT *, i, r_halfmass(i)
  END DO

!
!******************************************************************************!
!
! Find total kinetic, gravitational potential and total energy.
  ALLOCATE(kinetic_energy(1:snapnum))
  ALLOCATE(potential_energy(1:snapnum))
  ALLOCATE(total_energy(1:snapnum))
  write(6,*)"       Calculating cluster energies..."
! Loop over all snapshots
  DO i=1, snapnum
     CALL find_energy(i,nstars(i))
!!$     PRINT *, kinetic_energy(i),potential_energy(i),total_energy(i)
  END DO


!
! Write out the half mass radius and energy data
!
!TODO: calculate different planes for energy and write out all data;
! currently rhalf is calculated for all planes but just written out for 3D
  OPEN(4,file=TRIM(outarg)//'/macro',status='new')
  DO i=1,snapnum
     WRITE(4,40) i,kinetic_energy(i),potential_energy(i),total_energy(i),r_halfmass(i,4)
  END DO
40 FORMAT(1X,I4,3(2X,E9.3),2X,F7.3)
  CLOSE(4)

!
!******************************************************************************!
! Mass segratation
!******************************************************************************!
! Find the degree of mass segregation (lambda).

! Define the number of stars in the MST:
  nmst=10

! All lambda in all planes with all stars in cluster
  call reduce_cluster(i,nstars(i))

! All lambda in all planes with 5 pc field of view
  FoV_lim = 5
  call reduce_FoV(i,nstars(i))

!TODO: change reduce_half so it uses only rhalf at final snapshot
! All lambda in all planes where cluster within  2*r_half
  rfac = 2
  call reduce_rhalf(i,nstars(i))

! All lambda in all planes where cluster within  3*r_half
  rfac = 3
  call reduce_rhalf(i,nstars(i))


!******************************************************************************!
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
! Lambda given its own separate file, given there are 5 different calculations
  !OPEN(4,file=TRIM(outarg)//'/macro',status='new')
  !OPEN(5,file=TRIM(outarg)//'/lambda',status='new')
  !DO i=1,snapnum
     !WRITE(4,40) i,kinetic_energy(i),potential_energy(i),total_energy(i),r_halfmass(i,4)
     !WRITE(5,50) i,lambda(i),l_low(i),l_up(i),lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          !& lambda_tilde(i),l_low_tilde(i),l_up_tilde(i),lambda_star(i),l_low_star(i),l_up_star(i), &
          !& gamm(i),g_low(i),g_up(i)
  !END DO
!40 FORMAT(1X,I4,3(2X,E9.3),2X,F7.3)
!50 FORMAT(1X,I4,15(2X,F8.3))
  !CLOSE(4)
  !CLOSE(5)
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

  DEALLOCATE(r_halfmass)

  DEALLOCATE(kinetic_energy)
  DEALLOCATE(potential_energy)
  DEALLOCATE(total_energy)

END SUBROUTINE DEALLOCATE

