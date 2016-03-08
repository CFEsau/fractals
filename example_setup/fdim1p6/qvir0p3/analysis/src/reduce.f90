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
! ignore certain stars?
  logical :: ignore
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
! Make a directory for this simulation data
  CALL SYSTEM('mkdir -p '//TRIM(outarg))
 
 
!
!******************************************************************************!
!
! Find distance between each star and cluster centre of mass.
! r_com(:,1:3) gives x, y, z from com.
! col 4 gives 2D distance magnitude & 5 gives 3D.
  ALLOCATE(r_com(1:snapnum,1:nstars(1),1:5))
  r_com=0.

! Loop over all snapshots
! Calculate com in each case and populate the array
  DO i=1, snapnum
     CALL c_of_m(i,nstars(i))
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
! If stars beyond a certain distance should be ignore from following calculations,
! call the 'in_cluster' subroutine to populate a 'logical' array,
! setting 'true' for stars in cluster and 'false' if it has escaped.
  
  allocate(this_in2Dcluster(1:snapnum,1:nstars(1)))
  allocate(this_in3Dcluster(1:snapnum,1:nstars(1)))
! All stars start off in the cluster
  this_in2Dcluster=.true.
  this_in3Dcluster=.true.

! Should I ignore certain stars? (currently ignores outliers)
  ignore=.true.
  
  if (ignore) then
     open(10,file=TRIM(outarg)//'/2Dignored.txt')
     open(11,file=TRIM(outarg)//'/3Dignored.txt')
     do i=1,snapnum
        call in_cluster(i,nstars(i))
     end do
     close(10)
     close(11)
  end if

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
! Mass segratation
!******************************************************************************!
! Find the degree of mass segregation (lambda).

! Define the number of stars in the MST:
  nmst=10

! IDs of the most massive stars in the cluster:
  allocate(obj_mass(1:2,1:nmst))

!*****************
! Lambda:

  ALLOCATE(lambda(1:snapnum))
  ALLOCATE(l_low(1:snapnum))
  ALLOCATE(l_up(1:snapnum))
  lambda=0.
  l_up=0.
  l_low=0.
  open(20,file=TRIM(outarg)//'/obj_masses.txt')
  open(21,file=TRIM(outarg)//'/escaped.txt')
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_lambda(i,nstars(i))
     !if ((i==1).or.(i==snapnum)) then
        !PRINT *, 'Lambda:', i, lambda(i)
     !end if
  END DO
  close(20)
  close(21)

!
!*****************
! Lambda_bar:

  ALLOCATE(lambda_bar(1:snapnum))
  ALLOCATE(l_low_bar(1:snapnum))
  ALLOCATE(l_up_bar(1:snapnum))
  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.
  !open(20,file=TRIM(outarg)//'/obj_masses_lambar.txt')
  !open(21,file=TRIM(outarg)//'/escaped_lambar.txt')
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_lambda_bar(i,nstars(i))
  END DO
  !close(20)
  !close(21)

!
!*****************
! Lambda_tilde:

  ALLOCATE(lambda_tilde(1:snapnum))
  ALLOCATE(l_low_tilde(1:snapnum))
  ALLOCATE(l_up_tilde(1:snapnum))
  lambda_tilde=0.
  l_up_tilde=0.
  l_low_tilde=0.
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_lambda_tilde(i,nstars(i))
  END DO

!
!*****************
! Lambda_star:

  ALLOCATE(lambda_star(1:snapnum))
  ALLOCATE(l_low_star(1:snapnum))
  ALLOCATE(l_up_star(1:snapnum))
  lambda_star=0.
  l_up_star=0.
  l_low_star=0.
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_lambda_star(i,nstars(i))
  END DO

!
!*****************
! Gamma:

  ALLOCATE(gamm(1:snapnum))
  ALLOCATE(g_low(1:snapnum))
  ALLOCATE(g_up(1:snapnum))
  gamm=0.
  g_up=0.
  g_low=0.
! Loop over all snapshots
  DO i=1,snapnum
     CALL find_gamma(i,nstars(i))
  END DO


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
        WRITE(4,*) j, ids(i,j),t(i,j),m(i,j),r(i,j,1:3),v(i,j,1:3),this_in2Dcluster(i,j),this_in3Dcluster(i,j)
     END DO
     CLOSE(4)
  END DO
!
! Write out the half mass radius and energy data
! Lambda given its own separate file, given there are 5 different calculations
  OPEN(4,file=TRIM(outarg)//'/macro',status='new')
  OPEN(5,file=TRIM(outarg)//'/lambda',status='new')
  DO i=1,snapnum
     WRITE(4,40) i,kinetic_energy(i),potential_energy(i),total_energy(i),r_halfmass(i)
     WRITE(5,50) i,lambda(i),l_low(i),l_up(i),lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_tilde(i),l_low_tilde(i),l_up_tilde(i),lambda_star(i),l_low_star(i),l_up_star(i), &
          & gamm(i),g_low(i),g_up(i)
  END DO
40 FORMAT(1X,I4,3(2X,E9.3),2X,F7.3)
50 FORMAT(1X,I4,15(2X,F8.3))
  CLOSE(4)
  CLOSE(5)
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
  
  deallocate(this_in2Dcluster)
  deallocate(this_in3Dcluster)

  DEALLOCATE(kinetic_energy)
  DEALLOCATE(potential_energy)
  DEALLOCATE(total_energy)
  
  deallocate(obj_mass)
  DEALLOCATE(lambda)
  DEALLOCATE(l_up)
  DEALLOCATE(l_low)
  DEALLOCATE(lambda_bar)
  DEALLOCATE(l_up_bar)
  DEALLOCATE(l_low_bar)
  DEALLOCATE(lambda_tilde)
  DEALLOCATE(l_up_tilde)
  DEALLOCATE(l_low_tilde)
  DEALLOCATE(lambda_star)
  DEALLOCATE(l_up_star)
  DEALLOCATE(l_low_star)
  DEALLOCATE(gamm)
  DEALLOCATE(g_up)
  DEALLOCATE(g_low)
  
END SUBROUTINE DEALLOCATE

