!******************************************************************************!
!
! Star Wars Day 2017 - Claire Esau
! Bastardisation of code from Dan Griffiths 08/12/14
! and Richard Allison from aeons ago. Also added a bunch of new stuff.
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
  
  CALL GETARG(1,inarg) ! Get name & path of runfile from command line
  PRINT *, 'Input file: ',TRIM(inarg)

  CALL GETARG(2,outarg) ! Get path of output from command line
  PRINT *, 'Output directory: ',TRIM(outarg)
!
  CALL read_sl_out(inarg)
! Reads in the sl file & saves required data in arrays.
! Data saved: the overall properties of each individual snapshot and
! properties (mass, position, velocity) of the particles in each snapshot.
!
! Convert these masses, positions & velocities from N-body units
! to t in Myr, m in solar masses, r in pc, v in km/s
  tstar=tstar/tunit
  mstar=mstar/munit
  rstar=rstar*(2.255e-8/runit)
  vstar=vstar*(2.255e-8/runit)*tunit   ! N-body to pc/Myr
  vstar=vstar*pc*1.e-9/yr              ! pc/Myr to km/s
!
! Make a directory for this simulation data
  CALL SYSTEM('mkdir -p '//TRIM(outarg))
!
!
!******************************************************************************!
!
! Write out data from snapshots. One file for each snapshot.
!
  IF (writesnap) THEN
     CALL SYSTEM('mkdir -p '//TRIM(outarg)//'/snapshots')
     DO j=1,nsnaps
        IF (j<10) THEN
           WRITE(ofilen,'(i1)')j
           outfile='snap'//'000'//ofilen
        ELSE IF (j<100) THEN
           WRITE(ofilen,'(i2)')j
           outfile='snap'//'00'//ofilen
        ELSE IF (j<1000) THEN
           WRITE(ofilen,'(i3)')j
           outfile='snap'//'0'//ofilen
        ELSE
           WRITE(ofilen,'(i4)')j
           outfile='snap'//ofilen
        END IF
        OPEN(4,file=TRIM(outarg)//'/snapshots/'//outfile,status='replace')
        DO i=1,nstars(j)
           WRITE(4,104) i,idstar(i,j),tstar(i,j),mstar(i,j),&
                & rstar(1:3,i,j),vstar(1:3,i,j)
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
  ALLOCATE(incluster(1:nstars(1),1:nsnaps))

! x, y, z positions of centre of mass:
  ALLOCATE(com_cluster(1:3,1:nsnaps))

! x, y, z distances of each star from com:
  ALLOCATE(ri_com(1:3,1:nstars(1),1:nsnaps))

! halfmass radius of cluster:
  ALLOCATE(r_halfmass(1:nsnaps))

! Energy of each snapshot
  ALLOCATE(kinetic_energy(1:nsnaps))
  ALLOCATE(potential_energy(1:nsnaps))
  ALLOCATE(total_energy(1:nsnaps))


! Define the number of stars in the MST:
  nmst=10

! Find energy, c of m, rhalf, mass segregation for all stars in cluster:
  CALL reduce_cluster(nstars(1))

! Find energy, c of m, rhalf, mass segregation for stars
! within a field of view of FoV_lim pc:
  FoV_lim = 5
  CALL reduce_FoV(nstars(1))

! Find energy, c of m, rhalf, mass segregation for stars
! within rfac*r_halfmass(1:4,nsnaps) of all stars:
  rfac = 2
  CALL reduce_rhalf(nstars(1))
  rfac = 3
  CALL reduce_rhalf(nstars(1))
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
  
! Stellar information
  DEALLOCATE(idstar)
  DEALLOCATE(tstar)
  DEALLOCATE(mstar)
  DEALLOCATE(rstar)
  DEALLOCATE(vstar)

  DEALLOCATE(idmax)
  DEALLOCATE(tmax)
  DEALLOCATE(mmax)
  DEALLOCATE(rmax)
  DEALLOCATE(vmax)

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
