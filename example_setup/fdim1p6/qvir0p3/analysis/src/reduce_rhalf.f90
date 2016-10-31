SUBROUTINE reduce_rhalf(ni)

  USE sl_input_module
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ni
  LOGICAL :: dirExists
  CHARACTER(len=100) :: newDir, rfac_char
  INTEGER :: i,j

  limittype='rhalf'
! write value of rfac to string:
  WRITE(rfac_char,'(I0)') rfac

! create 'rhalf' directory for results:
  newDir = 'cluster_r'//TRIM(rfac_char)//'rhalf'
  newPath = TRIM(outarg)//'/'//TRIM(newDir)

  INQUIRE(file = TRIM(newPath)//'./', exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)

  IF (.NOT. dirExists) THEN
     WRITE(6,'(a)') "Creating new directory: '"//TRIM(newPath)//"'"
     CALL system('mkdir -p '//TRIM(newPath))
  END IF

  WRITE(6,'(a)') "Saving data in '"//TRIM(newPath)//"'..."


  DO projnum = 1,4
     IF (projnum==1) THEN
        proj='xy'
        WRITE(6,*)""
        WRITE(6,*)"   "//proj//":"
     ELSE IF (projnum==2) THEN
        proj='yz'
        WRITE(6,*)""
        WRITE(6,*)"   "//proj//":"
     ELSE IF (projnum==3)  THEN
        proj='xz'
        WRITE(6,*)""
        WRITE(6,*)"   "//proj//":"
     ELSE IF (projnum==4)  THEN
        proj='3D'
        WRITE(6,*)""
        WRITE(6,*)"   "//proj//":"
     END IF

!*************************************************!
! Stars within fac*rhalf
!*************************************************!
!
! Initially all stars are in the cluster; none have escaped
!
     incluster = .TRUE.

     OPEN(10,file=TRIM(newPath)//'/escaped_'//proj//'.dat')

     DO i=1,snapnum
        CALL in_cluster(i,nstars(i))
     END DO

     CLOSE(10)

!
!******************************************************************************!
!
! Find centre of mass of cluster.
     com_cluster=0.

! Find distance between each star and cluster centre of mass.
     ri_com=0.

! Loop over all snapshots
! Calculate com in each case and populate the array
     WRITE(6,*)"       Calculating centre of mass..."
     DO i=1, snapnum
        CALL c_of_m(i,nstars(i))
     END DO

!******************************************************************************!
!
! Find the Half mass radius.

     r_halfmass=0.

     WRITE(6,*)"       Calculating half-mass radius..."
! Loop over all snapshots
     DO i=1,snapnum
        CALL find_halfmass(i,nstars(i))
!!$     PRINT *, i, r_halfmass(i)
     END DO


!*******************************************
! Write out distance data
!
! Centre of mass and half-mass radius for each snapshot:
! output: i 2(xy yz xz) xyz
     OPEN(3,file=TRIM(newPath)//'/distances_'//proj//'.dat',status='new')
     DO i=1,snapnum
        WRITE(3,30) i,com_cluster(i,1),com_cluster(i,2),com_cluster(i,3), &
             & r_halfmass(i)
     END DO
30   FORMAT(1X,I4,3(2X,F7.4),2X,F7.3)
     CLOSE(3)

!
!*************************************************!
! Mass segregation
!*************************************************!
! Find the degree of mass segregation (lambda).

! All lambda in all planes for stars within rfac*rhalf:

     CALL lambda_setup(i,nstars(i))

  END DO
! End of 'projection' loop

!
!******************************************************************************!
!
! Find total kinetic, gravitational potential and total energy.
! (don't need if it's been called from 'reduce_cluster'
! as all stars are used for this calculation)

  !WRITE(6,*)""
  !WRITE(6,*)"   Calculating cluster energy..."
! Loop over all snapshots
  !DO i=1, snapnum
  !   CALL find_energy(i,nstars(i))
!!$     PRINT *, kinetic_energy(i),potential_energy(i),total_energy(i)
  !END DO
  !WRITE(6,*)"   ...done"
  !WRITE(6,*)""


!*******************************************
! Write out energy data
!
  !OPEN(4,file=TRIM(newPath)//'/energies.dat',status='new')
  !DO i=1,snapnum
  !   WRITE(4,40) i,kinetic_energy(i),potential_energy(i),total_energy(i)
  !END DO
!40 FORMAT(1X,I4,3(2X,E9.3))
  !CLOSE(4)


END SUBROUTINE reduce_rhalf
