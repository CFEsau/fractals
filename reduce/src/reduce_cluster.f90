SUBROUTINE reduce_cluster(ni)

  USE sl_input_module
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ni
  LOGICAL :: dirExists
  CHARACTER(len=100) :: newDir
  INTEGER :: i,j

! create 'cluster' directory for results:
  newDir = 'cluster_all'
  newPath = TRIM(outarg)//'/'//TRIM(newDir)

  INQUIRE(file = TRIM(newPath), exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)
  
  IF (.NOT. dirExists) THEN
     WRITE(6,'(a)') "Creating new directory: '"//TRIM(newPath)//"'"
     CALL system('mkdir -p '//TRIM(newPath))
  END IF
  
  WRITE(6,*)""
  WRITE(6,'(a)') "Saving data in '"//TRIM(newPath)//"'..."
  
  DO projnum = 1,4
     IF (projnum==1) THEN
        thisproj='xy'
        WRITE(6,*)""
        WRITE(6,*)"   "//thisproj//":"
     ELSE IF (projnum==2) THEN
        thisproj='yz'
        WRITE(6,*)""
        WRITE(6,*)"   "//thisproj//":"
     ELSE IF (projnum==3)  THEN
        thisproj='xz'
        WRITE(6,*)""
        WRITE(6,*)"   "//thisproj//":"
     ELSE IF (projnum==4)  THEN
        thisproj='3D'
        WRITE(6,*)""
        WRITE(6,*)"   "//thisproj//":"
     END IF

!*************************************************!
! Stars in cluster
!*************************************************!
!
! Initially all stars are in the cluster; none have escaped
!
     incluster = .TRUE.


! Don't need to call this here as all stars included
!     open(10,file=trim(newPath)//'/escaped_'//thisproj//'.dat')
!
!     do i=1,nsnaps
!       call in_cluster(i,nstars(i))
!     end do
!
!     close(10)

!
!******************************************************************************!
!
! Find centre of mass of cluster after ejections.
!
! Set centre of mass of cluster as (0,0,0). True when all in cluster.
     com_cluster=0.

! Find distance between each star and cluster centre of mass.
     ri_com=rstar

! Loop over all snapshots
! Calculate com in each case and populate the array
     WRITE(6,*)"       Calculating centre of mass..."
     DO i=1, nsnaps
        CALL c_of_m(i,nstars(i))
     END DO
     
!******************************************************************************!
!
! Find the Half mass radius.

     r_halfmass=0.

     WRITE(6,*)"       Calculating half-mass radius..."
! Loop over all snapshots
     DO i=1,nsnaps
        CALL find_halfmass(i,nstars(i))
!!$     PRINT *, i, r_halfmass(i)
     END DO


!*******************************************
! Write out distance data
!
! Centre of mass and half-mass radius:
! output: i com_x com_y com_z r1/2
     OPEN(3,file=TRIM(newPath)//'/c_of_m_'//thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(3,30) i,com_cluster(1,i),com_cluster(2,i),com_cluster(3,i), &
             & r_halfmass(i)
     END DO
30   FORMAT(1X,I4,3(2X,F7.4),2X,F7.3)
     CLOSE(3)

!
!*************************************************!
! Mass segregation
!*************************************************!
! Find the degree of mass segregation (lambda).

! All lambda in all planes with all stars in cluster:

     CALL lambda_setup

  END DO
! End of 'projection' loop

!
!******************************************************************************!
!
! Find total kinetic, gravitational potential and total energy.

  WRITE(6,*)""
  WRITE(6,*)"   Calculating cluster energy..."
! Loop over all snapshots
  DO i=1, nsnaps
     CALL find_energy(i,nstars(i))
!!$     PRINT *, kinetic_energy(i),potential_energy(i),total_energy(i)
  END DO
  WRITE(6,*)"   ...done"
  WRITE(6,*)""


!*******************************************
! Write out energy data
!
! Save in 'outarg' as this is the same for all cluster types
! (don't need to call from FoV / rhalf cluster types)
  OPEN(4,file=TRIM(outarg)//'/energies.dat',status='replace')
  DO i=1,nsnaps
     WRITE(4,40) i,kinetic_energy(i),potential_energy(i),total_energy(i)
  END DO
40 FORMAT(1X,I4,3(2X,E9.3))
  CLOSE(4)


END SUBROUTINE reduce_cluster
