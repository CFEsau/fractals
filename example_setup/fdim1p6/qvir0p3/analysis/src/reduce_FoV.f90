SUBROUTINE reduce_FoV(ni)

  USE sl_input_module
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ni
  LOGICAL :: dirExists, fileExists
  CHARACTER(len=100) :: newDir, FoV_char
  INTEGER :: i,j

  limittype='FoV'
! write value of FoV_lim (in pc) to string:
  WRITE(FoV_char,'(I0)') FoV_lim

! create 'FoV' directory for results:
  newDir = 'cluster_FoV'//TRIM(FoV_char)//'pc'
  newPath = TRIM(outarg)//'/'//TRIM(newDir)

  INQUIRE(file = TRIM(newPath), exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)

  IF (.NOT. dirExists) THEN
     WRITE(6,'(a)') "Creating new directory: '"//TRIM(newPath)//"'"
     CALL system('mkdir -p '//TRIM(newPath))
  END IF

  WRITE(6,'(a)') "Saving data in '"//TRIM(newPath)//"'..."

!also create directory for lambda data:
  INQUIRE(file = TRIM(newPath)//'/lambda', exist = dirExists)
   IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//TRIM(newPath)//'/lambda')
  END IF


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
! Stars within FoV
!*************************************************!
!
! Initially all stars are in the cluster; none have escaped
!
     incluster = .TRUE.

!
! Set centre of mass of cluster as (0,0,0). True when all in cluster.
     com_cluster=0.

! This means the distance between each star & cluster c of m
! is the distance between ri_com & centre of grid, i.e. r.
! (This assumes all stars are initially in FoV.)
     ri_com=r


! Is the star in the cluster?
     OPEN(10,file=TRIM(newPath)//'/escaped_'//proj//'.dat')

     DO i=1,snapnum
        CALL in_cluster(i,nstars(i))
     END DO

     CLOSE(10)

!
!******************************************************************************!
!
! Find new centre of mass of cluster after ejections.

! Loop over all snapshots
! Calculate com in each case and populate the array
     WRITE(6,*)"       Calculating new centre of mass..."
     DO i=1,snapnum
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
! output: i com_x com_y com_z r1/2
     OPEN(3,file=TRIM(newPath)//'/c_of_m_'//proj//'.dat',status='replace')
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

! All lambda in all planes for stars within 5 pc of centre:

     CALL lambda_setup

  END DO
! End of 'projection' loop

!
!******************************************************************************!
!
! Find total kinetic, gravitational potential and total energy.
! (don't need if it's been called from 'reduce_cluster'
! as all stars are used for this calculation)
! If these files don't exist, do them here:
  INQUIRE(file = TRIM(outarg)//'energies.dat', exist = fileExists)
   IF (.NOT. fileExists) THEN

      WRITE(6,*)""
      WRITE(6,*)"   Calculating cluster energy..."
! Loop over all snapshots
      DO i=1, snapnum
         CALL find_energy(i,nstars(i))
!!$     PRINT *, kinetic_energy(i),potential_energy(i),total_energy(i)
      END DO
      WRITE(6,*)"   ...done"
      WRITE(6,*)""


!*******************************************
! Write out energy data
!
      OPEN(4,file=TRIM(outarg)//'/energies.dat',status='new')
      DO i=1,snapnum
         WRITE(4,40) i,kinetic_energy(i),potential_energy(i),total_energy(i)
      END DO
40    FORMAT(1X,I4,3(2X,E9.3))
      CLOSE(4)
      
   END IF


END SUBROUTINE reduce_FoV
