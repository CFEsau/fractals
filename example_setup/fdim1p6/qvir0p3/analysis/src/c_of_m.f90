
SUBROUTINE c_of_m(snapshoti,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
! com x, y, z positions
  DOUBLE PRECISION :: com_x,com_y,com_z
  INTEGER :: i
  LOGICAL, DIMENSION(:), ALLOCATABLE :: i_incluster

! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(i_incluster(1:ni))

! Initialise variables
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  i_incluster(1:ni)=incluster(snapshoti,1:ni)
  totalmass=0.
  com_x=0.
  com_y=0.
  com_z=0.

! First find the centre of mass of the cluster.
! Loop over all stars

  DO i=1,ni
     IF (i_incluster(i)) THEN
        totalmass = totalmass + mi(i)
        com_x = com_x + (mi(i)*ri(i,1))
        com_y = com_y + (mi(i)*ri(i,2))
        com_z = com_z + (mi(i)*ri(i,3))
     END IF
  END DO

! divide all by total mass of cluster:
  com_x = com_x / totalmass
  com_y = com_y / totalmass
  com_z = com_z / totalmass

  com_cluster(snapshoti,1) = com_x
  com_cluster(snapshoti,2) = com_y
  com_cluster(snapshoti,3) = com_z

!Then find the distance between each star and the com of the cluster
  DO i=1,ni
!x, y, z distances:
     ri_com(snapshoti,i,1) = ri(i,1) - com_x
     ri_com(snapshoti,i,2) = ri(i,2) - com_y
     ri_com(snapshoti,i,3) = ri(i,3) - com_z
  END DO

! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(i_incluster)

END SUBROUTINE c_of_m
