
SUBROUTINE c_of_m(snapi,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapi,ni
! mi,ri = mass & position of stars
! (just m & r in 2D array, without snapi)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
! com x, y, z positions of cluster after ejections
  DOUBLE PRECISION :: com_x,com_y,com_z
  INTEGER :: i
  LOGICAL, DIMENSION(:), ALLOCATABLE :: i_incluster

! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:3,1:ni))
  ALLOCATE(i_incluster(1:ni))

! Initialise variables
  mi(1:ni)=mstar(1:ni,snapi)
  ri(1:3,1:ni)=rstar(1:3,1:ni,snapi)
  i_incluster(1:ni)=incluster(1:ni,snapi)
  totalmass=0.
  com_x=0.
  com_y=0.
  com_z=0.

! First find the centre of mass of the cluster.
! Loop over all stars

  DO i=1,ni
     IF (i_incluster(i)) THEN
        totalmass = totalmass + mi(i)
        com_x = com_x + (mi(i)*ri(1,i))
        com_y = com_y + (mi(i)*ri(2,i))
        com_z = com_z + (mi(i)*ri(3,i))
     END IF
  END DO

! divide all by total mass of cluster:
  com_x = com_x / totalmass
  com_y = com_y / totalmass
  com_z = com_z / totalmass

  com_cluster(1,snapi) = com_x
  com_cluster(2,snapi) = com_y
  com_cluster(3,snapi) = com_z

!Then find the distance between each star and the com of the cluster
  DO i=1,ni
!x, y, z distances:
     ri_com(1,i,snapi) = ri(1,i) - com_x
     ri_com(2,i,snapi) = ri(2,i) - com_y
     ri_com(3,i,snapi) = ri(3,i) - com_z
  END DO

! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(i_incluster)

END SUBROUTINE c_of_m
