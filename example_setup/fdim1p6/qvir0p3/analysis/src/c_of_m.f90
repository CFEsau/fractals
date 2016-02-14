
SUBROUTINE c_of_m(snapshoti,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
  DOUBLE PRECISION, DIMENSION(3) :: com_dist
  INTEGER :: i

! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  
! First find the com of the cluster.
! initialise variables
  totalmass=0.
  com_dist=0.
! Loop over all stars
  DO i=1,ni
     totalmass=totalmass+mi(i)
     com_dist(1:3)=com_dist(1:3)+(mi(i)*ri(i,1:3))
  END DO
  
  com_dist(1:3)=com_dist(1:3)/totalmass

!Then find the distance between each star and the com of the cluster
  DO i=1,ni
!x, y, z distances:
     r_com(snapshoti,i,1)=ri(i,1)-com_dist(1)
     r_com(snapshoti,i,2)=ri(i,2)-com_dist(2)
     r_com(snapshoti,i,3)=ri(i,3)-com_dist(3)
!Distance magnitude
!2D (observer), ignore z:
     r_com(snapshoti,i,4)=SQRT(r_com(snapshoti,i,1)**2 + &
          &         r_com(snapshoti,i,2)**2)
!3D:
     r_com(snapshoti,i,5)=SQRT(r_com(snapshoti,i,1)**2 + &
          &         r_com(snapshoti,i,2)**2 + r_com(snapshoti,i,3)**2)
     
  END DO

! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
END SUBROUTINE c_of_m
