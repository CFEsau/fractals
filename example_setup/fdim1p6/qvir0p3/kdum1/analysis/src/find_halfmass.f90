
SUBROUTINE find_halfmass(snapshoti,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
! com_dist = the centre of mass of the distribution
! totalmass = the total mass of the distribution
! rmag = array with the magnitude of the separation of all stars from the distribution com
! rlist = a list of all stars in the distribution
  DOUBLE PRECISION, DIMENSION(3) :: com_dist
  DOUBLE PRECISION :: totalmass, massi, massj
  REAL, DIMENSION(:), ALLOCATABLE :: rmag
  INTEGER, DIMENSION(:), ALLOCATABLE :: rlist
  INTEGER :: i,j,k
! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(rmag(1:ni))
  ALLOCATE(rlist(1:ni))
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
! First find the com of the cluster.
! initialise variables
  totalmass=0.
  com_dist=0.
  rlist=0
  rmag=0.
! Loop over all stars
  DO i=1,ni
     totalmass=totalmass+mi(i)
     com_dist(1:3)=com_dist(1:3)+(mi(i)*ri(i,1:3))
  END DO
  com_dist(1:3)=com_dist(1:3)/totalmass
!
! Next find the distance between each star and the com of the cluster
  DO i=1,ni
     rlist(i)=i
     rmag(i)=REAL(SQRT((ri(i,1)-com_dist(1))**2+(ri(i,2)-com_dist(2))**2+(ri(i,3)-com_dist(3))**2))
  END DO
!
! Sort the stars by distance from the com
  CALL heapsort(ni,rmag,rlist)
! The half mass radius is then the radius at which half of the stelalr mass in inside and half is outside.
! So basically keep adding up the masses of the stars, starting from the centre moving out, until
! half of the mass is inside the radius. There'll be a radius ri-1 where M < totalM/2 and ri > totalM/2.
! Use linear interpolation of these two values to find the halfmass radius.
  i=0
  j=0
  massi=0.
  massj=0.
  DO WHILE(massi<(0.5*totalmass))
     massj=massi
     j=i
     i=i+1
     massi=massi+mi(rlist(i))
  END DO
  ! Linear interpolation to find the halfmass radius
  r_halfmass(snapshoti)=rmag(rlist(j))+((rmag(rlist(i))-rmag(rlist(j)))*(((0.5*totalmass)-massj)/(massi-massj)))
!
  ! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(rmag)
  DEALLOCATE(rlist)
END SUBROUTINE find_halfmass
