SUBROUTINE find_halfmass(snapshoti,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri = mass & position of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
! ri_x, ri_y, ri_z = x, y, z distances of star i from c of m
  DOUBLE PRECISION :: ri_x, ri_y, ri_z
! rmag = distance magnitude of all stars from the distribution's com
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rmag
! rlist = a list of all stars in the distribution
  INTEGER, DIMENSION(:), ALLOCATABLE :: rlist
  DOUBLE PRECISION :: massi, massj
  INTEGER :: i,j,k

! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:3,1:ni))
  ALLOCATE(rmag(1:ni))
  ALLOCATE(rlist(1:ni))

! Assign values
  mi(1:ni)=mstar(1:ni,snapshoti)
  ri(1:3,1:ni)=rstar(1:3,1:ni,snapshoti)

! Find distance magnitude of star i from c of m:
  DO i = 1, ni
     ri_x = ri_com(1,i,snapshoti)
     ri_y = ri_com(2,i,snapshoti)
     ri_z = ri_com(3,i,snapshoti)
     IF (thisproj=='xy') THEN
        rmag(i) = SQRT(ri_x**2 + ri_y**2)
     ELSE IF (thisproj=='yz') THEN
        rmag(i) = SQRT(ri_y**2 + ri_z**2)
     ELSE IF (thisproj=='xz') THEN
        rmag(i) = SQRT(ri_x**2 + ri_z**2)
     ELSE IF (thisproj=='3D') THEN
        rmag(i) = SQRT(ri_x**2 + ri_y**2 + ri_z**2)
     END IF
  END DO

!***********************************
!* 3D half-mass-radius calculation *
!***********************************
  rlist=0
  
! assign IDs to c of m list in prep for heapsort
  DO i=1,ni
     rlist(i)=i
  END DO
!
! Sort the stars by distance from the com
  CALL heapsort(ni,rmag,rlist)
  
! The half mass radius is then the radius at which half of the stellar mass in inside and half is outside.
! Keep adding up the masses of the stars, starting from the centre moving out,
! until half of the mass is inside the radius.
! There'll be a radius ri-1 where M < totalM/2 and ri > totalM/2.
! Use linear interpolation of these two values to find the halfmass radius.

  i=0
  j=0
  massi=0.   ! cumulative stellar mass
  massj=0.

  DO WHILE(massi<(0.5*totalmass))
     massj=massi               ! massj is now the mass within the shell
     j=i
     i=i+1
     massi=massi+mi(rlist(i))  ! cumulative mass in shell plus next star
  END DO

! rhalf is somewhere between r(massj) and r(massi)
  r_halfmass(snapshoti) = rmag(rlist(j)) + &
       & ( (rmag(rlist(i)) - rmag(rlist(j))) * &
       & (((0.5*totalmass)-massj)/(massi-massj)) )

!
  ! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(rmag)
  DEALLOCATE(rlist)

END SUBROUTINE find_halfmass
