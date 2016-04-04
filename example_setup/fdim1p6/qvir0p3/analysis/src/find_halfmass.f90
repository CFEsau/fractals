
SUBROUTINE find_halfmass(snapshoti,ni)
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri, ri_com
  DOUBLE PRECISION :: massi, massj
! rmag = arrays with the magnitude of the separation of all
!        stars from the distribution's com in different planes
  double precision, DIMENSION(:), ALLOCATABLE :: rmag, rmagxy, rmagyz, rmagxz
! rlist = a list of all stars in the distribution
  INTEGER, DIMENSION(:), ALLOCATABLE :: rlist
  INTEGER :: i,j,k
! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(ri_com(1:ni,1:3))
  ALLOCATE(rmagxy(1:ni))
  ALLOCATE(rmagyz(1:ni))
  ALLOCATE(rmagxz(1:ni))
  ALLOCATE(rmag(1:ni))
  ALLOCATE(rlist(1:ni))
  
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  ri_com(1:ni,1:3)=r_com(snapshoti,1:ni,1:3)
! use a separate array for magnitude as we'll be doing a heapsort
! 2D distance magnitudes xy, yz, xz in columns 4:6
  rmagxy(1:ni)=r_com(snapshoti,1:ni,4)
  rmagyz(1:ni)=r_com(snapshoti,1:ni,5)
  rmagxz(1:ni)=r_com(snapshoti,1:ni,6)
! column 7 of r_com is 3D distance magnitude
  rmag(1:ni)=r_com(snapshoti,1:ni,7)

!***********************************
!* 3D half-mass-radius calculation *
!***********************************
  rlist=0
  
! assign IDs to c of m list in prep for heapsort
  do i=1,ni
     rlist(i)=i
  end do
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
  r_halfmass(snapshoti,4)=rmag(rlist(j))+((rmag(rlist(i))-rmag(rlist(j)))*(((0.5*totalmass)-massj)/(massi-massj)))

!************************************
!* 2D half-mass-radius calculations *
!************************************

!xy:
  rlist=0
  do i=1,ni
     rlist(i)=i
  end do
  CALL heapsort(ni,rmagxy,rlist)

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

  r_halfmass(snapshoti,1)=rmagxy(rlist(j))+((rmagxy(rlist(i))-rmagxy(rlist(j)))*(((0.5*totalmass)-massj)/(massi-massj)))


!yz:
  rlist=0
  do i=1,ni
     rlist(i)=i
  end do
  CALL heapsort(ni,rmagyz,rlist)

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

  r_halfmass(snapshoti,2)=rmagyz(rlist(j))+((rmagyz(rlist(i))-rmagyz(rlist(j)))*(((0.5*totalmass)-massj)/(massi-massj)))


!xz:
  rlist=0
  do i=1,ni
     rlist(i)=i
  end do
  CALL heapsort(ni,rmagxz,rlist)

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

  r_halfmass(snapshoti,3)=rmagxz(rlist(j))+((rmagxz(rlist(i))-rmagxz(rlist(j)))*(((0.5*totalmass)-massj)/(massi-massj)))
  
!
  ! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(ri_com)
  DEALLOCATE(rmagxy)
  DEALLOCATE(rmagyz)
  DEALLOCATE(rmagxz)
  DEALLOCATE(rmag)
  DEALLOCATE(rlist)
END SUBROUTINE find_halfmass
