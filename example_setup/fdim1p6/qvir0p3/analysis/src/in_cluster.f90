!If certain stars are to be ignored then make a new array
!of stars in each snapshot where nstars = n_in_cluster
SUBROUTINE in_cluster(snapshoti,ni)
  
  USE sl_input_module
  USE parameters_module
  USE constants_module

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: snapshoti,ni
! n_escaped = number of escaped stars in given snapshot & projection
  INTEGER :: n_escaped
! ri_x, ri_y, ri_z = distance in x, y, z of star i from c of m
  double precision :: ri_x, ri_y, ri_z
! rmag = distance magnitude between star i and centre of mass
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rmag
  INTEGER :: i,j
! number of stars in given projection of cluster:
  INTEGER :: n_proj
! masses are just for data output; nothing to do with calculations
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi

! Allocate memory for arrays
  ALLOCATE(rmag(1:ni))
  ALLOCATE(mi(1:ni))

  mi(1:ni)=m(snapshoti,1:ni)
! Populate distance magnitude arrays in observer planes & 3D
! between each star and centre of mass
  do i = 1, ni
     ri_x = ri_com(snapshoti,i,1)
     ri_y = ri_com(snapshoti,i,2)
     ri_z = ri_com(snapshoti,i,3)
     if (proj=='xy') then
        rmag(i) = sqrt(ri_x**2 + ri_y**2)
     else if (proj=='yz') then
        rmag(i) = sqrt(ri_y**2 + ri_z**2)
     else if (proj=='xz') then
        rmag(i) = sqrt(ri_x**2 + ri_z**2)
     else if (proj=='3D') then
        rmag(i) = sqrt(ri_x**2 + ri_y**2 + ri_z**2)
     end if
  end do

! initially all stars are in the cluster
  n_proj = ni

! Loop over all stars
  DO i=1,ni

!===============
! Field of view
!===============

! if rmag is greater than FoV_lim, it has left the cluster
     IF (limittype=='FoV' .and. rmag(i) > FoV_lim) THEN
        incluster(snapshoti,i) = .FALSE.
        n_proj = n_proj - 1
     END IF


!=================
! halfmass radius
!=================

! if rmag is greater than rfac*r_halfmass, it has left the cluster
     IF (limittype=='rhalf' .and. &
          & rmag(i) > rfac*rhalf_all(projnum)) THEN
        incluster(snapshoti,i) = .FALSE.
        n_proj = n_proj - 1
     END IF
  END DO


! ***********
! * Outputs *
! ***********
! Write out number of escaped stars if value
! is different from previous snapshot.

! Set an unphysical value for n_escaped in 1st snapshot
! so we are guaranteed to get a writeout:
  IF (snapshoti==1) THEN
     n_escaped=-1
  END IF


! if new stars escaped, write snapshot number
  IF (ni - n_proj /= n_escaped) THEN
     WRITE(10,'(2X,"Snapshot",I4,":",1X,I4,1X,"stars escaped")') &
          & snapshoti, ni-n_proj

     DO i=1,ni
! if star isn't in cluster, write distance magnitude & star mass
        IF (.NOT. incluster(snapshoti,i)) WRITE(10,80) i, rmag(i), mi(i)
     END DO

     IF (limittype=='rhalf') THEN
        WRITE(10,81) rfac, rfac * rhalf_all(projnum)
     END IF

! write number of stars escaped:
     n_escaped = ni - n_proj
! For some reason need this in or the above 'if' is ignored...
     WRITE(10,'(30X,I4,1X,I4)') ni-n_proj, n_escaped
     WRITE(10,*)""
  END IF

80   FORMAT(1X,I4,2X,F8.3,2X,F7.3)
81   FORMAT(4X,I2,'*r_halfmass:',F8.3)


  DEALLOCATE(rmag)
  DEALLOCATE(mi)

END SUBROUTINE in_cluster
