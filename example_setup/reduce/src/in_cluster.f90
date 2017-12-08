!If certain stars are to be ignored then make a new array
!of stars in each snapshot where nstars = n_in_cluster
SUBROUTINE in_cluster(snapi,ni)
  
  USE sl_input_module
  USE parameters_module
  USE constants_module

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: snapi,ni
! n_escaped = number of escaped stars in given snapshot & projection
  INTEGER :: n_escaped
! ri_x, ri_y, ri_z = distance in x, y, z of star i from cluster centre
  double precision :: ri_x, ri_y, ri_z
! rmag = distance magnitude between each star and centre of grid
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rmag
  INTEGER :: i,j
! number of stars in given projection of cluster:
  INTEGER :: n_proj
! masses are just for data output; nothing to do with calculations
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  logical :: writeescaped

! Allocate memory for arrays
  ALLOCATE(rmag(1:ni))
  ALLOCATE(mi(1:ni))

! Write out list of escaped stars?
! This file in full is large so use only if needed.
  writeescaped=.FALSE.

  mi(1:ni)=mstar(1:ni,snapi)
! Populate distance magnitude arrays in observer planes
! & 3D between each star and centre of the cluster
  do i = 1, ni
     ri_x = rstar(1,i,snapi)
     ri_y = rstar(2,i,snapi)
     ri_z = rstar(3,i,snapi)
     if (thisproj=='xy') then
        rmag(i) = sqrt(ri_x**2 + ri_y**2)
     else if (thisproj=='yz') then
        rmag(i) = sqrt(ri_y**2 + ri_z**2)
     else if (thisproj=='xz') then
        rmag(i) = sqrt(ri_x**2 + ri_z**2)
     else if (thisproj=='3D') then
        rmag(i) = sqrt(ri_x**2 + ri_y**2 + ri_z**2)
     end if
  end do

! The distribution is recentred on (0,0,0) for each snapshot
! and positions are written relative to this, so just need to
! find which stars fall in a circle of radius 5 pc.

! initially all stars are in the cluster
  n_proj = ni

! Loop over all stars
  DO i=1,ni

!===============
! Field of view
!===============

! if rmag is greater than FoV_lim, it has left the cluster
     IF (limittype=='FoV' .and. rmag(i) > FoV_lim) THEN
        incluster(i,snapi) = .FALSE.
        n_proj = n_proj - 1
     END IF
! This only sets star to F in current snapshot. Initially all T so
! if it re-enters, it should be kept as T in following snapshot.


!=================
! halfmass radius
!=================

! if rmag is greater than rfac*r_halfmass, it has left the cluster
     IF (limittype=='rhalf' .and. &
          & rmag(i) > rfac*rhalf_all(projnum)) THEN
        incluster(i,snapi) = .FALSE.
        n_proj = n_proj - 1
     END IF
  END DO
! This only sets star to F in current snapshot. Initially all T so
! if it re-enters, it should be kept as T in following snapshot.


! ***********
! * Outputs *
! ***********
! Write out number of escaped stars if value
! is different from previous snapshot.

! Set an unphysical value for n_escaped in 1st snapshot
! so we are guaranteed to get a writeout:
  IF (snapi==1) THEN
     n_escaped=-1
  END IF
! (This is actually physical if a star re-enters,
! but this won't happen in 1st snapshot)


! if new stars escaped & want this in output file, write snapshot number
  IF (ni-n_proj/=n_escaped) THEN
     n_escaped = ni-n_proj !update number of escaped stars
     
     WRITE(10,'(2X,"Snapshot",I4,":",1X,I4,1X,"stars escaped")') &
          & snapi, n_escaped

     IF (writeescaped) THEN
! if star isn't in cluster, write distance magnitude & star mass
        DO i=1,ni
           IF (.NOT. incluster(i,snapi)) WRITE(10,80) i, rmag(i), mi(i)
        END DO

        IF (limittype=='rhalf') THEN
           WRITE(10,81) rfac, rfac * rhalf_all(projnum)
        END IF
     END IF
  END IF

80   FORMAT(1X,I4,2X,F8.3,2X,F7.3)
81   FORMAT(4X,I2,'*r_halfmass:',F8.3)


  DEALLOCATE(rmag)
  DEALLOCATE(mi)

END SUBROUTINE in_cluster
