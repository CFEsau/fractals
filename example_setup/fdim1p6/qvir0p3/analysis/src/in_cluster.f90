!If certain stars are to be ignored then make a new array
!of stars in each snapshot where nstars = n_in_cluster
SUBROUTINE in_cluster(snapshoti,ni)
  
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! rmag = distance between star i and centre of mass
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rmag2D, rmag3D
  INTEGER :: i
! ignoring stars beyond a certain distance:
! r_lim is distance limit of allowed stars in pc
! fac is factor multiplied by r_halfmass for alternative distance
  real :: r_lim
  integer :: fac
  integer :: n_2Dcluster, n_3Dcluster !number of stars in the cluster
  
  ALLOCATE(rmag2D(1:ni))
  ALLOCATE(rmag3D(1:ni))
  rmag2D(1:ni)=r_com(snapshoti,1:ni,4)
  rmag3D(1:ni)=r_com(snapshoti,1:ni,5)

! initially all stars are in cluster
  n_2Dcluster = ni
  n_3Dcluster = ni
  !ignore stars when one r co-ord > this limit (it's outside the cluster)
  r_lim = 5.
  fac = 2

  write(10,*) 'Snapshot',snapshoti
! Loop over all stars
  DO i=1,ni


! if rmag is graeter than 5 pc & fac*r_halfmass, it has left the cluster
! (calculates both 2D and 3D distance magnitude)
     if ((rmag2D(i) > r_lim) .and. (rmag2D(i) > real(fac)*r_halfmass(snapshoti))) then
        this_in2Dcluster(snapshoti,i) = .false.
        n_2Dcluster = n_2Dcluster-1
        write(10,*) i, rmag3D(i)
     end if

     if ((rmag3D(i) > r_lim) .and. (rmag3D(i) > real(fac)*r_halfmass(snapshoti))) then
        this_in3Dcluster(snapshoti,i) = .false.
        n_3Dcluster = n_3Dcluster-1
        write(11,*) i, rmag3D(i)
     end if

  end do

!Write out number of ignored stars (if any have been ignored)
  if (ni-n_2Dcluster>0) write(10,*) fac,'* r_halfmass:',real(fac)*r_halfmass(snapshoti)
     write(11,*) 'Ignored stars:',ni-n_2Dcluster

  if (ni-n_3Dcluster>0) write(11,*) fac,'* r_halfmass:',real(fac)*r_halfmass(snapshoti)
     write(11,*) 'Ignored stars:',ni-n_3Dcluster

  DEALLOCATE(rmag2D)
  DEALLOCATE(rmag3D)
END SUBROUTINE in_cluster
