SUBROUTINE calc_lambda(obj_avl,ran_avl,lam_val,lam_all,minus,plus,avranmstlen)
! A subroutine to calculate final lambda values; various different 
! measures of lambda (e.g. lambda bar, lambda tilde) can be passed 
! in to save the need for repeating this chunk of code every time 
! we want a new lambda.

  USE parameters_module

  implicit none
  DOUBLE PRECISION, intent(in) :: obj_avl ! Average edge length (obj)
  DOUBLE PRECISION, dimension(1:nloop), intent(in) :: ran_avl ! " (ran)
  double precision, intent(out) :: lam_val,minus,plus
  double precision, intent(out) :: avranmstlen ! Average total length of
                                               ! random msts
  double precision, dimension(1:nloop) :: lam_all ! Lambda for each random MST
  integer :: i
  INTEGER, DIMENSION(1:nloop) :: listID        ! IDs of ran_avl list
  REAL :: ranup, ranlow          ! Upper & lower boundaries, 1 sigma
  
  do i = 1,nloop
     listID(i) = i
  END DO

!Sort nloop random MSTs in order of average edge length:
  CALL heapsort(nloop,ran_avl,listID)

! Find lambda for all random MSTs for lambda CDF plots:
  lam_all(1:nloop)=ran_avl(listID(1:nloop))/obj_avl

!Median of the average MST edge length:
  avranmstlen = ran_avl(listID(NINT(REAL(nloop)/2.)))
!1/6 and 5/6 boundaries for significance:
  ranlow = ran_avl(listID(NINT(REAL(nloop)/6.)))
  ranup = ran_avl(listID(NINT(5.*REAL(nloop)/6.)))

  lam_val = avranmstlen/obj_avl
  minus = ranlow/obj_avl
  plus = ranup/obj_avl

! check lambdas are ordered and match median value using this method:  
!  print*,
!  print*,lam_all(1)
!  print*,lam_all(2)
!  print*,lam_all(3)
!  print*,lam_all(4)
!  print*,lam_all(5)
!  print*,'  :'
!  print*,'  : median mst: ',ran_avl(listID(NINT(REAL(nloop)/2.)))*9.
!  print*,'  : object mst: ',obj_avl*9.
!  print*,'  : median lambda: ',ran_avl(listID(NINT(REAL(nloop)/2.)))/obj_avl
!  print*,'  : median lambda from lam_all: ',lam_all(nint(real(nloop/2.)))
!  print*,'  :'
!  print*,lam_all(nloop-4)
!  print*,lam_all(nloop-3)
!  print*,lam_all(nloop-2)
!  print*,lam_all(nloop-1)
!  print*,lam_all(nloop)

end subroutine calc_lambda
