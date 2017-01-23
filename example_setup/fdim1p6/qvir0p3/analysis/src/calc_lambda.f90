SUBROUTINE calc_lambda(obj_avl,ran_avl,lam_val,minus,plus,avranmstlen)
! A subroutine to calculate final lambda values; various different 
! measures of lambda (e.g. lambda bar, lambda tilde) can be passed 
! in to save the need for repeating this chunk of code every time 
! we want a new lambda.

  USE parameters_module

  implicit none
  DOUBLE PRECISION, intent(in) :: obj_avl !average edge length (obj)
  DOUBLE PRECISION, dimension(1:nloop), intent(in) :: ran_avl ! " (ran)
  double precision, intent(out) :: lam_val,minus,plus
  double precision, intent(out) :: avranmstlen !Average tot length of random trees
  integer :: i
  INTEGER, DIMENSION(1:nloop) :: listID !IDs of ran_avl list
  REAL :: ranup, ranlow !Upper & lower boundaries, 1 sigma
  
  do i = 1,nloop
     listID(i) = i
  END DO

!Sort nloop random MSTs in order of average edge length:
  CALL heapsort(nloop,ran_avl,listID)

!Median of the average MST edge length:
  avranmstlen = ran_avl(listID(NINT(REAL(nloop)/2.)))
!1/6 and 5/6 boundaries for significance:
  ranlow = ran_avl(listID(NINT(REAL(nloop)/6.)))
  ranup = ran_avl(listID(NINT(5.*REAL(nloop)/6.)))

  lam_val = avranmstlen/obj_avl
  minus = ranlow/obj_avl
  plus = ranup/obj_avl

end subroutine calc_lambda
