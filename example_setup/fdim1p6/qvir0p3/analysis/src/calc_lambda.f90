SUBROUTINE calc_lambda(obj_mst,rand_mst,lam_val,minus,plus)
! A subroutine to calculate final lambda values; various different 
! measures of lambda (e.g. lambda bar, lambda tilde) can be passed 
! in to save the need for repeating this chunk of code every time 
! we want a new lambda.

  USE parameters_module

  implicit none
  DOUBLE PRECISION, intent(in) :: obj_mst
  DOUBLE PRECISION, dimension(1:nloop), intent(in) :: rand_mst
  double precision, intent(out) :: lam_val,minus,plus
  integer :: i
  INTEGER, DIMENSION(1:nloop) :: listID !IDs of rand_mst list
  REAL :: avranmst !Average length of random trees
  REAL :: ranup, ranlow !Upper & lower boundaries, 1 sigma

  do i = 1,nloop
     listID(i) = i
  END DO

!Sort random MST lengths:
  CALL heapsort(nloop,rand_mst,listID)

  avranmst = rand_mst(listID(NINT(REAL(nloop)/2.))) !Median MST.
  ranlow = rand_mst(listID(NINT(REAL(nloop)/6.)))
  ranup = rand_mst(listID(NINT(5.*REAL(nloop)/6.)))

  lam_val = avranmst/obj_mst
  minus = ranlow/obj_mst
  plus = ranup/obj_mst

end subroutine calc_lambda
