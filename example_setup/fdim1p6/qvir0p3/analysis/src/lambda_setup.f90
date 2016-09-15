subroutine lambda_setup(snapshoti,ni)

  USE sl_input_module
  use constants_module
  use parameters_module
  implicit none
  INTEGER, INTENT(IN) :: snapshoti,ni
  integer :: i,j

!
!*************************************************!
! Mass segregation
!*************************************************!
!
! Find the degree of mass segregation (lambda).

  ALLOCATE(lambda(1:snapnum))
  ALLOCATE(l_low(1:snapnum))
  ALLOCATE(l_up(1:snapnum))

  ALLOCATE(lambda_bar(1:snapnum))
  ALLOCATE(l_low_bar(1:snapnum))
  ALLOCATE(l_up_bar(1:snapnum))

  ALLOCATE(lambda_til(1:snapnum))
  ALLOCATE(l_low_til(1:snapnum))
  ALLOCATE(l_up_til(1:snapnum))

  ALLOCATE(lambda_star(1:snapnum))
  ALLOCATE(l_low_star(1:snapnum))
  ALLOCATE(l_up_star(1:snapnum))

  ALLOCATE(gamm(1:snapnum))
  ALLOCATE(g_low(1:snapnum))
  ALLOCATE(g_up(1:snapnum))

! IDs of the most massive stars in the cluster:
  allocate(obj_mass(1:2,1:nmst))

  write(6,*)"       Calculating lambda..."

  lambda=0.
  l_up=0.
  l_low=0.

  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.

  lambda_til=0.
  l_up_til=0.
  l_low_til=0.

  lambda_star=0.
  l_up_star=0.
  l_low_star=0.

  gamm=0.
  g_up=0.
  g_low=0.

  open(50,file=trim(newPath)//'/objm_'//proj//'.dat')
  !open(51,file=trim(newPath)//'/objescaped_'//proj//'.dat')
!(only need this if you want to check distances of escaped object stars)


!Record any stars that fall on top of each other &
!need their separation changing for the mst
  open(4,file=TRIM(newPath)//'/sep_obj.dat',position='append')
  open(5,file=TRIM(newPath)//'/sep.dat',position='append')
  write(4,*) "**** ",proj," ****"
  write(5,*) "**** ",proj," ****"

  do i=1,snapnum
     call find_lambda_all(i,nstars(i))
  end do

  close(4)
  close(5)

  write(6,*)"       ...done"

  close(50)
  !close(51)

! Write out lambda data:
  OPEN(5,file=trim(newPath)//'/lambda_'//proj,status='new')
  DO i=1,snapnum
     WRITE(5,150) i,lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_til(i),l_low_til(i),l_up_til(i),lambda_star(i), &
          & l_low_star(i),l_up_star(i),gamm(i),g_low(i),g_up(i)
  END DO
  CLOSE(5)

!TODO: change writeout to vary depending on which lambda values are calculated
!e.g. if calclambda; values=values + lambda(i) + lambdalow(i) + lambdaup(i);
! form = form + 2X,F8.3; end if

150 FORMAT(1X,I4,12(2X,F9.3))

!===========================================
  deallocate(obj_mass)

  deALLOCATE(lambda)
  deALLOCATE(l_low)
  deALLOCATE(l_up)

  deALLOCATE(lambda_bar)
  deALLOCATE(l_low_bar)
  deALLOCATE(l_up_bar)

  deALLOCATE(lambda_til)
  deALLOCATE(l_low_til)
  deALLOCATE(l_up_til)

  deALLOCATE(lambda_star)
  deALLOCATE(l_low_star)
  deALLOCATE(l_up_star)

  deALLOCATE(gamm)
  deALLOCATE(g_low)
  deALLOCATE(g_up)

end subroutine lambda_setup
