subroutine reduce_rhalf(snapshoti,ni)

  USE sl_input_module
  use constants_module
  use parameters_module
  implicit none
  INTEGER, INTENT(IN) :: snapshoti,ni
  logical :: dirExists
  character(len=100) :: newDir,rfac_char
  integer :: i,j

  limittype='rhalf'
! write value of rfac to string:
  write(rfac_char,'(I0)') rfac

! create 'rhalf' directory for results:
  newDir = 'cluster_rhalf'//trim(rfac_char)//'r'
  inquire(file=trim(outarg)//'/'//trim(newDir)//'/.', exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)

  if (.not. dirExists) then
     write(6,'(a)') "Creating new directory: '"//trim(outarg)//'/'//trim(newDir)//"'"
     call system('mkdir -p '//trim(outarg)//'/'//trim(newDir))
  end if

  write(6,'(a)') "Saving data in '"//trim(outarg)//'/'//trim(newDir)//"'..."


!*************************************************!
! Stars in cluster (within fac*rhalf
!*************************************************!
!
  allocate(incluster(1:snapnum,1:nstars(1),1:4))
  incluster = .true.


  open(10,file=trim(outarg)//'/'//TRIM(newDir)//'/escapedxy.txt')
  open(11,file=trim(outarg)//'/'//TRIM(newDir)//'/escapedyz.txt')
  open(12,file=trim(outarg)//'/'//TRIM(newDir)//'/escapedxz.txt')
  open(13,file=trim(outarg)//'/'//TRIM(newDir)//'/escaped3D.txt')

  do i=1,snapnum
    call in_cluster(i,nstars(i))
  end do

  do j = 10, 13
    close(j)
 end do

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

  ALLOCATE(lambda_tilde(1:snapnum))
  ALLOCATE(l_low_tilde(1:snapnum))
  ALLOCATE(l_up_tilde(1:snapnum))

  ALLOCATE(lambda_star(1:snapnum))
  ALLOCATE(l_low_star(1:snapnum))
  ALLOCATE(l_up_star(1:snapnum))

  ALLOCATE(gamm(1:snapnum))
  ALLOCATE(g_low(1:snapnum))
  ALLOCATE(g_up(1:snapnum))

! IDs of the most massive stars in the cluster:
  allocate(obj_mass(1:2,1:nmst))

!======!
!  xy  !
!======!

  write(6,*)" ------"
  write(6,*)"   xy  "
  write(6,*)" ------"

  projection='xy'

  lambda=0.
  l_up=0.
  l_low=0.

  unit1 = 50
  unit2 = 51
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lam_xy.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lam_xy.txt')

  write(6,*)"       Calculating lambda..."
  do i=1,snapnum
     call find_lambda(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lambar_xy.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lambar_xy.txt')

  write(6,*)"       Calculating lambdabar..."
  do i=1,snapnum
     call find_lambda_bar(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_tilde=0.
  l_up_tilde=0.
  l_low_tilde=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamtilde_xy.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamtilde_xy.txt')

  write(6,*)"       Calculating lambdatilde..."
  do i=1,snapnum
     call find_lambda_tilde(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_star=0.
  l_up_star=0.
  l_low_star=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamstar_xy.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamstar_xy.txt')

  write(6,*)"       Calculating lambdastar..."
  do i=1,snapnum
     call find_lambda_star(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  gamm=0.
  g_up=0.
  g_low=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_gam_xy.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_gam_xy.txt')

  write(6,*)"       Calculating gamma..."
  do i=1,snapnum
     call find_gamma(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

! Write out lambda data:
  OPEN(5,file=TRIM(outarg)//'/'//trim(newDir)//'/lambda_xy',status='new')
  DO i=1,snapnum
     WRITE(5,150) i,lambda(i),l_low(i),l_up(i),lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_tilde(i),l_low_tilde(i),l_up_tilde(i),lambda_star(i), &
          & l_low_star(i),l_up_star(i),gamm(i),g_low(i),g_up(i)
  END DO
  CLOSE(5)


!======!
!  yz  !
!======!

  write(6,*)" ------"
  write(6,*)"   yz  "
  write(6,*)" ------"

  projection='yz'

  lambda=0.
  l_up=0.
  l_low=0.

  unit1 = 60
  unit2 = 61
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lam_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lam_yz.txt')

  write(6,*)"       Calculating lambda..."
  do i=1,snapnum
     call find_lambda(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lambar_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lambar_yz.txt')

  write(6,*)"       Calculating lambdabar..."
  do i=1,snapnum
     call find_lambda_bar(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_tilde=0.
  l_up_tilde=0.
  l_low_tilde=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamtilde_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamtilde_yz.txt')

  write(6,*)"       Calculating lambdatilde..."
  do i=1,snapnum
     call find_lambda_tilde(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_star=0.
  l_up_star=0.
  l_low_star=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamstar_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamstar_yz.txt')

  write(6,*)"       Calculating lambdastar..."
  do i=1,snapnum
     call find_lambda_star(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  gamm=0.
  g_up=0.
  g_low=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_gam_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_gam_yz.txt')

  write(6,*)"       Calculating gamma..."
  do i=1,snapnum
     call find_gamma(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

! Write out lambda data:
  OPEN(5,file=TRIM(outarg)//'/'//trim(newDir)//'/lambda_yz',status='new')
  DO i=1,snapnum
     WRITE(5,150) i,lambda(i),l_low(i),l_up(i),lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_tilde(i),l_low_tilde(i),l_up_tilde(i),lambda_star(i), &
          & l_low_star(i),l_up_star(i),gamm(i),g_low(i),g_up(i)
  END DO
  CLOSE(5)

!======!
!  xz  !
!======!

  write(6,*)" ------"
  write(6,*)"   xz  "
  write(6,*)" ------"

  projection='xz'

  lambda=0.
  l_up=0.
  l_low=0.

  unit1 = 70
  unit2 = 71
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lam_xz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lam_xz.txt')

  write(6,*)"       Calculating lambda..."
  do i=1,snapnum
     call find_lambda(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lambar_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lambar_yz.txt')

  write(6,*)"       Calculating lambdabar..."
  do i=1,snapnum
     call find_lambda_bar(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_tilde=0.
  l_up_tilde=0.
  l_low_tilde=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamtilde_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamtilde_yz.txt')

  write(6,*)"       Calculating lambdatilde..."
  do i=1,snapnum
     call find_lambda_tilde(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_star=0.
  l_up_star=0.
  l_low_star=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamstar_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamstar_yz.txt')

  write(6,*)"       Calculating lambdabar..."
  do i=1,snapnum
     call find_lambda_star(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  gamm=0.
  g_up=0.
  g_low=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_gam_yz.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_gam_yz.txt')

  write(6,*)"       Calculating gamma..."
  do i=1,snapnum
     call find_gamma(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

! Write out lambda data:
  OPEN(5,file=TRIM(outarg)//'/'//trim(newDir)//'/lambda_xz',status='new')
  DO i=1,snapnum
     WRITE(5,150) i,lambda(i),l_low(i),l_up(i),lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_tilde(i),l_low_tilde(i),l_up_tilde(i),lambda_star(i), &
          & l_low_star(i),l_up_star(i),gamm(i),g_low(i),g_up(i)
  END DO
  CLOSE(5)

!======!
!  3D  !
!======!

  write(6,*)" ------"
  write(6,*)"   3D  "
  write(6,*)" ------"

  projection='3D'

  lambda=0.
  l_up=0.
  l_low=0.

  unit1 = 80
  unit2 = 81
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lam_3D.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lam_3D.txt')

  write(6,*)"       Calculating lambda..."
  do i=1,snapnum
     call find_lambda(i,nstars(i))
  end do

  close(unit1)
  close(unit2)     

  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lambar_3D.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lambar_3D.txt')

  write(6,*)"       Calculating lambdabar..."
  do i=1,snapnum
     call find_lambda_bar(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_tilde=0.
  l_up_tilde=0.
  l_low_tilde=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamtilde_3D.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamtilde_3D.txt')

  write(6,*)"       Calculating lambdatilde..."
  do i=1,snapnum
     call find_lambda_tilde(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  lambda_star=0.
  l_up_star=0.
  l_low_star=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_lamstar_3D.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_lamstar_3D.txt')

  write(6,*)"       Calculating lambdastar..."
  do i=1,snapnum
     call find_lambda_star(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  gamm=0.
  g_up=0.
  g_low=0.

  unit1 = unit1 + 2
  unit2 = unit2 + 2
  open(unit1,file=trim(outarg)//'/'//trim(newDir)//'/objm_gam_3D.txt')
  open(unit2,file=trim(outarg)//'/'//trim(newDir)//'/objescaped_gam_3D.txt')

  write(6,*)"       Calculating lambda..."
  do i=1,snapnum
     call find_gamma(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

! Write out lambda data:
  OPEN(5,file=TRIM(outarg)//'/'//trim(newDir)//'/lambda_3D',status='new')
  DO i=1,snapnum
     WRITE(5,150) i,lambda(i),l_low(i),l_up(i),lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_tilde(i),l_low_tilde(i),l_up_tilde(i),lambda_star(i), &
          & l_low_star(i),l_up_star(i),gamm(i),g_low(i),g_up(i)
  END DO
  CLOSE(5)

! Write-out formats:

150 FORMAT(1X,I4,15(2X,F8.3))


!===========================================
  deallocate(incluster)
  deallocate(obj_mass)

  deALLOCATE(lambda)
  deALLOCATE(l_low)
  deALLOCATE(l_up)

  deALLOCATE(lambda_bar)
  deALLOCATE(l_low_bar)
  deALLOCATE(l_up_bar)

  deALLOCATE(lambda_tilde)
  deALLOCATE(l_low_tilde)
  deALLOCATE(l_up_tilde)

  deALLOCATE(lambda_star)
  deALLOCATE(l_low_star)
  deALLOCATE(l_up_star)

  deALLOCATE(gamm)
  deALLOCATE(g_low)
  deALLOCATE(g_up)

end subroutine reduce_rhalf
