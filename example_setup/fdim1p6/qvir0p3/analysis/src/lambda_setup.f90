subroutine lambda_setup

  USE sl_input_module
  use constants_module
  use parameters_module
  implicit none
  integer :: i,j

!
!*************************************************!
! Mass segregation
!*************************************************!
!
! Lengths of the edges of object stars:
  ALLOCATE(edgelengths(1:nmst-1))

! Find the degree of mass segregation (lambda).

  ALLOCATE(lambda(1:snapnum))
  ALLOCATE(l_low(1:snapnum))
  ALLOCATE(l_up(1:snapnum))
  ALLOCATE(lam_avranmst(1:snapnum))
  ALLOCATE(lam_objmst(1:snapnum))

  ALLOCATE(lambda_bar(1:snapnum))
  ALLOCATE(l_low_bar(1:snapnum))
  ALLOCATE(l_up_bar(1:snapnum))
  ALLOCATE(lbar_avranmst(1:snapnum))
  ALLOCATE(lbar_objmst(1:snapnum))

  ALLOCATE(lambda_til(1:snapnum))
  ALLOCATE(l_low_til(1:snapnum))
  ALLOCATE(l_up_til(1:snapnum))
  ALLOCATE(ltil_avranmst(1:snapnum))
  ALLOCATE(ltil_objmst(1:snapnum))

  ALLOCATE(lambda_star(1:snapnum))
  ALLOCATE(l_low_star(1:snapnum))
  ALLOCATE(l_up_star(1:snapnum))
  ALLOCATE(lstar_avranmst(1:snapnum))
  ALLOCATE(lstar_objmst(1:snapnum))

  ALLOCATE(gamm(1:snapnum))
  ALLOCATE(g_low(1:snapnum))
  ALLOCATE(g_up(1:snapnum))
  ALLOCATE(gam_avranmst(1:snapnum))
  ALLOCATE(gam_objmst(1:snapnum))

! IDs of the most massive stars in the cluster:
  allocate(obj_mass(1:2,1:nmst))

  write(6,*)"       Calculating lambda..."

  lambda=0.
  l_up=0.
  l_low=0.
  lam_avranmst=0.
  lam_objmst=0.

  lambda_bar=0.
  l_up_bar=0.
  l_low_bar=0.
  lbar_avranmst=0.
  lbar_objmst=0.

  lambda_til=0.
  l_up_til=0.
  l_low_til=0.
  ltil_avranmst=0.
  ltil_objmst=0.

  lambda_star=0.
  l_up_star=0.
  l_low_star=0.
  lstar_avranmst=0.
  lstar_objmst=0.

  gamm=0.
  g_up=0.
  g_low=0.
  gam_avranmst=0.
  gam_objmst=0.

  open(50,file=trim(newPath)//'/obj_m_'//proj//'.dat')
  !open(51,file=trim(newPath)//'/objescaped_lam_'//proj//'.dat')
!(only need this if you want to check distances of escaped object stars)
  OPEN(52,file=trim(newPath)//'/MSTedges_'//proj//'.dat',status='new')


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
  close(52)

! Write out lambda values with errors and
! numerator & denominator (MST lengths) used for each lambda:
  OPEN(4,file=TRIM(newPath)//'/MSTlengths_'//proj//'.dat',status='new')
  OPEN(5,file=trim(newPath)//'/lambda_'//proj//'.dat',status='new')
  
  DO i=1,snapnum
     WRITE(4,104) i,lam_avranmst(i),lam_objmst(i),lbar_avranmst(i), &
          & lbar_objmst(i),ltil_avranmst(i),ltil_objmst(i), &
          & lstar_avranmst(i),lstar_objmst(i),gam_avranmst(i),gam_objmst(i)
     
     WRITE(5,105) i,lambda_bar(i),l_low_bar(i),l_up_bar(i), &
          & lambda_til(i),l_low_til(i),l_up_til(i),lambda_star(i), &
          & l_low_star(i),l_up_star(i),gamm(i),g_low(i),g_up(i)
  END DO
  
104 FORMAT(1X,I4,10(2X,F9.4))
105 FORMAT(1X,I4,12(2X,F8.3))
  CLOSE(4)
  CLOSE(5)

!TODO: change writeout to vary depending on which lambda values are calculated
!e.g. if calclambda; values=values + lambda(i) + lambdalow(i) + lambdaup(i);
! form = form + 2X,F8.3; end if


!===========================================
  deallocate(obj_mass)
  deallocate(edgelengths)

  deALLOCATE(lambda)
  deALLOCATE(l_low)
  deALLOCATE(l_up)
  deallocate(lam_avranmst)
  deallocate(lam_objmst)

  deALLOCATE(lambda_bar)
  deALLOCATE(l_low_bar)
  deALLOCATE(l_up_bar)
  deallocate(lbar_avranmst)
  deallocate(lbar_objmst)

  deALLOCATE(lambda_til)
  deALLOCATE(l_low_til)
  deALLOCATE(l_up_til)
  deallocate(ltil_avranmst)
  deallocate(ltil_objmst)

  deALLOCATE(lambda_star)
  deALLOCATE(l_low_star)
  deALLOCATE(l_up_star)
  deallocate(lstar_avranmst)
  deallocate(lstar_objmst)

  deALLOCATE(gamm)
  deALLOCATE(g_low)
  deALLOCATE(g_up)
  deallocate(gam_avranmst)
  deallocate(gam_objmst)

end subroutine lambda_setup
