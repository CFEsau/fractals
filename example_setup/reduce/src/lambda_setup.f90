subroutine lambda_setup
! Allocate, zero, and deallocate various arrays needed for MST & lambda.
! File I/O for obj star masses, lambda values, edge lengths for CDFs, etc.

  USE sl_input_module
  use constants_module
  use parameters_module
  implicit none
  integer :: i,j
  character(len=6) :: filei !string version of 'i'
  integer :: nlam !number of lambda types
  integer :: lamunit, lambarunit, lamrmsunit, lamsmrunit !output units
  integer :: lamharunit, lamtilunit, lamstarunit, gamunit, lamlnunit
  integer :: lamNmedunit
  character(len=100) :: CDFpath
  LOGICAL :: dirExists
  !could just use 'lamunit' & increment but the have to keep track
  !of the order in which different lambda measures are calculated
  
!lambda: uses average random mst length & total obj length
!lambda bar: uses mean edge lengths (arithmetic mean)
!            Equivalent to lambda as 1/n cancels
!lambda rms: uses root mean square. p=2 in 'generalised mean' eq.
!lambda smr: uses square mean root (if that's a thing...). p=1/2.
!lambda har: generalised mean with p=-1
!lambda til: uses median edge length
!lambda star: uses mst made of median edge length, then add total MST
!gamma: uses geometric mean
!lambda ln: uses generalised f-mean, where f^(-1) is ln
!lambda N median: takes the 2 or 3 median points (2 if even number in MST)
!                 then finds the mean of these

!======================================================
! Set # of stars in MST, # of random MSTs, & # of CDFs
!======================================================
  
  nmst=10 !number of stars in the MST
  
!For large nmst it can take a while to build the MST, so reduce
!the number of times we do the loop - higher nmst MSTs are less 
!stochastic anyway...
  nloop = 1000
  IF(nmst >= 100) nloop = 50
  
  nCDF = 20 !Number of CDFs plotted for random MSTs; must be =< nloop
  CDFPath = TRIM(newPath)//'/CDFdata'
  !directory for CDF data (edge lengths in MST for different subsets)
!======================================================

! Lengths of the edges of object stars:
  ALLOCATE(edgelengths(1:nmst-1))

! IDs of the most massive stars in the cluster:
  allocate(obj_mass(1:2,1:nmst))

!
!*******************************************************!
! Mass segregation
!
!########################################  
!Set up different types types of average
!######################################## 

! Find the degree of mass segregation, lambda, using
! average MST edge lengths, with different types of average.
! Average edge lengths for object stars and
! random stars are saved for each snapshot.

  nlam=0 !number of different types of lambda
  
  findlam = .FALSE.
  findlambar = .TRUE.
  findlamrms = .TRUE.
  findlamsmr = .TRUE.
  findlamhar = .TRUE.
  findlamtil = .TRUE.
  findlamNmed = .TRUE.
  findlamstar = .FALSE.
  findgam = .TRUE.
  findlamln = .TRUE.
  
  IF (findlam) THEN
     nlam = nlam + 1
     lamunit = nlam+100
     
! final lambda value calculated (output)
     ALLOCATE(lambda(1:snapnum))
! +ve error bar
     ALLOCATE(l_low(1:snapnum))
! -ve error bar
     ALLOCATE(l_up(1:snapnum))
! median average edge length of random MTSs
     ALLOCATE(l_avranmst(1:snapnum))
! average edge length of object MST
     ALLOCATE(l_objmst(1:snapnum))

     lambda=0.
     l_up=0.
     l_low=0.
     l_avranmst=0.
     l_objmst=0.
  END IF

  IF (findlambar) THEN
     nlam = nlam + 1
     lambarunit = nlam+100
     
     ALLOCATE(lambda_bar(1:snapnum))
     ALLOCATE(l_low_bar(1:snapnum))
     ALLOCATE(l_up_bar(1:snapnum))
     ALLOCATE(lbar_avranmst(1:snapnum))
     ALLOCATE(lbar_objmst(1:snapnum))
     
     lambda_bar=0.
     l_up_bar=0.
     l_low_bar=0.
     lbar_avranmst=0.
     lbar_objmst=0.
  END IF

  IF (findlamrms) THEN
     nlam = nlam + 1
     lamrmsunit = nlam + 100
     
     ALLOCATE(lambda_rms(1:snapnum))
     ALLOCATE(l_low_rms(1:snapnum))
     ALLOCATE(l_up_rms(1:snapnum))
     ALLOCATE(lrms_avranmst(1:snapnum))
     ALLOCATE(lrms_objmst(1:snapnum))

     lambda_rms=0.
     l_up_rms=0.
     l_low_rms=0.
     lrms_avranmst=0.
     lrms_objmst=0.
  END IF

  IF (findlamsmr) THEN
     nlam = nlam + 1
     lamsmrunit = nlam + 100
     
     ALLOCATE(lambda_smr(1:snapnum))
     ALLOCATE(l_low_smr(1:snapnum))
     ALLOCATE(l_up_smr(1:snapnum))
     ALLOCATE(lsmr_avranmst(1:snapnum))
     ALLOCATE(lsmr_objmst(1:snapnum))

     lambda_smr=0.
     l_up_smr=0.
     l_low_smr=0.
     lsmr_avranmst=0.
     lsmr_objmst=0.
  END IF

  IF (findlamhar) THEN
     nlam = nlam + 1
     lamharunit = nlam + 100
     
     ALLOCATE(lambda_har(1:snapnum))
     ALLOCATE(l_low_har(1:snapnum))
     ALLOCATE(l_up_har(1:snapnum))
     ALLOCATE(lhar_avranmst(1:snapnum))
     ALLOCATE(lhar_objmst(1:snapnum))

     lambda_har=0.
     l_up_har=0.
     l_low_har=0.
     lhar_avranmst=0.
     lhar_objmst=0.
  END IF
  
  IF (findlamtil) THEN
     nlam = nlam + 1
     lamtilunit = nlam + 100
     
     ALLOCATE(lambda_til(1:snapnum))
     ALLOCATE(l_low_til(1:snapnum))
     ALLOCATE(l_up_til(1:snapnum))
     ALLOCATE(ltil_avranmst(1:snapnum))
     ALLOCATE(ltil_objmst(1:snapnum))
     
     lambda_til=0.
     l_up_til=0.
     l_low_til=0.
     ltil_avranmst=0.
     ltil_objmst=0.
  END IF
  
  IF (findlamNmed) THEN
     nlam = nlam + 1
     lamlnunit = nlam + 100
     
     ALLOCATE(lambda_Nmed(1:snapnum))
     ALLOCATE(l_low_Nmed(1:snapnum))
     ALLOCATE(l_up_Nmed(1:snapnum))
     ALLOCATE(lNmed_avranmst(1:snapnum))
     ALLOCATE(lNmed_objmst(1:snapnum))
     
     lambda_Nmed=0.
     l_up_Nmed=0.
     l_low_Nmed=0.
     lNmed_avranmst=0.
     lNmed_objmst=0.
  END IF
  
  IF (findlamstar) THEN
     nlam = nlam + 1
     lamstarunit = nlam + 100
     
     ALLOCATE(lambda_star(1:snapnum))
     ALLOCATE(l_low_star(1:snapnum))
     ALLOCATE(l_up_star(1:snapnum))
     ALLOCATE(lstar_avranmst(1:snapnum))
     ALLOCATE(lstar_objmst(1:snapnum))
     
     lambda_star=0.
     l_up_star=0.
     l_low_star=0.
     lstar_avranmst=0.
     lstar_objmst=0.
  END IF
  
  IF (findgam) THEN
     nlam = nlam + 1
     gamunit = nlam + 100
     
     ALLOCATE(lambda_gam(1:snapnum))
     ALLOCATE(l_low_gam(1:snapnum))
     ALLOCATE(l_up_gam(1:snapnum))
     ALLOCATE(lgam_avranmst(1:snapnum))
     ALLOCATE(lgam_objmst(1:snapnum))
     
     lambda_gam=0.
     l_up_gam=0.
     l_low_gam=0.
     lgam_avranmst=0.
     lgam_objmst=0.
  END IF
  
  IF (findlamln) THEN
     nlam = nlam + 1
     lamlnunit = nlam + 100
     
     ALLOCATE(lambda_ln(1:snapnum))
     ALLOCATE(l_low_ln(1:snapnum))
     ALLOCATE(l_up_ln(1:snapnum))
     ALLOCATE(lln_avranmst(1:snapnum))
     ALLOCATE(lln_objmst(1:snapnum))
     
     lambda_ln=0.
     l_up_ln=0.
     l_low_ln=0.
     lln_avranmst=0.
     lln_objmst=0.
  END IF
  
  IF(nlam==0) THEN
     write(6,*) 'All lambda methods set to FALSE; need at least one TRUE'
     stop
  END IF
  
!====================================
! NB if changing file units, remember to also change in find_lambda.
  
  write(6,*)"       Calculating lambda..."
  
  open(10,file=trim(newPath)//'/objm_'//thisproj//'.dat')
  !open(11,file=trim(newPath)//'/objescaped_lam_'//thisproj//'.dat')
!(only need this if you want to check distances of escaped object stars)

!Files for CDFs:  (only bother doing this in xy)
!  if(thisproj=='xy') then
  
  INQUIRE(file = TRIM(CDFPath), exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)
  IF (.NOT. dirExists) THEN
     WRITE(6,'(a)') "Creating new directory: '"//TRIM(CDFPath)//"'"
     CALL system('mkdir -p '//TRIM(CDFPath))
  END IF
  
! open file for MST of object stars
  OPEN(12,file=trim(CDFPath)//'/MSTedgeL_'//thisproj//'.dat',status='replace')
! open files for MSTs of randomly selected stars:
  do i=1,nCDF
     write(filei,'(I6)') i
     open(300+i,file=trim(CDFPath)//'/MSTedgeL_'//thisproj//'_'&
          & //trim(adjustl(filei))//'.dat',status='replace')
  end do
 ! end if

!Record any stars that fall on top of each other &
!need their separation changing for the mst.
!Using 'unit#' as these are both populated in same bit of code (mst.f90)
  unit1=4
  unit2=5
  open(unit1,file=TRIM(newPath)//'/sepchange_obj.dat',position='append')
  open(unit2,file=TRIM(newPath)//'/sepchange.dat',position='append')
  write(unit1,*) "**** ",thisproj," ****"
  write(unit2,*) "**** ",thisproj," ****"

  do i=1,snapnum
     call find_lambda(i,nstars(i))
  end do

  close(unit1)
  close(unit2)

  write(6,*)"       ...done"
  
  close(10)
  !close(11)

!close CDF files:
  !if(thisproj=='xy') then
     close(12)
     do i=1,nCDF
        close(300+i)
     end do
  !end if

  
!**********************
! Lambda outputs
!**********************

  IF (findlam) THEN
!lamunit
! Median random MST edge lengths &
! object edge lengths used for each lambda, and lambda values with errors:
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lambda_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamunit,99) i,l_avranmst(i),l_objmst(i),lambda(i),l_low(i),l_up(i)
     END DO

     CLOSE(lamunit)
  END IF
     
  IF (findlambar) THEN
     OPEN(lambarunit,file=TRIM(newPath)//'/lambda/MST_lambar_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lambarunit,99) i,lbar_avranmst(i),lbar_objmst(i), &
             & lambda_bar(i),l_low_bar(i),l_up_bar(i)
     END DO
     
     CLOSE(lambarunit)
  END IF
  
  IF (findlamrms) THEN
     OPEN(lamrmsunit,file=TRIM(newPath)//'/lambda/MST_lamrms_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamrmsunit,99) i,lrms_avranmst(i),lrms_objmst(i), &
             & lambda_rms(i),l_low_rms(i),l_up_rms(i)
     END DO
     
     CLOSE(lamrmsunit)
  END IF
  
  IF (findlamsmr) THEN
     OPEN(lamsmrunit,file=TRIM(newPath)//'/lambda/MST_lamsmr_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamsmrunit,99) i,lsmr_avranmst(i),lsmr_objmst(i), &
             & lambda_smr(i),l_low_smr(i),l_up_smr(i)
     END DO
     
     CLOSE(lamsmrunit)
  END IF
  
  IF (findlamhar) THEN
     OPEN(lamharunit,file=TRIM(newPath)//'/lambda/MST_lamhar_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamharunit,99) i,lhar_avranmst(i),lhar_objmst(i), &
             & lambda_har(i),l_low_har(i),l_up_har(i)
     END DO
     
     CLOSE(lamharunit)
  END IF
  
  IF (findlamtil) THEN
     OPEN(lamtilunit,file=TRIM(newPath)//'/lambda/MST_lamtil_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamtilunit,99) i,ltil_avranmst(i),ltil_objmst(i), &
             & lambda_til(i),l_low_til(i),l_up_til(i)
     END DO
     
     CLOSE(lamtilunit)
  END IF
  
  IF (findlamNmed) THEN
     OPEN(lamNmedunit,file=TRIM(newPath)//'/lambda/MST_lamNmed_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamNmedunit,99) i,lNmed_avranmst(i),lNmed_objmst(i), &
             & lambda_Nmed(i),l_low_Nmed(i),l_up_Nmed(i)
     END DO
     
     CLOSE(lamNmedunit)
  END IF
  
  IF (findlamstar) THEN
     OPEN(lamstarunit,file=TRIM(newPath)//'/lambda/MST_lamstar_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamstarunit,99) i,lstar_avranmst(i),lstar_objmst(i), &
             & lambda_star(i),l_low_star(i),l_up_star(i)
     END DO
     
     CLOSE(lamstarunit)
  END IF
  
  IF (findgam) THEN
     OPEN(gamunit,file=TRIM(newPath)//'/lambda/MST_gam_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(gamunit,99) i,lgam_avranmst(i),lgam_objmst(i), &
             & lambda_gam(i),l_low_gam(i),l_up_gam(i)
     END DO
     
     CLOSE(gamunit)
  END IF
  
  IF (findlamln) THEN
     OPEN(lamlnunit,file=TRIM(newPath)//'/lambda/MST_lamln_'//thisproj//'.dat',status='replace')
     DO i=1,snapnum
        WRITE(lamlnunit,99) i,lln_avranmst(i),lln_objmst(i), &
             & lambda_ln(i),l_low_ln(i),l_up_ln(i)
     END DO
     
     CLOSE(lamlnunit)
  END IF
  
99 FORMAT(1X,I4,2(2X,F10.5),3(2X,F9.3))


!===========================================
  deallocate(obj_mass)
  deallocate(edgelengths)
  
  IF (findlam) THEN
     deALLOCATE(lambda)
     deALLOCATE(l_low)
     deALLOCATE(l_up)
     deallocate(l_avranmst)
     deallocate(l_objmst)
  END IF
  
  IF (findlambar) THEN
     deALLOCATE(lambda_bar)
     deALLOCATE(l_low_bar)
     deALLOCATE(l_up_bar)
     deallocate(lbar_avranmst)
     deallocate(lbar_objmst)
  END IF
  
  IF (findlamrms) THEN
     deALLOCATE(lambda_rms)
     deALLOCATE(l_low_rms)
     deALLOCATE(l_up_rms)
     deallocate(lrms_avranmst)
     deallocate(lrms_objmst)
  END IF
  
  IF (findlamsmr) THEN
     deALLOCATE(lambda_smr)
     deALLOCATE(l_low_smr)
     deALLOCATE(l_up_smr)
     deallocate(lsmr_avranmst)
     deallocate(lsmr_objmst)
  END IF
  
  IF (findlamhar) THEN
     deALLOCATE(lambda_har)
     deALLOCATE(l_low_har)
     deALLOCATE(l_up_har)
     deallocate(lhar_avranmst)
     deallocate(lhar_objmst)
  END IF
  
  IF (findlamtil) THEN
     deALLOCATE(lambda_til)
     deALLOCATE(l_low_til)
     deALLOCATE(l_up_til)
     deallocate(ltil_avranmst)
     deallocate(ltil_objmst)
  END IF
  
  IF (findlamNmed) THEN
     deALLOCATE(lambda_Nmed)
     deALLOCATE(l_low_Nmed)
     deALLOCATE(l_up_Nmed)
     deallocate(lNmed_avranmst)
     deallocate(lNmed_objmst)
  END IF
  
  IF (findlamstar) THEN
     deALLOCATE(lambda_star)
     deALLOCATE(l_low_star)
     deALLOCATE(l_up_star)
     deallocate(lstar_avranmst)
     deallocate(lstar_objmst)
  END IF
  
  IF (findgam) THEN
     deALLOCATE(lambda_gam)
     deALLOCATE(l_low_gam)
     deALLOCATE(l_up_gam)
     deallocate(lgam_avranmst)
     deallocate(lgam_objmst)
  END IF
  
  IF (findlamln) THEN
     deALLOCATE(lambda_ln)
     deALLOCATE(l_low_ln)
     deALLOCATE(l_up_ln)
     deallocate(lln_avranmst)
     deallocate(lln_objmst)
  END IF
  
end subroutine lambda_setup
