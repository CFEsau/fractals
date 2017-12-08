subroutine lambda_setup
! Allocate, zero, and deallocate various arrays needed for MST & lambda.
! File I/O for obj star masses, lambda values, edge lengths for CDFs, etc.

  USE sl_input_module
  use constants_module
  use parameters_module
  implicit none
  integer :: i,j
  character(len=6) :: filei !string version of 'i'
  integer :: nedge       ! number of edges in MST (nmst-1)
  integer :: nlam        ! number of lambda types
  integer :: lamunit     ! output unit
  character(len=5) :: Nmedstr
  LOGICAL :: dirExists
  
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
  nmed=3 !when using x number of median points, findlamNmed
         !odd for odd nedge, even nmst & v/v)
  
!For large nmst it can take a while to build the MST, so reduce
!the number of times we do the loop - higher nmst MSTs are less 
!stochastic anyway...
  nloop = 1000
  IF(nmst >= 100) nloop = 50

! keep track of file units for i/o
  sepunit1=4     ! make a note of stars with separation 0. for object
  sepunit2=5     ! and random stars, as separation is changed for MST
  objmunit=10    ! file unit for masses of object stars
  cdfobjunit=11  ! file unit for edge lengths in object MST
  cdfranunit=300 ! file unit SHIFT for nCDF output files (u+1 to u+nranCDF)
  lamunit=12     ! format unit: 98

! create directory for lambda data:
  lampath = trim(newPath)//'/lambda'
  INQUIRE(file = TRIM(lampath), exist = dirExists)
   IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//TRIM(lampath))
  END IF
! create directory for MST positions/edge coordinates:
  INQUIRE(file = TRIM(lampath)//'/coords', exist = dirExists)
   IF (.NOT. dirExists) THEN
     CALL system('mkdir -p '//TRIM(lampath)//'/coords')
  END IF
  
  nCDF = 20 !Number of CDFs plotted for random MSTs; must be =< nloop
  
  ! create directory for CDF data (edge lengths in MST for different subsets)
  CDFPath = TRIM(newPath)//'/CDFdata'
  INQUIRE(file = TRIM(CDFPath), exist = dirExists)
!(Works for gfortran. For ifort: ...directory=newDir,exist...)
  IF (.NOT. dirExists) THEN
     WRITE(6,'(a)') "Creating new directory: '"//TRIM(CDFPath)//"'"
     CALL system('mkdir -p '//TRIM(CDFPath))
  END IF
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
     
! final lambda value calculated (output)
     ALLOCATE(lambda(1:nsnaps))
! +ve error bar
     ALLOCATE(l_low(1:nsnaps))
! -ve error bar
     ALLOCATE(l_up(1:nsnaps))
! median average edge length of random MTSs
     ALLOCATE(l_avranmst(1:nsnaps))
! average edge length of object MST
     ALLOCATE(l_objmst(1:nsnaps))
! lambda values for all the random MSTs that are generated
     allocate(l_allmsts(1:nloop,1:nsnaps))

     lambda=0.
     l_up=0.
     l_low=0.
     l_avranmst=0.
     l_objmst=0.
     l_allmsts=0.
  END IF

  IF (findlambar) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_bar(1:nsnaps))
     ALLOCATE(l_low_bar(1:nsnaps))
     ALLOCATE(l_up_bar(1:nsnaps))
     ALLOCATE(lbar_avranmst(1:nsnaps))
     ALLOCATE(lbar_objmst(1:nsnaps))
     allocate(lbar_allmsts(1:nloop,1:nsnaps))
     
     lambda_bar=0.
     l_up_bar=0.
     l_low_bar=0.
     lbar_avranmst=0.
     lbar_objmst=0.
     lbar_allmsts=0.
  END IF

  IF (findlamrms) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_rms(1:nsnaps))
     ALLOCATE(l_low_rms(1:nsnaps))
     ALLOCATE(l_up_rms(1:nsnaps))
     ALLOCATE(lrms_avranmst(1:nsnaps))
     ALLOCATE(lrms_objmst(1:nsnaps))
     allocate(lrms_allmsts(1:nloop,1:nsnaps))

     lambda_rms=0.
     l_up_rms=0.
     l_low_rms=0.
     lrms_avranmst=0.
     lrms_objmst=0.
     lrms_allmsts=0.
  END IF

  IF (findlamsmr) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_smr(1:nsnaps))
     ALLOCATE(l_low_smr(1:nsnaps))
     ALLOCATE(l_up_smr(1:nsnaps))
     ALLOCATE(lsmr_avranmst(1:nsnaps))
     ALLOCATE(lsmr_objmst(1:nsnaps))
     allocate(lsmr_allmsts(1:nloop,1:nsnaps))

     lambda_smr=0.
     l_up_smr=0.
     l_low_smr=0.
     lsmr_avranmst=0.
     lsmr_objmst=0.
     lsmr_allmsts=0.
  END IF

  IF (findlamhar) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_har(1:nsnaps))
     ALLOCATE(l_low_har(1:nsnaps))
     ALLOCATE(l_up_har(1:nsnaps))
     ALLOCATE(lhar_avranmst(1:nsnaps))
     ALLOCATE(lhar_objmst(1:nsnaps))
     allocate(lhar_allmsts(1:nloop,1:nsnaps))

     lambda_har=0.
     l_up_har=0.
     l_low_har=0.
     lhar_avranmst=0.
     lhar_objmst=0.
     lhar_allmsts=0.
  END IF
  
  IF (findlamtil) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_til(1:nsnaps))
     ALLOCATE(l_low_til(1:nsnaps))
     ALLOCATE(l_up_til(1:nsnaps))
     ALLOCATE(ltil_avranmst(1:nsnaps))
     ALLOCATE(ltil_objmst(1:nsnaps))
     allocate(ltil_allmsts(1:nloop,1:nsnaps))
     
     lambda_til=0.
     l_up_til=0.
     l_low_til=0.
     ltil_avranmst=0.
     ltil_objmst=0.
     ltil_allmsts=0.
  END IF
  
  IF (findlamNmed) THEN
     
     !Nmed checks:
     nedge=nmst-1
     !------------
     ! Even nmst: (odd edges)
     if (.not. MOD(nedge,2)==0) then
        if (MOD(Nmed,2)==0) stop 'Nmed must be odd'
        if (.not. Nmed .gt. 1) stop 'Nmed must be greater than 1'
        if (Nmed .gt. nedge) stop 'Nmed must be less than nedge'
     else
     !-----------
     ! Odd nmst: (even edges)
        if (.not. MOD(Nmed,2)==0) stop 'Nmed must be even'
        if (.not. Nmed .gt. 2) stop 'Nmed must be greater than 2'
        if (Nmed .gt. nedge) stop 'Nmed must be less than nedge'
     end if
     
     write(Nmedstr,'(I5)') Nmed
     
     nlam = nlam + 1
     
     ALLOCATE(lambda_Nmed(1:nsnaps))
     ALLOCATE(l_low_Nmed(1:nsnaps))
     ALLOCATE(l_up_Nmed(1:nsnaps))
     ALLOCATE(lNmed_avranmst(1:nsnaps))
     ALLOCATE(lNmed_objmst(1:nsnaps))
     allocate(lNmed_allmsts(1:nloop,1:nsnaps))
     
     lambda_Nmed=0.
     l_up_Nmed=0.
     l_low_Nmed=0.
     lNmed_avranmst=0.
     lNmed_objmst=0.
     lNmed_allmsts=0.
  END IF
  
  IF (findlamstar) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_star(1:nsnaps))
     ALLOCATE(l_low_star(1:nsnaps))
     ALLOCATE(l_up_star(1:nsnaps))
     ALLOCATE(lstar_avranmst(1:nsnaps))
     ALLOCATE(lstar_objmst(1:nsnaps))
     allocate(lstar_allmsts(1:nloop,1:nsnaps))
     
     lambda_star=0.
     l_up_star=0.
     l_low_star=0.
     lstar_avranmst=0.
     lstar_objmst=0.
     lstar_allmsts=0.
  END IF
  
  IF (findgam) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_gam(1:nsnaps))
     ALLOCATE(l_low_gam(1:nsnaps))
     ALLOCATE(l_up_gam(1:nsnaps))
     ALLOCATE(lgam_avranmst(1:nsnaps))
     ALLOCATE(lgam_objmst(1:nsnaps))
     allocate(lgam_allmsts(1:nloop,1:nsnaps))
     
     lambda_gam=0.
     l_up_gam=0.
     l_low_gam=0.
     lgam_avranmst=0.
     lgam_objmst=0.
     lgam_allmsts=0.
  END IF
  
  IF (findlamln) THEN
     nlam = nlam + 1
     
     ALLOCATE(lambda_ln(1:nsnaps))
     ALLOCATE(l_low_ln(1:nsnaps))
     ALLOCATE(l_up_ln(1:nsnaps))
     ALLOCATE(lln_avranmst(1:nsnaps))
     ALLOCATE(lln_objmst(1:nsnaps))
     allocate(lln_allmsts(1:nloop,1:nsnaps))
     
     lambda_ln=0.
     l_up_ln=0.
     l_low_ln=0.
     lln_avranmst=0.
     lln_objmst=0.
     lln_allmsts=0.
  END IF
  
  IF(nlam==0) THEN
     write(6,*) 'All lambda methods set to FALSE; need at least one TRUE'
     stop
  END IF
  
!====================================
! NB if changing file units, remember to also change in find_lambda.
  
  write(6,*)"       Calculating lambda..."
  
  open(objmunit,file=trim(newPath)//'/objm_'//thisproj//'.dat')
  !open(objescunit,file=trim(newPath)//'/objescaped_lam_'//thisproj//'.dat')
  !(only need this if you want to check distances of escaped object stars)

!Files for CDFs:  (only bother doing this in xy)
!  if(thisproj=='xy') then

! open file for MST of object stars
  OPEN(cdfobjunit,file=trim(CDFPath)//'/MSTedgeL_'//thisproj//'.dat',&
       & status='replace')
! open files for MSTs of randomly selected stars:
  do i=1,nCDF
     write(filei,'(I6)') i
     open(cdfranunit+i,file=trim(CDFPath)//'/MSTedgeL_'//thisproj//'_'&
          & //trim(adjustl(filei))//'.dat',status='replace')
  end do
 ! end if

!Record any stars that fall on top of each other &
!need their separation changing for the mst.
!Using 'unit#' as these are both populated in same bit of code (mst.f90)
  !sepunit1=4
  !sepunit2=5
  open(sepunit1,file=TRIM(newPath)//'/sepchange_obj.dat',position='append')
  open(sepunit2,file=TRIM(newPath)//'/sepchange.dat',position='append')
  write(sepunit1,*) "**** ",thisproj," ****"
  write(sepunit2,*) "**** ",thisproj," ****"

  do i=1,nsnaps
     call find_lambda(i,nstars(i))
  end do

  close(sepunit1)
  close(sepunit2)

  write(6,*)"       ...done"
  
  close(objmunit)
  !close(objescunit)

!close CDF files:
  !if(thisproj=='xy') then
  close(cdfobjunit)
  do i=1,nCDF
     close(cdfranunit+i)
  end do
  !end if

  
!**********************
! Lambda outputs
!**********************
  
  
  IF (findlam) THEN
!lamunit
! Median random MST edge lengths &
! object edge lengths used for each lambda, and lambda values with errors:
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lambda_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,l_avranmst(i),l_objmst(i),lambda(i),l_low(i),l_up(i)
     END DO
     CLOSE(lamunit)
     
! Lambda CDF: lambda bar values for every random MST at each snapshot
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lam_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,l_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF
  
  
! Lambda bar: average MST, object MST, lambda bar value, and errors
  IF (findlambar) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lambar_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lbar_avranmst(i),lbar_objmst(i), &
             & lambda_bar(i),l_low_bar(i),l_up_bar(i)
     END DO
     CLOSE(lamunit)
     
! Lambda bar CDF: lambda bar values for every random MST at each snapshot
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lambar_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lbar_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF
  
  
  IF (findlamrms) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lamrms_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lrms_avranmst(i),lrms_objmst(i), &
             & lambda_rms(i),l_low_rms(i),l_up_rms(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lamrms_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lrms_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findlamsmr) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lamsmr_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lsmr_avranmst(i),lsmr_objmst(i), &
             & lambda_smr(i),l_low_smr(i),l_up_smr(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lamsmr_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lsmr_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findlamhar) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lamhar_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lhar_avranmst(i),lhar_objmst(i), &
             & lambda_har(i),l_low_har(i),l_up_har(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lamhar_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lhar_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findlamtil) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lamtil_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,ltil_avranmst(i),ltil_objmst(i), &
             & lambda_til(i),l_low_til(i),l_up_til(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lamtil_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,ltil_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findlamNmed) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lam'//&
          & trim(adjustl(Nmedstr))//'med_'//thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lNmed_avranmst(i),lNmed_objmst(i), &
             & lambda_Nmed(i),l_low_Nmed(i),l_up_Nmed(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lam_'//&
          & trim(adjustl(Nmedstr))//'med_'//thisproj//'.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lNmed_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findlamstar) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lamstar_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lstar_avranmst(i),lstar_objmst(i), &
             & lambda_star(i),l_low_star(i),l_up_star(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lamstar_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lstar_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findgam) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_gam_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lgam_avranmst(i),lgam_objmst(i), &
             & lambda_gam(i),l_low_gam(i),l_up_gam(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_gam_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lgam_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF

  
  IF (findlamln) THEN
     
     OPEN(lamunit,file=TRIM(newPath)//'/lambda/MST_lamln_'//&
          thisproj//'.dat',status='replace')
     DO i=1,nsnaps
        WRITE(lamunit,98) i,lln_avranmst(i),lln_objmst(i), &
             & lambda_ln(i),l_low_ln(i),l_up_ln(i)
     END DO
     CLOSE(lamunit)
     
     open(lamunit,file=TRIM(CDFpath)//'/allMSTs_lamln_'//thisproj//&
          & '.dat',status='replace')
     do i=1,nsnaps
        write(lamunit,99) i,lln_allmsts(1:nloop,i) 
     end do
     close(lamunit)
     
  END IF
  
98 FORMAT(1X,I4,2(2X,F10.5),3(2X,F9.3))
99 format(1X, I4, 1000(2X,F10.5))


!===========================================
  deallocate(obj_mass)
  deallocate(edgelengths)
  
  IF (findlam) THEN
     deALLOCATE(lambda)
     deALLOCATE(l_low)
     deALLOCATE(l_up)
     deallocate(l_avranmst)
     deallocate(l_objmst)
     deallocate(l_allmsts)
  END IF
  
  IF (findlambar) THEN
     deALLOCATE(lambda_bar)
     deALLOCATE(l_low_bar)
     deALLOCATE(l_up_bar)
     deallocate(lbar_avranmst)
     deallocate(lbar_objmst)
     deallocate(lbar_allmsts)
  END IF
  
  IF (findlamrms) THEN
     deALLOCATE(lambda_rms)
     deALLOCATE(l_low_rms)
     deALLOCATE(l_up_rms)
     deallocate(lrms_avranmst)
     deallocate(lrms_objmst)
     deallocate(lrms_allmsts)
  END IF
  
  IF (findlamsmr) THEN
     deALLOCATE(lambda_smr)
     deALLOCATE(l_low_smr)
     deALLOCATE(l_up_smr)
     deallocate(lsmr_avranmst)
     deallocate(lsmr_objmst)
     deallocate(lsmr_allmsts)
  END IF
  
  IF (findlamhar) THEN
     deALLOCATE(lambda_har)
     deALLOCATE(l_low_har)
     deALLOCATE(l_up_har)
     deallocate(lhar_avranmst)
     deallocate(lhar_objmst)
     deallocate(lhar_allmsts)
  END IF
  
  IF (findlamtil) THEN
     deALLOCATE(lambda_til)
     deALLOCATE(l_low_til)
     deALLOCATE(l_up_til)
     deallocate(ltil_avranmst)
     deallocate(ltil_objmst)
     deallocate(ltil_allmsts)
  END IF
  
  IF (findlamNmed) THEN
     deALLOCATE(lambda_Nmed)
     deALLOCATE(l_low_Nmed)
     deALLOCATE(l_up_Nmed)
     deallocate(lNmed_avranmst)
     deallocate(lNmed_objmst)
     deallocate(lNmed_allmsts)
  END IF
  
  IF (findlamstar) THEN
     deALLOCATE(lambda_star)
     deALLOCATE(l_low_star)
     deALLOCATE(l_up_star)
     deallocate(lstar_avranmst)
     deallocate(lstar_objmst)
     deallocate(lstar_allmsts)
  END IF
  
  IF (findgam) THEN
     deALLOCATE(lambda_gam)
     deALLOCATE(l_low_gam)
     deALLOCATE(l_up_gam)
     deallocate(lgam_avranmst)
     deallocate(lgam_objmst)
     deallocate(lgam_allmsts)
  END IF
  
  IF (findlamln) THEN
     deALLOCATE(lambda_ln)
     deALLOCATE(l_low_ln)
     deALLOCATE(l_up_ln)
     deallocate(lln_avranmst)
     deallocate(lln_objmst)
     deallocate(lln_allmsts)
  END IF
  
end subroutine lambda_setup
