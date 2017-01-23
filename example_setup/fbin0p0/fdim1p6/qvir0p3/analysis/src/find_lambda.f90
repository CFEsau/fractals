!***Modifications***
!-Ignore stars a certain distance from the com / beyond n*half-mass-radius

!This subroutine finds the value of mass segregation and its significance
!for the set of OBJECT STARS that you provide it.
!It works by calculating the MST for the OBJECT STARS, and then returning
!a value produced from sampling MSTs from randomly selected stars.
!This subroutine can be treated more-or-less like a black box.
!
!This algorithm is presented in Allison et al, 2009, MNRAS, 395, 1449
!

SUBROUTINE find_lambda(snapi,ni)
  USE parameters_module

  IMPLICIT NONE

!-----------
! star data
!-----------
  INTEGER, INTENT(in) :: snapi, ni !Snapshot number and star number
! mi = mass of star i
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
! ri = position of star i
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
! ri_x, ri_y, ri_z = distance in x, y, z of star i from c of m
  DOUBLE PRECISION :: ri_x, ri_y, ri_z
! rmag = distance of star i from centre of mass
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rmag
  REAL :: ti !time of snapshot
  INTEGER, DIMENSION(:), ALLOCATABLE :: mlist !ID numbers for heapsort
 !i_incluster = is star i still in the cluster?
  LOGICAL, DIMENSION(:), ALLOCATABLE :: i_incluster

!-----------
! MST stuff
!-----------
! length = length of each MST edge
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: length
! length_list = IDs of 'length' entries (for heapsort)
  INTEGER, DIMENSION(:), ALLOCATABLE :: length_list
  double precision :: totallength !total length of mst (sum of 'length')
  double precision :: medianlength !median edge length in MST
! Length of object tree (e.g. most massive stars) for different lambda:
  DOUBLE PRECISION :: obj_mst, obj_mst_bar, obj_mst_til, obj_mst_star
  DOUBLE PRECISION :: obj_mst_gam
! obj_r = positions of obj stars
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: obj_r
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,y,z !split of obj_r
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: node !node connections for mst

!total lengths of random MSTs for different lambda:
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: l_ranmst, lbar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lrms_ranmst, lsmr_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lhar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ltil_ranmst, lstar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lgam_ranmst, lln_ranmst

  INTEGER, DIMENSION(:), ALLOCATABLE :: rand_list !IDs of _randmst list

! done = record which stars have been randomly selected
  LOGICAL, DIMENSION(:), ALLOCATABLE :: done
  REAL :: rand
  
  INTEGER :: i,j,k !generic counters

! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(rmag(1:ni))
  ALLOCATE(mlist(1:ni))
  ALLOCATE(done(1:ni))
!Allocate memory for arrays of length nmst:
  ALLOCATE(length(1:nmst-1))
  ALLOCATE(length_list(1:nmst-1))
  ALLOCATE(x(1:nmst))
  ALLOCATE(y(1:nmst))
  ALLOCATE(z(1:nmst))
  ALLOCATE(node(1:nmst-1,1:2))
  ALLOCATE(obj_r(1:nmst,1:3))

! And fill arrays 
  mi(1:ni)=m(snapi,1:ni)
  ri(1:ni,1:3)=r(snapi,1:ni,1:3)
  ti=t(snapi,1)

!====================================
!For when outliers are being ignored
!====================================
!
  ALLOCATE(i_incluster(1:ni))
! first entry of obj_mass stores IDs of original most massive stars,
! second stores new most-massive stars
  IF (snapi==1) obj_mass(1,1:nmst)=0
!=====================================

!Find distance of each star from centre of mass:
  DO i = 1, ni

     ri_x = ri_com(snapi,i,1)
     ri_y = ri_com(snapi,i,2)
     ri_z = ri_com(snapi,i,3)

     IF (proj=='xy') THEN
        ri(i,3) = 0.
        rmag(i) = SQRT(ri_x**2 + ri_y**2)
     ELSE IF (proj=='yz') THEN
        ri(i,1) = 0.
        rmag(i) = SQRT(ri_y**2 + ri_z**2)
     ELSE IF (proj=='xz') THEN
        ri(i,2) = 0.
        rmag(i) = SQRT(ri_x**2 + ri_z**2)
     ELSE IF (proj=='3D') THEN
        rmag(i) = SQRT(ri_x**2 + ri_y**2 + ri_z**2)
     END IF

     i_incluster(i) = incluster(snapi,i)
  END DO


!============================================================
!Sort stars in order of mass & make selection of 'obj' stars
!============================================================
!
  obj_r=0.
  mlist=0

!Assign IDs to stars for heapsort
  DO j=1,ni
     mlist(j)=j
  END DO

  CALL heapsort(ni,mi,mlist)

!select nmst most massive stars:
!(while checking whether star should be ignored)
  j=0 !j tracks with i but changes with each *attempted* iteration,
      ! even if i does not change

  DO i = 1,nmst
5    j=j+1
     k = mlist(ni+1-j) !heapsort orders from small to large - invert

     IF (.NOT. i_incluster(k)) THEN
        !WRITE(11,'(1X,I4,1X,I5,1X,I5,1X,F7.3,1X,F8.3,1X,A1)') &
             !& snapi, i, k, mi(k), rmag(k), i_incluster(k)

! if star has escaped cluster, leave it out of the list
        GOTO 5 !don't use 'cycle' as we don't want i to inrease
     END IF

     obj_r(i,1) = ri(k,1) !x position
     obj_r(i,2) = ri(k,2) !y position
     obj_r(i,3) = ri(k,3) !z position

! track masses selected after each snapshot for output file:
     obj_mass(2,i) = mi(mlist(ni+1-j))
  END DO

! write snapshot, time, and list of object masses to output file
! when IDs of most massive stars change
  DO i=1, nmst
     IF(obj_mass(1,i) /= obj_mass(2,i)) THEN
        WRITE(10,110) snapi, ti, obj_mass(2,1:nmst)
        EXIT
     END IF
  END DO
110   FORMAT(1X,I4,2X,E9.3,*(2X,F8.3))

  obj_mass(1,1:nmst)=obj_mass(2,1:nmst)


!          *********************************
!*****************************************************
!      Find the MST length for the OBJECT stars
!*****************************************************
!          *********************************

  length = 0.

!Split obj_r across x,y,z arrays. Just helps readability in mst subroutine.
  x(1:nmst)=obj_r(1:nmst,1)
  y(1:nmst)=obj_r(1:nmst,2)
  z(1:nmst)=obj_r(1:nmst,3)

!set unit for output to file (separation data)
  fileunit=unit1
  CALL mst(snapi,nmst,x,y,z,node,length)

!Assign IDs for heapsort to order MST edges
  DO j=1,nmst-1
     length_list(j)=j
  END DO
  CALL heapsort(nmst-1,length,length_list)
  
! Write out the edge lengths of the MST to plot a CDF:
  do i=1,nmst-1
     edgelengths(i)=length(length_list(i))
  end do
  write(12,112) snapi,edgelengths(1:nmst-1)
112 FORMAT(1X,I4,*(2X,F9.5))

!######################################
!Average edge lengths for object stars
!######################################

  totallength=0.
  !DO i = 1,nmst-1
  !   totallength = totallength + length(i)  !Add the edges of the mst
  !END DO                                    ! to find the total length
  totallength=SUM(length)

  medianlength=0.
  DO i = 1, nmst-1
     length_list(i) = i
  END DO
  CALL heapsort(nmst-1, length, length_list)
  medianlength=length(length_list( NINT(REAL(nmst-1)/2.) ) )

!Lambda MST:
! This is just the same as lambar without 1/n
  IF (findlam) l_objmst(snapi) = totallength


!Lambda bar MST:
!Use the mean MST edge length
  IF (findlambar) lbar_objmst(snapi) = totallength/(nmst-1)


!Generalised mean: ((1/n)*SUM(x**p))**(1/p). Arithmetic mean p=1
!rms p=2, smr p=1/2, harmonic p=-1
  
!Lambda rms MST:
  IF (findlamrms) THEN
     DO i=1,nmst-1
        lrms_objmst(snapi) = lrms_objmst(snapi) + length(i)**2.
     END DO
  END IF
  lrms_objmst(snapi) = ( (1./REAL(nmst-1))*lrms_objmst(snapi) )**0.5

  
!Lambda smr MST:
  IF (findlamsmr) THEN
     DO i=1,nmst-1
        lsmr_objmst(snapi) = lsmr_objmst(snapi) + (length(i))**0.5
     END DO
  END IF
  lsmr_objmst(snapi) = ( (1./REAL(nmst-1))*(lsmr_objmst(snapi)) )**2.

  
!Lambda har MST:
  IF (findlamhar) THEN
     DO i=1,nmst-1
        lhar_objmst(snapi) = lhar_objmst(snapi) + (length(i))**(-1.)
     END DO
  END IF
  lhar_objmst(snapi) = ( (1./REAL(nmst-1))*(lhar_objmst(snapi)) )**(-1.)
  
  
!Lambda tilde MST:
  ! Use the median edge length
  IF (findlamtil) ltil_objmst(snapi) = medianlength


!Lambda star MST:
! Find the median length in the tree, & find length of a tree made from these
  IF (findlamstar) THEN
     lstar_objmst(snapi) = (nmst-1)*medianlength
     
! Then add on the actual length of the tree
     lstar_objmst(snapi) = lstar_objmst(snapi)+totallength
  END IF


!Generalised f-mean: f**(-1) * ((1/n)*SUM(f(x)) where f is some function
!Gamma is a special case where f=ln, so f**(-1)=exp
!Also try ln - take log of summed exponent (f=exp, f**(-1)=ln)

!Gamma MST:
  IF (findgam) THEN
     DO i = 1,nmst-1
        lgam_objmst(snapi) = lgam_objmst(snapi) + LOG(length(i))  !Add edges to
     END DO                                                  !find total length
!Calculate the geometric mean
     lgam_objmst(snapi) = EXP( (1./REAL(nmst-1)) * lgam_objmst(snapi))
  END IF

  
!Lambda ln MST:
  IF (findlamln) THEN
     DO i = 1,nmst-1
        lln_objmst(snapi) = lln_objmst(snapi) + EXP(length(i))  !Add edges to
     END DO                                                !find total length
!Calculate the geometric mean
     lln_objmst(snapi) = LOG( (1./REAL(nmst-1)) * lln_objmst(snapi))
  END IF


!          *********************************
!*****************************************************
!    Find the MST length for the nloop random MSTs
!*****************************************************
!          *********************************

  length = 0.
  
  ALLOCATE(l_ranmst(1:nloop))
  ALLOCATE(lbar_ranmst(1:nloop))
  ALLOCATE(lrms_ranmst(1:nloop))
  ALLOCATE(lsmr_ranmst(1:nloop))
  ALLOCATE(lhar_ranmst(1:nloop))
  ALLOCATE(ltil_ranmst(1:nloop))
  ALLOCATE(lstar_ranmst(1:nloop))
  ALLOCATE(lgam_ranmst(1:nloop))
  ALLOCATE(lln_ranmst(1:nloop))

  l_ranmst = 0. !set lengths of random MSTs to 0.
  lbar_ranmst = 0.
  lrms_ranmst = 0.
  lsmr_ranmst = 0.
  lhar_ranmst = 0.
  ltil_ranmst = 0.
  lstar_ranmst = 0.
  lgam_ranmst = 0.
  lln_ranmst = 0.
  
  ALLOCATE(rand_list(1:nloop))

!set unit for output to file (separation data)
  fileunit=unit2

  DO j = 1,nloop          !Do nloop random MSTs
     x = 0. ; y = 0. ; z = 0.
     done = .FALSE.          
     DO i = 1,nmst        !Select nmst random stars
6       CALL RANDOM_NUMBER(rand) !random number between 0. and 1.
        k = NINT(rand*ni)
        IF (k == 0) GOTO 6    !There is no star with id=0
        IF (done(k)) GOTO 6    !Don't chose the same star twice
        IF (.NOT. i_incluster(k)) GOTO 6 !Don't use escaped stars
        done(k) = .TRUE.       !You're selected
        x(i) = ri(k,1)
        y(i) = ri(k,2)
        z(i) = ri(k,3)
     END DO

     CALL mst(snapi,nmst,x,y,z,node,length)

!CDFs of random stars:
     if(proj=='xy') then
        if (j.le.nCDF) then
! Find edge lengths for nCDF random MSTs
! Need to heapsort mst edge lengths for cdf data:
           DO k=1,nmst-1
              length_list(k)=k
           END DO
           CALL heapsort(nmst-1,length,length_list)
! Write out the edge lengths of the MST to plot a CDF:
           do k=1,nmst-1
              edgelengths(k)=length(length_list(k))
           end do
           write(300+j,300) snapi,edgelengths(1:nmst-1)
300        FORMAT(1X,I4,*(2X,F9.5))
        end if
     end if

!######################################
!Average edge lengths for random stars
!######################################

     totallength=0.
     !DO i = 1,nmst-1
     !   totallength = totallength + length(i)  !Add the edges of the mst
     !END DO                                    ! to find the total length
     totallength=SUM(length)

     medianlength=0.
     DO i = 1, nmst-1
        length_list(i) = i
     END DO
     CALL heapsort(nmst-1, length, length_list)
     medianlength=length(length_list( NINT(REAL(nmst-1)/2.) ) )

!Lambda MST:
! This is just the same as lambar without 1/n
     IF (findlam) l_ranmst(j) = totallength


!Lambda bar MST:
!Use the mean MST edge length
     IF (findlambar) lbar_ranmst(j) = totallength/(nmst-1)


!Generalised mean: ((1/n)*SUM(x**p))**(1/p). Arithmetic mean p=1
!rms p=2, smr p=1/2, harmonic p=-1
  
!Lambda rms MST:
     IF (findlamrms) THEN
        DO i=1,nmst-1
           lrms_ranmst(j) = lrms_ranmst(j) + (length(i))**2.
        END DO
     END IF
     lrms_ranmst(j) = ( (1./REAL(nmst-1))*(lrms_ranmst(j)) )**0.5

  
!Lambda smr MST:
     IF (findlamsmr) THEN
        DO i=1,nmst-1
           lsmr_ranmst(j) = lsmr_ranmst(j) + (length(i))**0.5
        END DO
     END IF
     lsmr_ranmst(j) = ( (1./REAL(nmst-1))*(lsmr_ranmst(j)) )**2.

  
!Lambda har MST:
     IF (findlamhar) THEN
        DO i=1,nmst-1
           lhar_ranmst(j) = lhar_ranmst(j) + (length(i))**(-1.)
        END DO
     END IF
     lhar_ranmst(j) = ( (1./REAL(nmst-1))*(lhar_ranmst(j)) )**(-1.)


!Lambda tilde MST:
! Use the median edge length
     IF (findlamtil) ltil_ranmst(j) = medianlength


!Lambda star MST:
! Find the median length in the tree, & find length of a tree made from these
     IF (findlamstar) THEN
        lstar_ranmst(j) = (nmst-1) * medianlength
        
! Then add on the actual length of the tree
        lstar_ranmst(j) = lstar_ranmst(j) + totallength
     END IF


!Generalised f-mean: f**(-1) * ((1/n)*SUM(f(x)) where f is some function
!Gamma is a special case where f=ln, so f**(-1)=exp
!Also try ln - take log of summed exponent (f=exp, f**(-1)=ln)


!Gamma MST:
     IF (findgam) THEN
        DO i = 1,nmst-1
           lgam_ranmst(j) = lgam_ranmst(j) + LOG(length(i))  !Add edges to
        END DO                                               !find total length
        lgam_ranmst(j) = EXP( (1./REAL(nmst-1)) * lgam_ranmst(j) )
     END IF

  
!Lambda ln MST:
     IF (findlamln) THEN
        DO i = 1,nmst-1
           lln_ranmst(j) = lln_ranmst(j) + EXP(length(i))  !Add edges to
        END DO                                                !find total length
!Calculate the geometric mean
        lln_ranmst(j) = LOG( (1./REAL(nmst-1)) * lln_ranmst(j))
     END IF
     
  END DO !end of nloop


!========================================================================
!Finally - find the average value of the random MSTs and the significance
!========================================================================
!
!The average is found using the median value - this means that small
!values don't skew the value.
!The signifincance is found using the 1/6 and 5/6 boundaries


!Calculate lambda:
  IF (findlam) THEN
     CALL calc_lambda(l_objmst(snapi), l_ranmst, lambda(snapi), &
          & l_low(snapi), l_up(snapi), l_avranmst(snapi))
  END IF


!Calculate lambda bar:
  IF (findlambar) THEN
     CALL calc_lambda(lbar_objmst(snapi), lbar_ranmst, lambda_bar(snapi), &
          & l_low_bar(snapi), l_up_bar(snapi), lbar_avranmst(snapi))
  END IF


!Calculate lambda rms:
  IF (findlamrms) THEN
     CALL calc_lambda(lrms_objmst(snapi), lrms_ranmst, lambda_rms(snapi), &
          & l_low_rms(snapi), l_up_rms(snapi), lrms_avranmst(snapi))
  END IF


!Calculate lambda smr:
  IF (findlamsmr) THEN
     CALL calc_lambda(lsmr_objmst(snapi), lsmr_ranmst, lambda_smr(snapi), &
          & l_low_smr(snapi), l_up_smr(snapi), lsmr_avranmst(snapi))
  END IF


!Calculate lambda har:
  IF (findlamhar) THEN
     CALL calc_lambda(lhar_objmst(snapi), lhar_ranmst, lambda_har(snapi), &
          & l_low_har(snapi), l_up_har(snapi), lhar_avranmst(snapi))
  END IF


!Calculate lambda tilde:
  IF (findlamtil) THEN
     CALL calc_lambda(ltil_objmst(snapi), ltil_ranmst, lambda_til(snapi), &
          & l_low_til(snapi), l_up_til(snapi), ltil_avranmst(snapi))
  END IF


!Calculate lambda star:
  IF (findlamstar) THEN
     CALL calc_lambda(lstar_objmst(snapi), lstar_ranmst, lambda_star(snapi), &
          & l_low_star(snapi), l_up_star(snapi), lstar_avranmst(snapi))
  END IF


!Calculate gamma:
  IF (findgam) THEN
     CALL calc_lambda(lgam_objmst(snapi), lgam_ranmst, lambda_gam(snapi), &
          & l_low_gam(snapi), l_up_gam(snapi), lgam_avranmst(snapi))
  END IF


!Calculate lambda ln:
  IF (findlamln) THEN
     CALL calc_lambda(lln_objmst(snapi), lln_ranmst, lambda_ln(snapi), &
          & l_low_ln(snapi), l_up_ln(snapi), lln_avranmst(snapi))
  END IF


!Deallocate arrays:
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(rmag)
  DEALLOCATE(mlist)
  DEALLOCATE(done)
  DEALLOCATE(i_incluster)
  DEALLOCATE(length)
  DEALLOCATE(length_list)
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(node)
  DEALLOCATE(obj_r)
  IF (findlam) DEALLOCATE(l_ranmst)
  IF (findlambar) DEALLOCATE(lbar_ranmst)
  IF (findlamrms) DEALLOCATE(lrms_ranmst)
  IF (findlamsmr) DEALLOCATE(lsmr_ranmst)
  IF (findlamhar) DEALLOCATE(lhar_ranmst)
  IF (findlamtil) DEALLOCATE(ltil_ranmst)
  IF (findlamstar) DEALLOCATE(lstar_ranmst)
  IF (findgam) DEALLOCATE(lgam_ranmst)
  IF (findlamln) DEALLOCATE(lln_ranmst)
  DEALLOCATE(rand_list)

END SUBROUTINE find_lambda
