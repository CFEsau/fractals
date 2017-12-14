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
  DOUBLE PRECISION :: mi(1:ni)     !mass of star i
  DOUBLE PRECISION :: ri(1:3,1:ni) !position of star i
! ri_x, ri_y, ri_z = distance in x, y, z of star i from cluster centre
  DOUBLE PRECISION :: ri_x, ri_y, ri_z
  DOUBLE PRECISION :: rmag(1:ni) !distance of star i from cluster centre
  REAL :: ti                     !time of snapshot
  INTEGER :: IDs(1:ni)           !ID numbers for heapsort
  LOGICAL :: i_incluster(1:ni)   !is star i in the cluster?
                                 !(for when outliers are being ignored)

!-----------
! MST stuff
!-----------
  INTEGER :: nedge !number of edge lengths
  INTEGER :: nint_nedged2, ceint_nedged2 !NINT(... & CEILING(  REAL(nedge)/2.)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: edgeL !length of each MST edge
  INTEGER, DIMENSION(:), ALLOCATABLE :: edgeLlist !IDs of 'edgeL' entries
  double precision, dimension(:,:), allocatable :: connections
  double precision :: totallength  !total length of mst (sum of 'edgeL')
  double precision :: medianlength !median edge length in MST
  DOUBLE PRECISION :: obj_mst !Length of object tree (e.g. massive stars)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: obj_r !positions of obj stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,y,z !split of obj_r

!total lengths of random MSTs for different lambda:
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: l_ranmst, lbar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lrms_ranmst, lsmr_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lhar_ranmst, ltil_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lNmed_ranmst, lstar_ranmst
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lgam_ranmst, lln_ranmst

  INTEGER, DIMENSION(:), ALLOCATABLE :: rand_list !IDs of _randmst list

  LOGICAL :: done(1:ni) !record which stars have been randomly selected
  REAL :: rand
  
  INTEGER :: i,j,k !generic counters
  character(len=4) :: snapchar !snapi as string with leading zeroes

!Frequently used expressions:
  nedge=nmst-1
  !for median calculations:
  nint_nedged2=nint(real(nedge)/2.) !nint to round to nearest int from ~x.0
  ceint_nedged2=ceiling(real(nedge)/2.) !ceint to round up from x.5)

!Allocate memory for arrays of length nmst:
  ALLOCATE(edgeL(1:nedge))
  ALLOCATE(edgeLlist(1:nedge))
! Coordinates of the edge connections: (xi,yi,zi,xj,yj,zj for each edge)
  ALLOCATE(connections(1:6,1:nedge))
  ALLOCATE(x(1:nmst))
  ALLOCATE(y(1:nmst))
  ALLOCATE(z(1:nmst))
  ALLOCATE(obj_r(1:3,1:nmst))
  
! And populate arrays 
  mi(1:ni)=mstar(1:ni,snapi)
  ri(1:3,1:ni)=rstar(1:3,1:ni,snapi)
  ti=tstar(1,snapi)
  
!====================================
!For when outliers are being ignored:
! first entry of obj_mass stores IDs of original most massive stars,
! second stores new most-massive stars
  IF (snapi==1) obj_mass(1,1:nmst)=0
!=====================================

!Find distance of each star from cluster centre (spatial distribution):
  DO i = 1, ni
     
     ri_x = rstar(1,i,snapi)
     ri_y = rstar(2,i,snapi)
     ri_z = rstar(3,i,snapi)

     IF (thisproj=='xy') THEN
        ri(3,i) = 0.          ! z coordinate is 0
        rmag(i) = SQRT(ri_x**2 + ri_y**2)
     ELSE IF (thisproj=='yz') THEN
        ri(1,i) = 0.          ! x coordinate is 0
        rmag(i) = SQRT(ri_y**2 + ri_z**2)
     ELSE IF (thisproj=='xz') THEN
        ri(2,i) = 0.          ! y coordinate is 0
        rmag(i) = SQRT(ri_x**2 + ri_z**2)
     ELSE IF (thisproj=='3D') THEN
        rmag(i) = SQRT(ri_x**2 + ri_y**2 + ri_z**2)
     END IF

     i_incluster(i) = incluster(i,snapi)
  END DO


!============================================================
!Sort stars in order of mass & make selection of 'obj' stars
!============================================================
!
  obj_r=0.
  IDs=0

!Assign IDs to stars for heapsort
  DO j=1,ni
     IDs(j)=j
  END DO

  CALL heapsort(ni,mi,IDs)

!select nmst most massive stars:
!(while checking whether star should be ignored)
  j=0 !j tracks with i but changes with each *attempted* iteration,
      ! even if i does not change

! write snapshot number as string for output directory
  write(snapchar,2) snapi
2 format(I4.4)
  
! write positions of object stars:
  open(1,file=trim(lampath)//'/coords/snap'//snapchar//'_objpositions_'&
       &//thisproj//'.dat' ,status='replace')
  
  DO i = 1,nmst
5    j=j+1
     k = IDs(ni+1-j) !heapsort orders from small to large - invert
     
     IF (.NOT. i_incluster(k)) THEN
        !WRITE(objescunit,'(1X,I4,1X,I5,1X,I5,1X,F7.3,1X,F8.3,1X,A1)') &
             !& snapi, i, k, mi(k), rmag(k), i_incluster(k)
        
! if star has escaped cluster, leave it out of the list
        GOTO 5 !don't use 'cycle' as we don't want i to inrease
     END IF
     
     obj_r(1,i) = ri(1,k) !x position
     obj_r(2,i) = ri(2,k) !y position
     obj_r(3,i) = ri(3,k) !z position
     write(1,'(3(2X,F10.4))') obj_r(1,i),obj_r(2,i),obj_r(3,i)
     
! track masses selected after each snapshot for output file:
     obj_mass(2,i) = mi(IDs(ni+1-j))
  END DO
  
  close(1) !close 'objpositions'
  
! write snapshot, time, and list of object masses to output file
! when IDs of most massive stars change
  DO i=1, nmst
     IF(obj_mass(1,i) /= obj_mass(2,i)) THEN
        WRITE(objmunit,110) snapi, ti, obj_mass(2,1:nmst)
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

  edgeL = 0.

!Split obj_r across x,y,z arrays. Just helps readability in mst subroutine.
  x(1:nmst)=obj_r(1,1:nmst)
  y(1:nmst)=obj_r(2,1:nmst)
  z(1:nmst)=obj_r(3,1:nmst)

!set unit for output to file (separation data, used in mst.f90)
  fileunit=sepunit1
  CALL mst(snapi,nmst,x,y,z,edgeL,connections)
  
! write out coordinates of all edge connections:
  open(1,file=trim(lampath)//'/coords/snap'//snapchar//'_objconnections_'&
       &//thisproj//'.dat',status='replace')
  do i=1,nedge
     write(1,'(6(2X,F10.4))') connections(1:6,i)
  end do
  close(1)
  
!Assign IDs for heapsort to order MST edges
  DO j=1,nedge
     edgeLlist(j)=j
  END DO
  CALL heapsort(nedge,edgeL,edgeLlist)
  
! Write out the edge lengths of the MST to plot a CDF:
  do i=1,nedge
     edgelengths(i)=edgeL(edgeLlist(i))
  end do
  write(cdfobjunit,112) snapi,edgelengths(1:nedge)
112 FORMAT(1X,I4,*(2X,F9.5))
  
  
!######################################
!Average edge lengths for object stars
!######################################

  totallength=0.
  !DO i = 1,nedge
  !   totallength = totallength + edgeL(i)  !Add the edges of the mst
  !END DO                                    ! to find the total length
  totallength=SUM(edgeL)

  !Median edge length:
  medianlength=0.
  DO i = 1, nedge
     edgeLlist(i) = i
  END DO
  CALL heapsort(nedge, edgeL, edgeLlist)
  
  if (MOD(nedge,2)==0) then !even nedge (odd nmst), take mean of median two
     medianlength=edgeL(edgeLlist(nint_nedged2)) &
          + edgeL(edgeLlist(nint_nedged2 + 1)) !nedge/2=x.0, ensure nearest
     medianlength=medianlength/2.
  else !.not. MOD(nedge,2)==0, odd nedge (even nmst), take median.
     medianlength=edgeL(edgeLlist(ceint_nedged2)) !nedge/2=x.5, round up
  end if
  
  
!Lambda MST:
! This is just the same as lambar without 1/n
  IF (findlam) l_objmst(snapi) = totallength


!Lambda bar MST:
!Use the mean MST edge length
  IF (findlambar) lbar_objmst(snapi) = totallength/(nedge)


!Generalised mean: ((1/n)*SUM(x**p))**(1/p). Arithmetic mean p=1
!rms p=2, smr p=1/2, harmonic p=-1
  
!Lambda rms MST:
  IF (findlamrms) THEN
     DO i=1,nedge
        lrms_objmst(snapi) = lrms_objmst(snapi) + edgeL(i)**2.
     END DO
  END IF
  lrms_objmst(snapi) = ( (1./REAL(nedge))*lrms_objmst(snapi) )**0.5

  
!Lambda smr MST:
  IF (findlamsmr) THEN
     DO i=1,nedge
        lsmr_objmst(snapi) = lsmr_objmst(snapi) + (edgeL(i))**0.5
     END DO
  END IF
  lsmr_objmst(snapi) = ( (1./REAL(nedge))*(lsmr_objmst(snapi)) )**2.

  
!Lambda har MST:
  IF (findlamhar) THEN
     DO i=1,nedge
        lhar_objmst(snapi) = lhar_objmst(snapi) + (edgeL(i))**(-1.)
     END DO
  END IF
  lhar_objmst(snapi) = ( (1./REAL(nedge))*(lhar_objmst(snapi)) )**(-1.)
  
  
!Lambda tilde MST:
  ! Use the median edge length
  IF (findlamtil) ltil_objmst(snapi) = medianlength
  
  
!Lambda N-median MST:
  IF (findlamNmed) THEN
     if (mod(Nmed,2)==0) then !even edge lengths (odd nmst), round to nint
        lNmed_objmst(snapi)=medianlength*2 !i=1. *2 for l1+l2 as /2 for median.
        do i=2,Nmed/2
           !take values either side of two median:
           lNmed_objmst(snapi)=lNmed_objmst(snapi) &
                + edgeL(edgeLlist(nint_nedged2 - (i-1))) &
                + edgeL(edgeLlist(nint_nedged2 + i))
        end do
     else !odd edge lengths (even nmst), round to ceiling
        lNmed_objmst(snapi)=medianlength !i=1
        do i=2,(Nmed+1)/2        !take values either side of median:
           lNmed_objmst(snapi)=lNmed_objmst(snapi) &
                + edgeL(edgeLlist(ceint_nedged2 - (i-1))) &
                + edgeL(edgeLlist(ceint_nedged2 + (i-1)))
        end do
     end if
     lNmed_objmst(snapi)=lNmed_objmst(snapi)/real(Nmed)
  END IF
  
  
!Lambda star MST:
! Find the median length in the tree, & find length of a tree made from these
  IF (findlamstar) THEN
     lstar_objmst(snapi) = (nedge)*medianlength
     
! Then add on the actual length of the tree
     lstar_objmst(snapi) = lstar_objmst(snapi)+totallength
  END IF


!Generalised f-mean: f**(-1) * ((1/n)*SUM(f(x)) where f is some function
!Gamma is a special case where f=ln, so f**(-1)=exp

!Gamma MST:
  IF (findgam) THEN
     DO i = 1,nedge
        lgam_objmst(snapi) = lgam_objmst(snapi) + LOG(edgeL(i))  !Add edges to
     END DO                                                  !find total length
!Calculate the geometric mean
     lgam_objmst(snapi) = EXP( (1./REAL(nedge)) * lgam_objmst(snapi))
  END IF

  
!Lambda ln MST:
  IF (findlamln) THEN
     DO i = 1,nedge
        lln_objmst(snapi) = lln_objmst(snapi) + EXP(edgeL(i))  !Add edges to
     END DO                                                !find total length
!Calculate the geometric mean
     lln_objmst(snapi) = LOG( (1./REAL(nedge)) * lln_objmst(snapi))
  END IF


!          *********************************
!*****************************************************
!    Find the MST length for the nloop random MSTs
!*****************************************************
!          *********************************

  edgeL = 0.
  
  !set lengths of random MSTs to 0:
  IF (findlam) then
     ALLOCATE(l_ranmst(1:nloop))
     l_ranmst = 0. 
  END IF
  IF (findlambar) THEN
     ALLOCATE(lbar_ranmst(1:nloop))
     lbar_ranmst = 0.
  END IF
  IF (findlamrms) THEN
     ALLOCATE(lrms_ranmst(1:nloop))
     lrms_ranmst = 0.
  END IF
  IF (findlamsmr) THEN
     ALLOCATE(lsmr_ranmst(1:nloop))
     lsmr_ranmst = 0.
  END IF
  IF (findlamhar) THEN
     ALLOCATE(lhar_ranmst(1:nloop))
     lhar_ranmst = 0.
  END IF
  IF (findlamtil) THEN
     ALLOCATE(ltil_ranmst(1:nloop))
     ltil_ranmst = 0.
  END IF
  IF (findlamNmed) THEN
     ALLOCATE(lNmed_ranmst(1:nloop))
     lNmed_ranmst = 0.
  END IF
  IF (findlamstar) THEN
     ALLOCATE(lstar_ranmst(1:nloop))
     lstar_ranmst = 0.
  END IF
  IF (findgam) THEN
     ALLOCATE(lgam_ranmst(1:nloop))
     lgam_ranmst = 0.
  END IF
  IF (findlamln) THEN
     ALLOCATE(lln_ranmst(1:nloop))
     lln_ranmst = 0.
  END IF
  
  ALLOCATE(rand_list(1:nloop))
  
!set unit for output to file (separation data, used in mst.f90)
  fileunit=sepunit2
  
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
        
        x(i) = ri(1,k)
        y(i) = ri(2,k)
        z(i) = ri(3,k)
     END DO
     
     CALL mst(snapi,nmst,x,y,z,edgeL,connections)
     
!######################################
!Average edge lengths for random stars
!######################################

     totallength=0.
     !DO i = 1,nedge
     !   totallength = totallength + edgeL(i)  !Add the edges of the mst
     !END DO                                    ! to find the total length
     totallength=SUM(edgeL)

     medianlength=0.
     DO i = 1, nedge
        edgeLlist(i) = i
     END DO
     CALL heapsort(nedge, edgeL, edgeLlist)
     
!CDFs of random stars:
     !if(thisproj=='3D') then
        if (j.le.nCDF) then
! Write out the edge lengths of the MST to plot a CDF:
           do k=1,nedge
              edgelengths(k)=edgeL(edgeLlist(k))
           end do
           write(cdfranunit+j,300) snapi,edgelengths(1:nedge)
300        FORMAT(1X,I4,*(2X,F9.5))
        end if
     !end if

        !Median edge length of this MST (depends on even/odd nedge):
        if (MOD(nedge,2)==0) then
           ! if even no. of edge lengths, take mean of median 2.
           medianlength=edgeL(edgeLlist(nint_nedged2 +1)) &
                + edgeL(edgeLlist(nint_nedged2))
           medianlength=medianlength/2.
        else !.not. MOD(nedge,2)==0
           !if odd no. of edge lengths (even nmst), take median edge length.
           !(Use CEILING as always need to round up from #.5)
           medianlength=edgeL(edgeLlist(ceint_nedged2))
        end if
     
!Lambda MST:
! This is just the same as lambar without 1/n
     IF (findlam) l_ranmst(j) = totallength


!Lambda bar MST:
!Use the mean MST edge length
     IF (findlambar) lbar_ranmst(j) = totallength/(nedge)


!Generalised mean: ((1/n)*SUM(x**p))**(1/p). Arithmetic mean p=1
!rms p=2, smr p=1/2, harmonic p=-1
  
!Lambda rms MST:
     IF (findlamrms) THEN
        DO i=1,nedge
           lrms_ranmst(j) = lrms_ranmst(j) + (edgeL(i))**2.
        END DO
     END IF
     lrms_ranmst(j) = ( (1./REAL(nedge))*(lrms_ranmst(j)) )**0.5

  
!Lambda smr MST:
     IF (findlamsmr) THEN
        DO i=1,nedge
           lsmr_ranmst(j) = lsmr_ranmst(j) + (edgeL(i))**0.5
        END DO
     END IF
     lsmr_ranmst(j) = ( (1./REAL(nedge))*(lsmr_ranmst(j)) )**2.

  
!Lambda har MST:
     IF (findlamhar) THEN
        DO i=1,nedge
           lhar_ranmst(j) = lhar_ranmst(j) + (edgeL(i))**(-1.)
        END DO
     END IF
     lhar_ranmst(j) = ( (1./REAL(nedge))*(lhar_ranmst(j)) )**(-1.)


!Lambda tilde MST:
! Use the median edge length
     IF (findlamtil) ltil_ranmst(j) = medianlength


!Lambda N-median MST:
!Find the N median points and take the mean of these
     IF (findlamNmed) THEN
        if (mod(Nmed,2)==0) then !even nedge (odd nmst)
           lNmed_ranmst(j) = medianlength*2 !i=1. *2 as /2 earlier.
           do i=2,Nmed/2
              !take values either side of two used for median:
              lNmed_ranmst(j)=lNmed_ranmst(j) &
                   + edgeL(edgeLlist(nint_nedged2 - (i-1))) &
                   + edgeL(edgeLlist(nint_nedged2 + i))
           end do
        else !.not. mod(Nmed,2)==0, odd nedge (even nmst)
           lNmed_ranmst(j) = medianlength !i=1
           do i=2,(Nmed+1)/2
              !take values either side of median:
              lNmed_ranmst(j)=lNmed_ranmst(j) &
                   + edgeL(edgeLlist(ceint_nedged2 - (i-1))) &
                   + edgeL(edgeLlist(ceint_nedged2 + (i-1)))
           end do
        end if
        lNmed_ranmst(j)=lNmed_ranmst(j)/real(Nmed)
     END IF


!Lambda star MST:
! Find the median length in the tree, & find length of a tree made from these
     IF (findlamstar) THEN
        lstar_ranmst(j) = (nedge) * medianlength
        
! Then add on the actual length of the tree
        lstar_ranmst(j) = lstar_ranmst(j) + totallength
     END IF


!Generalised f-mean: f**(-1) * ((1/n)*SUM(f(x)) where f is some function
!Gamma is a special case where f=ln, so f**(-1)=exp

!Gamma MST:
     IF (findgam) THEN
        DO i = 1,nedge
           lgam_ranmst(j) = lgam_ranmst(j) + LOG(edgeL(i))  !Add edges to
        END DO                                               !find total length
        lgam_ranmst(j) = EXP( (1./REAL(nedge)) * lgam_ranmst(j) )
     END IF

  
!Lambda ln MST:
     IF (findlamln) THEN
        DO i = 1,nedge
           lln_ranmst(j) = lln_ranmst(j) + EXP(edgeL(i))  !Add edges to
        END DO                                                !find total length
!Calculate the geometric mean
        lln_ranmst(j) = LOG( (1./REAL(nedge)) * lln_ranmst(j))
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
          & l_allmsts(1:nloop,snapi), l_low(snapi), l_up(snapi), &
          & l_avranmst(snapi))
  END IF
  
  
!Calculate lambda bar:
  IF (findlambar) THEN
     CALL calc_lambda(lbar_objmst(snapi), lbar_ranmst, lambda_bar(snapi), &
          & lbar_allmsts(1:nloop,snapi), l_low_bar(snapi), l_up_bar(snapi), &
          & lbar_avranmst(snapi))
  END IF
  
  
!Calculate lambda rms:
  IF (findlamrms) THEN
     CALL calc_lambda(lrms_objmst(snapi), lrms_ranmst, lambda_rms(snapi), &
          & lrms_allmsts(1:nloop,snapi), l_low_rms(snapi), l_up_rms(snapi), &
          & lrms_avranmst(snapi))
  END IF
  
  
!Calculate lambda smr:
  IF (findlamsmr) THEN
     CALL calc_lambda(lsmr_objmst(snapi), lsmr_ranmst, lambda_smr(snapi), &
          & lsmr_allmsts(1:nloop,snapi), l_low_smr(snapi), l_up_smr(snapi), &
          & lsmr_avranmst(snapi))
  END IF
  
  
!Calculate lambda har:
  IF (findlamhar) THEN
     CALL calc_lambda(lhar_objmst(snapi), lhar_ranmst, lambda_har(snapi), &
          & lhar_allmsts(1:nloop,snapi), l_low_har(snapi), l_up_har(snapi), &
          & lhar_avranmst(snapi))
  END IF
  
  
!Calculate lambda tilde:
  IF (findlamtil) THEN
     CALL calc_lambda(ltil_objmst(snapi), ltil_ranmst, lambda_til(snapi), &
          & ltil_allmsts(1:nloop,snapi), l_low_til(snapi), l_up_til(snapi), &
          & ltil_avranmst(snapi))
  END IF
  
  
!Calculate lambda N median:
  IF (findlamNmed) THEN
     CALL calc_lambda(lNmed_objmst(snapi), lNmed_ranmst, lambda_Nmed(snapi), &
          & lNmed_allmsts(1:nloop,snapi), l_low_Nmed(snapi), l_up_Nmed(snapi), &
          & lNmed_avranmst(snapi))
  END IF
  
  
!Calculate lambda star:
  IF (findlamstar) THEN
     CALL calc_lambda(lstar_objmst(snapi), lstar_ranmst, lambda_star(snapi), &
          & lstar_allmsts(1:nloop,snapi), l_low_star(snapi), l_up_star(snapi), &
          & lstar_avranmst(snapi))
  END IF
  
  
!Calculate gamma:
  IF (findgam) THEN
     CALL calc_lambda(lgam_objmst(snapi), lgam_ranmst, lambda_gam(snapi), &
          & lgam_allmsts(1:nloop,snapi), l_low_gam(snapi), l_up_gam(snapi), &
          & lgam_avranmst(snapi))
  END IF
  
  
!Calculate lambda ln:
  IF (findlamln) THEN
     CALL calc_lambda(lln_objmst(snapi), lln_ranmst, lambda_ln(snapi), &
          & lln_allmsts(1:nloop,snapi), l_low_ln(snapi), l_up_ln(snapi), &
          & lln_avranmst(snapi))
  END IF
  
  
!Deallocate arrays:
  DEALLOCATE(edgeL)
  DEALLOCATE(edgeLlist)
  deallocate(connections)
  DEALLOCATE(x)
  DEALLOCATE(y)
  DEALLOCATE(z)
  DEALLOCATE(obj_r)
  IF (findlam) DEALLOCATE(l_ranmst)
  IF (findlambar) DEALLOCATE(lbar_ranmst)
  IF (findlamrms) DEALLOCATE(lrms_ranmst)
  IF (findlamsmr) DEALLOCATE(lsmr_ranmst)
  IF (findlamhar) DEALLOCATE(lhar_ranmst)
  IF (findlamtil) DEALLOCATE(ltil_ranmst)
  IF (findlamNmed) DEALLOCATE(lNmed_ranmst)
  IF (findlamstar) DEALLOCATE(lstar_ranmst)
  IF (findgam) DEALLOCATE(lgam_ranmst)
  IF (findlamln) DEALLOCATE(lln_ranmst)
  DEALLOCATE(rand_list)

END SUBROUTINE find_lambda
