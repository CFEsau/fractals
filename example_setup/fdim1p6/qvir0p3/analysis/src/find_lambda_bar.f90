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

SUBROUTINE find_lambda_bar(snapshoti,ni)
  USE parameters_module

  IMPLICIT NONE

  INTEGER, INTENT(in) :: snapshoti, ni !Snapshot number and star number
! mi,ri= mass, position & distance from com of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi, rcom_i
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
  REAL :: ti !time of snapshot
  INTEGER, DIMENSION(:), ALLOCATABLE :: mlist !ID numbers for heapsort
  LOGICAL, DIMENSION(:), ALLOCATABLE :: i_incluster !has this star escaped?

  INTEGER :: ndim !2D or 3D analysis
  
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: length !length of each MST edge
  DOUBLE PRECISION :: obj_mst !Length of object tree (e.g. most massive stars)
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: obj_r !positions of obj stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,y,z !split of obj_r
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: node !node connections for mst

  INTEGER :: nloop !number of random MSTs in loop
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rand_mst !total lengths of random MSTs
  INTEGER, DIMENSION(:), ALLOCATABLE :: rand_list !IDs of rand_mst list
  REAL :: avranmst !Average length of random trees
  REAL ::ranup, ranlow !Upper & lower boundaries, 1 sigma

  LOGICAL, DIMENSION(:), ALLOCATABLE :: done !record which stars have been randomly selected
  REAL :: rand
  
  INTEGER :: i,j,k !generic counters

! Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(rcom_i(1:ni))
  ALLOCATE(mlist(1:ni))
  ALLOCATE(done(1:ni))
!Allocate memory for arrays of length nmst:
  ALLOCATE(length(1:nmst-1))
  ALLOCATE(x(1:nmst))
  ALLOCATE(y(1:nmst))
  ALLOCATE(z(1:nmst))
  ALLOCATE(node(1:nmst-1,1:2))
  ALLOCATE(obj_r(1:nmst,1:3))
  

! And fill arrays 
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  ti=t(snapshoti,1)

!====================================
!For when outliers are being ignored
!====================================
!
  ALLOCATE(i_incluster(1:ni))
! first entry of obj_mass stores IDs of original most massive stars,
! second stores new most-massive stars
  IF (snapshoti==1) obj_mass(1,1:nmst)=0
!=====================================


!2D or 3D analysis?
  ndim = 2
  IF (ndim==2) THEN
     rcom_i(1:ni)=r_com(snapshoti,1:ni,4)
     i_incluster(1:ni)=this_in2Dcluster(snapshoti,1:ni)
  ELSE IF (ndim==3) THEN
     rcom_i(1:ni)=r_com(snapshoti,1:ni,5)
     i_incluster(1:ni)=this_in3Dcluster(snapshoti,1:ni)
  END IF

!============================================================
!Sort stars in order of mass & make selection of 'obj' stars
!============================================================
!
  obj_r=0.
  mlist=0
 
!Assign IDs to stars 
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
     
! if star has escaped cluster, leave it out of the list
     IF (.NOT. i_incluster(k)) THEN
        !WRITE(21,*) i, k, mi(k), rcom_i(k), i_incluster(k)
        GOTO 5 !don't use 'cycle' as we don't want i to inrease
     END IF
     obj_r(i,1) = ri(k,1) !x position
     obj_r(i,2) = ri(k,2) !y position
!For 2D set z to 0
     IF (ndim == 2) THEN
        obj_r(i,3) = 0. !z position
     ELSE IF (ndim == 3) THEN
        obj_r(i,3) = ri(k,3) !z position 
     ELSE
        STOP 'ndim must be 2 or 3'
     END IF

! track masses selected after each snapshot for output file:
     obj_mass(2,i) = mi(mlist(ni+1-j))
  END DO

! write snapshot, time, and list of object masses to output file
! when IDs of most massive stars change
  !DO i=1, nmst
     !IF(obj_mass(1,i) /= obj_mass(2,i)) THEN
        !WRITE(20,99) snapshoti, ti, obj_mass(2,1:nmst)
        !EXIT
     !END IF
  !END DO
!99   FORMAT(1X,I4,2X,E9.3,*(2X,F8.3))

  obj_mass(1,1:nmst)=obj_mass(2,1:nmst)

!======================================================================
!Find the MST length for the OBJECT stars
!======================================================================
!
   length = 0.
   obj_mst = 0.

!Split obj_r across x,y,z arrays. Just helps readability in mst subroutine.
   x(1:nmst)=obj_r(1:nmst,1)
   y(1:nmst)=obj_r(1:nmst,2)
   z(1:nmst)=obj_r(1:nmst,3)

   CALL mst(nmst,x,y,z,node,length)
   
   DO i = 1,nmst-1
      obj_mst = obj_mst + length(i)  !Add the edges of the mst to find 
   END DO                            !the total length

!Use mean of MST instead of total length
  obj_mst = obj_mst / (nmst-1)

!======================================================================
!Then find nloop random MSTs
!======================================================================
!  
!For large nmst it can take a while to build the MST, so reduce
!the number of times we do the loop - higher nmst MSTs are less 
!stochstic anyway...

    nloop = 1000
    IF(nmst >= 100) nloop = 50
   
   ALLOCATE(rand_mst(1:nloop))   
   ALLOCATE(rand_list(1:nloop))
   rand_mst = 0. !set lengths of random MSTs to 0.

   DO j = 1,nloop          !Do nloop random MSTs
      x = 0. ; y = 0. ; z = 0.
      done = .FALSE.          
      DO i = 1,nmst        !Select nmst random stars
6        CALL RANDOM_NUMBER(rand) !random number between 0. and 1.
         k = NINT(rand*ni)
         IF (k == 0) GOTO 6    !There is no star with id=0
         IF (done(k)) GOTO 6    !Don't chose the same star twice
         IF (.NOT. i_incluster(k)) GOTO 6 !Don't use escaped stars
         done(k) = .TRUE.       !You're selected
         x(i) = ri(k,1)
         y(i) = ri(k,2)
         z(i) = ri(k,3)

!============================================
!For 2D set z to 0
         IF (ndim == 2) z(i) = 0.
!============================================
      END DO

      CALL mst(nmst,x,y,z,node,length)

      DO i = 1,nmst-1
         rand_mst(j) = rand_mst(j) + length(i)    !Add the edges of the mst
      END DO                                  !to find the total length

!Use mean of MST instead of total length
     rand_mst(j) = rand_mst(j) / (nmst-1)
   END DO

!======================================================================
!Finally - find the average value of the random MSTs and the significance
!======================================================================
!
!The average is found using the median value - this means that small
!values don't skew the value
!The signifincance is found using the 1/6 and 5/6 boundaries

   DO j = 1,nloop
      rand_list(j) = j
   END DO
   
   
!Sort random MST lengths:
   CALL heapsort(nloop,rand_mst,rand_list)
   
   avranmst = rand_mst(rand_list(NINT(REAL(nloop)/2.)))
   ranlow = rand_mst(rand_list(NINT(REAL(nloop)/6.)))
   ranup = rand_mst(rand_list(NINT(5.*REAL(nloop)/6.)))

   lambda_bar(snapshoti) = avranmst/obj_mst
   l_low_bar(snapshoti) = ranlow/obj_mst
   l_up_bar(snapshoti) = ranup/obj_mst

  ! if ((snapshoti==1).or. (snapshoti==900))then
  !    print*,"avranmst",avranmst
  !    print*,"ranlow",ranlow
  !    print*,"ranup",ranup
  !    print*,"obj_mst",obj_mst
  ! end if
!Deallocate arrays:
   DEALLOCATE(mi)
   DEALLOCATE(ri)
   DEALLOCATE(i_incluster)
   DEALLOCATE(mlist)
   DEALLOCATE(done)
   DEALLOCATE(length)
   DEALLOCATE(obj_r)
   DEALLOCATE(x)
   DEALLOCATE(y)
   DEALLOCATE(z)
   DEALLOCATE(node)

 END SUBROUTINE find_lambda_bar