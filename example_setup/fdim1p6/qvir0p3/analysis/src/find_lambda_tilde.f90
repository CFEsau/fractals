!This subroutine finds the value of mass segregation and its significance
!for the set of OBJECT STARS that you provide it.
!It works by calculating the MST for the OBJECT STARS, and then returning
!a value produced from sampling MSTs from randomly selected stars.
!This subroutine can be treated more-or-less like a black box.
!
! This is the version from Maschberger & Clarke, 2011 ,MNRAS, 416, 541
! That uses the MEDIAN value of the MST edge lengths

SUBROUTINE find_lambda_tilde(snapshoti,ni)
  use parameters_module

  IMPLICIT NONE

  INTEGER, INTENT(in) :: snapshoti, ni !Snapshot number and star number
! mi,ri= mass, position & distance from com of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi, rcom_i
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri
  REAL :: ti !time of snapshot
  INTEGER, DIMENSION(:), ALLOCATABLE :: mlist !ID numbers for heapsort
  !LOGICAL, DIMENSION(:), ALLOCATABLE :: i_incluster !has this star escaped?

  INTEGER :: ndim !2D or 3D analysis
  
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: length !length of each MST edge
  INTEGER, DIMENSION(:), ALLOCATABLE :: length_list !IDs of 'length' entries (for heapsort)
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
  
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(rcom_i(1:ni))
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
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  ti=t(snapshoti,1)
  
!2D or 3D analysis?
  ndim = 2
  IF (ndim==2) THEN
     rcom_i(1:ni)=r_com(snapshoti,1:ni,4)
  ELSE IF (ndim==3) THEN
     rcom_i(1:ni)=r_com(snapshoti,1:ni,5)
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
  DO i = 1,nmst
     k = mlist(ni+1-i) !heapsort orders from small to large - invert
     
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

  END DO

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

! Find the median MST length
  DO i = 1, nmst-1
     length_list(i) = i
  end do

  call heapsort(nmst-1, length, length_list)

  obj_mst = length(length_list( NINT(REAL(nmst-1)/2.) ) )


!======================================================================
!Then find nloop random star MSTs
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
6       CALL RANDOM_NUMBER(rand)
        k = NINT(rand*ni)
        IF (k == 0.) GOTO 6     !There is no star with id=0
        IF (done(k)) GOTO 6   !Don't chose the same star twice
        done(k) = .TRUE.        !You're selected
        x(i) = ri(k,1)
        y(i) = ri(k,2)
        z(i) = ri(k,3)

!============================================
!For 2D set z to 0
        IF (ndim == 2) z(i) = 0.
!============================================
     END DO
     
     CALL mst(nmst,x,y,z,node,length)

! Find the median MST length
     DO i = 1, nmst-1
        length_list(i) = i
     END DO

     CALL heapsort(nmst-1, length, length_list)
      
     rand_mst(j) = length(length_list( NINT(REAL(nmst-1)/2.) ) )     
    
  END DO


!======================================================================
!Finally - find the average value of the random MSTs and the significance
!======================================================================
!
!The average is found using the median value - this means that small
!values don't skew the value
!The signifincance is found using the 1/6 and 2/6 boundaries

  DO j = 1,nloop
     rand_list(j) = j
  END DO
  
!Sort random MST lengths:
  CALL heapsort(nloop,rand_mst,rand_list)
  
  avranmst = rand_mst(rand_list(NINT(REAL(nloop)/2.)))
  ranlow = rand_mst(rand_list(NINT(REAL(nloop)/6.)))
  ranup = rand_mst(rand_list(NINT(5.*REAL(nloop)/6.)))

  lambda_tilde(snapshoti) = avranmst/obj_mst
  l_low_tilde(snapshoti) = ranlow/obj_mst
  l_up_tilde(snapshoti) = ranup/obj_mst


END SUBROUTINE find_lambda_tilde
