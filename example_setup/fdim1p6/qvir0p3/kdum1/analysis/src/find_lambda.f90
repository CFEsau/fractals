!***Modifications***
!-Ignore stars a certain distance from the com

!This subroutine finds the value of mass segregation and its significance
!for the set of OBJECT STARS that you provide it.
!It works by calculating the MST for the OBJECT STARS, and then returning
!a value produced from sampling MSTs from randomly selected stars.
!This subroutine can be treated more-or-less like a black box.
!
!This algorithm is presented in Allison et al, 2009, MNRAS, 395, 1449
!

subroutine find_lambda(snapshoti,ni)
  use parameters_module
  implicit none
  integer, intent(in) :: snapshoti, ni !Snapshot number and star number
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi !masses of all stars
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri !positions of all stars
  integer, dimension(:), allocatable :: mlist !ID numbers for heapsort

  integer :: ndim !2D or 3D analysis
  
  integer :: nmst !Number of stars in each MST
  double precision, dimension(:), allocatable :: length !length of each MST edge
  double precision :: obj_mst !Length of object tree (e.g. most massive stars)
  double precision, dimension(:,:), allocatable :: obj_r !positions of obj stars
  double precision, DIMENSION(:), ALLOCATABLE :: x,y,z !split of obj_r
  integer, dimension(:,:), allocatable :: node !node connections for mst

  integer :: nloop !number of random MSTs in loop
  double precision, DIMENSION(:), ALLOCATABLE :: rand_mst !total lengths of random MSTs
  INTEGER, DIMENSION(:), ALLOCATABLE :: rand_list !IDs of rand_mst list
  real :: avranmst !Average length of random trees
  real ::ranup, ranlow !Upper & lower boundaries, 1 sigma

  LOGICAL, DIMENSION(:), ALLOCATABLE :: done !record which stars have been randomly selected
  REAL :: rand
  
  integer :: i,j,k !generic counters

!Allocate memory for arrays
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(mlist(1:ni))
  ALLOCATE(done(1:ni))
!And fill arrays 
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)

!Set number of stars in MST:
  nmst = 10
!and allocate memory for associated arrays:
  allocate(length(1:nmst-1))
  allocate(obj_r(1:nmst,1:3))
  allocate(x(1:nmst))
  allocate(y(1:nmst))
  allocate(z(1:nmst))
  allocate(node(1:nmst-1,1:2))
  
!2D or 3D analysis?
  ndim = 2

!======================================================================
!Sort stars in order of mass & make selection
!======================================================================
!

  obj_r=0.
  mlist=0
  
  do j=1,ni
     mlist(j)=j
  end do
  
  call heapsort(ni,mi,mlist)
  
!select nmst most massive stars:
  !(heapsort orders from smallest to largest, so invert)
  do i = 1,nmst
     k = mlist(ni+1-i) !heapsort orders from small to large - invert
     obj_r(i,1) = ri(k,1) !x position
     obj_r(i,2) = ri(k,2) !y position
!For 2D set z to 0
     IF (ndim == 2) then
        obj_r(i,3) = 0. !z position
     else if (ndim == 3) then
        obj_r(i,3) = ri(k,3) !z position 
     else
        stop 'ndim must be 2 or 3'
     end if
  end do
  
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
   !IF (z(1)==0) WRITE(6,*) '====***2D OBJECT MST***===='
   CALL mst(nmst,x,y,z,node,length)
   
   DO i = 1,nmst-1
      obj_mst = obj_mst + length(i)  !Add the edges of the mst to find 
   END DO                            !the total length

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
5        CALL RANDOM_NUMBER(rand) !random number between 0. and 1.
         k = NINT(rand*ni)
         IF (k == 0) GOTO 5    !There is no star with id=0
         IF (done(k)) GOTO 5    !Don't chose the same star twice
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

      DO i = 1,nmst - 1
         rand_mst(j) = rand_mst(j) + length(i)    !Add the edges of the mst
      END DO                                  !to find the total length
   END DO

   !IF (ndim == 2) WRITE(6,*) '====***2D RANDOM MST***===='

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
   
   avranmst = rand_mst(rand_list(NINT(REAL(nloop)/2.))) !Mean MST.
   ranlow = rand_mst(rand_list(NINT(REAL(nloop)/6.)))
   ranup = rand_mst(rand_list(NINT(5.*REAL(nloop)/6.)))

   lambda(snapshoti) = avranmst/obj_mst
   l_low(snapshoti) = ranlow/obj_mst
   l_up(snapshoti) = ranup/obj_mst

  ! if ((snapshoti==1).or. (snapshoti==900))then
  !    print*,"avranmst",avranmst
  !    print*,"ranlow",ranlow
  !    print*,"ranup",ranup
  !    print*,"obj_mst",obj_mst
  ! end if
!Deallocate arrays:
   DEALLOCATE(mi)
   DEALLOCATE(ri)
   DEALLOCATE(mlist)
   DEALLOCATE(done)
   deallocate(length)
   deallocate(obj_r)
   deallocate(x)
   deallocate(y)
   deallocate(z)
   deallocate(node)

 END SUBROUTINE find_lambda
