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

SUBROUTINE find_lambda_all(snapshoti,ni)
  USE parameters_module

  IMPLICIT NONE

!-----------
! star data
!-----------
  INTEGER, INTENT(in) :: snapshoti, ni !Snapshot number and star number
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
! Length of object tree (e.g. most massive stars) for different lambda:
  DOUBLE PRECISION :: obj_mst, obj_mst_bar, obj_mst_til, obj_mst_star
  DOUBLE PRECISION :: obj_mst_gam
! obj_r = positions of obj stars
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: obj_r
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x,y,z !split of obj_r
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: node !node connections for mst

!total lengths of random MSTs for different lambda:
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rand_mst, rand_mst_bar
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rand_mst_til, rand_mst_star
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rand_mst_gam

  INTEGER, DIMENSION(:), ALLOCATABLE :: rand_list !IDs of rand_mst list

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

!Find distance of each star from centre of mass:
  DO i = 1, ni

     ri_x = ri_com(snapshoti,i,1)
     ri_y = ri_com(snapshoti,i,2)
     ri_z = ri_com(snapshoti,i,3)

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

     i_incluster(i) = incluster(snapshoti,i)
  END DO


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

     IF (.NOT. i_incluster(k)) THEN
        !WRITE(51,'(1X,I4,1X,I5,1X,I5,1X,F7.3,1X,F8.3,1X,A1)') &
             !& snapshoti, i, k, mi(k), rmag(k), i_incluster(k)

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
        WRITE(50,99) snapshoti, ti, obj_mass(2,1:nmst)
        EXIT
     END IF
  END DO
99   FORMAT(1X,I4,2X,E9.3,*(2X,F8.3))

  obj_mass(1,1:nmst)=obj_mass(2,1:nmst)

!======================================================================
!Find the MST length for the OBJECT stars
!======================================================================
!
  length = 0.
  obj_mst = 0.
  obj_mst_bar = 0.
  obj_mst_til = 0.
  obj_mst_star = 0.
  obj_mst_gam = 0.

!Split obj_r across x,y,z arrays. Just helps readability in mst subroutine.
  x(1:nmst)=obj_r(1:nmst,1)
  y(1:nmst)=obj_r(1:nmst,2)
  z(1:nmst)=obj_r(1:nmst,3)

!set unit for output to file
  unit1=4
  CALL mst(nmst,x,y,z,node,length,snapshoti)

!Lambda MST:
  DO i = 1,nmst-1
     obj_mst = obj_mst + length(i)  !Add the edges of the mst to find 
  END DO                            !the total length

!Lambda bar MST:
  DO i = 1,nmst-1
     obj_mst_bar = obj_mst_bar + length(i)  !Add the edges then find mean
  END DO
  obj_mst_bar = obj_mst_bar/(nmst-1)


!Lambda tilde MST:
! Find the median vertex length
  DO i = 1, nmst-1
     length_list(i) = i
  END DO
  CALL heapsort(nmst-1, length, length_list)
  obj_mst_til = length(length_list( NINT(REAL(nmst-1)/2.) ) )


!Lambda star MST:
! Find the median vertex length
  DO i = 1, nmst-1
     length_list(i) = i
  END DO
  CALL heapsort(nmst-1, length, length_list)
! Find the median length in the tree, & find length of a tree made from these
  obj_mst_star = length(length_list( NINT(REAL(nmst-1)/2.) ) )
! Then add on the actual length of the tree
  DO i = 1, nmst-1
     obj_mst_star = obj_mst_star + length(i)  !Add the edges of the mst to
  END DO                      !find the total length


!Gamma MST:
  DO i = 1,nmst-1
     obj_mst_gam = obj_mst_gam + LOG(length(i))  !Add the edges of the mst to
  END DO                                 !find the total length
!Calculate the geometric mean
  obj_mst_gam = EXP( (1./REAL(nmst-1)) * obj_mst_gam)


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
  ALLOCATE(rand_mst_bar(1:nloop))
  ALLOCATE(rand_mst_til(1:nloop))
  ALLOCATE(rand_mst_star(1:nloop))
  ALLOCATE(rand_mst_gam(1:nloop))

  rand_mst = 0. !set lengths of random MSTs to 0.
  rand_mst_bar = 0.
  rand_mst_til = 0.
  rand_mst_star = 0.
  rand_mst_gam = 0.

  ALLOCATE(rand_list(1:nloop))


!set unit for output to file
  unit1=5

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

     CALL mst(nmst,x,y,z,node,length)

!Lambda MST:
     DO i = 1,nmst-1
        rand_mst(j) = rand_mst(j) + length(i)    !Add the edges of the mst
     END DO                                  !to find the total length


!Lambda bar MST:
     DO i = 1,nmst-1
        rand_mst_bar(j) = rand_mst_bar(j) + length(i) !Add the edges
     END DO   
     rand_mst_bar(j) = rand_mst_bar(j)/(nmst-1) !& find mean


!Lambda tilde MST:
! Find the median veretx length
     DO i = 1, nmst-1
        length_list(i) = i
     END DO
     CALL heapsort(nmst-1, length, length_list)
     rand_mst_til(j) = length(length_list( NINT(REAL(nmst-1)/2.) ) )


!Lambda star MST:
! Find the median veretx length
     DO i = 1, nmst-1
        length_list(i) = i
     END DO
     CALL heapsort(nmst-1, length, length_list)
! Find the median length in the tree, & find length of a tree made from these
     rand_mst_star(j) = (nmst-1) * length(length_list( NINT(REAL(nmst-1)/2.) ) )
! Then add on the actual length of the tree
     DO i = 1, nmst-1
        rand_mst_star(j) = rand_mst_star(j) + length(i)
     END DO


!Gamma MST:
     DO i = 1,nmst-1
        rand_mst_gam(j) = rand_mst_gam(j) + LOG(length(i))    !Add the edges of the mst
     END DO                                  !to find the total length
     rand_mst_gam(j) = EXP( (1./REAL(nmst-1)) * rand_mst_gam(j) )

  END DO


!======================================================================
!Finally - find the average value of the random MSTs and the significance
!======================================================================
!
!The average is found using the median value - this means that small
!values don't skew the value
!The signifincance is found using the 1/6 and 5/6 boundaries


!Calculate lambda:
  CALL calc_lambda(obj_mst,rand_mst,lambda(snapshoti), &
       & l_low(snapshoti),l_up(snapshoti))


!Calculate lambda bar:
  CALL calc_lambda(obj_mst_bar,rand_mst_bar,lambda_bar(snapshoti), &
       & l_low_bar(snapshoti),l_up_bar(snapshoti))


!Calculate lambda tilde:
  CALL calc_lambda(obj_mst_til,rand_mst_til,lambda_til(snapshoti), &
       & l_low_til(snapshoti),l_up_til(snapshoti))


!Calculate lambda star:
  CALL calc_lambda(obj_mst_star,rand_mst_star,lambda_star(snapshoti), &
       & l_low_star(snapshoti),l_up_star(snapshoti))


!Calculate gamma:
  CALL calc_lambda(obj_mst_gam,rand_mst_gam,gamm(snapshoti), &
       & g_low(snapshoti),g_up(snapshoti))


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
  DEALLOCATE(rand_mst)
  DEALLOCATE(rand_mst_bar)
  DEALLOCATE(rand_mst_til)
  DEALLOCATE(rand_mst_star)
  DEALLOCATE(rand_mst_gam)
  DEALLOCATE(rand_list)

END SUBROUTINE find_lambda_all
