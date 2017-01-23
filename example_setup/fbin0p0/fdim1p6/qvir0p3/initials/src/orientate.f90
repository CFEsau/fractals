   SUBROUTINE orientate(v1,v2,r1,r2,angA,angB) !,kdum also if called in main program 
!
   USE constants_module
   IMPLICIT NONE
!   INTEGER, intent(in) :: kdum        ! random number seed    - uncomment if called 
!   REAL :: rnum,ran2                  ! random number         - in main program
!   DOUBLE PRECISION :: pi,twopi       ! Constants: pi, 2*pi   ---------------------

   DOUBLE PRECISION :: angA,angB      ! matrix angles
   DOUBLE PRECISION :: v1p(1:3)       ! primary velocity components after 1st matrix 
   DOUBLE PRECISION :: v2p(1:3)       ! secondary velocity components after 1st matrix
   DOUBLE PRECISION :: v1pp(1:3)      ! primary velocity components after 2nd matrix
   DOUBLE PRECISION :: v2pp(1:3)      ! secondary velocity components after 2nd matrix
   DOUBLE PRECISION :: r1p(1:3)       ! primary separation components after 1st matrix
   DOUBLE PRECISION :: r2p(1:3)       ! secondary separation components after 1st matrix
   DOUBLE PRECISION :: r1pp(1:3)      ! primary separation components after 2nd matrix
   DOUBLE PRECISION :: r2pp(1:3)      ! secondary separation components after 2nd matrix
   DOUBLE PRECISION, INTENT(inout) :: v1(1:3) ! velocity vector components, primary star
   DOUBLE PRECISION, INTENT(inout) :: v2(1:3) ! velocity vector components, secondary star
   DOUBLE PRECISION, INTENT(inout) :: r1(1:3) ! separation vector components, primary star
   DOUBLE PRECISION, INTENT(inout) :: r2(1:3) ! separation vector components, secondary star   
!   EXTERNAL ran2  - uncomment if called in main program
!
! Set up pi and 2*pi if used as a main subroutine
!   pi=4.*ATAN(1.)               ! Gets pi
!   twopi=8.*ATAN(1.)            ! Gets 2*pi
!
!=============================================================================
!           CONVERSION OF UNITS - uncomment if called in main program
!=============================================================================
! Convert separation components from AU/pc (N-body) to metres (MKS)
! From pc:
!   r1=r1*pc
!   r2=r2*pc
! From AU:
!   r1=r1*au
!   r2=r2*au 
!
! Convert velocity components from kms^-1 to ms^-1
!   v1=v1*1000.
!   v2=v2*1000.
!
!==============================================================================
!                 Angles Used in Matrix Transformations
!==============================================================================
!
!  Un-comment random number in main program
!
! random number for angle A
!   rnum=ran2(kdum) ! rnum is random number
!   angA = rnum*2.*pi
!
! random number seed for angle B
!   rnum=ran2(kdum) ! rnum is random number
!   angB = rnum*2.*pi
!
! ===============================================================================
!                            MATRIX TRANSFORMATION
!================================================================================
! Matrix Rotation about the x-axis - we use this operator: 
!
!          (1    0       0  )
!  MatRx = (0  cos A  -sin A)
!          (0  sin A   cos A)
!
! This gives v' = MatRx(v):
   v1p(1)=v1(1)                              ! velocity' in x-direction \
   v1p(2)=v1(2)*COS(angA) - v1(3)*SIN(angA)  ! velocity' in y-direction  } Primary star
   v1p(3)=v1(2)*SIN(angA) + v1(3)*COS(angA)  ! velocity' in z-direction /
!
   v2p(1)=v2(1)                              ! velocity' in x-direction \
   v2p(2)=v2(2)*COS(angA) - v2(3)*SIN(angA)  ! velocity' in y-direction  } Secondary star
   v2p(3)=v2(2)*SIN(angA) + v2(3)*COS(angA)  ! velocity' in z-direction / 
!
! And r' = MatRx(r)
   r1p(1)=r1(1)                              ! separation' in x-direction \
   r1p(2)=r1(2)*COS(angA) - r1(3)*SIN(angA)  ! separation' in y-direction  } Primary star
   r1p(3)=r1(2)*SIN(angA) + r1(3)*COS(angA)  ! separation' in z-direction /
!
   r2p(1)=r2(1)                              ! separation' in x-direction \
   r2p(2)=r2(2)*COS(angA) - r2(3)*SIN(angA)  ! separation' in y-direction  } Secondary star
   r2p(3)=r2(2)*SIN(angA) + r2(3)*COS(angA)  ! separation' in z-direction /
!
! Matrix Rotation about the z-axis - we use this operator:
!
!          (cos B  -sin B  0)
!  MatRz = (sin B   cos B  0)
!          (  0       0    1)
!
! This gives v'' = MatRx(v'):
   v1pp(1)=v1p(1)*COS(angB) - v1p(2)*SIN(angB) ! velocity'' in x-direction \
   v1pp(2)=v1p(1)*SIN(angB) + v1p(2)*COS(angB) ! velocity'' in y-direction  } Primary star
   v1pp(3)=v1p(3)                              ! velocity'' in z-direction /
!
   v2pp(1)=v2p(1)*COS(angB) - v2p(2)*SIN(angB) ! velocity'' in x-direction \
   v2pp(2)=v2p(1)*SIN(angB) + v2p(2)*COS(angB) ! velocity'' in y-direction  } Secondary star
   v2pp(3)=v2p(3)                              ! velocity'' in z-direction /
!
! And r'' = MatRx(r')
   r1pp(1)=r1p(1)*COS(angB) - r1p(2)*SIN(angB) ! separation'' in x-direction \
   r1pp(2)=r1p(1)*SIN(angB) + r1p(2)*COS(angB) ! separation'' in y-direction  } Primary star
   r1pp(3)=r1p(3)                              ! separation'' in z-direction /
!
   r2pp(1)=r2p(1)*COS(angB) - r2p(2)*SIN(angB) ! separation'' in x-direction \
   r2pp(2)=r2p(1)*SIN(angB) + r2p(2)*COS(angB) ! separation'' in y-direction  } Secondary star
   r2pp(3)=r2p(3)                              ! separation'' in z-direction /
!
!============================================================================================
! Re-assign values (v = v'' and r = r''):
   v1(1:3)=v1pp(1:3)
   v2(1:3)=v2pp(1:3)
!
   r1(1:3)=r1pp(1:3)
   r2(1:3)=r2pp(1:3)
!
!=============================================================================
!              CONVERSION OF UNITS  - uncomment if called in main program
!=============================================================================
! Convert separation components to AU/pc
! In pc:
!   r1=r1/pc
!   r2=r2/pc
! In AU:
!   r1=r1/au
!   r2=r2/au
!
! Convert velocity components to kms^-1
!   v1=v1/1000.
!   v2=v2/1000.
!===========================================================================
   RETURN
   END SUBROUTINE orientate

