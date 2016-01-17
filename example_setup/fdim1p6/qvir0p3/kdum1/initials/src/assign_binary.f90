   SUBROUTINE assign_binary(rm1,rm2,re,rsma,rv1,rv2,rr1,rr2,kdum)
!
   USE constants_module
   IMPLICIT NONE
   REAL :: rnum,ran2                        ! random number
   DOUBLE PRECISION :: f,sinphi,phi         ! angles and operations on angles
   DOUBLE PRECISION :: mu                   ! G*(sum of component masses)
   DOUBLE PRECISION :: magvsq               ! square of velocity magnitude
   DOUBLE PRECISION :: e                    ! eccentricity 
   DOUBLE PRECISION :: sma                  ! semi-major axis              
   DOUBLE PRECISION :: m1,m2                ! primary and secondary masses
   DOUBLE PRECISION :: magr                 ! separation magnitude
   DOUBLE PRECISION :: magv                 ! velocity magnitude
   DOUBLE PRECISION :: r1(1:3)              ! separation components star 1
   DOUBLE PRECISION :: r2(1:3)              ! separation components star 2
   DOUBLE PRECISION :: v1(1:3)              ! velocity components star 1
   DOUBLE PRECISION :: v2(1:3)              ! velocity components star 2
   DOUBLE PRECISION :: COM(1:3)             ! centre of mass of system
   DOUBLE PRECISION :: COV(1:3)             ! centre of velocity of system
   DOUBLE PRECISION :: angA,angB            ! matrix transformation angles
   REAL, INTENT(in) :: re,rsma,rm1,rm2      ! real values into subroutine
   REAL, INTENT(out) :: rr1(3),rr2(3),rv1(3),rv2(3) ! real values out of subroutine
   INTEGER, INTENT(in) :: kdum              ! random number seed     
   EXTERNAL ran2

!   write(6,*) rm1, rm2, rsma, re
!   pause

! get into double precision
   e=DBLE(re)
   sma=DBLE(rsma)
   m1=DBLE(rm1)
   m2=DBLE(rm2)
!
! Set up pi and 2*pi
   pi=4.*ATAN(1.)                ! Gets pi
   twopi=8.*ATAN(1.)             ! Gets 2*pi
!
!=============================================================================
!                          CONVERSION OF UNITS
!=============================================================================
! Convert from AU/pc, msun (N-body units) to MKS units  
!   sma=sma*au          ! convert from AU to metres
!  ------------------- (assumes input semi-major axis is in AU)
   sma=sma*pc          ! if input is in pc
   m1=m1*msun          ! convert from solar masses to kilograms 
   m2=m2*msun          ! convert from solar masses to kilograms 
!=============================================================================
!
!=============================================================================
!                  2-DIMENSIONAL COMPONENTS FOR STAR
!=============================================================================
! random number for angle f
   rnum=ran2(kdum)                ! gets random number from seed kdum
!
! Select random orientation in two-dimensional space 
   f=rnum*twopi                   ! chooses random angle between 0 and 2*pi  
!
! Get the magnitude of the separation vector between the binary components
   magr=(sma*(1. - (e**2.)))/(1. + (e*COS(f))) ! separation from m1 - m2 (metres)
!
! Get phi - the angle between separation vector and velocity vector
   sinphi=(((sma**2.)*(1. - (e**2.)))/(magr*((2.*sma) - magr)))**0.5
   phi=ASIN(sinphi) 
!
! Determine the magnitude of the orbital velocity vector
   mu=G*(m1 + m2)                   ! G x (mass component 1 + mass component 2)
   magvsq=mu*((2./magr) - (1./sma)) ! square of orbital velocity magnitude
   magv=SQRT(magvsq)                ! orbital velocity magnitude
!
! Set components of separation vector for star 1 to zero 
   r1(1:3)=0.d0
! Calculate the components of the separation vector for star 2 
   r2(1)=magr*COS(f)                   ! x-component of separation
   r2(2)=magr*SIN(f)                   ! y-component of separation
   r2(3)=0.d0                          ! z-component of separation
!
! Set components of velocity vector for star 1 to zero
   v1(1:3)=0.d0 
! Calculate the components of the velocity vector for star 2
   IF (f<pi) THEN
     v2(1)=-magv*COS(phi - f)        ! velocity component in x direction (f<180)
     v2(2)=magv*SIN(phi - f)         ! velocity component in y direction (f<180)
     v2(3)=0.d0                      ! velocity component in z direction (f<180)
   ELSE
     v2(1)=-magv*COS(twopi - f - phi)! velocity component in x direction (f>180)
     v2(2)=magv*SIN(twopi - f - phi) ! velocity component in y direction (f>180)
     v2(3)=0.d0                      ! velocity component in z direction (f>180)
   END IF
!
!============================================================================
!          CENTRE OF MASS (C.O.M.), CENTRE OF VELOCITY (C.O.V.)
!============================================================================
!
   COM(1:3)=(r1(1:3)*m1 + r2(1:3)*m2)/(m1 + m2)
   COV(1:3)=(v1(1:3)*m1 + v2(1:3)*m2)/(m1 + m2) 
!
!   write(6,*) COM,'       ',COV
!   write(6,*) v1,'        ',v2
   r1=r1 - COM
   r2=r2 - COM
   v1=v1 - COV  
   v2=v2 - COV
!   write(6,*) ''
!   write(6,*) r1,'        ',r2
!   write(6,*) v1,'        ',v2
!   write(6,*) ''
!
!=============================================================================
!                   RANDOMLY ORIENTATE SYSTEM AXES
!=============================================================================
! **** These angles are calculated here in case we use orientate to rotate eg the 
! entire cluster ******
!
! random number to get angle A
   rnum=ran2(kdum)          
!   angA = rnum*2.*pi
   angA=ACOS(rnum*2. - 1.)
!
! random number seed for angle B
   rnum=ran2(kdum) 
   angB = rnum*2.*pi
!
   CALL orientate(v1,v2,r1,r2,angA,angB)
!   write(6,*) 'stellar',r1,'        ',r2
!   write(6,*) 'stellar',v1,'        ',v2
!   write(6,*) ''
!   pause
!
!=============================================================================
!                          CONVERSION OF UNITS
!=============================================================================
! Convert masses back to N-body units (Solar Masses)
   m1=m1/msun  
   m2=m2/msun       
!
! Convert separation components and magnitude to AU/pc 
! In pc:
   magr=magr/pc
   r1=r1/pc
   r2=r2/pc
! In AU:
!   magr=magr/au
!   r2(1)=r2(1)/au
!   r2(2)=r2(2)/au
!   r2(3)=r2(3)/au
!
! Convert velocity components and magnitude to kms^-1
   magv=magv/1000.
   v1=v1/1000.
   v2=v2/1000.
!============================================================================
! put r1,r2,v1,v2 into real*4
   rr1=REAL(r1)
   rr2=REAL(r2)
   rv1=REAL(v1)
   rv2=REAL(v2)
!
   RETURN
   END SUBROUTINE assign_binary
