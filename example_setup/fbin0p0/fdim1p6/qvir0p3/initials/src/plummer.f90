!********************************************************!
! Updtaed 01/12/14 DG
! 
! Initial Conditions Generator for KIRA N-Body Program
!
! Subroutine: Plummer Sphere Creator
! Nicked from Richard Parker (annotations added)
! 
! Algorithm taken from Aarseth, Henon and Wielen, 1974 
! (A&A 37, 183)
!
!********************************************************!
!
! Inputs:
! - n = Number of systems
! - seed = random number seed
!
! Outputs:
! - r = x,y,x positions for all systems
! - v = vx,vy,vz velocities for all systems
!
SUBROUTINE plummer(n,r,v,seed)
  IMPLICIT NONE
  ! Declare inputs
  INTEGER, INTENT(in) :: n, seed
  ! Declare output arrays
  REAL, INTENT(out) :: r(1:3,1:n),v(1:3,1:n)
  ! x(i) = Random numbers
  ! ran2 = random number generator
  ! twopi = constant 2*pi
  ! ri = radial distance
  ! vi = velocity magnitude
  ! gq = g(q) from Aarseth
  ! i = iterator
  REAL :: x(7)
  INTEGER :: i
  REAL :: ran2, twopi, ri, vi, gq
  ! Random number generator RAN2
  EXTERNAL ran2
  ! Define constants
  twopi = 8.*ATAN(1.)
!
!********************************************************!
  ! Loop over all stars
  DO i=1,n
     ! Find r:
     ! Get a normalised random number
48   x(1) = RAN2(seed)
     ! Reject number if very small
     IF(x(1)<1.0e-10) GO TO 48
     ! Calculate radius corresponding to number
     ri = ((x(1)**(-2./3.))-1.)**(-0.5)
     ! If ri is very distant, i.e. greater than 14.6154 (for some reason)
     IF(ri>14.6154) GO TO 48
     ! Calculate x,y,z:
     ! Get 2 normalised random numbers
     x(2) = RAN2(seed)
     x(3) = RAN2(seed)
     ! Get z, then get x and y from z
     r(3,i) = (1.-2.*x(2))*ri
     r(1,i) = SQRT(ri**2-r(3,i)**2)*COS(twopi*x(3))
     r(2,i) = SQRT(ri**2-r(3,i)**2)*SIN(twopi*x(3))
!
!********************************************************!
     ! Find v:
     ! - von Neumann's rejection technique:
67   x(4) = RAN2(seed)
     x(5) = RAN2(seed)
     ! find g(x(4))
     gq = (x(4)**2)*((1.-x(4)**2)**3.5)
     ! If 0.1*x(5)<g(q) then q=x(4); else reject the two numbers 
     IF(0.1*x(5)>gq) GO TO 67
     ! v = q*v(esc)
     vi = x(4)*SQRT(2.)/((1+ri**2)**0.25)
     ! Then find vx,vy,vz in the same way as x,y,z
     ! Get 2 normalised random numbers
     x(6) = RAN2(seed)
     x(7) = RAN2(seed)
     ! Get vz, then get vx and vy from vz
     v(3,i) = (1.-2.*x(6))*vi
     v(1,i) = SQRT(vi**2-v(3,i)**2)*COS(twopi*x(7))
     v(2,i) = SQRT(vi**2-v(3,i)**2)*SIN(twopi*x(7))
!
!********************************************************!
  END DO
  RETURN
END SUBROUTINE plummer


