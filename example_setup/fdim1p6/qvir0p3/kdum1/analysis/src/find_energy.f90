
SUBROUTINE find_energy(snapshoti,ni)
  ! Find total kinetic, gravitational potential and total energy.
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi,ri,vi = mass, position and velocity of stars
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri, vi
! ke = the total kinetic energy
! rij = the magnitude of the separation between two stars
! epoti = the gravitational potential energy between two stars
! epot = the total gravitational potential energy
! etot = the total energy
  DOUBLE PRECISION :: ekin,rij,epoti,epot,etot
  INTEGER :: i,j,k
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:ni,1:3))
  ALLOCATE(vi(1:ni,1:3))
  mi(1:ni)=m(snapshoti,1:ni)
  ri(1:ni,1:3)=r(snapshoti,1:ni,1:3)
  vi(1:ni,1:3)=v(snapshoti,1:ni,1:3)
! First convert to SI units
  mi=mi*msun
  ri=ri*pc
  vi=vi*1000
  ekin=0.
! Loop over all stars
  DO i=1,ni
! Find kinetic energy
     ekin=ekin+(0.5*mi(i)*(vi(i,1)**2+vi(i,2)**2+vi(i,3)**2))
  END DO
  epot=0.
  rij=0.
! Loop over all stars
  DO i=1,ni-1
     epoti=0.d0
     DO j=i+1,ni
! Find separation between two stars
        rij=(ri(i,1)-ri(j,1))**2 + (ri(i,2)-ri(j,2))**2 + & 
             &              (ri(i,3)-ri(j,3))**2
! Find gravitational potential energy between two stars
        epoti=epoti + DBLE(mi(j)/SQRT(rij))
     END DO
! Find total gravitational potential energy
     epot=epot + G*DBLE(mi(i))*epoti
  END DO
  etot=ekin-epot
! add energies to array
  kinetic_energy(snapshoti)=ekin
  potential_energy(snapshoti)=-epot
  total_energy(snapshoti)=etot
! Deallocate arrays
  DEALLOCATE(mi)
  DEALLOCATE(ri)
  DEALLOCATE(vi)
!!$  Print *, ekin, -epot, etot
END SUBROUTINE find_energy

