SUBROUTINE find_energy(snapshoti,ni)
  ! Find total kinetic, gravitational potential and total energy.
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! mi = mass of star i
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi
! ri, vi = position & velocity of star i
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ri, vi
  !logical, dimension(:), allocatable :: i_incluster !has this star escaped?
! ke = the total kinetic energy
! rij = the magnitude of the separation between two stars
! epoti = the gravitational potential energy between two stars
! epot = the total gravitational potential energy
! etot = the total energy
  DOUBLE PRECISION :: ekin,rij,epoti,epot,etot
  INTEGER :: i,j,k
  
  ALLOCATE(mi(1:ni))
  ALLOCATE(ri(1:3,1:ni))
  ALLOCATE(vi(1:3,1:ni))
  !allocate(i_incluster(1:ni))
  mi(1:ni)=m(1:ni,snapshoti)
  ri(1:3,1:ni)=r(1:3,1:ni,snapshoti)
  vi(1:3,1:ni)=v(1:3,1:ni,snapshoti)
  !i_incluster(1:ni)=incluster_FoV(4,1:ni,snapshoti) !in 3D cluster

! First convert to SI units
  mi=mi*msun
  ri=ri*pc
  vi=vi*1000
  ekin=0.
! Loop over all stars
  DO i=1,ni
     ! Find kinetic energy
     !if (.not. i_incluster(i)) cycle
     ekin = ekin + (0.5*mi(i)*(vi(1,i)**2+vi(2,i)**2+vi(3,i)**2))
  END DO
  
  epot=0.
  rij=0.
! Loop over all stars
  DO i=1,ni-1
     !if (.not. i_incluster(i)) cycle
     epoti=0.d0
     DO j=i+1,ni
        !if (.not. i_incluster(j)) cycle
! Find separation between two stars
        rij=(ri(1,i)-ri(1,j))**2 + (ri(2,i)-ri(2,j))**2 + & 
             &              (ri(3,i)-ri(3,j))**2
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
  !deallocate(i_incluster)
!!$  Print *, ekin, -epot, etot
END SUBROUTINE find_energy
