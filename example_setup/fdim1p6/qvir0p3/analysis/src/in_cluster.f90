!If certain stars are to be ignored then make a new array
!of stars in each snapshot where nstars = n_in_cluster
SUBROUTINE in_cluster(snapshoti,ni)
  
  USE constants_module
  USE parameters_module
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: snapshoti,ni
! n_escaped(1:3) gives number of escaped stars
! in xy, yz and xz in this snapshot.
! n_escaped(4) gives number of escaped stars in 3D in this snapshot.
  INTEGER, DIMENSION(1:4) :: n_escaped
! rmag = distance between star i and centre of mass
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rmag, rmagxy, rmagyz, rmagxz
  INTEGER :: i,j
! number of stars in 2D projections and 3D cluster:
  INTEGER :: n_xy, n_yz, n_xz, n_3D
! masses are just for data output; nothing to do with calculations
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mi

!Allocate memory for arrays
  ALLOCATE(rmag(1:ni))
  ALLOCATE(rmagxy(1:ni))
  ALLOCATE(rmagyz(1:ni))
  ALLOCATE(rmagxz(1:ni))
  ALLOCATE(mi(1:ni))
!Populate arrays
  rmagxy(1:ni)=r_com(snapshoti,1:ni,4)
  rmagyz(1:ni)=r_com(snapshoti,1:ni,5)
  rmagxz(1:ni)=r_com(snapshoti,1:ni,6)
  rmag(1:ni)=r_com(snapshoti,1:ni,7)
  mi(1:ni)=m(snapshoti,1:ni)

! initially all stars are in the cluster
  n_xy = ni
  n_yz = ni
  n_xz = ni
  n_3D = ni

  !do j=10, 13
  !   write(j,*) 'Snapshot',snapshoti
  !end do

! Loop over all stars
  DO i=1,ni

!===============
! Field of view
!===============

! if rmag is greater than FoV_lim, it has left the cluster
     IF (limittype=='FoV') THEN
        IF (rmagxy(i) > FoV_lim) THEN
           incluster(snapshoti,i,1) = .FALSE.
           n_xy = n_xy - 1
           !write(10,*) i, rmagxy(i)
        END IF
        IF (rmagyz(i) > FoV_lim) THEN
           incluster(snapshoti,i,2) = .FALSE.
           n_yz = n_yz - 1
           !write(11,*) i, rmagyz(i)
        END IF
        IF (rmagxz(i) > FoV_lim) THEN
           incluster(snapshoti,i,3) = .FALSE.
           n_xz = n_xz - 1
           !write(12,*) i, rmagxz(i)
        END IF
        IF (rmag(i) > FoV_lim) THEN
           incluster(snapshoti,i,4) = .FALSE.
           n_3D = n_3D - 1
           !write(13,*) i, rmag(i)
        END IF
     END IF


!=================
! halfmass radius
!=================

! if rmag is greater than rfac*r_halfmass, it has left the cluster

     IF (limittype=='rhalf') THEN
        IF (rmagxy(i) > rfac*r_halfmass(snapshoti,1)) THEN
           incluster(snapshoti,i,1) = .FALSE.
           n_xy = n_xy - 1
           !write(10,*) i, rmagxy(i)
        END IF
        IF (rmagyz(i) > rfac*r_halfmass(snapshoti,2)) THEN
           incluster(snapshoti,i,2) = .FALSE.
           n_yz = n_yz - 1
           !write(11,*) i, rmagyz(i)
        END IF
        IF (rmagxz(i) > rfac*r_halfmass(snapshoti,3)) THEN
           incluster(snapshoti,i,3) = .FALSE.
           n_xz = n_xz - 1
           !write(12,*) i, rmagxz(i)
        END IF
        IF (rmag(i) > rfac*r_halfmass(snapshoti,4)) THEN
           incluster(snapshoti,i,4) = .FALSE.
           n_3D = n_3D - 1
           !write(13,*) i, rmag(i)
        END IF
     END IF
  END DO


! ***********
! * Outputs *
! ***********
! Write out number of escaped stars if value
! is different from previous snapshot.
! Set an unphysical value for n_escaped in 1st snapshot
! so we are guaranteed to get a writeout:
  IF (snapshoti==1) THEN
     n_escaped=-1
  END IF

!xy:
! if new stars escaped, write snapshot number
  IF (ni - n_xy /= n_escaped(1)) THEN
     WRITE(10,'(2X,"Snapshot",I4,":",1X,I4,1X,"stars escaped")') &
          & snapshoti, ni-n_xy

     DO i=1,ni
! if star isn't in cluster, write distance magnitude & star mass
        IF (.NOT. incluster(snapshoti,i,1)) WRITE(10,80) i, rmagxy(i), mi(i)
     END DO

     IF (limittype=='rhalf') THEN
        WRITE(10,81) rfac, rfac * r_halfmass(snapshoti,1)
     END IF

! write number of stars escaped:
     n_escaped(1) = ni - n_xy
! For some reason need this in or the above 'if' is ignored...
     WRITE(10,'(30X,I4,1X,I4)') ni-n_xy, n_escaped(1)
        WRITE(10,*)""
  END IF

!yz:
! if new stars escaped, write snapshot number
  IF (ni - n_yz /= n_escaped(2)) THEN
     WRITE(11,'(2X,"Snapshot",I4,":",1X,I4,1X,"stars escaped")') &
          & snapshoti, ni-n_yz
     DO i=1,ni
! if star isn't in cluster, write distance magnitude & star mass
        IF (.NOT. incluster(snapshoti,i,2)) WRITE(11,80) i, rmagyz(i), mi(i)
     END DO

     IF (limittype=='rhalf') THEN
        WRITE(11,81) rfac, rfac * r_halfmass(snapshoti,2)
        WRITE(11,*)""
     END IF
! write number of stars escaped:
     n_escaped(2) = ni - n_yz
  END IF

!xz:
! if new stars escaped, write snapshot number
  IF (ni - n_xz /= n_escaped(3)) THEN
     WRITE(12,'(2X,"Snapshot",I4,":",1X,I4,1X,"stars escaped")') &
          & snapshoti, ni-n_xz
     DO i=1,ni
! if star isn't in cluster, write distance magnitude & star mass
        IF (.NOT. incluster(snapshoti,i,3)) WRITE(12,80) i, rmagxz(i), mi(i)
     END DO

     IF (limittype=='rhalf') THEN
        WRITE(12,81) rfac, rfac * r_halfmass(snapshoti,3)
     END IF
! write number of stars escaped:
     n_escaped(3) = ni - n_xz
        WRITE(12,*)""
  END IF

!3D:
! if new stars escaped, write snapshot number
  IF (ni - n_3D /= n_escaped(4)) THEN
     WRITE(13,'(2X,"Snapshot",I4,":",1X,I4,1X,"stars escaped")') &
          & snapshoti, ni-n_3D

     DO i=1,ni
! if star isn't in cluster, write distance magnitude & star mass
        IF (.NOT. incluster(snapshoti,i,4)) WRITE(13,80) i, rmag(i), mi(i)
     END DO

     IF (limittype=='rhalf') THEN
        WRITE(13,81) rfac, rfac * r_halfmass(snapshoti,4)
     END IF

! write number of stars escaped:x
     n_escaped(4) = ni - n_3D
! For some reason need this in or the above 'if' is ignored...
     WRITE(13,'(30X,I4,1X,I4)') ni-n_xy, n_escaped(4)
     WRITE(13,*)""
  END IF

80   FORMAT(1X,I4,2X,F8.3,2X,F7.3)
81   FORMAT(4X,I2,'*r_halfmass:',F8.3)


  DEALLOCATE(rmag)
  DEALLOCATE(rmagxy)
  DEALLOCATE(rmagyz)
  DEALLOCATE(rmagxz)
  DEALLOCATE(mi)

END SUBROUTINE in_cluster
