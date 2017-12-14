!***************************************************!
! 
! particletree.f90
!
! This contains the subroutines for analysing the
! particle information read in from read_sl_out.
!
!***************************************************!

SUBROUTINE particletree(snapi)
!
! This subroutine reconstructs the hierarchical tree. It searches the starlab
! hierarchy so that all the particles 'know' which other particles they are related
! to. This allows the user to obtain information on a single star in a hierarchical
! system and then relate that to the information describing the system itself.
!
! 1 T     e.g. in this tree there are 5 particles (npart=5) from j=1 to j=5. There
! 2   B     is 'system' info for the triple and the binary within this triple,
! 3     S   and also 'star' info for the single star and the stars that form the
! 4     S   binary. All of this info is preserved using the following code...
! 5   S
  use sl_input_module
  use parameters_module
  implicit none
  integer, intent(in) :: snapi        ! This snapshot number
  integer :: nsing                    ! number of single stars
  integer :: parent                   ! Particles particle ID
  integer :: ncounter                 ! counts number of stars in each multiple
  integer, dimension(npart) :: childof,sibling
  integer :: nfalse
  integer :: i,k                      ! Iterators
  logical, dimension(npart) :: done
  !    Total mass, velocity, and position of particle (sum of particle + parent):
  double precision, dimension(npart) :: totm
  double precision, dimension(3,npart) :: totv,totr


  done=.FALSE.             ! Set a variable 'done' to false
  parent=0                 ! Zero parent/sibling/child integers
  childof=0
  
  do
     ncounter=0             ! Set ncouncter to be zero
     nfalse=0               ! Set number of 'false' occurrences to be zero
     do i=1,npart
        if (done(i)) cycle   ! if value of i has been assigned a parent/child, cycle
        if (ncounter==0) then  ! if ncounter is zero, go into 'if' block
           if (multN(i)>1) ncounter=multN(i) ! If multiple, set ncounter to be the number
                                           ! of constituents in the multiple particle.
           parent=i                        ! Set parent to be i
           done(i)=.true.                  ! Set flag 'done' to be true
        else
           if (multN(i)==1) ncounter=ncounter-1 ! if single star, -1 from total in system
           childof(i)=parent         ! set this star's parent to be the last multiple system
        end if
        if (.not. done(i)) nfalse=nfalse+1
     end do
     if (nfalse==0) exit          ! exit loop when all parents/children have been assigned
  end do
  
! This loop assigns siblings to particles
  do i=1,npart
     if (childof(i)==0) cycle      ! if it's not a child then there are no siblings.
     if (sibling(i) .gt. 0) cycle  ! sibling already been assigned, move to next particle
     do k=1,npart
        if (i==k) cycle            ! sibling cannot be sibling to itself, so cycle
        if (childof(k)==childof(i)) then
           sibling(k)=i
           sibling(i)=k
           exit                    ! exit when siblings are assigned
        end if
     end do
  end do

! Display parents, siblings and children in a table
!  write(6,*) 'particle number |',' number of sub-particles |',' parent |',' sibling'
!  do i=1,npart
!     write(6,*) '     ',i,'             ',multN(i),'        ',childof(i),'      ,'sibling(i)
!  end do

  
  totr=0
  totv=0
  totm=0
  do i=1,npart     
     if (childof(i)==0) then    ! If particle has no parents (top-level particle)
        totr(1:3,i)=rmax(1:3,i,snapi)    ! then its r, v and m are right already
        totv(1:3,i)=vmax(1:3,i,snapi)
        totm(i)=mmax(i,snapi)
     else                    ! It it has a parent then add those values
        parent=childof(i)         ! ID of i's parent
        do
           ! sum particle i's values  with parent particle parent's values:
           totr(1:3,i)=rmax(1:3,i,snapi) + rmax(1:3,parent,snapi)
           totv(1:3,i)=vmax(1:3,i,snapi) + vmax(1:3,parent,snapi)
           totm(i)=mmax(i,snapi)
           parent=childof(parent)       ! ID of the parent's parent
           if (parent==0) exit     ! Exit if top level
        end do
     end if
  end do
  
  ! Make all objects from this point onwards single stars:
  nsing=0
  do i=1,npart
     if (multN(i)==1) then
        nsing=nsing+1
        idstar(nsing,snapi)=idmax(i,snapi)
        rstar(1:3,nsing,snapi)=totr(1:3,i)   ! position of single star
        vstar(1:3,nsing,snapi)=totv(1:3,i)   ! velocity of single star
        mstar(nsing,snapi)=totm(i)           ! mass of single star
        tstar(nsing,snapi)=tmax(i,snapi)
     end if
  end do

  if (nsing/=numstars) then
     write(6,*)'nsing ',nsing,' not equal to numstars ',numstars; stop
  end if

end subroutine particletree
