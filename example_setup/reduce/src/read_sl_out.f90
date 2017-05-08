!******************************************************************************!
!
! 08/12/14 Dan Griffiths
! 
! read_sl_out.f90
!
!******************************************************************************!

SUBROUTINE read_sl_out(slName)

! Declare modules used
! input_module has most of the variables needed for file input, 
! plus descriptions of those variables.
  USE sl_input_module
  IMPLICIT NONE
!
! Declare variables
! Iterators
  INTEGER :: i,j,k,l
! slName = the inputted name of the sl file
  CHARACTER*150, INTENT(in) :: slName
!
!******************************************************************************!
!
! Allocate memory for arrays
! Numbers of stars/multiples
  ALLOCATE(nMult(1:3000))
  ALLOCATE(nstars(1:3000))
! Multiple information
  ALLOCATE(mult_nstars(1:3000,1:3000))
  ALLOCATE(mult_ids(1:3,1:3000,1:3000))
  ALLOCATE(mult_t(1:3000,1:3000))
  ALLOCATE(mult_m(1:3000,1:3000))
  ALLOCATE(mult_r(1:3,1:3000,1:3000))
  ALLOCATE(mult_v(1:3,1:3000,1:3000))
! Initialise variables
  snapnum=0 ! Number of snapshots
  i=0 ! Iterator
  top_energy=0 ! Total energy according to top info
!
! Stellar information
  ALLOCATE(ids(1:3000,1:3000))
  ALLOCATE(star_t(1:3000,1:3000))
  ALLOCATE(star_m(1:3000,1:3000))
  ALLOCATE(star_r(1:3,1:3000,1:3000))
  ALLOCATE(star_v(1:3,1:3000,1:3000))
! Initialise arrays
  nMult=0
  nstars=0
l=0
!
! Open the runfile
  OPEN(2,file=TRIM(slName),status='old')
! Open an output file for snapshot top info
!!$  OPEN(3,file='out_topinfo',status='new')
!
!******************************************************************************!
! 
! READING SL FILE
!
!******************************************************************************!
  DO
     READ(2,*,END=1) descriptor                 ! Reads line
     IF (descriptor=='(Particle') THEN          ! If line is an opening 'Particle'
        READ(2,*) descriptor                    ! reads next line
        IF (descriptor=='name') THEN            ! If next line is 'name'
           BACKSPACE 2                          ! Rewinds file one line
           READ(2,*) descriptor,equals,partname ! Reads next line again
           IF (partname=='root') THEN
              snapnum=snapnum + 1 ! new snapshot if 'name = root'
!
! It's a new snapshot, so it will start by reading the top info
!              
!******************************************************************************!
! 
! READING TOPINFO
!
! part of this is just a sanity check, to make sure KIRA doesn't calculate a change
! in total energy or something like that. So the output file out_topinfo is for
! diagnostics.
!
!******************************************************************************!
!
              CALL read_topinfo
              l=l+1
              !Print *,l
!
!******************************************************************************!
!
! END OF READING TOPINFO
!
! READING MAIN SNAPSHOT INFO
!
!******************************************************************************!
!
! Now to take care of the systems themselves
! If description does equal 'name', but that name isn't 'root'...
! The difference between a single star and a multiple is that a single star has an
!'i =' line before its 'name =' line. Therefore the  single stars can be taken care
! of seperately underneath. Here we take care of multiples.
! 
           ELSE IF(partname/='root') THEN ! Multiple
!
              BACKSPACE 2 ! Rewinds file one line
              nMult(snapnum) = nMult(snapnum)+1 ! The number of multiples in the snapshot increases by 1
              CALL read_multinfo
!
           END IF
!
!******************************************************************************!
!
! SINGLE STARS
!
!******************************************************************************!
!
! If descriptor does not equal 'name', but instead equals 'i'...
        ELSE IF(descriptor=='i') THEN
           nstars(snapnum)=nstars(snapnum)+1
           IF(l==1005) Print *, 'hello'
           BACKSPACE 2
           IF(l==1005) Print *, 'hello'
           CALL read_singinfo
        END IF
     END IF
  END DO ! End of file
1 PRINT *, "End of sl file"
  PRINT *, "Snapshots: ", snapnum
  CLOSE(2)
!!$  CLOSE(3)

END SUBROUTINE
