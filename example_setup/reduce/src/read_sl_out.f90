!******************************************************************************!
!
! 09/05/17 Claire Esau
! 
! read_sl_out.f90
!
!******************************************************************************!

SUBROUTINE read_sl_out(slName)

! ===== Declarations ==========================================
  USE sl_input_module    ! variables needed for file input
  USE parameters_module  ! m, t, r, v, etc
  IMPLICIT NONE
  integer :: snapi                     ! snapshot counter
  INTEGER :: j                         ! iterator
  INTEGER :: snapreq                   ! required snapshot number
  CHARACTER*150, INTENT(in) :: slName  ! name of sl file (input file)
  INTEGER :: inparticle  ! >0 if we're in an instance of '(Particle' in sl file
! ===== End of declarations ===================================
!
  call count_slfile(slName) ! find numstars & nsnaps
  
  nmax=numstars*3     ! guess number of 'particles' in starlab tree
  
! Allocate memory for arrays
!
! Stellar information
  ALLOCATE(tmax(1:nmax,1:nsnaps))
  ALLOCATE(mmax(1:nmax,1:nsnaps))
  ALLOCATE(rmax(1:3,1:nmax,1:nsnaps))
  ALLOCATE(vmax(1:3,1:nmax,1:nsnaps))
  ALLOCATE(idmax(1:nmax,1:nsnaps))    ! star ids for each snapshot
  ALLOCATE(multN(1:nmax))             ! number of stars in each 'particle'
  tmax=0.
  mmax=0.
  rmax=0.
  vmax=0.
  idmax = 0  ! zero particle id array for all particles and snapshots

  
  ALLOCATE(nMult(1:nsnaps))        ! number of multiples in each snapshot
  ALLOCATE(nstars(1:nsnaps))       ! number of stars in each snapshot
  nMult=0
  nstars=0
  snapi=0     ! Snapshot counter
!
!
! Data arrays for single stars:
  allocate(tstar(1:numstars,1:nsnaps))
  allocate(mstar(1:numstars,1:nsnaps))
  allocate(rstar(1:3,1:numstars,1:nsnaps))
  allocate(vstar(1:3,1:numstars,1:nsnaps))
  ALLOCATE(idstar(1:numstars,1:nsnaps))
  tstar=0.
  mstar=0.
  rstar=0.
  vstar=0.
!
! Open the runfile
  OPEN(2,file=TRIM(slName),status='old')
! Open an output file for snapshot top info
!!$  OPEN(3,file='out_topinfo',status='new')
!
!**************************************!
! 
! READING SL FILE
!
!**************************************!
  inparticle=0
  DO
     READ(2,*,END=1) descriptor              ! Read next line until end of file
     IF (descriptor=='(Particle') THEN       ! If line is an opening 'Particle'
        inparticle=inparticle+1              ! Count '(Particle' instances
        READ(2,*) descriptor                 ! reads next line
        IF (descriptor=='name') THEN            ! If next line is 'name'
           BACKSPACE 2                          ! Rewinds file one line
           READ(2,*) descriptor,equals,partname ! Reads line again
           IF (partname=='root') THEN
              snapi=snapi + 1               ! new snapshot
              j = 0                         ! zero particle counter j & particle
              multN = 0                     ! family counter (new family tree)
!
              CALL read_topinfo    ! read info from top of snapshot output
! This subroutine reads the 'top' or 'system' information from the Starlab
! output. Currently nothing is required from this but it is called in case
! the user requires it and so the later output is read correctly.
              ! WRITE(6,*) "top information extracted."
!
!              
!***************************************************!
! READING MAIN SNAPSHOT INFO
!***************************************************!
           ELSE ! If name/='root' & is after '(Particle' then it's a multiple
                ! A single has an 'i =' line before the 'name =' line.
!
              j=j+1                                ! increment particle counter
              READ(2,*) descriptor,equals,multN(j) ! 'N' stars in this multiple
              ! WRITE(6,*) 'multN = ',multN(j)!
              idmax(j,snapi)=0
              
!             === Extract multiple particle information ===
              CALL read_multinfo(j,snapi)
!             This subroutine extracts infrmation from 'multiple particles'
!             and places this info into the respective arrays mmax(j,snapi),
!             rmax(1:3,j,snapi), vmax(1:3,j,snapi) & tmax(j,snapi)
!             WRITE(6,*) "Multiple particle information extraced"
!             WRITE(6,*) ""
           END IF                               ! end of root/multiple 'if'
!
        ELSE IF(descriptor=='i') THEN           ! descriptor 'i=', not 'name='
           nstars(snapi)=nstars(snapi)+1        ! add 1 to number of stars
           BACKSPACE 2
           j=j+1                                ! increment particle counter
           READ(2,*) descriptor, equals, idmax(j,snapi) ! read line for star ID
           CALL read_singinfo(j,snapi)
           multN(j)=1
        END IF        ! End of descriptor = name/i. Still in (Particle 'if'.
        

     ELSE IF (descriptor==')Particle') then
        inparticle=inparticle-1       ! We have left an instance of '(Particle'
     end if
!
!==================================
!    Make particle family tree
!==================================
!     
     if (inparticle==0) then  ! If we are in one instance of '(Particle'
        !if (partname .ne. 'root') then
! After reading info then we are at the end of the snapshot. Do hierarchy.
           npart=j                    ! Set number of particles from counter j
           call particletree(snapi)   ! Find which particles are related
        !end if
     end if  
     
  END DO ! End of reading file
1 PRINT *, "End of sl file"
  !PRINT *, "Snapshots: ", snapi
  CLOSE(2)
!!$  CLOSE(3)

END SUBROUTINE
