!*************************************************!
! 
! read_sl_sub.f90
!
! This contains the subroutines used for read_sl
!
!*************************************************!

subroutine count_slfile(slName)
  use sl_input_module
  
  CHARACTER*150, INTENT(in) :: slName  ! name of sl file (input file)
  
! Find the number of stars in this simulation. Take 
! number of stars in 1st snapshot as none are added later
  open(2,file=TRIM(slName),status='old') ! Open file
  do
     read(2,*) descriptor
     if (descriptor=='N') then   ! read N= at top of file for numstars
        backspace 2
        read(2,*) descriptor,equals,numstars
        exit
     end if
  end do
  close(2) ! Close file
  PRINT *, "Number of stars: ", numstars

! Count the number of snapshots
  nsnaps=0     ! Set the number of snapshots to zero
  open(2,file=TRIM(slName),status='old') ! Open file
  do
     READ(2,*,end=20) descriptor    ! Read each line until end of file (20)
     if (descriptor=='name') then   ! If line is 'name'
        backspace 2                 ! rewind one line
        read(2,*) descriptor,equals,partname      ! read next line
        if (partname=='root') nsnaps=nsnaps+1 ! If 'root',= a new snapshot
     end if
  END DO
20  close(2)    ! File has ended - close.
  PRINT *, "Number of snapshots: ", nsnaps
  
end subroutine count_slfile




SUBROUTINE read_topinfo
! This subroutine reads the 'top' or 'system' info from the Starlab
! output. Currently nothing is required from here but it is called
! in case the user requires it and so the output is read correctly.
!
  USE sl_input_module
  USE parameters_module
  IMPLICIT NONE
!
! Start reading file
  READ(2,*) descriptor                 ! Reads line
  IF (descriptor=='N') THEN            ! If line starts with 'N ='
     BACKSPACE 2                       ! Rewinds file one line
     READ(2,*) descriptor,equals,Ntop  ! Reads in number of stars
  END IF
!
  DO WHILE(descriptor/='(Particle')   ! Topinfo ends when next
                                      ! '(Particle' line is read.
     READ(2,*) descriptor             ! Read next line

! Example loop to extract information from 'Log' at the top of the output.
! Currently nothing is required from here but the below commented out loop
! allows the user to obtain information if they wish.

!     IF (descriptor=='(Log') THEN
!        DO WHILE(descriptor/=')Log')
!           READ(2,*) descriptor
!           ! === Example of pulling data from loop ===
!           IF (descriptor=='initial_mass') THEN
!              BACKSPACE 2
!              READ(2,*) descriptor,equals,value
!              WRITE(6,*) descriptor,equals,value
!        END DO
!     END IF
     
! This loop runs through the 'Dynamics' information and extracts
! system time, mass, position, velocity, and energy. It can be
! extended to extract other information if the user requires.
!
     IF (descriptor=='(Dynamics') THEN
        DO WHILE(descriptor/=')Dynamics')
           READ(2,*) descriptor                       ! Reads next line
           IF (descriptor=='system_time') THEN        ! If line is 'system time'
              BACKSPACE 2                             ! rewind file one line
              READ(2,*) descriptor,equals,system_time ! read in system_time
           !END IF
           ELSE IF (descriptor=='m') THEN             ! If line is 'm'
              BACKSPACE 2                             ! rewind file one line
              READ (2,*) descriptor,equals,topm       ! read in cluster mass
           !END IF
           ELSE IF (descriptor=='r') THEN             ! If line is 'r'
              BACKSPACE 2                             ! rewind file one line
              READ (2,*) descriptor,equals,topr1,topr2,topr3 !cluster position
           !END IF
           ELSE IF (descriptor=='v') THEN             ! If line is 'v'
              BACKSPACE 2                             ! rewind file one line
              READ (2,*) descriptor,equals,topv1,topv2,topv3 !cluster velocity
           !END IF
           ELSE IF (descriptor=='t') THEN             ! If line is 't'
              BACKSPACE 2                             ! rewind file one line
              READ (2,*) descriptor,equals,topt       ! read in snapshot time
           !END IF
           ELSE IF (descriptor=='total_energy') THEN ! If line is 'total energy'
              BACKSPACE 2                            ! rewind file one line
              READ(2,*) descriptor,equals,top_energy ! read in total_energy
           END IF
        END DO                        ! ')Dynamics' has been read. End of loop.

! This loop reads in Hydro information if required.
! It will need editing to extract any information.
!     ELSE IF (descriptor=='(Hydro')
!        DO WHILE(descriptor/=')Hydro')
!           READ(2,*) descriptor
!        END DO                           ! ')Hydro' has been read. End of loop.
        
! This loop reads Star information if required
     ELSE IF (descriptor=='(Star') THEN
        DO WHILE(descriptor/=')Star')
           READ(2,*) descriptor                        ! Reads next line
           IF (descriptor=='mass_scale') THEN          ! N-body mass scale
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,munit        ! Reads in mass_scale
           !END IF
           ELSE IF (descriptor=='size_scale') THEN     ! N-body length scale
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,runit        ! Reads in size_scale
           !END IF
           ELSE IF (descriptor=='time_scale') THEN     ! N-body time scale
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,tunit        ! Reads time_scale
           END IF
        END DO                            ! ')Star' has been read. End of loop.
     
     ELSE IF (descriptor=='(Particle') THEN
        BACKSPACE 2
     END IF
  END DO             ! The next '(Particle' line has been read so topinfo ends.
  
END SUBROUTINE read_topinfo




SUBROUTINE read_multinfo(j,snapi)
!**********************************************************************!
! 
! READING MULTIPLE INFO
!
! Extract info from 'multiple particles' & place into respective arrays
!  - mmax(j,snapi), rmax(1:3,j,snapi), vmax(1:3,j,snapi) and tmax(j,snapi)
!
!**********************************************************************!
!
  USE sl_input_module
  USE parameters_module
  IMPLICIT NONE
  integer,intent(in) :: j,snapi    ! import particle & snapshot counters
!
!
  DO WHILE(descriptor/='(Particle') ! Loop until multiple info is over & new particle read in
     READ(2,*) descriptor           ! Reads in line
!
     IF(descriptor=='t') THEN
        BACKSPACE 2                 ! Rewinds file one line
        READ(2,*) descriptor,equals,tmax(j,snapi) ! Reads in time in system particle i
     !END IF
     ELSE IF(descriptor=='m') THEN
        BACKSPACE 2                 ! Rewinds file one line
        READ(2,*) descriptor,equals,mmax(j,snapi) ! Reads mass of system particle i
     !END IF
     ELSE IF(descriptor=='r') THEN
        BACKSPACE 2                 ! Rewinds file one line
        READ(2,*) descriptor,equals,rmax(1,j,snapi),&
             & rmax(2,j,snapi),rmax(3,j,snapi) ! Reads in position of com of system particle i
     !END IF
     ELSE IF(descriptor=='v') THEN
        BACKSPACE 2                 ! Rewinds file one line
        READ(2,*) descriptor,equals,vmax(1,j,snapi)&
             & ,vmax(2,j,snapi),vmax(3,j,snapi) ! Reads in velocity of cov of system particle i

     ELSE IF(descriptor=='(Particle') then
        BACKSPACE 2           ! Go back to beginning of (Particle line
     END IF                   ! ready for return to 'read_sl_out'.
  END DO
  
END SUBROUTINE read_multinfo




SUBROUTINE read_singinfo(j,snapi)
  USE sl_input_module
  USE parameters_module
  IMPLICIT NONE
  integer,intent(in) :: j,snapi    ! import particle & snapshot counters
  
  !READ(2,*) descriptor,equals,ids(j,snapi) ! Reads in ids of stars
  DO WHILE(descriptor/=')Particle') ! Loop over particle info
     READ(2,*) descriptor           ! Reads next line
!
     IF(descriptor=='t') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,tmax(j,snapi) ! Reads in star time
     !END IF
     ELSE IF(descriptor=='m') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,mmax(j,snapi) ! Reads  in star mass
     !END IF
     ELSE IF(descriptor=='r') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,rmax(1,j,snapi),&
             & rmax(2,j,snapi),rmax(3,j,snapi) ! Reads in star position
     !END IF
     ELSE IF(descriptor=='v') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,vmax(1,j,snapi),&
             & vmax(2,j,snapi),vmax(3,j,snapi) ! Reads star velocity

     ELSE IF(descriptor==')Particle') then
        BACKSPACE 2           ! Go back to beginning of (Particle line
     END IF                   ! ready for return to 'read_sl_out'.
  END DO                      ! End of star info
END SUBROUTINE read_singinfo
