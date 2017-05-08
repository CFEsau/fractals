!******************************************************************************!
! 
! read_sl_sub.f90
!
! This contains three subroutines used for read_sl
!
!******************************************************************************!

SUBROUTINE read_topinfo
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
  USE sl_input_module
  IMPLICIT NONE
!
! Start reading file
  READ(2,*) descriptor                ! Reads line
  IF (descriptor=='N') THEN           ! If line starts with 'N ='
     BACKSPACE 2                      ! Rewinds file one line
     READ(2,*) descriptor,equals,nTop ! Reads in nTop (Number of stars in topinfo)
  END IF
!
  DO WHILE(descriptor/='(Particle')   ! Topinfo ends when the next '(Particle' line is read
     READ(2,*) descriptor             ! Reads next line
!
     IF (descriptor=='(Dynamics') THEN                 ! If we're in the 'Dynamics' part of topinfo
        DO WHILE(descriptor/=')Dynamics')              ! Loop until we're out of the 'Dynamics' part of topinfo
           READ(2,*) descriptor                        ! Reads next line
           IF (descriptor=='system_time') THEN        
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,system_time  ! Reads in system_time
           END IF
           IF (descriptor=='total_energy') THEN 
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,top_energy   ! Reads in total_energy
           END IF
        END DO                                         ! End loop
     END IF
!
     IF (descriptor=='(Star') THEN                     ! If we're in the 'Star' part of topinfo
        DO WHILE(descriptor/=')Star')                  ! Loop until we're out of the 'Dynamics' part of topinfo
           READ(2,*) descriptor                        ! Reads next line
           IF (descriptor=='mass_scale') THEN 
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,mass_scale   ! Reads in mass_scale
           END IF
           IF (descriptor=='size_scale') THEN 
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,size_scale   ! Reads in size_scale
           END IF
           IF (descriptor=='time_scale') THEN 
              BACKSPACE 2                              ! Rewinds file one line
              READ(2,*) descriptor,equals,time_scale   ! Reads time_scale
!
! Print these variables to the terminal and to output file out_topinfo.
              !PRINT *, snapnum,nTop,system_time,top_energy,mass_scale,size_scale,time_scale
!!$              WRITE(3,*) snapnum,nTop,system_time,top_energy,mass_scale,size_scale,time_scale 
           END IF
        END DO                                         ! End loop
     END IF
     IF (descriptor=='(Particle') BACKSPACE 2
  END DO ! Ends loop: The next '(Particle' line has been read so topinfo ends.
END SUBROUTINE read_topinfo

SUBROUTINE read_multinfo
!******************************************************************************!
! 
! READING MULTIPLE INFO
!
! This will hopefully take care of binaries, triples, quadruples etc.
!
!******************************************************************************!
!
  USE sl_input_module
  IMPLICIT NONE
!
! First, need to find out how many stars are in the multiple.
  READ(2,*) descriptor,equals,partname ! Reads one line
  READ(2,*) descriptor ! Reads in line, hopefully 'N =' line
  IF(descriptor=='N') THEN
     BACKSPACE 2                                                     ! Rewinds file one line
     READ(2,*) descriptor,equals,mult_nstars(nMult(snapnum),snapnum) ! Reads in nMult (number of stars in system)  
     BACKSPACE 2
     BACKSPACE 2 ! Rewinds file two lines to 'name =' line.
  END IF

  IF(mult_nstars(nMult(snapnum),snapnum)==3) CALL read_tripleinfo
  IF(mult_nstars(nMult(snapnum),snapnum)==2) CALL read_binaryinfo

END SUBROUTINE
  
SUBROUTINE read_tripleinfo
!******************************************************************************!
! 
! READING TRIPLE INFO
!
! This will hopefully take care of triples
!
!******************************************************************************!
!
  USE sl_input_module
  IMPLICIT NONE
  REAL :: bin_t,bin_m
  REAL, DIMENSION(3) :: bin_r,bin_v
  INTEGER, DIMENSION(2) :: bin_id
  INTEGER :: i
  CHARACTER*2 :: ptest
  bin_t=0
  bin_m=0
  bin_r=0
  bin_v=0
  bin_id=0
! ! Find IDs of star1, star2 and star3
  READ(2,*) descriptor,equals,cstar(1),cstar(2),cstar(3) ! Reads in line with IDs of three stars
! The next few lines format the IDs and stick them into array mult_ids
  ! Need to know here whether there are two opening parentheses on the first one, or
  ! two closing parentheses on the last one.
  ptest=cstar(1)
  IF(ptest=="((") THEN
     WRITE(cstar_len(1),'(I1)') LEN_TRIM(cstar(1))-2
     WRITE(cstar_len(2),'(I1)') LEN_TRIM(cstar(2))-1
     WRITE(cstar_len(3),'(I1)') LEN_TRIM(cstar(3))-1
     READ(cstar(1),"(A1,A1,I"//TRIM(cstar_len(1))//")") parop,parop,mult_ids(1,nMult(snapnum),snapnum)
     READ(cstar(2),"(I"//TRIM(cstar_len(2))//",A1)") mult_ids(2,nMult(snapnum),snapnum),parcl
     READ(cstar(3),"(I"//TRIM(cstar_len(3))//",A1)") mult_ids(3,nMult(snapnum),snapnum),parcl
  ELSE
     WRITE(cstar_len(1),'(I1)') LEN_TRIM(cstar(1))-1
     WRITE(cstar_len(2),'(I1)') LEN_TRIM(cstar(2))-1
     WRITE(cstar_len(3),'(I1)') LEN_TRIM(cstar(3))-2
     READ(cstar(1),"(A1,I"//TRIM(cstar_len(1))//")") parop,mult_ids(1,nMult(snapnum),snapnum)
     READ(cstar(2),"(A1,I"//TRIM(cstar_len(2))//")") parop,mult_ids(2,nMult(snapunm),snapnum)
     READ(cstar(3),"(I"//TRIM(cstar_len(3))//",A1,A1)") mult_ids(3,nMult(snapnum),snapnum),parcl,parcl
  END IF
!
  DO WHILE(descriptor/='(Particle')                                ! Loops until the triple info is over
     READ(2,*) descriptor                                          ! Reads in line
!
     IF(descriptor=='t') THEN
        BACKSPACE 2                                                ! Rewinds file one line
        READ(2,*) descriptor,equals,mult_t(nMult(snapnum),snapnum) ! Reads in mult_t (time in system) 
     END IF
!
     IF(descriptor=='m') THEN
        BACKSPACE 2                                                ! Rewinds file one line
        READ(2,*) descriptor,equals,mult_m(nMult(snapnum),snapnum) ! Reads mult_m (mass of system)
     END IF
!
     IF(descriptor=='r') THEN
        BACKSPACE 2                                                              ! Rewinds file one line
        READ(2,*) descriptor,equals,mult_r(1,nMult(snapunm),snapnum),&
             & mult_r(2,nMult(snapnum),snapnum),mult_r(3,nMult(snapnum),snapnum) ! Reads in mult_r (position of com of system)
     END IF
!
     IF(descriptor=='v') THEN
        BACKSPACE 2                                                               ! Rewinds file one line
        READ(2,*) descriptor,equals,mult_v(1,nMult(snapnum),snapnum)&
             & ,mult_v(2,nMult(snapnum),snapnum),mult_v(3,nMult(snapnum),snapnum) ! Reads in mult_v (velocity of cov of system)
        !PRINT *, mult_ids(1:3,nMult(snapnum),snapnum), mult_r(1:3,nMult(snapnum),snapnum), mult_v(1:3,nMult(snapnum),snapnum)
     END IF
     !IF (descriptor=='(Particle') BACKSPACE 2
  END DO ! Ends loop because data for that triple is over.
! So the triple part is over, but there'll now be a binary part underneath!
! This is because the binary is the child of the triple. I need to make sure that
! the overall number of multiples doesn't go up (because it's a triple not a binary)
! but I need:
! com of star1= r(star1) + com(binary) + com(triple)
! com of star2= r(star2) + com(binary) + com(triple)
! com of star3= r(star3) + com(triple)
!
  IF((ptest=="((")) THEN
! Find IDs of star1 and star2
     READ(2,*) descriptor,equals,cstar(1),cstar(2) ! Reads in line with IDs of both stars
     !PRINT *, ptest
! The next few lines format the IDs and stick them into array mult_ids
     WRITE(cstar_len(1),'(I1)') LEN_TRIM(cstar(1))-1
     WRITE(cstar_len(2),'(I1)') LEN_TRIM(cstar(2))-1
     READ(cstar(1),"(A1,I"//TRIM(cstar_len(1))//")") parop,bin_id(1)
     READ(cstar(2),"(I"//TRIM(cstar_len(2))//",A1)") bin_id(2),parcl
!
     DO WHILE(descriptor/='(Particle')                                ! Loops until the binary info is over.
        READ(2,*) descriptor                                          ! Reads in line
!
        IF(descriptor=='t') THEN
           BACKSPACE 2                                                ! Rewinds file one line
           READ(2,*) descriptor,equals,bin_t ! Reads in bin_t (time in system) 
        END IF
!
        IF(descriptor=='m') THEN
           BACKSPACE 2                                                ! Rewinds file one line
           READ(2,*) descriptor,equals,bin_m ! Reads bin_m (mass of system)
        END IF
!
        IF(descriptor=='r') THEN
           BACKSPACE 2                                                              ! Rewinds file one line
           READ(2,*) descriptor,equals,bin_r(1),&
                & bin_r(2),bin_r(3) ! Reads in bin_r (position of com of system)
        END IF
!
        IF(descriptor=='v') THEN
           BACKSPACE 2                                                               ! Rewinds file one line
           READ(2,*) descriptor,equals,bin_v(1)&
                & ,bin_v(2),bin_v(3) ! Reads in bin_v (velocity of cov of system)
        !PRINT *, bin_id(1:2), bin_r(1:3), bin_v(1:3)
        END IF
     !IF (descriptor=='(Particle') BACKSPACE 2
     END DO ! Ends loop because data for that binary info is over.
!
! Now the three stars in the triple system...
     i=0
     DO WHILE(i<mult_nstars(nMult(snapnum),snapnum)) ! This loop cycles through the stars in the multiple system.
        READ(2,*) descriptor                         ! Reads line
!
        IF(descriptor=='i') THEN ! It's a star if it has an id number.
           i=i+1
           nstars(snapnum)=nstars(snapnum)+1 ! Add to the overall number of stars.
           BACKSPACE 2                                              ! Rewinds file one line
           READ(2,*) descriptor,equals,ids(nstars(snapnum),snapnum) ! Reads in id number of star
!
           DO WHILE(descriptor/=')Particle') ! Loop. Star info ends at next ')Particle' line
              READ(2,*) descriptor           ! Reads line
!
              IF(descriptor=='t') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_t(nstars(snapnum),snapnum) ! Reads in star time
              END IF
!
              IF(descriptor=='m') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_m(nstars(snapnum),snapnum) ! Reads in star mass
              END IF
!
              IF(descriptor=='r') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_r(1,nstars(snapnum),snapnum),&
                      & star_r(2,nstars(snapnum),snapnum),star_r(3,nstars(snapnum),snapnum) ! Reads in position relative to com
              !PRINT *, ids(nstars(snapnum),snapnum), star_r(1:3,nstars(snapnum),snapnum)
                 IF(ids(nstars(snapnum),snapnum)==bin_id(1).OR. ids(nstars(snapnum),snapnum)==bin_id(2)) THEN
                    star_r(1:3,nstars(snapnum),snapnum)=star_r(1:3,nstars(snapnum),snapnum)+&
                         & mult_r(1:3,nMult(snapnum),snapnum) + bin_r(1:3) ! star position = [position relative to com] + [position of com]
                 ELSE 
                    star_r(1:3,nstars(snapnum),snapnum)=star_r(1:3,nstars(snapnum),snapnum)+&
                         & mult_r(1:3,nMult(snapnum),snapnum)
                 END IF
!
              END IF
!
              IF(descriptor=='v') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_v(1,nstars(snapnum),nstars),&
                      & star_v(2,nstars(snapnum),snapnum),star_v(3,nstars(snapnum),snapnum) ! Reads  in velocity relative to cov
              !PRINT *, ids(nstars(snapnum),snapnum),star_v(1:3,nstars(snapnum),snapnum)
                 IF(ids(nstars(snapnum),snapnum)==bin_id(1).OR. ids(nstars(snapnum),snapnum)==bin_id(2)) THEN
                    star_v(1:3,nstars(snapnum),snapnum)=star_v(1:3,nstars(snapnum),snapnum)+&
                         & mult_v(1:3,nMult(snapnum),snapnum) + bin_v(1:3) ! star velocity = [velocity relative to cov] + [velocity of cov]
                 ELSE
                    star_v(1:3,nMult(snapnum),snapnum)=star_v(1:3,nMult(snapnum),snapnum)+&
                         & mult_v(1:3,nMult(snapnum),snapnum)
                 END IF
              !PRINT *, ids(nstars(snapnum),snapnum),mult_r(1:3,nMult(snapnum),snapnum),mult_v(1:3,nMult(snapnum),snapnum)
              END IF
           END DO
        END IF
     END DO
     ELSE
! This part is for when the file has the first single star in the multiple listed before the binary is. Very annoying KIRA!
! First star...
        READ(2,*) descriptor                         ! Reads line
        IF(descriptor=='i') THEN ! It's a star if it has an id number.
           nstars(snapnum)=nstars(snapnum)+1 ! Add to the overall number of stars.
           BACKSPACE 2                                              ! Rewinds file one line
           READ(2,*) descriptor,equals,ids(nstars(snapnum),snapnum) ! Reads in id number of star
!
           DO WHILE(descriptor/=')Particle') ! Loop. Star info ends at next ')Particle' line
              READ(2,*) descriptor           ! Reads line
!
              IF(descriptor=='t') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_t(nstars(snapnum),snapnum) ! Reads in star time
              END IF
!
              IF(descriptor=='m') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_m(nstars(snapnum),snapnum) ! Reads in star mass
              END IF
!
              IF(descriptor=='r') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_r(1,nstars(snapnum),nstars),&
                      & star_r(2,nstars(snapnum),snapnum),star_r(3,nstars(snapnum),snapnum) ! Reads in position relative to com
              !PRINT *, ids(nstars(snapnum),snapnum), star_r(1:3,nstars(snapnum),snapnum)
                 star_r(1:3,nstars(snapnum),snapnum)=star_r(1:3,nstars(snapnum),snapnum)+&
                      & mult_r(1:3,nMult(snapnum),snapnum)
!
              END IF
!
              IF(descriptor=='v') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_v(1,nstars(snapnum),nstars),&
                      & star_v(2,nstars(snapnum),snapnum),star_v(3,nstars(snapnum),snapnum) ! Reads  in velocity relative to cov
              !PRINT *, ids(nstars(snapnum),snapnum),star_v(1:3,nstars(snapnum),snapnum)
                 star_v(1:3,nstars(snapnum),snapnum)=star_v(1:3,nstars(snapnum),snapnum)+&
                      & mult_v(1:3,nMult(snapnum),snapnum)
              !PRINT *, ids(nstars(snapnum),snapnum),mult_r(1:3,nMult(snapnum),snapnum),mult_v(1:3,nMult(snapnum),snapnum)
              END IF
           END DO
        END IF
! Now binary...
! Find IDs of star1 and star2
        READ(2,*) descriptor ! Reads in line, hopefully '(Particle' line
        READ(2,*) descriptor,equals,cstar(1),cstar(2) ! Reads in line with IDs of both stars
        !PRINT *, descriptor,equals,cstar(1),cstar(2)
! The next few lines format the IDs and stick them into array mult_ids
        WRITE(cstar_len(1),'(I1)') LEN_TRIM(cstar(1))-1
        WRITE(cstar_len(2),'(I1)') LEN_TRIM(cstar(2))-1
        READ(cstar(1),"(A1,I"//TRIM(cstar_len(1))//")") parop,bin_id(1)
        READ(cstar(2),"(I"//TRIM(cstar_len(2))//",A1)") bin_id(2),parcl
!
        DO WHILE(descriptor/='(Particle')                                ! Loops until the binary info is over.
           READ(2,*) descriptor                                          ! Reads in line
!
           IF(descriptor=='t') THEN
              BACKSPACE 2                                                ! Rewinds file one line
              READ(2,*) descriptor,equals,bin_t ! Reads in bin_t (time in system) 
           END IF
!
           IF(descriptor=='m') THEN
              BACKSPACE 2                                                ! Rewinds file one line
              READ(2,*) descriptor,equals,bin_m ! Reads bin_m (mass of system)
           END IF
!
           IF(descriptor=='r') THEN
              BACKSPACE 2                                                              ! Rewinds file one line
              READ(2,*) descriptor,equals,bin_r(1),&
                   & bin_r(2),bin_r(3) ! Reads in bin_r (position of com of system)
           END IF
!
           IF(descriptor=='v') THEN
              BACKSPACE 2                                                               ! Rewinds file one line
              READ(2,*) descriptor,equals,bin_v(1)&
                   & ,bin_v(2),bin_v(3) ! Reads in bin_v (velocity of cov of system)
              !PRINT *, bin_id(1:2), bin_r(1:3), bin_v(1:3)
           END IF
           !IF (descriptor=='(Particle') BACKSPACE 2
        END DO ! Ends loop because data for that binary info is over.
! now final two stars...
        i=0
        DO WHILE(i<2) ! This loop cycles through the stars in the multiple system.
           READ(2,*) descriptor                         ! Reads line
!
           IF(descriptor=='i') THEN ! It's a star if it has an id number.
              i=i+1
              nstars(snapnum)=nstars(snapnum)+1 ! Add to the overall number of stars.
              BACKSPACE 2                                              ! Rewinds file one line
              READ(2,*) descriptor,equals,ids(nstars(snapnum),snapnum) ! Reads in id number of star
!
              DO WHILE(descriptor/=')Particle') ! Loop. Star info ends at next ')Particle' line
                 READ(2,*) descriptor           ! Reads line
!
                 IF(descriptor=='t') THEN
                    BACKSPACE 2                ! Rewinds file one line
                    READ(2,*) descriptor,equals,star_t(nstars(snapnum),snapnum) ! Reads in star time
                 END IF
!
                 IF(descriptor=='m') THEN
                    BACKSPACE 2                ! Rewinds file one line
                    READ(2,*) descriptor,equals,star_m(nstars(snapnum),snapnum) ! Reads in star mass
                 END IF
!
                 IF(descriptor=='r') THEN
                    BACKSPACE 2                ! Rewinds file one line
                    READ(2,*) descriptor,equals,star_r(1,nstars(snapnum),nstars),&
                         & star_r(2,nstars(snapnum),snapnum),star_r(3,nstars(snapnum),snapnum) ! Reads in position relative to com
              !PRINT *, ids(nstars(snapnum),snapnum), star_r(1:3,nstars(snapnum),snapnum)
                    IF(ids(nstars(snapnum),snapnum)==bin_id(1).OR. ids(nstars(snapnum),snapnum)==bin_id(2)) THEN
                       star_r(1:3,nstars(snapnum),snapnum)=star_r(1:3,nstars(snapnum),snapnum)+&
                            & mult_r(1:3,nMult(snapnum),snapnum) + bin_r(1:3) ! star position = [position relative to com] + [position of com]
                    ELSE 
                       star_r(1:3,nstars(snapnum),snapnum)=star_r(1:3,nstars(snapnum),snapnum)+&
                            & mult_r(1:3,nMult(snapnum),snapnum)
                    END IF
!
                 END IF
!
                 IF(descriptor=='v') THEN
                    BACKSPACE 2                ! Rewinds file one line
                    READ(2,*) descriptor,equals,star_v(1,nstars(snapnum),nstars),&
                         & star_v(2,nstars(snapnum),snapnum),star_v(3,nstars(snapnum),snapnum) ! Reads  in velocity relative to cov
                    !PRINT *, ids(nstars(snapnum),snapnum),star_v(1:3,nstars(snapnum),snapnum)
                    IF(ids(nstars(snapnum),snapnum)==bin_id(1).OR. ids(nstars(snapnum),snapnum)==bin_id(2)) THEN
                       star_v(1:3,nstars(snapnum),snapnum)=star_v(1:3,nstars(snapnum),snapnum)+&
                            & mult_v(1:3,nMult(snapnum),snapnum) + bin_v(1:3) ! star velocity = [velocity relative to cov] + [velocity of cov]
                    ELSE
                       star_v(1:3,nstars(snapnum),snapnum)=star_v(1:3,nstars(snapnum),snapnum)+&
                            & mult_v(1:3,nMult(snapnum),snapnum)
                    END IF
                    !PRINT *, ids(nstars(snapnum),snapnum),mult_r(1:3,nMult(snapnum),snapnum),mult_v(1:3,nMult(snapnum),snapnum)
                 END IF
              END DO
           END IF
        END DO
     END IF
END SUBROUTINE

SUBROUTINE read_binaryinfo
  USE sl_input_module
  IMPLICIT NONE
  INTEGER :: i
!
! So, for a binary:
  IF(mult_nstars(nMult(snapnum),snapnum)==2) THEN
! Find IDs of star1 and star2
     READ(2,*) descriptor,equals,cstar(1),cstar(2) ! Reads in line with IDs of both stars
! The next few lines format the IDs and stick them into array mult_ids
     WRITE(cstar_len(1),'(I1)') LEN_TRIM(cstar(1))-1
     WRITE(cstar_len(2),'(I1)') LEN_TRIM(cstar(2))-1
     READ(cstar(1),"(A1,I"//TRIM(cstar_len(1))//")") parop,mult_ids(1,nMult(snapnum),snapnum)
     READ(cstar(2),"(I"//TRIM(cstar_len(2))//",A1)") mult_ids(2,nMult(snapnum),snapnum),parcl
!
     DO WHILE(descriptor/=')Dynamics')                                ! Loops until the binary info is over.
        READ(2,*) descriptor                                          ! Reads in line
!
        IF(descriptor=='t') THEN
           BACKSPACE 2                                                ! Rewinds file one line
           READ(2,*) descriptor,equals,mult_t(nMult(snapnum),snapnum) ! Reads in mult_t (time in system) 
        END IF
!
        IF(descriptor=='m') THEN
           BACKSPACE 2                                                ! Rewinds file one line
           READ(2,*) descriptor,equals,mult_m(nMult(snapnum),snapnum) ! Reads mult_m (mass of system)
        END IF
!
        IF(descriptor=='r') THEN
           BACKSPACE 2                                                              ! Rewinds file one line
           READ(2,*) descriptor,equals,mult_r(1,nMult(snapnum),snapnum),&
                & mult_r(2,nMult(snapnum),snapnum),mult_r(3,nMult(snapnum),snapnum) ! Reads in mult_r (position of com of system)
        END IF
!
        IF(descriptor=='v') THEN
           BACKSPACE 2                                                               ! Rewinds file one line
           READ(2,*) descriptor,equals,mult_v(1,nMult(snapnum),snapnum)&
                & ,mult_v(2,nMult(snapnum),snapnum),mult_v(3,nMult(snapnum),snapnum) ! Reads in mult_v (velocity of cov of system)
        END IF
     END DO ! Ends loop because data for that binary info is over.
     i=0
     DO WHILE(i<mult_nstars(nMult(snapnum),snapnum)) ! This loop cycles through the stars in the multiple system.
        READ(2,*) descriptor                         ! Reads line
!
        IF(descriptor=='i') THEN ! It's a star if it has an id number.
           i=i+1
           nstars(snapnum)=nstars(snapnum)+1 ! Add to the overall number of stars.
           BACKSPACE 2                                              ! Rewinds file one line
           READ(2,*) descriptor,equals,ids(nstars(snapnum),snapnum) ! Reads in id number of star
!
           DO WHILE(descriptor/=')Particle') ! Loop. Star info ends at next ')Particle' line
              READ(2,*) descriptor           ! Reads line

              IF(descriptor=='t') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_t(nstars(snapnum),snapnum) ! Reads in star time
              END IF
!
              IF(descriptor=='m') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_m(nstars(snapnum),snapnum) ! Reads in star mass
              END IF
!
              IF(descriptor=='r') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_r(1,nstars(snapnum),nstars),&
                      & star_r(2,nstars(snapnum),snapnum),star_r(3,nstars(snapnum),snapnum)
                 star_r(1:3,nstars(snapnum),snapnum)=star_r(1:3,nstars(snapnum),snapnum)+&
                      & mult_r(1:3,nMult(snapnum),snapnum)
              END IF
!
              IF(descriptor=='v') THEN
                 BACKSPACE 2                ! Rewinds file one line
                 READ(2,*) descriptor,equals,star_v(1,nstars(snapnum),nstars),&
                      & star_v(2,nstars(snapnum),snapnum),star_v(3,nstars(snapnum),snapnum)
                 star_v(1:3,nstars(snapnum),snapnum)=star_v(1:3,nstars(snapnum),snapnum)+&
                      & mult_v(1:3,nMult(snapnum),snapnum)
              END IF
           END DO
        END IF
     END DO
  END IF ! end of binary star
END SUBROUTINE
  
SUBROUTINE read_singinfo
  USE sl_input_module
  IMPLICIT NONE
  READ(2,*) descriptor,equals,ids(nstars(snapnum),snapnum) ! Reads in star id
  DO WHILE(descriptor/=')Particle') ! Loop : Star info ends at next ')Particle' line.
     READ(2,*) descriptor           ! Reads next line
!
     IF(descriptor=='t') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,star_t(nstars(snapnum),snapnum) ! Reads in star time
     END IF
!
     IF(descriptor=='m') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,star_m(nstars(snapnum),snapnum) ! Reads  in star mass
     END IF
!
     IF(descriptor=='r') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,star_r(1,nstars(snapnum),nstars),&
             & star_r(2,nstars(snapnum),snapnum),star_r(3,nstars(snapnum),snapnum) ! Reads in star position
     END IF
!
     IF(descriptor=='v') THEN
        BACKSPACE 2                ! Rewinds file one line
        READ(2,*) descriptor,equals,star_v(1,nstars(snapnum),nstars),&
             & star_v(2,nstars(snapnum),snapnum),star_v(3,nstars(snapnum),snapnum) ! Reads star velocity
     END IF
  END DO ! End of star info
END SUBROUTINE
