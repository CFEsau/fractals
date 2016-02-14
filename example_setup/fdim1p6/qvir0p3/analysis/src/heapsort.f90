   SUBROUTINE heapsort(psort,measureof,pwhichhas)
! does a heapsort (by APW)
   IMPLICIT NONE
   INTEGER, INTENT(IN)  :: psort              ! number of values to be sorted.
   double precision, INTENT(IN)  :: measureof(1:psort) ! values to be sorted.
   INTEGER, INTENT(OUT) :: pwhichhas(1:psort) ! identifier of value.
   INTEGER              :: rank               ! rank of value.
   INTEGER              :: ranknow            ! dummy rank.
   INTEGER              :: ranktest           ! dummy rank.
!
   DO rank=2,psort                ! THIS DO-LOOP BUILDS THE BINARY HEAP
     ranknow=rank
1    IF (ranknow==1) CYCLE
     ranktest=ranknow/2
     IF (measureof(pwhichhas(ranktest))>=measureof(pwhichhas(ranknow))) CYCLE
     CALL swapi(pwhichhas(ranknow),pwhichhas(ranktest))
     ranknow=ranktest
     GOTO 1
   END DO
!
   DO rank=psort,2,-1             ! AND THIS DO-LOOP INVERTS THE BINARY HEAP
     CALL swapi(pwhichhas(rank),pwhichhas(1))
     ranknow=1
2    ranktest=2*ranknow
     IF (ranktest>=rank) CYCLE
     IF ((measureof(pwhichhas(ranktest+1))>measureof(pwhichhas(ranktest))) &
                            & .AND.(ranktest+1<rank)) ranktest=ranktest+1
     IF (measureof(pwhichhas(ranktest))<=measureof(pwhichhas(ranknow))) CYCLE
     CALL swapi(pwhichhas(ranknow),pwhichhas(ranktest))
     ranknow=ranktest
     GOTO 2
   END DO
!
   RETURN
   END SUBROUTINE heapsort
!
! ===========================================================================
! ===========================================================================
!
   SUBROUTINE swapi(item1,item2)
!
   IMPLICIT NONE
   INTEGER :: item0,item1,item2
!
   item0=item1; item1=item2; item2=item0
!
   RETURN
   END SUBROUTINE swapi

