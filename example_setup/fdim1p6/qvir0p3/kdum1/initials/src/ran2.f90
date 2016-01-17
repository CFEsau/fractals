      FUNCTION ran2(idum)
!  random number generator
      PARAMETER (m=714025,ia=1366,ic=150889,rm=1./m)
      COMMON/rand/  iy,iff,ir(97) 
!
      IF (idum.LT.0.OR.iff.EQ.0) THEN
          iff = 1
          idum = MOD(ic-idum,m)
          DO 11 j = 1,97
              idum = MOD(ia*idum+ic,m)
              ir(j) = idum
   11     CONTINUE
          idum = MOD(ia*idum+ic,m)
          iy = idum
      END IF
      j = 1 + (97*iy)/m
      IF (j.GT.97.OR.j.LT.1) WRITE (6,12)  j, idum
   12 FORMAT (/,'  troubles in ran2   j idum ',2i12)
      iy = ir(j)
      ran2 = iy*rm
      idum = MOD(ia*idum+ic,m)
      ir(j) = idum
!
      RETURN
      END FUNCTION ran2
