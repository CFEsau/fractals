SUBROUTINE maschberger_imf(alpha3,beta,mu,mlo,mup,kdum,tmmass)

  IMPLICIT NONE
  REAL, INTENT(in) :: alpha3,beta,mu,mlo,mup
  INTEGER, INTENT(in) :: kdum
  REAL, INTENT(out) :: tmmass
  REAL :: omb,oma3
  REAL :: Gm,Gmlo,Gmup  
  REAL :: ran2
  EXTERNAL ran2

  oma3=1. - alpha3
  omb=1. - beta
  Gmlo=(1. + (mlo/mu)**(oma3))**(omb)
  Gmup=(1. + (mup/mu)**(oma3))**(omb)

  Gm=ran2(kdum)*(Gmup - Gmlo) + Gmlo
  tmmass=mu*(Gm**(1./omb) - 1.)**(1./oma3)

  RETURN
END SUBROUTINE maschberger_imf
