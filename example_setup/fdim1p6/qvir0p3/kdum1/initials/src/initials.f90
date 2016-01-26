!******************************************************************************!
!
! 21/01/16 Claire Esau
!
! Stolen from Dan, modified binary pairing method
!
!******************************************************************************!
!
PROGRAM initials
  USE input_module
  USE properties_module
  USE output_module
  USE constants_module
  IMPLICIT NONE
! Define variables:
! Iterators i, j and k
! kdum = a random seed
  INTEGER :: i,j,k,kdum
! Random number generator RAN2
  REAL :: ran2, random, test
  EXTERNAL ran2
!
!******************************************************************************!
! 
! Declare constants
  pi=4.*ATAN(1.)
  twopi=8.*ATAN(1.)
!
!******************************************************************************!
!
! File input: Read in initial parameters from input file.
  !WRITE(6,*) 'Read in parameters - choose input file:'
  !READ(5,*) infilename
  infilename='input.dat'
  WRITE(6,*) 'Chosen input file is: ',infilename
! Open input file
  OPEN(1,file=infilename,status='old')
  READ(1,*) nsys           ! number of systems (not stars)
  READ(1,'(a)') distrib    ! mass distribution
  READ(1,*) fdim           ! fractal dimension
  READ(1,*) runit          ! unit scaling. For Plummer sphere, halfmass radius (in pc) = approximately 1.3x the Plummer radius
  READ(1,*) tend           ! simulation length (Myr)
  READ(1,*) tout           ! snapshot interval (Myr)
  READ(1,'(a)') filestem   ! used to make filenames
  READ(1,*) fbinary        ! fraction of systems in binaries
  READ(1,'(a)') pairing    ! binary pairing method - ratio or imf
  CLOSE(1)
!
! Loop over as many simulations as you wish to create.
  DO k=0,0
! File Output: Creates the names of each output file.
  IF (k<10) THEN 
     WRITE(id,'(i1)')k
     outfilestem='out'//filestem//'0'//id
     outfilename=outfilestem//'.dat'
     icfilestem='ic'//filestem//'0'//id
     icfilename=icfilestem//'.sl'
     runfilestem='run'//filestem//'0'//id
     runfilename=runfilestem//'.sl'
     scriptname='scr_'//filestem//'0'//id
     restartfilestem='RS_'//filestem//'0'//id
     restartfilename=restartfilestem//'.sl'
  ELSE
     WRITE(id,'(i2)')k
     outfilestem='out'//filestem//id
     outfilename=outfilestem//'.dat'
     icfilestem='ic'//filestem//id
     icfilename=icfilestem//'.sl'
     runfilestem='run'//filestem//id
     runfilename=runfilestem//'.sl'
     scriptname='scr_'//filestem//id
     restartfilestem='RS_'//filestem//id
     restartfilename=restartfilestem//'.sl'
  END IF
!
!******************************************************************************!
!
! Assume that stars will be in triple systems at most
  nmax=nsys*3
! Manually set seed to be able to reproduce simulations
  kdum=1
! Allocate memory space for arrays
  ALLOCATE(rsys(1:3,1:nmax))
  ALLOCATE(vsys(1:3,1:nmax))
  ALLOCATE(msys(1:nmax))
  ALLOCATE(sinfo(1:nmax))
  ALLOCATE(r(1:3,1:nmax))
  ALLOCATE(v(1:3,1:nmax))
  ALLOCATE(m(1:nmax))
! Initialise arrays
  rsys=0.
  vsys=0.
  msys=0.
  r=0.
  v=0.
  m=0.
! The default number of stars in a system is 1
  sinfo=1
! Initialise variables
! The number of stars is initialised at 0 & calculated from nsys
  nstars=0
  numsingle=0
  numbinary=0
  mtotsys=0.
  mtotstars=0.
! Set system as sub-virial (0.5 for virialised)
  qvir=0.3
!
!******************************************************************************!
! 
! Open outfile: This file will give data on the macroscopic properties of the simulation
  OPEN(2,file='./output/'//outfilename,status='new')
  WRITE(2,*) 'Number of systems: ',nsys
  WRITE(2,*) 'Simulation length (Myr): ',tend
  WRITE(2,*) 'Snapshot interval (Myr): ',tout
  WRITE(2,*) 'Random number seed: ',kdum
  WRITE(2,*) 'output filename: ',outfilename 
  WRITE(2,*) 'Virial ratio: ',qvir 
!
!******************************************************************************!
!
! Finding the masses of each system!
!
  alpha=2.3
  beta=1.4
  mLower=0.1
  mUpper=50
  mu=0.2
  i=0
  j=0
  
  numsingle=INT((1.-fbinary)*DBLE(nsys))
  numbinary=INT(2.*fbinary*DBLE(nsys))
  nstars=numsingle+numbinary
  
  WRITE(2,*)'Binary fraction:',fbinary
  write(2,*)'Pairing method: ',pairing
  WRITE(2,*)'Number of systems:',nsys
  WRITE(2,*)'Number of single stars:',numsingle
  WRITE(2,*)'Number of stars in binaries:',numbinary
  WRITE(2,*)'Total number of stars:',nstars
  
  sinfo(1:numsingle)=1
  IF (nstars>numsingle) sinfo(numsingle+1:nstars)=2
  DO WHILE(i<nstars)
!i counts number of stars. j counts number of systems
     i=i+1
     j=j+1
     IF (sinfo(i)==1) THEN
        CALL maschberger_imf(alpha,beta,mu,mLower,mUpper,kdum,m(i))
        msys(j)=m(i)
        mtotstars=mtotstars+DBLE(m(i))
        mtotsys=mtotsys+DBLE(msys(j))
     ELSE IF (sinfo(i)==2) THEN
!call masch for m1:
        CALL maschberger_imf(alpha,beta,mu,mLower,mUpper,kdum,m(i))
!Let's assume that pairing will either be in upper case, lower case, or
!first character upper...
        IF (pairing=='IMF' .OR. pairing=='imf' .OR. pairing=='IMF') THEN
!add another star and call masch for m2:
           i=i+1
           CALL maschberger_imf(alpha,beta,mu,mLower,mUpper,kdum,m(i))  

        ELSE IF (pairing=='ratio' .OR. pairing=='Ratio' .OR. pairing=='RATIO') THEN
!pick random mass ratio between e.g. 0.1 & 1 for mass of secondary:
!m(2)=(RAN2(kdum*(1.-0.1))+0.1) ---> (...kdum*0.9)+0.1
           i=i+1
           m(i)=(RAN2(kdum)*0.9)+0.1
        ELSE
           STOP 'No pairing method selected!'
        END IF
!set system mass as sum of m1 and m2:
        msys(j)=m(i-1)+m(i)
        mtotstars=mtotstars+DBLE(m(i-1))+DBLE(m(i))
        mtotsys=mtotsys+DBLE(msys(j))
     ELSE
        STOP 'SINFO DOES NOT EQUAL 1 OR 2!'
     END IF
  END DO

  !check i=nstars and j=nsys
!  PRINT *,nstars, i, nsys, j
  IF(i/=nstars)STOP 'i is not equal to number of stars'
  IF(j/=nsys)STOP 'j is not equal to number of systems'
     
!
!******************************************************************************!
!
!
! Find the positions and velocities of systems and fill arrays
! rsys and vsys. For systems consisting of a single star,this
! is the position and velocity of that star. For multiple systems, 
! this is the centre of mass and the centre of velocity for that system.
!
! If distrib=p then we create a Plummer sphere.
  IF (distrib=='p' .OR. distrib=='P') THEN
     WRITE(2,*)'Distribution: Plummer sphere'
! Calling the plummer subroutine creates a plummer sphere of 
! plummer radius 1
     CALL plummer(nsys,rsys,vsys,kdum)
! as this gives a system of plummer radius 1, scale to runit here
     rsys=rsys*runit
! otherwise, if distrib=f then we make a fractal distribution with fractal dimension fdim.
  ELSE IF (distrib=='f' .OR. distrib=='F') THEN
     WRITE(2,*) 'Fractal cluster, fractal dimension:',fdim
     WRITE(2,*) 'Radius of fractal cluster',runit
     CALL makefractal(nsys,rsys,vsys,fdim,kdum)
! scale to runit to make the size of the cluster here 
     rsys=rsys*runit
  END IF
  WRITE(2,*) 'r and v values set up.'
  WRITE(2,*) 'Radius (pc): ', runit

!******************************************************************************!

  nstars=0
  numsingle=0
  numbinary=0
  mtotsys=0.
  mtotstars=0.

  DO j=1,nsys
     IF(sinfo(j)==1) THEN
! Increase the total number of stars by one
        nstars=nstars+1
! The mass of the star = the mass of the system
        m(nstars)=msys(j)
! same for position and velocity
        r(1:3,nstars)=rsys(1:3,j)
        v(1:3,nstars)=vsys(1:3,j)
! Add the system mass to the total system mass
        mtotsys=mtotsys + DBLE(msys(j))
! Add the star mass to the overall stellar mass
        mtotstars=mtotstars + DBLE(m(nstars))

! Add one to the number of single-star systems
        numsingle=numsingle+1
        
     ELSE ! Otherwise, make a binary

!Find period & hence separation.
!Period/separation distributions are found in the literature
!For binary separation between 1000-5000 AU:
!        asep=(RAN2(kdum)*4000.)+1000.
!Constant separation (in AU):
        asep=50.
        
        ! Then the period (in seconds) is given by:
        asep=asep*au
        P=SQRT(((asep**3.)*(twopi**2.))/((m(nstars+1) + m(nstars+2))*msun*G))

! Need to find the eccentricity of the system.
! Set this to zero for now.
        ecc=0

! Convert period to separation
        asep=(((P/twopi)**2.)*G*((m(nstars+1) + m(nstars+2))*msun))**(1./3.)

! Add 1 to the overall binary number
        numbinary=numbinary + 1
! Add the system mass to the total system mass
        mtotsys=mtotsys + DBLE(msys(i))
! Add the star masses to the overall stellar mass
        mtotstars=mtotstars + DBLE(m(nstars+1)) + DBLE(m(nstars+2))
!
! subroutine assign_binary uses the stellar masses and the eccentricity &
! separation of the binary to find the positions & velocities of
! the stars in the binary relative to the binary's centres of mass & velocity.
! Separation needs to be a Real and in parsecs.
        asep=asep/pc
        rsep=REAL(asep)
        CALL assign_binary(m(nstars+1),m(nstars+2),ecc,rsep,v(1:3,nstars+1), &
             &                           v(1:3,nstars+2),r(1:3,nstars+1),r(1:3,nstars+2),kdum)

! Add 2 to the total number of stars
        nstars=nstars + 2
! Binary is created!
     END IF
  END DO
! All systems should now be specified
  !PRINT *, nstars, nsys, mtotsys, mtotstars, numsingle, numbinary
!
!******************************************************************************! 
!
! Write new properties to output file
  WRITE(2,*) 'Total number of systems: ',nsys
  WRITE(2,*) 'Total number of stars: ',nstars
  WRITE(2,*) 'Total number of binaries: ',numbinary
  WRITE(2,*) 'Total mass calculated from from systems (mSun): ',mtotsys
  WRITE(2,*) 'Total mass calculated from from stars (mSun): ',mtotstars
!
!******************************************************************************! 
!
! Get the kinetic/potential/total energy:
! Kinetic energy:
  ekin=0.d0
  DO i=1,nsys
     ekin=ekin + msys(i)*(vsys(1,i)*vsys(1,i) + vsys(2,i)*vsys(2,i) + &
          &                      vsys(3,i)*vsys(3,i))
  END DO
  rkin=REAL(ekin)*0.5
! Get units right
  ekin=ekin*0.5*msun
! Potential energy:
   epot=0.d0
   DO i=1,nsys-1
      epoti=0.d0
      DO j=i+1,nsys
         rij=(rsys(1,i)-rsys(1,j))**2 + (rsys(2,i)-rsys(2,j))**2 + & 
              &              (rsys(3,i)-rsys(3,j))**2 + 1.e-10
         epoti=epoti + DBLE(msys(j)/SQRT(rij))
      END DO
      epot=epot + DBLE(msys(i))*epoti
   END DO
   rpot=REAL(epot)
! Get units right
   epot=epot*G*msun*msun/pc
! So far, we don't know what units the velocities are in.
! We can use ekin, epot and the virial ratio qvir to convert
! the velocities into units of km/s.
   DO i=1,nsys
      vsys(1:3,i)=vsys(1:3,i)*SQRT(epot*qvir/ekin)/1.e3
   END DO
!
!******************************************************************************! 
!
! Set the Centre of mass and the centre of velocity of the distribution to zero
!
! Initialise:
   com=0.d0
   cov=0.d0
   mtotsys=0.d0
! Add up the centres of mass and velocity of all the systems
   DO i=1,nsys
      com(1:3)=com(1:3) + DBLE(msys(i)*rsys(1:3,i))
      cov(1:3)=cov(1:3) + DBLE(msys(i)*vsys(1:3,i))
      mtotsys=mtotsys + DBLE(msys(i))
   END DO
! Calculate the centre of mass + velocity of the distribution
   com=com/mtotsys
   cov=cov/mtotsys
! Take this away from the positions and velocities of the systems
! This centres the distribution on zero.
   DO i=1,nsys
      rsys(1:3,i)=rsys(1:3,i) - REAL(com(1:3))
      vsys(1:3,i)=vsys(1:3,i) - REAL(cov(1:3))
   END DO
   com=0.d0
   cov=0.d0
   mtotsys=0.d0
! Find new com/cov
   DO i=1,nsys
      com(1:3)=com(1:3) + DBLE(msys(i)*rsys(1:3,i))
      cov(1:3)=cov(1:3) + DBLE(msys(i)*vsys(1:3,i))
      mtotsys=mtotsys + DBLE(msys(i))
   END DO
   com=com/mtotsys
   cov=cov/mtotsys
!
! TODO: look at the above code in more detail
!
!******************************************************************************! 
! PREPARE DATA FOR WRITING-OUT 
!
! KIRA requires the input data to be in the correct units and format.
!
!******************************************************************************! 
!
! Get N-body units:
!
! Total mass of the distribution.
   mtotstars=0.d0
   DO i=1,nstars
      mtotstars=mtotstars + DBLE(m(i))
   END DO
   munit=REAL(mtotstars)
   ! Time in Myr
   tunit=REAL(SQRT((DBLE(runit)*pc)**3./(G*DBLE(munit)*msun))/Myr)
   ! Velocity in km/s
   vunit=REAL(DBLE(runit)*pc/(DBLE(tunit)*Myr*1.d3))
!
! Convert the masses, positions and velocities of the stars and the
! systems into N-body units.
! This assumes m in msun, r in pc, v in km/s
   m=m/munit
   r=r/runit
   v=v/vunit
   vsys=vsys/vunit
   rsys=rsys/runit
   msys=msys/munit
!
!******************************************************************************! 
!
! Write new properties to output file
  WRITE(2,*) 'Total kinetic energy: ',ekin
  WRITE(2,*) 'Total potential energy: ',epot
  WRITE(2,*) 'Total energy: ',ekin-epot
  WRITE(2,*) 'Centre of mass ',com
  WRITE(2,*) 'Centre of velocity ',cov
  WRITE(2,*) 'munit: ', munit
  WRITE(2,*) 'runit: ', runit
  WRITE(2,*) 'vunit: ', vunit
  WRITE(2,*) 'tunit: ', tunit
  WRITE(2,*) '1 Myr is', 1./tunit, ' n-body units'
!
!******************************************************************************! 
!
! Specify the close encounter distance:
!
! Close encounter distance in AU:
   cenc=1.e3
! Write to output file
   WRITE(2,*) 'Close encounter distance (AU): ',cenc
! Convert to N-body units:
   cenc=cenc*au/pc
   ncenc=REAL(DBLE(cenc)/DBLE(runit))
! Convert to a string:
   WRITE(ncencchar,'(f9.5)') ncenc
!
!******************************************************************************! 
!
! KIRA command line
!
! Convert the end time from Myr to N-body units:
   ntend=tend/tunit
! Convert to a string:
   WRITE(ntendchar,'(f9.4)') ntend
! Convert the snapshot interval from Myr to N-body units:
   ntout=tout/tunit
! Convert to a string:
   WRITE(ntoutchar,'(f9.4)') ntout
! Fill the kira command string:
! -t = end time
! -D = snapshot interval
! -a = accuracy parameter
! -y = close encounter distance
! -u = disable unperturbed multiples
! -O = save (and overwrite) extra snapshot at each output
   kiracomm = 'kira -t '//ntendchar//' -D '//ntoutchar//' -a 0.01 -y'//ncencchar//' -u -O '//restartfilename//&
        &' < '//icfilename//' > '//runfilename//' 2> '//runfilename(4:12)//'err'
!
!******************************************************************************! 
!
! Write command string to output file
   WRITE(2,*) 'Kira command line: ',kiracomm
! Write script
   OPEN(7,file='./output/'//scriptname,status='new')
   WRITE(7,'(a)') '#!/bin/sh'
   WRITE(7,'(a)') '#$ -cwd'
   WRITE(7,'(a)') '#$ -q long.q'
   !WRITE(7,'(a)') 'mkdir ./'//filestem
   !WRITE(7,'(a)') 'cd ./'//filestem
   WRITE(7,'(a)') kiracomm
   CLOSE(7) 
!
!******************************************************************************! 
!
! Fake starlab format:
      OPEN(8,file='./output/'//icfilename,status='new')
!
! write header information
        WRITE(8,'(a)') '(Particle'
        WRITE(8,*) '  name = root'
        WRITE(8,*) '  N =', nstars
        WRITE(8,'(a)') '(Log'
        WRITE(8,'(a)') ')Log'
!
        WRITE(8,'(a)') '(Dynamics'
        WRITE(8,*) '  system_time  =  0'
        WRITE(8,*) '  m  =  1'
        WRITE(8,*) '  r  =  0 0 0'
        WRITE(8,*) '  v  =  0 0 0'
        WRITE(8,*) '  com_time = 0'
        WRITE(8,*) '  com_pos =',REAL(com(1)), REAL(com(2)), REAL(com(3))
        WRITE(8,*) '  com_vel =',REAL(cov(1)), REAL(cov(2)), REAL(cov(3))
        WRITE(8,*) '  total_energy =', rkin-rpot
        WRITE(8,'(a)') ')Dynamics'
!
        WRITE(8,'(a)') '(Hydro'
        WRITE(8,'(a)') ')Hydro'
        WRITE(8,'(a)') '(Star'
        WRITE(8,'(a18,2x,f17.15)') '  mass_scale     =', 1./munit
        WRITE(8,'(a18,2x,e17.11)') '  size_scale     =', 2.255e-8/runit
        WRITE(8,'(a18,2x,e17.11)') '  time_scale     =', 1./tunit
        WRITE(8,'(a)') ')Star'
!
        nstars=0
        DO i=1,nsys
! write system data
          IF (sinfo(i)==1) THEN
            nstars=nstars + 1
            WRITE(8,'(a)') '(Particle'
            WRITE(8,*) '  i =', i
            WRITE(8,*) '  N = 1'
            WRITE(8,'(a)') '(Log'
            WRITE(8,'(a)') ')Log'
            WRITE(8,'(a)') '(Dynamics'
            WRITE(8,*) '  m  =',  msys(i)
            WRITE(8,*) '  r  =',  rsys(1:3,i)
            WRITE(8,*) '  v  =',  vsys(1:3,i)
            WRITE(8,'(a)') ')Dynamics'
            WRITE(8,'(a)') '(Hydro'
            WRITE(8,'(a)') ')Hydro'
            WRITE(8,'(a)') '(Star'
            WRITE(8,'(a)') ')Star'
            WRITE(8,'(a)') ')Particle'
         ELSE
            WRITE(8,'(a)') '(Particle'
            IF (i<=9) THEN
               WRITE(8,10) '   name = (',i,',',i + 100000,')'
            ELSE IF (i>9 .AND. i<=99) THEN
               WRITE(8,11) '   name = (',i,',',i + 100000,')'
            ELSE IF (i>99 .AND. i<=999) THEN
               WRITE(8,12) '   name = (',i,',',i + 100000,')'
            ELSE IF (i>999 .AND. i<=9999) THEN
               WRITE(8,13) '   name = (',i,',',i + 100000,')'
            ELSE IF (i>9999 .AND. i<=99999) THEN
               WRITE(8,14) '   name = (',i,',',i + 100000,')'
            END IF
            WRITE(8,*) '  N =', sinfo(i)
            WRITE(8,'(a)') '(Log'
            WRITE(8,'(a)') ')Log'
            WRITE(8,'(a)') '(Dynamics'
            WRITE(8,*) '  m  =',  msys(i)
            WRITE(8,*) '  r  =',  rsys(1:3,i)
            WRITE(8,*) '  v  =',  vsys(1:3,i)
            WRITE(8,'(a)') ')Dynamics'
            WRITE(8,'(a)') '(Hydro'
            WRITE(8,'(a)') ')Hydro'
            WRITE(8,'(a)') '(Star'
            WRITE(8,'(a)') ')Star'

            DO j=1,sinfo(i)
               WRITE(8,'(a)') '(Particle'
               IF (j==1) THEN
                  WRITE(8,*) '  i =', i
                  IF (i>99999) STOP 'number of stars too big for names'
               ELSE
                  WRITE(8,*) '  i =', i + 100000
               END IF
            WRITE(8,*) '  N = 1'
            WRITE(8,'(a)') '(Log'
            WRITE(8,'(a)') ')Log'
            WRITE(8,'(a)') '(Dynamics'
            WRITE(8,*) '  m  =',  m(nstars+j)
            WRITE(8,*) '  r  =',  r(1:3,nstars+j)
            WRITE(8,*) '  v  =',  v(1:3,nstars+j)
            WRITE(8,'(a)') ')Dynamics'
            WRITE(8,'(a)') '(Hydro'
            WRITE(8,'(a)') ')Hydro'
            WRITE(8,'(a)') '(Star'
            WRITE(8,'(a)') ')Star'
            WRITE(8,'(a)') ')Particle'
           END DO
          WRITE(8,'(a)') ')Particle'
          nstars=nstars + sinfo(i)
          END IF
        END DO
!
        WRITE(8,'(a)') ')Particle'
!
      CLOSE(8)
!
! For formatting the above output:
10    FORMAT(A11,I1,A1,I6,A1)
11    FORMAT(A11,I2,A1,I6,A1)
12    FORMAT(A11,I3,A1,I6,A1)
13    FORMAT(A11,I4,A1,I6,A1)
14    FORMAT(A11,I5,A1,I6,A1)
!
!******************************************************************************! 
!
! Deallocate array memory
  DEALLOCATE(rsys)
  DEALLOCATE(vsys)
  DEALLOCATE(msys)
  DEALLOCATE(sinfo)
  DEALLOCATE(r)
  DEALLOCATE(v)
  DEALLOCATE(m)
  CLOSE(2)
END DO
END PROGRAM initials

