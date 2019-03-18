****************************************************************
*** function dshmjave: average over the solid angle          ***
*** delta (sr) of the function dshmj(cospsi0)                ***
***                                                          ***
*** dshmj is the line of sight integral which enters in the  ***
*** computation of the gamma-ray and neutrino fluxes from    ***
*** pair annihilations of wimps in the halo.                 ***
***                                                          ***
*** see definition in e.g. bergstrom et al.,                 ***
***   phys. rev. d59 (1999) 043506                           ***
*** in case of the many unresolved clump scenario the term   ***
*** fdelta is factorized out                                 *** 
***                                                          ***  
*** it is valid for a spherical dark matter halo             ***
*** psi0 is the angle between direction of observation       ***
*** and the direction of the galactic center; cospsi0 is     ***
*** its cosine.                                              ***
***                                                          ***
*** author: piero ullio (piero@tapir.caltech.edu)            ***
*** date: 00-07-13                                           *** 
*** modified: 09-08-08 pat scott (pat@fysik.su.se)           ***
****************************************************************


      real*8 function dshmjave(cospsi0in,deltain)
      implicit none
      include 'dshmcom.h'
      include 'dsucmh.h'
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsrel,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      real*8 cospsi0in,deltain,par,sd,logrmin,logrmax,inter,locinf,
     & locsup,checkres
      integer k
      real*8 cospsi0,delta,rr
      common/dshmjavecom/cospsi0,delta,rr
      external dshmjavepar1
      real*8 epsabs
      integer hmjaveset,ndigits
      common/hmjavesetcom/epsabs,hmjaveset,ndigits
      if(hmjaveset.ne.123456) then
         ndigits=5
         epsabs=1.d-10          !numerical accuracy
         hmjaveset=123456
      endif
      limit=5000
      cospsi0=cospsi0in
      delta=deltain
      par=0.d0
      if (halotype .eq. 'ucmh') then
        logrmax=dlog10(Rh_ucmh)
        if (rhcut .lt. Rh_ucmh) then
          sd=1.d1*rhcut
        else
          sd=0.1d0*Rh_ucmh
        endif
      else
        sd=1.d-5
        logrmax=dlog10(r_0)
      endif
      epsrel=10.d0**(-ndigits)
      logrmin=dlog10(sd)
      inter=(logrmax-logrmin)/6.d0
      par=0.d0
      do k=0,5
        locinf=r_0-10.d0**(logrmax-inter*k)
        locsup=r_0-10.d0**(logrmax-inter*(k+1))
 10     call dqagse(dshmjavepar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        checkres=par+result
        if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
          goto 10
        elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
        endif  
        par=par+result
      enddo
      locinf=r_0-sd
      locsup=r_0+sd 
 11   call dqagse(dshmjavepar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      checkres=par+result
      if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
        epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
        goto 11
      elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
        epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
      endif  
      par=par+result
      logrmin=dlog10(sd)
      if (halotype .eq. 'ucmh') then
        logrmax=dlog10(20.d0)
      else
        logrmax=dlog10(100.d0)
      endif
      inter=(logrmax-logrmin)/7.d0
      do k=0,6
        locinf=r_0+10.d0**(logrmin+inter*k)
        locsup=r_0+10.d0**(logrmin+inter*(k+1)) 
 12     call dqagse(dshmjavepar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        checkres=par+result
        if(dabs(checkres)/10.d0**ndigits.lt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2) !too low accuracy
          goto 12
        elseif(dabs(checkres)/10.d0**(ndigits+2).gt.epsabs) then
          epsabs=dabs(checkres)/10.d0**(ndigits+2)  !too high accuracy 
        endif  
        par=par+result
      enddo
      dshmjave=par/8.5d0/delta*(rho0/0.3d0)**2
      return
      end









