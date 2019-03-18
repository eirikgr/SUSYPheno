****************************************************************
*** function dshmj: line of sight integral which enters in the ***
*** computation of the gamma-ray and neutrino fluxes from    ***
*** pair annihilations of wimps in the halo.                 ***
***                                                          ***
*** see definition in e.g. bergstrom et al.,                 ***
***   phys. rev d59 (1999) 043506                            ***
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

      real*8 function dshmj(cospsi0in)
      implicit none
      include 'dshmcom.h'
      include 'dsucmh.h'
      integer ier,iord,last,limit,neval
      real*8 abserr,alist,blist,elist,epsabs,epsrel,result,rlist
      dimension alist(5000),blist(5000),elist(5000),iord(5000),
     & rlist(5000)
      real*8 cospsi0in,par,sd,logrmin,logrmax,inter,locinf,
     & locsup
      integer k
      real*8 cospsi0
      common/dshmjcom/cospsi0
      external dshmjpar1
      cospsi0=cospsi0in
      epsabs=1.d-10     !numerical accuracy
      epsrel=1.d-10
      limit=5000
      if (halotype .eq. 'ucmh') then
        logrmax=dlog10(Rh_ucmh)
        if (rhcut .lt. Rh_ucmh) then
          sd=rhcut
        else
          sd=0.1d0*Rh_ucmh
        endif
      else
        sd=1.d-5
        logrmax=dlog10(r_0)
      endif
      logrmin=dlog10(sd)
      inter=(logrmax-logrmin)/6.d0
      par=0.d0
      do k=0,5
        locinf=r_0-10.d0**(logrmax-inter*k)
        locsup=r_0-10.d0**(logrmax-inter*(k+1)) 
        call dqagse(dshmjpar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        par=par+result
      enddo
      locinf=r_0-sd
      locsup=r_0+sd 
      call dqagse(dshmjpar1,locinf,locsup,epsabs,epsrel,limit,
     &  result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
      logrmin=dlog10(sd)
      if (halotype .eq. 'ucmh') then
        logrmax=dlog10(20.d0)
      else
        logrmax=dlog10(100.d0)
      endif
      par=par+result
      inter=(logrmax-logrmin)/7.d0
      do k=0,6
        locinf=r_0+10.d0**(logrmin+inter*k)
        locsup=r_0+10.d0**(logrmin+inter*(k+1)) 
        call dqagse(dshmjpar1,locinf,locsup,epsabs,epsrel,limit,
     &    result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
        par=par+result
      enddo
      dshmj=par/8.5d0*(rho0/0.3d0)**2
      return
      end


