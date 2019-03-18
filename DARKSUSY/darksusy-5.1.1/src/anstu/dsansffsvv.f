*********************************************************
*** subroutine dsansffsvv                             ***
*** fermion + fermion -> gauge boson + gauge boson in ***
*** s-channel scalar exchange (index k)               ***
*** 1 - arrow in, 2 - arrow out, k intermediate       ***
*** this code is computer generated by reduce         ***
*** and gentran.                                      ***
*** author: joakim edsjo, edsjo@physto.se             ***
*********************************************************

      subroutine dsansffsvv(p,costheta,kp1,kp2,kpk,kp3,kp4)
      implicit none

      include 'dsmssm.h'
      include 'dsandiacom.h'
      integer kp1,kp2,kpk,kp3,kp4
      real*8 p,costheta
      complex*16 dh,gv2,ga2

      call dsankinvar(p,costheta,kp1,kp2,kpk,kp3,kp4)
      if (s.lt.(mp3+mp4)**2) return
      dh=1.0d0/dcmplx(mk**2-s,-width(kpk)*mk)
      gv2=gl(kpk,kp3,kp4)*
     &  (gl(kpk,kp2,kp1)+gr(kpk,kp2,kp1))
      ga2=gl(kpk,kp3,kp4)*
     &  (gl(kpk,kp2,kp1)-gr(kpk,kp2,kp1))

      if (kp3.eq.kgamma) then
        if (kp4.eq.kgamma) then

      aa(0,0,-1,-1)=aa(0,0,-1,-1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(0,0,1,1)=aa(0,0,1,1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(1,0,-1,-1)=aa(1,0,-1,-1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)
      aa(1,0,1,1)=aa(1,0,1,1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)

        else

      aa(0,0,-1,-1)=aa(0,0,-1,-1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(0,0,1,1)=aa(0,0,1,1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(1,0,-1,-1)=aa(1,0,-1,-1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)
      aa(1,0,1,1)=aa(1,0,1,1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)

        endif
      else
        if (kp4.eq.kgamma) then

      aa(0,0,-1,-1)=aa(0,0,-1,-1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(0,0,1,1)=aa(0,0,1,1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(1,0,-1,-1)=aa(1,0,-1,-1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)
      aa(1,0,1,1)=aa(1,0,1,1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)

        else

      aa(0,0,-1,-1)=aa(0,0,-1,-1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(0,0,0,0)=aa(0,0,0,0)-(dsqrt(2.0d0)*(e3*e4+kk**2)*dh*epl*ga2*
     . mw)/(mp3*mp4)
      aa(0,0,1,1)=aa(0,0,1,1)+dsqrt(2.0d0)*dh*epl*ga2*mw
      aa(1,0,-1,-1)=aa(1,0,-1,-1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)
      aa(1,0,0,0)=aa(1,0,0,0)+(dsqrt(2.0d0)*(e3*e4+kk**2)*dh*gv2*mw*
     . ppl)/(mp3*mp4)
      aa(1,0,1,1)=aa(1,0,1,1)-(dsqrt(2.0d0)*dh*gv2*mw*ppl)

        endif
      endif

      return

      end


