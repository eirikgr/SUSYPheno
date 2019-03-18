!N=16 beta calculator for feeder scaling. 
!
!Author:  Pat Scott
!          pscott@imperial.ac.uk
!Date: 2012
! 
!Input:  M     M_3 non-Gaussianty parameter (dimensionless)
!        nu    ratio delta/sigma of density contrast to mass variance (dimensionless)
!Output: beta  (dimensionless)

      double precision function dsucmh_BetaF16(M,v)

      implicit none
      include 'dsmpconst.h'

      double precision, intent(IN) :: M,v

      dsucmh_BetaF16 = (M**1.3333333333333333*v*(-3. + v**2))/(4.*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (M*(-1. + v**2))/(3.*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (7.*M**2*v*(15. - 10.*v**2 + v**4))/(36.*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (M**1.6666666666666667*(3. - 6.*v**2 + v**4))/(5.*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (167.*M**2.6666666666666665*v*(-105. + 105.*v**2 - 21.*v**4 + v**6))/(960.*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (31.*M**2.3333333333333335*(-15. + 45.*v**2 - 15.*v**4 + v**6))/(168.*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (7969.*M**3.3333333333333335*v*(945. - 1260.*v**2 + 378.*v**4 - 36.*v**6 + v**8))/
     &   (50400.*exp(v**2/2.)*Sqrt(2.*pi)) + (67.*M**3*(105. - 420.*v**2 + 210.*v**4 - 28.*v**6 + v**8))/
     &   (405.*exp(v**2/2.)*Sqrt(2.*pi)) + (635347.*M**4*v*
     &     (-10395. + 17325.*v**2 - 6930.*v**4 + 990.*v**6 - 55.*v**8 + v**10))/(4.35456e6*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (67259.*M**3.6666666666666665*(-945. + 4725.*v**2 - 3150.*v**4 + 630.*v**6 - 45.*v**8 + v**10))/
     &   (443520.*exp(v**2/2.)*Sqrt(2.*pi)) + (76070549.*M**4.666666666666667*v*
     &     (135135. - 270270.*v**2 + 135135.*v**4 - 25740.*v**6 + 2145.*v**8 - 78.*v**10 + v**12))/
     &   (5.588352e8*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (474311.*M**4.333333333333333*(10395. - 62370.*v**2 + 51975.*v**4 - 13860.*v**6 + 1485.*v**8 - 66.*v**10 + v**12))/
     &   (3.3696e6*exp(v**2/2.)*Sqrt(2.*pi)) + (12759542113.*M**5.333333333333333*v*
     &     (-2027025. + 4729725.*v**2 - 2837835.*v**4 + 675675.*v**6 - 75075.*v**8 + 4095.*v**10 - 105.*v**12 + v**14))/
     &   (9.96323328e10*exp(v**2/2.)*Sqrt(2.*pi)) + 
     &  (947790673.*M**5*(-135135. + 945945.*v**2 - 945945.*v**4 + 315315.*v**6 - 45045.*v**8 + 3003.*v**10 - 91.*v**12 + 
     &       v**14))/(7.185024e9*exp(v**2/2.)*Sqrt(2.*pi)) + Erfc(v/Sqrt(2.))

      end function dsucmh_BetaF16

