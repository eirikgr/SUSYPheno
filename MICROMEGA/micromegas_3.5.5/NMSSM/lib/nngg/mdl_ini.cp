*     LanHEP output produced at Tue May 18 00:00:25 2010
*     Model named 'NMSSM-Ug'

C#include "DgetNEWv0.F"
C#include "CgetNEWv0.F"

      subroutine ModelDefaults
      implicit none
      integer slharead
      real*8 slhaVal0,slhaVal1,slhaVal2
      integer slhaValExists0,slhaValExists1,slhaValExists2
      integer err
	  double precision findValW	
      real*8 Q,MSt_1,MSt_2        

#include "model.h"
 
      Q=0
      MST_1=slhaVal1('MASS',Q,1000006)
      MST_2=slhaVal1('MASS',Q,2000006)
      Q=sqrt(MST_1*MST_2)

      wtop = 0D0
      EE = 0.31344D0
      SW = 0.481D0
      s12 = 0.221D0
      s23 = 0.041D0
      s13 = 0.0035D0
      GG = 1.238D0
	  
      Q=0
      MST_1=slhaVal1('MASS',Q,1000006)
      MST_2=slhaVal1('MASS',Q,2000006)
      Q=sqrt(MST_1*MST_2)
	  
      MZ = 91.187D0
	  
      Me = 0.000511D0
      Mm = 0.1057D0
      Mu = 0.046D0
      Md = 0.046D0
      Mc = 1.42D0
      Ms = 0.2D0
      Ml = 1.777D0

      Mt = slhaVal1('SMINPUTS',Q,6)
      Mb = 2.7D0
      wt = 1.7524D0
      wZ = 2.4944D0
      wW = 2.08895D0
      tb = 3D0
      hL = 0.1D0
      hK = 0.1D0
      hLs = 150D0
      hKs = 0D0
      mue = 210D0
      MG2 = 1000D0
      MG3 = 1000D0
      MG1 = 500D0
      Ml2 = 1000D0
      Ml3 = 1000D0
      Mr2 = 1000D0
      Mr3 = 1000D0
      Al = 1000D0
      Ae = 0D0
      Am = 0D0
      Mq2 = 1000D0
      Mq3 = 1000D0
      Mu2 = 1000D0
      Mu3 = 1000D0
      Md2 = 1000D0
      Md3 = 1000D0
      At = 1000D0
      Ab = 1000D0
      Au = 0D0
      Ad = 0D0
      As = 0D0
      Ac = 0D0
      dlh2 = 0D0
      wh1 = 0D0
      wh2 = 0D0
      wh3 = 0D0
      wha = 0D0
      whb = 0D0
      wHc = 0D0
      wC1 = 0D0
      wC2 = 0D0
      wNE2 = 0D0
      wNE3 = 0D0
      wNE4 = 0D0
      wSG = 0D0
      wNE5 = 0D0
      wSe1 = 0D0
      wSe2 = 0D0
      wSm1 = 0D0
      wSm2 = 0D0
      wSl1 = 0D0
      wSl2 = 0D0
      wSne = 0D0
      wSnm = 0D0
      wSnl = 0D0
      wSu1 = 0D0
      wSu2 = 0D0
      wSd1 = 0D0
      wSd2 = 0D0
      wSc1 = 0D0
      wSc2 = 0D0
      wSs1 = 0D0
      wSs2 = 0D0
      wSt1 = 0D0
      wSt2 = 0D0
      wSb1 = 0D0
      wSb2 = 0D0
	  
      if(slhaValExists1('MINPAR',3).GT.0) then
        tb=  slhaVal1('MINPAR',Q,3)
      else 
        tb= slhaVal1('HMIX',Q,2) 
      endif

      mue=slhaVal1('HMIX', Q,1)

      MG1=slhaVal1('MSOFT',Q,1)
      MG2=slhaVal1('MSOFT',Q,2)
      MG3=slhaVal1('MSOFT',Q,3)
      Ml2=slhaVal1('MSOFT',Q,32)
c	  Ml1=Ml2
      Ml3=slhaVal1('MSOFT',Q,33)  
      Mr2=slhaVal1('MSOFT',Q,35)
c	  Mr1=Mr2
      Mr3=slhaVal1('MSOFT',Q,36)
      Mq2=slhaVal1('MSOFT',Q,42) 
c	  Mq1=Mq2
      Mq3=slhaVal1('MSOFT',Q,43)
      Mu2=slhaVal1('MSOFT',Q,45)  
c	  Mu1=Mu2
      Mu3=slhaVal1('MSOFT',Q,46)      
      Md2=slhaVal1('MSOFT',Q,48)
c	  Md1=Md2
      Md3=slhaVal1('MSOFT',Q,49)

      Al= slhaVal2('Ae',   Q,3,3)
      Am=Al
      if(slhaValExists2('Ae',2,2).ne.0) Am=slhaVal2('Ae',Q,2,2)
      At= slhaVal2('Au',   Q,3,3)
      Ab= slhaVal2('Ad',   Q,3,3)

      MH3=slhaVal1('MASS', Q,36)


	  
c      tb = 5D0
      hL = slhaVal1('EXTPAR',Q,61)
      hK = slhaVal1('EXTPAR',Q,62)
      hLs = slhaVal1('EXTPAR',Q,63)
      hKs = slhaVal1('EXTPAR',Q,64)

      end



      subroutine ModelConstIni(*)
      implicit none

#include "model.h"

      double precision higs1,higs2,higp1,higp2,neu5,neutd

      call ModelDefaults
      CW = sqrt(1D0-SW**2)
      S2W = 2D0*SW*CW
      C2W = CW**2-SW**2
      MW = MZ*CW
      c12 = sqrt(1D0-s12**2)
      c23 = sqrt(1D0-s23**2)
      c13 = sqrt(1D0-s13**2)
      Vud = c12
      Vus = s12
      Vcs = Vud
      Vcd = -Vus
      sb = tb/sqrt(1D0+tb**2)
      cb = sqrt(1D0-sb**2)
      t2b = 2D0*tb/(1D0-tb**2)
      c2b = cb**2-sb**2
      s2b = 4D0*sb*cb/2D0
      xvev = mue/hL
      ls1 = Me*Ae/Sqrt2*EE/SW/MW/cb
      ls2 = Mm*Am/Sqrt2*EE/SW/MW/cb
      ls3 = Ml*Al/Sqrt2*EE/SW/MW/cb
      MSne = sqrt(c2b*MW**2/2D0+CW**2*Ml2**2)/CW
      MSnm = sqrt(c2b*MW**2/2D0+CW**2*Ml2**2)/CW
      MSnl = sqrt(c2b*MW**2/2D0+CW**2*Ml3**2)/CW
      MSeLL = Ml2**2+Me**2-MW**2/CW**2*(1D0/2D0-SW**2)*c2b
      MSeRR = Mr2**2+Me**2-MW**2/CW**2*SW**2*c2b
      MSeLR = -Me*mue*tb
      MSe1 = sqrt((MSeLL+MSeRR-sqrt((MSeLL-MSeRR)**2+4D0
     &      *MSeLR**2))/2D0)
      MSe2 = sqrt((MSeLL+MSeRR+sqrt((MSeLL-MSeRR)**2+4D0
     &      *MSeLR**2))/2D0)
      MSeth = atan2(-2D0*MSeLR, -MSeLL+MSeRR)/2D0
      Zl11 = cos(MSeth)
      Zl14 = -sin(MSeth)
      Zl41 = -Zl14
      Zl44 = Zl11
      MSmuLL = Ml2**2+Mm**2-MW**2/CW**2*(1D0/2D0-SW**2)*c2b
      MSmuLR = MW*SW*Sqrt2*cb*ls2/EE-Mm*mue*tb
      MSmuRR = Mr2**2+Mm**2-MW**2/CW**2*SW**2*c2b
      MSm1 = sqrt((MSmuLL+MSmuRR-sqrt((MSmuLL-MSmuRR)**2
     &      +4D0*MSmuLR**2))/2D0)
      MSm2 = sqrt((MSmuLL+MSmuRR+sqrt((MSmuLL-MSmuRR)**2
     &      +4D0*MSmuLR**2))/2D0)
      MSmuth = atan2(-2D0*MSmuLR, -MSmuLL+MSmuRR)/2D0
      Zl22 = cos(MSmuth)
      Zl25 = -sin(MSmuth)
      Zl52 = -Zl25
      Zl55 = Zl22
      MSlLL = Ml3**2+Ml**2-MW**2/CW**2*(1D0/2D0-SW**2)*c2b
      MSlLR = MW*SW*Sqrt2*cb*ls3/EE-Ml*mue*tb
      MSlRR = Mr3**2+Ml**2-MW**2/CW**2*SW**2*c2b
      MSl1 = sqrt((MSlLL+MSlRR-sqrt((MSlLL-MSlRR)**2+4D0
     &      *MSlLR**2))/2D0)
      MSl2 = sqrt((MSlLL+MSlRR+sqrt((MSlLL-MSlRR)**2+4D0
     &      *MSlLR**2))/2D0)
      MSlth = atan2(-2D0*MSlLR, -MSlLL+MSlRR)/2D0
      Zl33 = cos(MSlth)
      Zl36 = -sin(MSlth)
      Zl63 = -Zl36
      Zl66 = Zl33
      us1 = -Mu*Au*EE/SW/MW/Sqrt2/sb
      ds1 = Md*Ad*EE/SW/MW/Sqrt2/cb
      ds2 = Ms*As*EE/SW/MW/Sqrt2/cb
      us2 = -Mc*Ac*EE/SW/MW/Sqrt2/sb
      ds3 = Mb*Ab*EE/SW/MW/Sqrt2/cb
      us3 = -At*Mt*EE/SW/MW/Sqrt2/sb
      MSuLL = Mq2**2+Mu**2+MW**2/CW**2*(1D0/2D0-2D0/3D0*SW
     &      **2)*c2b
      MSuLR = -MW*SW*Sqrt2*sb*us1/EE-Mu*mue/tb
      MSuRR = Mu2**2+Mu**2+2D0/3D0*MW**2/CW**2*SW**2*c2b
      MSu1 = sqrt((MSuLL+MSuRR-sqrt((MSuLL-MSuRR)**2+4D0
     &      *MSuLR**2))/2D0)
      MSu2 = sqrt((MSuLL+MSuRR+sqrt((MSuLL-MSuRR)**2+4D0
     &      *MSuLR**2))/2D0)
      MSuth = atan2(-2D0*MSuLR, -MSuLL+MSuRR)/2D0
      Zu11 = cos(MSuth)
      Zu14 = -sin(MSuth)
      Zu41 = -Zu14
      Zu44 = Zu11
      MSdLL = Mq2**2+Md**2-MW**2/CW**2*(1D0/2D0-1D0/3D0*SW
     &      **2)*c2b
      MSdLR = MW*SW*Sqrt2*cb*ds1/EE-Md*mue*tb
      MSdRR = Md2**2+Md**2-1D0/3D0*MW**2/CW**2*SW**2*c2b
      MSd1 = sqrt((MSdLL+MSdRR-sqrt((MSdLL-MSdRR)**2+4D0
     &      *MSdLR**2))/2D0)
      MSd2 = sqrt((MSdLL+MSdRR+sqrt((MSdLL-MSdRR)**2+4D0
     &      *MSdLR**2))/2D0)
      MSdth = atan2(-2D0*MSdLR, -MSdLL+MSdRR)/2D0
      Zd11 = cos(MSdth)
      Zd14 = -sin(MSdth)
      Zd41 = -Zd14
      Zd44 = Zd11
      MScLL = Mq2**2+Mc**2+MW**2/CW**2*(1D0/2D0-2D0/3D0*SW
     &      **2)*c2b
      MScLR = -MW*SW*Sqrt2*sb*us2/EE-Mc*mue/tb
      MScRR = Mu2**2+Mc**2+2D0/3D0*MW**2/CW**2*SW**2*c2b
      MSc1 = sqrt((MScLL+MScRR-sqrt((MScLL-MScRR)**2+4D0
     &      *MScLR**2))/2D0)
      MSc2 = sqrt((MScLL+MScRR+sqrt((MScLL-MScRR)**2+4D0
     &      *MScLR**2))/2D0)
      MScth = atan2(-2D0*MScLR, -MScLL+MScRR)/2D0
      Zu22 = cos(MScth)
      Zu25 = -sin(MScth)
      Zu52 = -Zu25
      Zu55 = Zu22
      MSsLL = Mq2**2+Ms**2-MW**2/CW**2*(1D0/2D0-1D0/3D0*SW
     &      **2)*c2b
      MSsLR = MW*SW*Sqrt2*cb*ds2/EE-Ms*mue*tb
      MSsRR = Md2**2+Ms**2-1D0/3D0*MW**2/CW**2*SW**2*c2b
      MSs1 = sqrt((MSsLL+MSsRR-sqrt((MSsLL-MSsRR)**2+4D0
     &      *MSsLR**2))/2D0)
      MSs2 = sqrt((MSsLL+MSsRR+sqrt((MSsLL-MSsRR)**2+4D0
     &      *MSsLR**2))/2D0)
      MSsth = atan2(-2D0*MSsLR, -MSsLL+MSsRR)/2D0
      Zd22 = cos(MSsth)
      Zd25 = -sin(MSsth)
      Zd52 = -Zd25
      Zd55 = Zd22
      MStLL = Mq3**2+Mt**2+MW**2/CW**2*(1D0/2D0-2D0/3D0*SW
     &      **2)*c2b
      MStLR = -MW*SW*Sqrt2*sb*us3/EE-Mt*mue/tb
      MStRR = Mu3**2+Mt**2+2D0/3D0*MW**2/CW**2*SW**2*c2b
      MSt1 = sqrt((MStLL+MStRR-sqrt((MStLL-MStRR)**2+4D0
     &      *MStLR**2))/2D0)
      MSt2 = sqrt((MStLL+MStRR+sqrt((MStLL-MStRR)**2+4D0
     &      *MStLR**2))/2D0)
      MStth = atan2(-2D0*MStLR, -MStLL+MStRR)/2D0
      Zu33 = cos(MStth)
      Zu36 = -sin(MStth)
      Zu63 = -Zu36
      Zu66 = Zu33
      MSbLL = Mq3**2+Mb**2-MW**2/CW**2*(1D0/2D0-1D0/3D0*SW
     &      **2)*c2b
      MSbLR = MW*SW*Sqrt2*cb*ds3/EE-Mb*mue*tb
      MSbRR = Md3**2+Mb**2-1D0/3D0*MW**2/CW**2*SW**2*c2b
      MSb1 = sqrt((MSbLL+MSbRR-sqrt((MSbLL-MSbRR)**2+4D0
     &      *MSbLR**2))/2D0)
      MSb2 = sqrt((MSbLL+MSbRR+sqrt((MSbLL-MSbRR)**2+4D0
     &      *MSbLR**2))/2D0)
      MSbth = atan2(-2D0*MSbLR, -MSbLL+MSbRR)/2D0
      Zd33 = cos(MSbth)
      Zd36 = -sin(MSbth)
      Zd63 = -Zd36
      Zd66 = Zd33
      ntkhs = higs1(SW, MW, EE, sb, cb, hK, hL, xvev, hKs
     &      , hLs, dlh2)
      Mh1 = higs2(0D0, 1D0, ntkhs)
      Mh2 = higs2(0D0, 2D0, ntkhs)
      Mh3 = higs2(0D0, 3D0, ntkhs)
      Zh11 = higs2(1D0, 1D0, ntkhs)
      Zh12 = higs2(1D0, 2D0, ntkhs)
      Zh13 = higs2(1D0, 3D0, ntkhs)
      Zh21 = higs2(2D0, 1D0, ntkhs)
      Zh22 = higs2(2D0, 2D0, ntkhs)
      Zh23 = higs2(2D0, 3D0, ntkhs)
      Zh31 = higs2(3D0, 1D0, ntkhs)
      Zh32 = higs2(3D0, 2D0, ntkhs)
      Zh33 = higs2(3D0, 3D0, ntkhs)
      ntkhp = higp1(SW, MW, EE, sb, cb, hK, hL, xvev, hKs
     &      , hLs)
      Mha = higp2(0D0, 1D0, ntkhs)
      Mhb = higp2(0D0, 2D0, ntkhp)
      Za11 = -higp2(1D0, 1D0, ntkhp)
      Za12 = -higp2(1D0, 2D0, ntkhp)
      Za13 = -higp2(1D0, 3D0, ntkhp)
      Za21 = -higp2(2D0, 1D0, ntkhp)
      Za22 = -higp2(2D0, 2D0, ntkhp)
      Za23 = -higp2(2D0, 3D0, ntkhp)
      Za31 = -higp2(3D0, 1D0, ntkhp)
      Za32 = -higp2(3D0, 2D0, ntkhp)
      Za33 = -higp2(3D0, 3D0, ntkhp)
      MHc = sqrt(MW**2*(1D0-2D0/EE**2*SW**2*hL**2)+2D0*hL
     &      *xvev/s2b*(hLs+hK*xvev))
      Zcsx = (MG2**2+mue**2)/2D0+MW**2
      Zctx = (MG2**2-mue**2)**2/4D0+MW**4*c2b**2+MW**2*(MG2
     &      **2+mue**2+2D0*mue*MG2*s2b)
      MC01 = sqrt(Zcsx-sqrt(Zctx))
      MC02 = sqrt(Zcsx+sqrt(Zctx))
      Zcc2u = -(MG2**2-mue**2-2D0*MW**2*c2b)/(MC02**2-MC01
     &      **2)
      Zcc2v = -(MG2**2-mue**2+2D0*MW**2*c2b)/(MC02**2-MC01
     &      **2)
      Zcsigu = -(MG2*cb+mue*sb)/sqrt((MG2*cb+mue*sb)**2)
      Zcsigv = -(MG2*sb+mue*cb)/sqrt((MG2*sb+mue*cb)**2)
      Zccu = sqrt((1D0+Zcc2u)/2D0)
      Zcsu = sqrt((1D0-Zcc2u)/2D0)*Zcsigu
      Zccv = sqrt((1D0+Zcc2v)/2D0)
      Zcsv = sqrt((1D0-Zcc2v)/2D0)*Zcsigv
      Zcsig1 = (mue*MG2-MW**2*s2b)/sqrt((mue*MG2-MW**2*s2b
     &      )**2)
      Zcsig2 = MW*Sqrt2*(cb*Zcsu*Zccv+sb*Zccu*Zcsv)+mue*Zcsv
     &      *Zcsu+MG2*Zccu*Zccv
      Zcsig3 = Zcsig1*Zcsig2/sqrt(Zcsig2**2)
      Zcsig4 = -MW*Sqrt2*(cb*Zccu*Zcsv+sb*Zcsu*Zccv)+mue
     &      *Zccv*Zccu+MG2*Zcsu*Zcsv
      MC1 = Zcsig2**2/sqrt(Zcsig2**2)
      MC2 = Zcsig4**2/sqrt(Zcsig4**2)
      Zm11 = Zccu*Zcsig1*Zcsig3
      Zm12 = -Zcsu*Zcsig3
      Zm21 = Zcsu*Zcsig1*Zcsig3
      Zm22 = Zccu*Zcsig3
      Zp11 = Zccv
      Zp12 = -Zcsv
      Zp21 = Zcsv
      Zp22 = Zccv
      neutk = neu5(MG1, MG2, SW, MW, EE, sb, cb, hK, hL, xvev
     &      )
      Zn11r = neutd(neutk, 1D0, 1D0, 1D0)
      Zn11i = -neutd(neutk, 2D0, 1D0, 1D0)
      Zn12r = neutd(neutk, 1D0, 1D0, 2D0)
      Zn12i = -neutd(neutk, 2D0, 1D0, 2D0)
      Zn13r = neutd(neutk, 1D0, 1D0, 3D0)
      Zn13i = -neutd(neutk, 2D0, 1D0, 3D0)
      Zn14r = neutd(neutk, 1D0, 1D0, 4D0)
      Zn14i = -neutd(neutk, 2D0, 1D0, 4D0)
      Zn15r = neutd(neutk, 1D0, 1D0, 5D0)
      Zn15i = -neutd(neutk, 2D0, 1D0, 5D0)
      Zn21r = neutd(neutk, 1D0, 2D0, 1D0)
      Zn21i = -neutd(neutk, 2D0, 2D0, 1D0)
      Zn22r = neutd(neutk, 1D0, 2D0, 2D0)
      Zn22i = -neutd(neutk, 2D0, 2D0, 2D0)
      Zn23r = neutd(neutk, 1D0, 2D0, 3D0)
      Zn23i = -neutd(neutk, 2D0, 2D0, 3D0)
      Zn24r = neutd(neutk, 1D0, 2D0, 4D0)
      Zn24i = -neutd(neutk, 2D0, 2D0, 4D0)
      Zn25r = neutd(neutk, 1D0, 2D0, 5D0)
      Zn25i = -neutd(neutk, 2D0, 2D0, 5D0)
      Zn31r = neutd(neutk, 1D0, 3D0, 1D0)
      Zn31i = -neutd(neutk, 2D0, 3D0, 1D0)
      Zn32r = neutd(neutk, 1D0, 3D0, 2D0)
      Zn32i = -neutd(neutk, 2D0, 3D0, 2D0)
      Zn33r = neutd(neutk, 1D0, 3D0, 3D0)
      Zn33i = -neutd(neutk, 2D0, 3D0, 3D0)
      Zn34r = neutd(neutk, 1D0, 3D0, 4D0)
      Zn34i = -neutd(neutk, 2D0, 3D0, 4D0)
      Zn35r = neutd(neutk, 1D0, 3D0, 5D0)
      Zn35i = -neutd(neutk, 2D0, 3D0, 5D0)
      Zn41r = neutd(neutk, 1D0, 4D0, 1D0)
      Zn41i = -neutd(neutk, 2D0, 4D0, 1D0)
      Zn42r = neutd(neutk, 1D0, 4D0, 2D0)
      Zn42i = -neutd(neutk, 2D0, 4D0, 2D0)
      Zn43r = neutd(neutk, 1D0, 4D0, 3D0)
      Zn43i = -neutd(neutk, 2D0, 4D0, 3D0)
      Zn44r = neutd(neutk, 1D0, 4D0, 4D0)
      Zn44i = -neutd(neutk, 2D0, 4D0, 4D0)
      Zn45r = neutd(neutk, 1D0, 4D0, 5D0)
      Zn45i = -neutd(neutk, 2D0, 4D0, 5D0)
      Zn51r = neutd(neutk, 1D0, 5D0, 1D0)
      Zn51i = -neutd(neutk, 2D0, 5D0, 1D0)
      Zn52r = neutd(neutk, 1D0, 5D0, 2D0)
      Zn52i = -neutd(neutk, 2D0, 5D0, 2D0)
      Zn53r = neutd(neutk, 1D0, 5D0, 3D0)
      Zn53i = -neutd(neutk, 2D0, 5D0, 3D0)
      Zn54r = neutd(neutk, 1D0, 5D0, 4D0)
      Zn54i = -neutd(neutk, 2D0, 5D0, 4D0)
      Zn55r = neutd(neutk, 1D0, 5D0, 5D0)
      Zn55i = -neutd(neutk, 2D0, 5D0, 5D0)
      MNE1 = neutd(neutk, 3D0, 1D0, 1D0)
      MNE2 = neutd(neutk, 3D0, 2D0, 2D0)
      MNE3 = neutd(neutk, 3D0, 3D0, 3D0)
      MNE4 = neutd(neutk, 3D0, 4D0, 4D0)
      MNE5 = neutd(neutk, 3D0, 5D0, 5D0)
      Zn11 = Zn11r+cI*Zn11i
      Zn12 = Zn12r+cI*Zn12i
      Zn13 = Zn13r+cI*Zn13i
      Zn14 = Zn14r+cI*Zn14i
      Zn15 = Zn15r+cI*Zn15i
      Zn21 = Zn21r+cI*Zn21i
      Zn22 = Zn22r+cI*Zn22i
      Zn23 = Zn23r+cI*Zn23i
      Zn24 = Zn24r+cI*Zn24i
      Zn25 = Zn25r+cI*Zn25i
      Zn31 = Zn31r+cI*Zn31i
      Zn32 = Zn32r+cI*Zn32i
      Zn33 = Zn33r+cI*Zn33i
      Zn34 = Zn34r+cI*Zn34i
      Zn35 = Zn35r+cI*Zn35i
      Zn41 = Zn41r+cI*Zn41i
      Zn42 = Zn42r+cI*Zn42i
      Zn43 = Zn43r+cI*Zn43i
      Zn44 = Zn44r+cI*Zn44i
      Zn45 = Zn45r+cI*Zn45i
      Zn51 = Zn51r+cI*Zn51i
      Zn52 = Zn52r+cI*Zn52i
      Zn53 = Zn53r+cI*Zn53i
      Zn54 = Zn54r+cI*Zn54i
      Zn55 = Zn55r+cI*Zn55i
      MSG = MG3
      Mesq = Me**2
      Mmsq = Mm**2
      Mlsq = Ml**2
      Musq = Mu**2
      Mdsq = Md**2
      Mcsq = Mc**2
      Mssq = Ms**2
      Mbsq = Mb**2
      Mtsq = Mt**2
      MSe1sq = MSe1**2
      MSm1sq = MSm1**2
      MSl1sq = MSl1**2
      MSu1sq = MSu1**2
      MSd1sq = MSd1**2
      MSc1sq = MSc1**2
      MSs1sq = MSs1**2
      MSb1sq = MSb1**2
      MSt1sq = MSt1**2
      MSe2sq = MSe2**2
      MSm2sq = MSm2**2
      MSl2sq = MSl2**2
      MSu2sq = MSu2**2
      MSd2sq = MSd2**2
      MSc2sq = MSc2**2
      MSs2sq = MSs2**2
      MSb2sq = MSb2**2
      MSt2sq = MSt2**2
      MSnesq = MSne**2
      MSnmsq = MSnm**2
      MSnlsq = MSnl**2

      call aaini01
      call aaini02
      call mtrini
      end

      subroutine aaini01
      implicit none
#include "model.h"

      AAABR(1) = MW*(CW*Za13*cb-CW*Za23*sb+SW**2/CW*Za13
     &      *cb-SW**2/CW*Za23*sb)
      AAABR(2) = MW/CW
      AAABR(3) = MW*(CW*Za11*cb-CW*Za21*sb+SW**2/CW*Za11
     &      *cb-SW**2/CW*Za21*sb)
      AAABR(4) = MW*(CW*Za12*cb-CW*Za22*sb+SW**2/CW*Za12
     &      *cb-SW**2/CW*Za22*sb)
      AAABR(5) = CW**2+SW**2
      AAABR(6) = MW**2*(CW**2+2D0*SW**2+SW**4/CW**2)
      AAABR(7) = EE*MW
      AAABR(8) = EE*MW/SW
      AAABR(9) = EE*MW/SW*(Zh12*cb+Zh22*sb)
      AAABR(10) = EE*MW/SW*(Zh13*cb+Zh23*sb)
      AAABR(11) = EE*MW/SW*(Zh11*cb+Zh21*sb)
      AAABR(12) = EE*MW*(SW/CW-CW/SW)
      AAABR(13) = EE/CW*MW/SW
      AAABR(14) = EE/CW**2*MW/SW*(Zh11*cb+Zh21*sb)
      AAABR(15) = EE/CW**2*MW/SW*(Zh13*cb+Zh23*sb)
      AAABR(16) = EE/CW**2*MW/SW*(Zh12*cb+Zh22*sb)
      AAABR(17) = CW*EE/SW
      AAABR(18) = 2D0*Sqrt2*Zh32*hL**2*xvev+EE*MW/SW*Zh12
     &      *cb+3D0*EE*MW/SW*Zh22*sb-2D0*EE*MW/SW*Zh22*sb
     &      **3-EE/CW**2*MW*SW*Zh12*cb+EE/CW**2*MW*SW*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zh22*sb**3-4D0*MW/EE*SW
     &      *Zh22*hL**2*sb+4D0*MW/EE*SW*Zh22*hL**2*sb**3
     &      +2D0*EE*MW/SW*Zh12*cb*sb**2+2D0*Sqrt2*Zh32*cb
     &      *hL*hLs*sb+2D0*EE/CW**2*MW*SW*Zh12*cb*sb**2-4D0
     &      *MW/EE*SW*Zh12*cb*hL**2*sb**2+4D0*Sqrt2*Zh32
     &      *cb*hK*hL*sb*xvev
      AAABR(19) = 2D0*Sqrt2*Zh33*hL**2*xvev+EE*MW/SW*Zh13
     &      *cb+3D0*EE*MW/SW*Zh23*sb-2D0*EE*MW/SW*Zh23*sb
     &      **3-EE/CW**2*MW*SW*Zh13*cb+EE/CW**2*MW*SW*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zh23*sb**3-4D0*MW/EE*SW
     &      *Zh23*hL**2*sb+4D0*MW/EE*SW*Zh23*hL**2*sb**3
     &      +2D0*EE*MW/SW*Zh13*cb*sb**2+2D0*Sqrt2*Zh33*cb
     &      *hL*hLs*sb+2D0*EE/CW**2*MW*SW*Zh13*cb*sb**2-4D0
     &      *MW/EE*SW*Zh13*cb*hL**2*sb**2+4D0*Sqrt2*Zh33
     &      *cb*hK*hL*sb*xvev
      AAABR(20) = 2D0*Sqrt2*Zh31*hL**2*xvev+EE*MW/SW*Zh11
     &      *cb+3D0*EE*MW/SW*Zh21*sb-2D0*EE*MW/SW*Zh21*sb
     &      **3-EE/CW**2*MW*SW*Zh11*cb+EE/CW**2*MW*SW*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zh21*sb**3-4D0*MW/EE*SW
     &      *Zh21*hL**2*sb+4D0*MW/EE*SW*Zh21*hL**2*sb**3
     &      +2D0*EE*MW/SW*Zh11*cb*sb**2+2D0*Sqrt2*Zh31*cb
     &      *hL*hLs*sb+2D0*EE/CW**2*MW*SW*Zh11*cb*sb**2-4D0
     &      *MW/EE*SW*Zh11*cb*hL**2*sb**2+4D0*Sqrt2*Zh31
     &      *cb*hK*hL*sb*xvev
      AAABR(21) = 2D0*Zl41*ls1*sb-2D0*EE*MW/SW*Sqrt2*Zl11
     &      /cb*sb+2D0*EE*MW/SW*Sqrt2*Zl11/cb*sb**3+EE/MW
     &      *Me/SW*Sqrt2*Zl41*hL*xvev+EE/MW*Me**2/SW*Sqrt2
     &      *Zl11/cb*sb
      AAABR(22) = 2D0*Zl52*ls2*sb-2D0*EE*MW/SW*Sqrt2*Zl22
     &      /cb*sb+2D0*EE*MW/SW*Sqrt2*Zl22/cb*sb**3+EE/MW
     &      *Mm/SW*Sqrt2*Zl52*hL*xvev+EE/MW*Mm**2/SW*Sqrt2
     &      *Zl22/cb*sb
      AAABR(23) = 2D0*Zl63*ls3*sb-2D0*EE*MW/SW*Sqrt2*Zl33
     &      /cb*sb+2D0*EE*MW/SW*Sqrt2*Zl33/cb*sb**3+EE/MW
     &      *Ml/SW*Sqrt2*Zl63*hL*xvev+EE/MW*Ml**2/SW*Sqrt2
     &      *Zl33/cb*sb
      AAABR(24) = 2D0*Zl55*ls2*sb-2D0*EE*MW/SW*Sqrt2*Zl25
     &      /cb*sb+2D0*EE*MW/SW*Sqrt2*Zl25/cb*sb**3+EE/MW
     &      *Mm/SW*Sqrt2*Zl55*hL*xvev+EE/MW*Mm**2/SW*Sqrt2
     &      *Zl25/cb*sb
      AAABR(25) = 2D0*Zl66*ls3*sb-2D0*EE*MW/SW*Sqrt2*Zl36
     &      /cb*sb+2D0*EE*MW/SW*Sqrt2*Zl36/cb*sb**3+EE/MW
     &      *Ml/SW*Sqrt2*Zl66*hL*xvev+EE/MW*Ml**2/SW*Sqrt2
     &      *Zl36/cb*sb
      AAABR(26) = 2D0*Zl44*ls1*sb-2D0*EE*MW/SW*Sqrt2*Zl14
     &      /cb*sb+2D0*EE*MW/SW*Sqrt2*Zl14/cb*sb**3+EE/MW
     &      *Me/SW*Sqrt2*Zl44*hL*xvev+EE/MW*Me**2/SW*Sqrt2
     &      *Zl14/cb*sb
      AAABR(27) = Vud*(Zd11*Zu41*cb*us1-Zd41*Zu11*ds1*sb
     &      +EE*MW/SW*Sqrt2*Zd11*Zu11*cb*sb)
      AAABR(28) = Vcd*(2D0*Zd11*Zu52*cb*us2-2D0*Zd41*Zu22
     &      *ds1*sb+2D0*EE*MW/SW*Sqrt2*Zd11*Zu22*cb*sb-EE
     &      /MW*Mc/SW*Sqrt2*Zd11*Zu52*hL*xvev-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd11*Zu22*cb/sb)
      AAABR(29) = Vus*(2D0*Zd22*Zu41/cb*us1-2D0*Zd52*Zu11
     &      *ds2*sb-2D0*Zd22*Zu41/cb*sb**2*us1+2D0*EE*MW
     &      /SW*Sqrt2*Zd22*Zu11/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu11/cb*sb**3-EE/MW*Ms/SW*Sqrt2*Zd52*Zu11
     &      *hL*xvev-EE/MW*Ms**2/SW*Sqrt2*Zd22*Zu11/cb*sb)
      AAABR(30) = 2D0*Zd33*Zu63/cb*us3-2D0*Zd63*Zu33*ds3
     &      *sb-2D0*Zd33*Zu63/cb*sb**2*us3+2D0*EE*MW/SW*Sqrt2
     &      *Zd33*Zu33/cb*sb-2D0*EE*MW/SW*Sqrt2*Zd33*Zu33
     &      /cb*sb**3-EE/MW*Mb/SW*Sqrt2*Zd63*Zu33*hL*xvev
     &      -EE/MW*Mb**2/SW*Sqrt2*Zd33*Zu33/cb*sb-EE/MW*Mt
     &      /SW*Sqrt2*Zd33*Zu63*hL*xvev-EE/MW*Mt**2/SW*Sqrt2
     &      *Zd33*Zu33/cb/sb+EE/MW*Mt**2/SW*Sqrt2*Zd33*Zu33
     &      /cb*sb-EE/MW*Mb*Mt/SW*Sqrt2*Zd63*Zu63/cb/sb
      AAABR(31) = Vcs*(2D0*Zd22*Zu52/cb*us2-2D0*Zd52*Zu22
     &      *ds2*sb-2D0*Zd22*Zu52/cb*sb**2*us2+2D0*EE*MW
     &      /SW*Sqrt2*Zd22*Zu22/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu22/cb*sb**3-EE/MW*Mc/SW*Sqrt2*Zd22*Zu52
     &      *hL*xvev-EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu22/cb/sb
     &      +EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu22/cb*sb-EE/MW*Ms
     &      /SW*Sqrt2*Zd52*Zu22*hL*xvev-EE/MW*Ms**2/SW*Sqrt2
     &      *Zd22*Zu22/cb*sb-EE/MW*Mc*Ms/SW*Sqrt2*Zd52*Zu52
     &      /cb/sb)
      AAABR(32) = 2D0*Zd36*Zu63/cb*us3-2D0*Zd66*Zu33*ds3
     &      *sb-2D0*Zd36*Zu63/cb*sb**2*us3+2D0*EE*MW/SW*Sqrt2
     &      *Zd36*Zu33/cb*sb-2D0*EE*MW/SW*Sqrt2*Zd36*Zu33
     &      /cb*sb**3-EE/MW*Mb/SW*Sqrt2*Zd66*Zu33*hL*xvev
     &      -EE/MW*Mb**2/SW*Sqrt2*Zd36*Zu33/cb*sb-EE/MW*Mt
     &      /SW*Sqrt2*Zd36*Zu63*hL*xvev-EE/MW*Mt**2/SW*Sqrt2
     &      *Zd36*Zu33/cb/sb+EE/MW*Mt**2/SW*Sqrt2*Zd36*Zu33
     &      /cb*sb-EE/MW*Mb*Mt/SW*Sqrt2*Zd66*Zu63/cb/sb
      AAABR(33) = Vcd*(2D0*Zd14*Zu52*cb*us2-2D0*Zd44*Zu22
     &      *ds1*sb+2D0*EE*MW/SW*Sqrt2*Zd14*Zu22*cb*sb-EE
     &      /MW*Mc/SW*Sqrt2*Zd14*Zu52*hL*xvev-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd14*Zu22*cb/sb)
      AAABR(34) = Vcs*(2D0*Zd25*Zu52/cb*us2-2D0*Zd55*Zu22
     &      *ds2*sb-2D0*Zd25*Zu52/cb*sb**2*us2+2D0*EE*MW
     &      /SW*Sqrt2*Zd25*Zu22/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu22/cb*sb**3-EE/MW*Mc/SW*Sqrt2*Zd25*Zu52
     &      *hL*xvev-EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu22/cb/sb
     &      +EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu22/cb*sb-EE/MW*Ms
     &      /SW*Sqrt2*Zd55*Zu22*hL*xvev-EE/MW*Ms**2/SW*Sqrt2
     &      *Zd25*Zu22/cb*sb-EE/MW*Mc*Ms/SW*Sqrt2*Zd55*Zu52
     &      /cb/sb)
      AAABR(35) = Vud*(Zd14*Zu41*cb*us1-Zd44*Zu11*ds1*sb
     &      +EE*MW/SW*Sqrt2*Zd14*Zu11*cb*sb)
      AAABR(36) = Vus*(2D0*Zd25*Zu41/cb*us1-2D0*Zd55*Zu11
     &      *ds2*sb-2D0*Zd25*Zu41/cb*sb**2*us1+2D0*EE*MW
     &      /SW*Sqrt2*Zd25*Zu11/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu11/cb*sb**3-EE/MW*Ms/SW*Sqrt2*Zd55*Zu11
     &      *hL*xvev-EE/MW*Ms**2/SW*Sqrt2*Zd25*Zu11/cb*sb)
      AAABR(37) = 2D0*Zd33*Zu66/cb*us3-2D0*Zd63*Zu36*ds3
     &      *sb-2D0*Zd33*Zu66/cb*sb**2*us3+2D0*EE*MW/SW*Sqrt2
     &      *Zd33*Zu36/cb*sb-2D0*EE*MW/SW*Sqrt2*Zd33*Zu36
     &      /cb*sb**3-EE/MW*Mb/SW*Sqrt2*Zd63*Zu36*hL*xvev
     &      -EE/MW*Mb**2/SW*Sqrt2*Zd33*Zu36/cb*sb-EE/MW*Mt
     &      /SW*Sqrt2*Zd33*Zu66*hL*xvev-EE/MW*Mt**2/SW*Sqrt2
     &      *Zd33*Zu36/cb/sb+EE/MW*Mt**2/SW*Sqrt2*Zd33*Zu36
     &      /cb*sb-EE/MW*Mb*Mt/SW*Sqrt2*Zd63*Zu66/cb/sb
      AAABR(38) = Vcd*(2D0*Zd11*Zu55*cb*us2-2D0*Zd41*Zu25
     &      *ds1*sb+2D0*EE*MW/SW*Sqrt2*Zd11*Zu25*cb*sb-EE
     &      /MW*Mc/SW*Sqrt2*Zd11*Zu55*hL*xvev-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd11*Zu25*cb/sb)
      AAABR(39) = Vcs*(2D0*Zd22*Zu55/cb*us2-2D0*Zd52*Zu25
     &      *ds2*sb-2D0*Zd22*Zu55/cb*sb**2*us2+2D0*EE*MW
     &      /SW*Sqrt2*Zd22*Zu25/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu25/cb*sb**3-EE/MW*Mc/SW*Sqrt2*Zd22*Zu55
     &      *hL*xvev-EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu25/cb/sb
     &      +EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu25/cb*sb-EE/MW*Ms
     &      /SW*Sqrt2*Zd52*Zu25*hL*xvev-EE/MW*Ms**2/SW*Sqrt2
     &      *Zd22*Zu25/cb*sb-EE/MW*Mc*Ms/SW*Sqrt2*Zd52*Zu55
     &      /cb/sb)
      AAABR(40) = Vud*(Zd11*Zu44*cb*us1-Zd41*Zu14*ds1*sb
     &      +EE*MW/SW*Sqrt2*Zd11*Zu14*cb*sb)
      AAABR(41) = Vus*(2D0*Zd22*Zu44/cb*us1-2D0*Zd52*Zu14
     &      *ds2*sb-2D0*Zd22*Zu44/cb*sb**2*us1+2D0*EE*MW
     &      /SW*Sqrt2*Zd22*Zu14/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu14/cb*sb**3-EE/MW*Ms/SW*Sqrt2*Zd52*Zu14
     &      *hL*xvev-EE/MW*Ms**2/SW*Sqrt2*Zd22*Zu14/cb*sb)
      AAABR(42) = 2D0*Zd36*Zu66/cb*us3-2D0*Zd66*Zu36*ds3
     &      *sb-2D0*Zd36*Zu66/cb*sb**2*us3+2D0*EE*MW/SW*Sqrt2
     &      *Zd36*Zu36/cb*sb-2D0*EE*MW/SW*Sqrt2*Zd36*Zu36
     &      /cb*sb**3-EE/MW*Mb/SW*Sqrt2*Zd66*Zu36*hL*xvev
     &      -EE/MW*Mb**2/SW*Sqrt2*Zd36*Zu36/cb*sb-EE/MW*Mt
     &      /SW*Sqrt2*Zd36*Zu66*hL*xvev-EE/MW*Mt**2/SW*Sqrt2
     &      *Zd36*Zu36/cb/sb+EE/MW*Mt**2/SW*Sqrt2*Zd36*Zu36
     &      /cb*sb-EE/MW*Mb*Mt/SW*Sqrt2*Zd66*Zu66/cb/sb
      AAABR(43) = Vud*(Zd14*Zu44*cb*us1-Zd44*Zu14*ds1*sb
     &      +EE*MW/SW*Sqrt2*Zd14*Zu14*cb*sb)
      AAABR(44) = Vus*(2D0*Zd25*Zu44/cb*us1-2D0*Zd55*Zu14
     &      *ds2*sb-2D0*Zd25*Zu44/cb*sb**2*us1+2D0*EE*MW
     &      /SW*Sqrt2*Zd25*Zu14/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu14/cb*sb**3-EE/MW*Ms/SW*Sqrt2*Zd55*Zu14
     &      *hL*xvev-EE/MW*Ms**2/SW*Sqrt2*Zd25*Zu14/cb*sb)
      AAABR(45) = Vcd*(2D0*Zd14*Zu55*cb*us2-2D0*Zd44*Zu25
     &      *ds1*sb+2D0*EE*MW/SW*Sqrt2*Zd14*Zu25*cb*sb-EE
     &      /MW*Mc/SW*Sqrt2*Zd14*Zu55*hL*xvev-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd14*Zu25*cb/sb)
      AAABR(46) = Vcs*(2D0*Zd25*Zu55/cb*us2-2D0*Zd55*Zu25
     &      *ds2*sb-2D0*Zd25*Zu55/cb*sb**2*us2+2D0*EE*MW
     &      /SW*Sqrt2*Zd25*Zu25/cb*sb-2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu25/cb*sb**3-EE/MW*Mc/SW*Sqrt2*Zd25*Zu55
     &      *hL*xvev-EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu25/cb/sb
     &      +EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu25/cb*sb-EE/MW*Ms
     &      /SW*Sqrt2*Zd55*Zu25*hL*xvev-EE/MW*Ms**2/SW*Sqrt2
     &      *Zd25*Zu25/cb*sb-EE/MW*Mc*Ms/SW*Sqrt2*Zd55*Zu55
     &      /cb/sb)
      AAABR(47) = Sqrt2*Za33*hL*hLs-EE*MW/SW*Za13*sb-EE*MW
     &      /SW*Za23*cb-2D0*Sqrt2*Za33*hK*hL*xvev+2D0*MW
     &      /EE*SW*Za13*hL**2*sb+2D0*MW/EE*SW*Za23*cb*hL**2
      AAABR(48) = Sqrt2*Za31*hL*hLs-EE*MW/SW*Za11*sb-EE*MW
     &      /SW*Za21*cb-2D0*Sqrt2*Za31*hK*hL*xvev+2D0*MW
     &      /EE*SW*Za11*hL**2*sb+2D0*MW/EE*SW*Za21*cb*hL**2
      AAABR(49) = Sqrt2*Za32*hL*hLs-EE*MW/SW*Za12*sb-EE*MW
     &      /SW*Za22*cb-2D0*Sqrt2*Za32*hK*hL*xvev+2D0*MW
     &      /EE*SW*Za12*hL**2*sb+2D0*MW/EE*SW*Za22*cb*hL**2
      AAABR(50) = Sqrt2*Zh31*hL*hLs+EE*MW/SW*Zh11*sb-2D0
     &      *EE*MW/SW*Zh11*sb**3+EE*MW/SW*Zh21*cb+2D0*Sqrt2
     &      *Zh31*hK*hL*xvev-2D0*Sqrt2*Zh31*hL*hLs*sb**2
     &      +2D0*EE/CW**2*MW*SW*Zh11*sb-2D0*EE/CW**2*MW*SW
     &      *Zh11*sb**3-2D0*MW/EE*SW*Zh11*hL**2*sb+4D0*MW
     &      /EE*SW*Zh11*hL**2*sb**3-2D0*MW/EE*SW*Zh21*cb
     &      *hL**2-2D0*EE*MW/SW*Zh21*cb*sb**2-4D0*Sqrt2*Zh31
     &      *hK*hL*sb**2*xvev-2D0*EE/CW**2*MW*SW*Zh21*cb
     &      *sb**2+4D0*MW/EE*SW*Zh21*cb*hL**2*sb**2
      AAABR(51) = Sqrt2*Zh32*hL*hLs+EE*MW/SW*Zh12*sb-2D0
     &      *EE*MW/SW*Zh12*sb**3+EE*MW/SW*Zh22*cb+2D0*Sqrt2
     &      *Zh32*hK*hL*xvev-2D0*Sqrt2*Zh32*hL*hLs*sb**2
     &      +2D0*EE/CW**2*MW*SW*Zh12*sb-2D0*EE/CW**2*MW*SW
     &      *Zh12*sb**3-2D0*MW/EE*SW*Zh12*hL**2*sb+4D0*MW
     &      /EE*SW*Zh12*hL**2*sb**3-2D0*MW/EE*SW*Zh22*cb
     &      *hL**2-2D0*EE*MW/SW*Zh22*cb*sb**2-4D0*Sqrt2*Zh32
     &      *hK*hL*sb**2*xvev-2D0*EE/CW**2*MW*SW*Zh22*cb
     &      *sb**2+4D0*MW/EE*SW*Zh22*cb*hL**2*sb**2
      AAABR(52) = Sqrt2*Zh33*hL*hLs+EE*MW/SW*Zh13*sb-2D0
     &      *EE*MW/SW*Zh13*sb**3+EE*MW/SW*Zh23*cb+2D0*Sqrt2
     &      *Zh33*hK*hL*xvev-2D0*Sqrt2*Zh33*hL*hLs*sb**2
     &      +2D0*EE/CW**2*MW*SW*Zh13*sb-2D0*EE/CW**2*MW*SW
     &      *Zh13*sb**3-2D0*MW/EE*SW*Zh13*hL**2*sb+4D0*MW
     &      /EE*SW*Zh13*hL**2*sb**3-2D0*MW/EE*SW*Zh23*cb
     &      *hL**2-2D0*EE*MW/SW*Zh23*cb*sb**2-4D0*Sqrt2*Zh33
     &      *hK*hL*sb**2*xvev-2D0*EE/CW**2*MW*SW*Zh23*cb
     &      *sb**2+4D0*MW/EE*SW*Zh23*cb*hL**2*sb**2
      AAABR(53) = 2D0*Zl52/cb*ls2-2D0*Zl52/cb*ls2*sb**2-EE
     &      *MW/SW*Sqrt2*Zl22+EE/MW*Mm**2/SW*Sqrt2*Zl22+2D0
     &      *EE*MW/SW*Sqrt2*Zl22*sb**2-EE/MW*Mm/SW*Sqrt2
     &      *Zl52/cb*hL*sb*xvev
      AAABR(54) = 2D0*Zl41/cb*ls1-2D0*Zl41/cb*ls1*sb**2-EE
     &      *MW/SW*Sqrt2*Zl11+EE/MW*Me**2/SW*Sqrt2*Zl11+2D0
     &      *EE*MW/SW*Sqrt2*Zl11*sb**2-EE/MW*Me/SW*Sqrt2
     &      *Zl41/cb*hL*sb*xvev
      AAABR(55) = 2D0*Zl63/cb*ls3-2D0*Zl63/cb*ls3*sb**2-EE
     &      *MW/SW*Sqrt2*Zl33+EE/MW*Ml**2/SW*Sqrt2*Zl33+2D0
     &      *EE*MW/SW*Sqrt2*Zl33*sb**2-EE/MW*Ml/SW*Sqrt2
     &      *Zl63/cb*hL*sb*xvev
      AAABR(56) = Sqrt2*Za13*Zl22*Zl55*ls2-Sqrt2*Za13*Zl25
     &      *Zl52*ls2+Mm*Sqrt2*Za33*Zl22*Zl55/cb*hL*sb-Mm
     &      *Sqrt2*Za33*Zl25*Zl52/cb*hL*sb+EE/MW*Mm/SW*Za23
     &      *Zl22*Zl55/cb*hL*xvev-EE/MW*Mm/SW*Za23*Zl25*Zl52
     &      /cb*hL*xvev
      AAABR(57) = Sqrt2*Za13*Zl11*Zl44*ls1-Sqrt2*Za13*Zl14
     &      *Zl41*ls1+Me*Sqrt2*Za33*Zl11*Zl44/cb*hL*sb-Me
     &      *Sqrt2*Za33*Zl14*Zl41/cb*hL*sb+EE/MW*Me/SW*Za23
     &      *Zl11*Zl44/cb*hL*xvev-EE/MW*Me/SW*Za23*Zl14*Zl41
     &      /cb*hL*xvev
      AAABR(58) = Sqrt2*Za13*Zl33*Zl66*ls3-Sqrt2*Za13*Zl36
     &      *Zl63*ls3+Ml*Sqrt2*Za33*Zl33*Zl66/cb*hL*sb-Ml
     &      *Sqrt2*Za33*Zl36*Zl63/cb*hL*sb+EE/MW*Ml/SW*Za23
     &      *Zl33*Zl66/cb*hL*xvev-EE/MW*Ml/SW*Za23*Zl36*Zl63
     &      /cb*hL*xvev
      AAABR(59) = Sqrt2*Za11*Zl33*Zl66*ls3-Sqrt2*Za11*Zl36
     &      *Zl63*ls3+Ml*Sqrt2*Za31*Zl33*Zl66/cb*hL*sb-Ml
     &      *Sqrt2*Za31*Zl36*Zl63/cb*hL*sb+EE/MW*Ml/SW*Za21
     &      *Zl33*Zl66/cb*hL*xvev-EE/MW*Ml/SW*Za21*Zl36*Zl63
     &      /cb*hL*xvev
      AAABR(60) = Sqrt2*Za12*Zl22*Zl55*ls2-Sqrt2*Za12*Zl25
     &      *Zl52*ls2+Mm*Sqrt2*Za32*Zl22*Zl55/cb*hL*sb-Mm
     &      *Sqrt2*Za32*Zl25*Zl52/cb*hL*sb+EE/MW*Mm/SW*Za22
     &      *Zl22*Zl55/cb*hL*xvev-EE/MW*Mm/SW*Za22*Zl25*Zl52
     &      /cb*hL*xvev
      AAABR(61) = Sqrt2*Za11*Zl22*Zl55*ls2-Sqrt2*Za11*Zl25
     &      *Zl52*ls2+Mm*Sqrt2*Za31*Zl22*Zl55/cb*hL*sb-Mm
     &      *Sqrt2*Za31*Zl25*Zl52/cb*hL*sb+EE/MW*Mm/SW*Za21
     &      *Zl22*Zl55/cb*hL*xvev-EE/MW*Mm/SW*Za21*Zl25*Zl52
     &      /cb*hL*xvev
      AAABR(62) = Sqrt2*Za11*Zl11*Zl44*ls1-Sqrt2*Za11*Zl14
     &      *Zl41*ls1+Me*Sqrt2*Za31*Zl11*Zl44/cb*hL*sb-Me
     &      *Sqrt2*Za31*Zl14*Zl41/cb*hL*sb+EE/MW*Me/SW*Za21
     &      *Zl11*Zl44/cb*hL*xvev-EE/MW*Me/SW*Za21*Zl14*Zl41
     &      /cb*hL*xvev
      AAABR(63) = Sqrt2*Za12*Zl33*Zl66*ls3-Sqrt2*Za12*Zl36
     &      *Zl63*ls3+Ml*Sqrt2*Za32*Zl33*Zl66/cb*hL*sb-Ml
     &      *Sqrt2*Za32*Zl36*Zl63/cb*hL*sb+EE/MW*Ml/SW*Za22
     &      *Zl33*Zl66/cb*hL*xvev-EE/MW*Ml/SW*Za22*Zl36*Zl63
     &      /cb*hL*xvev
      AAABR(64) = Sqrt2*Za12*Zl11*Zl44*ls1-Sqrt2*Za12*Zl14
     &      *Zl41*ls1+Me*Sqrt2*Za32*Zl11*Zl44/cb*hL*sb-Me
     &      *Sqrt2*Za32*Zl14*Zl41/cb*hL*sb+EE/MW*Me/SW*Za22
     &      *Zl11*Zl44/cb*hL*xvev-EE/MW*Me/SW*Za22*Zl14*Zl41
     &      /cb*hL*xvev
      AAABR(65) = 2D0*Sqrt2*Zh11*Zl11*Zl41*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl11**2/cb+EE
     &      *MW/SW*Zh21*Zl11**2*sb+EE/CW**2*MW*SW*Zh11*Zl11
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl41**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl11**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl41**2*sb+EE*MW/SW*Zh11*Zl11**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl11**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl41**2/cb*sb**2-2D0*Me
     &      *Sqrt2*Zh31*Zl11*Zl41/cb*hL*sb-2D0*EE/MW*Me/SW
     &      *Zh21*Zl11*Zl41/cb*hL*xvev
      AAABR(66) = 2D0*Sqrt2*Zh11*Zl33*Zl63*ls3+2D0*EE/MW
     &      *Ml**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl33**2/cb+EE
     &      *MW/SW*Zh21*Zl33**2*sb+EE/CW**2*MW*SW*Zh11*Zl33
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl63**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl33**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl63**2*sb+EE*MW/SW*Zh11*Zl33**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl33**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl63**2/cb*sb**2-2D0*Ml
     &      *Sqrt2*Zh31*Zl33*Zl63/cb*hL*sb-2D0*EE/MW*Ml/SW
     &      *Zh21*Zl33*Zl63/cb*hL*xvev
      AAABR(67) = 2D0*Sqrt2*Zh12*Zl11*Zl41*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl11**2/cb+EE
     &      *MW/SW*Zh22*Zl11**2*sb+EE/CW**2*MW*SW*Zh12*Zl11
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl41**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl11**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl41**2*sb+EE*MW/SW*Zh12*Zl11**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl11**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl41**2/cb*sb**2-2D0*Me
     &      *Sqrt2*Zh32*Zl11*Zl41/cb*hL*sb-2D0*EE/MW*Me/SW
     &      *Zh22*Zl11*Zl41/cb*hL*xvev
      AAABR(68) = 2D0*Sqrt2*Zh12*Zl33*Zl63*ls3+2D0*EE/MW
     &      *Ml**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl33**2/cb+EE
     &      *MW/SW*Zh22*Zl33**2*sb+EE/CW**2*MW*SW*Zh12*Zl33
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl63**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl33**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl63**2*sb+EE*MW/SW*Zh12*Zl33**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl33**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl63**2/cb*sb**2-2D0*Ml
     &      *Sqrt2*Zh32*Zl33*Zl63/cb*hL*sb-2D0*EE/MW*Ml/SW
     &      *Zh22*Zl33*Zl63/cb*hL*xvev
      AAABR(69) = 2D0*Sqrt2*Zh11*Zl22*Zl52*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl22**2/cb+EE
     &      *MW/SW*Zh21*Zl22**2*sb+EE/CW**2*MW*SW*Zh11*Zl22
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl52**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl22**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl52**2*sb+EE*MW/SW*Zh11*Zl22**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl22**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl52**2/cb*sb**2-2D0*Mm
     &      *Sqrt2*Zh31*Zl22*Zl52/cb*hL*sb-2D0*EE/MW*Mm/SW
     &      *Zh21*Zl22*Zl52/cb*hL*xvev
      AAABR(70) = 2D0*Sqrt2*Zh13*Zl22*Zl52*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl22**2/cb+EE
     &      *MW/SW*Zh23*Zl22**2*sb+EE/CW**2*MW*SW*Zh13*Zl22
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl52**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl22**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl52**2*sb+EE*MW/SW*Zh13*Zl22**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl22**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl52**2/cb*sb**2-2D0*Mm
     &      *Sqrt2*Zh33*Zl22*Zl52/cb*hL*sb-2D0*EE/MW*Mm/SW
     &      *Zh23*Zl22*Zl52/cb*hL*xvev
      AAABR(71) = 2D0*Sqrt2*Zh13*Zl11*Zl41*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl11**2/cb+EE
     &      *MW/SW*Zh23*Zl11**2*sb+EE/CW**2*MW*SW*Zh13*Zl11
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl41**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl11**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl41**2*sb+EE*MW/SW*Zh13*Zl11**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl11**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl41**2/cb*sb**2-2D0*Me
     &      *Sqrt2*Zh33*Zl11*Zl41/cb*hL*sb-2D0*EE/MW*Me/SW
     &      *Zh23*Zl11*Zl41/cb*hL*xvev
      AAABR(72) = 2D0*Sqrt2*Zh12*Zl22*Zl52*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl22**2/cb+EE
     &      *MW/SW*Zh22*Zl22**2*sb+EE/CW**2*MW*SW*Zh12*Zl22
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl52**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl22**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl52**2*sb+EE*MW/SW*Zh12*Zl22**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl22**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl52**2/cb*sb**2-2D0*Mm
     &      *Sqrt2*Zh32*Zl22*Zl52/cb*hL*sb-2D0*EE/MW*Mm/SW
     &      *Zh22*Zl22*Zl52/cb*hL*xvev
      AAABR(73) = 2D0*Sqrt2*Zh13*Zl33*Zl63*ls3+2D0*EE/MW
     &      *Ml**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl33**2/cb+EE
     &      *MW/SW*Zh23*Zl33**2*sb+EE/CW**2*MW*SW*Zh13*Zl33
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl63**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl33**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl63**2*sb+EE*MW/SW*Zh13*Zl33**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl33**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl63**2/cb*sb**2-2D0*Ml
     &      *Sqrt2*Zh33*Zl33*Zl63/cb*hL*sb-2D0*EE/MW*Ml/SW
     &      *Zh23*Zl33*Zl63/cb*hL*xvev
      AAABR(74) = Sqrt2*Zh11*Zl22*Zl55*ls2+Sqrt2*Zh11*Zl25
     &      *Zl52*ls2-EE*MW/SW*Zh11*Zl22*Zl25/cb+EE*MW/SW
     &      *Zh21*Zl22*Zl25*sb+EE/CW**2*MW*SW*Zh11*Zl22*Zl25
     &      /cb-2D0*EE/CW**2*MW*SW*Zh11*Zl52*Zl55/cb-EE/CW
     &      **2*MW*SW*Zh21*Zl22*Zl25*sb+2D0*EE/CW**2*MW*SW
     &      *Zh21*Zl52*Zl55*sb+EE*MW/SW*Zh11*Zl22*Zl25/cb
     &      *sb**2-Mm*Sqrt2*Zh31*Zl22*Zl55/cb*hL*sb-Mm*Sqrt2
     &      *Zh31*Zl25*Zl52/cb*hL*sb-EE/CW**2*MW*SW*Zh11
     &      *Zl22*Zl25/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl52
     &      *Zl55/cb*sb**2-EE/MW*Mm/SW*Zh21*Zl22*Zl55/cb
     &      *hL*xvev-EE/MW*Mm/SW*Zh21*Zl25*Zl52/cb*hL*xvev
      AAABR(75) = Sqrt2*Zh11*Zl33*Zl66*ls3+Sqrt2*Zh11*Zl36
     &      *Zl63*ls3-EE*MW/SW*Zh11*Zl33*Zl36/cb+EE*MW/SW
     &      *Zh21*Zl33*Zl36*sb+EE/CW**2*MW*SW*Zh11*Zl33*Zl36
     &      /cb-2D0*EE/CW**2*MW*SW*Zh11*Zl63*Zl66/cb-EE/CW
     &      **2*MW*SW*Zh21*Zl33*Zl36*sb+2D0*EE/CW**2*MW*SW
     &      *Zh21*Zl63*Zl66*sb+EE*MW/SW*Zh11*Zl33*Zl36/cb
     &      *sb**2-Ml*Sqrt2*Zh31*Zl33*Zl66/cb*hL*sb-Ml*Sqrt2
     &      *Zh31*Zl36*Zl63/cb*hL*sb-EE/CW**2*MW*SW*Zh11
     &      *Zl33*Zl36/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl63
     &      *Zl66/cb*sb**2-EE/MW*Ml/SW*Zh21*Zl33*Zl66/cb
     &      *hL*xvev-EE/MW*Ml/SW*Zh21*Zl36*Zl63/cb*hL*xvev
      AAABR(76) = Sqrt2*Zh13*Zl33*Zl66*ls3+Sqrt2*Zh13*Zl36
     &      *Zl63*ls3-EE*MW/SW*Zh13*Zl33*Zl36/cb+EE*MW/SW
     &      *Zh23*Zl33*Zl36*sb+EE/CW**2*MW*SW*Zh13*Zl33*Zl36
     &      /cb-2D0*EE/CW**2*MW*SW*Zh13*Zl63*Zl66/cb-EE/CW
     &      **2*MW*SW*Zh23*Zl33*Zl36*sb+2D0*EE/CW**2*MW*SW
     &      *Zh23*Zl63*Zl66*sb+EE*MW/SW*Zh13*Zl33*Zl36/cb
     &      *sb**2-Ml*Sqrt2*Zh33*Zl33*Zl66/cb*hL*sb-Ml*Sqrt2
     &      *Zh33*Zl36*Zl63/cb*hL*sb-EE/CW**2*MW*SW*Zh13
     &      *Zl33*Zl36/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh13*Zl63
     &      *Zl66/cb*sb**2-EE/MW*Ml/SW*Zh23*Zl33*Zl66/cb
     &      *hL*xvev-EE/MW*Ml/SW*Zh23*Zl36*Zl63/cb*hL*xvev
      AAABR(77) = Sqrt2*Zh13*Zl22*Zl55*ls2+Sqrt2*Zh13*Zl25
     &      *Zl52*ls2-EE*MW/SW*Zh13*Zl22*Zl25/cb+EE*MW/SW
     &      *Zh23*Zl22*Zl25*sb+EE/CW**2*MW*SW*Zh13*Zl22*Zl25
     &      /cb-2D0*EE/CW**2*MW*SW*Zh13*Zl52*Zl55/cb-EE/CW
     &      **2*MW*SW*Zh23*Zl22*Zl25*sb+2D0*EE/CW**2*MW*SW
     &      *Zh23*Zl52*Zl55*sb+EE*MW/SW*Zh13*Zl22*Zl25/cb
     &      *sb**2-Mm*Sqrt2*Zh33*Zl22*Zl55/cb*hL*sb-Mm*Sqrt2
     &      *Zh33*Zl25*Zl52/cb*hL*sb-EE/CW**2*MW*SW*Zh13
     &      *Zl22*Zl25/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh13*Zl52
     &      *Zl55/cb*sb**2-EE/MW*Mm/SW*Zh23*Zl22*Zl55/cb
     &      *hL*xvev-EE/MW*Mm/SW*Zh23*Zl25*Zl52/cb*hL*xvev
      AAABR(78) = Sqrt2*Zh13*Zl11*Zl44*ls1+Sqrt2*Zh13*Zl14
     &      *Zl41*ls1-EE*MW/SW*Zh13*Zl11*Zl14/cb+EE*MW/SW
     &      *Zh23*Zl11*Zl14*sb+EE/CW**2*MW*SW*Zh13*Zl11*Zl14
     &      /cb-2D0*EE/CW**2*MW*SW*Zh13*Zl41*Zl44/cb-EE/CW
     &      **2*MW*SW*Zh23*Zl11*Zl14*sb+2D0*EE/CW**2*MW*SW
     &      *Zh23*Zl41*Zl44*sb+EE*MW/SW*Zh13*Zl11*Zl14/cb
     &      *sb**2-Me*Sqrt2*Zh33*Zl11*Zl44/cb*hL*sb-Me*Sqrt2
     &      *Zh33*Zl14*Zl41/cb*hL*sb-EE/CW**2*MW*SW*Zh13
     &      *Zl11*Zl14/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh13*Zl41
     &      *Zl44/cb*sb**2-EE/MW*Me/SW*Zh23*Zl11*Zl44/cb
     &      *hL*xvev-EE/MW*Me/SW*Zh23*Zl14*Zl41/cb*hL*xvev
      AAABR(79) = Sqrt2*Zh11*Zl11*Zl44*ls1+Sqrt2*Zh11*Zl14
     &      *Zl41*ls1-EE*MW/SW*Zh11*Zl11*Zl14/cb+EE*MW/SW
     &      *Zh21*Zl11*Zl14*sb+EE/CW**2*MW*SW*Zh11*Zl11*Zl14
     &      /cb-2D0*EE/CW**2*MW*SW*Zh11*Zl41*Zl44/cb-EE/CW
     &      **2*MW*SW*Zh21*Zl11*Zl14*sb+2D0*EE/CW**2*MW*SW
     &      *Zh21*Zl41*Zl44*sb+EE*MW/SW*Zh11*Zl11*Zl14/cb
     &      *sb**2-Me*Sqrt2*Zh31*Zl11*Zl44/cb*hL*sb-Me*Sqrt2
     &      *Zh31*Zl14*Zl41/cb*hL*sb-EE/CW**2*MW*SW*Zh11
     &      *Zl11*Zl14/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh11*Zl41
     &      *Zl44/cb*sb**2-EE/MW*Me/SW*Zh21*Zl11*Zl44/cb
     &      *hL*xvev-EE/MW*Me/SW*Zh21*Zl14*Zl41/cb*hL*xvev
      AAABR(80) = Sqrt2*Zh12*Zl22*Zl55*ls2+Sqrt2*Zh12*Zl25
     &      *Zl52*ls2-EE*MW/SW*Zh12*Zl22*Zl25/cb+EE*MW/SW
     &      *Zh22*Zl22*Zl25*sb+EE/CW**2*MW*SW*Zh12*Zl22*Zl25
     &      /cb-2D0*EE/CW**2*MW*SW*Zh12*Zl52*Zl55/cb-EE/CW
     &      **2*MW*SW*Zh22*Zl22*Zl25*sb+2D0*EE/CW**2*MW*SW
     &      *Zh22*Zl52*Zl55*sb+EE*MW/SW*Zh12*Zl22*Zl25/cb
     &      *sb**2-Mm*Sqrt2*Zh32*Zl22*Zl55/cb*hL*sb-Mm*Sqrt2
     &      *Zh32*Zl25*Zl52/cb*hL*sb-EE/CW**2*MW*SW*Zh12
     &      *Zl22*Zl25/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh12*Zl52
     &      *Zl55/cb*sb**2-EE/MW*Mm/SW*Zh22*Zl22*Zl55/cb
     &      *hL*xvev-EE/MW*Mm/SW*Zh22*Zl25*Zl52/cb*hL*xvev
      AAABR(81) = Sqrt2*Zh12*Zl11*Zl44*ls1+Sqrt2*Zh12*Zl14
     &      *Zl41*ls1-EE*MW/SW*Zh12*Zl11*Zl14/cb+EE*MW/SW
     &      *Zh22*Zl11*Zl14*sb+EE/CW**2*MW*SW*Zh12*Zl11*Zl14
     &      /cb-2D0*EE/CW**2*MW*SW*Zh12*Zl41*Zl44/cb-EE/CW
     &      **2*MW*SW*Zh22*Zl11*Zl14*sb+2D0*EE/CW**2*MW*SW
     &      *Zh22*Zl41*Zl44*sb+EE*MW/SW*Zh12*Zl11*Zl14/cb
     &      *sb**2-Me*Sqrt2*Zh32*Zl11*Zl44/cb*hL*sb-Me*Sqrt2
     &      *Zh32*Zl14*Zl41/cb*hL*sb-EE/CW**2*MW*SW*Zh12
     &      *Zl11*Zl14/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh12*Zl41
     &      *Zl44/cb*sb**2-EE/MW*Me/SW*Zh22*Zl11*Zl44/cb
     &      *hL*xvev-EE/MW*Me/SW*Zh22*Zl14*Zl41/cb*hL*xvev
      AAABR(82) = Sqrt2*Zh12*Zl33*Zl66*ls3+Sqrt2*Zh12*Zl36
     &      *Zl63*ls3-EE*MW/SW*Zh12*Zl33*Zl36/cb+EE*MW/SW
     &      *Zh22*Zl33*Zl36*sb+EE/CW**2*MW*SW*Zh12*Zl33*Zl36
     &      /cb-2D0*EE/CW**2*MW*SW*Zh12*Zl63*Zl66/cb-EE/CW
     &      **2*MW*SW*Zh22*Zl33*Zl36*sb+2D0*EE/CW**2*MW*SW
     &      *Zh22*Zl63*Zl66*sb+EE*MW/SW*Zh12*Zl33*Zl36/cb
     &      *sb**2-Ml*Sqrt2*Zh32*Zl33*Zl66/cb*hL*sb-Ml*Sqrt2
     &      *Zh32*Zl36*Zl63/cb*hL*sb-EE/CW**2*MW*SW*Zh12
     &      *Zl33*Zl36/cb*sb**2+2D0*EE/CW**2*MW*SW*Zh12*Zl63
     &      *Zl66/cb*sb**2-EE/MW*Ml/SW*Zh22*Zl33*Zl66/cb
     &      *hL*xvev-EE/MW*Ml/SW*Zh22*Zl36*Zl63/cb*hL*xvev
      AAABR(83) = 2D0*Zl44/cb*ls1-2D0*Zl44/cb*ls1*sb**2-EE
     &      *MW/SW*Sqrt2*Zl14+EE/MW*Me**2/SW*Sqrt2*Zl14+2D0
     &      *EE*MW/SW*Sqrt2*Zl14*sb**2-EE/MW*Me/SW*Sqrt2
     &      *Zl44/cb*hL*sb*xvev
      AAABR(84) = 2D0*Zl55/cb*ls2-2D0*Zl55/cb*ls2*sb**2-EE
     &      *MW/SW*Sqrt2*Zl25+EE/MW*Mm**2/SW*Sqrt2*Zl25+2D0
     &      *EE*MW/SW*Sqrt2*Zl25*sb**2-EE/MW*Mm/SW*Sqrt2
     &      *Zl55/cb*hL*sb*xvev
      AAABR(85) = 2D0*Zl66/cb*ls3-2D0*Zl66/cb*ls3*sb**2-EE
     &      *MW/SW*Sqrt2*Zl36+EE/MW*Ml**2/SW*Sqrt2*Zl36+2D0
     &      *EE*MW/SW*Sqrt2*Zl36*sb**2-EE/MW*Ml/SW*Sqrt2
     &      *Zl66/cb*hL*sb*xvev
      AAABR(86) = 2D0*Sqrt2*Zh12*Zl36*Zl66*ls3+2D0*EE/MW
     &      *Ml**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl36**2/cb+EE
     &      *MW/SW*Zh22*Zl36**2*sb+EE/CW**2*MW*SW*Zh12*Zl36
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl66**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl36**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl66**2*sb+EE*MW/SW*Zh12*Zl36**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl36**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl66**2/cb*sb**2-2D0*Ml
     &      *Sqrt2*Zh32*Zl36*Zl66/cb*hL*sb-2D0*EE/MW*Ml/SW
     &      *Zh22*Zl36*Zl66/cb*hL*xvev
      AAABR(87) = 2D0*Sqrt2*Zh11*Zl25*Zl55*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl25**2/cb+EE
     &      *MW/SW*Zh21*Zl25**2*sb+EE/CW**2*MW*SW*Zh11*Zl25
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl55**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl25**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl55**2*sb+EE*MW/SW*Zh11*Zl25**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl25**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl55**2/cb*sb**2-2D0*Mm
     &      *Sqrt2*Zh31*Zl25*Zl55/cb*hL*sb-2D0*EE/MW*Mm/SW
     &      *Zh21*Zl25*Zl55/cb*hL*xvev
      AAABR(88) = 2D0*Sqrt2*Zh13*Zl36*Zl66*ls3+2D0*EE/MW
     &      *Ml**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl36**2/cb+EE
     &      *MW/SW*Zh23*Zl36**2*sb+EE/CW**2*MW*SW*Zh13*Zl36
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl66**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl36**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl66**2*sb+EE*MW/SW*Zh13*Zl36**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl36**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl66**2/cb*sb**2-2D0*Ml
     &      *Sqrt2*Zh33*Zl36*Zl66/cb*hL*sb-2D0*EE/MW*Ml/SW
     &      *Zh23*Zl36*Zl66/cb*hL*xvev
      AAABR(89) = 2D0*Sqrt2*Zh13*Zl14*Zl44*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl14**2/cb+EE
     &      *MW/SW*Zh23*Zl14**2*sb+EE/CW**2*MW*SW*Zh13*Zl14
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl44**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl14**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl44**2*sb+EE*MW/SW*Zh13*Zl14**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl14**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl44**2/cb*sb**2-2D0*Me
     &      *Sqrt2*Zh33*Zl14*Zl44/cb*hL*sb-2D0*EE/MW*Me/SW
     &      *Zh23*Zl14*Zl44/cb*hL*xvev
      AAABR(90) = 2D0*Sqrt2*Zh13*Zl25*Zl55*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh13/cb-EE*MW/SW*Zh13*Zl25**2/cb+EE
     &      *MW/SW*Zh23*Zl25**2*sb+EE/CW**2*MW*SW*Zh13*Zl25
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh13*Zl55**2/cb-EE
     &      /CW**2*MW*SW*Zh23*Zl25**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh23*Zl55**2*sb+EE*MW/SW*Zh13*Zl25**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh13*Zl25**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh13*Zl55**2/cb*sb**2-2D0*Mm
     &      *Sqrt2*Zh33*Zl25*Zl55/cb*hL*sb-2D0*EE/MW*Mm/SW
     &      *Zh23*Zl25*Zl55/cb*hL*xvev
      AAABR(91) = 2D0*Sqrt2*Zh12*Zl14*Zl44*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl14**2/cb+EE
     &      *MW/SW*Zh22*Zl14**2*sb+EE/CW**2*MW*SW*Zh12*Zl14
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl44**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl14**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl44**2*sb+EE*MW/SW*Zh12*Zl14**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl14**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl44**2/cb*sb**2-2D0*Me
     &      *Sqrt2*Zh32*Zl14*Zl44/cb*hL*sb-2D0*EE/MW*Me/SW
     &      *Zh22*Zl14*Zl44/cb*hL*xvev
      AAABR(92) = 2D0*Sqrt2*Zh11*Zl36*Zl66*ls3+2D0*EE/MW
     &      *Ml**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl36**2/cb+EE
     &      *MW/SW*Zh21*Zl36**2*sb+EE/CW**2*MW*SW*Zh11*Zl36
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl66**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl36**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl66**2*sb+EE*MW/SW*Zh11*Zl36**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl36**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl66**2/cb*sb**2-2D0*Ml
     &      *Sqrt2*Zh31*Zl36*Zl66/cb*hL*sb-2D0*EE/MW*Ml/SW
     &      *Zh21*Zl36*Zl66/cb*hL*xvev
      AAABR(93) = 2D0*Sqrt2*Zh12*Zl25*Zl55*ls2+2D0*EE/MW
     &      *Mm**2/SW*Zh12/cb-EE*MW/SW*Zh12*Zl25**2/cb+EE
     &      *MW/SW*Zh22*Zl25**2*sb+EE/CW**2*MW*SW*Zh12*Zl25
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh12*Zl55**2/cb-EE
     &      /CW**2*MW*SW*Zh22*Zl25**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh22*Zl55**2*sb+EE*MW/SW*Zh12*Zl25**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh12*Zl25**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh12*Zl55**2/cb*sb**2-2D0*Mm
     &      *Sqrt2*Zh32*Zl25*Zl55/cb*hL*sb-2D0*EE/MW*Mm/SW
     &      *Zh22*Zl25*Zl55/cb*hL*xvev
      AAABR(94) = 2D0*Sqrt2*Zh11*Zl14*Zl44*ls1+2D0*EE/MW
     &      *Me**2/SW*Zh11/cb-EE*MW/SW*Zh11*Zl14**2/cb+EE
     &      *MW/SW*Zh21*Zl14**2*sb+EE/CW**2*MW*SW*Zh11*Zl14
     &      **2/cb-2D0*EE/CW**2*MW*SW*Zh11*Zl44**2/cb-EE
     &      /CW**2*MW*SW*Zh21*Zl14**2*sb+2D0*EE/CW**2*MW
     &      *SW*Zh21*Zl44**2*sb+EE*MW/SW*Zh11*Zl14**2/cb
     &      *sb**2-EE/CW**2*MW*SW*Zh11*Zl14**2/cb*sb**2+2D0
     &      *EE/CW**2*MW*SW*Zh11*Zl44**2/cb*sb**2-2D0*Me
     &      *Sqrt2*Zh31*Zl14*Zl44/cb*hL*sb-2D0*EE/MW*Me/SW
     &      *Zh21*Zl14*Zl44/cb*hL*xvev
      AAABR(95) = EE*MW*(Zh13/SW*cb-Zh23/SW*sb+SW/CW**2*Zh13
     &      *cb-SW/CW**2*Zh23*sb)
      AAABR(96) = EE*MW*(Zh11/SW*cb-Zh21/SW*sb+SW/CW**2*Zh11
     &      *cb-SW/CW**2*Zh21*sb)
      AAABR(97) = EE*MW*(Zh12/SW*cb-Zh22/SW*sb+SW/CW**2*Zh12
     &      *cb-SW/CW**2*Zh22*sb)
      AAABR(98) = Vus*(2D0*Zd22*Zu41*sb*us1+2D0*Zd52*Zu11
     &      /cb*ds2-2D0*Zd52*Zu11/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd22*Zu11+EE/MW*Ms**2/SW*Sqrt2*Zd22*Zu11+2D0
     &      *EE*MW/SW*Sqrt2*Zd22*Zu11*sb**2-EE/MW*Ms/SW*Sqrt2
     &      *Zd52*Zu11/cb*hL*sb*xvev)
      AAABR(99) = Vcs*(2D0*Zd22*Zu52*sb*us2+2D0*Zd52*Zu22
     &      /cb*ds2-2D0*Zd52*Zu22/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd22*Zu22-EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu22+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd22*Zu22+2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu22*sb**2+EE/MW*Mc/SW*Sqrt2*Zd22*Zu52
     &      /cb*hL/sb*xvev-EE/MW*Mc/SW*Sqrt2*Zd22*Zu52/cb
     &      *hL*sb*xvev-EE/MW*Ms/SW*Sqrt2*Zd52*Zu22/cb*hL
     &      *sb*xvev)
      AAABR(100) = Vud*(2D0*Zd11*Zu41*sb*us1+2D0*Zd41*Zu11
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd11*Zu11+2D0*EE*MW/SW
     &      *Sqrt2*Zd11*Zu11*sb**2)
      AAABR(101) = Vcd*(2D0*Zd11*Zu52*sb*us2+2D0*Zd41*Zu22
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd11*Zu22-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd11*Zu22+2D0*EE*MW/SW*Sqrt2*Zd11*Zu22
     &      *sb**2+EE/MW*Mc/SW*Sqrt2*Zd11*Zu52*cb*hL/sb*xvev
     &      )
      AAABR(102) = 2D0*Zd33*Zu63*sb*us3+2D0*Zd63*Zu33/cb
     &      *ds3-2D0*Zd63*Zu33/cb*ds3*sb**2-EE*MW/SW*Sqrt2
     &      *Zd33*Zu33+EE/MW*Mb**2/SW*Sqrt2*Zd33*Zu33-EE
     &      /MW*Mt**2/SW*Sqrt2*Zd33*Zu33+2D0*EE*MW/SW*Sqrt2
     &      *Zd33*Zu33*sb**2-EE/MW*Mb/SW*Sqrt2*Zd63*Zu33
     &      /cb*hL*sb*xvev+EE/MW*Mt/SW*Sqrt2*Zd33*Zu63/cb
     &      *hL/sb*xvev-EE/MW*Mt/SW*Sqrt2*Zd33*Zu63/cb*hL
     &      *sb*xvev
      AAABR(103) = 2D0*Zd33*Zu66*sb*us3+2D0*Zd63*Zu36/cb
     &      *ds3-2D0*Zd63*Zu36/cb*ds3*sb**2-EE*MW/SW*Sqrt2
     &      *Zd33*Zu36+EE/MW*Mb**2/SW*Sqrt2*Zd33*Zu36-EE
     &      /MW*Mt**2/SW*Sqrt2*Zd33*Zu36+2D0*EE*MW/SW*Sqrt2
     &      *Zd33*Zu36*sb**2-EE/MW*Mb/SW*Sqrt2*Zd63*Zu36
     &      /cb*hL*sb*xvev+EE/MW*Mt/SW*Sqrt2*Zd33*Zu66/cb
     &      *hL/sb*xvev-EE/MW*Mt/SW*Sqrt2*Zd33*Zu66/cb*hL
     &      *sb*xvev
      AAABR(104) = Vud*(2D0*Zd11*Zu44*sb*us1+2D0*Zd41*Zu14
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd11*Zu14+2D0*EE*MW/SW
     &      *Sqrt2*Zd11*Zu14*sb**2)
      AAABR(105) = Vus*(2D0*Zd22*Zu44*sb*us1+2D0*Zd52*Zu14
     &      /cb*ds2-2D0*Zd52*Zu14/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd22*Zu14+EE/MW*Ms**2/SW*Sqrt2*Zd22*Zu14+2D0
     &      *EE*MW/SW*Sqrt2*Zd22*Zu14*sb**2-EE/MW*Ms/SW*Sqrt2
     &      *Zd52*Zu14/cb*hL*sb*xvev)
      AAABR(106) = Vcs*(2D0*Zd22*Zu55*sb*us2+2D0*Zd52*Zu25
     &      /cb*ds2-2D0*Zd52*Zu25/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd22*Zu25-EE/MW*Mc**2/SW*Sqrt2*Zd22*Zu25+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd22*Zu25+2D0*EE*MW/SW*Sqrt2
     &      *Zd22*Zu25*sb**2+EE/MW*Mc/SW*Sqrt2*Zd22*Zu55
     &      /cb*hL/sb*xvev-EE/MW*Mc/SW*Sqrt2*Zd22*Zu55/cb
     &      *hL*sb*xvev-EE/MW*Ms/SW*Sqrt2*Zd52*Zu25/cb*hL
     &      *sb*xvev)
      AAABR(107) = Vcd*(2D0*Zd11*Zu55*sb*us2+2D0*Zd41*Zu25
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd11*Zu25-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd11*Zu25+2D0*EE*MW/SW*Sqrt2*Zd11*Zu25
     &      *sb**2+EE/MW*Mc/SW*Sqrt2*Zd11*Zu55*cb*hL/sb*xvev
     &      )
      AAABR(108) = Sqrt2*Za13*Zd22*Zd55*ds2-Sqrt2*Za13*Zd25
     &      *Zd52*ds2+Ms*Sqrt2*Za33*Zd22*Zd55/cb*hL*sb-Ms
     &      *Sqrt2*Za33*Zd25*Zd52/cb*hL*sb+EE/MW*Ms/SW*Za23
     &      *Zd22*Zd55/cb*hL*xvev-EE/MW*Ms/SW*Za23*Zd25*Zd52
     &      /cb*hL*xvev
      AAABR(109) = Sqrt2*Za13*Zd33*Zd66*ds3-Sqrt2*Za13*Zd36
     &      *Zd63*ds3+Mb*Sqrt2*Za33*Zd33*Zd66/cb*hL*sb-Mb
     &      *Sqrt2*Za33*Zd36*Zd63/cb*hL*sb+EE/MW*Mb/SW*Za23
     &      *Zd33*Zd66/cb*hL*xvev-EE/MW*Mb/SW*Za23*Zd36*Zd63
     &      /cb*hL*xvev
      AAABR(110) = Sqrt2*Za13*ds1*(Zd11*Zd44-Zd14*Zd41)
      AAABR(111) = Sqrt2*Za12*Zd33*Zd66*ds3-Sqrt2*Za12*Zd36
     &      *Zd63*ds3+Mb*Sqrt2*Za32*Zd33*Zd66/cb*hL*sb-Mb
     &      *Sqrt2*Za32*Zd36*Zd63/cb*hL*sb+EE/MW*Mb/SW*Za22
     &      *Zd33*Zd66/cb*hL*xvev-EE/MW*Mb/SW*Za22*Zd36*Zd63
     &      /cb*hL*xvev
      AAABR(112) = Sqrt2*Za12*ds1*(Zd11*Zd44-Zd14*Zd41)
      AAABR(113) = Sqrt2*Za11*Zd33*Zd66*ds3-Sqrt2*Za11*Zd36
     &      *Zd63*ds3+Mb*Sqrt2*Za31*Zd33*Zd66/cb*hL*sb-Mb
     &      *Sqrt2*Za31*Zd36*Zd63/cb*hL*sb+EE/MW*Mb/SW*Za21
     &      *Zd33*Zd66/cb*hL*xvev-EE/MW*Mb/SW*Za21*Zd36*Zd63
     &      /cb*hL*xvev
      AAABR(114) = Sqrt2*Za12*Zd22*Zd55*ds2-Sqrt2*Za12*Zd25
     &      *Zd52*ds2+Ms*Sqrt2*Za32*Zd22*Zd55/cb*hL*sb-Ms
     &      *Sqrt2*Za32*Zd25*Zd52/cb*hL*sb+EE/MW*Ms/SW*Za22
     &      *Zd22*Zd55/cb*hL*xvev-EE/MW*Ms/SW*Za22*Zd25*Zd52
     &      /cb*hL*xvev
      AAABR(115) = Sqrt2*Za11*ds1*(Zd11*Zd44-Zd14*Zd41)
      AAABR(116) = Sqrt2*Za11*Zd22*Zd55*ds2-Sqrt2*Za11*Zd25
     &      *Zd52*ds2+Ms*Sqrt2*Za31*Zd22*Zd55/cb*hL*sb-Ms
     &      *Sqrt2*Za31*Zd25*Zd52/cb*hL*sb+EE/MW*Ms/SW*Za21
     &      *Zd22*Zd55/cb*hL*xvev-EE/MW*Ms/SW*Za21*Zd25*Zd52
     &      /cb*hL*xvev
      AAABR(117) = 6D0*Sqrt2*Zd22*Zd52*Zh11*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh11/cb-3D0*EE*MW/SW*Zd22**2*Zh11/cb
     &      +3D0*EE*MW/SW*Zd22**2*Zh21*sb-EE/CW**2*MW*SW
     &      *Zd22**2*Zh11/cb+EE/CW**2*MW*SW*Zd22**2*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zd52**2*Zh11/cb+2D0*EE
     &      /CW**2*MW*SW*Zd52**2*Zh21*sb+3D0*EE*MW/SW*Zd22
     &      **2*Zh11/cb*sb**2+EE/CW**2*MW*SW*Zd22**2*Zh11
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd52**2*Zh11/cb
     &      *sb**2-6D0*Ms*Sqrt2*Zd22*Zd52*Zh31/cb*hL*sb-6D0
     &      *EE/MW*Ms/SW*Zd22*Zd52*Zh21/cb*hL*xvev
      AAABR(118) = 6D0*Sqrt2*Zd11*Zd41*Zh12*ds1-3D0*EE*MW
     &      /SW*Zd11**2*Zh12*cb+3D0*EE*MW/SW*Zd11**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zd11**2*Zh12*cb+EE/CW**2*MW
     &      *SW*Zd11**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zd41**2
     &      *Zh12*cb+2D0*EE/CW**2*MW*SW*Zd41**2*Zh22*sb
      AAABR(119) = 6D0*Sqrt2*Zd22*Zd52*Zh13*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh13/cb-3D0*EE*MW/SW*Zd22**2*Zh13/cb
     &      +3D0*EE*MW/SW*Zd22**2*Zh23*sb-EE/CW**2*MW*SW
     &      *Zd22**2*Zh13/cb+EE/CW**2*MW*SW*Zd22**2*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zd52**2*Zh13/cb+2D0*EE
     &      /CW**2*MW*SW*Zd52**2*Zh23*sb+3D0*EE*MW/SW*Zd22
     &      **2*Zh13/cb*sb**2+EE/CW**2*MW*SW*Zd22**2*Zh13
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd52**2*Zh13/cb
     &      *sb**2-6D0*Ms*Sqrt2*Zd22*Zd52*Zh33/cb*hL*sb-6D0
     &      *EE/MW*Ms/SW*Zd22*Zd52*Zh23/cb*hL*xvev
      AAABR(120) = 6D0*Sqrt2*Zd33*Zd63*Zh11*ds3+6D0*EE/MW
     &      *Mb**2/SW*Zh11/cb-3D0*EE*MW/SW*Zd33**2*Zh11/cb
     &      +3D0*EE*MW/SW*Zd33**2*Zh21*sb-EE/CW**2*MW*SW
     &      *Zd33**2*Zh11/cb+EE/CW**2*MW*SW*Zd33**2*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zd63**2*Zh11/cb+2D0*EE
     &      /CW**2*MW*SW*Zd63**2*Zh21*sb+3D0*EE*MW/SW*Zd33
     &      **2*Zh11/cb*sb**2+EE/CW**2*MW*SW*Zd33**2*Zh11
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd63**2*Zh11/cb
     &      *sb**2-6D0*Mb*Sqrt2*Zd33*Zd63*Zh31/cb*hL*sb-6D0
     &      *EE/MW*Mb/SW*Zd33*Zd63*Zh21/cb*hL*xvev
      AAABR(121) = 6D0*Sqrt2*Zd11*Zd41*Zh13*ds1-3D0*EE*MW
     &      /SW*Zd11**2*Zh13*cb+3D0*EE*MW/SW*Zd11**2*Zh23
     &      *sb-EE/CW**2*MW*SW*Zd11**2*Zh13*cb+EE/CW**2*MW
     &      *SW*Zd11**2*Zh23*sb-2D0*EE/CW**2*MW*SW*Zd41**2
     &      *Zh13*cb+2D0*EE/CW**2*MW*SW*Zd41**2*Zh23*sb
      AAABR(122) = 6D0*Sqrt2*Zd33*Zd63*Zh13*ds3+6D0*EE/MW
     &      *Mb**2/SW*Zh13/cb-3D0*EE*MW/SW*Zd33**2*Zh13/cb
     &      +3D0*EE*MW/SW*Zd33**2*Zh23*sb-EE/CW**2*MW*SW
     &      *Zd33**2*Zh13/cb+EE/CW**2*MW*SW*Zd33**2*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zd63**2*Zh13/cb+2D0*EE
     &      /CW**2*MW*SW*Zd63**2*Zh23*sb+3D0*EE*MW/SW*Zd33
     &      **2*Zh13/cb*sb**2+EE/CW**2*MW*SW*Zd33**2*Zh13
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd63**2*Zh13/cb
     &      *sb**2-6D0*Mb*Sqrt2*Zd33*Zd63*Zh33/cb*hL*sb-6D0
     &      *EE/MW*Mb/SW*Zd33*Zd63*Zh23/cb*hL*xvev
      AAABR(123) = 6D0*Sqrt2*Zd33*Zd63*Zh12*ds3+6D0*EE/MW
     &      *Mb**2/SW*Zh12/cb-3D0*EE*MW/SW*Zd33**2*Zh12/cb
     &      +3D0*EE*MW/SW*Zd33**2*Zh22*sb-EE/CW**2*MW*SW
     &      *Zd33**2*Zh12/cb+EE/CW**2*MW*SW*Zd33**2*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zd63**2*Zh12/cb+2D0*EE
     &      /CW**2*MW*SW*Zd63**2*Zh22*sb+3D0*EE*MW/SW*Zd33
     &      **2*Zh12/cb*sb**2+EE/CW**2*MW*SW*Zd33**2*Zh12
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd63**2*Zh12/cb
     &      *sb**2-6D0*Mb*Sqrt2*Zd33*Zd63*Zh32/cb*hL*sb-6D0
     &      *EE/MW*Mb/SW*Zd33*Zd63*Zh22/cb*hL*xvev
      AAABR(124) = 6D0*Sqrt2*Zd22*Zd52*Zh12*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh12/cb-3D0*EE*MW/SW*Zd22**2*Zh12/cb
     &      +3D0*EE*MW/SW*Zd22**2*Zh22*sb-EE/CW**2*MW*SW
     &      *Zd22**2*Zh12/cb+EE/CW**2*MW*SW*Zd22**2*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zd52**2*Zh12/cb+2D0*EE
     &      /CW**2*MW*SW*Zd52**2*Zh22*sb+3D0*EE*MW/SW*Zd22
     &      **2*Zh12/cb*sb**2+EE/CW**2*MW*SW*Zd22**2*Zh12
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd52**2*Zh12/cb
     &      *sb**2-6D0*Ms*Sqrt2*Zd22*Zd52*Zh32/cb*hL*sb-6D0
     &      *EE/MW*Ms/SW*Zd22*Zd52*Zh22/cb*hL*xvev
      AAABR(125) = 6D0*Sqrt2*Zd11*Zd41*Zh11*ds1-3D0*EE*MW
     &      /SW*Zd11**2*Zh11*cb+3D0*EE*MW/SW*Zd11**2*Zh21
     &      *sb-EE/CW**2*MW*SW*Zd11**2*Zh11*cb+EE/CW**2*MW
     &      *SW*Zd11**2*Zh21*sb-2D0*EE/CW**2*MW*SW*Zd41**2
     &      *Zh11*cb+2D0*EE/CW**2*MW*SW*Zd41**2*Zh21*sb
      AAABR(126) = 3D0*Sqrt2*Zd33*Zd66*Zh11*ds3+3D0*Sqrt2
     &      *Zd36*Zd63*Zh11*ds3-3D0*EE*MW/SW*Zd33*Zd36*Zh11
     &      /cb+3D0*EE*MW/SW*Zd33*Zd36*Zh21*sb-EE/CW**2*MW
     &      *SW*Zd33*Zd36*Zh11/cb+EE/CW**2*MW*SW*Zd33*Zd36
     &      *Zh21*sb-2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh11/cb
     &      +2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh21*sb+3D0*EE
     &      *MW/SW*Zd33*Zd36*Zh11/cb*sb**2-3D0*Mb*Sqrt2*Zd33
     &      *Zd66*Zh31/cb*hL*sb-3D0*Mb*Sqrt2*Zd36*Zd63*Zh31
     &      /cb*hL*sb+EE/CW**2*MW*SW*Zd33*Zd36*Zh11/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh11/cb*sb**2
     &      -3D0*EE/MW*Mb/SW*Zd33*Zd66*Zh21/cb*hL*xvev-3D0
     &      *EE/MW*Mb/SW*Zd36*Zd63*Zh21/cb*hL*xvev
      AAABR(127) = 3D0*Sqrt2*Zd22*Zd55*Zh11*ds2+3D0*Sqrt2
     &      *Zd25*Zd52*Zh11*ds2-3D0*EE*MW/SW*Zd22*Zd25*Zh11
     &      /cb+3D0*EE*MW/SW*Zd22*Zd25*Zh21*sb-EE/CW**2*MW
     &      *SW*Zd22*Zd25*Zh11/cb+EE/CW**2*MW*SW*Zd22*Zd25
     &      *Zh21*sb-2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh11/cb
     &      +2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh21*sb+3D0*EE
     &      *MW/SW*Zd22*Zd25*Zh11/cb*sb**2-3D0*Ms*Sqrt2*Zd22
     &      *Zd55*Zh31/cb*hL*sb-3D0*Ms*Sqrt2*Zd25*Zd52*Zh31
     &      /cb*hL*sb+EE/CW**2*MW*SW*Zd22*Zd25*Zh11/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh11/cb*sb**2
     &      -3D0*EE/MW*Ms/SW*Zd22*Zd55*Zh21/cb*hL*xvev-3D0
     &      *EE/MW*Ms/SW*Zd25*Zd52*Zh21/cb*hL*xvev
      AAABR(128) = 3D0*Sqrt2*Zd22*Zd55*Zh13*ds2+3D0*Sqrt2
     &      *Zd25*Zd52*Zh13*ds2-3D0*EE*MW/SW*Zd22*Zd25*Zh13
     &      /cb+3D0*EE*MW/SW*Zd22*Zd25*Zh23*sb-EE/CW**2*MW
     &      *SW*Zd22*Zd25*Zh13/cb+EE/CW**2*MW*SW*Zd22*Zd25
     &      *Zh23*sb-2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh13/cb
     &      +2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh23*sb+3D0*EE
     &      *MW/SW*Zd22*Zd25*Zh13/cb*sb**2-3D0*Ms*Sqrt2*Zd22
     &      *Zd55*Zh33/cb*hL*sb-3D0*Ms*Sqrt2*Zd25*Zd52*Zh33
     &      /cb*hL*sb+EE/CW**2*MW*SW*Zd22*Zd25*Zh13/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh13/cb*sb**2
     &      -3D0*EE/MW*Ms/SW*Zd22*Zd55*Zh23/cb*hL*xvev-3D0
     &      *EE/MW*Ms/SW*Zd25*Zd52*Zh23/cb*hL*xvev
      AAABR(129) = 3D0*Sqrt2*Zd33*Zd66*Zh13*ds3+3D0*Sqrt2
     &      *Zd36*Zd63*Zh13*ds3-3D0*EE*MW/SW*Zd33*Zd36*Zh13
     &      /cb+3D0*EE*MW/SW*Zd33*Zd36*Zh23*sb-EE/CW**2*MW
     &      *SW*Zd33*Zd36*Zh13/cb+EE/CW**2*MW*SW*Zd33*Zd36
     &      *Zh23*sb-2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh13/cb
     &      +2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh23*sb+3D0*EE
     &      *MW/SW*Zd33*Zd36*Zh13/cb*sb**2-3D0*Mb*Sqrt2*Zd33
     &      *Zd66*Zh33/cb*hL*sb-3D0*Mb*Sqrt2*Zd36*Zd63*Zh33
     &      /cb*hL*sb+EE/CW**2*MW*SW*Zd33*Zd36*Zh13/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh13/cb*sb**2
     &      -3D0*EE/MW*Mb/SW*Zd33*Zd66*Zh23/cb*hL*xvev-3D0
     &      *EE/MW*Mb/SW*Zd36*Zd63*Zh23/cb*hL*xvev
      AAABR(130) = 3D0*Sqrt2*Zd11*Zd44*Zh13*ds1+3D0*Sqrt2
     &      *Zd14*Zd41*Zh13*ds1-3D0*EE*MW/SW*Zd11*Zd14*Zh13
     &      *cb+3D0*EE*MW/SW*Zd11*Zd14*Zh23*sb-EE/CW**2*MW
     &      *SW*Zd11*Zd14*Zh13*cb+EE/CW**2*MW*SW*Zd11*Zd14
     &      *Zh23*sb-2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh13*cb
     &      +2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh23*sb
      AAABR(131) = 3D0*Sqrt2*Zd22*Zd55*Zh12*ds2+3D0*Sqrt2
     &      *Zd25*Zd52*Zh12*ds2-3D0*EE*MW/SW*Zd22*Zd25*Zh12
     &      /cb+3D0*EE*MW/SW*Zd22*Zd25*Zh22*sb-EE/CW**2*MW
     &      *SW*Zd22*Zd25*Zh12/cb+EE/CW**2*MW*SW*Zd22*Zd25
     &      *Zh22*sb-2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh12/cb
     &      +2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh22*sb+3D0*EE
     &      *MW/SW*Zd22*Zd25*Zh12/cb*sb**2-3D0*Ms*Sqrt2*Zd22
     &      *Zd55*Zh32/cb*hL*sb-3D0*Ms*Sqrt2*Zd25*Zd52*Zh32
     &      /cb*hL*sb+EE/CW**2*MW*SW*Zd22*Zd25*Zh12/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd52*Zd55*Zh12/cb*sb**2
     &      -3D0*EE/MW*Ms/SW*Zd22*Zd55*Zh22/cb*hL*xvev-3D0
     &      *EE/MW*Ms/SW*Zd25*Zd52*Zh22/cb*hL*xvev
      AAABR(132) = 3D0*Sqrt2*Zd33*Zd66*Zh12*ds3+3D0*Sqrt2
     &      *Zd36*Zd63*Zh12*ds3-3D0*EE*MW/SW*Zd33*Zd36*Zh12
     &      /cb+3D0*EE*MW/SW*Zd33*Zd36*Zh22*sb-EE/CW**2*MW
     &      *SW*Zd33*Zd36*Zh12/cb+EE/CW**2*MW*SW*Zd33*Zd36
     &      *Zh22*sb-2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh12/cb
     &      +2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh22*sb+3D0*EE
     &      *MW/SW*Zd33*Zd36*Zh12/cb*sb**2-3D0*Mb*Sqrt2*Zd33
     &      *Zd66*Zh32/cb*hL*sb-3D0*Mb*Sqrt2*Zd36*Zd63*Zh32
     &      /cb*hL*sb+EE/CW**2*MW*SW*Zd33*Zd36*Zh12/cb*sb
     &      **2+2D0*EE/CW**2*MW*SW*Zd63*Zd66*Zh12/cb*sb**2
     &      -3D0*EE/MW*Mb/SW*Zd33*Zd66*Zh22/cb*hL*xvev-3D0
     &      *EE/MW*Mb/SW*Zd36*Zd63*Zh22/cb*hL*xvev
      AAABR(133) = 3D0*Sqrt2*Zd11*Zd44*Zh12*ds1+3D0*Sqrt2
     &      *Zd14*Zd41*Zh12*ds1-3D0*EE*MW/SW*Zd11*Zd14*Zh12
     &      *cb+3D0*EE*MW/SW*Zd11*Zd14*Zh22*sb-EE/CW**2*MW
     &      *SW*Zd11*Zd14*Zh12*cb+EE/CW**2*MW*SW*Zd11*Zd14
     &      *Zh22*sb-2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh12*cb
     &      +2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh22*sb
      AAABR(134) = 3D0*Sqrt2*Zd11*Zd44*Zh11*ds1+3D0*Sqrt2
     &      *Zd14*Zd41*Zh11*ds1-3D0*EE*MW/SW*Zd11*Zd14*Zh11
     &      *cb+3D0*EE*MW/SW*Zd11*Zd14*Zh21*sb-EE/CW**2*MW
     &      *SW*Zd11*Zd14*Zh11*cb+EE/CW**2*MW*SW*Zd11*Zd14
     &      *Zh21*sb-2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh11*cb
     &      +2D0*EE/CW**2*MW*SW*Zd41*Zd44*Zh21*sb
      AAABR(135) = Vus*(2D0*Zd25*Zu41*sb*us1+2D0*Zd55*Zu11
     &      /cb*ds2-2D0*Zd55*Zu11/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd25*Zu11+EE/MW*Ms**2/SW*Sqrt2*Zd25*Zu11+2D0
     &      *EE*MW/SW*Sqrt2*Zd25*Zu11*sb**2-EE/MW*Ms/SW*Sqrt2
     &      *Zd55*Zu11/cb*hL*sb*xvev)
      AAABR(136) = Vud*(2D0*Zd14*Zu41*sb*us1+2D0*Zd44*Zu11
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd14*Zu11+2D0*EE*MW/SW
     &      *Sqrt2*Zd14*Zu11*sb**2)
      AAABR(137) = Vcs*(2D0*Zd25*Zu52*sb*us2+2D0*Zd55*Zu22
     &      /cb*ds2-2D0*Zd55*Zu22/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd25*Zu22-EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu22+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd25*Zu22+2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu22*sb**2+EE/MW*Mc/SW*Sqrt2*Zd25*Zu52
     &      /cb*hL/sb*xvev-EE/MW*Mc/SW*Sqrt2*Zd25*Zu52/cb
     &      *hL*sb*xvev-EE/MW*Ms/SW*Sqrt2*Zd55*Zu22/cb*hL
     &      *sb*xvev)
      AAABR(138) = 2D0*Zd36*Zu63*sb*us3+2D0*Zd66*Zu33/cb
     &      *ds3-2D0*Zd66*Zu33/cb*ds3*sb**2-EE*MW/SW*Sqrt2
     &      *Zd36*Zu33+EE/MW*Mb**2/SW*Sqrt2*Zd36*Zu33-EE
     &      /MW*Mt**2/SW*Sqrt2*Zd36*Zu33+2D0*EE*MW/SW*Sqrt2
     &      *Zd36*Zu33*sb**2-EE/MW*Mb/SW*Sqrt2*Zd66*Zu33
     &      /cb*hL*sb*xvev+EE/MW*Mt/SW*Sqrt2*Zd36*Zu63/cb
     &      *hL/sb*xvev-EE/MW*Mt/SW*Sqrt2*Zd36*Zu63/cb*hL
     &      *sb*xvev
      AAABR(139) = Vcd*(2D0*Zd14*Zu52*sb*us2+2D0*Zd44*Zu22
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd14*Zu22-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd14*Zu22+2D0*EE*MW/SW*Sqrt2*Zd14*Zu22
     &      *sb**2+EE/MW*Mc/SW*Sqrt2*Zd14*Zu52*cb*hL/sb*xvev
     &      )
      AAABR(140) = 2D0*Zd36*Zu66*sb*us3+2D0*Zd66*Zu36/cb
     &      *ds3-2D0*Zd66*Zu36/cb*ds3*sb**2-EE*MW/SW*Sqrt2
     &      *Zd36*Zu36+EE/MW*Mb**2/SW*Sqrt2*Zd36*Zu36-EE
     &      /MW*Mt**2/SW*Sqrt2*Zd36*Zu36+2D0*EE*MW/SW*Sqrt2
     &      *Zd36*Zu36*sb**2-EE/MW*Mb/SW*Sqrt2*Zd66*Zu36
     &      /cb*hL*sb*xvev+EE/MW*Mt/SW*Sqrt2*Zd36*Zu66/cb
     &      *hL/sb*xvev-EE/MW*Mt/SW*Sqrt2*Zd36*Zu66/cb*hL
     &      *sb*xvev
      AAABR(141) = Vud*(2D0*Zd14*Zu44*sb*us1+2D0*Zd44*Zu14
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd14*Zu14+2D0*EE*MW/SW
     &      *Sqrt2*Zd14*Zu14*sb**2)
      AAABR(142) = Vus*(2D0*Zd25*Zu44*sb*us1+2D0*Zd55*Zu14
     &      /cb*ds2-2D0*Zd55*Zu14/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd25*Zu14+EE/MW*Ms**2/SW*Sqrt2*Zd25*Zu14+2D0
     &      *EE*MW/SW*Sqrt2*Zd25*Zu14*sb**2-EE/MW*Ms/SW*Sqrt2
     &      *Zd55*Zu14/cb*hL*sb*xvev)
      AAABR(143) = Vcs*(2D0*Zd25*Zu55*sb*us2+2D0*Zd55*Zu25
     &      /cb*ds2-2D0*Zd55*Zu25/cb*ds2*sb**2-EE*MW/SW*Sqrt2
     &      *Zd25*Zu25-EE/MW*Mc**2/SW*Sqrt2*Zd25*Zu25+EE
     &      /MW*Ms**2/SW*Sqrt2*Zd25*Zu25+2D0*EE*MW/SW*Sqrt2
     &      *Zd25*Zu25*sb**2+EE/MW*Mc/SW*Sqrt2*Zd25*Zu55
     &      /cb*hL/sb*xvev-EE/MW*Mc/SW*Sqrt2*Zd25*Zu55/cb
     &      *hL*sb*xvev-EE/MW*Ms/SW*Sqrt2*Zd55*Zu25/cb*hL
     &      *sb*xvev)
      AAABR(144) = Vcd*(2D0*Zd14*Zu55*sb*us2+2D0*Zd44*Zu25
     &      *cb*ds1-EE*MW/SW*Sqrt2*Zd14*Zu25-EE/MW*Mc**2
     &      /SW*Sqrt2*Zd14*Zu25+2D0*EE*MW/SW*Sqrt2*Zd14*Zu25
     &      *sb**2+EE/MW*Mc/SW*Sqrt2*Zd14*Zu55*cb*hL/sb*xvev
     &      )
      AAABR(145) = 6D0*Sqrt2*Zd36*Zd66*Zh13*ds3+6D0*EE/MW
     &      *Mb**2/SW*Zh13/cb-3D0*EE*MW/SW*Zd36**2*Zh13/cb
     &      +3D0*EE*MW/SW*Zd36**2*Zh23*sb-EE/CW**2*MW*SW
     &      *Zd36**2*Zh13/cb+EE/CW**2*MW*SW*Zd36**2*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zd66**2*Zh13/cb+2D0*EE
     &      /CW**2*MW*SW*Zd66**2*Zh23*sb+3D0*EE*MW/SW*Zd36
     &      **2*Zh13/cb*sb**2+EE/CW**2*MW*SW*Zd36**2*Zh13
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd66**2*Zh13/cb
     &      *sb**2-6D0*Mb*Sqrt2*Zd36*Zd66*Zh33/cb*hL*sb-6D0
     &      *EE/MW*Mb/SW*Zd36*Zd66*Zh23/cb*hL*xvev
      AAABR(146) = 6D0*Sqrt2*Zd25*Zd55*Zh11*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh11/cb-3D0*EE*MW/SW*Zd25**2*Zh11/cb
     &      +3D0*EE*MW/SW*Zd25**2*Zh21*sb-EE/CW**2*MW*SW
     &      *Zd25**2*Zh11/cb+EE/CW**2*MW*SW*Zd25**2*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zd55**2*Zh11/cb+2D0*EE
     &      /CW**2*MW*SW*Zd55**2*Zh21*sb+3D0*EE*MW/SW*Zd25
     &      **2*Zh11/cb*sb**2+EE/CW**2*MW*SW*Zd25**2*Zh11
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd55**2*Zh11/cb
     &      *sb**2-6D0*Ms*Sqrt2*Zd25*Zd55*Zh31/cb*hL*sb-6D0
     &      *EE/MW*Ms/SW*Zd25*Zd55*Zh21/cb*hL*xvev
      AAABR(147) = 6D0*Sqrt2*Zd14*Zd44*Zh13*ds1-3D0*EE*MW
     &      /SW*Zd14**2*Zh13*cb+3D0*EE*MW/SW*Zd14**2*Zh23
     &      *sb-EE/CW**2*MW*SW*Zd14**2*Zh13*cb+EE/CW**2*MW
     &      *SW*Zd14**2*Zh23*sb-2D0*EE/CW**2*MW*SW*Zd44**2
     &      *Zh13*cb+2D0*EE/CW**2*MW*SW*Zd44**2*Zh23*sb
      AAABR(148) = 6D0*Sqrt2*Zd14*Zd44*Zh12*ds1-3D0*EE*MW
     &      /SW*Zd14**2*Zh12*cb+3D0*EE*MW/SW*Zd14**2*Zh22
     &      *sb-EE/CW**2*MW*SW*Zd14**2*Zh12*cb+EE/CW**2*MW
     &      *SW*Zd14**2*Zh22*sb-2D0*EE/CW**2*MW*SW*Zd44**2
     &      *Zh12*cb+2D0*EE/CW**2*MW*SW*Zd44**2*Zh22*sb
      AAABR(149) = 6D0*Sqrt2*Zd14*Zd44*Zh11*ds1-3D0*EE*MW
     &      /SW*Zd14**2*Zh11*cb+3D0*EE*MW/SW*Zd14**2*Zh21
     &      *sb-EE/CW**2*MW*SW*Zd14**2*Zh11*cb+EE/CW**2*MW
     &      *SW*Zd14**2*Zh21*sb-2D0*EE/CW**2*MW*SW*Zd44**2
     &      *Zh11*cb+2D0*EE/CW**2*MW*SW*Zd44**2*Zh21*sb
      AAABR(150) = 6D0*Sqrt2*Zd25*Zd55*Zh12*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh12/cb-3D0*EE*MW/SW*Zd25**2*Zh12/cb
     &      +3D0*EE*MW/SW*Zd25**2*Zh22*sb-EE/CW**2*MW*SW
     &      *Zd25**2*Zh12/cb+EE/CW**2*MW*SW*Zd25**2*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zd55**2*Zh12/cb+2D0*EE
     &      /CW**2*MW*SW*Zd55**2*Zh22*sb+3D0*EE*MW/SW*Zd25
     &      **2*Zh12/cb*sb**2+EE/CW**2*MW*SW*Zd25**2*Zh12
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd55**2*Zh12/cb
     &      *sb**2-6D0*Ms*Sqrt2*Zd25*Zd55*Zh32/cb*hL*sb-6D0
     &      *EE/MW*Ms/SW*Zd25*Zd55*Zh22/cb*hL*xvev
      AAABR(151) = 6D0*Sqrt2*Zd25*Zd55*Zh13*ds2+6D0*EE/MW
     &      *Ms**2/SW*Zh13/cb-3D0*EE*MW/SW*Zd25**2*Zh13/cb
     &      +3D0*EE*MW/SW*Zd25**2*Zh23*sb-EE/CW**2*MW*SW
     &      *Zd25**2*Zh13/cb+EE/CW**2*MW*SW*Zd25**2*Zh23
     &      *sb-2D0*EE/CW**2*MW*SW*Zd55**2*Zh13/cb+2D0*EE
     &      /CW**2*MW*SW*Zd55**2*Zh23*sb+3D0*EE*MW/SW*Zd25
     &      **2*Zh13/cb*sb**2+EE/CW**2*MW*SW*Zd25**2*Zh13
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd55**2*Zh13/cb
     &      *sb**2-6D0*Ms*Sqrt2*Zd25*Zd55*Zh33/cb*hL*sb-6D0
     &      *EE/MW*Ms/SW*Zd25*Zd55*Zh23/cb*hL*xvev
      AAABR(152) = 6D0*Sqrt2*Zd36*Zd66*Zh12*ds3+6D0*EE/MW
     &      *Mb**2/SW*Zh12/cb-3D0*EE*MW/SW*Zd36**2*Zh12/cb
     &      +3D0*EE*MW/SW*Zd36**2*Zh22*sb-EE/CW**2*MW*SW
     &      *Zd36**2*Zh12/cb+EE/CW**2*MW*SW*Zd36**2*Zh22
     &      *sb-2D0*EE/CW**2*MW*SW*Zd66**2*Zh12/cb+2D0*EE
     &      /CW**2*MW*SW*Zd66**2*Zh22*sb+3D0*EE*MW/SW*Zd36
     &      **2*Zh12/cb*sb**2+EE/CW**2*MW*SW*Zd36**2*Zh12
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd66**2*Zh12/cb
     &      *sb**2-6D0*Mb*Sqrt2*Zd36*Zd66*Zh32/cb*hL*sb-6D0
     &      *EE/MW*Mb/SW*Zd36*Zd66*Zh22/cb*hL*xvev
      AAABR(153) = 6D0*Sqrt2*Zd36*Zd66*Zh11*ds3+6D0*EE/MW
     &      *Mb**2/SW*Zh11/cb-3D0*EE*MW/SW*Zd36**2*Zh11/cb
     &      +3D0*EE*MW/SW*Zd36**2*Zh21*sb-EE/CW**2*MW*SW
     &      *Zd36**2*Zh11/cb+EE/CW**2*MW*SW*Zd36**2*Zh21
     &      *sb-2D0*EE/CW**2*MW*SW*Zd66**2*Zh11/cb+2D0*EE
     &      /CW**2*MW*SW*Zd66**2*Zh21*sb+3D0*EE*MW/SW*Zd36
     &      **2*Zh11/cb*sb**2+EE/CW**2*MW*SW*Zd36**2*Zh11
     &      /cb*sb**2+2D0*EE/CW**2*MW*SW*Zd66**2*Zh11/cb
     &      *sb**2-6D0*Mb*Sqrt2*Zd36*Zd66*Zh31/cb*hL*sb-6D0
     &      *EE/MW*Mb/SW*Zd36*Zd66*Zh21/cb*hL*xvev
      AAABR(154) = Sqrt2*Za23*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(155) = Sqrt2*Za23*Zu22*Zu55*us2-Sqrt2*Za23*Zu25
     &      *Zu52*us2-Mc*Sqrt2*Za33*Zu22*Zu55*cb*hL/sb+Mc
     &      *Sqrt2*Za33*Zu25*Zu52*cb*hL/sb-EE/MW*Mc/SW*Za13
     &      *Zu22*Zu55*hL/sb*xvev+EE/MW*Mc/SW*Za13*Zu25*Zu52
     &      *hL/sb*xvev
      AAABR(156) = Sqrt2*Za23*Zu33*Zu66*us3-Sqrt2*Za23*Zu36
     &      *Zu63*us3-Mt*Sqrt2*Za33*Zu33*Zu66*cb*hL/sb+Mt
     &      *Sqrt2*Za33*Zu36*Zu63*cb*hL/sb-EE/MW*Mt/SW*Za13
     &      *Zu33*Zu66*hL/sb*xvev+EE/MW*Mt/SW*Za13*Zu36*Zu63
     &      *hL/sb*xvev
      AAABR(157) = Sqrt2*Za21*Zu22*Zu55*us2-Sqrt2*Za21*Zu25
     &      *Zu52*us2-Mc*Sqrt2*Za31*Zu22*Zu55*cb*hL/sb+Mc
     &      *Sqrt2*Za31*Zu25*Zu52*cb*hL/sb-EE/MW*Mc/SW*Za11
     &      *Zu22*Zu55*hL/sb*xvev+EE/MW*Mc/SW*Za11*Zu25*Zu52
     &      *hL/sb*xvev
      AAABR(158) = Sqrt2*Za21*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(159) = Sqrt2*Za22*Zu33*Zu66*us3-Sqrt2*Za22*Zu36
     &      *Zu63*us3-Mt*Sqrt2*Za32*Zu33*Zu66*cb*hL/sb+Mt
     &      *Sqrt2*Za32*Zu36*Zu63*cb*hL/sb-EE/MW*Mt/SW*Za12
     &      *Zu33*Zu66*hL/sb*xvev+EE/MW*Mt/SW*Za12*Zu36*Zu63
     &      *hL/sb*xvev
      AAABR(160) = Sqrt2*Za22*Zu22*Zu55*us2-Sqrt2*Za22*Zu25
     &      *Zu52*us2-Mc*Sqrt2*Za32*Zu22*Zu55*cb*hL/sb+Mc
     &      *Sqrt2*Za32*Zu25*Zu52*cb*hL/sb-EE/MW*Mc/SW*Za12
     &      *Zu22*Zu55*hL/sb*xvev+EE/MW*Mc/SW*Za12*Zu25*Zu52
     &      *hL/sb*xvev
      AAABR(161) = Sqrt2*Za21*Zu33*Zu66*us3-Sqrt2*Za21*Zu36
     &      *Zu63*us3-Mt*Sqrt2*Za31*Zu33*Zu66*cb*hL/sb+Mt
     &      *Sqrt2*Za31*Zu36*Zu63*cb*hL/sb-EE/MW*Mt/SW*Za11
     &      *Zu33*Zu66*hL/sb*xvev+EE/MW*Mt/SW*Za11*Zu36*Zu63
     &      *hL/sb*xvev
      AAABR(162) = Sqrt2*Za22*us1*(Zu11*Zu44-Zu14*Zu41)
      AAABR(163) = 6D0*Sqrt2*Zh21*Zu11*Zu41*us1-3D0*EE*MW
     &      /SW*Zh11*Zu11**2*cb+3D0*EE*MW/SW*Zh21*Zu11**2
     &      *sb+EE/CW**2*MW*SW*Zh11*Zu11**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh11*Zu41**2*cb-EE/CW**2*MW*SW*Zh21
     &      *Zu11**2*sb+4D0*EE/CW**2*MW*SW*Zh21*Zu41**2*sb
      AAABR(164) = 6D0*Sqrt2*Zh21*Zu33*Zu63*us3-6D0*EE/MW
     &      *Mt**2/SW*Zh21/sb-3D0*EE*MW/SW*Zh11*Zu33**2*cb
     &      +3D0*EE*MW/SW*Zh21*Zu33**2*sb+EE/CW**2*MW*SW
     &      *Zh11*Zu33**2*cb-4D0*EE/CW**2*MW*SW*Zh11*Zu63
     &      **2*cb-EE/CW**2*MW*SW*Zh21*Zu33**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh21*Zu63**2*sb+6D0*Mt*Sqrt2*Zh31
     &      *Zu33*Zu63*cb*hL/sb+6D0*EE/MW*Mt/SW*Zh11*Zu33
     &      *Zu63*hL/sb*xvev
      AAABR(165) = 6D0*Sqrt2*Zh23*Zu33*Zu63*us3-6D0*EE/MW
     &      *Mt**2/SW*Zh23/sb-3D0*EE*MW/SW*Zh13*Zu33**2*cb
     &      +3D0*EE*MW/SW*Zh23*Zu33**2*sb+EE/CW**2*MW*SW
     &      *Zh13*Zu33**2*cb-4D0*EE/CW**2*MW*SW*Zh13*Zu63
     &      **2*cb-EE/CW**2*MW*SW*Zh23*Zu33**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh23*Zu63**2*sb+6D0*Mt*Sqrt2*Zh33
     &      *Zu33*Zu63*cb*hL/sb+6D0*EE/MW*Mt/SW*Zh13*Zu33
     &      *Zu63*hL/sb*xvev
      AAABR(166) = 6D0*Sqrt2*Zh22*Zu11*Zu41*us1-3D0*EE*MW
     &      /SW*Zh12*Zu11**2*cb+3D0*EE*MW/SW*Zh22*Zu11**2
     &      *sb+EE/CW**2*MW*SW*Zh12*Zu11**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh12*Zu41**2*cb-EE/CW**2*MW*SW*Zh22
     &      *Zu11**2*sb+4D0*EE/CW**2*MW*SW*Zh22*Zu41**2*sb
      AAABR(167) = 6D0*Sqrt2*Zh22*Zu33*Zu63*us3-6D0*EE/MW
     &      *Mt**2/SW*Zh22/sb-3D0*EE*MW/SW*Zh12*Zu33**2*cb
     &      +3D0*EE*MW/SW*Zh22*Zu33**2*sb+EE/CW**2*MW*SW
     &      *Zh12*Zu33**2*cb-4D0*EE/CW**2*MW*SW*Zh12*Zu63
     &      **2*cb-EE/CW**2*MW*SW*Zh22*Zu33**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh22*Zu63**2*sb+6D0*Mt*Sqrt2*Zh32
     &      *Zu33*Zu63*cb*hL/sb+6D0*EE/MW*Mt/SW*Zh12*Zu33
     &      *Zu63*hL/sb*xvev
      AAABR(168) = 6D0*Sqrt2*Zh23*Zu22*Zu52*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh23/sb-3D0*EE*MW/SW*Zh13*Zu22**2*cb
     &      +3D0*EE*MW/SW*Zh23*Zu22**2*sb+EE/CW**2*MW*SW
     &      *Zh13*Zu22**2*cb-4D0*EE/CW**2*MW*SW*Zh13*Zu52
     &      **2*cb-EE/CW**2*MW*SW*Zh23*Zu22**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh23*Zu52**2*sb+6D0*Mc*Sqrt2*Zh33
     &      *Zu22*Zu52*cb*hL/sb+6D0*EE/MW*Mc/SW*Zh13*Zu22
     &      *Zu52*hL/sb*xvev
      AAABR(169) = 6D0*Sqrt2*Zh21*Zu22*Zu52*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh21/sb-3D0*EE*MW/SW*Zh11*Zu22**2*cb
     &      +3D0*EE*MW/SW*Zh21*Zu22**2*sb+EE/CW**2*MW*SW
     &      *Zh11*Zu22**2*cb-4D0*EE/CW**2*MW*SW*Zh11*Zu52
     &      **2*cb-EE/CW**2*MW*SW*Zh21*Zu22**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh21*Zu52**2*sb+6D0*Mc*Sqrt2*Zh31
     &      *Zu22*Zu52*cb*hL/sb+6D0*EE/MW*Mc/SW*Zh11*Zu22
     &      *Zu52*hL/sb*xvev
      AAABR(170) = 6D0*Sqrt2*Zh23*Zu11*Zu41*us1-3D0*EE*MW
     &      /SW*Zh13*Zu11**2*cb+3D0*EE*MW/SW*Zh23*Zu11**2
     &      *sb+EE/CW**2*MW*SW*Zh13*Zu11**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh13*Zu41**2*cb-EE/CW**2*MW*SW*Zh23
     &      *Zu11**2*sb+4D0*EE/CW**2*MW*SW*Zh23*Zu41**2*sb
      AAABR(171) = 6D0*Sqrt2*Zh22*Zu22*Zu52*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh22/sb-3D0*EE*MW/SW*Zh12*Zu22**2*cb
     &      +3D0*EE*MW/SW*Zh22*Zu22**2*sb+EE/CW**2*MW*SW
     &      *Zh12*Zu22**2*cb-4D0*EE/CW**2*MW*SW*Zh12*Zu52
     &      **2*cb-EE/CW**2*MW*SW*Zh22*Zu22**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh22*Zu52**2*sb+6D0*Mc*Sqrt2*Zh32
     &      *Zu22*Zu52*cb*hL/sb+6D0*EE/MW*Mc/SW*Zh12*Zu22
     &      *Zu52*hL/sb*xvev
      AAABR(172) = 3D0*Sqrt2*Zh23*Zu33*Zu66*us3+3D0*Sqrt2
     &      *Zh23*Zu36*Zu63*us3-3D0*EE*MW/SW*Zh13*Zu33*Zu36
     &      *cb+3D0*EE*MW/SW*Zh23*Zu33*Zu36*sb+EE/CW**2*MW
     &      *SW*Zh13*Zu33*Zu36*cb-4D0*EE/CW**2*MW*SW*Zh13
     &      *Zu63*Zu66*cb-EE/CW**2*MW*SW*Zh23*Zu33*Zu36*sb
     &      +4D0*EE/CW**2*MW*SW*Zh23*Zu63*Zu66*sb+3D0*Mt
     &      *Sqrt2*Zh33*Zu33*Zu66*cb*hL/sb+3D0*Mt*Sqrt2*Zh33
     &      *Zu36*Zu63*cb*hL/sb+3D0*EE/MW*Mt/SW*Zh13*Zu33
     &      *Zu66*hL/sb*xvev+3D0*EE/MW*Mt/SW*Zh13*Zu36*Zu63
     &      *hL/sb*xvev
      AAABR(173) = 3D0*Sqrt2*Zh21*Zu11*Zu44*us1+3D0*Sqrt2
     &      *Zh21*Zu14*Zu41*us1-3D0*EE*MW/SW*Zh11*Zu11*Zu14
     &      *cb+3D0*EE*MW/SW*Zh21*Zu11*Zu14*sb+EE/CW**2*MW
     &      *SW*Zh11*Zu11*Zu14*cb-4D0*EE/CW**2*MW*SW*Zh11
     &      *Zu41*Zu44*cb-EE/CW**2*MW*SW*Zh21*Zu11*Zu14*sb
     &      +4D0*EE/CW**2*MW*SW*Zh21*Zu41*Zu44*sb
      AAABR(174) = 3D0*Sqrt2*Zh23*Zu11*Zu44*us1+3D0*Sqrt2
     &      *Zh23*Zu14*Zu41*us1-3D0*EE*MW/SW*Zh13*Zu11*Zu14
     &      *cb+3D0*EE*MW/SW*Zh23*Zu11*Zu14*sb+EE/CW**2*MW
     &      *SW*Zh13*Zu11*Zu14*cb-4D0*EE/CW**2*MW*SW*Zh13
     &      *Zu41*Zu44*cb-EE/CW**2*MW*SW*Zh23*Zu11*Zu14*sb
     &      +4D0*EE/CW**2*MW*SW*Zh23*Zu41*Zu44*sb
      AAABR(175) = 3D0*Sqrt2*Zh22*Zu22*Zu55*us2+3D0*Sqrt2
     &      *Zh22*Zu25*Zu52*us2-3D0*EE*MW/SW*Zh12*Zu22*Zu25
     &      *cb+3D0*EE*MW/SW*Zh22*Zu22*Zu25*sb+EE/CW**2*MW
     &      *SW*Zh12*Zu22*Zu25*cb-4D0*EE/CW**2*MW*SW*Zh12
     &      *Zu52*Zu55*cb-EE/CW**2*MW*SW*Zh22*Zu22*Zu25*sb
     &      +4D0*EE/CW**2*MW*SW*Zh22*Zu52*Zu55*sb+3D0*Mc
     &      *Sqrt2*Zh32*Zu22*Zu55*cb*hL/sb+3D0*Mc*Sqrt2*Zh32
     &      *Zu25*Zu52*cb*hL/sb+3D0*EE/MW*Mc/SW*Zh12*Zu22
     &      *Zu55*hL/sb*xvev+3D0*EE/MW*Mc/SW*Zh12*Zu25*Zu52
     &      *hL/sb*xvev
      AAABR(176) = 3D0*Sqrt2*Zh22*Zu11*Zu44*us1+3D0*Sqrt2
     &      *Zh22*Zu14*Zu41*us1-3D0*EE*MW/SW*Zh12*Zu11*Zu14
     &      *cb+3D0*EE*MW/SW*Zh22*Zu11*Zu14*sb+EE/CW**2*MW
     &      *SW*Zh12*Zu11*Zu14*cb-4D0*EE/CW**2*MW*SW*Zh12
     &      *Zu41*Zu44*cb-EE/CW**2*MW*SW*Zh22*Zu11*Zu14*sb
     &      +4D0*EE/CW**2*MW*SW*Zh22*Zu41*Zu44*sb
      AAABR(177) = 3D0*Sqrt2*Zh21*Zu33*Zu66*us3+3D0*Sqrt2
     &      *Zh21*Zu36*Zu63*us3-3D0*EE*MW/SW*Zh11*Zu33*Zu36
     &      *cb+3D0*EE*MW/SW*Zh21*Zu33*Zu36*sb+EE/CW**2*MW
     &      *SW*Zh11*Zu33*Zu36*cb-4D0*EE/CW**2*MW*SW*Zh11
     &      *Zu63*Zu66*cb-EE/CW**2*MW*SW*Zh21*Zu33*Zu36*sb
     &      +4D0*EE/CW**2*MW*SW*Zh21*Zu63*Zu66*sb+3D0*Mt
     &      *Sqrt2*Zh31*Zu33*Zu66*cb*hL/sb+3D0*Mt*Sqrt2*Zh31
     &      *Zu36*Zu63*cb*hL/sb+3D0*EE/MW*Mt/SW*Zh11*Zu33
     &      *Zu66*hL/sb*xvev+3D0*EE/MW*Mt/SW*Zh11*Zu36*Zu63
     &      *hL/sb*xvev
      AAABR(178) = 3D0*Sqrt2*Zh21*Zu22*Zu55*us2+3D0*Sqrt2
     &      *Zh21*Zu25*Zu52*us2-3D0*EE*MW/SW*Zh11*Zu22*Zu25
     &      *cb+3D0*EE*MW/SW*Zh21*Zu22*Zu25*sb+EE/CW**2*MW
     &      *SW*Zh11*Zu22*Zu25*cb-4D0*EE/CW**2*MW*SW*Zh11
     &      *Zu52*Zu55*cb-EE/CW**2*MW*SW*Zh21*Zu22*Zu25*sb
     &      +4D0*EE/CW**2*MW*SW*Zh21*Zu52*Zu55*sb+3D0*Mc
     &      *Sqrt2*Zh31*Zu22*Zu55*cb*hL/sb+3D0*Mc*Sqrt2*Zh31
     &      *Zu25*Zu52*cb*hL/sb+3D0*EE/MW*Mc/SW*Zh11*Zu22
     &      *Zu55*hL/sb*xvev+3D0*EE/MW*Mc/SW*Zh11*Zu25*Zu52
     &      *hL/sb*xvev
      AAABR(179) = 3D0*Sqrt2*Zh22*Zu33*Zu66*us3+3D0*Sqrt2
     &      *Zh22*Zu36*Zu63*us3-3D0*EE*MW/SW*Zh12*Zu33*Zu36
     &      *cb+3D0*EE*MW/SW*Zh22*Zu33*Zu36*sb+EE/CW**2*MW
     &      *SW*Zh12*Zu33*Zu36*cb-4D0*EE/CW**2*MW*SW*Zh12
     &      *Zu63*Zu66*cb-EE/CW**2*MW*SW*Zh22*Zu33*Zu36*sb
     &      +4D0*EE/CW**2*MW*SW*Zh22*Zu63*Zu66*sb+3D0*Mt
     &      *Sqrt2*Zh32*Zu33*Zu66*cb*hL/sb+3D0*Mt*Sqrt2*Zh32
     &      *Zu36*Zu63*cb*hL/sb+3D0*EE/MW*Mt/SW*Zh12*Zu33
     &      *Zu66*hL/sb*xvev+3D0*EE/MW*Mt/SW*Zh12*Zu36*Zu63
     &      *hL/sb*xvev
      AAABR(180) = 3D0*Sqrt2*Zh23*Zu22*Zu55*us2+3D0*Sqrt2
     &      *Zh23*Zu25*Zu52*us2-3D0*EE*MW/SW*Zh13*Zu22*Zu25
     &      *cb+3D0*EE*MW/SW*Zh23*Zu22*Zu25*sb+EE/CW**2*MW
     &      *SW*Zh13*Zu22*Zu25*cb-4D0*EE/CW**2*MW*SW*Zh13
     &      *Zu52*Zu55*cb-EE/CW**2*MW*SW*Zh23*Zu22*Zu25*sb
     &      +4D0*EE/CW**2*MW*SW*Zh23*Zu52*Zu55*sb+3D0*Mc
     &      *Sqrt2*Zh33*Zu22*Zu55*cb*hL/sb+3D0*Mc*Sqrt2*Zh33
     &      *Zu25*Zu52*cb*hL/sb+3D0*EE/MW*Mc/SW*Zh13*Zu22
     &      *Zu55*hL/sb*xvev+3D0*EE/MW*Mc/SW*Zh13*Zu25*Zu52
     &      *hL/sb*xvev
      AAABR(181) = 6D0*Sqrt2*Zh22*Zu14*Zu44*us1-3D0*EE*MW
     &      /SW*Zh12*Zu14**2*cb+3D0*EE*MW/SW*Zh22*Zu14**2
     &      *sb+EE/CW**2*MW*SW*Zh12*Zu14**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh12*Zu44**2*cb-EE/CW**2*MW*SW*Zh22
     &      *Zu14**2*sb+4D0*EE/CW**2*MW*SW*Zh22*Zu44**2*sb
      AAABR(182) = 6D0*Sqrt2*Zh21*Zu25*Zu55*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh21/sb-3D0*EE*MW/SW*Zh11*Zu25**2*cb
     &      +3D0*EE*MW/SW*Zh21*Zu25**2*sb+EE/CW**2*MW*SW
     &      *Zh11*Zu25**2*cb-4D0*EE/CW**2*MW*SW*Zh11*Zu55
     &      **2*cb-EE/CW**2*MW*SW*Zh21*Zu25**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh21*Zu55**2*sb+6D0*Mc*Sqrt2*Zh31
     &      *Zu25*Zu55*cb*hL/sb+6D0*EE/MW*Mc/SW*Zh11*Zu25
     &      *Zu55*hL/sb*xvev
      AAABR(183) = 6D0*Sqrt2*Zh21*Zu14*Zu44*us1-3D0*EE*MW
     &      /SW*Zh11*Zu14**2*cb+3D0*EE*MW/SW*Zh21*Zu14**2
     &      *sb+EE/CW**2*MW*SW*Zh11*Zu14**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh11*Zu44**2*cb-EE/CW**2*MW*SW*Zh21
     &      *Zu14**2*sb+4D0*EE/CW**2*MW*SW*Zh21*Zu44**2*sb
      AAABR(184) = 6D0*Sqrt2*Zh23*Zu25*Zu55*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh23/sb-3D0*EE*MW/SW*Zh13*Zu25**2*cb
     &      +3D0*EE*MW/SW*Zh23*Zu25**2*sb+EE/CW**2*MW*SW
     &      *Zh13*Zu25**2*cb-4D0*EE/CW**2*MW*SW*Zh13*Zu55
     &      **2*cb-EE/CW**2*MW*SW*Zh23*Zu25**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh23*Zu55**2*sb+6D0*Mc*Sqrt2*Zh33
     &      *Zu25*Zu55*cb*hL/sb+6D0*EE/MW*Mc/SW*Zh13*Zu25
     &      *Zu55*hL/sb*xvev
      AAABR(185) = 6D0*Sqrt2*Zh23*Zu14*Zu44*us1-3D0*EE*MW
     &      /SW*Zh13*Zu14**2*cb+3D0*EE*MW/SW*Zh23*Zu14**2
     &      *sb+EE/CW**2*MW*SW*Zh13*Zu14**2*cb-4D0*EE/CW
     &      **2*MW*SW*Zh13*Zu44**2*cb-EE/CW**2*MW*SW*Zh23
     &      *Zu14**2*sb+4D0*EE/CW**2*MW*SW*Zh23*Zu44**2*sb
      AAABR(186) = 6D0*Sqrt2*Zh22*Zu36*Zu66*us3-6D0*EE/MW
     &      *Mt**2/SW*Zh22/sb-3D0*EE*MW/SW*Zh12*Zu36**2*cb
     &      +3D0*EE*MW/SW*Zh22*Zu36**2*sb+EE/CW**2*MW*SW
     &      *Zh12*Zu36**2*cb-4D0*EE/CW**2*MW*SW*Zh12*Zu66
     &      **2*cb-EE/CW**2*MW*SW*Zh22*Zu36**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh22*Zu66**2*sb+6D0*Mt*Sqrt2*Zh32
     &      *Zu36*Zu66*cb*hL/sb+6D0*EE/MW*Mt/SW*Zh12*Zu36
     &      *Zu66*hL/sb*xvev
      AAABR(187) = 6D0*Sqrt2*Zh21*Zu36*Zu66*us3-6D0*EE/MW
     &      *Mt**2/SW*Zh21/sb-3D0*EE*MW/SW*Zh11*Zu36**2*cb
     &      +3D0*EE*MW/SW*Zh21*Zu36**2*sb+EE/CW**2*MW*SW
     &      *Zh11*Zu36**2*cb-4D0*EE/CW**2*MW*SW*Zh11*Zu66
     &      **2*cb-EE/CW**2*MW*SW*Zh21*Zu36**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh21*Zu66**2*sb+6D0*Mt*Sqrt2*Zh31
     &      *Zu36*Zu66*cb*hL/sb+6D0*EE/MW*Mt/SW*Zh11*Zu36
     &      *Zu66*hL/sb*xvev
      AAABR(188) = 6D0*Sqrt2*Zh22*Zu25*Zu55*us2-6D0*EE/MW
     &      *Mc**2/SW*Zh22/sb-3D0*EE*MW/SW*Zh12*Zu25**2*cb
     &      +3D0*EE*MW/SW*Zh22*Zu25**2*sb+EE/CW**2*MW*SW
     &      *Zh12*Zu25**2*cb-4D0*EE/CW**2*MW*SW*Zh12*Zu55
     &      **2*cb-EE/CW**2*MW*SW*Zh22*Zu25**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh22*Zu55**2*sb+6D0*Mc*Sqrt2*Zh32
     &      *Zu25*Zu55*cb*hL/sb+6D0*EE/MW*Mc/SW*Zh12*Zu25
     &      *Zu55*hL/sb*xvev
      AAABR(189) = 6D0*Sqrt2*Zh23*Zu36*Zu66*us3-6D0*EE/MW
     &      *Mt**2/SW*Zh23/sb-3D0*EE*MW/SW*Zh13*Zu36**2*cb
     &      +3D0*EE*MW/SW*Zh23*Zu36**2*sb+EE/CW**2*MW*SW
     &      *Zh13*Zu36**2*cb-4D0*EE/CW**2*MW*SW*Zh13*Zu66
     &      **2*cb-EE/CW**2*MW*SW*Zh23*Zu36**2*sb+4D0*EE
     &      /CW**2*MW*SW*Zh23*Zu66**2*sb+6D0*Mt*Sqrt2*Zh33
     &      *Zu36*Zu66*cb*hL/sb+6D0*EE/MW*Mt/SW*Zh13*Zu36
     &      *Zu66*hL/sb*xvev
      AAABR(190) = 2D0*Sqrt2*Zh32*hL**2*xvev+EE*MW/SW*Zh12
     &      *cb-EE*MW/SW*Zh22*sb+2D0*EE*MW/SW*Zh22*sb**3
     &      +EE/CW**2*MW*SW*Zh12*cb-EE/CW**2*MW*SW*Zh22*sb
     &      +2D0*EE/CW**2*MW*SW*Zh22*sb**3+4D0*MW/EE*SW*Zh22
     &      *hL**2*sb-4D0*MW/EE*SW*Zh22*hL**2*sb**3-2D0*EE
     &      *MW/SW*Zh12*cb*sb**2-2D0*Sqrt2*Zh32*cb*hL*hLs
     &      *sb-2D0*EE/CW**2*MW*SW*Zh12*cb*sb**2+4D0*MW/EE
     &      *SW*Zh12*cb*hL**2*sb**2-4D0*Sqrt2*Zh32*cb*hK
     &      *hL*sb*xvev
      AAABR(191) = 2D0*Sqrt2*Zh33*hL**2*xvev+EE*MW/SW*Zh13
     &      *cb-EE*MW/SW*Zh23*sb+2D0*EE*MW/SW*Zh23*sb**3
     &      +EE/CW**2*MW*SW*Zh13*cb-EE/CW**2*MW*SW*Zh23*sb
     &      +2D0*EE/CW**2*MW*SW*Zh23*sb**3+4D0*MW/EE*SW*Zh23
     &      *hL**2*sb-4D0*MW/EE*SW*Zh23*hL**2*sb**3-2D0*EE
     &      *MW/SW*Zh13*cb*sb**2-2D0*Sqrt2*Zh33*cb*hL*hLs
     &      *sb-2D0*EE/CW**2*MW*SW*Zh13*cb*sb**2+4D0*MW/EE
     &      *SW*Zh13*cb*hL**2*sb**2-4D0*Sqrt2*Zh33*cb*hK
     &      *hL*sb*xvev
      AAABR(192) = 2D0*Sqrt2*Zh31*hL**2*xvev+EE*MW/SW*Zh11
     &      *cb-EE*MW/SW*Zh21*sb+2D0*EE*MW/SW*Zh21*sb**3
     &      +EE/CW**2*MW*SW*Zh11*cb-EE/CW**2*MW*SW*Zh21*sb
     &      +2D0*EE/CW**2*MW*SW*Zh21*sb**3+4D0*MW/EE*SW*Zh21
     &      *hL**2*sb-4D0*MW/EE*SW*Zh21*hL**2*sb**3-2D0*EE
     &      *MW/SW*Zh11*cb*sb**2-2D0*Sqrt2*Zh31*cb*hL*hLs
     &      *sb-2D0*EE/CW**2*MW*SW*Zh11*cb*sb**2+4D0*MW/EE
     &      *SW*Zh11*cb*hL**2*sb**2-4D0*Sqrt2*Zh31*cb*hK
     &      *hL*sb*xvev
      AAABR(193) = 2D0*Sqrt2*Za13**2*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Za23**2*Zh33*hL**2*xvev-2D0*Sqrt2*Za33
     &      **2*Zh33*hK*hKs+4D0*Sqrt2*Za33**2*Zh33*hK**2
     &      *xvev+EE*MW/SW*Za13**2*Zh13*cb-EE*MW/SW*Za13
     &      **2*Zh23*sb-EE*MW/SW*Za23**2*Zh13*cb+EE*MW/SW
     &      *Za23**2*Zh23*sb+2D0*Sqrt2*Za13*Za23*Zh33*hL
     &      *hLs+2D0*Sqrt2*Za13*Za33*Zh23*hL*hLs+2D0*Sqrt2
     &      *Za23*Za33*Zh13*hL*hLs+EE/CW**2*MW*SW*Za13**2
     &      *Zh13*cb-EE/CW**2*MW*SW*Za13**2*Zh23*sb-EE/CW
     &      **2*MW*SW*Za23**2*Zh13*cb+EE/CW**2*MW*SW*Za23
     &      **2*Zh23*sb+4D0*MW/EE*SW*Za13**2*Zh23*hL**2*sb
     &      +4D0*MW/EE*SW*Za23**2*Zh13*cb*hL**2+4D0*MW/EE
     &      *SW*Za33**2*Zh13*cb*hL**2+4D0*MW/EE*SW*Za33**2
     &      *Zh23*hL**2*sb+4D0*Sqrt2*Za13*Za23*Zh33*hK*hL
     &      *xvev-4D0*Sqrt2*Za13*Za33*Zh23*hK*hL*xvev-4D0
     &      *Sqrt2*Za23*Za33*Zh13*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za33**2*Zh13*hK*hL*sb+4D0*MW/EE*SW*Za33**2*Zh23
     &      *cb*hK*hL-8D0*MW/EE*SW*Za13*Za33*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Za23*Za33*Zh33*cb*hK*hL
      AAABR(194) = 2D0*Sqrt2*Za13**2*Zh32*hL**2*xvev+2D0
     &      *Sqrt2*Za23**2*Zh32*hL**2*xvev-2D0*Sqrt2*Za33
     &      **2*Zh32*hK*hKs+4D0*Sqrt2*Za33**2*Zh32*hK**2
     &      *xvev+EE*MW/SW*Za13**2*Zh12*cb-EE*MW/SW*Za13
     &      **2*Zh22*sb-EE*MW/SW*Za23**2*Zh12*cb+EE*MW/SW
     &      *Za23**2*Zh22*sb+2D0*Sqrt2*Za13*Za23*Zh32*hL
     &      *hLs+2D0*Sqrt2*Za13*Za33*Zh22*hL*hLs+2D0*Sqrt2
     &      *Za23*Za33*Zh12*hL*hLs+EE/CW**2*MW*SW*Za13**2
     &      *Zh12*cb-EE/CW**2*MW*SW*Za13**2*Zh22*sb-EE/CW
     &      **2*MW*SW*Za23**2*Zh12*cb+EE/CW**2*MW*SW*Za23
     &      **2*Zh22*sb+4D0*MW/EE*SW*Za13**2*Zh22*hL**2*sb
     &      +4D0*MW/EE*SW*Za23**2*Zh12*cb*hL**2+4D0*MW/EE
     &      *SW*Za33**2*Zh12*cb*hL**2+4D0*MW/EE*SW*Za33**2
     &      *Zh22*hL**2*sb+4D0*Sqrt2*Za13*Za23*Zh32*hK*hL
     &      *xvev-4D0*Sqrt2*Za13*Za33*Zh22*hK*hL*xvev-4D0
     &      *Sqrt2*Za23*Za33*Zh12*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za33**2*Zh12*hK*hL*sb+4D0*MW/EE*SW*Za33**2*Zh22
     &      *cb*hK*hL-8D0*MW/EE*SW*Za13*Za33*Zh32*hK*hL*sb
     &      -8D0*MW/EE*SW*Za23*Za33*Zh32*cb*hK*hL
      AAABR(195) = 2D0*Sqrt2*Za13**2*Zh31*hL**2*xvev+2D0
     &      *Sqrt2*Za23**2*Zh31*hL**2*xvev-2D0*Sqrt2*Za33
     &      **2*Zh31*hK*hKs+4D0*Sqrt2*Za33**2*Zh31*hK**2
     &      *xvev+EE*MW/SW*Za13**2*Zh11*cb-EE*MW/SW*Za13
     &      **2*Zh21*sb-EE*MW/SW*Za23**2*Zh11*cb+EE*MW/SW
     &      *Za23**2*Zh21*sb+2D0*Sqrt2*Za13*Za23*Zh31*hL
     &      *hLs+2D0*Sqrt2*Za13*Za33*Zh21*hL*hLs+2D0*Sqrt2
     &      *Za23*Za33*Zh11*hL*hLs+EE/CW**2*MW*SW*Za13**2
     &      *Zh11*cb-EE/CW**2*MW*SW*Za13**2*Zh21*sb-EE/CW
     &      **2*MW*SW*Za23**2*Zh11*cb+EE/CW**2*MW*SW*Za23
     &      **2*Zh21*sb+4D0*MW/EE*SW*Za13**2*Zh21*hL**2*sb
     &      +4D0*MW/EE*SW*Za23**2*Zh11*cb*hL**2+4D0*MW/EE
     &      *SW*Za33**2*Zh11*cb*hL**2+4D0*MW/EE*SW*Za33**2
     &      *Zh21*hL**2*sb+4D0*Sqrt2*Za13*Za23*Zh31*hK*hL
     &      *xvev-4D0*Sqrt2*Za13*Za33*Zh21*hK*hL*xvev-4D0
     &      *Sqrt2*Za23*Za33*Zh11*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za33**2*Zh11*hK*hL*sb+4D0*MW/EE*SW*Za33**2*Zh21
     &      *cb*hK*hL-8D0*MW/EE*SW*Za13*Za33*Zh31*hK*hL*sb
     &      -8D0*MW/EE*SW*Za23*Za33*Zh31*cb*hK*hL
      AAABR(196) = 2D0*Sqrt2*Za12*Za13*Zh32*hL**2*xvev+Sqrt2
     &      *Za12*Za23*Zh32*hL*hLs+Sqrt2*Za12*Za33*Zh22*hL
     &      *hLs+Sqrt2*Za13*Za22*Zh32*hL*hLs+Sqrt2*Za13*Za32
     &      *Zh22*hL*hLs+2D0*Sqrt2*Za22*Za23*Zh32*hL**2*xvev
     &      +Sqrt2*Za22*Za33*Zh12*hL*hLs+Sqrt2*Za23*Za32
     &      *Zh12*hL*hLs-2D0*Sqrt2*Za32*Za33*Zh32*hK*hKs
     &      +4D0*Sqrt2*Za32*Za33*Zh32*hK**2*xvev+EE*MW/SW
     &      *Za12*Za13*Zh12*cb-EE*MW/SW*Za12*Za13*Zh22*sb
     &      -EE*MW/SW*Za22*Za23*Zh12*cb+EE*MW/SW*Za22*Za23
     &      *Zh22*sb+2D0*Sqrt2*Za12*Za23*Zh32*hK*hL*xvev
     &      -2D0*Sqrt2*Za12*Za33*Zh22*hK*hL*xvev+2D0*Sqrt2
     &      *Za13*Za22*Zh32*hK*hL*xvev-2D0*Sqrt2*Za13*Za32
     &      *Zh22*hK*hL*xvev-2D0*Sqrt2*Za22*Za33*Zh12*hK
     &      *hL*xvev-2D0*Sqrt2*Za23*Za32*Zh12*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za12*Za13*Zh12*cb-EE/CW**2*MW
     &      *SW*Za12*Za13*Zh22*sb-EE/CW**2*MW*SW*Za22*Za23
     &      *Zh12*cb+EE/CW**2*MW*SW*Za22*Za23*Zh22*sb+4D0
     &      *MW/EE*SW*Za12*Za13*Zh22*hL**2*sb+4D0*MW/EE*SW
     &      *Za22*Za23*Zh12*cb*hL**2+4D0*MW/EE*SW*Za32*Za33
     &      *Zh12*cb*hL**2+4D0*MW/EE*SW*Za32*Za33*Zh22*hL
     &      **2*sb-4D0*MW/EE*SW*Za12*Za33*Zh32*hK*hL*sb-4D0
     &      *MW/EE*SW*Za13*Za32*Zh32*hK*hL*sb-4D0*MW/EE*SW
     &      *Za22*Za33*Zh32*cb*hK*hL-4D0*MW/EE*SW*Za23*Za32
     &      *Zh32*cb*hK*hL+4D0*MW/EE*SW*Za32*Za33*Zh12*hK
     &      *hL*sb+4D0*MW/EE*SW*Za32*Za33*Zh22*cb*hK*hL
      AAABR(197) = 2D0*Sqrt2*Za11*Za13*Zh31*hL**2*xvev+Sqrt2
     &      *Za11*Za23*Zh31*hL*hLs+Sqrt2*Za11*Za33*Zh21*hL
     &      *hLs+Sqrt2*Za13*Za21*Zh31*hL*hLs+Sqrt2*Za13*Za31
     &      *Zh21*hL*hLs+2D0*Sqrt2*Za21*Za23*Zh31*hL**2*xvev
     &      +Sqrt2*Za21*Za33*Zh11*hL*hLs+Sqrt2*Za23*Za31
     &      *Zh11*hL*hLs-2D0*Sqrt2*Za31*Za33*Zh31*hK*hKs
     &      +4D0*Sqrt2*Za31*Za33*Zh31*hK**2*xvev+EE*MW/SW
     &      *Za11*Za13*Zh11*cb-EE*MW/SW*Za11*Za13*Zh21*sb
     &      -EE*MW/SW*Za21*Za23*Zh11*cb+EE*MW/SW*Za21*Za23
     &      *Zh21*sb+2D0*Sqrt2*Za11*Za23*Zh31*hK*hL*xvev
     &      -2D0*Sqrt2*Za11*Za33*Zh21*hK*hL*xvev+2D0*Sqrt2
     &      *Za13*Za21*Zh31*hK*hL*xvev-2D0*Sqrt2*Za13*Za31
     &      *Zh21*hK*hL*xvev-2D0*Sqrt2*Za21*Za33*Zh11*hK
     &      *hL*xvev-2D0*Sqrt2*Za23*Za31*Zh11*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za11*Za13*Zh11*cb-EE/CW**2*MW
     &      *SW*Za11*Za13*Zh21*sb-EE/CW**2*MW*SW*Za21*Za23
     &      *Zh11*cb+EE/CW**2*MW*SW*Za21*Za23*Zh21*sb+4D0
     &      *MW/EE*SW*Za11*Za13*Zh21*hL**2*sb+4D0*MW/EE*SW
     &      *Za21*Za23*Zh11*cb*hL**2+4D0*MW/EE*SW*Za31*Za33
     &      *Zh11*cb*hL**2+4D0*MW/EE*SW*Za31*Za33*Zh21*hL
     &      **2*sb-4D0*MW/EE*SW*Za11*Za33*Zh31*hK*hL*sb-4D0
     &      *MW/EE*SW*Za13*Za31*Zh31*hK*hL*sb-4D0*MW/EE*SW
     &      *Za21*Za33*Zh31*cb*hK*hL-4D0*MW/EE*SW*Za23*Za31
     &      *Zh31*cb*hK*hL+4D0*MW/EE*SW*Za31*Za33*Zh11*hK
     &      *hL*sb+4D0*MW/EE*SW*Za31*Za33*Zh21*cb*hK*hL
      AAABR(198) = 2D0*Sqrt2*Za12*Za13*Zh33*hL**2*xvev+Sqrt2
     &      *Za12*Za23*Zh33*hL*hLs+Sqrt2*Za12*Za33*Zh23*hL
     &      *hLs+Sqrt2*Za13*Za22*Zh33*hL*hLs+Sqrt2*Za13*Za32
     &      *Zh23*hL*hLs+2D0*Sqrt2*Za22*Za23*Zh33*hL**2*xvev
     &      +Sqrt2*Za22*Za33*Zh13*hL*hLs+Sqrt2*Za23*Za32
     &      *Zh13*hL*hLs-2D0*Sqrt2*Za32*Za33*Zh33*hK*hKs
     &      +4D0*Sqrt2*Za32*Za33*Zh33*hK**2*xvev+EE*MW/SW
     &      *Za12*Za13*Zh13*cb-EE*MW/SW*Za12*Za13*Zh23*sb
     &      -EE*MW/SW*Za22*Za23*Zh13*cb+EE*MW/SW*Za22*Za23
     &      *Zh23*sb+2D0*Sqrt2*Za12*Za23*Zh33*hK*hL*xvev
     &      -2D0*Sqrt2*Za12*Za33*Zh23*hK*hL*xvev+2D0*Sqrt2
     &      *Za13*Za22*Zh33*hK*hL*xvev-2D0*Sqrt2*Za13*Za32
     &      *Zh23*hK*hL*xvev-2D0*Sqrt2*Za22*Za33*Zh13*hK
     &      *hL*xvev-2D0*Sqrt2*Za23*Za32*Zh13*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za12*Za13*Zh13*cb-EE/CW**2*MW
     &      *SW*Za12*Za13*Zh23*sb-EE/CW**2*MW*SW*Za22*Za23
     &      *Zh13*cb+EE/CW**2*MW*SW*Za22*Za23*Zh23*sb+4D0
     &      *MW/EE*SW*Za12*Za13*Zh23*hL**2*sb+4D0*MW/EE*SW
     &      *Za22*Za23*Zh13*cb*hL**2+4D0*MW/EE*SW*Za32*Za33
     &      *Zh13*cb*hL**2+4D0*MW/EE*SW*Za32*Za33*Zh23*hL
     &      **2*sb-4D0*MW/EE*SW*Za12*Za33*Zh33*hK*hL*sb-4D0
     &      *MW/EE*SW*Za13*Za32*Zh33*hK*hL*sb-4D0*MW/EE*SW
     &      *Za22*Za33*Zh33*cb*hK*hL-4D0*MW/EE*SW*Za23*Za32
     &      *Zh33*cb*hK*hL+4D0*MW/EE*SW*Za32*Za33*Zh13*hK
     &      *hL*sb+4D0*MW/EE*SW*Za32*Za33*Zh23*cb*hK*hL
      AAABR(199) = 2D0*Sqrt2*Za11*Za13*Zh32*hL**2*xvev+Sqrt2
     &      *Za11*Za23*Zh32*hL*hLs+Sqrt2*Za11*Za33*Zh22*hL
     &      *hLs+Sqrt2*Za13*Za21*Zh32*hL*hLs+Sqrt2*Za13*Za31
     &      *Zh22*hL*hLs+2D0*Sqrt2*Za21*Za23*Zh32*hL**2*xvev
     &      +Sqrt2*Za21*Za33*Zh12*hL*hLs+Sqrt2*Za23*Za31
     &      *Zh12*hL*hLs-2D0*Sqrt2*Za31*Za33*Zh32*hK*hKs
     &      +4D0*Sqrt2*Za31*Za33*Zh32*hK**2*xvev+EE*MW/SW
     &      *Za11*Za13*Zh12*cb-EE*MW/SW*Za11*Za13*Zh22*sb
     &      -EE*MW/SW*Za21*Za23*Zh12*cb+EE*MW/SW*Za21*Za23
     &      *Zh22*sb+2D0*Sqrt2*Za11*Za23*Zh32*hK*hL*xvev
     &      -2D0*Sqrt2*Za11*Za33*Zh22*hK*hL*xvev+2D0*Sqrt2
     &      *Za13*Za21*Zh32*hK*hL*xvev-2D0*Sqrt2*Za13*Za31
     &      *Zh22*hK*hL*xvev-2D0*Sqrt2*Za21*Za33*Zh12*hK
     &      *hL*xvev-2D0*Sqrt2*Za23*Za31*Zh12*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za11*Za13*Zh12*cb-EE/CW**2*MW
     &      *SW*Za11*Za13*Zh22*sb-EE/CW**2*MW*SW*Za21*Za23
     &      *Zh12*cb+EE/CW**2*MW*SW*Za21*Za23*Zh22*sb+4D0
     &      *MW/EE*SW*Za11*Za13*Zh22*hL**2*sb+4D0*MW/EE*SW
     &      *Za21*Za23*Zh12*cb*hL**2+4D0*MW/EE*SW*Za31*Za33
     &      *Zh12*cb*hL**2+4D0*MW/EE*SW*Za31*Za33*Zh22*hL
     &      **2*sb-4D0*MW/EE*SW*Za11*Za33*Zh32*hK*hL*sb-4D0
     &      *MW/EE*SW*Za13*Za31*Zh32*hK*hL*sb-4D0*MW/EE*SW
     &      *Za21*Za33*Zh32*cb*hK*hL-4D0*MW/EE*SW*Za23*Za31
     &      *Zh32*cb*hK*hL+4D0*MW/EE*SW*Za31*Za33*Zh12*hK
     &      *hL*sb+4D0*MW/EE*SW*Za31*Za33*Zh22*cb*hK*hL
      AAABR(200) = 2D0*Sqrt2*Za12*Za13*Zh31*hL**2*xvev+Sqrt2
     &      *Za12*Za23*Zh31*hL*hLs+Sqrt2*Za12*Za33*Zh21*hL
     &      *hLs+Sqrt2*Za13*Za22*Zh31*hL*hLs+Sqrt2*Za13*Za32
     &      *Zh21*hL*hLs+2D0*Sqrt2*Za22*Za23*Zh31*hL**2*xvev
     &      +Sqrt2*Za22*Za33*Zh11*hL*hLs+Sqrt2*Za23*Za32
     &      *Zh11*hL*hLs-2D0*Sqrt2*Za32*Za33*Zh31*hK*hKs
     &      +4D0*Sqrt2*Za32*Za33*Zh31*hK**2*xvev+EE*MW/SW
     &      *Za12*Za13*Zh11*cb-EE*MW/SW*Za12*Za13*Zh21*sb
     &      -EE*MW/SW*Za22*Za23*Zh11*cb+EE*MW/SW*Za22*Za23
     &      *Zh21*sb+2D0*Sqrt2*Za12*Za23*Zh31*hK*hL*xvev
     &      -2D0*Sqrt2*Za12*Za33*Zh21*hK*hL*xvev+2D0*Sqrt2
     &      *Za13*Za22*Zh31*hK*hL*xvev-2D0*Sqrt2*Za13*Za32
     &      *Zh21*hK*hL*xvev-2D0*Sqrt2*Za22*Za33*Zh11*hK
     &      *hL*xvev-2D0*Sqrt2*Za23*Za32*Zh11*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za12*Za13*Zh11*cb-EE/CW**2*MW
     &      *SW*Za12*Za13*Zh21*sb-EE/CW**2*MW*SW*Za22*Za23
     &      *Zh11*cb+EE/CW**2*MW*SW*Za22*Za23*Zh21*sb+4D0
     &      *MW/EE*SW*Za12*Za13*Zh21*hL**2*sb+4D0*MW/EE*SW
     &      *Za22*Za23*Zh11*cb*hL**2+4D0*MW/EE*SW*Za32*Za33
     &      *Zh11*cb*hL**2+4D0*MW/EE*SW*Za32*Za33*Zh21*hL
     &      **2*sb-4D0*MW/EE*SW*Za12*Za33*Zh31*hK*hL*sb-4D0
     &      *MW/EE*SW*Za13*Za32*Zh31*hK*hL*sb-4D0*MW/EE*SW
     &      *Za22*Za33*Zh31*cb*hK*hL-4D0*MW/EE*SW*Za23*Za32
     &      *Zh31*cb*hK*hL+4D0*MW/EE*SW*Za32*Za33*Zh11*hK
     &      *hL*sb+4D0*MW/EE*SW*Za32*Za33*Zh21*cb*hK*hL
      AAABR(201) = 2D0*Sqrt2*Za11*Za13*Zh33*hL**2*xvev+Sqrt2
     &      *Za11*Za23*Zh33*hL*hLs+Sqrt2*Za11*Za33*Zh23*hL
     &      *hLs+Sqrt2*Za13*Za21*Zh33*hL*hLs+Sqrt2*Za13*Za31
     &      *Zh23*hL*hLs+2D0*Sqrt2*Za21*Za23*Zh33*hL**2*xvev
     &      +Sqrt2*Za21*Za33*Zh13*hL*hLs+Sqrt2*Za23*Za31
     &      *Zh13*hL*hLs-2D0*Sqrt2*Za31*Za33*Zh33*hK*hKs
     &      +4D0*Sqrt2*Za31*Za33*Zh33*hK**2*xvev+EE*MW/SW
     &      *Za11*Za13*Zh13*cb-EE*MW/SW*Za11*Za13*Zh23*sb
     &      -EE*MW/SW*Za21*Za23*Zh13*cb+EE*MW/SW*Za21*Za23
     &      *Zh23*sb+2D0*Sqrt2*Za11*Za23*Zh33*hK*hL*xvev
     &      -2D0*Sqrt2*Za11*Za33*Zh23*hK*hL*xvev+2D0*Sqrt2
     &      *Za13*Za21*Zh33*hK*hL*xvev-2D0*Sqrt2*Za13*Za31
     &      *Zh23*hK*hL*xvev-2D0*Sqrt2*Za21*Za33*Zh13*hK
     &      *hL*xvev-2D0*Sqrt2*Za23*Za31*Zh13*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za11*Za13*Zh13*cb-EE/CW**2*MW
     &      *SW*Za11*Za13*Zh23*sb-EE/CW**2*MW*SW*Za21*Za23
     &      *Zh13*cb+EE/CW**2*MW*SW*Za21*Za23*Zh23*sb+4D0
     &      *MW/EE*SW*Za11*Za13*Zh23*hL**2*sb+4D0*MW/EE*SW
     &      *Za21*Za23*Zh13*cb*hL**2+4D0*MW/EE*SW*Za31*Za33
     &      *Zh13*cb*hL**2+4D0*MW/EE*SW*Za31*Za33*Zh23*hL
     &      **2*sb-4D0*MW/EE*SW*Za11*Za33*Zh33*hK*hL*sb-4D0
     &      *MW/EE*SW*Za13*Za31*Zh33*hK*hL*sb-4D0*MW/EE*SW
     &      *Za21*Za33*Zh33*cb*hK*hL-4D0*MW/EE*SW*Za23*Za31
     &      *Zh33*cb*hK*hL+4D0*MW/EE*SW*Za31*Za33*Zh13*hK
     &      *hL*sb+4D0*MW/EE*SW*Za31*Za33*Zh23*cb*hK*hL
      AAABR(202) = 2D0*Sqrt2*Za11**2*Zh32*hL**2*xvev+2D0
     &      *Sqrt2*Za21**2*Zh32*hL**2*xvev-2D0*Sqrt2*Za31
     &      **2*Zh32*hK*hKs+4D0*Sqrt2*Za31**2*Zh32*hK**2
     &      *xvev+EE*MW/SW*Za11**2*Zh12*cb-EE*MW/SW*Za11
     &      **2*Zh22*sb-EE*MW/SW*Za21**2*Zh12*cb+EE*MW/SW
     &      *Za21**2*Zh22*sb+2D0*Sqrt2*Za11*Za21*Zh32*hL
     &      *hLs+2D0*Sqrt2*Za11*Za31*Zh22*hL*hLs+2D0*Sqrt2
     &      *Za21*Za31*Zh12*hL*hLs+EE/CW**2*MW*SW*Za11**2
     &      *Zh12*cb-EE/CW**2*MW*SW*Za11**2*Zh22*sb-EE/CW
     &      **2*MW*SW*Za21**2*Zh12*cb+EE/CW**2*MW*SW*Za21
     &      **2*Zh22*sb+4D0*MW/EE*SW*Za11**2*Zh22*hL**2*sb
     &      +4D0*MW/EE*SW*Za21**2*Zh12*cb*hL**2+4D0*MW/EE
     &      *SW*Za31**2*Zh12*cb*hL**2+4D0*MW/EE*SW*Za31**2
     &      *Zh22*hL**2*sb+4D0*Sqrt2*Za11*Za21*Zh32*hK*hL
     &      *xvev-4D0*Sqrt2*Za11*Za31*Zh22*hK*hL*xvev-4D0
     &      *Sqrt2*Za21*Za31*Zh12*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za31**2*Zh12*hK*hL*sb+4D0*MW/EE*SW*Za31**2*Zh22
     &      *cb*hK*hL-8D0*MW/EE*SW*Za11*Za31*Zh32*hK*hL*sb
     &      -8D0*MW/EE*SW*Za21*Za31*Zh32*cb*hK*hL
      AAABR(203) = 2D0*Sqrt2*Za11**2*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Za21**2*Zh33*hL**2*xvev-2D0*Sqrt2*Za31
     &      **2*Zh33*hK*hKs+4D0*Sqrt2*Za31**2*Zh33*hK**2
     &      *xvev+EE*MW/SW*Za11**2*Zh13*cb-EE*MW/SW*Za11
     &      **2*Zh23*sb-EE*MW/SW*Za21**2*Zh13*cb+EE*MW/SW
     &      *Za21**2*Zh23*sb+2D0*Sqrt2*Za11*Za21*Zh33*hL
     &      *hLs+2D0*Sqrt2*Za11*Za31*Zh23*hL*hLs+2D0*Sqrt2
     &      *Za21*Za31*Zh13*hL*hLs+EE/CW**2*MW*SW*Za11**2
     &      *Zh13*cb-EE/CW**2*MW*SW*Za11**2*Zh23*sb-EE/CW
     &      **2*MW*SW*Za21**2*Zh13*cb+EE/CW**2*MW*SW*Za21
     &      **2*Zh23*sb+4D0*MW/EE*SW*Za11**2*Zh23*hL**2*sb
     &      +4D0*MW/EE*SW*Za21**2*Zh13*cb*hL**2+4D0*MW/EE
     &      *SW*Za31**2*Zh13*cb*hL**2+4D0*MW/EE*SW*Za31**2
     &      *Zh23*hL**2*sb+4D0*Sqrt2*Za11*Za21*Zh33*hK*hL
     &      *xvev-4D0*Sqrt2*Za11*Za31*Zh23*hK*hL*xvev-4D0
     &      *Sqrt2*Za21*Za31*Zh13*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za31**2*Zh13*hK*hL*sb+4D0*MW/EE*SW*Za31**2*Zh23
     &      *cb*hK*hL-8D0*MW/EE*SW*Za11*Za31*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Za21*Za31*Zh33*cb*hK*hL
      AAABR(204) = 2D0*Sqrt2*Za12**2*Zh31*hL**2*xvev+2D0
     &      *Sqrt2*Za22**2*Zh31*hL**2*xvev-2D0*Sqrt2*Za32
     &      **2*Zh31*hK*hKs+4D0*Sqrt2*Za32**2*Zh31*hK**2
     &      *xvev+EE*MW/SW*Za12**2*Zh11*cb-EE*MW/SW*Za12
     &      **2*Zh21*sb-EE*MW/SW*Za22**2*Zh11*cb+EE*MW/SW
     &      *Za22**2*Zh21*sb+2D0*Sqrt2*Za12*Za22*Zh31*hL
     &      *hLs+2D0*Sqrt2*Za12*Za32*Zh21*hL*hLs+2D0*Sqrt2
     &      *Za22*Za32*Zh11*hL*hLs+EE/CW**2*MW*SW*Za12**2
     &      *Zh11*cb-EE/CW**2*MW*SW*Za12**2*Zh21*sb-EE/CW
     &      **2*MW*SW*Za22**2*Zh11*cb+EE/CW**2*MW*SW*Za22
     &      **2*Zh21*sb+4D0*MW/EE*SW*Za12**2*Zh21*hL**2*sb
     &      +4D0*MW/EE*SW*Za22**2*Zh11*cb*hL**2+4D0*MW/EE
     &      *SW*Za32**2*Zh11*cb*hL**2+4D0*MW/EE*SW*Za32**2
     &      *Zh21*hL**2*sb+4D0*Sqrt2*Za12*Za22*Zh31*hK*hL
     &      *xvev-4D0*Sqrt2*Za12*Za32*Zh21*hK*hL*xvev-4D0
     &      *Sqrt2*Za22*Za32*Zh11*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za32**2*Zh11*hK*hL*sb+4D0*MW/EE*SW*Za32**2*Zh21
     &      *cb*hK*hL-8D0*MW/EE*SW*Za12*Za32*Zh31*hK*hL*sb
     &      -8D0*MW/EE*SW*Za22*Za32*Zh31*cb*hK*hL
      AAABR(205) = 2D0*Sqrt2*Za11*Za12*Zh32*hL**2*xvev+Sqrt2
     &      *Za11*Za22*Zh32*hL*hLs+Sqrt2*Za11*Za32*Zh22*hL
     &      *hLs+Sqrt2*Za12*Za21*Zh32*hL*hLs+Sqrt2*Za12*Za31
     &      *Zh22*hL*hLs+2D0*Sqrt2*Za21*Za22*Zh32*hL**2*xvev
     &      +Sqrt2*Za21*Za32*Zh12*hL*hLs+Sqrt2*Za22*Za31
     &      *Zh12*hL*hLs-2D0*Sqrt2*Za31*Za32*Zh32*hK*hKs
     &      +4D0*Sqrt2*Za31*Za32*Zh32*hK**2*xvev+EE*MW/SW
     &      *Za11*Za12*Zh12*cb-EE*MW/SW*Za11*Za12*Zh22*sb
     &      -EE*MW/SW*Za21*Za22*Zh12*cb+EE*MW/SW*Za21*Za22
     &      *Zh22*sb+2D0*Sqrt2*Za11*Za22*Zh32*hK*hL*xvev
     &      -2D0*Sqrt2*Za11*Za32*Zh22*hK*hL*xvev+2D0*Sqrt2
     &      *Za12*Za21*Zh32*hK*hL*xvev-2D0*Sqrt2*Za12*Za31
     &      *Zh22*hK*hL*xvev-2D0*Sqrt2*Za21*Za32*Zh12*hK
     &      *hL*xvev-2D0*Sqrt2*Za22*Za31*Zh12*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za11*Za12*Zh12*cb-EE/CW**2*MW
     &      *SW*Za11*Za12*Zh22*sb-EE/CW**2*MW*SW*Za21*Za22
     &      *Zh12*cb+EE/CW**2*MW*SW*Za21*Za22*Zh22*sb+4D0
     &      *MW/EE*SW*Za11*Za12*Zh22*hL**2*sb+4D0*MW/EE*SW
     &      *Za21*Za22*Zh12*cb*hL**2+4D0*MW/EE*SW*Za31*Za32
     &      *Zh12*cb*hL**2+4D0*MW/EE*SW*Za31*Za32*Zh22*hL
     &      **2*sb-4D0*MW/EE*SW*Za11*Za32*Zh32*hK*hL*sb-4D0
     &      *MW/EE*SW*Za12*Za31*Zh32*hK*hL*sb-4D0*MW/EE*SW
     &      *Za21*Za32*Zh32*cb*hK*hL-4D0*MW/EE*SW*Za22*Za31
     &      *Zh32*cb*hK*hL+4D0*MW/EE*SW*Za31*Za32*Zh12*hK
     &      *hL*sb+4D0*MW/EE*SW*Za31*Za32*Zh22*cb*hK*hL
      AAABR(206) = 2D0*Sqrt2*Za12**2*Zh32*hL**2*xvev+2D0
     &      *Sqrt2*Za22**2*Zh32*hL**2*xvev-2D0*Sqrt2*Za32
     &      **2*Zh32*hK*hKs+4D0*Sqrt2*Za32**2*Zh32*hK**2
     &      *xvev+EE*MW/SW*Za12**2*Zh12*cb-EE*MW/SW*Za12
     &      **2*Zh22*sb-EE*MW/SW*Za22**2*Zh12*cb+EE*MW/SW
     &      *Za22**2*Zh22*sb+2D0*Sqrt2*Za12*Za22*Zh32*hL
     &      *hLs+2D0*Sqrt2*Za12*Za32*Zh22*hL*hLs+2D0*Sqrt2
     &      *Za22*Za32*Zh12*hL*hLs+EE/CW**2*MW*SW*Za12**2
     &      *Zh12*cb-EE/CW**2*MW*SW*Za12**2*Zh22*sb-EE/CW
     &      **2*MW*SW*Za22**2*Zh12*cb+EE/CW**2*MW*SW*Za22
     &      **2*Zh22*sb+4D0*MW/EE*SW*Za12**2*Zh22*hL**2*sb
     &      +4D0*MW/EE*SW*Za22**2*Zh12*cb*hL**2+4D0*MW/EE
     &      *SW*Za32**2*Zh12*cb*hL**2+4D0*MW/EE*SW*Za32**2
     &      *Zh22*hL**2*sb+4D0*Sqrt2*Za12*Za22*Zh32*hK*hL
     &      *xvev-4D0*Sqrt2*Za12*Za32*Zh22*hK*hL*xvev-4D0
     &      *Sqrt2*Za22*Za32*Zh12*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za32**2*Zh12*hK*hL*sb+4D0*MW/EE*SW*Za32**2*Zh22
     &      *cb*hK*hL-8D0*MW/EE*SW*Za12*Za32*Zh32*hK*hL*sb
     &      -8D0*MW/EE*SW*Za22*Za32*Zh32*cb*hK*hL
      AAABR(207) = 2D0*Sqrt2*Za11*Za12*Zh33*hL**2*xvev+Sqrt2
     &      *Za11*Za22*Zh33*hL*hLs+Sqrt2*Za11*Za32*Zh23*hL
     &      *hLs+Sqrt2*Za12*Za21*Zh33*hL*hLs+Sqrt2*Za12*Za31
     &      *Zh23*hL*hLs+2D0*Sqrt2*Za21*Za22*Zh33*hL**2*xvev
     &      +Sqrt2*Za21*Za32*Zh13*hL*hLs+Sqrt2*Za22*Za31
     &      *Zh13*hL*hLs-2D0*Sqrt2*Za31*Za32*Zh33*hK*hKs
     &      +4D0*Sqrt2*Za31*Za32*Zh33*hK**2*xvev+EE*MW/SW
     &      *Za11*Za12*Zh13*cb-EE*MW/SW*Za11*Za12*Zh23*sb
     &      -EE*MW/SW*Za21*Za22*Zh13*cb+EE*MW/SW*Za21*Za22
     &      *Zh23*sb+2D0*Sqrt2*Za11*Za22*Zh33*hK*hL*xvev
     &      -2D0*Sqrt2*Za11*Za32*Zh23*hK*hL*xvev+2D0*Sqrt2
     &      *Za12*Za21*Zh33*hK*hL*xvev-2D0*Sqrt2*Za12*Za31
     &      *Zh23*hK*hL*xvev-2D0*Sqrt2*Za21*Za32*Zh13*hK
     &      *hL*xvev-2D0*Sqrt2*Za22*Za31*Zh13*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za11*Za12*Zh13*cb-EE/CW**2*MW
     &      *SW*Za11*Za12*Zh23*sb-EE/CW**2*MW*SW*Za21*Za22
     &      *Zh13*cb+EE/CW**2*MW*SW*Za21*Za22*Zh23*sb+4D0
     &      *MW/EE*SW*Za11*Za12*Zh23*hL**2*sb+4D0*MW/EE*SW
     &      *Za21*Za22*Zh13*cb*hL**2+4D0*MW/EE*SW*Za31*Za32
     &      *Zh13*cb*hL**2+4D0*MW/EE*SW*Za31*Za32*Zh23*hL
     &      **2*sb-4D0*MW/EE*SW*Za11*Za32*Zh33*hK*hL*sb-4D0
     &      *MW/EE*SW*Za12*Za31*Zh33*hK*hL*sb-4D0*MW/EE*SW
     &      *Za21*Za32*Zh33*cb*hK*hL-4D0*MW/EE*SW*Za22*Za31
     &      *Zh33*cb*hK*hL+4D0*MW/EE*SW*Za31*Za32*Zh13*hK
     &      *hL*sb+4D0*MW/EE*SW*Za31*Za32*Zh23*cb*hK*hL
      AAABR(208) = 2D0*Sqrt2*Za12**2*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Za22**2*Zh33*hL**2*xvev-2D0*Sqrt2*Za32
     &      **2*Zh33*hK*hKs+4D0*Sqrt2*Za32**2*Zh33*hK**2
     &      *xvev+EE*MW/SW*Za12**2*Zh13*cb-EE*MW/SW*Za12
     &      **2*Zh23*sb-EE*MW/SW*Za22**2*Zh13*cb+EE*MW/SW
     &      *Za22**2*Zh23*sb+2D0*Sqrt2*Za12*Za22*Zh33*hL
     &      *hLs+2D0*Sqrt2*Za12*Za32*Zh23*hL*hLs+2D0*Sqrt2
     &      *Za22*Za32*Zh13*hL*hLs+EE/CW**2*MW*SW*Za12**2
     &      *Zh13*cb-EE/CW**2*MW*SW*Za12**2*Zh23*sb-EE/CW
     &      **2*MW*SW*Za22**2*Zh13*cb+EE/CW**2*MW*SW*Za22
     &      **2*Zh23*sb+4D0*MW/EE*SW*Za12**2*Zh23*hL**2*sb
     &      +4D0*MW/EE*SW*Za22**2*Zh13*cb*hL**2+4D0*MW/EE
     &      *SW*Za32**2*Zh13*cb*hL**2+4D0*MW/EE*SW*Za32**2
     &      *Zh23*hL**2*sb+4D0*Sqrt2*Za12*Za22*Zh33*hK*hL
     &      *xvev-4D0*Sqrt2*Za12*Za32*Zh23*hK*hL*xvev-4D0
     &      *Sqrt2*Za22*Za32*Zh13*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za32**2*Zh13*hK*hL*sb+4D0*MW/EE*SW*Za32**2*Zh23
     &      *cb*hK*hL-8D0*MW/EE*SW*Za12*Za32*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Za22*Za32*Zh33*cb*hK*hL
      AAABR(209) = 2D0*Sqrt2*Za11*Za12*Zh31*hL**2*xvev+Sqrt2
     &      *Za11*Za22*Zh31*hL*hLs+Sqrt2*Za11*Za32*Zh21*hL
     &      *hLs+Sqrt2*Za12*Za21*Zh31*hL*hLs+Sqrt2*Za12*Za31
     &      *Zh21*hL*hLs+2D0*Sqrt2*Za21*Za22*Zh31*hL**2*xvev
     &      +Sqrt2*Za21*Za32*Zh11*hL*hLs+Sqrt2*Za22*Za31
     &      *Zh11*hL*hLs-2D0*Sqrt2*Za31*Za32*Zh31*hK*hKs
     &      +4D0*Sqrt2*Za31*Za32*Zh31*hK**2*xvev+EE*MW/SW
     &      *Za11*Za12*Zh11*cb-EE*MW/SW*Za11*Za12*Zh21*sb
     &      -EE*MW/SW*Za21*Za22*Zh11*cb+EE*MW/SW*Za21*Za22
     &      *Zh21*sb+2D0*Sqrt2*Za11*Za22*Zh31*hK*hL*xvev
     &      -2D0*Sqrt2*Za11*Za32*Zh21*hK*hL*xvev+2D0*Sqrt2
     &      *Za12*Za21*Zh31*hK*hL*xvev-2D0*Sqrt2*Za12*Za31
     &      *Zh21*hK*hL*xvev-2D0*Sqrt2*Za21*Za32*Zh11*hK
     &      *hL*xvev-2D0*Sqrt2*Za22*Za31*Zh11*hK*hL*xvev
     &      +EE/CW**2*MW*SW*Za11*Za12*Zh11*cb-EE/CW**2*MW
     &      *SW*Za11*Za12*Zh21*sb-EE/CW**2*MW*SW*Za21*Za22
     &      *Zh11*cb+EE/CW**2*MW*SW*Za21*Za22*Zh21*sb+4D0
     &      *MW/EE*SW*Za11*Za12*Zh21*hL**2*sb+4D0*MW/EE*SW
     &      *Za21*Za22*Zh11*cb*hL**2+4D0*MW/EE*SW*Za31*Za32
     &      *Zh11*cb*hL**2+4D0*MW/EE*SW*Za31*Za32*Zh21*hL
     &      **2*sb-4D0*MW/EE*SW*Za11*Za32*Zh31*hK*hL*sb-4D0
     &      *MW/EE*SW*Za12*Za31*Zh31*hK*hL*sb-4D0*MW/EE*SW
     &      *Za21*Za32*Zh31*cb*hK*hL-4D0*MW/EE*SW*Za22*Za31
     &      *Zh31*cb*hK*hL+4D0*MW/EE*SW*Za31*Za32*Zh11*hK
     &      *hL*sb+4D0*MW/EE*SW*Za31*Za32*Zh21*cb*hK*hL
      AAABR(210) = 2D0*Sqrt2*Za11**2*Zh31*hL**2*xvev+2D0
     &      *Sqrt2*Za21**2*Zh31*hL**2*xvev-2D0*Sqrt2*Za31
     &      **2*Zh31*hK*hKs+4D0*Sqrt2*Za31**2*Zh31*hK**2
     &      *xvev+EE*MW/SW*Za11**2*Zh11*cb-EE*MW/SW*Za11
     &      **2*Zh21*sb-EE*MW/SW*Za21**2*Zh11*cb+EE*MW/SW
     &      *Za21**2*Zh21*sb+2D0*Sqrt2*Za11*Za21*Zh31*hL
     &      *hLs+2D0*Sqrt2*Za11*Za31*Zh21*hL*hLs+2D0*Sqrt2
     &      *Za21*Za31*Zh11*hL*hLs+EE/CW**2*MW*SW*Za11**2
     &      *Zh11*cb-EE/CW**2*MW*SW*Za11**2*Zh21*sb-EE/CW
     &      **2*MW*SW*Za21**2*Zh11*cb+EE/CW**2*MW*SW*Za21
     &      **2*Zh21*sb+4D0*MW/EE*SW*Za11**2*Zh21*hL**2*sb
     &      +4D0*MW/EE*SW*Za21**2*Zh11*cb*hL**2+4D0*MW/EE
     &      *SW*Za31**2*Zh11*cb*hL**2+4D0*MW/EE*SW*Za31**2
     &      *Zh21*hL**2*sb+4D0*Sqrt2*Za11*Za21*Zh31*hK*hL
     &      *xvev-4D0*Sqrt2*Za11*Za31*Zh21*hK*hL*xvev-4D0
     &      *Sqrt2*Za21*Za31*Zh11*hK*hL*xvev+4D0*MW/EE*SW
     &      *Za31**2*Zh11*hK*hL*sb+4D0*MW/EE*SW*Za31**2*Zh21
     &      *cb*hK*hL-8D0*MW/EE*SW*Za11*Za31*Zh31*hK*hL*sb
     &      -8D0*MW/EE*SW*Za21*Za31*Zh31*cb*hK*hL
      AAABR(211) = 2D0*Sqrt2*Zh11**2*Zh32*hL**2*xvev+2D0
     &      *Sqrt2*Zh21**2*Zh32*hL**2*xvev+2D0*Sqrt2*Zh31
     &      **2*Zh32*hK*hKs+12D0*Sqrt2*Zh31**2*Zh32*hK**2
     &      *xvev+3D0*EE*MW/SW*Zh11**2*Zh12*cb-EE*MW/SW*Zh11
     &      **2*Zh22*sb-EE*MW/SW*Zh12*Zh21**2*cb+3D0*EE*MW
     &      /SW*Zh21**2*Zh22*sb+4D0*Sqrt2*Zh11*Zh12*Zh31
     &      *hL**2*xvev-2D0*Sqrt2*Zh11*Zh21*Zh32*hL*hLs-2D0
     &      *Sqrt2*Zh11*Zh22*Zh31*hL*hLs-2D0*Sqrt2*Zh12*Zh21
     &      *Zh31*hL*hLs+4D0*Sqrt2*Zh21*Zh22*Zh31*hL**2*xvev
     &      +3D0*EE/CW**2*MW*SW*Zh11**2*Zh12*cb-EE/CW**2
     &      *MW*SW*Zh11**2*Zh22*sb-EE/CW**2*MW*SW*Zh12*Zh21
     &      **2*cb+3D0*EE/CW**2*MW*SW*Zh21**2*Zh22*sb+4D0
     &      *MW/EE*SW*Zh11**2*Zh22*hL**2*sb+4D0*MW/EE*SW
     &      *Zh12*Zh21**2*cb*hL**2+4D0*MW/EE*SW*Zh12*Zh31
     &      **2*cb*hL**2+4D0*MW/EE*SW*Zh22*Zh31**2*hL**2
     &      *sb-2D0*EE*MW/SW*Zh11*Zh12*Zh21*sb-2D0*EE*MW
     &      /SW*Zh11*Zh21*Zh22*cb-4D0*Sqrt2*Zh11*Zh21*Zh32
     &      *hK*hL*xvev-4D0*Sqrt2*Zh11*Zh22*Zh31*hK*hL*xvev
     &      -4D0*Sqrt2*Zh12*Zh21*Zh31*hK*hL*xvev-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zh12*Zh21*sb-2D0*EE/CW**2*MW*SW
     &      *Zh11*Zh21*Zh22*cb+8D0*MW/EE*SW*Zh11*Zh12*Zh21
     &      *hL**2*sb+8D0*MW/EE*SW*Zh11*Zh21*Zh22*cb*hL**2
     &      +8D0*MW/EE*SW*Zh11*Zh31*Zh32*cb*hL**2-4D0*MW
     &      /EE*SW*Zh12*Zh31**2*hK*hL*sb+8D0*MW/EE*SW*Zh21
     &      *Zh31*Zh32*hL**2*sb-4D0*MW/EE*SW*Zh22*Zh31**2
     &      *cb*hK*hL-8D0*MW/EE*SW*Zh11*Zh31*Zh32*hK*hL*sb
     &      -8D0*MW/EE*SW*Zh21*Zh31*Zh32*cb*hK*hL
      AAABR(212) = 2D0*Sqrt2*Zh12**2*Zh31*hL**2*xvev+2D0
     &      *Sqrt2*Zh22**2*Zh31*hL**2*xvev+2D0*Sqrt2*Zh31
     &      *Zh32**2*hK*hKs+12D0*Sqrt2*Zh31*Zh32**2*hK**2
     &      *xvev+3D0*EE*MW/SW*Zh11*Zh12**2*cb-EE*MW/SW*Zh11
     &      *Zh22**2*cb-EE*MW/SW*Zh12**2*Zh21*sb+3D0*EE*MW
     &      /SW*Zh21*Zh22**2*sb+4D0*Sqrt2*Zh11*Zh12*Zh32
     &      *hL**2*xvev-2D0*Sqrt2*Zh11*Zh22*Zh32*hL*hLs-2D0
     &      *Sqrt2*Zh12*Zh21*Zh32*hL*hLs-2D0*Sqrt2*Zh12*Zh22
     &      *Zh31*hL*hLs+4D0*Sqrt2*Zh21*Zh22*Zh32*hL**2*xvev
     &      +3D0*EE/CW**2*MW*SW*Zh11*Zh12**2*cb-EE/CW**2
     &      *MW*SW*Zh11*Zh22**2*cb-EE/CW**2*MW*SW*Zh12**2
     &      *Zh21*sb+3D0*EE/CW**2*MW*SW*Zh21*Zh22**2*sb+4D0
     &      *MW/EE*SW*Zh11*Zh22**2*cb*hL**2+4D0*MW/EE*SW
     &      *Zh11*Zh32**2*cb*hL**2+4D0*MW/EE*SW*Zh12**2*Zh21
     &      *hL**2*sb+4D0*MW/EE*SW*Zh21*Zh32**2*hL**2*sb
     &      -2D0*EE*MW/SW*Zh11*Zh12*Zh22*sb-2D0*EE*MW/SW
     &      *Zh12*Zh21*Zh22*cb-4D0*Sqrt2*Zh11*Zh22*Zh32*hK
     &      *hL*xvev-4D0*Sqrt2*Zh12*Zh21*Zh32*hK*hL*xvev
     &      -4D0*Sqrt2*Zh12*Zh22*Zh31*hK*hL*xvev-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zh12*Zh22*sb-2D0*EE/CW**2*MW*SW
     &      *Zh12*Zh21*Zh22*cb+8D0*MW/EE*SW*Zh11*Zh12*Zh22
     &      *hL**2*sb-4D0*MW/EE*SW*Zh11*Zh32**2*hK*hL*sb
     &      +8D0*MW/EE*SW*Zh12*Zh21*Zh22*cb*hL**2+8D0*MW
     &      /EE*SW*Zh12*Zh31*Zh32*cb*hL**2-4D0*MW/EE*SW*Zh21
     &      *Zh32**2*cb*hK*hL+8D0*MW/EE*SW*Zh22*Zh31*Zh32
     &      *hL**2*sb-8D0*MW/EE*SW*Zh12*Zh31*Zh32*hK*hL*sb
     &      -8D0*MW/EE*SW*Zh22*Zh31*Zh32*cb*hK*hL
      AAABR(213) = 2D0*Sqrt2*Zh11*Zh12*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Zh11*Zh13*Zh32*hL**2*xvev-Sqrt2*Zh11*Zh22
     &      *Zh33*hL*hLs-Sqrt2*Zh11*Zh23*Zh32*hL*hLs+2D0
     &      *Sqrt2*Zh12*Zh13*Zh31*hL**2*xvev-Sqrt2*Zh12*Zh21
     &      *Zh33*hL*hLs-Sqrt2*Zh12*Zh23*Zh31*hL*hLs-Sqrt2
     &      *Zh13*Zh21*Zh32*hL*hLs-Sqrt2*Zh13*Zh22*Zh31*hL
     &      *hLs+2D0*Sqrt2*Zh21*Zh22*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Zh21*Zh23*Zh32*hL**2*xvev+2D0*Sqrt2*Zh22
     &      *Zh23*Zh31*hL**2*xvev+2D0*Sqrt2*Zh31*Zh32*Zh33
     &      *hK*hKs+12D0*Sqrt2*Zh31*Zh32*Zh33*hK**2*xvev
     &      +3D0*EE*MW/SW*Zh11*Zh12*Zh13*cb-EE*MW/SW*Zh11
     &      *Zh12*Zh23*sb-EE*MW/SW*Zh11*Zh13*Zh22*sb-EE*MW
     &      /SW*Zh11*Zh22*Zh23*cb-EE*MW/SW*Zh12*Zh13*Zh21
     &      *sb-EE*MW/SW*Zh12*Zh21*Zh23*cb-EE*MW/SW*Zh13
     &      *Zh21*Zh22*cb+3D0*EE*MW/SW*Zh21*Zh22*Zh23*sb
     &      -2D0*Sqrt2*Zh11*Zh22*Zh33*hK*hL*xvev-2D0*Sqrt2
     &      *Zh11*Zh23*Zh32*hK*hL*xvev-2D0*Sqrt2*Zh12*Zh21
     &      *Zh33*hK*hL*xvev-2D0*Sqrt2*Zh12*Zh23*Zh31*hK
     &      *hL*xvev-2D0*Sqrt2*Zh13*Zh21*Zh32*hK*hL*xvev
     &      -2D0*Sqrt2*Zh13*Zh22*Zh31*hK*hL*xvev+3D0*EE/CW
     &      **2*MW*SW*Zh11*Zh12*Zh13*cb-EE/CW**2*MW*SW*Zh11
     &      *Zh12*Zh23*sb-EE/CW**2*MW*SW*Zh11*Zh13*Zh22*sb
     &      -EE/CW**2*MW*SW*Zh11*Zh22*Zh23*cb-EE/CW**2*MW
     &      *SW*Zh12*Zh13*Zh21*sb-EE/CW**2*MW*SW*Zh12*Zh21
     &      *Zh23*cb-EE/CW**2*MW*SW*Zh13*Zh21*Zh22*cb+3D0
     &      *EE/CW**2*MW*SW*Zh21*Zh22*Zh23*sb+4D0*MW/EE*SW
     &      *Zh11*Zh12*Zh23*hL**2*sb+4D0*MW/EE*SW*Zh11*Zh13
     &      *Zh22*hL**2*sb+4D0*MW/EE*SW*Zh11*Zh22*Zh23*cb
     &      *hL**2+4D0*MW/EE*SW*Zh11*Zh32*Zh33*cb*hL**2+4D0
     &      *MW/EE*SW*Zh12*Zh13*Zh21*hL**2*sb+4D0*MW/EE*SW
     &      *Zh12*Zh21*Zh23*cb*hL**2+4D0*MW/EE*SW*Zh12*Zh31
     &      *Zh33*cb*hL**2+4D0*MW/EE*SW*Zh13*Zh21*Zh22*cb
     &      *hL**2+4D0*MW/EE*SW*Zh13*Zh31*Zh32*cb*hL**2+4D0
     &      *MW/EE*SW*Zh21*Zh32*Zh33*hL**2*sb+4D0*MW/EE*SW
     &      *Zh22*Zh31*Zh33*hL**2*sb+4D0*MW/EE*SW*Zh23*Zh31
     &      *Zh32*hL**2*sb-4D0*MW/EE*SW*Zh11*Zh32*Zh33*hK
     &      *hL*sb-4D0*MW/EE*SW*Zh12*Zh31*Zh33*hK*hL*sb-4D0
     &      *MW/EE*SW*Zh13*Zh31*Zh32*hK*hL*sb-4D0*MW/EE*SW
     &      *Zh21*Zh32*Zh33*cb*hK*hL-4D0*MW/EE*SW*Zh22*Zh31
     &      *Zh33*cb*hK*hL-4D0*MW/EE*SW*Zh23*Zh31*Zh32*cb
     &      *hK*hL
      AAABR(214) = 2D0*Sqrt2*Zh13**2*Zh31*hL**2*xvev+2D0
     &      *Sqrt2*Zh23**2*Zh31*hL**2*xvev+2D0*Sqrt2*Zh31
     &      *Zh33**2*hK*hKs+12D0*Sqrt2*Zh31*Zh33**2*hK**2
     &      *xvev+3D0*EE*MW/SW*Zh11*Zh13**2*cb-EE*MW/SW*Zh11
     &      *Zh23**2*cb-EE*MW/SW*Zh13**2*Zh21*sb+3D0*EE*MW
     &      /SW*Zh21*Zh23**2*sb+4D0*Sqrt2*Zh11*Zh13*Zh33
     &      *hL**2*xvev-2D0*Sqrt2*Zh11*Zh23*Zh33*hL*hLs-2D0
     &      *Sqrt2*Zh13*Zh21*Zh33*hL*hLs-2D0*Sqrt2*Zh13*Zh23
     &      *Zh31*hL*hLs+4D0*Sqrt2*Zh21*Zh23*Zh33*hL**2*xvev
     &      +3D0*EE/CW**2*MW*SW*Zh11*Zh13**2*cb-EE/CW**2
     &      *MW*SW*Zh11*Zh23**2*cb-EE/CW**2*MW*SW*Zh13**2
     &      *Zh21*sb+3D0*EE/CW**2*MW*SW*Zh21*Zh23**2*sb+4D0
     &      *MW/EE*SW*Zh11*Zh23**2*cb*hL**2+4D0*MW/EE*SW
     &      *Zh11*Zh33**2*cb*hL**2+4D0*MW/EE*SW*Zh13**2*Zh21
     &      *hL**2*sb+4D0*MW/EE*SW*Zh21*Zh33**2*hL**2*sb
     &      -2D0*EE*MW/SW*Zh11*Zh13*Zh23*sb-2D0*EE*MW/SW
     &      *Zh13*Zh21*Zh23*cb-4D0*Sqrt2*Zh11*Zh23*Zh33*hK
     &      *hL*xvev-4D0*Sqrt2*Zh13*Zh21*Zh33*hK*hL*xvev
     &      -4D0*Sqrt2*Zh13*Zh23*Zh31*hK*hL*xvev-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zh13*Zh23*sb-2D0*EE/CW**2*MW*SW
     &      *Zh13*Zh21*Zh23*cb+8D0*MW/EE*SW*Zh11*Zh13*Zh23
     &      *hL**2*sb-4D0*MW/EE*SW*Zh11*Zh33**2*hK*hL*sb
     &      +8D0*MW/EE*SW*Zh13*Zh21*Zh23*cb*hL**2+8D0*MW
     &      /EE*SW*Zh13*Zh31*Zh33*cb*hL**2-4D0*MW/EE*SW*Zh21
     &      *Zh33**2*cb*hK*hL+8D0*MW/EE*SW*Zh23*Zh31*Zh33
     &      *hL**2*sb-8D0*MW/EE*SW*Zh13*Zh31*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Zh23*Zh31*Zh33*cb*hK*hL
      AAABR(215) = 2D0*Sqrt2*Zh13**2*Zh32*hL**2*xvev+2D0
     &      *Sqrt2*Zh23**2*Zh32*hL**2*xvev+2D0*Sqrt2*Zh32
     &      *Zh33**2*hK*hKs+12D0*Sqrt2*Zh32*Zh33**2*hK**2
     &      *xvev+3D0*EE*MW/SW*Zh12*Zh13**2*cb-EE*MW/SW*Zh12
     &      *Zh23**2*cb-EE*MW/SW*Zh13**2*Zh22*sb+3D0*EE*MW
     &      /SW*Zh22*Zh23**2*sb+4D0*Sqrt2*Zh12*Zh13*Zh33
     &      *hL**2*xvev-2D0*Sqrt2*Zh12*Zh23*Zh33*hL*hLs-2D0
     &      *Sqrt2*Zh13*Zh22*Zh33*hL*hLs-2D0*Sqrt2*Zh13*Zh23
     &      *Zh32*hL*hLs+4D0*Sqrt2*Zh22*Zh23*Zh33*hL**2*xvev
     &      +3D0*EE/CW**2*MW*SW*Zh12*Zh13**2*cb-EE/CW**2
     &      *MW*SW*Zh12*Zh23**2*cb-EE/CW**2*MW*SW*Zh13**2
     &      *Zh22*sb+3D0*EE/CW**2*MW*SW*Zh22*Zh23**2*sb+4D0
     &      *MW/EE*SW*Zh12*Zh23**2*cb*hL**2+4D0*MW/EE*SW
     &      *Zh12*Zh33**2*cb*hL**2+4D0*MW/EE*SW*Zh13**2*Zh22
     &      *hL**2*sb+4D0*MW/EE*SW*Zh22*Zh33**2*hL**2*sb
     &      -2D0*EE*MW/SW*Zh12*Zh13*Zh23*sb-2D0*EE*MW/SW
     &      *Zh13*Zh22*Zh23*cb-4D0*Sqrt2*Zh12*Zh23*Zh33*hK
     &      *hL*xvev-4D0*Sqrt2*Zh13*Zh22*Zh33*hK*hL*xvev
     &      -4D0*Sqrt2*Zh13*Zh23*Zh32*hK*hL*xvev-2D0*EE/CW
     &      **2*MW*SW*Zh12*Zh13*Zh23*sb-2D0*EE/CW**2*MW*SW
     &      *Zh13*Zh22*Zh23*cb+8D0*MW/EE*SW*Zh12*Zh13*Zh23
     &      *hL**2*sb-4D0*MW/EE*SW*Zh12*Zh33**2*hK*hL*sb
     &      +8D0*MW/EE*SW*Zh13*Zh22*Zh23*cb*hL**2+8D0*MW
     &      /EE*SW*Zh13*Zh32*Zh33*cb*hL**2-4D0*MW/EE*SW*Zh22
     &      *Zh33**2*cb*hK*hL+8D0*MW/EE*SW*Zh23*Zh32*Zh33
     &      *hL**2*sb-8D0*MW/EE*SW*Zh13*Zh32*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Zh23*Zh32*Zh33*cb*hK*hL
      AAABR(216) = 2D0*Sqrt2*Zh33**3*hK*hKs+12D0*Sqrt2*Zh33
     &      **3*hK**2*xvev+3D0*EE*MW/SW*Zh13**3*cb+3D0*EE
     &      *MW/SW*Zh23**3*sb+6D0*Sqrt2*Zh13**2*Zh33*hL**2
     &      *xvev+6D0*Sqrt2*Zh23**2*Zh33*hL**2*xvev+3D0*EE
     &      /CW**2*MW*SW*Zh13**3*cb+3D0*EE/CW**2*MW*SW*Zh23
     &      **3*sb-3D0*EE*MW/SW*Zh13*Zh23**2*cb-3D0*EE*MW
     &      /SW*Zh13**2*Zh23*sb-6D0*Sqrt2*Zh13*Zh23*Zh33
     &      *hL*hLs-3D0*EE/CW**2*MW*SW*Zh13*Zh23**2*cb-3D0
     &      *EE/CW**2*MW*SW*Zh13**2*Zh23*sb+12D0*MW/EE*SW
     &      *Zh13*Zh23**2*cb*hL**2+12D0*MW/EE*SW*Zh13*Zh33
     &      **2*cb*hL**2+12D0*MW/EE*SW*Zh13**2*Zh23*hL**2
     &      *sb+12D0*MW/EE*SW*Zh23*Zh33**2*hL**2*sb-12D0
     &      *Sqrt2*Zh13*Zh23*Zh33*hK*hL*xvev-12D0*MW/EE*SW
     &      *Zh13*Zh33**2*hK*hL*sb-12D0*MW/EE*SW*Zh23*Zh33
     &      **2*cb*hK*hL
      AAABR(217) = 2D0*Sqrt2*Zh31**3*hK*hKs+12D0*Sqrt2*Zh31
     &      **3*hK**2*xvev+3D0*EE*MW/SW*Zh11**3*cb+3D0*EE
     &      *MW/SW*Zh21**3*sb+6D0*Sqrt2*Zh11**2*Zh31*hL**2
     &      *xvev+6D0*Sqrt2*Zh21**2*Zh31*hL**2*xvev+3D0*EE
     &      /CW**2*MW*SW*Zh11**3*cb+3D0*EE/CW**2*MW*SW*Zh21
     &      **3*sb-3D0*EE*MW/SW*Zh11*Zh21**2*cb-3D0*EE*MW
     &      /SW*Zh11**2*Zh21*sb-6D0*Sqrt2*Zh11*Zh21*Zh31
     &      *hL*hLs-3D0*EE/CW**2*MW*SW*Zh11*Zh21**2*cb-3D0
     &      *EE/CW**2*MW*SW*Zh11**2*Zh21*sb+12D0*MW/EE*SW
     &      *Zh11*Zh21**2*cb*hL**2+12D0*MW/EE*SW*Zh11*Zh31
     &      **2*cb*hL**2+12D0*MW/EE*SW*Zh11**2*Zh21*hL**2
     &      *sb+12D0*MW/EE*SW*Zh21*Zh31**2*hL**2*sb-12D0
     &      *Sqrt2*Zh11*Zh21*Zh31*hK*hL*xvev-12D0*MW/EE*SW
     &      *Zh11*Zh31**2*hK*hL*sb-12D0*MW/EE*SW*Zh21*Zh31
     &      **2*cb*hK*hL
      AAABR(218) = 2D0*Sqrt2*Zh12**2*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Zh22**2*Zh33*hL**2*xvev+2D0*Sqrt2*Zh32
     &      **2*Zh33*hK*hKs+12D0*Sqrt2*Zh32**2*Zh33*hK**2
     &      *xvev+3D0*EE*MW/SW*Zh12**2*Zh13*cb-EE*MW/SW*Zh12
     &      **2*Zh23*sb-EE*MW/SW*Zh13*Zh22**2*cb+3D0*EE*MW
     &      /SW*Zh22**2*Zh23*sb+4D0*Sqrt2*Zh12*Zh13*Zh32
     &      *hL**2*xvev-2D0*Sqrt2*Zh12*Zh22*Zh33*hL*hLs-2D0
     &      *Sqrt2*Zh12*Zh23*Zh32*hL*hLs-2D0*Sqrt2*Zh13*Zh22
     &      *Zh32*hL*hLs+4D0*Sqrt2*Zh22*Zh23*Zh32*hL**2*xvev
     &      +3D0*EE/CW**2*MW*SW*Zh12**2*Zh13*cb-EE/CW**2
     &      *MW*SW*Zh12**2*Zh23*sb-EE/CW**2*MW*SW*Zh13*Zh22
     &      **2*cb+3D0*EE/CW**2*MW*SW*Zh22**2*Zh23*sb+4D0
     &      *MW/EE*SW*Zh12**2*Zh23*hL**2*sb+4D0*MW/EE*SW
     &      *Zh13*Zh22**2*cb*hL**2+4D0*MW/EE*SW*Zh13*Zh32
     &      **2*cb*hL**2+4D0*MW/EE*SW*Zh23*Zh32**2*hL**2
     &      *sb-2D0*EE*MW/SW*Zh12*Zh13*Zh22*sb-2D0*EE*MW
     &      /SW*Zh12*Zh22*Zh23*cb-4D0*Sqrt2*Zh12*Zh22*Zh33
     &      *hK*hL*xvev-4D0*Sqrt2*Zh12*Zh23*Zh32*hK*hL*xvev
     &      -4D0*Sqrt2*Zh13*Zh22*Zh32*hK*hL*xvev-2D0*EE/CW
     &      **2*MW*SW*Zh12*Zh13*Zh22*sb-2D0*EE/CW**2*MW*SW
     &      *Zh12*Zh22*Zh23*cb+8D0*MW/EE*SW*Zh12*Zh13*Zh22
     &      *hL**2*sb+8D0*MW/EE*SW*Zh12*Zh22*Zh23*cb*hL**2
     &      +8D0*MW/EE*SW*Zh12*Zh32*Zh33*cb*hL**2-4D0*MW
     &      /EE*SW*Zh13*Zh32**2*hK*hL*sb+8D0*MW/EE*SW*Zh22
     &      *Zh32*Zh33*hL**2*sb-4D0*MW/EE*SW*Zh23*Zh32**2
     &      *cb*hK*hL-8D0*MW/EE*SW*Zh12*Zh32*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Zh22*Zh32*Zh33*cb*hK*hL
      AAABR(219) = 2D0*Sqrt2*Zh11**2*Zh33*hL**2*xvev+2D0
     &      *Sqrt2*Zh21**2*Zh33*hL**2*xvev+2D0*Sqrt2*Zh31
     &      **2*Zh33*hK*hKs+12D0*Sqrt2*Zh31**2*Zh33*hK**2
     &      *xvev+3D0*EE*MW/SW*Zh11**2*Zh13*cb-EE*MW/SW*Zh11
     &      **2*Zh23*sb-EE*MW/SW*Zh13*Zh21**2*cb+3D0*EE*MW
     &      /SW*Zh21**2*Zh23*sb+4D0*Sqrt2*Zh11*Zh13*Zh31
     &      *hL**2*xvev-2D0*Sqrt2*Zh11*Zh21*Zh33*hL*hLs-2D0
     &      *Sqrt2*Zh11*Zh23*Zh31*hL*hLs-2D0*Sqrt2*Zh13*Zh21
     &      *Zh31*hL*hLs+4D0*Sqrt2*Zh21*Zh23*Zh31*hL**2*xvev
     &      +3D0*EE/CW**2*MW*SW*Zh11**2*Zh13*cb-EE/CW**2
     &      *MW*SW*Zh11**2*Zh23*sb-EE/CW**2*MW*SW*Zh13*Zh21
     &      **2*cb+3D0*EE/CW**2*MW*SW*Zh21**2*Zh23*sb+4D0
     &      *MW/EE*SW*Zh11**2*Zh23*hL**2*sb+4D0*MW/EE*SW
     &      *Zh13*Zh21**2*cb*hL**2+4D0*MW/EE*SW*Zh13*Zh31
     &      **2*cb*hL**2+4D0*MW/EE*SW*Zh23*Zh31**2*hL**2
     &      *sb-2D0*EE*MW/SW*Zh11*Zh13*Zh21*sb-2D0*EE*MW
     &      /SW*Zh11*Zh21*Zh23*cb-4D0*Sqrt2*Zh11*Zh21*Zh33
     &      *hK*hL*xvev-4D0*Sqrt2*Zh11*Zh23*Zh31*hK*hL*xvev
     &      -4D0*Sqrt2*Zh13*Zh21*Zh31*hK*hL*xvev-2D0*EE/CW
     &      **2*MW*SW*Zh11*Zh13*Zh21*sb-2D0*EE/CW**2*MW*SW
     &      *Zh11*Zh21*Zh23*cb+8D0*MW/EE*SW*Zh11*Zh13*Zh21
     &      *hL**2*sb+8D0*MW/EE*SW*Zh11*Zh21*Zh23*cb*hL**2
     &      +8D0*MW/EE*SW*Zh11*Zh31*Zh33*cb*hL**2-4D0*MW
     &      /EE*SW*Zh13*Zh31**2*hK*hL*sb+8D0*MW/EE*SW*Zh21
     &      *Zh31*Zh33*hL**2*sb-4D0*MW/EE*SW*Zh23*Zh31**2
     &      *cb*hK*hL-8D0*MW/EE*SW*Zh11*Zh31*Zh33*hK*hL*sb
     &      -8D0*MW/EE*SW*Zh21*Zh31*Zh33*cb*hK*hL
      AAABR(220) = 2D0*Sqrt2*Zh32**3*hK*hKs+12D0*Sqrt2*Zh32
     &      **3*hK**2*xvev+3D0*EE*MW/SW*Zh12**3*cb+3D0*EE
     &      *MW/SW*Zh22**3*sb+6D0*Sqrt2*Zh12**2*Zh32*hL**2
     &      *xvev+6D0*Sqrt2*Zh22**2*Zh32*hL**2*xvev+3D0*EE
     &      /CW**2*MW*SW*Zh12**3*cb+3D0*EE/CW**2*MW*SW*Zh22
     &      **3*sb-3D0*EE*MW/SW*Zh12*Zh22**2*cb-3D0*EE*MW
     &      /SW*Zh12**2*Zh22*sb-6D0*Sqrt2*Zh12*Zh22*Zh32
     &      *hL*hLs-3D0*EE/CW**2*MW*SW*Zh12*Zh22**2*cb-3D0
     &      *EE/CW**2*MW*SW*Zh12**2*Zh22*sb+12D0*MW/EE*SW
     &      *Zh12*Zh22**2*cb*hL**2+12D0*MW/EE*SW*Zh12*Zh32
     &      **2*cb*hL**2+12D0*MW/EE*SW*Zh12**2*Zh22*hL**2
     &      *sb+12D0*MW/EE*SW*Zh22*Zh32**2*hL**2*sb-12D0
     &      *Sqrt2*Zh12*Zh22*Zh32*hK*hL*xvev-12D0*MW/EE*SW
     &      *Zh12*Zh32**2*hK*hL*sb-12D0*MW/EE*SW*Zh22*Zh32
     &      **2*cb*hK*hL
      AAABR(221) = EE*(SW/CW-CW/SW)
      AAABR(222) = EE/SW*(Za13*sb+Za23*cb)
      AAABR(223) = EE/SW*(Za12*sb+Za22*cb)
      AAABR(224) = EE/SW*(Za11*sb+Za21*cb)
      AAABR(225) = EE/SW*(Zh13*sb-Zh23*cb)
      AAABR(226) = EE/SW*(Zh12*sb-Zh22*cb)
      AAABR(227) = EE/SW*(Zh11*sb-Zh21*cb)
      AAABR(228) = EE*(SW/CW*Zl33**2+2D0*SW/CW*Zl63**2-CW
     &      /SW*Zl33**2)
      AAABR(229) = EE*(SW/CW*Zl11**2+2D0*SW/CW*Zl41**2-CW
     &      /SW*Zl11**2)
      AAABR(230) = EE*(SW/CW*Zl22**2+2D0*SW/CW*Zl52**2-CW
     &      /SW*Zl22**2)
      AAABR(231) = EE*(SW/CW*Zl22*Zl25+2D0*SW/CW*Zl52*Zl55
     &      -CW/SW*Zl22*Zl25)
      AAABR(232) = EE*(SW/CW*Zl11*Zl14+2D0*SW/CW*Zl41*Zl44
     &      -CW/SW*Zl11*Zl14)
      AAABR(233) = EE*(SW/CW*Zl33*Zl36+2D0*SW/CW*Zl63*Zl66
     &      -CW/SW*Zl33*Zl36)
      AAABR(234) = EE/SW*Sqrt2*Zl33
      AAABR(235) = EE/SW*Sqrt2*Zl11
      AAABR(236) = EE/SW*Sqrt2*Zl22
      AAABR(237) = EE*(SW/CW*Zl14**2+2D0*SW/CW*Zl44**2-CW
     &      /SW*Zl14**2)
      AAABR(238) = EE*(SW/CW*Zl36**2+2D0*SW/CW*Zl66**2-CW
     &      /SW*Zl36**2)
      AAABR(239) = EE*(SW/CW*Zl25**2+2D0*SW/CW*Zl55**2-CW
     &      /SW*Zl25**2)
      AAABR(240) = EE/SW*Sqrt2*Zl36
      AAABR(241) = EE/SW*Sqrt2*Zl25
      AAABR(242) = EE/SW*Sqrt2*Zl14
      AAABR(243) = EE*(SW/CW+CW/SW)
      AAABR(244) = EE*(SW/CW*Zd33**2-2D0*SW/CW*Zd63**2+3D0
     &      *CW/SW*Zd33**2)
      AAABR(245) = EE*(SW/CW*Zd22**2-2D0*SW/CW*Zd52**2+3D0
     &      *CW/SW*Zd22**2)
      AAABR(246) = EE*(SW/CW*Zd11**2-2D0*SW/CW*Zd41**2+3D0
     &      *CW/SW*Zd11**2)
      AAABR(247) = EE*(SW/CW*Zd22*Zd25-2D0*SW/CW*Zd52*Zd55
     &      +3D0*CW/SW*Zd22*Zd25)
      AAABR(248) = EE*(SW/CW*Zd33*Zd36-2D0*SW/CW*Zd63*Zd66
     &      +3D0*CW/SW*Zd33*Zd36)
      AAABR(249) = EE*(SW/CW*Zd11*Zd14-2D0*SW/CW*Zd41*Zd44
     &      +3D0*CW/SW*Zd11*Zd14)
      AAABR(250) = EE/SW*Sqrt2*Vus*Zd22*Zu11
      AAABR(251) = EE/SW*Sqrt2*Vcs*Zd22*Zu22
      AAABR(252) = EE/SW*Sqrt2*Zd33*Zu33
      AAABR(253) = EE/SW*Sqrt2*Vcd*Zd11*Zu22
      AAABR(254) = EE/SW*Sqrt2*Vud*Zd11*Zu11
      AAABR(255) = EE/SW*Sqrt2*Vcd*Zd11*Zu25
      AAABR(256) = EE/SW*Sqrt2*Vcs*Zd22*Zu25
      AAABR(257) = EE/SW*Sqrt2*Vud*Zd11*Zu14
      AAABR(258) = EE/SW*Sqrt2*Zd33*Zu36
      AAABR(259) = EE/SW*Sqrt2*Vus*Zd22*Zu14
      AAABR(260) = EE*(SW/CW*Zd25**2-2D0*SW/CW*Zd55**2+3D0
     &      *CW/SW*Zd25**2)
      AAABR(261) = EE*(SW/CW*Zd36**2-2D0*SW/CW*Zd66**2+3D0
     &      *CW/SW*Zd36**2)
      AAABR(262) = EE*(SW/CW*Zd14**2-2D0*SW/CW*Zd44**2+3D0
     &      *CW/SW*Zd14**2)
      AAABR(263) = EE/SW*Sqrt2*Vcs*Zd25*Zu22
      AAABR(264) = EE/SW*Sqrt2*Vud*Zd14*Zu11
      AAABR(265) = EE/SW*Sqrt2*Vcd*Zd14*Zu22
      AAABR(266) = EE/SW*Sqrt2*Vus*Zd25*Zu11
      AAABR(267) = EE/SW*Sqrt2*Zd36*Zu33
      AAABR(268) = EE/SW*Sqrt2*Vcs*Zd25*Zu25
      AAABR(269) = EE/SW*Sqrt2*Vcd*Zd14*Zu25
      AAABR(270) = EE/SW*Sqrt2*Zd36*Zu36
      AAABR(271) = EE/SW*Sqrt2*Vus*Zd25*Zu14
      AAABR(272) = EE/SW*Sqrt2*Vud*Zd14*Zu14
      AAABR(273) = EE*(SW/CW*Zu11**2+4D0*SW/CW*Zu41**2-3D0
     &      *CW/SW*Zu11**2)
      AAABR(274) = EE*(SW/CW*Zu33**2+4D0*SW/CW*Zu63**2-3D0
     &      *CW/SW*Zu33**2)
      AAABR(275) = EE*(SW/CW*Zu22**2+4D0*SW/CW*Zu52**2-3D0
     &      *CW/SW*Zu22**2)
      AAABR(276) = EE*(SW/CW*Zu11*Zu14+4D0*SW/CW*Zu41*Zu44
     &      -3D0*CW/SW*Zu11*Zu14)
      AAABR(277) = EE*(SW/CW*Zu22*Zu25+4D0*SW/CW*Zu52*Zu55
     &      -3D0*CW/SW*Zu22*Zu25)
      AAABR(278) = EE*(SW/CW*Zu33*Zu36+4D0*SW/CW*Zu63*Zu66
     &      -3D0*CW/SW*Zu33*Zu36)
      AAABR(279) = EE*(SW/CW*Zu36**2+4D0*SW/CW*Zu66**2-3D0
     &      *CW/SW*Zu36**2)
      AAABR(280) = EE*(SW/CW*Zu14**2+4D0*SW/CW*Zu44**2-3D0
     &      *CW/SW*Zu14**2)
      AAABR(281) = EE*(SW/CW*Zu25**2+4D0*SW/CW*Zu55**2-3D0
     &      *CW/SW*Zu25**2)
      AAABR(282) = EE/SW*(Za13*cb-Za23*sb)
      AAABR(283) = EE/SW*(Za12*cb-Za22*sb)
      AAABR(284) = EE/SW*(Za11*cb-Za21*sb)
      AAABR(285) = EE/SW*(Zh11*cb+Zh21*sb)
      AAABR(286) = EE/SW*(Zh12*cb+Zh22*sb)
      AAABR(287) = EE/SW*(Zh13*cb+Zh23*sb)
      AAABR(288) = EE*(SW/CW*Za13*Zh12-SW/CW*Za23*Zh22+CW
     &      /SW*Za13*Zh12-CW/SW*Za23*Zh22)
      AAABR(289) = EE*(SW/CW*Za13*Zh13-SW/CW*Za23*Zh23+CW
     &      /SW*Za13*Zh13-CW/SW*Za23*Zh23)
      AAABR(290) = EE*(SW/CW*Za13*Zh11-SW/CW*Za23*Zh21+CW
     &      /SW*Za13*Zh11-CW/SW*Za23*Zh21)
      AAABR(291) = EE*(SW/CW*Za11*Zh11-SW/CW*Za21*Zh21+CW
     &      /SW*Za11*Zh11-CW/SW*Za21*Zh21)
      AAABR(292) = EE*(SW/CW*Za11*Zh13-SW/CW*Za21*Zh23+CW
     &      /SW*Za11*Zh13-CW/SW*Za21*Zh23)
      AAABR(293) = EE*(SW/CW*Za12*Zh12-SW/CW*Za22*Zh22+CW
     &      /SW*Za12*Zh12-CW/SW*Za22*Zh22)
      AAABR(294) = EE*(SW/CW*Za12*Zh13-SW/CW*Za22*Zh23+CW
     &      /SW*Za12*Zh13-CW/SW*Za22*Zh23)
      AAABR(295) = EE*(SW/CW*Za12*Zh11-SW/CW*Za22*Zh21+CW
     &      /SW*Za12*Zh11-CW/SW*Za22*Zh21)
      AAABR(296) = EE*(SW/CW*Za11*Zh12-SW/CW*Za21*Zh22+CW
     &      /SW*Za11*Zh12-CW/SW*Za21*Zh22)
      AAABR(297) = EE/CW*MW*SW
      AAABR(298) = EE*MW*(2D0*SW*Zh11*cb+2D0*SW*Zh21*sb+SW
     &      **3/CW**2*Zh11*cb+SW**3/CW**2*Zh21*sb+CW**2/SW
     &      *Zh11*cb+CW**2/SW*Zh21*sb)
      AAABR(299) = EE*MW*(2D0*SW*Zh12*cb+2D0*SW*Zh22*sb+SW
     &      **3/CW**2*Zh12*cb+SW**3/CW**2*Zh22*sb+CW**2/SW
     &      *Zh12*cb+CW**2/SW*Zh22*sb)
      AAABR(300) = EE*MW*(2D0*SW*Zh13*cb+2D0*SW*Zh23*sb+SW
     &      **3/CW**2*Zh13*cb+SW**3/CW**2*Zh23*sb+CW**2/SW
     &      *Zh13*cb+CW**2/SW*Zh23*sb)
      AAABR(301) = Sqrt2*(Za33*Zm22*Zp22*hL-EE/SW*Za13*Zm22
     &      *Zp12-EE/SW*Za23*Zm12*Zp22)
      AAABR(302) = Sqrt2*(Za33*Zm21*Zp22*hL-EE/SW*Za13*Zm21
     &      *Zp12-EE/SW*Za23*Zm11*Zp22)
      AAABR(303) = Sqrt2*(Za33*Zm22*Zp21*hL-EE/SW*Za13*Zm22
     &      *Zp11-EE/SW*Za23*Zm12*Zp21)
      AAABR(304) = Sqrt2*(Za33*Zm21*Zp21*hL-EE/SW*Za13*Zm21
     &      *Zp11-EE/SW*Za23*Zm11*Zp21)
      AAABR(305) = Sqrt2*(Za32*Zm21*Zp21*hL-EE/SW*Za12*Zm21
     &      *Zp11-EE/SW*Za22*Zm11*Zp21)
      AAABR(306) = Sqrt2*(Za32*Zm22*Zp21*hL-EE/SW*Za12*Zm22
     &      *Zp11-EE/SW*Za22*Zm12*Zp21)
      AAABR(307) = Sqrt2*(Za32*Zm21*Zp22*hL-EE/SW*Za12*Zm21
     &      *Zp12-EE/SW*Za22*Zm11*Zp22)
      AAABR(308) = Sqrt2*(Za31*Zm22*Zp22*hL-EE/SW*Za11*Zm22
     &      *Zp12-EE/SW*Za21*Zm12*Zp22)
      AAABR(309) = Sqrt2*(Za32*Zm22*Zp22*hL-EE/SW*Za12*Zm22
     &      *Zp12-EE/SW*Za22*Zm12*Zp22)
      AAABR(310) = Sqrt2*(Za31*Zm21*Zp22*hL-EE/SW*Za11*Zm21
     &      *Zp12-EE/SW*Za21*Zm11*Zp22)
      AAABR(311) = Sqrt2*(Za31*Zm22*Zp21*hL-EE/SW*Za11*Zm22
     &      *Zp11-EE/SW*Za21*Zm12*Zp21)
      AAABR(312) = Sqrt2*(Za31*Zm21*Zp21*hL-EE/SW*Za11*Zm21
     &      *Zp11-EE/SW*Za21*Zm11*Zp21)
      AAABR(313) = Sqrt2*(Zh31*Zm22*Zp22*hL+EE/SW*Zh11*Zm22
     &      *Zp12+EE/SW*Zh21*Zm12*Zp22)
      AAABR(314) = Sqrt2*(Zh32*Zm22*Zp22*hL+EE/SW*Zh12*Zm22
     &      *Zp12+EE/SW*Zh22*Zm12*Zp22)
      AAABR(315) = Sqrt2*(Zh31*Zm21*Zp21*hL+EE/SW*Zh11*Zm21
     &      *Zp11+EE/SW*Zh21*Zm11*Zp21)
      AAABR(316) = Sqrt2*(Zh32*Zm22*Zp21*hL+EE/SW*Zh12*Zm22
     &      *Zp11+EE/SW*Zh22*Zm12*Zp21)
      AAABR(317) = Sqrt2*(Zh32*Zm21*Zp22*hL+EE/SW*Zh12*Zm21
     &      *Zp12+EE/SW*Zh22*Zm11*Zp22)
      AAABR(318) = Sqrt2*(Zh33*Zm22*Zp22*hL+EE/SW*Zh13*Zm22
     &      *Zp12+EE/SW*Zh23*Zm12*Zp22)
      AAABR(319) = Sqrt2*(Zh32*Zm21*Zp21*hL+EE/SW*Zh12*Zm21
     &      *Zp11+EE/SW*Zh22*Zm11*Zp21)
      AAABR(320) = Sqrt2*(Zh31*Zm22*Zp21*hL+EE/SW*Zh11*Zm22
     &      *Zp11+EE/SW*Zh21*Zm12*Zp21)
      AAABR(321) = Sqrt2*(Zh31*Zm21*Zp22*hL+EE/SW*Zh11*Zm21
     &      *Zp12+EE/SW*Zh21*Zm11*Zp22)
      AAABR(322) = Sqrt2*(Zh33*Zm22*Zp21*hL+EE/SW*Zh13*Zm22
     &      *Zp11+EE/SW*Zh23*Zm12*Zp21)
      AAABR(323) = Sqrt2*(Zh33*Zm21*Zp22*hL+EE/SW*Zh13*Zm21
     &      *Zp12+EE/SW*Zh23*Zm11*Zp22)
      AAABR(324) = Sqrt2*(Zh33*Zm21*Zp21*hL+EE/SW*Zh13*Zm21
     &      *Zp11+EE/SW*Zh23*Zm11*Zp21)
      AAABR(325) = EE/SW*(2D0*Zl22*Zm12-Mm/MW*Sqrt2*Zl52
     &      *Zm22/cb)
      AAABR(326) = EE/SW*(2D0*Zl11*Zm12-Me/MW*Sqrt2*Zl41
     &      *Zm22/cb)
      AAABR(327) = EE/SW*(2D0*Zl11*Zm11-Me/MW*Sqrt2*Zl41
     &      *Zm21/cb)
      AAABR(328) = EE/SW*(2D0*Zl22*Zm11-Mm/MW*Sqrt2*Zl52
     &      *Zm21/cb)
      AAABR(329) = EE/SW*(2D0*Zl33*Zm11-Ml/MW*Sqrt2*Zl63
     &      *Zm21/cb)
      AAABR(330) = EE/SW*(2D0*Zl33*Zm12-Ml/MW*Sqrt2*Zl63
     &      *Zm22/cb)
      AAABR(331) = EE/SW*(2D0*Zl25*Zm12-Mm/MW*Sqrt2*Zl55
     &      *Zm22/cb)
      AAABR(332) = EE/SW*(2D0*Zl25*Zm11-Mm/MW*Sqrt2*Zl55
     &      *Zm21/cb)
      AAABR(333) = EE/SW*(2D0*Zl14*Zm12-Me/MW*Sqrt2*Zl44
     &      *Zm22/cb)
      AAABR(334) = EE/SW*(2D0*Zl36*Zm11-Ml/MW*Sqrt2*Zl66
     &      *Zm21/cb)
      AAABR(335) = EE/SW*(2D0*Zl14*Zm11-Me/MW*Sqrt2*Zl44
     &      *Zm21/cb)
      AAABR(336) = EE/SW*(2D0*Zl36*Zm12-Ml/MW*Sqrt2*Zl66
     &      *Zm22/cb)
      AAABC(1) = 2D0*Zm21*Zn52*cb*hL+2D0*EE/SW*Zm11*Zn32
     &      *sb-EE/CW*Sqrt2*Zm21*Zn12*sb-EE/SW*Sqrt2*Zm21
     &      *Zn22*sb
      AAABC(2) = 2D0*dconjg(Zn52)*Zp21*hL*sb+2D0*EE/SW*dconjg(Zn42
     &      )*Zp11*cb+EE/CW*Sqrt2*Zn12*Zp21*cb+EE/SW*Sqrt2
     &      *Zn22*Zp21*cb
      AAABC(3) = 2D0*Zm21*Zn54*cb*hL+2D0*EE/SW*Zm11*Zn34
     &      *sb-EE/CW*Sqrt2*Zm21*Zn14*sb-EE/SW*Sqrt2*Zm21
     &      *Zn24*sb
      AAABC(4) = 2D0*dconjg(Zn54)*Zp21*hL*sb+2D0*EE/SW*dconjg(Zn44
     &      )*Zp11*cb+EE/CW*Sqrt2*Zn14*Zp21*cb+EE/SW*Sqrt2
     &      *Zn24*Zp21*cb
      AAABC(5) = 2D0*Zm21*Zn53*cb*hL+2D0*EE/SW*Zm11*Zn33
     &      *sb-EE/CW*Sqrt2*Zm21*Zn13*sb-EE/SW*Sqrt2*Zm21
     &      *Zn23*sb
      AAABC(6) = 2D0*dconjg(Zn53)*Zp21*hL*sb+2D0*EE/SW*dconjg(Zn43
     &      )*Zp11*cb+EE/CW*Sqrt2*Zn13*Zp21*cb+EE/SW*Sqrt2
     &      *Zn23*Zp21*cb
      AAABC(7) = 2D0*Zm21*Zn51*cb*hL+2D0*EE/SW*Zm11*Zn31
     &      *sb-EE/CW*Sqrt2*Zm21*Zn11*sb-EE/SW*Sqrt2*Zm21
     &      *Zn21*sb
      AAABC(8) = 2D0*dconjg(Zn51)*Zp21*hL*sb+2D0*EE/SW*dconjg(Zn41
     &      )*Zp11*cb+EE/CW*Sqrt2*Zn11*Zp21*cb+EE/SW*Sqrt2
     &      *Zn21*Zp21*cb
      AAABC(9) = 2D0*Zm22*Zn53*cb*hL+2D0*EE/SW*Zm12*Zn33
     &      *sb-EE/CW*Sqrt2*Zm22*Zn13*sb-EE/SW*Sqrt2*Zm22
     &      *Zn23*sb
      AAABC(10) = 2D0*dconjg(Zn53)*Zp22*hL*sb+2D0*EE/SW*dconjg(Zn43
     &      )*Zp12*cb+EE/CW*Sqrt2*Zn13*Zp22*cb+EE/SW*Sqrt2
     &      *Zn23*Zp22*cb
      AAABC(11) = 2D0*Zm22*Zn52*cb*hL+2D0*EE/SW*Zm12*Zn32
     &      *sb-EE/CW*Sqrt2*Zm22*Zn12*sb-EE/SW*Sqrt2*Zm22
     &      *Zn22*sb
      AAABC(12) = 2D0*dconjg(Zn52)*Zp22*hL*sb+2D0*EE/SW*dconjg(Zn42
     &      )*Zp12*cb+EE/CW*Sqrt2*Zn12*Zp22*cb+EE/SW*Sqrt2
     &      *Zn22*Zp22*cb
      AAABC(13) = 2D0*Zm22*Zn51*cb*hL+2D0*EE/SW*Zm12*Zn31
     &      *sb-EE/CW*Sqrt2*Zm22*Zn11*sb-EE/SW*Sqrt2*Zm22
     &      *Zn21*sb
      AAABC(14) = 2D0*dconjg(Zn51)*Zp22*hL*sb+2D0*EE/SW*dconjg(Zn41
     &      )*Zp12*cb+EE/CW*Sqrt2*Zn11*Zp22*cb+EE/SW*Sqrt2
     &      *Zn21*Zp22*cb
      AAABC(15) = 2D0*Zm22*Zn54*cb*hL+2D0*EE/SW*Zm12*Zn34
     &      *sb-EE/CW*Sqrt2*Zm22*Zn14*sb-EE/SW*Sqrt2*Zm22
     &      *Zn24*sb
      AAABC(16) = 2D0*dconjg(Zn54)*Zp22*hL*sb+2D0*EE/SW*dconjg(Zn44
     &      )*Zp12*cb+EE/CW*Sqrt2*Zn14*Zp22*cb+EE/SW*Sqrt2
     &      *Zn24*Zp22*cb
      AAABC(17) = 2D0*Zm22*Zn53*hL*sb-2D0*EE/SW*Zm12*Zn33
     &      *cb+EE/CW*Sqrt2*Zm22*Zn13*cb+EE/SW*Sqrt2*Zm22
     &      *Zn23*cb
      AAABC(18) = 2D0*dconjg(Zn53)*Zp22*cb*hL-2D0*EE/SW*dconjg(Zn43
     &      )*Zp12*sb-EE/CW*Sqrt2*Zn13*Zp22*sb-EE/SW*Sqrt2
     &      *Zn23*Zp22*sb
      AAABC(19) = 2D0*Zm22*Zn52*hL*sb-2D0*EE/SW*Zm12*Zn32
     &      *cb+EE/CW*Sqrt2*Zm22*Zn12*cb+EE/SW*Sqrt2*Zm22
     &      *Zn22*cb
      AAABC(20) = 2D0*dconjg(Zn52)*Zp22*cb*hL-2D0*EE/SW*dconjg(Zn42
     &      )*Zp12*sb-EE/CW*Sqrt2*Zn12*Zp22*sb-EE/SW*Sqrt2
     &      *Zn22*Zp22*sb
      AAABC(21) = 2D0*Zm21*Zn54*hL*sb-2D0*EE/SW*Zm11*Zn34
     &      *cb+EE/CW*Sqrt2*Zm21*Zn14*cb+EE/SW*Sqrt2*Zm21
     &      *Zn24*cb
      AAABC(22) = 2D0*dconjg(Zn54)*Zp21*cb*hL-2D0*EE/SW*dconjg(Zn44
     &      )*Zp11*sb-EE/CW*Sqrt2*Zn14*Zp21*sb-EE/SW*Sqrt2
     &      *Zn24*Zp21*sb
      AAABC(23) = 2D0*Zm21*Zn53*hL*sb-2D0*EE/SW*Zm11*Zn33
     &      *cb+EE/CW*Sqrt2*Zm21*Zn13*cb+EE/SW*Sqrt2*Zm21
     &      *Zn23*cb
      AAABC(24) = 2D0*dconjg(Zn53)*Zp21*cb*hL-2D0*EE/SW*dconjg(Zn43
     &      )*Zp11*sb-EE/CW*Sqrt2*Zn13*Zp21*sb-EE/SW*Sqrt2
     &      *Zn23*Zp21*sb
      AAABC(25) = 2D0*Zm22*Zn51*hL*sb-2D0*EE/SW*Zm12*Zn31
     &      *cb+EE/CW*Sqrt2*Zm22*Zn11*cb+EE/SW*Sqrt2*Zm22
     &      *Zn21*cb
      AAABC(26) = 2D0*dconjg(Zn51)*Zp22*cb*hL-2D0*EE/SW*dconjg(Zn41
     &      )*Zp12*sb-EE/CW*Sqrt2*Zn11*Zp22*sb-EE/SW*Sqrt2
     &      *Zn21*Zp22*sb
      AAABC(27) = 2D0*Zm21*Zn51*hL*sb-2D0*EE/SW*Zm11*Zn31
     &      *cb+EE/CW*Sqrt2*Zm21*Zn11*cb+EE/SW*Sqrt2*Zm21
     &      *Zn21*cb
      AAABC(28) = 2D0*dconjg(Zn51)*Zp21*cb*hL-2D0*EE/SW*dconjg(Zn41
     &      )*Zp11*sb-EE/CW*Sqrt2*Zn11*Zp21*sb-EE/SW*Sqrt2
     &      *Zn21*Zp21*sb
      AAABC(29) = 2D0*Zm21*Zn52*hL*sb-2D0*EE/SW*Zm11*Zn32
     &      *cb+EE/CW*Sqrt2*Zm21*Zn12*cb+EE/SW*Sqrt2*Zm21
     &      *Zn22*cb
      AAABC(30) = 2D0*dconjg(Zn52)*Zp21*cb*hL-2D0*EE/SW*dconjg(Zn42
     &      )*Zp11*sb-EE/CW*Sqrt2*Zn12*Zp21*sb-EE/SW*Sqrt2
     &      *Zn22*Zp21*sb
      AAABC(31) = 2D0*Zm22*Zn54*hL*sb-2D0*EE/SW*Zm12*Zn34
     &      *cb+EE/CW*Sqrt2*Zm22*Zn14*cb+EE/SW*Sqrt2*Zm22
     &      *Zn24*cb
      AAABC(32) = 2D0*dconjg(Zn54)*Zp22*cb*hL-2D0*EE/SW*dconjg(Zn44
     &      )*Zp12*sb-EE/CW*Sqrt2*Zn14*Zp22*sb-EE/SW*Sqrt2
     &      *Zn24*Zp22*sb
      AAABR(337) = EE/SW*Vcd*Zd11*Zm11
      AAABR(338) = EE/MW*Mc/SW*Sqrt2*Vcd*Zd11*Zp21/sb
      AAABR(339) = EE/SW*(2D0*Zd33*Zm11-Mb/MW*Sqrt2*Zd63
     &      *Zm21/cb)
      AAABR(340) = EE/MW*Mt/SW*Sqrt2*Zd33*Zp21/sb
      AAABR(341) = EE/SW*Vus*(2D0*Zd22*Zm12-Ms/MW*Sqrt2*Zd52
     &      *Zm22/cb)
      AAABR(342) = EE/SW*Vud*Zd11*Zm12
      AAABR(343) = EE/SW*Vcs*(2D0*Zd22*Zm11-Ms/MW*Sqrt2*Zd52
     &      *Zm21/cb)
      AAABR(344) = EE/MW*Mc/SW*Sqrt2*Vcs*Zd22*Zp21/sb
      AAABR(345) = EE/SW*Vud*Zd11*Zm11
      AAABR(346) = EE/SW*Vcs*(2D0*Zd22*Zm12-Ms/MW*Sqrt2*Zd52
     &      *Zm22/cb)
      AAABR(347) = EE/MW*Mc/SW*Sqrt2*Vcs*Zd22*Zp22/sb
      AAABR(348) = EE/SW*Vus*(2D0*Zd22*Zm11-Ms/MW*Sqrt2*Zd52
     &      *Zm21/cb)
      AAABR(349) = EE/SW*(2D0*Zd33*Zm12-Mb/MW*Sqrt2*Zd63
     &      *Zm22/cb)
      AAABR(350) = EE/MW*Mt/SW*Sqrt2*Zd33*Zp22/sb
      AAABR(351) = EE/SW*Vcd*Zd11*Zm12
      AAABR(352) = EE/MW*Mc/SW*Sqrt2*Vcd*Zd11*Zp22/sb
      AAABR(353) = EE/SW*Vud*Zd14*Zm12
      AAABR(354) = EE/SW*Vcd*Zd14*Zm11
      AAABR(355) = EE/MW*Mc/SW*Sqrt2*Vcd*Zd14*Zp21/sb
      AAABR(356) = EE/SW*(2D0*Zd36*Zm12-Mb/MW*Sqrt2*Zd66
     &      *Zm22/cb)
      AAABR(357) = EE/MW*Mt/SW*Sqrt2*Zd36*Zp22/sb
      AAABR(358) = EE/SW*(2D0*Zd36*Zm11-Mb/MW*Sqrt2*Zd66
     &      *Zm21/cb)
      AAABR(359) = EE/MW*Mt/SW*Sqrt2*Zd36*Zp21/sb
      AAABR(360) = EE/SW*Vus*(2D0*Zd25*Zm11-Ms/MW*Sqrt2*Zd55
     &      *Zm21/cb)
      AAABR(361) = EE/SW*Vus*(2D0*Zd25*Zm12-Ms/MW*Sqrt2*Zd55
     &      *Zm22/cb)
      AAABR(362) = EE/SW*Vud*Zd14*Zm11
      AAABR(363) = EE/SW*Vcs*(2D0*Zd25*Zm11-Ms/MW*Sqrt2*Zd55
     &      *Zm21/cb)
      AAABR(364) = EE/MW*Mc/SW*Sqrt2*Vcs*Zd25*Zp21/sb
      AAABR(365) = EE/SW*Vcs*(2D0*Zd25*Zm12-Ms/MW*Sqrt2*Zd55
     &      *Zm22/cb)
      AAABR(366) = EE/MW*Mc/SW*Sqrt2*Vcs*Zd25*Zp22/sb
      AAABR(367) = EE/SW*Vcd*Zd14*Zm12
      AAABR(368) = EE/MW*Mc/SW*Sqrt2*Vcd*Zd14*Zp22/sb
      AAABC(33) = 2D0*Zm21*Zn55*cb*hL+2D0*EE/SW*Zm11*Zn35
     &      *sb-EE/CW*Sqrt2*Zm21*Zn15*sb-EE/SW*Sqrt2*Zm21
     &      *Zn25*sb
      AAABC(34) = 2D0*dconjg(Zn55)*Zp21*hL*sb+2D0*EE/SW*dconjg(Zn45
     &      )*Zp11*cb+EE/CW*Sqrt2*Zn15*Zp21*cb+EE/SW*Sqrt2
     &      *Zn25*Zp21*cb
      AAABC(35) = 2D0*Zm22*Zn55*cb*hL+2D0*EE/SW*Zm12*Zn35
     &      *sb-EE/CW*Sqrt2*Zm22*Zn15*sb-EE/SW*Sqrt2*Zm22
     &      *Zn25*sb
      AAABC(36) = 2D0*dconjg(Zn55)*Zp22*hL*sb+2D0*EE/SW*dconjg(Zn45
     &      )*Zp12*cb+EE/CW*Sqrt2*Zn15*Zp22*cb+EE/SW*Sqrt2
     &      *Zn25*Zp22*cb
      AAABC(37) = 2D0*Zm22*Zn55*hL*sb-2D0*EE/SW*Zm12*Zn35
     &      *cb+EE/CW*Sqrt2*Zm22*Zn15*cb+EE/SW*Sqrt2*Zm22
     &      *Zn25*cb
      AAABC(38) = 2D0*dconjg(Zn55)*Zp22*cb*hL-2D0*EE/SW*dconjg(Zn45
     &      )*Zp12*sb-EE/CW*Sqrt2*Zn15*Zp22*sb-EE/SW*Sqrt2
     &      *Zn25*Zp22*sb
      AAABC(39) = 2D0*Zm21*Zn55*hL*sb-2D0*EE/SW*Zm11*Zn35
     &      *cb+EE/CW*Sqrt2*Zm21*Zn15*cb+EE/SW*Sqrt2*Zm21
     &      *Zn25*cb
      AAABC(40) = 2D0*dconjg(Zn55)*Zp21*cb*hL-2D0*EE/SW*dconjg(Zn45
     &      )*Zp11*sb-EE/CW*Sqrt2*Zn15*Zp21*sb-EE/SW*Sqrt2
     &      *Zn25*Zp21*sb
      AAABR(369) = EE/MW*Mm/SW*Za13/cb
      AAABR(370) = EE/MW*Me/SW*Za13/cb
      AAABR(371) = EE/MW*Ml/SW*Za13/cb
      AAABR(372) = EE/MW*Mm/SW*Za12/cb
      AAABR(373) = EE/MW*Ml/SW*Za11/cb
      AAABR(374) = EE/MW*Me/SW*Za11/cb
      AAABR(375) = EE/MW*Mm/SW*Za11/cb
      AAABR(376) = EE/MW*Me/SW*Za12/cb
      AAABR(377) = EE/MW*Ml/SW*Za12/cb
      AAABR(378) = EE/MW*Mm/SW*Zh11/cb
      AAABR(379) = EE/MW*Ml/SW*Zh13/cb
      AAABR(380) = EE/MW*Mm/SW*Zh13/cb
      AAABR(381) = EE/MW*Me/SW*Zh13/cb
      AAABR(382) = EE/MW*Mm/SW*Zh12/cb
      AAABR(383) = EE/MW*Me/SW*Zh11/cb
      AAABR(384) = EE/MW*Me/SW*Zh12/cb
      AAABR(385) = EE/MW*Ml/SW*Zh11/cb
      AAABR(386) = EE/MW*Ml/SW*Zh12/cb
      AAABR(387) = EE/MW*Me/SW*Sqrt2*tb
      AAABR(388) = EE/MW*Mm/SW*Sqrt2*tb
      AAABR(389) = EE/MW*Ml/SW*Sqrt2*tb
      AAABR(390) = EE/MW*Ml/SW*Sqrt2
      AAABR(391) = EE/MW*Me/SW*Sqrt2
      AAABR(392) = EE/MW*Mm/SW*Sqrt2
      AAABC(41) = EE*Sqrt2*(Zl11/CW*Zn12+Zl11/SW*Zn22-Me
     &      /MW/SW*Zl41*dconjg(Zn32)/cb)
      AAABC(42) = EE*Sqrt2*(2D0*Zl41/CW*Zn12+Me/MW/SW*Zl11
     &      *Zn32/cb)
      AAABC(43) = EE*Sqrt2*(Zl11/CW*Zn11+Zl11/SW*Zn21-Me
     &      /MW/SW*Zl41*dconjg(Zn31)/cb)
      AAABC(44) = EE*Sqrt2*(2D0*Zl41/CW*Zn11+Me/MW/SW*Zl11
     &      *Zn31/cb)
      AAABC(45) = EE*Sqrt2*(Zl22/CW*Zn14+Zl22/SW*Zn24-Mm
     &      /MW/SW*Zl52*dconjg(Zn34)/cb)
      AAABC(46) = EE*Sqrt2*(2D0*Zl52/CW*Zn14+Mm/MW/SW*Zl22
     &      *Zn34/cb)
      AAABC(47) = EE*Sqrt2*(Zl11/CW*Zn13+Zl11/SW*Zn23-Me
     &      /MW/SW*Zl41*dconjg(Zn33)/cb)
      AAABC(48) = EE*Sqrt2*(2D0*Zl41/CW*Zn13+Me/MW/SW*Zl11
     &      *Zn33/cb)
      AAABC(49) = EE*Sqrt2*(Zl33/CW*Zn14+Zl33/SW*Zn24-Ml
     &      /MW/SW*Zl63*dconjg(Zn34)/cb)
      AAABC(50) = EE*Sqrt2*(2D0*Zl63/CW*Zn14+Ml/MW/SW*Zl33
     &      *Zn34/cb)
      AAABC(51) = EE*Sqrt2*(Zl33/CW*Zn11+Zl33/SW*Zn21-Ml
     &      /MW/SW*Zl63*dconjg(Zn31)/cb)
      AAABC(52) = EE*Sqrt2*(2D0*Zl63/CW*Zn11+Ml/MW/SW*Zl33
     &      *Zn31/cb)
      AAABC(53) = EE*Sqrt2*(Zl22/CW*Zn13+Zl22/SW*Zn23-Mm
     &      /MW/SW*Zl52*dconjg(Zn33)/cb)
      AAABC(54) = EE*Sqrt2*(2D0*Zl52/CW*Zn13+Mm/MW/SW*Zl22
     &      *Zn33/cb)
      AAABC(55) = EE*Sqrt2*(Zl33/CW*Zn12+Zl33/SW*Zn22-Ml
     &      /MW/SW*Zl63*dconjg(Zn32)/cb)
      AAABC(56) = EE*Sqrt2*(2D0*Zl63/CW*Zn12+Ml/MW/SW*Zl33
     &      *Zn32/cb)
      AAABC(57) = EE*Sqrt2*(Zl33/CW*Zn13+Zl33/SW*Zn23-Ml
     &      /MW/SW*Zl63*dconjg(Zn33)/cb)
      AAABC(58) = EE*Sqrt2*(2D0*Zl63/CW*Zn13+Ml/MW/SW*Zl33
     &      *Zn33/cb)
      AAABC(59) = EE*Sqrt2*(Zl11/CW*Zn14+Zl11/SW*Zn24-Me
     &      /MW/SW*Zl41*dconjg(Zn34)/cb)
      AAABC(60) = EE*Sqrt2*(2D0*Zl41/CW*Zn14+Me/MW/SW*Zl11
     &      *Zn34/cb)
      AAABC(61) = EE*Sqrt2*(Zl22/CW*Zn12+Zl22/SW*Zn22-Mm
     &      /MW/SW*Zl52*dconjg(Zn32)/cb)
      AAABC(62) = EE*Sqrt2*(2D0*Zl52/CW*Zn12+Mm/MW/SW*Zl22
     &      *Zn32/cb)
      AAABC(63) = EE*Sqrt2*(Zl22/CW*Zn11+Zl22/SW*Zn21-Mm
     &      /MW/SW*Zl52*dconjg(Zn31)/cb)
      AAABC(64) = EE*Sqrt2*(2D0*Zl52/CW*Zn11+Mm/MW/SW*Zl22
     &      *Zn31/cb)
      AAABC(65) = EE*Sqrt2*(Zl14/CW*Zn13+Zl14/SW*Zn23-Me
     &      /MW/SW*Zl44*dconjg(Zn33)/cb)
      AAABC(66) = EE*Sqrt2*(2D0*Zl44/CW*Zn13+Me/MW/SW*Zl14
     &      *Zn33/cb)
      AAABC(67) = EE*Sqrt2*(Zl36/CW*Zn14+Zl36/SW*Zn24-Ml
     &      /MW/SW*Zl66*dconjg(Zn34)/cb)
      AAABC(68) = EE*Sqrt2*(2D0*Zl66/CW*Zn14+Ml/MW/SW*Zl36
     &      *Zn34/cb)
      AAABC(69) = EE*Sqrt2*(Zl36/CW*Zn13+Zl36/SW*Zn23-Ml
     &      /MW/SW*Zl66*dconjg(Zn33)/cb)
      AAABC(70) = EE*Sqrt2*(2D0*Zl66/CW*Zn13+Ml/MW/SW*Zl36
     &      *Zn33/cb)
      AAABC(71) = EE*Sqrt2*(Zl14/CW*Zn14+Zl14/SW*Zn24-Me
     &      /MW/SW*Zl44*dconjg(Zn34)/cb)
      AAABC(72) = EE*Sqrt2*(2D0*Zl44/CW*Zn14+Me/MW/SW*Zl14
     &      *Zn34/cb)
      AAABC(73) = EE*Sqrt2*(Zl25/CW*Zn11+Zl25/SW*Zn21-Mm
     &      /MW/SW*Zl55*dconjg(Zn31)/cb)
      AAABC(74) = EE*Sqrt2*(2D0*Zl55/CW*Zn11+Mm/MW/SW*Zl25
     &      *Zn31/cb)
      AAABC(75) = EE*Sqrt2*(Zl25/CW*Zn14+Zl25/SW*Zn24-Mm
     &      /MW/SW*Zl55*dconjg(Zn34)/cb)
      AAABC(76) = EE*Sqrt2*(2D0*Zl55/CW*Zn14+Mm/MW/SW*Zl25
     &      *Zn34/cb)
      AAABC(77) = EE*Sqrt2*(Zl25/CW*Zn13+Zl25/SW*Zn23-Mm
     &      /MW/SW*Zl55*dconjg(Zn33)/cb)
      AAABC(78) = EE*Sqrt2*(2D0*Zl55/CW*Zn13+Mm/MW/SW*Zl25
     &      *Zn33/cb)
      AAABC(79) = EE*Sqrt2*(Zl25/CW*Zn12+Zl25/SW*Zn22-Mm
     &      /MW/SW*Zl55*dconjg(Zn32)/cb)
      AAABC(80) = EE*Sqrt2*(2D0*Zl55/CW*Zn12+Mm/MW/SW*Zl25
     &      *Zn32/cb)
      AAABC(81) = EE*Sqrt2*(Zl14/CW*Zn12+Zl14/SW*Zn22-Me
     &      /MW/SW*Zl44*dconjg(Zn32)/cb)
      AAABC(82) = EE*Sqrt2*(2D0*Zl44/CW*Zn12+Me/MW/SW*Zl14
     &      *Zn32/cb)
      AAABC(83) = EE*Sqrt2*(Zl36/CW*Zn12+Zl36/SW*Zn22-Ml
     &      /MW/SW*Zl66*dconjg(Zn32)/cb)
      AAABC(84) = EE*Sqrt2*(2D0*Zl66/CW*Zn12+Ml/MW/SW*Zl36
     &      *Zn32/cb)
      AAABC(85) = EE*Sqrt2*(Zl14/CW*Zn11+Zl14/SW*Zn21-Me
     &      /MW/SW*Zl44*dconjg(Zn31)/cb)
      AAABC(86) = EE*Sqrt2*(2D0*Zl44/CW*Zn11+Me/MW/SW*Zl14
     &      *Zn31/cb)
      AAABC(87) = EE*Sqrt2*(Zl36/CW*Zn11+Zl36/SW*Zn21-Ml
     &      /MW/SW*Zl66*dconjg(Zn31)/cb)
      AAABC(88) = EE*Sqrt2*(2D0*Zl66/CW*Zn11+Ml/MW/SW*Zl36
     &      *Zn31/cb)
      AAABR(393) = EE/SW*Zp11
      AAABR(394) = EE/MW*Me/SW*Sqrt2*Zm21/cb
      AAABR(395) = EE/MW*Ml/SW*Sqrt2*Zm21/cb
      AAABR(396) = EE/MW*Mm/SW*Sqrt2*Zm21/cb
      AAABR(397) = EE/SW*Zp12
      AAABR(398) = EE/MW*Ml/SW*Sqrt2*Zm22/cb
      AAABR(399) = EE/MW*Me/SW*Sqrt2*Zm22/cb
      AAABR(400) = EE/MW*Mm/SW*Sqrt2*Zm22/cb
      AAABC(89) = EE*Sqrt2*(Zl22/CW*Zn15+Zl22/SW*Zn25-Mm
     &      /MW/SW*Zl52*dconjg(Zn35)/cb)
      AAABC(90) = EE*Sqrt2*(2D0*Zl52/CW*Zn15+Mm/MW/SW*Zl22
     &      *Zn35/cb)
      AAABC(91) = EE*Sqrt2*(Zl11/CW*Zn15+Zl11/SW*Zn25-Me
     &      /MW/SW*Zl41*dconjg(Zn35)/cb)
      AAABC(92) = EE*Sqrt2*(2D0*Zl41/CW*Zn15+Me/MW/SW*Zl11
     &      *Zn35/cb)
      AAABC(93) = EE*Sqrt2*(Zl33/CW*Zn15+Zl33/SW*Zn25-Ml
     &      /MW/SW*Zl63*dconjg(Zn35)/cb)
      AAABC(94) = EE*Sqrt2*(2D0*Zl63/CW*Zn15+Ml/MW/SW*Zl33
     &      *Zn35/cb)
      AAABC(95) = EE*Sqrt2*(Zl14/CW*Zn15+Zl14/SW*Zn25-Me
     &      /MW/SW*Zl44*dconjg(Zn35)/cb)
      AAABC(96) = EE*Sqrt2*(2D0*Zl44/CW*Zn15+Me/MW/SW*Zl14
     &      *Zn35/cb)
      AAABC(97) = EE*Sqrt2*(Zl25/CW*Zn15+Zl25/SW*Zn25-Mm
     &      /MW/SW*Zl55*dconjg(Zn35)/cb)
      AAABC(98) = EE*Sqrt2*(2D0*Zl55/CW*Zn15+Mm/MW/SW*Zl25
     &      *Zn35/cb)
      AAABC(99) = EE*Sqrt2*(Zl36/CW*Zn15+Zl36/SW*Zn25-Ml
     &      /MW/SW*Zl66*dconjg(Zn35)/cb)
      AAABC(100) = EE*Sqrt2*(2D0*Zl66/CW*Zn15+Ml/MW/SW*Zl36
     &      *Zn35/cb)
      AAABC(101) = EE*Sqrt2*(Zn14/CW-Zn24/SW)
      AAABC(102) = EE*Sqrt2*(Zn13/CW-Zn23/SW)
      AAABC(103) = EE*Sqrt2*(Zn12/CW-Zn22/SW)
      AAABC(104) = EE*Sqrt2*(Zn11/CW-Zn21/SW)
      AAABC(105) = EE*Sqrt2*(Zn15/CW-Zn25/SW)
      AAABC(106) = EE*Sqrt2*(Zd33/CW*Zn11-3D0*Zd33/SW*Zn21
     &      +3D0*Mb/MW/SW*Zd63*dconjg(Zn31)/cb)
      AAABC(107) = EE*Sqrt2*(2D0*Zd63/CW*Zn11+3D0*Mb/MW/SW
     &      *Zd33*Zn31/cb)
      AAABC(108) = EE/CW*Sqrt2*Zd41*Zn12
      AAABC(109) = EE*Sqrt2*Zd11*(Zn12/CW-3D0*Zn22/SW)
      AAABC(110) = EE/CW*Sqrt2*Zd41*Zn13
      AAABC(111) = EE*Sqrt2*Zd11*(Zn13/CW-3D0*Zn23/SW)
      AAABC(112) = EE*Sqrt2*(Zd33/CW*Zn14-3D0*Zd33/SW*Zn24
     &      +3D0*Mb/MW/SW*Zd63*dconjg(Zn34)/cb)
      AAABC(113) = EE*Sqrt2*(2D0*Zd63/CW*Zn14+3D0*Mb/MW/SW
     &      *Zd33*Zn34/cb)
      AAABC(114) = EE/CW*Sqrt2*Zd41*Zn11
      AAABC(115) = EE*Sqrt2*Zd11*(Zn11/CW-3D0*Zn21/SW)
      AAABC(116) = EE/CW*Sqrt2*Zd41*Zn14
      AAABC(117) = EE*Sqrt2*Zd11*(Zn14/CW-3D0*Zn24/SW)
      AAABC(118) = EE*Sqrt2*(2D0*Zd52/CW*Zn13+3D0*Ms/MW/SW
     &      *Zd22*Zn33/cb)
      AAABC(119) = EE*Sqrt2*(Zd22/CW*Zn13-3D0*Zd22/SW*Zn23
     &      +3D0*Ms/MW/SW*Zd52*dconjg(Zn33)/cb)
      AAABC(120) = EE*Sqrt2*(Zd33/CW*Zn13-3D0*Zd33/SW*Zn23
     &      +3D0*Mb/MW/SW*Zd63*dconjg(Zn33)/cb)
      AAABC(121) = EE*Sqrt2*(2D0*Zd63/CW*Zn13+3D0*Mb/MW/SW
     &      *Zd33*Zn33/cb)
      AAABC(122) = EE*Sqrt2*(2D0*Zd52/CW*Zn11+3D0*Ms/MW/SW
     &      *Zd22*Zn31/cb)
      AAABC(123) = EE*Sqrt2*(Zd22/CW*Zn11-3D0*Zd22/SW*Zn21
     &      +3D0*Ms/MW/SW*Zd52*dconjg(Zn31)/cb)
      AAABC(124) = EE*Sqrt2*(Zd33/CW*Zn12-3D0*Zd33/SW*Zn22
     &      +3D0*Mb/MW/SW*Zd63*dconjg(Zn32)/cb)
      AAABC(125) = EE*Sqrt2*(2D0*Zd63/CW*Zn12+3D0*Mb/MW/SW
     &      *Zd33*Zn32/cb)
      AAABC(126) = EE*Sqrt2*(2D0*Zd52/CW*Zn14+3D0*Ms/MW/SW
     &      *Zd22*Zn34/cb)
      AAABC(127) = EE*Sqrt2*(Zd22/CW*Zn14-3D0*Zd22/SW*Zn24
     &      +3D0*Ms/MW/SW*Zd52*dconjg(Zn34)/cb)
      AAABC(128) = EE*Sqrt2*(2D0*Zd52/CW*Zn12+3D0*Ms/MW/SW
     &      *Zd22*Zn32/cb)
      AAABC(129) = EE*Sqrt2*(Zd22/CW*Zn12-3D0*Zd22/SW*Zn22
     &      +3D0*Ms/MW/SW*Zd52*dconjg(Zn32)/cb)
      AAABC(130) = EE*Sqrt2*(2D0*Zd55/CW*Zn12+3D0*Ms/MW/SW
     &      *Zd25*Zn32/cb)
      AAABC(131) = EE*Sqrt2*(Zd25/CW*Zn12-3D0*Zd25/SW*Zn22
     &      +3D0*Ms/MW/SW*Zd55*dconjg(Zn32)/cb)
      AAABC(132) = EE*Sqrt2*(Zd36/CW*Zn14-3D0*Zd36/SW*Zn24
     &      +3D0*Mb/MW/SW*Zd66*dconjg(Zn34)/cb)
      AAABC(133) = EE*Sqrt2*(2D0*Zd66/CW*Zn14+3D0*Mb/MW/SW
     &      *Zd36*Zn34/cb)
      AAABC(134) = EE/CW*Sqrt2*Zd44*Zn13
      AAABC(135) = EE*Sqrt2*Zd14*(Zn13/CW-3D0*Zn23/SW)
      AAABC(136) = EE*Sqrt2*(2D0*Zd55/CW*Zn14+3D0*Ms/MW/SW
     &      *Zd25*Zn34/cb)
      AAABC(137) = EE*Sqrt2*(Zd25/CW*Zn14-3D0*Zd25/SW*Zn24
     &      +3D0*Ms/MW/SW*Zd55*dconjg(Zn34)/cb)
      AAABC(138) = EE*Sqrt2*(Zd36/CW*Zn11-3D0*Zd36/SW*Zn21
     &      +3D0*Mb/MW/SW*Zd66*dconjg(Zn31)/cb)
      AAABC(139) = EE*Sqrt2*(2D0*Zd66/CW*Zn11+3D0*Mb/MW/SW
     &      *Zd36*Zn31/cb)
      AAABC(140) = EE*Sqrt2*(Zd36/CW*Zn12-3D0*Zd36/SW*Zn22
     &      +3D0*Mb/MW/SW*Zd66*dconjg(Zn32)/cb)
      AAABC(141) = EE*Sqrt2*(2D0*Zd66/CW*Zn12+3D0*Mb/MW/SW
     &      *Zd36*Zn32/cb)
      AAABC(142) = EE/CW*Sqrt2*Zd44*Zn11
      AAABC(143) = EE*Sqrt2*Zd14*(Zn11/CW-3D0*Zn21/SW)
      AAABC(144) = EE/CW*Sqrt2*Zd44*Zn12
      AAABC(145) = EE*Sqrt2*Zd14*(Zn12/CW-3D0*Zn22/SW)
      AAABC(146) = EE*Sqrt2*(2D0*Zd55/CW*Zn11+3D0*Ms/MW/SW
     &      *Zd25*Zn31/cb)
      AAABC(147) = EE*Sqrt2*(Zd25/CW*Zn11-3D0*Zd25/SW*Zn21
     &      +3D0*Ms/MW/SW*Zd55*dconjg(Zn31)/cb)
      AAABC(148) = EE*Sqrt2*(2D0*Zd55/CW*Zn13+3D0*Ms/MW/SW
     &      *Zd25*Zn33/cb)
      AAABC(149) = EE*Sqrt2*(Zd25/CW*Zn13-3D0*Zd25/SW*Zn23
     &      +3D0*Ms/MW/SW*Zd55*dconjg(Zn33)/cb)
      AAABC(150) = EE*Sqrt2*(Zd36/CW*Zn13-3D0*Zd36/SW*Zn23
     &      +3D0*Mb/MW/SW*Zd66*dconjg(Zn33)/cb)
      AAABC(151) = EE*Sqrt2*(2D0*Zd66/CW*Zn13+3D0*Mb/MW/SW
     &      *Zd36*Zn33/cb)
      AAABC(152) = EE/CW*Sqrt2*Zd44*Zn14
      AAABC(153) = EE*Sqrt2*Zd14*(Zn14/CW-3D0*Zn24/SW)
      AAABR(401) = EE/MW*Ms/SW*Za13/cb
      AAABR(402) = EE/MW*Mb/SW*Za13/cb
      AAABR(403) = EE/MW*Ms/SW*Za12/cb
      AAABR(404) = EE/MW*Mb/SW*Za12/cb
      AAABR(405) = EE/MW*Mb/SW*Za11/cb
      AAABR(406) = EE/MW*Ms/SW*Za11/cb
      AAABR(407) = EE/MW*Ms/SW*Zh11/cb
      AAABR(408) = EE/MW*Ms/SW*Zh13/cb
      AAABR(409) = EE/MW*Ms/SW*Zh12/cb
      AAABR(410) = EE/MW*Mb/SW*Zh11/cb
      AAABR(411) = EE/MW*Mb/SW*Zh13/cb
      AAABR(412) = EE/MW*Mb/SW*Zh12/cb
      AAABR(413) = EE/MW*Mc/SW*Sqrt2*Vcs/cb*(1D0/sb-sb)
      AAABR(414) = EE/MW*Ms/SW*Sqrt2*Vcs/cb*sb
      AAABR(415) = EE/MW*Ms/SW*Sqrt2*Vus*tb
      AAABR(416) = EE/MW*Mt/SW*Sqrt2/cb*(1D0/sb-sb)
      AAABR(417) = EE/MW*Mb/SW*Sqrt2/cb*sb
      AAABR(418) = EE/MW*Mc/SW*Sqrt2*Vcd/tb
      AAABR(419) = EE/MW*Ms/SW*Sqrt2*Vus
      AAABR(420) = EE/MW*Mc/SW*Sqrt2*Vcs
      AAABR(421) = EE/MW*Ms/SW*Sqrt2*Vcs
      AAABR(422) = EE/MW*Mt/SW*Sqrt2
      AAABR(423) = EE/MW*Mb/SW*Sqrt2
      AAABR(424) = EE/MW*Mc/SW*Sqrt2*Vcd
      AAABR(425) = EE/SW*Vcs*(2D0*Zp11*Zu22-Mc/MW*Sqrt2*Zp21
     &      *Zu52/sb)
      AAABR(426) = EE/MW*Ms/SW*Sqrt2*Vcs*Zm21*Zu22/cb
      AAABR(427) = EE/SW*(2D0*Zp11*Zu33-Mt/MW*Sqrt2*Zp21
     &      *Zu63/sb)
      AAABR(428) = EE/MW*Mb/SW*Sqrt2*Zm21*Zu33/cb
      AAABR(429) = EE/SW*Vcd*(2D0*Zp11*Zu22-Mc/MW*Sqrt2*Zp21
     &      *Zu52/sb)
      AAABR(430) = EE/SW*Vud*Zp11*Zu11
      AAABR(431) = EE/SW*Vus*Zp11*Zu11
      AAABR(432) = EE/MW*Ms/SW*Sqrt2*Vus*Zm21*Zu11/cb
      AAABR(433) = EE/SW*Vus*Zp11*Zu14
      AAABR(434) = EE/MW*Ms/SW*Sqrt2*Vus*Zm21*Zu14/cb
      AAABR(435) = EE/SW*(2D0*Zp11*Zu36-Mt/MW*Sqrt2*Zp21
     &      *Zu66/sb)
      AAABR(436) = EE/MW*Mb/SW*Sqrt2*Zm21*Zu36/cb
      AAABR(437) = EE/SW*Vcs*(2D0*Zp11*Zu25-Mc/MW*Sqrt2*Zp21
     &      *Zu55/sb)
      AAABR(438) = EE/MW*Ms/SW*Sqrt2*Vcs*Zm21*Zu25/cb
      AAABR(439) = EE/SW*Vud*Zp11*Zu14
      AAABR(440) = EE/SW*Vcd*(2D0*Zp11*Zu25-Mc/MW*Sqrt2*Zp21
     &      *Zu55/sb)
      AAABR(441) = EE/SW*(2D0*Zp12*Zu33-Mt/MW*Sqrt2*Zp22
     &      *Zu63/sb)
      AAABR(442) = EE/MW*Mb/SW*Sqrt2*Zm22*Zu33/cb
      AAABR(443) = EE/SW*Vud*Zp12*Zu11
      AAABR(444) = EE/SW*Vcs*(2D0*Zp12*Zu22-Mc/MW*Sqrt2*Zp22
     &      *Zu52/sb)
      AAABR(445) = EE/MW*Ms/SW*Sqrt2*Vcs*Zm22*Zu22/cb
      AAABR(446) = EE/SW*Vcd*(2D0*Zp12*Zu22-Mc/MW*Sqrt2*Zp22
     &      *Zu52/sb)
      AAABR(447) = EE/SW*Vus*Zp12*Zu11
      AAABR(448) = EE/MW*Ms/SW*Sqrt2*Vus*Zm22*Zu11/cb
      AAABR(449) = EE/SW*(2D0*Zp12*Zu36-Mt/MW*Sqrt2*Zp22
     &      *Zu66/sb)
      AAABR(450) = EE/MW*Mb/SW*Sqrt2*Zm22*Zu36/cb
      AAABR(451) = EE/SW*Vcs*(2D0*Zp12*Zu25-Mc/MW*Sqrt2*Zp22
     &      *Zu55/sb)
      AAABR(452) = EE/MW*Ms/SW*Sqrt2*Vcs*Zm22*Zu25/cb
      AAABR(453) = EE/SW*Vud*Zp12*Zu14
      AAABR(454) = EE/SW*Vus*Zp12*Zu14
      AAABR(455) = EE/MW*Ms/SW*Sqrt2*Vus*Zm22*Zu14/cb
      AAABR(456) = EE/SW*Vcd*(2D0*Zp12*Zu25-Mc/MW*Sqrt2*Zp22
     &      *Zu55/sb)
      AAABR(457) = GG*Sqrt2*Zd41
      AAABR(458) = GG*Sqrt2*Zd11
      AAABR(459) = GG*Sqrt2*Zd52
      AAABR(460) = GG*Sqrt2*Zd22
      AAABR(461) = GG*Sqrt2*Zd33
      AAABR(462) = GG*Sqrt2*Zd63
      AAABR(463) = GG*Sqrt2*Zd44
      AAABR(464) = GG*Sqrt2*Zd14
      AAABR(465) = GG*Sqrt2*Zd36
      AAABR(466) = GG*Sqrt2*Zd66
      AAABR(467) = GG*Sqrt2*Zd55
      AAABR(468) = GG*Sqrt2*Zd25
      AAABC(154) = EE/CW*Sqrt2*Zd41*Zn15
      AAABC(155) = EE*Sqrt2*Zd11*(Zn15/CW-3D0*Zn25/SW)
      AAABC(156) = EE*Sqrt2*(2D0*Zd52/CW*Zn15+3D0*Ms/MW/SW
     &      *Zd22*Zn35/cb)
      AAABC(157) = EE*Sqrt2*(Zd22/CW*Zn15-3D0*Zd22/SW*Zn25
     &      +3D0*Ms/MW/SW*Zd52*dconjg(Zn35)/cb)
      AAABC(158) = EE*Sqrt2*(Zd33/CW*Zn15-3D0*Zd33/SW*Zn25
     &      +3D0*Mb/MW/SW*Zd63*dconjg(Zn35)/cb)
      AAABC(159) = EE*Sqrt2*(2D0*Zd63/CW*Zn15+3D0*Mb/MW/SW
     &      *Zd33*Zn35/cb)
      AAABC(160) = EE*Sqrt2*(Zd36/CW*Zn15-3D0*Zd36/SW*Zn25
     &      +3D0*Mb/MW/SW*Zd66*dconjg(Zn35)/cb)
      AAABC(161) = EE*Sqrt2*(2D0*Zd66/CW*Zn15+3D0*Mb/MW/SW
     &      *Zd36*Zn35/cb)
      AAABC(162) = EE/CW*Sqrt2*Zd44*Zn15
      AAABC(163) = EE*Sqrt2*Zd14*(Zn15/CW-3D0*Zn25/SW)
      AAABC(164) = EE*Sqrt2*(2D0*Zd55/CW*Zn15+3D0*Ms/MW/SW
     &      *Zd25*Zn35/cb)
      AAABC(165) = EE*Sqrt2*(Zd25/CW*Zn15-3D0*Zd25/SW*Zn25
     &      +3D0*Ms/MW/SW*Zd55*dconjg(Zn35)/cb)
      AAABC(166) = EE*Sqrt2*(Zn11/CW*Zu33+3D0*Zn21/SW*Zu33
     &      +3D0*Mt/MW/SW*dconjg(Zn41)*Zu63/sb)
      AAABC(167) = EE*Sqrt2*(4D0*Zn11/CW*Zu63-3D0*Mt/MW/SW
     &      *Zn41*Zu33/sb)
      AAABC(168) = EE*Sqrt2*(Zn11/CW*Zu22+3D0*Zn21/SW*Zu22
     &      +3D0*Mc/MW/SW*dconjg(Zn41)*Zu52/sb)
      AAABC(169) = EE*Sqrt2*(4D0*Zn11/CW*Zu52-3D0*Mc/MW/SW
     &      *Zn41*Zu22/sb)
      AAABC(170) = EE*Sqrt2*(Zn12/CW*Zu33+3D0*Zn22/SW*Zu33
     &      +3D0*Mt/MW/SW*dconjg(Zn42)*Zu63/sb)
      AAABC(171) = EE*Sqrt2*(4D0*Zn12/CW*Zu63-3D0*Mt/MW/SW
     &      *Zn42*Zu33/sb)
      AAABC(172) = EE*Sqrt2*Zu11*(Zn14/CW+3D0*Zn24/SW)
      AAABC(173) = EE/CW*Sqrt2*Zn14*Zu41
      AAABC(174) = EE*Sqrt2*Zu11*(Zn13/CW+3D0*Zn23/SW)
      AAABC(175) = EE/CW*Sqrt2*Zn13*Zu41
      AAABC(176) = EE*Sqrt2*Zu11*(Zn11/CW+3D0*Zn21/SW)
      AAABC(177) = EE/CW*Sqrt2*Zn11*Zu41
      AAABC(178) = EE*Sqrt2*(Zn14/CW*Zu22+3D0*Zn24/SW*Zu22
     &      +3D0*Mc/MW/SW*dconjg(Zn44)*Zu52/sb)
      AAABC(179) = EE*Sqrt2*(4D0*Zn14/CW*Zu52-3D0*Mc/MW/SW
     &      *Zn44*Zu22/sb)
      AAABC(180) = EE*Sqrt2*Zu11*(Zn12/CW+3D0*Zn22/SW)
      AAABC(181) = EE/CW*Sqrt2*Zn12*Zu41
      AAABC(182) = EE*Sqrt2*(Zn12/CW*Zu22+3D0*Zn22/SW*Zu22
     &      +3D0*Mc/MW/SW*dconjg(Zn42)*Zu52/sb)
      AAABC(183) = EE*Sqrt2*(4D0*Zn12/CW*Zu52-3D0*Mc/MW/SW
     &      *Zn42*Zu22/sb)
      AAABC(184) = EE*Sqrt2*(Zn13/CW*Zu33+3D0*Zn23/SW*Zu33
     &      +3D0*Mt/MW/SW*dconjg(Zn43)*Zu63/sb)
      AAABC(185) = EE*Sqrt2*(4D0*Zn13/CW*Zu63-3D0*Mt/MW/SW
     &      *Zn43*Zu33/sb)
      AAABC(186) = EE*Sqrt2*(Zn14/CW*Zu33+3D0*Zn24/SW*Zu33
     &      +3D0*Mt/MW/SW*dconjg(Zn44)*Zu63/sb)
      AAABC(187) = EE*Sqrt2*(4D0*Zn14/CW*Zu63-3D0*Mt/MW/SW
     &      *Zn44*Zu33/sb)
      AAABC(188) = EE*Sqrt2*(Zn13/CW*Zu22+3D0*Zn23/SW*Zu22
     &      +3D0*Mc/MW/SW*dconjg(Zn43)*Zu52/sb)
      AAABC(189) = EE*Sqrt2*(4D0*Zn13/CW*Zu52-3D0*Mc/MW/SW
     &      *Zn43*Zu22/sb)
      AAABC(190) = EE*Sqrt2*Zu14*(Zn13/CW+3D0*Zn23/SW)
      AAABC(191) = EE/CW*Sqrt2*Zn13*Zu44
      AAABC(192) = EE*Sqrt2*(Zn13/CW*Zu25+3D0*Zn23/SW*Zu25
     &      +3D0*Mc/MW/SW*dconjg(Zn43)*Zu55/sb)
      AAABC(193) = EE*Sqrt2*(4D0*Zn13/CW*Zu55-3D0*Mc/MW/SW
     &      *Zn43*Zu25/sb)
      AAABC(194) = EE*Sqrt2*Zu14*(Zn12/CW+3D0*Zn22/SW)
      AAABC(195) = EE/CW*Sqrt2*Zn12*Zu44
      AAABC(196) = EE*Sqrt2*(Zn14/CW*Zu25+3D0*Zn24/SW*Zu25
     &      +3D0*Mc/MW/SW*dconjg(Zn44)*Zu55/sb)
      AAABC(197) = EE*Sqrt2*(4D0*Zn14/CW*Zu55-3D0*Mc/MW/SW
     &      *Zn44*Zu25/sb)
      AAABC(198) = EE*Sqrt2*(Zn13/CW*Zu36+3D0*Zn23/SW*Zu36
     &      +3D0*Mt/MW/SW*dconjg(Zn43)*Zu66/sb)
      AAABC(199) = EE*Sqrt2*(4D0*Zn13/CW*Zu66-3D0*Mt/MW/SW
     &      *Zn43*Zu36/sb)
      AAABC(200) = EE*Sqrt2*Zu14*(Zn14/CW+3D0*Zn24/SW)
      AAABC(201) = EE/CW*Sqrt2*Zn14*Zu44
      AAABC(202) = EE*Sqrt2*(Zn11/CW*Zu36+3D0*Zn21/SW*Zu36
     &      +3D0*Mt/MW/SW*dconjg(Zn41)*Zu66/sb)
      AAABC(203) = EE*Sqrt2*(4D0*Zn11/CW*Zu66-3D0*Mt/MW/SW
     &      *Zn41*Zu36/sb)
      AAABC(204) = EE*Sqrt2*(Zn12/CW*Zu36+3D0*Zn22/SW*Zu36
     &      +3D0*Mt/MW/SW*dconjg(Zn42)*Zu66/sb)
      AAABC(205) = EE*Sqrt2*(4D0*Zn12/CW*Zu66-3D0*Mt/MW/SW
     &      *Zn42*Zu36/sb)
      AAABC(206) = EE*Sqrt2*(Zn11/CW*Zu25+3D0*Zn21/SW*Zu25
     &      +3D0*Mc/MW/SW*dconjg(Zn41)*Zu55/sb)
      AAABC(207) = EE*Sqrt2*(4D0*Zn11/CW*Zu55-3D0*Mc/MW/SW
     &      *Zn41*Zu25/sb)
      AAABC(208) = EE*Sqrt2*Zu14*(Zn11/CW+3D0*Zn21/SW)
      AAABC(209) = EE/CW*Sqrt2*Zn11*Zu44
      AAABC(210) = EE*Sqrt2*(Zn12/CW*Zu25+3D0*Zn22/SW*Zu25
     &      +3D0*Mc/MW/SW*dconjg(Zn42)*Zu55/sb)
      AAABC(211) = EE*Sqrt2*(4D0*Zn12/CW*Zu55-3D0*Mc/MW/SW
     &      *Zn42*Zu25/sb)
      AAABC(212) = EE*Sqrt2*(Zn14/CW*Zu36+3D0*Zn24/SW*Zu36
     &      +3D0*Mt/MW/SW*dconjg(Zn44)*Zu66/sb)
      AAABC(213) = EE*Sqrt2*(4D0*Zn14/CW*Zu66-3D0*Mt/MW/SW
     &      *Zn44*Zu36/sb)
      AAABR(469) = EE/MW*Mt/SW*Za23/sb
      AAABR(470) = EE/MW*Mc/SW*Za23/sb
      AAABR(471) = EE/MW*Mt/SW*Za21/sb
      AAABR(472) = EE/MW*Mt/SW*Za22/sb
      AAABR(473) = EE/MW*Mc/SW*Za21/sb
      AAABR(474) = EE/MW*Mc/SW*Za22/sb
      AAABR(475) = EE/MW*Mt/SW*Zh23/sb
      AAABR(476) = EE/MW*Mc/SW*Zh21/sb
      AAABR(477) = EE/MW*Mc/SW*Zh22/sb
      AAABR(478) = EE/MW*Mt/SW*Zh22/sb
      AAABR(479) = EE/MW*Mc/SW*Zh23/sb
      AAABR(480) = EE/MW*Mt/SW*Zh21/sb
      AAABR(481) = GG*Sqrt2*Zu11
      AAABR(482) = GG*Sqrt2*Zu41
      AAABR(483) = GG*Sqrt2*Zu33
      AAABR(484) = GG*Sqrt2*Zu63
      AAABR(485) = GG*Sqrt2*Zu22
      AAABR(486) = GG*Sqrt2*Zu52
      AAABR(487) = GG*Sqrt2*Zu36
      AAABR(488) = GG*Sqrt2*Zu66
      AAABR(489) = GG*Sqrt2*Zu25
      AAABR(490) = GG*Sqrt2*Zu55
      AAABR(491) = GG*Sqrt2*Zu14
      AAABR(492) = GG*Sqrt2*Zu44
      AAABC(214) = EE*Sqrt2*Zu11*(Zn15/CW+3D0*Zn25/SW)
      AAABC(215) = EE/CW*Sqrt2*Zn15*Zu41
      AAABC(216) = EE*Sqrt2*(Zn15/CW*Zu33+3D0*Zn25/SW*Zu33
     &      +3D0*Mt/MW/SW*dconjg(Zn45)*Zu63/sb)
      AAABC(217) = EE*Sqrt2*(4D0*Zn15/CW*Zu63-3D0*Mt/MW/SW
     &      *Zn45*Zu33/sb)
      AAABC(218) = EE*Sqrt2*(Zn15/CW*Zu22+3D0*Zn25/SW*Zu22
     &      +3D0*Mc/MW/SW*dconjg(Zn45)*Zu52/sb)
      AAABC(219) = EE*Sqrt2*(4D0*Zn15/CW*Zu52-3D0*Mc/MW/SW
     &      *Zn45*Zu22/sb)
      AAABC(220) = EE*Sqrt2*(Zn15/CW*Zu36+3D0*Zn25/SW*Zu36
     &      +3D0*Mt/MW/SW*dconjg(Zn45)*Zu66/sb)
      AAABC(221) = EE*Sqrt2*(4D0*Zn15/CW*Zu66-3D0*Mt/MW/SW
     &      *Zn45*Zu36/sb)
      AAABC(222) = EE*Sqrt2*(Zn15/CW*Zu25+3D0*Zn25/SW*Zu25
     &      +3D0*Mc/MW/SW*dconjg(Zn45)*Zu55/sb)
      AAABC(223) = EE*Sqrt2*(4D0*Zn15/CW*Zu55-3D0*Mc/MW/SW
     &      *Zn45*Zu25/sb)
      AAABC(224) = EE*Sqrt2*Zu14*(Zn15/CW+3D0*Zn25/SW)
      AAABC(225) = EE/CW*Sqrt2*Zn15*Zu44
      AAABC(226) = EE*Sqrt2*(Zd33/CW*Zn13-3D0*Zd33/SW*Zn23
     &      +3D0*Mb/MW/SW*Zd63*Zn33/cb)
      AAABC(227) = EE*Sqrt2*(2D0*Zd63/CW*Zn13+3D0*Mb/MW/SW
     &      *Zd33*dconjg(Zn33)/cb)
      AAABC(228) = EE*Sqrt2*(Zd33/CW*Zn14-3D0*Zd33/SW*Zn24
     &      +3D0*Mb/MW/SW*Zd63*Zn34/cb)
      AAABC(229) = EE*Sqrt2*(2D0*Zd63/CW*Zn14+3D0*Mb/MW/SW
     &      *Zd33*dconjg(Zn34)/cb)
      AAABC(230) = EE*Sqrt2*(Zd33/CW*Zn11-3D0*Zd33/SW*Zn21
     &      +3D0*Mb/MW/SW*Zd63*Zn31/cb)
      AAABC(231) = EE*Sqrt2*(2D0*Zd63/CW*Zn11+3D0*Mb/MW/SW
     &      *Zd33*dconjg(Zn31)/cb)
      AAABC(232) = EE*Sqrt2*(Zd33/CW*Zn12-3D0*Zd33/SW*Zn22
     &      +3D0*Mb/MW/SW*Zd63*Zn32/cb)
      AAABC(233) = EE*Sqrt2*(2D0*Zd63/CW*Zn12+3D0*Mb/MW/SW
     &      *Zd33*dconjg(Zn32)/cb)
      AAABC(234) = EE*Sqrt2*(Zd36/CW*Zn11-3D0*Zd36/SW*Zn21
     &      +3D0*Mb/MW/SW*Zd66*Zn31/cb)
      AAABC(235) = EE*Sqrt2*(2D0*Zd66/CW*Zn11+3D0*Mb/MW/SW
     &      *Zd36*dconjg(Zn31)/cb)
      AAABC(236) = EE*Sqrt2*(Zd36/CW*Zn12-3D0*Zd36/SW*Zn22
     &      +3D0*Mb/MW/SW*Zd66*Zn32/cb)
      AAABC(237) = EE*Sqrt2*(2D0*Zd66/CW*Zn12+3D0*Mb/MW/SW
     &      *Zd36*dconjg(Zn32)/cb)
      AAABC(238) = EE*Sqrt2*(Zd36/CW*Zn14-3D0*Zd36/SW*Zn24
     &      +3D0*Mb/MW/SW*Zd66*Zn34/cb)
      AAABC(239) = EE*Sqrt2*(2D0*Zd66/CW*Zn14+3D0*Mb/MW/SW
     &      *Zd36*dconjg(Zn34)/cb)
      AAABC(240) = EE*Sqrt2*(Zd36/CW*Zn13-3D0*Zd36/SW*Zn23
     &      +3D0*Mb/MW/SW*Zd66*Zn33/cb)
      AAABC(241) = EE*Sqrt2*(2D0*Zd66/CW*Zn13+3D0*Mb/MW/SW
     &      *Zd36*dconjg(Zn33)/cb)
      AAABC(242) = EE*Sqrt2*(Zd33/CW*Zn15-3D0*Zd33/SW*Zn25
     &      +3D0*Mb/MW/SW*Zd63*Zn35/cb)
      AAABC(243) = EE*Sqrt2*(2D0*Zd63/CW*Zn15+3D0*Mb/MW/SW
     &      *Zd33*dconjg(Zn35)/cb)
      AAABC(244) = EE*Sqrt2*(Zd36/CW*Zn15-3D0*Zd36/SW*Zn25
     &      +3D0*Mb/MW/SW*Zd66*Zn35/cb)
      AAABC(245) = EE*Sqrt2*(2D0*Zd66/CW*Zn15+3D0*Mb/MW/SW
     &      *Zd36*dconjg(Zn35)/cb)
      AAABC(246) = EE*Sqrt2*(Zn11/CW*Zu22+3D0*Zn21/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn41*Zu52/sb)
      AAABC(247) = EE*Sqrt2*(4D0*Zn11/CW*Zu52-3D0*Mc/MW/SW
     &      *dconjg(Zn41)*Zu22/sb)
      AAABC(248) = EE*Sqrt2*(Zn12/CW*Zu22+3D0*Zn22/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn42*Zu52/sb)
      AAABC(249) = EE*Sqrt2*(4D0*Zn12/CW*Zu52-3D0*Mc/MW/SW
     &      *dconjg(Zn42)*Zu22/sb)
      AAABC(250) = EE*Sqrt2*(Zn13/CW*Zu22+3D0*Zn23/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn43*Zu52/sb)
      AAABC(251) = EE*Sqrt2*(4D0*Zn13/CW*Zu52-3D0*Mc/MW/SW
     &      *dconjg(Zn43)*Zu22/sb)
      AAABC(252) = EE*Sqrt2*(Zn14/CW*Zu22+3D0*Zn24/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn44*Zu52/sb)
      AAABC(253) = EE*Sqrt2*(4D0*Zn14/CW*Zu52-3D0*Mc/MW/SW
     &      *dconjg(Zn44)*Zu22/sb)
      AAABC(254) = EE*Sqrt2*(Zn11/CW*Zu25+3D0*Zn21/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn41*Zu55/sb)
      AAABC(255) = EE*Sqrt2*(4D0*Zn11/CW*Zu55-3D0*Mc/MW/SW
     &      *dconjg(Zn41)*Zu25/sb)
      AAABC(256) = EE*Sqrt2*(Zn12/CW*Zu25+3D0*Zn22/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn42*Zu55/sb)
      AAABC(257) = EE*Sqrt2*(4D0*Zn12/CW*Zu55-3D0*Mc/MW/SW
     &      *dconjg(Zn42)*Zu25/sb)
      AAABC(258) = EE*Sqrt2*(Zn13/CW*Zu25+3D0*Zn23/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn43*Zu55/sb)
      AAABC(259) = EE*Sqrt2*(4D0*Zn13/CW*Zu55-3D0*Mc/MW/SW
     &      *dconjg(Zn43)*Zu25/sb)
      AAABC(260) = EE*Sqrt2*(Zn14/CW*Zu25+3D0*Zn24/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn44*Zu55/sb)
      AAABC(261) = EE*Sqrt2*(4D0*Zn14/CW*Zu55-3D0*Mc/MW/SW
     &      *dconjg(Zn44)*Zu25/sb)
      AAABC(262) = EE*Sqrt2*(Zn15/CW*Zu22+3D0*Zn25/SW*Zu22
     &      +3D0*Mc/MW/SW*Zn45*Zu52/sb)
      AAABC(263) = EE*Sqrt2*(4D0*Zn15/CW*Zu52-3D0*Mc/MW/SW
     &      *dconjg(Zn45)*Zu22/sb)
      AAABC(264) = EE*Sqrt2*(Zn15/CW*Zu25+3D0*Zn25/SW*Zu25
     &      +3D0*Mc/MW/SW*Zn45*Zu55/sb)
      AAABC(265) = EE*Sqrt2*(4D0*Zn15/CW*Zu55-3D0*Mc/MW/SW
     &      *dconjg(Zn45)*Zu25/sb)
      AAABC(266) = EE*Sqrt2*(Zl11/CW*Zn13+Zl11/SW*Zn23-Me
     &      /MW/SW*Zl41*Zn33/cb)
      AAABC(267) = EE*Sqrt2*(2D0*Zl41/CW*Zn13+Me/MW/SW*Zl11
     &      *dconjg(Zn33)/cb)
      AAABC(268) = EE*Sqrt2*(Zl11/CW*Zn11+Zl11/SW*Zn21-Me
     &      /MW/SW*Zl41*Zn31/cb)
      AAABC(269) = EE*Sqrt2*(2D0*Zl41/CW*Zn11+Me/MW/SW*Zl11
     &      *dconjg(Zn31)/cb)
      AAABC(270) = EE*Sqrt2*(Zl11/CW*Zn14+Zl11/SW*Zn24-Me
     &      /MW/SW*Zl41*Zn34/cb)
      AAABC(271) = EE*Sqrt2*(2D0*Zl41/CW*Zn14+Me/MW/SW*Zl11
     &      *dconjg(Zn34)/cb)
      AAABC(272) = EE*Sqrt2*(Zl11/CW*Zn12+Zl11/SW*Zn22-Me
     &      /MW/SW*Zl41*Zn32/cb)
      AAABC(273) = EE*Sqrt2*(2D0*Zl41/CW*Zn12+Me/MW/SW*Zl11
     &      *dconjg(Zn32)/cb)
      AAABC(274) = EE*Sqrt2*(Zl14/CW*Zn13+Zl14/SW*Zn23-Me
     &      /MW/SW*Zl44*Zn33/cb)
      AAABC(275) = EE*Sqrt2*(2D0*Zl44/CW*Zn13+Me/MW/SW*Zl14
     &      *dconjg(Zn33)/cb)
      AAABC(276) = EE*Sqrt2*(Zl14/CW*Zn11+Zl14/SW*Zn21-Me
     &      /MW/SW*Zl44*Zn31/cb)
      AAABC(277) = EE*Sqrt2*(2D0*Zl44/CW*Zn11+Me/MW/SW*Zl14
     &      *dconjg(Zn31)/cb)
      AAABC(278) = EE*Sqrt2*(Zl14/CW*Zn12+Zl14/SW*Zn22-Me
     &      /MW/SW*Zl44*Zn32/cb)
      AAABC(279) = EE*Sqrt2*(2D0*Zl44/CW*Zn12+Me/MW/SW*Zl14
     &      *dconjg(Zn32)/cb)
      AAABC(280) = EE*Sqrt2*(Zl14/CW*Zn14+Zl14/SW*Zn24-Me
     &      /MW/SW*Zl44*Zn34/cb)
      AAABC(281) = EE*Sqrt2*(2D0*Zl44/CW*Zn14+Me/MW/SW*Zl14
     &      *dconjg(Zn34)/cb)
      AAABC(282) = EE*Sqrt2*(Zl11/CW*Zn15+Zl11/SW*Zn25-Me
     &      /MW/SW*Zl41*Zn35/cb)
      AAABC(283) = EE*Sqrt2*(2D0*Zl41/CW*Zn15+Me/MW/SW*Zl11
     &      *dconjg(Zn35)/cb)
      AAABC(284) = EE*Sqrt2*(Zl14/CW*Zn15+Zl14/SW*Zn25-Me
     &      /MW/SW*Zl44*Zn35/cb)
      AAABC(285) = EE*Sqrt2*(2D0*Zl44/CW*Zn15+Me/MW/SW*Zl14
     &      *dconjg(Zn35)/cb)
      AAABC(286) = EE*Sqrt2*(Zl33/CW*Zn14+Zl33/SW*Zn24-Ml
     &      /MW/SW*Zl63*Zn34/cb)
      AAABC(287) = EE*Sqrt2*(2D0*Zl63/CW*Zn14+Ml/MW/SW*Zl33
     &      *dconjg(Zn34)/cb)
      AAABC(288) = EE*Sqrt2*(Zl33/CW*Zn12+Zl33/SW*Zn22-Ml
     &      /MW/SW*Zl63*Zn32/cb)
      AAABC(289) = EE*Sqrt2*(2D0*Zl63/CW*Zn12+Ml/MW/SW*Zl33
     &      *dconjg(Zn32)/cb)
      AAABC(290) = EE*Sqrt2*(Zl33/CW*Zn11+Zl33/SW*Zn21-Ml
     &      /MW/SW*Zl63*Zn31/cb)
      AAABC(291) = EE*Sqrt2*(2D0*Zl63/CW*Zn11+Ml/MW/SW*Zl33
     &      *dconjg(Zn31)/cb)
      AAABC(292) = EE*Sqrt2*(Zl33/CW*Zn13+Zl33/SW*Zn23-Ml
     &      /MW/SW*Zl63*Zn33/cb)
      AAABC(293) = EE*Sqrt2*(2D0*Zl63/CW*Zn13+Ml/MW/SW*Zl33
     &      *dconjg(Zn33)/cb)
      AAABC(294) = EE*Sqrt2*(Zl36/CW*Zn13+Zl36/SW*Zn23-Ml
     &      /MW/SW*Zl66*Zn33/cb)
      AAABC(295) = EE*Sqrt2*(2D0*Zl66/CW*Zn13+Ml/MW/SW*Zl36
     &      *dconjg(Zn33)/cb)
      AAABC(296) = EE*Sqrt2*(Zl36/CW*Zn11+Zl36/SW*Zn21-Ml
     &      /MW/SW*Zl66*Zn31/cb)
      AAABC(297) = EE*Sqrt2*(2D0*Zl66/CW*Zn11+Ml/MW/SW*Zl36
     &      *dconjg(Zn31)/cb)
      AAABC(298) = EE*Sqrt2*(Zl36/CW*Zn14+Zl36/SW*Zn24-Ml
     &      /MW/SW*Zl66*Zn34/cb)
      AAABC(299) = EE*Sqrt2*(2D0*Zl66/CW*Zn14+Ml/MW/SW*Zl36
     &      *dconjg(Zn34)/cb)
      AAABC(300) = EE*Sqrt2*(Zl36/CW*Zn12+Zl36/SW*Zn22-Ml
     &      /MW/SW*Zl66*Zn32/cb)
      AAABC(301) = EE*Sqrt2*(2D0*Zl66/CW*Zn12+Ml/MW/SW*Zl36
     &      *dconjg(Zn32)/cb)
      AAABC(302) = EE*Sqrt2*(Zl33/CW*Zn15+Zl33/SW*Zn25-Ml
     &      /MW/SW*Zl63*Zn35/cb)
      AAABC(303) = EE*Sqrt2*(2D0*Zl63/CW*Zn15+Ml/MW/SW*Zl33
     &      *dconjg(Zn35)/cb)
      AAABC(304) = EE*Sqrt2*(Zl36/CW*Zn15+Zl36/SW*Zn25-Ml
     &      /MW/SW*Zl66*Zn35/cb)
      AAABC(305) = EE*Sqrt2*(2D0*Zl66/CW*Zn15+Ml/MW/SW*Zl36
     &      *dconjg(Zn35)/cb)
      AAABC(306) = EE*Sqrt2*(Zl22/CW*Zn13+Zl22/SW*Zn23-Mm
     &      /MW/SW*Zl52*Zn33/cb)
      AAABC(307) = EE*Sqrt2*(2D0*Zl52/CW*Zn13+Mm/MW/SW*Zl22
     &      *dconjg(Zn33)/cb)
      AAABC(308) = EE*Sqrt2*(Zl22/CW*Zn12+Zl22/SW*Zn22-Mm
     &      /MW/SW*Zl52*Zn32/cb)
      AAABC(309) = EE*Sqrt2*(2D0*Zl52/CW*Zn12+Mm/MW/SW*Zl22
     &      *dconjg(Zn32)/cb)
      AAABC(310) = EE*Sqrt2*(Zl22/CW*Zn11+Zl22/SW*Zn21-Mm
     &      /MW/SW*Zl52*Zn31/cb)
      AAABC(311) = EE*Sqrt2*(2D0*Zl52/CW*Zn11+Mm/MW/SW*Zl22
     &      *dconjg(Zn31)/cb)
      AAABC(312) = EE*Sqrt2*(Zl22/CW*Zn14+Zl22/SW*Zn24-Mm
     &      /MW/SW*Zl52*Zn34/cb)
      AAABC(313) = EE*Sqrt2*(2D0*Zl52/CW*Zn14+Mm/MW/SW*Zl22
     &      *dconjg(Zn34)/cb)
      AAABC(314) = EE*Sqrt2*(Zl25/CW*Zn12+Zl25/SW*Zn22-Mm
     &      /MW/SW*Zl55*Zn32/cb)
      AAABC(315) = EE*Sqrt2*(2D0*Zl55/CW*Zn12+Mm/MW/SW*Zl25
     &      *dconjg(Zn32)/cb)
      AAABC(316) = EE*Sqrt2*(Zl25/CW*Zn13+Zl25/SW*Zn23-Mm
     &      /MW/SW*Zl55*Zn33/cb)
      AAABC(317) = EE*Sqrt2*(2D0*Zl55/CW*Zn13+Mm/MW/SW*Zl25
     &      *dconjg(Zn33)/cb)
      AAABC(318) = EE*Sqrt2*(Zl25/CW*Zn14+Zl25/SW*Zn24-Mm
     &      /MW/SW*Zl55*Zn34/cb)
      AAABC(319) = EE*Sqrt2*(2D0*Zl55/CW*Zn14+Mm/MW/SW*Zl25
     &      *dconjg(Zn34)/cb)
      AAABC(320) = EE*Sqrt2*(Zl25/CW*Zn11+Zl25/SW*Zn21-Mm
     &      /MW/SW*Zl55*Zn31/cb)
      AAABC(321) = EE*Sqrt2*(2D0*Zl55/CW*Zn11+Mm/MW/SW*Zl25
     &      *dconjg(Zn31)/cb)
      AAABC(322) = EE*Sqrt2*(Zl22/CW*Zn15+Zl22/SW*Zn25-Mm
     &      /MW/SW*Zl52*Zn35/cb)
      AAABC(323) = EE*Sqrt2*(2D0*Zl52/CW*Zn15+Mm/MW/SW*Zl22
     &      *dconjg(Zn35)/cb)
      AAABC(324) = EE*Sqrt2*(Zl25/CW*Zn15+Zl25/SW*Zn25-Mm
     &      /MW/SW*Zl55*Zn35/cb)
      AAABC(325) = EE*Sqrt2*(2D0*Zl55/CW*Zn15+Mm/MW/SW*Zl25
     &      *dconjg(Zn35)/cb)
      AAABC(326) = EE*Sqrt2*(2D0*Zd52/CW*Zn13+3D0*Ms/MW/SW
     &      *Zd22*dconjg(Zn33)/cb)
      AAABC(327) = EE*Sqrt2*(Zd22/CW*Zn13-3D0*Zd22/SW*Zn23
     &      +3D0*Ms/MW/SW*Zd52*Zn33/cb)
      AAABC(328) = EE*Sqrt2*(2D0*Zd52/CW*Zn14+3D0*Ms/MW/SW
     &      *Zd22*dconjg(Zn34)/cb)
      AAABC(329) = EE*Sqrt2*(Zd22/CW*Zn14-3D0*Zd22/SW*Zn24
     &      +3D0*Ms/MW/SW*Zd52*Zn34/cb)
      AAABC(330) = EE*Sqrt2*(2D0*Zd52/CW*Zn12+3D0*Ms/MW/SW
     &      *Zd22*dconjg(Zn32)/cb)
      AAABC(331) = EE*Sqrt2*(Zd22/CW*Zn12-3D0*Zd22/SW*Zn22
     &      +3D0*Ms/MW/SW*Zd52*Zn32/cb)
      AAABC(332) = EE*Sqrt2*(2D0*Zd52/CW*Zn11+3D0*Ms/MW/SW
     &      *Zd22*dconjg(Zn31)/cb)
      AAABC(333) = EE*Sqrt2*(Zd22/CW*Zn11-3D0*Zd22/SW*Zn21
     &      +3D0*Ms/MW/SW*Zd52*Zn31/cb)
      AAABC(334) = EE*Sqrt2*(2D0*Zd55/CW*Zn14+3D0*Ms/MW/SW
     &      *Zd25*dconjg(Zn34)/cb)
      AAABC(335) = EE*Sqrt2*(Zd25/CW*Zn14-3D0*Zd25/SW*Zn24
     &      +3D0*Ms/MW/SW*Zd55*Zn34/cb)
      AAABC(336) = EE*Sqrt2*(2D0*Zd55/CW*Zn13+3D0*Ms/MW/SW
     &      *Zd25*dconjg(Zn33)/cb)
      AAABC(337) = EE*Sqrt2*(Zd25/CW*Zn13-3D0*Zd25/SW*Zn23
     &      +3D0*Ms/MW/SW*Zd55*Zn33/cb)
      AAABC(338) = EE*Sqrt2*(2D0*Zd55/CW*Zn12+3D0*Ms/MW/SW
     &      *Zd25*dconjg(Zn32)/cb)
      AAABC(339) = EE*Sqrt2*(Zd25/CW*Zn12-3D0*Zd25/SW*Zn22
     &      +3D0*Ms/MW/SW*Zd55*Zn32/cb)
      AAABC(340) = EE*Sqrt2*(2D0*Zd55/CW*Zn11+3D0*Ms/MW/SW
     &      *Zd25*dconjg(Zn31)/cb)
      AAABC(341) = EE*Sqrt2*(Zd25/CW*Zn11-3D0*Zd25/SW*Zn21
     &      +3D0*Ms/MW/SW*Zd55*Zn31/cb)
      AAABC(342) = EE*Sqrt2*(2D0*Zd52/CW*Zn15+3D0*Ms/MW/SW
     &      *Zd22*dconjg(Zn35)/cb)
      AAABC(343) = EE*Sqrt2*(Zd22/CW*Zn15-3D0*Zd22/SW*Zn25
     &      +3D0*Ms/MW/SW*Zd52*Zn35/cb)
      AAABC(344) = EE*Sqrt2*(2D0*Zd55/CW*Zn15+3D0*Ms/MW/SW
     &      *Zd25*dconjg(Zn35)/cb)
      AAABC(345) = EE*Sqrt2*(Zd25/CW*Zn15-3D0*Zd25/SW*Zn25
     &      +3D0*Ms/MW/SW*Zd55*Zn35/cb)
      AAABC(346) = EE*Sqrt2*(Zn14/CW*Zu33+3D0*Zn24/SW*Zu33
     &      +3D0*Mt/MW/SW*Zn44*Zu63/sb)
      AAABC(347) = EE*Sqrt2*(4D0*Zn14/CW*Zu63-3D0*Mt/MW/SW
     &      *dconjg(Zn44)*Zu33/sb)
      AAABC(348) = EE*Sqrt2*(Zn13/CW*Zu33+3D0*Zn23/SW*Zu33
     &      +3D0*Mt/MW/SW*Zn43*Zu63/sb)
      AAABC(349) = EE*Sqrt2*(4D0*Zn13/CW*Zu63-3D0*Mt/MW/SW
     &      *dconjg(Zn43)*Zu33/sb)
      AAABC(350) = EE*Sqrt2*(Zn11/CW*Zu33+3D0*Zn21/SW*Zu33
     &      +3D0*Mt/MW/SW*Zn41*Zu63/sb)
      AAABC(351) = EE*Sqrt2*(4D0*Zn11/CW*Zu63-3D0*Mt/MW/SW
     &      *dconjg(Zn41)*Zu33/sb)
      AAABC(352) = EE*Sqrt2*(Zn12/CW*Zu33+3D0*Zn22/SW*Zu33
     &      +3D0*Mt/MW/SW*Zn42*Zu63/sb)
      AAABC(353) = EE*Sqrt2*(4D0*Zn12/CW*Zu63-3D0*Mt/MW/SW
     &      *dconjg(Zn42)*Zu33/sb)
      AAABC(354) = EE*Sqrt2*(Zn13/CW*Zu36+3D0*Zn23/SW*Zu36
     &      +3D0*Mt/MW/SW*Zn43*Zu66/sb)
      AAABC(355) = EE*Sqrt2*(4D0*Zn13/CW*Zu66-3D0*Mt/MW/SW
     &      *dconjg(Zn43)*Zu36/sb)
      AAABC(356) = EE*Sqrt2*(Zn14/CW*Zu36+3D0*Zn24/SW*Zu36
     &      +3D0*Mt/MW/SW*Zn44*Zu66/sb)
      AAABC(357) = EE*Sqrt2*(4D0*Zn14/CW*Zu66-3D0*Mt/MW/SW
     &      *dconjg(Zn44)*Zu36/sb)
      AAABC(358) = EE*Sqrt2*(Zn12/CW*Zu36+3D0*Zn22/SW*Zu36
     &      +3D0*Mt/MW/SW*Zn42*Zu66/sb)
      AAABC(359) = EE*Sqrt2*(4D0*Zn12/CW*Zu66-3D0*Mt/MW/SW
     &      *dconjg(Zn42)*Zu36/sb)
      AAABC(360) = EE*Sqrt2*(Zn11/CW*Zu36+3D0*Zn21/SW*Zu36
     &      +3D0*Mt/MW/SW*Zn41*Zu66/sb)
      AAABC(361) = EE*Sqrt2*(4D0*Zn11/CW*Zu66-3D0*Mt/MW/SW
     &      *dconjg(Zn41)*Zu36/sb)
      AAABC(362) = EE*Sqrt2*(Zn15/CW*Zu33+3D0*Zn25/SW*Zu33
     &      +3D0*Mt/MW/SW*Zn45*Zu63/sb)
      AAABC(363) = EE*Sqrt2*(4D0*Zn15/CW*Zu63-3D0*Mt/MW/SW
     &      *dconjg(Zn45)*Zu33/sb)
      AAABC(364) = EE*Sqrt2*(Zn15/CW*Zu36+3D0*Zn25/SW*Zu36
     &      +3D0*Mt/MW/SW*Zn45*Zu66/sb)
      AAABC(365) = EE*Sqrt2*(4D0*Zn15/CW*Zu66-3D0*Mt/MW/SW
     &      *dconjg(Zn45)*Zu36/sb)
      AAABC(366) = 2D0*Zm21*dconjg(Zn53)*cb*hL+2D0*EE/SW
     &      *Zm11*dconjg(Zn33)*sb-EE/CW*Sqrt2*Zm21*Zn13*sb
     &      -EE/SW*Sqrt2*Zm21*Zn23*sb
      AAABC(367) = 2D0*Zn53*Zp21*hL*sb+2D0*EE/SW*Zn43*Zp11
     &      *cb+EE/CW*Sqrt2*Zn13*Zp21*cb+EE/SW*Sqrt2*Zn23
     &      *Zp21*cb
      AAABC(368) = 2D0*Zm21*dconjg(Zn52)*cb*hL+2D0*EE/SW
     &      *Zm11*dconjg(Zn32)*sb-EE/CW*Sqrt2*Zm21*Zn12*sb
     &      -EE/SW*Sqrt2*Zm21*Zn22*sb
      AAABC(369) = 2D0*Zn52*Zp21*hL*sb+2D0*EE/SW*Zn42*Zp11
     &      *cb+EE/CW*Sqrt2*Zn12*Zp21*cb+EE/SW*Sqrt2*Zn22
     &      *Zp21*cb
      AAABC(370) = 2D0*Zm21*dconjg(Zn51)*cb*hL+2D0*EE/SW
     &      *Zm11*dconjg(Zn31)*sb-EE/CW*Sqrt2*Zm21*Zn11*sb
     &      -EE/SW*Sqrt2*Zm21*Zn21*sb
      AAABC(371) = 2D0*Zn51*Zp21*hL*sb+2D0*EE/SW*Zn41*Zp11
     &      *cb+EE/CW*Sqrt2*Zn11*Zp21*cb+EE/SW*Sqrt2*Zn21
     &      *Zp21*cb
      AAABC(372) = 2D0*Zm21*dconjg(Zn54)*cb*hL+2D0*EE/SW
     &      *Zm11*dconjg(Zn34)*sb-EE/CW*Sqrt2*Zm21*Zn14*sb
     &      -EE/SW*Sqrt2*Zm21*Zn24*sb
      AAABC(373) = 2D0*Zn54*Zp21*hL*sb+2D0*EE/SW*Zn44*Zp11
     &      *cb+EE/CW*Sqrt2*Zn14*Zp21*cb+EE/SW*Sqrt2*Zn24
     &      *Zp21*cb
      AAABC(374) = 2D0*Zm21*dconjg(Zn52)*hL*sb-2D0*EE/SW
     &      *Zm11*dconjg(Zn32)*cb+EE/CW*Sqrt2*Zm21*Zn12*cb
     &      +EE/SW*Sqrt2*Zm21*Zn22*cb
      AAABC(375) = 2D0*Zn52*Zp21*cb*hL-2D0*EE/SW*Zn42*Zp11
     &      *sb-EE/CW*Sqrt2*Zn12*Zp21*sb-EE/SW*Sqrt2*Zn22
     &      *Zp21*sb
      AAABC(376) = 2D0*Zm21*dconjg(Zn51)*hL*sb-2D0*EE/SW
     &      *Zm11*dconjg(Zn31)*cb+EE/CW*Sqrt2*Zm21*Zn11*cb
     &      +EE/SW*Sqrt2*Zm21*Zn21*cb
      AAABC(377) = 2D0*Zn51*Zp21*cb*hL-2D0*EE/SW*Zn41*Zp11
     &      *sb-EE/CW*Sqrt2*Zn11*Zp21*sb-EE/SW*Sqrt2*Zn21
     &      *Zp21*sb
      AAABC(378) = 2D0*Zm21*dconjg(Zn54)*hL*sb-2D0*EE/SW
     &      *Zm11*dconjg(Zn34)*cb+EE/CW*Sqrt2*Zm21*Zn14*cb
     &      +EE/SW*Sqrt2*Zm21*Zn24*cb
      AAABC(379) = 2D0*Zn54*Zp21*cb*hL-2D0*EE/SW*Zn44*Zp11
     &      *sb-EE/CW*Sqrt2*Zn14*Zp21*sb-EE/SW*Sqrt2*Zn24
     &      *Zp21*sb
      AAABC(380) = 2D0*Zm21*dconjg(Zn53)*hL*sb-2D0*EE/SW
     &      *Zm11*dconjg(Zn33)*cb+EE/CW*Sqrt2*Zm21*Zn13*cb
     &      +EE/SW*Sqrt2*Zm21*Zn23*cb
      AAABC(381) = 2D0*Zn53*Zp21*cb*hL-2D0*EE/SW*Zn43*Zp11
     &      *sb-EE/CW*Sqrt2*Zn13*Zp21*sb-EE/SW*Sqrt2*Zn23
     &      *Zp21*sb
      AAABC(382) = 2D0*Zm21*dconjg(Zn55)*cb*hL+2D0*EE/SW
     &      *Zm11*dconjg(Zn35)*sb-EE/CW*Sqrt2*Zm21*Zn15*sb
     &      -EE/SW*Sqrt2*Zm21*Zn25*sb
      AAABC(383) = 2D0*Zn55*Zp21*hL*sb+2D0*EE/SW*Zn45*Zp11
     &      *cb+EE/CW*Sqrt2*Zn15*Zp21*cb+EE/SW*Sqrt2*Zn25
     &      *Zp21*cb
      AAABC(384) = 2D0*Zm21*dconjg(Zn55)*hL*sb-2D0*EE/SW
     &      *Zm11*dconjg(Zn35)*cb+EE/CW*Sqrt2*Zm21*Zn15*cb
     &      +EE/SW*Sqrt2*Zm21*Zn25*cb
      AAABC(385) = 2D0*Zn55*Zp21*cb*hL-2D0*EE/SW*Zn45*Zp11
     &      *sb-EE/CW*Sqrt2*Zn15*Zp21*sb-EE/SW*Sqrt2*Zn25
     &      *Zp21*sb
      AAABC(386) = 2D0*Zm22*dconjg(Zn53)*cb*hL+2D0*EE/SW
     &      *Zm12*dconjg(Zn33)*sb-EE/CW*Sqrt2*Zm22*Zn13*sb
     &      -EE/SW*Sqrt2*Zm22*Zn23*sb
      AAABC(387) = 2D0*Zn53*Zp22*hL*sb+2D0*EE/SW*Zn43*Zp12
     &      *cb+EE/CW*Sqrt2*Zn13*Zp22*cb+EE/SW*Sqrt2*Zn23
     &      *Zp22*cb
      AAABC(388) = 2D0*Zm22*dconjg(Zn51)*cb*hL+2D0*EE/SW
     &      *Zm12*dconjg(Zn31)*sb-EE/CW*Sqrt2*Zm22*Zn11*sb
     &      -EE/SW*Sqrt2*Zm22*Zn21*sb
      AAABC(389) = 2D0*Zn51*Zp22*hL*sb+2D0*EE/SW*Zn41*Zp12
     &      *cb+EE/CW*Sqrt2*Zn11*Zp22*cb+EE/SW*Sqrt2*Zn21
     &      *Zp22*cb
      AAABC(390) = 2D0*Zm22*dconjg(Zn52)*cb*hL+2D0*EE/SW
     &      *Zm12*dconjg(Zn32)*sb-EE/CW*Sqrt2*Zm22*Zn12*sb
     &      -EE/SW*Sqrt2*Zm22*Zn22*sb
      AAABC(391) = 2D0*Zn52*Zp22*hL*sb+2D0*EE/SW*Zn42*Zp12
     &      *cb+EE/CW*Sqrt2*Zn12*Zp22*cb+EE/SW*Sqrt2*Zn22
     &      *Zp22*cb
      AAABC(392) = 2D0*Zm22*dconjg(Zn54)*cb*hL+2D0*EE/SW
     &      *Zm12*dconjg(Zn34)*sb-EE/CW*Sqrt2*Zm22*Zn14*sb
     &      -EE/SW*Sqrt2*Zm22*Zn24*sb
      AAABC(393) = 2D0*Zn54*Zp22*hL*sb+2D0*EE/SW*Zn44*Zp12
     &      *cb+EE/CW*Sqrt2*Zn14*Zp22*cb+EE/SW*Sqrt2*Zn24
     &      *Zp22*cb
      AAABC(394) = 2D0*Zm22*dconjg(Zn53)*hL*sb-2D0*EE/SW
     &      *Zm12*dconjg(Zn33)*cb+EE/CW*Sqrt2*Zm22*Zn13*cb
     &      +EE/SW*Sqrt2*Zm22*Zn23*cb
      AAABC(395) = 2D0*Zn53*Zp22*cb*hL-2D0*EE/SW*Zn43*Zp12
     &      *sb-EE/CW*Sqrt2*Zn13*Zp22*sb-EE/SW*Sqrt2*Zn23
     &      *Zp22*sb
      AAABC(396) = 2D0*Zm22*dconjg(Zn54)*hL*sb-2D0*EE/SW
     &      *Zm12*dconjg(Zn34)*cb+EE/CW*Sqrt2*Zm22*Zn14*cb
     &      +EE/SW*Sqrt2*Zm22*Zn24*cb
      AAABC(397) = 2D0*Zn54*Zp22*cb*hL-2D0*EE/SW*Zn44*Zp12
     &      *sb-EE/CW*Sqrt2*Zn14*Zp22*sb-EE/SW*Sqrt2*Zn24
     &      *Zp22*sb
      AAABC(398) = 2D0*Zm22*dconjg(Zn51)*hL*sb-2D0*EE/SW
     &      *Zm12*dconjg(Zn31)*cb+EE/CW*Sqrt2*Zm22*Zn11*cb
     &      +EE/SW*Sqrt2*Zm22*Zn21*cb
      AAABC(399) = 2D0*Zn51*Zp22*cb*hL-2D0*EE/SW*Zn41*Zp12
     &      *sb-EE/CW*Sqrt2*Zn11*Zp22*sb-EE/SW*Sqrt2*Zn21
     &      *Zp22*sb
      AAABC(400) = 2D0*Zm22*dconjg(Zn52)*hL*sb-2D0*EE/SW
     &      *Zm12*dconjg(Zn32)*cb+EE/CW*Sqrt2*Zm22*Zn12*cb
     &      +EE/SW*Sqrt2*Zm22*Zn22*cb
      AAABC(401) = 2D0*Zn52*Zp22*cb*hL-2D0*EE/SW*Zn42*Zp12
     &      *sb-EE/CW*Sqrt2*Zn12*Zp22*sb-EE/SW*Sqrt2*Zn22
     &      *Zp22*sb
      AAABC(402) = 2D0*Zm22*dconjg(Zn55)*cb*hL+2D0*EE/SW
     &      *Zm12*dconjg(Zn35)*sb-EE/CW*Sqrt2*Zm22*Zn15*sb
     &      -EE/SW*Sqrt2*Zm22*Zn25*sb
      AAABC(403) = 2D0*Zn55*Zp22*hL*sb+2D0*EE/SW*Zn45*Zp12
     &      *cb+EE/CW*Sqrt2*Zn15*Zp22*cb+EE/SW*Sqrt2*Zn25
     &      *Zp22*cb
      AAABC(404) = 2D0*Zm22*dconjg(Zn55)*hL*sb-2D0*EE/SW
     &      *Zm12*dconjg(Zn35)*cb+EE/CW*Sqrt2*Zm22*Zn15*cb
     &      +EE/SW*Sqrt2*Zm22*Zn25*cb
      AAABC(405) = 2D0*Zn55*Zp22*cb*hL-2D0*EE/SW*Zn45*Zp12
     &      *sb-EE/CW*Sqrt2*Zn15*Zp22*sb-EE/SW*Sqrt2*Zn25
     &      *Zp22*sb
      AAABC(406) = Sqrt2*Za33*Zn51**2*hK+EE/CW*Za13*Zn11
     &      *Zn31-EE/CW*Za23*Zn11*Zn41-EE/SW*Za13*Zn21*Zn31
     &      +EE/SW*Za23*Zn21*Zn41-Sqrt2*Za13*Zn41*Zn51*hL
     &      -Sqrt2*Za23*Zn31*Zn51*hL-Sqrt2*Za33*Zn31*Zn41*hL
      AAABC(407) = Sqrt2*Za33*dconjg(Zn51)**2*hK+EE/CW*Za13
     &      *Zn11*dconjg(Zn31)-EE/CW*Za23*Zn11*dconjg(Zn41
     &      )-EE/SW*Za13*Zn21*dconjg(Zn31)+EE/SW*Za23*Zn21
     &      *dconjg(Zn41)-Sqrt2*Za13*dconjg(Zn41)*dconjg(Zn51
     &      )*hL-Sqrt2*Za23*dconjg(Zn31)*dconjg(Zn51)*hL
     &      -Sqrt2*Za33*dconjg(Zn31)*dconjg(Zn41)*hL
      AAABC(408) = EE/CW*Za13*Zn11*Zn34+EE/CW*Za13*Zn14*Zn31
     &      -EE/CW*Za23*Zn11*Zn44-EE/CW*Za23*Zn14*Zn41-EE
     &      /SW*Za13*Zn21*Zn34-EE/SW*Za13*Zn24*Zn31+EE/SW
     &      *Za23*Zn21*Zn44+EE/SW*Za23*Zn24*Zn41-Sqrt2*Za13
     &      *Zn41*Zn54*hL-Sqrt2*Za13*Zn44*Zn51*hL-Sqrt2*Za23
     &      *Zn31*Zn54*hL-Sqrt2*Za23*Zn34*Zn51*hL-Sqrt2*Za33
     &      *Zn31*Zn44*hL-Sqrt2*Za33*Zn34*Zn41*hL+2D0*Sqrt2
     &      *Za33*Zn51*Zn54*hK
      AAABC(409) = EE/CW*Za13*Zn11*dconjg(Zn34)+EE/CW*Za13
     &      *Zn14*dconjg(Zn31)-EE/CW*Za23*Zn11*dconjg(Zn44
     &      )-EE/CW*Za23*Zn14*dconjg(Zn41)-EE/SW*Za13*Zn21
     &      *dconjg(Zn34)-EE/SW*Za13*Zn24*dconjg(Zn31)+EE
     &      /SW*Za23*Zn21*dconjg(Zn44)+EE/SW*Za23*Zn24*dconjg(Zn41
     &      )-Sqrt2*Za13*dconjg(Zn41)*dconjg(Zn54)*hL-Sqrt2
     &      *Za13*dconjg(Zn44)*dconjg(Zn51)*hL-Sqrt2*Za23
     &      *dconjg(Zn31)*dconjg(Zn54)*hL-Sqrt2*Za23*dconjg(Zn34
     &      )*dconjg(Zn51)*hL-Sqrt2*Za33*dconjg(Zn31)*dconjg(Zn44
     &      )*hL-Sqrt2*Za33*dconjg(Zn34)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn51)*dconjg(Zn54)*hK
      AAABC(410) = EE/CW*Za13*Zn11*Zn32+EE/CW*Za13*Zn12*Zn31
     &      -EE/CW*Za23*Zn11*Zn42-EE/CW*Za23*Zn12*Zn41-EE
     &      /SW*Za13*Zn21*Zn32-EE/SW*Za13*Zn22*Zn31+EE/SW
     &      *Za23*Zn21*Zn42+EE/SW*Za23*Zn22*Zn41-Sqrt2*Za13
     &      *Zn41*Zn52*hL-Sqrt2*Za13*Zn42*Zn51*hL-Sqrt2*Za23
     &      *Zn31*Zn52*hL-Sqrt2*Za23*Zn32*Zn51*hL-Sqrt2*Za33
     &      *Zn31*Zn42*hL-Sqrt2*Za33*Zn32*Zn41*hL+2D0*Sqrt2
     &      *Za33*Zn51*Zn52*hK
      AAABC(411) = EE/CW*Za13*Zn11*dconjg(Zn32)+EE/CW*Za13
     &      *Zn12*dconjg(Zn31)-EE/CW*Za23*Zn11*dconjg(Zn42
     &      )-EE/CW*Za23*Zn12*dconjg(Zn41)-EE/SW*Za13*Zn21
     &      *dconjg(Zn32)-EE/SW*Za13*Zn22*dconjg(Zn31)+EE
     &      /SW*Za23*Zn21*dconjg(Zn42)+EE/SW*Za23*Zn22*dconjg(Zn41
     &      )-Sqrt2*Za13*dconjg(Zn41)*dconjg(Zn52)*hL-Sqrt2
     &      *Za13*dconjg(Zn42)*dconjg(Zn51)*hL-Sqrt2*Za23
     &      *dconjg(Zn31)*dconjg(Zn52)*hL-Sqrt2*Za23*dconjg(Zn32
     &      )*dconjg(Zn51)*hL-Sqrt2*Za33*dconjg(Zn31)*dconjg(Zn42
     &      )*hL-Sqrt2*Za33*dconjg(Zn32)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn51)*dconjg(Zn52)*hK
      AAABC(412) = EE/CW*Za13*Zn11*Zn33+EE/CW*Za13*Zn13*Zn31
     &      -EE/CW*Za23*Zn11*Zn43-EE/CW*Za23*Zn13*Zn41-EE
     &      /SW*Za13*Zn21*Zn33-EE/SW*Za13*Zn23*Zn31+EE/SW
     &      *Za23*Zn21*Zn43+EE/SW*Za23*Zn23*Zn41-Sqrt2*Za13
     &      *Zn41*Zn53*hL-Sqrt2*Za13*Zn43*Zn51*hL-Sqrt2*Za23
     &      *Zn31*Zn53*hL-Sqrt2*Za23*Zn33*Zn51*hL-Sqrt2*Za33
     &      *Zn31*Zn43*hL-Sqrt2*Za33*Zn33*Zn41*hL+2D0*Sqrt2
     &      *Za33*Zn51*Zn53*hK
      AAABC(413) = EE/CW*Za13*Zn11*dconjg(Zn33)+EE/CW*Za13
     &      *Zn13*dconjg(Zn31)-EE/CW*Za23*Zn11*dconjg(Zn43
     &      )-EE/CW*Za23*Zn13*dconjg(Zn41)-EE/SW*Za13*Zn21
     &      *dconjg(Zn33)-EE/SW*Za13*Zn23*dconjg(Zn31)+EE
     &      /SW*Za23*Zn21*dconjg(Zn43)+EE/SW*Za23*Zn23*dconjg(Zn41
     &      )-Sqrt2*Za13*dconjg(Zn41)*dconjg(Zn53)*hL-Sqrt2
     &      *Za13*dconjg(Zn43)*dconjg(Zn51)*hL-Sqrt2*Za23
     &      *dconjg(Zn31)*dconjg(Zn53)*hL-Sqrt2*Za23*dconjg(Zn33
     &      )*dconjg(Zn51)*hL-Sqrt2*Za33*dconjg(Zn31)*dconjg(Zn43
     &      )*hL-Sqrt2*Za33*dconjg(Zn33)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn51)*dconjg(Zn53)*hK
      AAABC(414) = EE/CW*Za12*Zn11*Zn32+EE/CW*Za12*Zn12*Zn31
     &      -EE/CW*Za22*Zn11*Zn42-EE/CW*Za22*Zn12*Zn41-EE
     &      /SW*Za12*Zn21*Zn32-EE/SW*Za12*Zn22*Zn31+EE/SW
     &      *Za22*Zn21*Zn42+EE/SW*Za22*Zn22*Zn41-Sqrt2*Za12
     &      *Zn41*Zn52*hL-Sqrt2*Za12*Zn42*Zn51*hL-Sqrt2*Za22
     &      *Zn31*Zn52*hL-Sqrt2*Za22*Zn32*Zn51*hL-Sqrt2*Za32
     &      *Zn31*Zn42*hL-Sqrt2*Za32*Zn32*Zn41*hL+2D0*Sqrt2
     &      *Za32*Zn51*Zn52*hK
      AAABC(415) = EE/CW*Za12*Zn11*dconjg(Zn32)+EE/CW*Za12
     &      *Zn12*dconjg(Zn31)-EE/CW*Za22*Zn11*dconjg(Zn42
     &      )-EE/CW*Za22*Zn12*dconjg(Zn41)-EE/SW*Za12*Zn21
     &      *dconjg(Zn32)-EE/SW*Za12*Zn22*dconjg(Zn31)+EE
     &      /SW*Za22*Zn21*dconjg(Zn42)+EE/SW*Za22*Zn22*dconjg(Zn41
     &      )-Sqrt2*Za12*dconjg(Zn41)*dconjg(Zn52)*hL-Sqrt2
     &      *Za12*dconjg(Zn42)*dconjg(Zn51)*hL-Sqrt2*Za22
     &      *dconjg(Zn31)*dconjg(Zn52)*hL-Sqrt2*Za22*dconjg(Zn32
     &      )*dconjg(Zn51)*hL-Sqrt2*Za32*dconjg(Zn31)*dconjg(Zn42
     &      )*hL-Sqrt2*Za32*dconjg(Zn32)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn51)*dconjg(Zn52)*hK
      AAABC(416) = EE/CW*Za11*Zn11*Zn33+EE/CW*Za11*Zn13*Zn31
     &      -EE/CW*Za21*Zn11*Zn43-EE/CW*Za21*Zn13*Zn41-EE
     &      /SW*Za11*Zn21*Zn33-EE/SW*Za11*Zn23*Zn31+EE/SW
     &      *Za21*Zn21*Zn43+EE/SW*Za21*Zn23*Zn41-Sqrt2*Za11
     &      *Zn41*Zn53*hL-Sqrt2*Za11*Zn43*Zn51*hL-Sqrt2*Za21
     &      *Zn31*Zn53*hL-Sqrt2*Za21*Zn33*Zn51*hL-Sqrt2*Za31
     &      *Zn31*Zn43*hL-Sqrt2*Za31*Zn33*Zn41*hL+2D0*Sqrt2
     &      *Za31*Zn51*Zn53*hK
      AAABC(417) = EE/CW*Za11*Zn11*dconjg(Zn33)+EE/CW*Za11
     &      *Zn13*dconjg(Zn31)-EE/CW*Za21*Zn11*dconjg(Zn43
     &      )-EE/CW*Za21*Zn13*dconjg(Zn41)-EE/SW*Za11*Zn21
     &      *dconjg(Zn33)-EE/SW*Za11*Zn23*dconjg(Zn31)+EE
     &      /SW*Za21*Zn21*dconjg(Zn43)+EE/SW*Za21*Zn23*dconjg(Zn41
     &      )-Sqrt2*Za11*dconjg(Zn41)*dconjg(Zn53)*hL-Sqrt2
     &      *Za11*dconjg(Zn43)*dconjg(Zn51)*hL-Sqrt2*Za21
     &      *dconjg(Zn31)*dconjg(Zn53)*hL-Sqrt2*Za21*dconjg(Zn33
     &      )*dconjg(Zn51)*hL-Sqrt2*Za31*dconjg(Zn31)*dconjg(Zn43
     &      )*hL-Sqrt2*Za31*dconjg(Zn33)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn51)*dconjg(Zn53)*hK
      AAABC(418) = Sqrt2*Za31*Zn51**2*hK+EE/CW*Za11*Zn11
     &      *Zn31-EE/CW*Za21*Zn11*Zn41-EE/SW*Za11*Zn21*Zn31
     &      +EE/SW*Za21*Zn21*Zn41-Sqrt2*Za11*Zn41*Zn51*hL
     &      -Sqrt2*Za21*Zn31*Zn51*hL-Sqrt2*Za31*Zn31*Zn41*hL
      AAABC(419) = Sqrt2*Za31*dconjg(Zn51)**2*hK+EE/CW*Za11
     &      *Zn11*dconjg(Zn31)-EE/CW*Za21*Zn11*dconjg(Zn41
     &      )-EE/SW*Za11*Zn21*dconjg(Zn31)+EE/SW*Za21*Zn21
     &      *dconjg(Zn41)-Sqrt2*Za11*dconjg(Zn41)*dconjg(Zn51
     &      )*hL-Sqrt2*Za21*dconjg(Zn31)*dconjg(Zn51)*hL
     &      -Sqrt2*Za31*dconjg(Zn31)*dconjg(Zn41)*hL
      AAABC(420) = EE/CW*Za11*Zn11*Zn34+EE/CW*Za11*Zn14*Zn31
     &      -EE/CW*Za21*Zn11*Zn44-EE/CW*Za21*Zn14*Zn41-EE
     &      /SW*Za11*Zn21*Zn34-EE/SW*Za11*Zn24*Zn31+EE/SW
     &      *Za21*Zn21*Zn44+EE/SW*Za21*Zn24*Zn41-Sqrt2*Za11
     &      *Zn41*Zn54*hL-Sqrt2*Za11*Zn44*Zn51*hL-Sqrt2*Za21
     &      *Zn31*Zn54*hL-Sqrt2*Za21*Zn34*Zn51*hL-Sqrt2*Za31
     &      *Zn31*Zn44*hL-Sqrt2*Za31*Zn34*Zn41*hL+2D0*Sqrt2
     &      *Za31*Zn51*Zn54*hK
      AAABC(421) = EE/CW*Za11*Zn11*dconjg(Zn34)+EE/CW*Za11
     &      *Zn14*dconjg(Zn31)-EE/CW*Za21*Zn11*dconjg(Zn44
     &      )-EE/CW*Za21*Zn14*dconjg(Zn41)-EE/SW*Za11*Zn21
     &      *dconjg(Zn34)-EE/SW*Za11*Zn24*dconjg(Zn31)+EE
     &      /SW*Za21*Zn21*dconjg(Zn44)+EE/SW*Za21*Zn24*dconjg(Zn41
     &      )-Sqrt2*Za11*dconjg(Zn41)*dconjg(Zn54)*hL-Sqrt2
     &      *Za11*dconjg(Zn44)*dconjg(Zn51)*hL-Sqrt2*Za21
     &      *dconjg(Zn31)*dconjg(Zn54)*hL-Sqrt2*Za21*dconjg(Zn34
     &      )*dconjg(Zn51)*hL-Sqrt2*Za31*dconjg(Zn31)*dconjg(Zn44
     &      )*hL-Sqrt2*Za31*dconjg(Zn34)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn51)*dconjg(Zn54)*hK
      AAABC(422) = Sqrt2*Za32*Zn51**2*hK+EE/CW*Za12*Zn11
     &      *Zn31-EE/CW*Za22*Zn11*Zn41-EE/SW*Za12*Zn21*Zn31
     &      +EE/SW*Za22*Zn21*Zn41-Sqrt2*Za12*Zn41*Zn51*hL
     &      -Sqrt2*Za22*Zn31*Zn51*hL-Sqrt2*Za32*Zn31*Zn41*hL
      AAABC(423) = Sqrt2*Za32*dconjg(Zn51)**2*hK+EE/CW*Za12
     &      *Zn11*dconjg(Zn31)-EE/CW*Za22*Zn11*dconjg(Zn41
     &      )-EE/SW*Za12*Zn21*dconjg(Zn31)+EE/SW*Za22*Zn21
     &      *dconjg(Zn41)-Sqrt2*Za12*dconjg(Zn41)*dconjg(Zn51
     &      )*hL-Sqrt2*Za22*dconjg(Zn31)*dconjg(Zn51)*hL
     &      -Sqrt2*Za32*dconjg(Zn31)*dconjg(Zn41)*hL
      AAABC(424) = EE/CW*Za12*Zn11*Zn33+EE/CW*Za12*Zn13*Zn31
     &      -EE/CW*Za22*Zn11*Zn43-EE/CW*Za22*Zn13*Zn41-EE
     &      /SW*Za12*Zn21*Zn33-EE/SW*Za12*Zn23*Zn31+EE/SW
     &      *Za22*Zn21*Zn43+EE/SW*Za22*Zn23*Zn41-Sqrt2*Za12
     &      *Zn41*Zn53*hL-Sqrt2*Za12*Zn43*Zn51*hL-Sqrt2*Za22
     &      *Zn31*Zn53*hL-Sqrt2*Za22*Zn33*Zn51*hL-Sqrt2*Za32
     &      *Zn31*Zn43*hL-Sqrt2*Za32*Zn33*Zn41*hL+2D0*Sqrt2
     &      *Za32*Zn51*Zn53*hK
      AAABC(425) = EE/CW*Za12*Zn11*dconjg(Zn33)+EE/CW*Za12
     &      *Zn13*dconjg(Zn31)-EE/CW*Za22*Zn11*dconjg(Zn43
     &      )-EE/CW*Za22*Zn13*dconjg(Zn41)-EE/SW*Za12*Zn21
     &      *dconjg(Zn33)-EE/SW*Za12*Zn23*dconjg(Zn31)+EE
     &      /SW*Za22*Zn21*dconjg(Zn43)+EE/SW*Za22*Zn23*dconjg(Zn41
     &      )-Sqrt2*Za12*dconjg(Zn41)*dconjg(Zn53)*hL-Sqrt2
     &      *Za12*dconjg(Zn43)*dconjg(Zn51)*hL-Sqrt2*Za22
     &      *dconjg(Zn31)*dconjg(Zn53)*hL-Sqrt2*Za22*dconjg(Zn33
     &      )*dconjg(Zn51)*hL-Sqrt2*Za32*dconjg(Zn31)*dconjg(Zn43
     &      )*hL-Sqrt2*Za32*dconjg(Zn33)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn51)*dconjg(Zn53)*hK
      AAABC(426) = EE/CW*Za12*Zn11*Zn34+EE/CW*Za12*Zn14*Zn31
     &      -EE/CW*Za22*Zn11*Zn44-EE/CW*Za22*Zn14*Zn41-EE
     &      /SW*Za12*Zn21*Zn34-EE/SW*Za12*Zn24*Zn31+EE/SW
     &      *Za22*Zn21*Zn44+EE/SW*Za22*Zn24*Zn41-Sqrt2*Za12
     &      *Zn41*Zn54*hL-Sqrt2*Za12*Zn44*Zn51*hL-Sqrt2*Za22
     &      *Zn31*Zn54*hL-Sqrt2*Za22*Zn34*Zn51*hL-Sqrt2*Za32
     &      *Zn31*Zn44*hL-Sqrt2*Za32*Zn34*Zn41*hL+2D0*Sqrt2
     &      *Za32*Zn51*Zn54*hK
      AAABC(427) = EE/CW*Za12*Zn11*dconjg(Zn34)+EE/CW*Za12
     &      *Zn14*dconjg(Zn31)-EE/CW*Za22*Zn11*dconjg(Zn44
     &      )-EE/CW*Za22*Zn14*dconjg(Zn41)-EE/SW*Za12*Zn21
     &      *dconjg(Zn34)-EE/SW*Za12*Zn24*dconjg(Zn31)+EE
     &      /SW*Za22*Zn21*dconjg(Zn44)+EE/SW*Za22*Zn24*dconjg(Zn41
     &      )-Sqrt2*Za12*dconjg(Zn41)*dconjg(Zn54)*hL-Sqrt2
     &      *Za12*dconjg(Zn44)*dconjg(Zn51)*hL-Sqrt2*Za22
     &      *dconjg(Zn31)*dconjg(Zn54)*hL-Sqrt2*Za22*dconjg(Zn34
     &      )*dconjg(Zn51)*hL-Sqrt2*Za32*dconjg(Zn31)*dconjg(Zn44
     &      )*hL-Sqrt2*Za32*dconjg(Zn34)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn51)*dconjg(Zn54)*hK
      AAABC(428) = EE/CW*Za11*Zn11*Zn32+EE/CW*Za11*Zn12*Zn31
     &      -EE/CW*Za21*Zn11*Zn42-EE/CW*Za21*Zn12*Zn41-EE
     &      /SW*Za11*Zn21*Zn32-EE/SW*Za11*Zn22*Zn31+EE/SW
     &      *Za21*Zn21*Zn42+EE/SW*Za21*Zn22*Zn41-Sqrt2*Za11
     &      *Zn41*Zn52*hL-Sqrt2*Za11*Zn42*Zn51*hL-Sqrt2*Za21
     &      *Zn31*Zn52*hL-Sqrt2*Za21*Zn32*Zn51*hL-Sqrt2*Za31
     &      *Zn31*Zn42*hL-Sqrt2*Za31*Zn32*Zn41*hL+2D0*Sqrt2
     &      *Za31*Zn51*Zn52*hK
      AAABC(429) = EE/CW*Za11*Zn11*dconjg(Zn32)+EE/CW*Za11
     &      *Zn12*dconjg(Zn31)-EE/CW*Za21*Zn11*dconjg(Zn42
     &      )-EE/CW*Za21*Zn12*dconjg(Zn41)-EE/SW*Za11*Zn21
     &      *dconjg(Zn32)-EE/SW*Za11*Zn22*dconjg(Zn31)+EE
     &      /SW*Za21*Zn21*dconjg(Zn42)+EE/SW*Za21*Zn22*dconjg(Zn41
     &      )-Sqrt2*Za11*dconjg(Zn41)*dconjg(Zn52)*hL-Sqrt2
     &      *Za11*dconjg(Zn42)*dconjg(Zn51)*hL-Sqrt2*Za21
     &      *dconjg(Zn31)*dconjg(Zn52)*hL-Sqrt2*Za21*dconjg(Zn32
     &      )*dconjg(Zn51)*hL-Sqrt2*Za31*dconjg(Zn31)*dconjg(Zn42
     &      )*hL-Sqrt2*Za31*dconjg(Zn32)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn51)*dconjg(Zn52)*hK
      AAABC(430) = EE/CW*Zh13*Zn11*Zn33+EE/CW*Zh13*Zn13*Zn31
     &      -EE/CW*Zh23*Zn11*Zn43-EE/CW*Zh23*Zn13*Zn41-EE
     &      /SW*Zh13*Zn21*Zn33-EE/SW*Zh13*Zn23*Zn31+EE/SW
     &      *Zh23*Zn21*Zn43+EE/SW*Zh23*Zn23*Zn41+Sqrt2*Zh13
     &      *Zn41*Zn53*hL+Sqrt2*Zh13*Zn43*Zn51*hL+Sqrt2*Zh23
     &      *Zn31*Zn53*hL+Sqrt2*Zh23*Zn33*Zn51*hL+Sqrt2*Zh33
     &      *Zn31*Zn43*hL+Sqrt2*Zh33*Zn33*Zn41*hL-2D0*Sqrt2
     &      *Zh33*Zn51*Zn53*hK
      AAABC(431) = EE/CW*Zh13*Zn11*dconjg(Zn33)+EE/CW*Zh13
     &      *Zn13*dconjg(Zn31)-EE/CW*Zh23*Zn11*dconjg(Zn43
     &      )-EE/CW*Zh23*Zn13*dconjg(Zn41)-EE/SW*Zh13*Zn21
     &      *dconjg(Zn33)-EE/SW*Zh13*Zn23*dconjg(Zn31)+EE
     &      /SW*Zh23*Zn21*dconjg(Zn43)+EE/SW*Zh23*Zn23*dconjg(Zn41
     &      )+Sqrt2*Zh13*dconjg(Zn41)*dconjg(Zn53)*hL+Sqrt2
     &      *Zh13*dconjg(Zn43)*dconjg(Zn51)*hL+Sqrt2*Zh23
     &      *dconjg(Zn31)*dconjg(Zn53)*hL+Sqrt2*Zh23*dconjg(Zn33
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh33*dconjg(Zn31)*dconjg(Zn43
     &      )*hL+Sqrt2*Zh33*dconjg(Zn33)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn51)*dconjg(Zn53)*hK
      AAABC(432) = EE/CW*Zh11*Zn11*Zn33+EE/CW*Zh11*Zn13*Zn31
     &      -EE/CW*Zh21*Zn11*Zn43-EE/CW*Zh21*Zn13*Zn41-EE
     &      /SW*Zh11*Zn21*Zn33-EE/SW*Zh11*Zn23*Zn31+EE/SW
     &      *Zh21*Zn21*Zn43+EE/SW*Zh21*Zn23*Zn41+Sqrt2*Zh11
     &      *Zn41*Zn53*hL+Sqrt2*Zh11*Zn43*Zn51*hL+Sqrt2*Zh21
     &      *Zn31*Zn53*hL+Sqrt2*Zh21*Zn33*Zn51*hL+Sqrt2*Zh31
     &      *Zn31*Zn43*hL+Sqrt2*Zh31*Zn33*Zn41*hL-2D0*Sqrt2
     &      *Zh31*Zn51*Zn53*hK
      AAABC(433) = EE/CW*Zh11*Zn11*dconjg(Zn33)+EE/CW*Zh11
     &      *Zn13*dconjg(Zn31)-EE/CW*Zh21*Zn11*dconjg(Zn43
     &      )-EE/CW*Zh21*Zn13*dconjg(Zn41)-EE/SW*Zh11*Zn21
     &      *dconjg(Zn33)-EE/SW*Zh11*Zn23*dconjg(Zn31)+EE
     &      /SW*Zh21*Zn21*dconjg(Zn43)+EE/SW*Zh21*Zn23*dconjg(Zn41
     &      )+Sqrt2*Zh11*dconjg(Zn41)*dconjg(Zn53)*hL+Sqrt2
     &      *Zh11*dconjg(Zn43)*dconjg(Zn51)*hL+Sqrt2*Zh21
     &      *dconjg(Zn31)*dconjg(Zn53)*hL+Sqrt2*Zh21*dconjg(Zn33
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh31*dconjg(Zn31)*dconjg(Zn43
     &      )*hL+Sqrt2*Zh31*dconjg(Zn33)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn51)*dconjg(Zn53)*hK
      AAABC(434) = EE/CW*Zh12*Zn11*Zn34+EE/CW*Zh12*Zn14*Zn31
     &      -EE/CW*Zh22*Zn11*Zn44-EE/CW*Zh22*Zn14*Zn41-EE
     &      /SW*Zh12*Zn21*Zn34-EE/SW*Zh12*Zn24*Zn31+EE/SW
     &      *Zh22*Zn21*Zn44+EE/SW*Zh22*Zn24*Zn41+Sqrt2*Zh12
     &      *Zn41*Zn54*hL+Sqrt2*Zh12*Zn44*Zn51*hL+Sqrt2*Zh22
     &      *Zn31*Zn54*hL+Sqrt2*Zh22*Zn34*Zn51*hL+Sqrt2*Zh32
     &      *Zn31*Zn44*hL+Sqrt2*Zh32*Zn34*Zn41*hL-2D0*Sqrt2
     &      *Zh32*Zn51*Zn54*hK
      AAABC(435) = EE/CW*Zh12*Zn11*dconjg(Zn34)+EE/CW*Zh12
     &      *Zn14*dconjg(Zn31)-EE/CW*Zh22*Zn11*dconjg(Zn44
     &      )-EE/CW*Zh22*Zn14*dconjg(Zn41)-EE/SW*Zh12*Zn21
     &      *dconjg(Zn34)-EE/SW*Zh12*Zn24*dconjg(Zn31)+EE
     &      /SW*Zh22*Zn21*dconjg(Zn44)+EE/SW*Zh22*Zn24*dconjg(Zn41
     &      )+Sqrt2*Zh12*dconjg(Zn41)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh12*dconjg(Zn44)*dconjg(Zn51)*hL+Sqrt2*Zh22
     &      *dconjg(Zn31)*dconjg(Zn54)*hL+Sqrt2*Zh22*dconjg(Zn34
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh32*dconjg(Zn31)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh32*dconjg(Zn34)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn51)*dconjg(Zn54)*hK
      AAABC(436) = EE/CW*Zh12*Zn11*Zn32+EE/CW*Zh12*Zn12*Zn31
     &      -EE/CW*Zh22*Zn11*Zn42-EE/CW*Zh22*Zn12*Zn41-EE
     &      /SW*Zh12*Zn21*Zn32-EE/SW*Zh12*Zn22*Zn31+EE/SW
     &      *Zh22*Zn21*Zn42+EE/SW*Zh22*Zn22*Zn41+Sqrt2*Zh12
     &      *Zn41*Zn52*hL+Sqrt2*Zh12*Zn42*Zn51*hL+Sqrt2*Zh22
     &      *Zn31*Zn52*hL+Sqrt2*Zh22*Zn32*Zn51*hL+Sqrt2*Zh32
     &      *Zn31*Zn42*hL+Sqrt2*Zh32*Zn32*Zn41*hL-2D0*Sqrt2
     &      *Zh32*Zn51*Zn52*hK
      AAABC(437) = EE/CW*Zh12*Zn11*dconjg(Zn32)+EE/CW*Zh12
     &      *Zn12*dconjg(Zn31)-EE/CW*Zh22*Zn11*dconjg(Zn42
     &      )-EE/CW*Zh22*Zn12*dconjg(Zn41)-EE/SW*Zh12*Zn21
     &      *dconjg(Zn32)-EE/SW*Zh12*Zn22*dconjg(Zn31)+EE
     &      /SW*Zh22*Zn21*dconjg(Zn42)+EE/SW*Zh22*Zn22*dconjg(Zn41
     &      )+Sqrt2*Zh12*dconjg(Zn41)*dconjg(Zn52)*hL+Sqrt2
     &      *Zh12*dconjg(Zn42)*dconjg(Zn51)*hL+Sqrt2*Zh22
     &      *dconjg(Zn31)*dconjg(Zn52)*hL+Sqrt2*Zh22*dconjg(Zn32
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh32*dconjg(Zn31)*dconjg(Zn42
     &      )*hL+Sqrt2*Zh32*dconjg(Zn32)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn51)*dconjg(Zn52)*hK
      AAABC(438) = EE/CW*Zh13*Zn11*Zn34+EE/CW*Zh13*Zn14*Zn31
     &      -EE/CW*Zh23*Zn11*Zn44-EE/CW*Zh23*Zn14*Zn41-EE
     &      /SW*Zh13*Zn21*Zn34-EE/SW*Zh13*Zn24*Zn31+EE/SW
     &      *Zh23*Zn21*Zn44+EE/SW*Zh23*Zn24*Zn41+Sqrt2*Zh13
     &      *Zn41*Zn54*hL+Sqrt2*Zh13*Zn44*Zn51*hL+Sqrt2*Zh23
     &      *Zn31*Zn54*hL+Sqrt2*Zh23*Zn34*Zn51*hL+Sqrt2*Zh33
     &      *Zn31*Zn44*hL+Sqrt2*Zh33*Zn34*Zn41*hL-2D0*Sqrt2
     &      *Zh33*Zn51*Zn54*hK
      AAABC(439) = EE/CW*Zh13*Zn11*dconjg(Zn34)+EE/CW*Zh13
     &      *Zn14*dconjg(Zn31)-EE/CW*Zh23*Zn11*dconjg(Zn44
     &      )-EE/CW*Zh23*Zn14*dconjg(Zn41)-EE/SW*Zh13*Zn21
     &      *dconjg(Zn34)-EE/SW*Zh13*Zn24*dconjg(Zn31)+EE
     &      /SW*Zh23*Zn21*dconjg(Zn44)+EE/SW*Zh23*Zn24*dconjg(Zn41
     &      )+Sqrt2*Zh13*dconjg(Zn41)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh13*dconjg(Zn44)*dconjg(Zn51)*hL+Sqrt2*Zh23
     &      *dconjg(Zn31)*dconjg(Zn54)*hL+Sqrt2*Zh23*dconjg(Zn34
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh33*dconjg(Zn31)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh33*dconjg(Zn34)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn51)*dconjg(Zn54)*hK
      AAABC(440) = EE/CW*Zh11*Zn11*Zn32+EE/CW*Zh11*Zn12*Zn31
     &      -EE/CW*Zh21*Zn11*Zn42-EE/CW*Zh21*Zn12*Zn41-EE
     &      /SW*Zh11*Zn21*Zn32-EE/SW*Zh11*Zn22*Zn31+EE/SW
     &      *Zh21*Zn21*Zn42+EE/SW*Zh21*Zn22*Zn41+Sqrt2*Zh11
     &      *Zn41*Zn52*hL+Sqrt2*Zh11*Zn42*Zn51*hL+Sqrt2*Zh21
     &      *Zn31*Zn52*hL+Sqrt2*Zh21*Zn32*Zn51*hL+Sqrt2*Zh31
     &      *Zn31*Zn42*hL+Sqrt2*Zh31*Zn32*Zn41*hL-2D0*Sqrt2
     &      *Zh31*Zn51*Zn52*hK
      AAABC(441) = EE/CW*Zh11*Zn11*dconjg(Zn32)+EE/CW*Zh11
     &      *Zn12*dconjg(Zn31)-EE/CW*Zh21*Zn11*dconjg(Zn42
     &      )-EE/CW*Zh21*Zn12*dconjg(Zn41)-EE/SW*Zh11*Zn21
     &      *dconjg(Zn32)-EE/SW*Zh11*Zn22*dconjg(Zn31)+EE
     &      /SW*Zh21*Zn21*dconjg(Zn42)+EE/SW*Zh21*Zn22*dconjg(Zn41
     &      )+Sqrt2*Zh11*dconjg(Zn41)*dconjg(Zn52)*hL+Sqrt2
     &      *Zh11*dconjg(Zn42)*dconjg(Zn51)*hL+Sqrt2*Zh21
     &      *dconjg(Zn31)*dconjg(Zn52)*hL+Sqrt2*Zh21*dconjg(Zn32
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh31*dconjg(Zn31)*dconjg(Zn42
     &      )*hL+Sqrt2*Zh31*dconjg(Zn32)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn51)*dconjg(Zn52)*hK
      AAABC(442) = Sqrt2*Zh33*Zn51**2*hK-EE/CW*Zh13*Zn11
     &      *Zn31+EE/CW*Zh23*Zn11*Zn41+EE/SW*Zh13*Zn21*Zn31
     &      -EE/SW*Zh23*Zn21*Zn41-Sqrt2*Zh13*Zn41*Zn51*hL
     &      -Sqrt2*Zh23*Zn31*Zn51*hL-Sqrt2*Zh33*Zn31*Zn41*hL
      AAABC(443) = Sqrt2*Zh33*dconjg(Zn51)**2*hK-EE/CW*Zh13
     &      *Zn11*dconjg(Zn31)+EE/CW*Zh23*Zn11*dconjg(Zn41
     &      )+EE/SW*Zh13*Zn21*dconjg(Zn31)-EE/SW*Zh23*Zn21
     &      *dconjg(Zn41)-Sqrt2*Zh13*dconjg(Zn41)*dconjg(Zn51
     &      )*hL-Sqrt2*Zh23*dconjg(Zn31)*dconjg(Zn51)*hL
     &      -Sqrt2*Zh33*dconjg(Zn31)*dconjg(Zn41)*hL
      AAABC(444) = EE/CW*Zh13*Zn11*Zn32+EE/CW*Zh13*Zn12*Zn31
     &      -EE/CW*Zh23*Zn11*Zn42-EE/CW*Zh23*Zn12*Zn41-EE
     &      /SW*Zh13*Zn21*Zn32-EE/SW*Zh13*Zn22*Zn31+EE/SW
     &      *Zh23*Zn21*Zn42+EE/SW*Zh23*Zn22*Zn41+Sqrt2*Zh13
     &      *Zn41*Zn52*hL+Sqrt2*Zh13*Zn42*Zn51*hL+Sqrt2*Zh23
     &      *Zn31*Zn52*hL+Sqrt2*Zh23*Zn32*Zn51*hL+Sqrt2*Zh33
     &      *Zn31*Zn42*hL+Sqrt2*Zh33*Zn32*Zn41*hL-2D0*Sqrt2
     &      *Zh33*Zn51*Zn52*hK
      AAABC(445) = EE/CW*Zh13*Zn11*dconjg(Zn32)+EE/CW*Zh13
     &      *Zn12*dconjg(Zn31)-EE/CW*Zh23*Zn11*dconjg(Zn42
     &      )-EE/CW*Zh23*Zn12*dconjg(Zn41)-EE/SW*Zh13*Zn21
     &      *dconjg(Zn32)-EE/SW*Zh13*Zn22*dconjg(Zn31)+EE
     &      /SW*Zh23*Zn21*dconjg(Zn42)+EE/SW*Zh23*Zn22*dconjg(Zn41
     &      )+Sqrt2*Zh13*dconjg(Zn41)*dconjg(Zn52)*hL+Sqrt2
     &      *Zh13*dconjg(Zn42)*dconjg(Zn51)*hL+Sqrt2*Zh23
     &      *dconjg(Zn31)*dconjg(Zn52)*hL+Sqrt2*Zh23*dconjg(Zn32
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh33*dconjg(Zn31)*dconjg(Zn42
     &      )*hL+Sqrt2*Zh33*dconjg(Zn32)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn51)*dconjg(Zn52)*hK
      AAABC(446) = EE/CW*Zh12*Zn11*Zn33+EE/CW*Zh12*Zn13*Zn31
     &      -EE/CW*Zh22*Zn11*Zn43-EE/CW*Zh22*Zn13*Zn41-EE
     &      /SW*Zh12*Zn21*Zn33-EE/SW*Zh12*Zn23*Zn31+EE/SW
     &      *Zh22*Zn21*Zn43+EE/SW*Zh22*Zn23*Zn41+Sqrt2*Zh12
     &      *Zn41*Zn53*hL+Sqrt2*Zh12*Zn43*Zn51*hL+Sqrt2*Zh22
     &      *Zn31*Zn53*hL+Sqrt2*Zh22*Zn33*Zn51*hL+Sqrt2*Zh32
     &      *Zn31*Zn43*hL+Sqrt2*Zh32*Zn33*Zn41*hL-2D0*Sqrt2
     &      *Zh32*Zn51*Zn53*hK
      AAABC(447) = EE/CW*Zh12*Zn11*dconjg(Zn33)+EE/CW*Zh12
     &      *Zn13*dconjg(Zn31)-EE/CW*Zh22*Zn11*dconjg(Zn43
     &      )-EE/CW*Zh22*Zn13*dconjg(Zn41)-EE/SW*Zh12*Zn21
     &      *dconjg(Zn33)-EE/SW*Zh12*Zn23*dconjg(Zn31)+EE
     &      /SW*Zh22*Zn21*dconjg(Zn43)+EE/SW*Zh22*Zn23*dconjg(Zn41
     &      )+Sqrt2*Zh12*dconjg(Zn41)*dconjg(Zn53)*hL+Sqrt2
     &      *Zh12*dconjg(Zn43)*dconjg(Zn51)*hL+Sqrt2*Zh22
     &      *dconjg(Zn31)*dconjg(Zn53)*hL+Sqrt2*Zh22*dconjg(Zn33
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh32*dconjg(Zn31)*dconjg(Zn43
     &      )*hL+Sqrt2*Zh32*dconjg(Zn33)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn51)*dconjg(Zn53)*hK
      AAABC(448) = Sqrt2*Zh31*Zn51**2*hK-EE/CW*Zh11*Zn11
     &      *Zn31+EE/CW*Zh21*Zn11*Zn41+EE/SW*Zh11*Zn21*Zn31
     &      -EE/SW*Zh21*Zn21*Zn41-Sqrt2*Zh11*Zn41*Zn51*hL
     &      -Sqrt2*Zh21*Zn31*Zn51*hL-Sqrt2*Zh31*Zn31*Zn41*hL
      AAABC(449) = Sqrt2*Zh31*dconjg(Zn51)**2*hK-EE/CW*Zh11
     &      *Zn11*dconjg(Zn31)+EE/CW*Zh21*Zn11*dconjg(Zn41
     &      )+EE/SW*Zh11*Zn21*dconjg(Zn31)-EE/SW*Zh21*Zn21
     &      *dconjg(Zn41)-Sqrt2*Zh11*dconjg(Zn41)*dconjg(Zn51
     &      )*hL-Sqrt2*Zh21*dconjg(Zn31)*dconjg(Zn51)*hL
     &      -Sqrt2*Zh31*dconjg(Zn31)*dconjg(Zn41)*hL
      AAABC(450) = EE/CW*Zh11*Zn11*Zn34+EE/CW*Zh11*Zn14*Zn31
     &      -EE/CW*Zh21*Zn11*Zn44-EE/CW*Zh21*Zn14*Zn41-EE
     &      /SW*Zh11*Zn21*Zn34-EE/SW*Zh11*Zn24*Zn31+EE/SW
     &      *Zh21*Zn21*Zn44+EE/SW*Zh21*Zn24*Zn41+Sqrt2*Zh11
     &      *Zn41*Zn54*hL+Sqrt2*Zh11*Zn44*Zn51*hL+Sqrt2*Zh21
     &      *Zn31*Zn54*hL+Sqrt2*Zh21*Zn34*Zn51*hL+Sqrt2*Zh31
     &      *Zn31*Zn44*hL+Sqrt2*Zh31*Zn34*Zn41*hL-2D0*Sqrt2
     &      *Zh31*Zn51*Zn54*hK
      AAABC(451) = EE/CW*Zh11*Zn11*dconjg(Zn34)+EE/CW*Zh11
     &      *Zn14*dconjg(Zn31)-EE/CW*Zh21*Zn11*dconjg(Zn44
     &      )-EE/CW*Zh21*Zn14*dconjg(Zn41)-EE/SW*Zh11*Zn21
     &      *dconjg(Zn34)-EE/SW*Zh11*Zn24*dconjg(Zn31)+EE
     &      /SW*Zh21*Zn21*dconjg(Zn44)+EE/SW*Zh21*Zn24*dconjg(Zn41
     &      )+Sqrt2*Zh11*dconjg(Zn41)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh11*dconjg(Zn44)*dconjg(Zn51)*hL+Sqrt2*Zh21
     &      *dconjg(Zn31)*dconjg(Zn54)*hL+Sqrt2*Zh21*dconjg(Zn34
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh31*dconjg(Zn31)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh31*dconjg(Zn34)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn51)*dconjg(Zn54)*hK
      AAABC(452) = Sqrt2*Zh32*Zn51**2*hK-EE/CW*Zh12*Zn11
     &      *Zn31+EE/CW*Zh22*Zn11*Zn41+EE/SW*Zh12*Zn21*Zn31
     &      -EE/SW*Zh22*Zn21*Zn41-Sqrt2*Zh12*Zn41*Zn51*hL
     &      -Sqrt2*Zh22*Zn31*Zn51*hL-Sqrt2*Zh32*Zn31*Zn41*hL
      AAABC(453) = Sqrt2*Zh32*dconjg(Zn51)**2*hK-EE/CW*Zh12
     &      *Zn11*dconjg(Zn31)+EE/CW*Zh22*Zn11*dconjg(Zn41
     &      )+EE/SW*Zh12*Zn21*dconjg(Zn31)-EE/SW*Zh22*Zn21
     &      *dconjg(Zn41)-Sqrt2*Zh12*dconjg(Zn41)*dconjg(Zn51
     &      )*hL-Sqrt2*Zh22*dconjg(Zn31)*dconjg(Zn51)*hL
     &      -Sqrt2*Zh32*dconjg(Zn31)*dconjg(Zn41)*hL
      AAABC(454) = EE/CW*Za13*Zn11*Zn35+EE/CW*Za13*Zn15*Zn31
     &      -EE/CW*Za23*Zn11*Zn45-EE/CW*Za23*Zn15*Zn41-EE
     &      /SW*Za13*Zn21*Zn35-EE/SW*Za13*Zn25*Zn31+EE/SW
     &      *Za23*Zn21*Zn45+EE/SW*Za23*Zn25*Zn41-Sqrt2*Za13
     &      *Zn41*Zn55*hL-Sqrt2*Za13*Zn45*Zn51*hL-Sqrt2*Za23
     &      *Zn31*Zn55*hL-Sqrt2*Za23*Zn35*Zn51*hL-Sqrt2*Za33
     &      *Zn31*Zn45*hL-Sqrt2*Za33*Zn35*Zn41*hL+2D0*Sqrt2
     &      *Za33*Zn51*Zn55*hK
      AAABC(455) = EE/CW*Za13*Zn11*dconjg(Zn35)+EE/CW*Za13
     &      *Zn15*dconjg(Zn31)-EE/CW*Za23*Zn11*dconjg(Zn45
     &      )-EE/CW*Za23*Zn15*dconjg(Zn41)-EE/SW*Za13*Zn21
     &      *dconjg(Zn35)-EE/SW*Za13*Zn25*dconjg(Zn31)+EE
     &      /SW*Za23*Zn21*dconjg(Zn45)+EE/SW*Za23*Zn25*dconjg(Zn41
     &      )-Sqrt2*Za13*dconjg(Zn41)*dconjg(Zn55)*hL-Sqrt2
     &      *Za13*dconjg(Zn45)*dconjg(Zn51)*hL-Sqrt2*Za23
     &      *dconjg(Zn31)*dconjg(Zn55)*hL-Sqrt2*Za23*dconjg(Zn35
     &      )*dconjg(Zn51)*hL-Sqrt2*Za33*dconjg(Zn31)*dconjg(Zn45
     &      )*hL-Sqrt2*Za33*dconjg(Zn35)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn51)*dconjg(Zn55)*hK
      AAABC(456) = EE/CW*Za11*Zn11*Zn35+EE/CW*Za11*Zn15*Zn31
     &      -EE/CW*Za21*Zn11*Zn45-EE/CW*Za21*Zn15*Zn41-EE
     &      /SW*Za11*Zn21*Zn35-EE/SW*Za11*Zn25*Zn31+EE/SW
     &      *Za21*Zn21*Zn45+EE/SW*Za21*Zn25*Zn41-Sqrt2*Za11
     &      *Zn41*Zn55*hL-Sqrt2*Za11*Zn45*Zn51*hL-Sqrt2*Za21
     &      *Zn31*Zn55*hL-Sqrt2*Za21*Zn35*Zn51*hL-Sqrt2*Za31
     &      *Zn31*Zn45*hL-Sqrt2*Za31*Zn35*Zn41*hL+2D0*Sqrt2
     &      *Za31*Zn51*Zn55*hK
      AAABC(457) = EE/CW*Za11*Zn11*dconjg(Zn35)+EE/CW*Za11
     &      *Zn15*dconjg(Zn31)-EE/CW*Za21*Zn11*dconjg(Zn45
     &      )-EE/CW*Za21*Zn15*dconjg(Zn41)-EE/SW*Za11*Zn21
     &      *dconjg(Zn35)-EE/SW*Za11*Zn25*dconjg(Zn31)+EE
     &      /SW*Za21*Zn21*dconjg(Zn45)+EE/SW*Za21*Zn25*dconjg(Zn41
     &      )-Sqrt2*Za11*dconjg(Zn41)*dconjg(Zn55)*hL-Sqrt2
     &      *Za11*dconjg(Zn45)*dconjg(Zn51)*hL-Sqrt2*Za21
     &      *dconjg(Zn31)*dconjg(Zn55)*hL-Sqrt2*Za21*dconjg(Zn35
     &      )*dconjg(Zn51)*hL-Sqrt2*Za31*dconjg(Zn31)*dconjg(Zn45
     &      )*hL-Sqrt2*Za31*dconjg(Zn35)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn51)*dconjg(Zn55)*hK
      AAABC(458) = EE/CW*Za12*Zn11*Zn35+EE/CW*Za12*Zn15*Zn31
     &      -EE/CW*Za22*Zn11*Zn45-EE/CW*Za22*Zn15*Zn41-EE
     &      /SW*Za12*Zn21*Zn35-EE/SW*Za12*Zn25*Zn31+EE/SW
     &      *Za22*Zn21*Zn45+EE/SW*Za22*Zn25*Zn41-Sqrt2*Za12
     &      *Zn41*Zn55*hL-Sqrt2*Za12*Zn45*Zn51*hL-Sqrt2*Za22
     &      *Zn31*Zn55*hL-Sqrt2*Za22*Zn35*Zn51*hL-Sqrt2*Za32
     &      *Zn31*Zn45*hL-Sqrt2*Za32*Zn35*Zn41*hL+2D0*Sqrt2
     &      *Za32*Zn51*Zn55*hK
      AAABC(459) = EE/CW*Za12*Zn11*dconjg(Zn35)+EE/CW*Za12
     &      *Zn15*dconjg(Zn31)-EE/CW*Za22*Zn11*dconjg(Zn45
     &      )-EE/CW*Za22*Zn15*dconjg(Zn41)-EE/SW*Za12*Zn21
     &      *dconjg(Zn35)-EE/SW*Za12*Zn25*dconjg(Zn31)+EE
     &      /SW*Za22*Zn21*dconjg(Zn45)+EE/SW*Za22*Zn25*dconjg(Zn41
     &      )-Sqrt2*Za12*dconjg(Zn41)*dconjg(Zn55)*hL-Sqrt2
     &      *Za12*dconjg(Zn45)*dconjg(Zn51)*hL-Sqrt2*Za22
     &      *dconjg(Zn31)*dconjg(Zn55)*hL-Sqrt2*Za22*dconjg(Zn35
     &      )*dconjg(Zn51)*hL-Sqrt2*Za32*dconjg(Zn31)*dconjg(Zn45
     &      )*hL-Sqrt2*Za32*dconjg(Zn35)*dconjg(Zn41)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn51)*dconjg(Zn55)*hK
      AAABC(460) = EE/CW*Zh13*Zn11*Zn35+EE/CW*Zh13*Zn15*Zn31
     &      -EE/CW*Zh23*Zn11*Zn45-EE/CW*Zh23*Zn15*Zn41-EE
     &      /SW*Zh13*Zn21*Zn35-EE/SW*Zh13*Zn25*Zn31+EE/SW
     &      *Zh23*Zn21*Zn45+EE/SW*Zh23*Zn25*Zn41+Sqrt2*Zh13
     &      *Zn41*Zn55*hL+Sqrt2*Zh13*Zn45*Zn51*hL+Sqrt2*Zh23
     &      *Zn31*Zn55*hL+Sqrt2*Zh23*Zn35*Zn51*hL+Sqrt2*Zh33
     &      *Zn31*Zn45*hL+Sqrt2*Zh33*Zn35*Zn41*hL-2D0*Sqrt2
     &      *Zh33*Zn51*Zn55*hK
      AAABC(461) = EE/CW*Zh13*Zn11*dconjg(Zn35)+EE/CW*Zh13
     &      *Zn15*dconjg(Zn31)-EE/CW*Zh23*Zn11*dconjg(Zn45
     &      )-EE/CW*Zh23*Zn15*dconjg(Zn41)-EE/SW*Zh13*Zn21
     &      *dconjg(Zn35)-EE/SW*Zh13*Zn25*dconjg(Zn31)+EE
     &      /SW*Zh23*Zn21*dconjg(Zn45)+EE/SW*Zh23*Zn25*dconjg(Zn41
     &      )+Sqrt2*Zh13*dconjg(Zn41)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh13*dconjg(Zn45)*dconjg(Zn51)*hL+Sqrt2*Zh23
     &      *dconjg(Zn31)*dconjg(Zn55)*hL+Sqrt2*Zh23*dconjg(Zn35
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh33*dconjg(Zn31)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh33*dconjg(Zn35)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn51)*dconjg(Zn55)*hK
      AAABC(462) = EE/CW*Zh12*Zn11*Zn35+EE/CW*Zh12*Zn15*Zn31
     &      -EE/CW*Zh22*Zn11*Zn45-EE/CW*Zh22*Zn15*Zn41-EE
     &      /SW*Zh12*Zn21*Zn35-EE/SW*Zh12*Zn25*Zn31+EE/SW
     &      *Zh22*Zn21*Zn45+EE/SW*Zh22*Zn25*Zn41+Sqrt2*Zh12
     &      *Zn41*Zn55*hL+Sqrt2*Zh12*Zn45*Zn51*hL+Sqrt2*Zh22
     &      *Zn31*Zn55*hL+Sqrt2*Zh22*Zn35*Zn51*hL+Sqrt2*Zh32
     &      *Zn31*Zn45*hL+Sqrt2*Zh32*Zn35*Zn41*hL-2D0*Sqrt2
     &      *Zh32*Zn51*Zn55*hK
      AAABC(463) = EE/CW*Zh12*Zn11*dconjg(Zn35)+EE/CW*Zh12
     &      *Zn15*dconjg(Zn31)-EE/CW*Zh22*Zn11*dconjg(Zn45
     &      )-EE/CW*Zh22*Zn15*dconjg(Zn41)-EE/SW*Zh12*Zn21
     &      *dconjg(Zn35)-EE/SW*Zh12*Zn25*dconjg(Zn31)+EE
     &      /SW*Zh22*Zn21*dconjg(Zn45)+EE/SW*Zh22*Zn25*dconjg(Zn41
     &      )+Sqrt2*Zh12*dconjg(Zn41)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh12*dconjg(Zn45)*dconjg(Zn51)*hL+Sqrt2*Zh22
     &      *dconjg(Zn31)*dconjg(Zn55)*hL+Sqrt2*Zh22*dconjg(Zn35
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh32*dconjg(Zn31)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh32*dconjg(Zn35)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn51)*dconjg(Zn55)*hK
      AAABC(464) = EE/CW*Zh11*Zn11*Zn35+EE/CW*Zh11*Zn15*Zn31
     &      -EE/CW*Zh21*Zn11*Zn45-EE/CW*Zh21*Zn15*Zn41-EE
     &      /SW*Zh11*Zn21*Zn35-EE/SW*Zh11*Zn25*Zn31+EE/SW
     &      *Zh21*Zn21*Zn45+EE/SW*Zh21*Zn25*Zn41+Sqrt2*Zh11
     &      *Zn41*Zn55*hL+Sqrt2*Zh11*Zn45*Zn51*hL+Sqrt2*Zh21
     &      *Zn31*Zn55*hL+Sqrt2*Zh21*Zn35*Zn51*hL+Sqrt2*Zh31
     &      *Zn31*Zn45*hL+Sqrt2*Zh31*Zn35*Zn41*hL-2D0*Sqrt2
     &      *Zh31*Zn51*Zn55*hK
      AAABC(465) = EE/CW*Zh11*Zn11*dconjg(Zn35)+EE/CW*Zh11
     &      *Zn15*dconjg(Zn31)-EE/CW*Zh21*Zn11*dconjg(Zn45
     &      )-EE/CW*Zh21*Zn15*dconjg(Zn41)-EE/SW*Zh11*Zn21
     &      *dconjg(Zn35)-EE/SW*Zh11*Zn25*dconjg(Zn31)+EE
     &      /SW*Zh21*Zn21*dconjg(Zn45)+EE/SW*Zh21*Zn25*dconjg(Zn41
     &      )+Sqrt2*Zh11*dconjg(Zn41)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh11*dconjg(Zn45)*dconjg(Zn51)*hL+Sqrt2*Zh21
     &      *dconjg(Zn31)*dconjg(Zn55)*hL+Sqrt2*Zh21*dconjg(Zn35
     &      )*dconjg(Zn51)*hL+Sqrt2*Zh31*dconjg(Zn31)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh31*dconjg(Zn35)*dconjg(Zn41)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn51)*dconjg(Zn55)*hK
      AAABC(466) = EE/CW*Za13*Zn12*Zn33+EE/CW*Za13*Zn13*Zn32
     &      -EE/CW*Za23*Zn12*Zn43-EE/CW*Za23*Zn13*Zn42-EE
     &      /SW*Za13*Zn22*Zn33-EE/SW*Za13*Zn23*Zn32+EE/SW
     &      *Za23*Zn22*Zn43+EE/SW*Za23*Zn23*Zn42-Sqrt2*Za13
     &      *Zn42*Zn53*hL-Sqrt2*Za13*Zn43*Zn52*hL-Sqrt2*Za23
     &      *Zn32*Zn53*hL-Sqrt2*Za23*Zn33*Zn52*hL-Sqrt2*Za33
     &      *Zn32*Zn43*hL-Sqrt2*Za33*Zn33*Zn42*hL+2D0*Sqrt2
     &      *Za33*Zn52*Zn53*hK
      AAABC(467) = EE/CW*Za13*Zn12*dconjg(Zn33)+EE/CW*Za13
     &      *Zn13*dconjg(Zn32)-EE/CW*Za23*Zn12*dconjg(Zn43
     &      )-EE/CW*Za23*Zn13*dconjg(Zn42)-EE/SW*Za13*Zn22
     &      *dconjg(Zn33)-EE/SW*Za13*Zn23*dconjg(Zn32)+EE
     &      /SW*Za23*Zn22*dconjg(Zn43)+EE/SW*Za23*Zn23*dconjg(Zn42
     &      )-Sqrt2*Za13*dconjg(Zn42)*dconjg(Zn53)*hL-Sqrt2
     &      *Za13*dconjg(Zn43)*dconjg(Zn52)*hL-Sqrt2*Za23
     &      *dconjg(Zn32)*dconjg(Zn53)*hL-Sqrt2*Za23*dconjg(Zn33
     &      )*dconjg(Zn52)*hL-Sqrt2*Za33*dconjg(Zn32)*dconjg(Zn43
     &      )*hL-Sqrt2*Za33*dconjg(Zn33)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn52)*dconjg(Zn53)*hK
      AAABC(468) = EE/CW*Za13*Zn12*Zn34+EE/CW*Za13*Zn14*Zn32
     &      -EE/CW*Za23*Zn12*Zn44-EE/CW*Za23*Zn14*Zn42-EE
     &      /SW*Za13*Zn22*Zn34-EE/SW*Za13*Zn24*Zn32+EE/SW
     &      *Za23*Zn22*Zn44+EE/SW*Za23*Zn24*Zn42-Sqrt2*Za13
     &      *Zn42*Zn54*hL-Sqrt2*Za13*Zn44*Zn52*hL-Sqrt2*Za23
     &      *Zn32*Zn54*hL-Sqrt2*Za23*Zn34*Zn52*hL-Sqrt2*Za33
     &      *Zn32*Zn44*hL-Sqrt2*Za33*Zn34*Zn42*hL+2D0*Sqrt2
     &      *Za33*Zn52*Zn54*hK
      AAABC(469) = EE/CW*Za13*Zn12*dconjg(Zn34)+EE/CW*Za13
     &      *Zn14*dconjg(Zn32)-EE/CW*Za23*Zn12*dconjg(Zn44
     &      )-EE/CW*Za23*Zn14*dconjg(Zn42)-EE/SW*Za13*Zn22
     &      *dconjg(Zn34)-EE/SW*Za13*Zn24*dconjg(Zn32)+EE
     &      /SW*Za23*Zn22*dconjg(Zn44)+EE/SW*Za23*Zn24*dconjg(Zn42
     &      )-Sqrt2*Za13*dconjg(Zn42)*dconjg(Zn54)*hL-Sqrt2
     &      *Za13*dconjg(Zn44)*dconjg(Zn52)*hL-Sqrt2*Za23
     &      *dconjg(Zn32)*dconjg(Zn54)*hL-Sqrt2*Za23*dconjg(Zn34
     &      )*dconjg(Zn52)*hL-Sqrt2*Za33*dconjg(Zn32)*dconjg(Zn44
     &      )*hL-Sqrt2*Za33*dconjg(Zn34)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn52)*dconjg(Zn54)*hK
      AAABC(470) = Sqrt2*Za33*Zn52**2*hK+EE/CW*Za13*Zn12
     &      *Zn32-EE/CW*Za23*Zn12*Zn42-EE/SW*Za13*Zn22*Zn32
     &      +EE/SW*Za23*Zn22*Zn42-Sqrt2*Za13*Zn42*Zn52*hL
     &      -Sqrt2*Za23*Zn32*Zn52*hL-Sqrt2*Za33*Zn32*Zn42*hL
      AAABC(471) = Sqrt2*Za33*dconjg(Zn52)**2*hK+EE/CW*Za13
     &      *Zn12*dconjg(Zn32)-EE/CW*Za23*Zn12*dconjg(Zn42
     &      )-EE/SW*Za13*Zn22*dconjg(Zn32)+EE/SW*Za23*Zn22
     &      *dconjg(Zn42)-Sqrt2*Za13*dconjg(Zn42)*dconjg(Zn52
     &      )*hL-Sqrt2*Za23*dconjg(Zn32)*dconjg(Zn52)*hL
     &      -Sqrt2*Za33*dconjg(Zn32)*dconjg(Zn42)*hL
      AAABC(472) = EE/CW*Za11*Zn12*Zn33+EE/CW*Za11*Zn13*Zn32
     &      -EE/CW*Za21*Zn12*Zn43-EE/CW*Za21*Zn13*Zn42-EE
     &      /SW*Za11*Zn22*Zn33-EE/SW*Za11*Zn23*Zn32+EE/SW
     &      *Za21*Zn22*Zn43+EE/SW*Za21*Zn23*Zn42-Sqrt2*Za11
     &      *Zn42*Zn53*hL-Sqrt2*Za11*Zn43*Zn52*hL-Sqrt2*Za21
     &      *Zn32*Zn53*hL-Sqrt2*Za21*Zn33*Zn52*hL-Sqrt2*Za31
     &      *Zn32*Zn43*hL-Sqrt2*Za31*Zn33*Zn42*hL+2D0*Sqrt2
     &      *Za31*Zn52*Zn53*hK
      AAABC(473) = EE/CW*Za11*Zn12*dconjg(Zn33)+EE/CW*Za11
     &      *Zn13*dconjg(Zn32)-EE/CW*Za21*Zn12*dconjg(Zn43
     &      )-EE/CW*Za21*Zn13*dconjg(Zn42)-EE/SW*Za11*Zn22
     &      *dconjg(Zn33)-EE/SW*Za11*Zn23*dconjg(Zn32)+EE
     &      /SW*Za21*Zn22*dconjg(Zn43)+EE/SW*Za21*Zn23*dconjg(Zn42
     &      )-Sqrt2*Za11*dconjg(Zn42)*dconjg(Zn53)*hL-Sqrt2
     &      *Za11*dconjg(Zn43)*dconjg(Zn52)*hL-Sqrt2*Za21
     &      *dconjg(Zn32)*dconjg(Zn53)*hL-Sqrt2*Za21*dconjg(Zn33
     &      )*dconjg(Zn52)*hL-Sqrt2*Za31*dconjg(Zn32)*dconjg(Zn43
     &      )*hL-Sqrt2*Za31*dconjg(Zn33)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn52)*dconjg(Zn53)*hK
      AAABC(474) = EE/CW*Za12*Zn12*Zn34+EE/CW*Za12*Zn14*Zn32
     &      -EE/CW*Za22*Zn12*Zn44-EE/CW*Za22*Zn14*Zn42-EE
     &      /SW*Za12*Zn22*Zn34-EE/SW*Za12*Zn24*Zn32+EE/SW
     &      *Za22*Zn22*Zn44+EE/SW*Za22*Zn24*Zn42-Sqrt2*Za12
     &      *Zn42*Zn54*hL-Sqrt2*Za12*Zn44*Zn52*hL-Sqrt2*Za22
     &      *Zn32*Zn54*hL-Sqrt2*Za22*Zn34*Zn52*hL-Sqrt2*Za32
     &      *Zn32*Zn44*hL-Sqrt2*Za32*Zn34*Zn42*hL+2D0*Sqrt2
     &      *Za32*Zn52*Zn54*hK
      AAABC(475) = EE/CW*Za12*Zn12*dconjg(Zn34)+EE/CW*Za12
     &      *Zn14*dconjg(Zn32)-EE/CW*Za22*Zn12*dconjg(Zn44
     &      )-EE/CW*Za22*Zn14*dconjg(Zn42)-EE/SW*Za12*Zn22
     &      *dconjg(Zn34)-EE/SW*Za12*Zn24*dconjg(Zn32)+EE
     &      /SW*Za22*Zn22*dconjg(Zn44)+EE/SW*Za22*Zn24*dconjg(Zn42
     &      )-Sqrt2*Za12*dconjg(Zn42)*dconjg(Zn54)*hL-Sqrt2
     &      *Za12*dconjg(Zn44)*dconjg(Zn52)*hL-Sqrt2*Za22
     &      *dconjg(Zn32)*dconjg(Zn54)*hL-Sqrt2*Za22*dconjg(Zn34
     &      )*dconjg(Zn52)*hL-Sqrt2*Za32*dconjg(Zn32)*dconjg(Zn44
     &      )*hL-Sqrt2*Za32*dconjg(Zn34)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn52)*dconjg(Zn54)*hK
      AAABC(476) = Sqrt2*Za32*Zn52**2*hK+EE/CW*Za12*Zn12
     &      *Zn32-EE/CW*Za22*Zn12*Zn42-EE/SW*Za12*Zn22*Zn32
     &      +EE/SW*Za22*Zn22*Zn42-Sqrt2*Za12*Zn42*Zn52*hL
     &      -Sqrt2*Za22*Zn32*Zn52*hL-Sqrt2*Za32*Zn32*Zn42*hL
      AAABC(477) = Sqrt2*Za32*dconjg(Zn52)**2*hK+EE/CW*Za12
     &      *Zn12*dconjg(Zn32)-EE/CW*Za22*Zn12*dconjg(Zn42
     &      )-EE/SW*Za12*Zn22*dconjg(Zn32)+EE/SW*Za22*Zn22
     &      *dconjg(Zn42)-Sqrt2*Za12*dconjg(Zn42)*dconjg(Zn52
     &      )*hL-Sqrt2*Za22*dconjg(Zn32)*dconjg(Zn52)*hL
     &      -Sqrt2*Za32*dconjg(Zn32)*dconjg(Zn42)*hL
      AAABC(478) = Sqrt2*Za31*Zn52**2*hK+EE/CW*Za11*Zn12
     &      *Zn32-EE/CW*Za21*Zn12*Zn42-EE/SW*Za11*Zn22*Zn32
     &      +EE/SW*Za21*Zn22*Zn42-Sqrt2*Za11*Zn42*Zn52*hL
     &      -Sqrt2*Za21*Zn32*Zn52*hL-Sqrt2*Za31*Zn32*Zn42*hL
      AAABC(479) = Sqrt2*Za31*dconjg(Zn52)**2*hK+EE/CW*Za11
     &      *Zn12*dconjg(Zn32)-EE/CW*Za21*Zn12*dconjg(Zn42
     &      )-EE/SW*Za11*Zn22*dconjg(Zn32)+EE/SW*Za21*Zn22
     &      *dconjg(Zn42)-Sqrt2*Za11*dconjg(Zn42)*dconjg(Zn52
     &      )*hL-Sqrt2*Za21*dconjg(Zn32)*dconjg(Zn52)*hL
     &      -Sqrt2*Za31*dconjg(Zn32)*dconjg(Zn42)*hL
      AAABC(480) = EE/CW*Za12*Zn12*Zn33+EE/CW*Za12*Zn13*Zn32
     &      -EE/CW*Za22*Zn12*Zn43-EE/CW*Za22*Zn13*Zn42-EE
     &      /SW*Za12*Zn22*Zn33-EE/SW*Za12*Zn23*Zn32+EE/SW
     &      *Za22*Zn22*Zn43+EE/SW*Za22*Zn23*Zn42-Sqrt2*Za12
     &      *Zn42*Zn53*hL-Sqrt2*Za12*Zn43*Zn52*hL-Sqrt2*Za22
     &      *Zn32*Zn53*hL-Sqrt2*Za22*Zn33*Zn52*hL-Sqrt2*Za32
     &      *Zn32*Zn43*hL-Sqrt2*Za32*Zn33*Zn42*hL+2D0*Sqrt2
     &      *Za32*Zn52*Zn53*hK
      AAABC(481) = EE/CW*Za12*Zn12*dconjg(Zn33)+EE/CW*Za12
     &      *Zn13*dconjg(Zn32)-EE/CW*Za22*Zn12*dconjg(Zn43
     &      )-EE/CW*Za22*Zn13*dconjg(Zn42)-EE/SW*Za12*Zn22
     &      *dconjg(Zn33)-EE/SW*Za12*Zn23*dconjg(Zn32)+EE
     &      /SW*Za22*Zn22*dconjg(Zn43)+EE/SW*Za22*Zn23*dconjg(Zn42
     &      )-Sqrt2*Za12*dconjg(Zn42)*dconjg(Zn53)*hL-Sqrt2
     &      *Za12*dconjg(Zn43)*dconjg(Zn52)*hL-Sqrt2*Za22
     &      *dconjg(Zn32)*dconjg(Zn53)*hL-Sqrt2*Za22*dconjg(Zn33
     &      )*dconjg(Zn52)*hL-Sqrt2*Za32*dconjg(Zn32)*dconjg(Zn43
     &      )*hL-Sqrt2*Za32*dconjg(Zn33)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn52)*dconjg(Zn53)*hK
      AAABC(482) = EE/CW*Za11*Zn12*Zn34+EE/CW*Za11*Zn14*Zn32
     &      -EE/CW*Za21*Zn12*Zn44-EE/CW*Za21*Zn14*Zn42-EE
     &      /SW*Za11*Zn22*Zn34-EE/SW*Za11*Zn24*Zn32+EE/SW
     &      *Za21*Zn22*Zn44+EE/SW*Za21*Zn24*Zn42-Sqrt2*Za11
     &      *Zn42*Zn54*hL-Sqrt2*Za11*Zn44*Zn52*hL-Sqrt2*Za21
     &      *Zn32*Zn54*hL-Sqrt2*Za21*Zn34*Zn52*hL-Sqrt2*Za31
     &      *Zn32*Zn44*hL-Sqrt2*Za31*Zn34*Zn42*hL+2D0*Sqrt2
     &      *Za31*Zn52*Zn54*hK
      AAABC(483) = EE/CW*Za11*Zn12*dconjg(Zn34)+EE/CW*Za11
     &      *Zn14*dconjg(Zn32)-EE/CW*Za21*Zn12*dconjg(Zn44
     &      )-EE/CW*Za21*Zn14*dconjg(Zn42)-EE/SW*Za11*Zn22
     &      *dconjg(Zn34)-EE/SW*Za11*Zn24*dconjg(Zn32)+EE
     &      /SW*Za21*Zn22*dconjg(Zn44)+EE/SW*Za21*Zn24*dconjg(Zn42
     &      )-Sqrt2*Za11*dconjg(Zn42)*dconjg(Zn54)*hL-Sqrt2
     &      *Za11*dconjg(Zn44)*dconjg(Zn52)*hL-Sqrt2*Za21
     &      *dconjg(Zn32)*dconjg(Zn54)*hL-Sqrt2*Za21*dconjg(Zn34
     &      )*dconjg(Zn52)*hL-Sqrt2*Za31*dconjg(Zn32)*dconjg(Zn44
     &      )*hL-Sqrt2*Za31*dconjg(Zn34)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn52)*dconjg(Zn54)*hK
      AAABC(484) = EE/CW*Zh13*Zn12*Zn34+EE/CW*Zh13*Zn14*Zn32
     &      -EE/CW*Zh23*Zn12*Zn44-EE/CW*Zh23*Zn14*Zn42-EE
     &      /SW*Zh13*Zn22*Zn34-EE/SW*Zh13*Zn24*Zn32+EE/SW
     &      *Zh23*Zn22*Zn44+EE/SW*Zh23*Zn24*Zn42+Sqrt2*Zh13
     &      *Zn42*Zn54*hL+Sqrt2*Zh13*Zn44*Zn52*hL+Sqrt2*Zh23
     &      *Zn32*Zn54*hL+Sqrt2*Zh23*Zn34*Zn52*hL+Sqrt2*Zh33
     &      *Zn32*Zn44*hL+Sqrt2*Zh33*Zn34*Zn42*hL-2D0*Sqrt2
     &      *Zh33*Zn52*Zn54*hK
      AAABC(485) = EE/CW*Zh13*Zn12*dconjg(Zn34)+EE/CW*Zh13
     &      *Zn14*dconjg(Zn32)-EE/CW*Zh23*Zn12*dconjg(Zn44
     &      )-EE/CW*Zh23*Zn14*dconjg(Zn42)-EE/SW*Zh13*Zn22
     &      *dconjg(Zn34)-EE/SW*Zh13*Zn24*dconjg(Zn32)+EE
     &      /SW*Zh23*Zn22*dconjg(Zn44)+EE/SW*Zh23*Zn24*dconjg(Zn42
     &      )+Sqrt2*Zh13*dconjg(Zn42)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh13*dconjg(Zn44)*dconjg(Zn52)*hL+Sqrt2*Zh23
     &      *dconjg(Zn32)*dconjg(Zn54)*hL+Sqrt2*Zh23*dconjg(Zn34
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh33*dconjg(Zn32)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh33*dconjg(Zn34)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn52)*dconjg(Zn54)*hK
      AAABC(486) = EE/CW*Zh12*Zn12*Zn33+EE/CW*Zh12*Zn13*Zn32
     &      -EE/CW*Zh22*Zn12*Zn43-EE/CW*Zh22*Zn13*Zn42-EE
     &      /SW*Zh12*Zn22*Zn33-EE/SW*Zh12*Zn23*Zn32+EE/SW
     &      *Zh22*Zn22*Zn43+EE/SW*Zh22*Zn23*Zn42+Sqrt2*Zh12
     &      *Zn42*Zn53*hL+Sqrt2*Zh12*Zn43*Zn52*hL+Sqrt2*Zh22
     &      *Zn32*Zn53*hL+Sqrt2*Zh22*Zn33*Zn52*hL+Sqrt2*Zh32
     &      *Zn32*Zn43*hL+Sqrt2*Zh32*Zn33*Zn42*hL-2D0*Sqrt2
     &      *Zh32*Zn52*Zn53*hK
      AAABC(487) = EE/CW*Zh12*Zn12*dconjg(Zn33)+EE/CW*Zh12
     &      *Zn13*dconjg(Zn32)-EE/CW*Zh22*Zn12*dconjg(Zn43
     &      )-EE/CW*Zh22*Zn13*dconjg(Zn42)-EE/SW*Zh12*Zn22
     &      *dconjg(Zn33)-EE/SW*Zh12*Zn23*dconjg(Zn32)+EE
     &      /SW*Zh22*Zn22*dconjg(Zn43)+EE/SW*Zh22*Zn23*dconjg(Zn42
     &      )+Sqrt2*Zh12*dconjg(Zn42)*dconjg(Zn53)*hL+Sqrt2
     &      *Zh12*dconjg(Zn43)*dconjg(Zn52)*hL+Sqrt2*Zh22
     &      *dconjg(Zn32)*dconjg(Zn53)*hL+Sqrt2*Zh22*dconjg(Zn33
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh32*dconjg(Zn32)*dconjg(Zn43
     &      )*hL+Sqrt2*Zh32*dconjg(Zn33)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn52)*dconjg(Zn53)*hK
      AAABC(488) = EE/CW*Zh11*Zn12*Zn34+EE/CW*Zh11*Zn14*Zn32
     &      -EE/CW*Zh21*Zn12*Zn44-EE/CW*Zh21*Zn14*Zn42-EE
     &      /SW*Zh11*Zn22*Zn34-EE/SW*Zh11*Zn24*Zn32+EE/SW
     &      *Zh21*Zn22*Zn44+EE/SW*Zh21*Zn24*Zn42+Sqrt2*Zh11
     &      *Zn42*Zn54*hL+Sqrt2*Zh11*Zn44*Zn52*hL+Sqrt2*Zh21
     &      *Zn32*Zn54*hL+Sqrt2*Zh21*Zn34*Zn52*hL+Sqrt2*Zh31
     &      *Zn32*Zn44*hL+Sqrt2*Zh31*Zn34*Zn42*hL-2D0*Sqrt2
     &      *Zh31*Zn52*Zn54*hK
      AAABC(489) = EE/CW*Zh11*Zn12*dconjg(Zn34)+EE/CW*Zh11
     &      *Zn14*dconjg(Zn32)-EE/CW*Zh21*Zn12*dconjg(Zn44
     &      )-EE/CW*Zh21*Zn14*dconjg(Zn42)-EE/SW*Zh11*Zn22
     &      *dconjg(Zn34)-EE/SW*Zh11*Zn24*dconjg(Zn32)+EE
     &      /SW*Zh21*Zn22*dconjg(Zn44)+EE/SW*Zh21*Zn24*dconjg(Zn42
     &      )+Sqrt2*Zh11*dconjg(Zn42)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh11*dconjg(Zn44)*dconjg(Zn52)*hL+Sqrt2*Zh21
     &      *dconjg(Zn32)*dconjg(Zn54)*hL+Sqrt2*Zh21*dconjg(Zn34
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh31*dconjg(Zn32)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh31*dconjg(Zn34)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn52)*dconjg(Zn54)*hK
      AAABC(490) = EE/CW*Zh12*Zn12*Zn34+EE/CW*Zh12*Zn14*Zn32
     &      -EE/CW*Zh22*Zn12*Zn44-EE/CW*Zh22*Zn14*Zn42-EE
     &      /SW*Zh12*Zn22*Zn34-EE/SW*Zh12*Zn24*Zn32+EE/SW
     &      *Zh22*Zn22*Zn44+EE/SW*Zh22*Zn24*Zn42+Sqrt2*Zh12
     &      *Zn42*Zn54*hL+Sqrt2*Zh12*Zn44*Zn52*hL+Sqrt2*Zh22
     &      *Zn32*Zn54*hL+Sqrt2*Zh22*Zn34*Zn52*hL+Sqrt2*Zh32
     &      *Zn32*Zn44*hL+Sqrt2*Zh32*Zn34*Zn42*hL-2D0*Sqrt2
     &      *Zh32*Zn52*Zn54*hK
      AAABC(491) = EE/CW*Zh12*Zn12*dconjg(Zn34)+EE/CW*Zh12
     &      *Zn14*dconjg(Zn32)-EE/CW*Zh22*Zn12*dconjg(Zn44
     &      )-EE/CW*Zh22*Zn14*dconjg(Zn42)-EE/SW*Zh12*Zn22
     &      *dconjg(Zn34)-EE/SW*Zh12*Zn24*dconjg(Zn32)+EE
     &      /SW*Zh22*Zn22*dconjg(Zn44)+EE/SW*Zh22*Zn24*dconjg(Zn42
     &      )+Sqrt2*Zh12*dconjg(Zn42)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh12*dconjg(Zn44)*dconjg(Zn52)*hL+Sqrt2*Zh22
     &      *dconjg(Zn32)*dconjg(Zn54)*hL+Sqrt2*Zh22*dconjg(Zn34
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh32*dconjg(Zn32)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh32*dconjg(Zn34)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn52)*dconjg(Zn54)*hK
      AAABC(492) = EE/CW*Zh11*Zn12*Zn33+EE/CW*Zh11*Zn13*Zn32
     &      -EE/CW*Zh21*Zn12*Zn43-EE/CW*Zh21*Zn13*Zn42-EE
     &      /SW*Zh11*Zn22*Zn33-EE/SW*Zh11*Zn23*Zn32+EE/SW
     &      *Zh21*Zn22*Zn43+EE/SW*Zh21*Zn23*Zn42+Sqrt2*Zh11
     &      *Zn42*Zn53*hL+Sqrt2*Zh11*Zn43*Zn52*hL+Sqrt2*Zh21
     &      *Zn32*Zn53*hL+Sqrt2*Zh21*Zn33*Zn52*hL+Sqrt2*Zh31
     &      *Zn32*Zn43*hL+Sqrt2*Zh31*Zn33*Zn42*hL-2D0*Sqrt2
     &      *Zh31*Zn52*Zn53*hK
      AAABC(493) = EE/CW*Zh11*Zn12*dconjg(Zn33)+EE/CW*Zh11
     &      *Zn13*dconjg(Zn32)-EE/CW*Zh21*Zn12*dconjg(Zn43
     &      )-EE/CW*Zh21*Zn13*dconjg(Zn42)-EE/SW*Zh11*Zn22
     &      *dconjg(Zn33)-EE/SW*Zh11*Zn23*dconjg(Zn32)+EE
     &      /SW*Zh21*Zn22*dconjg(Zn43)+EE/SW*Zh21*Zn23*dconjg(Zn42
     &      )+Sqrt2*Zh11*dconjg(Zn42)*dconjg(Zn53)*hL+Sqrt2
     &      *Zh11*dconjg(Zn43)*dconjg(Zn52)*hL+Sqrt2*Zh21
     &      *dconjg(Zn32)*dconjg(Zn53)*hL+Sqrt2*Zh21*dconjg(Zn33
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh31*dconjg(Zn32)*dconjg(Zn43
     &      )*hL+Sqrt2*Zh31*dconjg(Zn33)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn52)*dconjg(Zn53)*hK
      AAABC(494) = Sqrt2*Zh33*Zn52**2*hK-EE/CW*Zh13*Zn12
     &      *Zn32+EE/CW*Zh23*Zn12*Zn42+EE/SW*Zh13*Zn22*Zn32
     &      -EE/SW*Zh23*Zn22*Zn42-Sqrt2*Zh13*Zn42*Zn52*hL
     &      -Sqrt2*Zh23*Zn32*Zn52*hL-Sqrt2*Zh33*Zn32*Zn42*hL
      AAABC(495) = Sqrt2*Zh33*dconjg(Zn52)**2*hK-EE/CW*Zh13
     &      *Zn12*dconjg(Zn32)+EE/CW*Zh23*Zn12*dconjg(Zn42
     &      )+EE/SW*Zh13*Zn22*dconjg(Zn32)-EE/SW*Zh23*Zn22
     &      *dconjg(Zn42)-Sqrt2*Zh13*dconjg(Zn42)*dconjg(Zn52
     &      )*hL-Sqrt2*Zh23*dconjg(Zn32)*dconjg(Zn52)*hL
     &      -Sqrt2*Zh33*dconjg(Zn32)*dconjg(Zn42)*hL
      AAABC(496) = Sqrt2*Zh32*Zn52**2*hK-EE/CW*Zh12*Zn12
     &      *Zn32+EE/CW*Zh22*Zn12*Zn42+EE/SW*Zh12*Zn22*Zn32
     &      -EE/SW*Zh22*Zn22*Zn42-Sqrt2*Zh12*Zn42*Zn52*hL
     &      -Sqrt2*Zh22*Zn32*Zn52*hL-Sqrt2*Zh32*Zn32*Zn42*hL
      AAABC(497) = Sqrt2*Zh32*dconjg(Zn52)**2*hK-EE/CW*Zh12
     &      *Zn12*dconjg(Zn32)+EE/CW*Zh22*Zn12*dconjg(Zn42
     &      )+EE/SW*Zh12*Zn22*dconjg(Zn32)-EE/SW*Zh22*Zn22
     &      *dconjg(Zn42)-Sqrt2*Zh12*dconjg(Zn42)*dconjg(Zn52
     &      )*hL-Sqrt2*Zh22*dconjg(Zn32)*dconjg(Zn52)*hL
     &      -Sqrt2*Zh32*dconjg(Zn32)*dconjg(Zn42)*hL
      AAABC(498) = Sqrt2*Zh31*Zn52**2*hK-EE/CW*Zh11*Zn12
     &      *Zn32+EE/CW*Zh21*Zn12*Zn42+EE/SW*Zh11*Zn22*Zn32
     &      -EE/SW*Zh21*Zn22*Zn42-Sqrt2*Zh11*Zn42*Zn52*hL
     &      -Sqrt2*Zh21*Zn32*Zn52*hL-Sqrt2*Zh31*Zn32*Zn42*hL
      AAABC(499) = Sqrt2*Zh31*dconjg(Zn52)**2*hK-EE/CW*Zh11
     &      *Zn12*dconjg(Zn32)+EE/CW*Zh21*Zn12*dconjg(Zn42
     &      )+EE/SW*Zh11*Zn22*dconjg(Zn32)-EE/SW*Zh21*Zn22
     &      *dconjg(Zn42)-Sqrt2*Zh11*dconjg(Zn42)*dconjg(Zn52
     &      )*hL-Sqrt2*Zh21*dconjg(Zn32)*dconjg(Zn52)*hL
     &      -Sqrt2*Zh31*dconjg(Zn32)*dconjg(Zn42)*hL
      AAABC(500) = EE/CW*Zh13*Zn12*Zn33+EE/CW*Zh13*Zn13*Zn32
     &      -EE/CW*Zh23*Zn12*Zn43-EE/CW*Zh23*Zn13*Zn42-EE
     &      /SW*Zh13*Zn22*Zn33-EE/SW*Zh13*Zn23*Zn32+EE/SW
     &      *Zh23*Zn22*Zn43+EE/SW*Zh23*Zn23*Zn42+Sqrt2*Zh13
     &      *Zn42*Zn53*hL+Sqrt2*Zh13*Zn43*Zn52*hL+Sqrt2*Zh23
     &      *Zn32*Zn53*hL+Sqrt2*Zh23*Zn33*Zn52*hL+Sqrt2*Zh33
     &      *Zn32*Zn43*hL+Sqrt2*Zh33*Zn33*Zn42*hL-2D0*Sqrt2
     &      *Zh33*Zn52*Zn53*hK
      AAABC(501) = EE/CW*Zh13*Zn12*dconjg(Zn33)+EE/CW*Zh13
     &      *Zn13*dconjg(Zn32)-EE/CW*Zh23*Zn12*dconjg(Zn43
     &      )-EE/CW*Zh23*Zn13*dconjg(Zn42)-EE/SW*Zh13*Zn22
     &      *dconjg(Zn33)-EE/SW*Zh13*Zn23*dconjg(Zn32)+EE
     &      /SW*Zh23*Zn22*dconjg(Zn43)+EE/SW*Zh23*Zn23*dconjg(Zn42
     &      )+Sqrt2*Zh13*dconjg(Zn42)*dconjg(Zn53)*hL+Sqrt2
     &      *Zh13*dconjg(Zn43)*dconjg(Zn52)*hL+Sqrt2*Zh23
     &      *dconjg(Zn32)*dconjg(Zn53)*hL+Sqrt2*Zh23*dconjg(Zn33
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh33*dconjg(Zn32)*dconjg(Zn43
     &      )*hL+Sqrt2*Zh33*dconjg(Zn33)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn52)*dconjg(Zn53)*hK
      AAABC(502) = EE/CW*Za13*Zn12*Zn35+EE/CW*Za13*Zn15*Zn32
     &      -EE/CW*Za23*Zn12*Zn45-EE/CW*Za23*Zn15*Zn42-EE
     &      /SW*Za13*Zn22*Zn35-EE/SW*Za13*Zn25*Zn32+EE/SW
     &      *Za23*Zn22*Zn45+EE/SW*Za23*Zn25*Zn42-Sqrt2*Za13
     &      *Zn42*Zn55*hL-Sqrt2*Za13*Zn45*Zn52*hL-Sqrt2*Za23
     &      *Zn32*Zn55*hL-Sqrt2*Za23*Zn35*Zn52*hL-Sqrt2*Za33
     &      *Zn32*Zn45*hL-Sqrt2*Za33*Zn35*Zn42*hL+2D0*Sqrt2
     &      *Za33*Zn52*Zn55*hK
      AAABC(503) = EE/CW*Za13*Zn12*dconjg(Zn35)+EE/CW*Za13
     &      *Zn15*dconjg(Zn32)-EE/CW*Za23*Zn12*dconjg(Zn45
     &      )-EE/CW*Za23*Zn15*dconjg(Zn42)-EE/SW*Za13*Zn22
     &      *dconjg(Zn35)-EE/SW*Za13*Zn25*dconjg(Zn32)+EE
     &      /SW*Za23*Zn22*dconjg(Zn45)+EE/SW*Za23*Zn25*dconjg(Zn42
     &      )-Sqrt2*Za13*dconjg(Zn42)*dconjg(Zn55)*hL-Sqrt2
     &      *Za13*dconjg(Zn45)*dconjg(Zn52)*hL-Sqrt2*Za23
     &      *dconjg(Zn32)*dconjg(Zn55)*hL-Sqrt2*Za23*dconjg(Zn35
     &      )*dconjg(Zn52)*hL-Sqrt2*Za33*dconjg(Zn32)*dconjg(Zn45
     &      )*hL-Sqrt2*Za33*dconjg(Zn35)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn52)*dconjg(Zn55)*hK
      AAABC(504) = EE/CW*Za11*Zn12*Zn35+EE/CW*Za11*Zn15*Zn32
     &      -EE/CW*Za21*Zn12*Zn45-EE/CW*Za21*Zn15*Zn42-EE
     &      /SW*Za11*Zn22*Zn35-EE/SW*Za11*Zn25*Zn32+EE/SW
     &      *Za21*Zn22*Zn45+EE/SW*Za21*Zn25*Zn42-Sqrt2*Za11
     &      *Zn42*Zn55*hL-Sqrt2*Za11*Zn45*Zn52*hL-Sqrt2*Za21
     &      *Zn32*Zn55*hL-Sqrt2*Za21*Zn35*Zn52*hL-Sqrt2*Za31
     &      *Zn32*Zn45*hL-Sqrt2*Za31*Zn35*Zn42*hL+2D0*Sqrt2
     &      *Za31*Zn52*Zn55*hK
      AAABC(505) = EE/CW*Za11*Zn12*dconjg(Zn35)+EE/CW*Za11
     &      *Zn15*dconjg(Zn32)-EE/CW*Za21*Zn12*dconjg(Zn45
     &      )-EE/CW*Za21*Zn15*dconjg(Zn42)-EE/SW*Za11*Zn22
     &      *dconjg(Zn35)-EE/SW*Za11*Zn25*dconjg(Zn32)+EE
     &      /SW*Za21*Zn22*dconjg(Zn45)+EE/SW*Za21*Zn25*dconjg(Zn42
     &      )-Sqrt2*Za11*dconjg(Zn42)*dconjg(Zn55)*hL-Sqrt2
     &      *Za11*dconjg(Zn45)*dconjg(Zn52)*hL-Sqrt2*Za21
     &      *dconjg(Zn32)*dconjg(Zn55)*hL-Sqrt2*Za21*dconjg(Zn35
     &      )*dconjg(Zn52)*hL-Sqrt2*Za31*dconjg(Zn32)*dconjg(Zn45
     &      )*hL-Sqrt2*Za31*dconjg(Zn35)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn52)*dconjg(Zn55)*hK
      AAABC(506) = EE/CW*Za12*Zn12*Zn35+EE/CW*Za12*Zn15*Zn32
     &      -EE/CW*Za22*Zn12*Zn45-EE/CW*Za22*Zn15*Zn42-EE
     &      /SW*Za12*Zn22*Zn35-EE/SW*Za12*Zn25*Zn32+EE/SW
     &      *Za22*Zn22*Zn45+EE/SW*Za22*Zn25*Zn42-Sqrt2*Za12
     &      *Zn42*Zn55*hL-Sqrt2*Za12*Zn45*Zn52*hL-Sqrt2*Za22
     &      *Zn32*Zn55*hL-Sqrt2*Za22*Zn35*Zn52*hL-Sqrt2*Za32
     &      *Zn32*Zn45*hL-Sqrt2*Za32*Zn35*Zn42*hL+2D0*Sqrt2
     &      *Za32*Zn52*Zn55*hK
      AAABC(507) = EE/CW*Za12*Zn12*dconjg(Zn35)+EE/CW*Za12
     &      *Zn15*dconjg(Zn32)-EE/CW*Za22*Zn12*dconjg(Zn45
     &      )-EE/CW*Za22*Zn15*dconjg(Zn42)-EE/SW*Za12*Zn22
     &      *dconjg(Zn35)-EE/SW*Za12*Zn25*dconjg(Zn32)+EE
     &      /SW*Za22*Zn22*dconjg(Zn45)+EE/SW*Za22*Zn25*dconjg(Zn42
     &      )-Sqrt2*Za12*dconjg(Zn42)*dconjg(Zn55)*hL-Sqrt2
     &      *Za12*dconjg(Zn45)*dconjg(Zn52)*hL-Sqrt2*Za22
     &      *dconjg(Zn32)*dconjg(Zn55)*hL-Sqrt2*Za22*dconjg(Zn35
     &      )*dconjg(Zn52)*hL-Sqrt2*Za32*dconjg(Zn32)*dconjg(Zn45
     &      )*hL-Sqrt2*Za32*dconjg(Zn35)*dconjg(Zn42)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn52)*dconjg(Zn55)*hK
      AAABC(508) = EE/CW*Zh13*Zn12*Zn35+EE/CW*Zh13*Zn15*Zn32
     &      -EE/CW*Zh23*Zn12*Zn45-EE/CW*Zh23*Zn15*Zn42-EE
     &      /SW*Zh13*Zn22*Zn35-EE/SW*Zh13*Zn25*Zn32+EE/SW
     &      *Zh23*Zn22*Zn45+EE/SW*Zh23*Zn25*Zn42+Sqrt2*Zh13
     &      *Zn42*Zn55*hL+Sqrt2*Zh13*Zn45*Zn52*hL+Sqrt2*Zh23
     &      *Zn32*Zn55*hL+Sqrt2*Zh23*Zn35*Zn52*hL+Sqrt2*Zh33
     &      *Zn32*Zn45*hL+Sqrt2*Zh33*Zn35*Zn42*hL-2D0*Sqrt2
     &      *Zh33*Zn52*Zn55*hK
      end

      subroutine aaini02
      implicit none
#include "model.h"

      AAABC(509) = EE/CW*Zh13*Zn12*dconjg(Zn35)+EE/CW*Zh13
     &      *Zn15*dconjg(Zn32)-EE/CW*Zh23*Zn12*dconjg(Zn45
     &      )-EE/CW*Zh23*Zn15*dconjg(Zn42)-EE/SW*Zh13*Zn22
     &      *dconjg(Zn35)-EE/SW*Zh13*Zn25*dconjg(Zn32)+EE
     &      /SW*Zh23*Zn22*dconjg(Zn45)+EE/SW*Zh23*Zn25*dconjg(Zn42
     &      )+Sqrt2*Zh13*dconjg(Zn42)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh13*dconjg(Zn45)*dconjg(Zn52)*hL+Sqrt2*Zh23
     &      *dconjg(Zn32)*dconjg(Zn55)*hL+Sqrt2*Zh23*dconjg(Zn35
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh33*dconjg(Zn32)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh33*dconjg(Zn35)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn52)*dconjg(Zn55)*hK
      AAABC(510) = EE/CW*Zh12*Zn12*Zn35+EE/CW*Zh12*Zn15*Zn32
     &      -EE/CW*Zh22*Zn12*Zn45-EE/CW*Zh22*Zn15*Zn42-EE
     &      /SW*Zh12*Zn22*Zn35-EE/SW*Zh12*Zn25*Zn32+EE/SW
     &      *Zh22*Zn22*Zn45+EE/SW*Zh22*Zn25*Zn42+Sqrt2*Zh12
     &      *Zn42*Zn55*hL+Sqrt2*Zh12*Zn45*Zn52*hL+Sqrt2*Zh22
     &      *Zn32*Zn55*hL+Sqrt2*Zh22*Zn35*Zn52*hL+Sqrt2*Zh32
     &      *Zn32*Zn45*hL+Sqrt2*Zh32*Zn35*Zn42*hL-2D0*Sqrt2
     &      *Zh32*Zn52*Zn55*hK
      AAABC(511) = EE/CW*Zh12*Zn12*dconjg(Zn35)+EE/CW*Zh12
     &      *Zn15*dconjg(Zn32)-EE/CW*Zh22*Zn12*dconjg(Zn45
     &      )-EE/CW*Zh22*Zn15*dconjg(Zn42)-EE/SW*Zh12*Zn22
     &      *dconjg(Zn35)-EE/SW*Zh12*Zn25*dconjg(Zn32)+EE
     &      /SW*Zh22*Zn22*dconjg(Zn45)+EE/SW*Zh22*Zn25*dconjg(Zn42
     &      )+Sqrt2*Zh12*dconjg(Zn42)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh12*dconjg(Zn45)*dconjg(Zn52)*hL+Sqrt2*Zh22
     &      *dconjg(Zn32)*dconjg(Zn55)*hL+Sqrt2*Zh22*dconjg(Zn35
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh32*dconjg(Zn32)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh32*dconjg(Zn35)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn52)*dconjg(Zn55)*hK
      AAABC(512) = EE/CW*Zh11*Zn12*Zn35+EE/CW*Zh11*Zn15*Zn32
     &      -EE/CW*Zh21*Zn12*Zn45-EE/CW*Zh21*Zn15*Zn42-EE
     &      /SW*Zh11*Zn22*Zn35-EE/SW*Zh11*Zn25*Zn32+EE/SW
     &      *Zh21*Zn22*Zn45+EE/SW*Zh21*Zn25*Zn42+Sqrt2*Zh11
     &      *Zn42*Zn55*hL+Sqrt2*Zh11*Zn45*Zn52*hL+Sqrt2*Zh21
     &      *Zn32*Zn55*hL+Sqrt2*Zh21*Zn35*Zn52*hL+Sqrt2*Zh31
     &      *Zn32*Zn45*hL+Sqrt2*Zh31*Zn35*Zn42*hL-2D0*Sqrt2
     &      *Zh31*Zn52*Zn55*hK
      AAABC(513) = EE/CW*Zh11*Zn12*dconjg(Zn35)+EE/CW*Zh11
     &      *Zn15*dconjg(Zn32)-EE/CW*Zh21*Zn12*dconjg(Zn45
     &      )-EE/CW*Zh21*Zn15*dconjg(Zn42)-EE/SW*Zh11*Zn22
     &      *dconjg(Zn35)-EE/SW*Zh11*Zn25*dconjg(Zn32)+EE
     &      /SW*Zh21*Zn22*dconjg(Zn45)+EE/SW*Zh21*Zn25*dconjg(Zn42
     &      )+Sqrt2*Zh11*dconjg(Zn42)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh11*dconjg(Zn45)*dconjg(Zn52)*hL+Sqrt2*Zh21
     &      *dconjg(Zn32)*dconjg(Zn55)*hL+Sqrt2*Zh21*dconjg(Zn35
     &      )*dconjg(Zn52)*hL+Sqrt2*Zh31*dconjg(Zn32)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh31*dconjg(Zn35)*dconjg(Zn42)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn52)*dconjg(Zn55)*hK
      AAABC(514) = Sqrt2*Za33*Zn53**2*hK+EE/CW*Za13*Zn13
     &      *Zn33-EE/CW*Za23*Zn13*Zn43-EE/SW*Za13*Zn23*Zn33
     &      +EE/SW*Za23*Zn23*Zn43-Sqrt2*Za13*Zn43*Zn53*hL
     &      -Sqrt2*Za23*Zn33*Zn53*hL-Sqrt2*Za33*Zn33*Zn43*hL
      AAABC(515) = Sqrt2*Za33*dconjg(Zn53)**2*hK+EE/CW*Za13
     &      *Zn13*dconjg(Zn33)-EE/CW*Za23*Zn13*dconjg(Zn43
     &      )-EE/SW*Za13*Zn23*dconjg(Zn33)+EE/SW*Za23*Zn23
     &      *dconjg(Zn43)-Sqrt2*Za13*dconjg(Zn43)*dconjg(Zn53
     &      )*hL-Sqrt2*Za23*dconjg(Zn33)*dconjg(Zn53)*hL
     &      -Sqrt2*Za33*dconjg(Zn33)*dconjg(Zn43)*hL
      AAABC(516) = EE/CW*Za13*Zn13*Zn34+EE/CW*Za13*Zn14*Zn33
     &      -EE/CW*Za23*Zn13*Zn44-EE/CW*Za23*Zn14*Zn43-EE
     &      /SW*Za13*Zn23*Zn34-EE/SW*Za13*Zn24*Zn33+EE/SW
     &      *Za23*Zn23*Zn44+EE/SW*Za23*Zn24*Zn43-Sqrt2*Za13
     &      *Zn43*Zn54*hL-Sqrt2*Za13*Zn44*Zn53*hL-Sqrt2*Za23
     &      *Zn33*Zn54*hL-Sqrt2*Za23*Zn34*Zn53*hL-Sqrt2*Za33
     &      *Zn33*Zn44*hL-Sqrt2*Za33*Zn34*Zn43*hL+2D0*Sqrt2
     &      *Za33*Zn53*Zn54*hK
      AAABC(517) = EE/CW*Za13*Zn13*dconjg(Zn34)+EE/CW*Za13
     &      *Zn14*dconjg(Zn33)-EE/CW*Za23*Zn13*dconjg(Zn44
     &      )-EE/CW*Za23*Zn14*dconjg(Zn43)-EE/SW*Za13*Zn23
     &      *dconjg(Zn34)-EE/SW*Za13*Zn24*dconjg(Zn33)+EE
     &      /SW*Za23*Zn23*dconjg(Zn44)+EE/SW*Za23*Zn24*dconjg(Zn43
     &      )-Sqrt2*Za13*dconjg(Zn43)*dconjg(Zn54)*hL-Sqrt2
     &      *Za13*dconjg(Zn44)*dconjg(Zn53)*hL-Sqrt2*Za23
     &      *dconjg(Zn33)*dconjg(Zn54)*hL-Sqrt2*Za23*dconjg(Zn34
     &      )*dconjg(Zn53)*hL-Sqrt2*Za33*dconjg(Zn33)*dconjg(Zn44
     &      )*hL-Sqrt2*Za33*dconjg(Zn34)*dconjg(Zn43)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn53)*dconjg(Zn54)*hK
      AAABC(518) = Sqrt2*Za31*Zn53**2*hK+EE/CW*Za11*Zn13
     &      *Zn33-EE/CW*Za21*Zn13*Zn43-EE/SW*Za11*Zn23*Zn33
     &      +EE/SW*Za21*Zn23*Zn43-Sqrt2*Za11*Zn43*Zn53*hL
     &      -Sqrt2*Za21*Zn33*Zn53*hL-Sqrt2*Za31*Zn33*Zn43*hL
      AAABC(519) = Sqrt2*Za31*dconjg(Zn53)**2*hK+EE/CW*Za11
     &      *Zn13*dconjg(Zn33)-EE/CW*Za21*Zn13*dconjg(Zn43
     &      )-EE/SW*Za11*Zn23*dconjg(Zn33)+EE/SW*Za21*Zn23
     &      *dconjg(Zn43)-Sqrt2*Za11*dconjg(Zn43)*dconjg(Zn53
     &      )*hL-Sqrt2*Za21*dconjg(Zn33)*dconjg(Zn53)*hL
     &      -Sqrt2*Za31*dconjg(Zn33)*dconjg(Zn43)*hL
      AAABC(520) = EE/CW*Za11*Zn13*Zn34+EE/CW*Za11*Zn14*Zn33
     &      -EE/CW*Za21*Zn13*Zn44-EE/CW*Za21*Zn14*Zn43-EE
     &      /SW*Za11*Zn23*Zn34-EE/SW*Za11*Zn24*Zn33+EE/SW
     &      *Za21*Zn23*Zn44+EE/SW*Za21*Zn24*Zn43-Sqrt2*Za11
     &      *Zn43*Zn54*hL-Sqrt2*Za11*Zn44*Zn53*hL-Sqrt2*Za21
     &      *Zn33*Zn54*hL-Sqrt2*Za21*Zn34*Zn53*hL-Sqrt2*Za31
     &      *Zn33*Zn44*hL-Sqrt2*Za31*Zn34*Zn43*hL+2D0*Sqrt2
     &      *Za31*Zn53*Zn54*hK
      AAABC(521) = EE/CW*Za11*Zn13*dconjg(Zn34)+EE/CW*Za11
     &      *Zn14*dconjg(Zn33)-EE/CW*Za21*Zn13*dconjg(Zn44
     &      )-EE/CW*Za21*Zn14*dconjg(Zn43)-EE/SW*Za11*Zn23
     &      *dconjg(Zn34)-EE/SW*Za11*Zn24*dconjg(Zn33)+EE
     &      /SW*Za21*Zn23*dconjg(Zn44)+EE/SW*Za21*Zn24*dconjg(Zn43
     &      )-Sqrt2*Za11*dconjg(Zn43)*dconjg(Zn54)*hL-Sqrt2
     &      *Za11*dconjg(Zn44)*dconjg(Zn53)*hL-Sqrt2*Za21
     &      *dconjg(Zn33)*dconjg(Zn54)*hL-Sqrt2*Za21*dconjg(Zn34
     &      )*dconjg(Zn53)*hL-Sqrt2*Za31*dconjg(Zn33)*dconjg(Zn44
     &      )*hL-Sqrt2*Za31*dconjg(Zn34)*dconjg(Zn43)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn53)*dconjg(Zn54)*hK
      AAABC(522) = Sqrt2*Za32*Zn53**2*hK+EE/CW*Za12*Zn13
     &      *Zn33-EE/CW*Za22*Zn13*Zn43-EE/SW*Za12*Zn23*Zn33
     &      +EE/SW*Za22*Zn23*Zn43-Sqrt2*Za12*Zn43*Zn53*hL
     &      -Sqrt2*Za22*Zn33*Zn53*hL-Sqrt2*Za32*Zn33*Zn43*hL
      AAABC(523) = Sqrt2*Za32*dconjg(Zn53)**2*hK+EE/CW*Za12
     &      *Zn13*dconjg(Zn33)-EE/CW*Za22*Zn13*dconjg(Zn43
     &      )-EE/SW*Za12*Zn23*dconjg(Zn33)+EE/SW*Za22*Zn23
     &      *dconjg(Zn43)-Sqrt2*Za12*dconjg(Zn43)*dconjg(Zn53
     &      )*hL-Sqrt2*Za22*dconjg(Zn33)*dconjg(Zn53)*hL
     &      -Sqrt2*Za32*dconjg(Zn33)*dconjg(Zn43)*hL
      AAABC(524) = EE/CW*Za12*Zn13*Zn34+EE/CW*Za12*Zn14*Zn33
     &      -EE/CW*Za22*Zn13*Zn44-EE/CW*Za22*Zn14*Zn43-EE
     &      /SW*Za12*Zn23*Zn34-EE/SW*Za12*Zn24*Zn33+EE/SW
     &      *Za22*Zn23*Zn44+EE/SW*Za22*Zn24*Zn43-Sqrt2*Za12
     &      *Zn43*Zn54*hL-Sqrt2*Za12*Zn44*Zn53*hL-Sqrt2*Za22
     &      *Zn33*Zn54*hL-Sqrt2*Za22*Zn34*Zn53*hL-Sqrt2*Za32
     &      *Zn33*Zn44*hL-Sqrt2*Za32*Zn34*Zn43*hL+2D0*Sqrt2
     &      *Za32*Zn53*Zn54*hK
      AAABC(525) = EE/CW*Za12*Zn13*dconjg(Zn34)+EE/CW*Za12
     &      *Zn14*dconjg(Zn33)-EE/CW*Za22*Zn13*dconjg(Zn44
     &      )-EE/CW*Za22*Zn14*dconjg(Zn43)-EE/SW*Za12*Zn23
     &      *dconjg(Zn34)-EE/SW*Za12*Zn24*dconjg(Zn33)+EE
     &      /SW*Za22*Zn23*dconjg(Zn44)+EE/SW*Za22*Zn24*dconjg(Zn43
     &      )-Sqrt2*Za12*dconjg(Zn43)*dconjg(Zn54)*hL-Sqrt2
     &      *Za12*dconjg(Zn44)*dconjg(Zn53)*hL-Sqrt2*Za22
     &      *dconjg(Zn33)*dconjg(Zn54)*hL-Sqrt2*Za22*dconjg(Zn34
     &      )*dconjg(Zn53)*hL-Sqrt2*Za32*dconjg(Zn33)*dconjg(Zn44
     &      )*hL-Sqrt2*Za32*dconjg(Zn34)*dconjg(Zn43)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn53)*dconjg(Zn54)*hK
      AAABC(526) = Sqrt2*Zh32*Zn53**2*hK-EE/CW*Zh12*Zn13
     &      *Zn33+EE/CW*Zh22*Zn13*Zn43+EE/SW*Zh12*Zn23*Zn33
     &      -EE/SW*Zh22*Zn23*Zn43-Sqrt2*Zh12*Zn43*Zn53*hL
     &      -Sqrt2*Zh22*Zn33*Zn53*hL-Sqrt2*Zh32*Zn33*Zn43*hL
      AAABC(527) = Sqrt2*Zh32*dconjg(Zn53)**2*hK-EE/CW*Zh12
     &      *Zn13*dconjg(Zn33)+EE/CW*Zh22*Zn13*dconjg(Zn43
     &      )+EE/SW*Zh12*Zn23*dconjg(Zn33)-EE/SW*Zh22*Zn23
     &      *dconjg(Zn43)-Sqrt2*Zh12*dconjg(Zn43)*dconjg(Zn53
     &      )*hL-Sqrt2*Zh22*dconjg(Zn33)*dconjg(Zn53)*hL
     &      -Sqrt2*Zh32*dconjg(Zn33)*dconjg(Zn43)*hL
      AAABC(528) = Sqrt2*Zh31*Zn53**2*hK-EE/CW*Zh11*Zn13
     &      *Zn33+EE/CW*Zh21*Zn13*Zn43+EE/SW*Zh11*Zn23*Zn33
     &      -EE/SW*Zh21*Zn23*Zn43-Sqrt2*Zh11*Zn43*Zn53*hL
     &      -Sqrt2*Zh21*Zn33*Zn53*hL-Sqrt2*Zh31*Zn33*Zn43*hL
      AAABC(529) = Sqrt2*Zh31*dconjg(Zn53)**2*hK-EE/CW*Zh11
     &      *Zn13*dconjg(Zn33)+EE/CW*Zh21*Zn13*dconjg(Zn43
     &      )+EE/SW*Zh11*Zn23*dconjg(Zn33)-EE/SW*Zh21*Zn23
     &      *dconjg(Zn43)-Sqrt2*Zh11*dconjg(Zn43)*dconjg(Zn53
     &      )*hL-Sqrt2*Zh21*dconjg(Zn33)*dconjg(Zn53)*hL
     &      -Sqrt2*Zh31*dconjg(Zn33)*dconjg(Zn43)*hL
      AAABC(530) = EE/CW*Zh13*Zn13*Zn34+EE/CW*Zh13*Zn14*Zn33
     &      -EE/CW*Zh23*Zn13*Zn44-EE/CW*Zh23*Zn14*Zn43-EE
     &      /SW*Zh13*Zn23*Zn34-EE/SW*Zh13*Zn24*Zn33+EE/SW
     &      *Zh23*Zn23*Zn44+EE/SW*Zh23*Zn24*Zn43+Sqrt2*Zh13
     &      *Zn43*Zn54*hL+Sqrt2*Zh13*Zn44*Zn53*hL+Sqrt2*Zh23
     &      *Zn33*Zn54*hL+Sqrt2*Zh23*Zn34*Zn53*hL+Sqrt2*Zh33
     &      *Zn33*Zn44*hL+Sqrt2*Zh33*Zn34*Zn43*hL-2D0*Sqrt2
     &      *Zh33*Zn53*Zn54*hK
      AAABC(531) = EE/CW*Zh13*Zn13*dconjg(Zn34)+EE/CW*Zh13
     &      *Zn14*dconjg(Zn33)-EE/CW*Zh23*Zn13*dconjg(Zn44
     &      )-EE/CW*Zh23*Zn14*dconjg(Zn43)-EE/SW*Zh13*Zn23
     &      *dconjg(Zn34)-EE/SW*Zh13*Zn24*dconjg(Zn33)+EE
     &      /SW*Zh23*Zn23*dconjg(Zn44)+EE/SW*Zh23*Zn24*dconjg(Zn43
     &      )+Sqrt2*Zh13*dconjg(Zn43)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh13*dconjg(Zn44)*dconjg(Zn53)*hL+Sqrt2*Zh23
     &      *dconjg(Zn33)*dconjg(Zn54)*hL+Sqrt2*Zh23*dconjg(Zn34
     &      )*dconjg(Zn53)*hL+Sqrt2*Zh33*dconjg(Zn33)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh33*dconjg(Zn34)*dconjg(Zn43)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn53)*dconjg(Zn54)*hK
      AAABC(532) = EE/CW*Zh11*Zn13*Zn34+EE/CW*Zh11*Zn14*Zn33
     &      -EE/CW*Zh21*Zn13*Zn44-EE/CW*Zh21*Zn14*Zn43-EE
     &      /SW*Zh11*Zn23*Zn34-EE/SW*Zh11*Zn24*Zn33+EE/SW
     &      *Zh21*Zn23*Zn44+EE/SW*Zh21*Zn24*Zn43+Sqrt2*Zh11
     &      *Zn43*Zn54*hL+Sqrt2*Zh11*Zn44*Zn53*hL+Sqrt2*Zh21
     &      *Zn33*Zn54*hL+Sqrt2*Zh21*Zn34*Zn53*hL+Sqrt2*Zh31
     &      *Zn33*Zn44*hL+Sqrt2*Zh31*Zn34*Zn43*hL-2D0*Sqrt2
     &      *Zh31*Zn53*Zn54*hK
      AAABC(533) = EE/CW*Zh11*Zn13*dconjg(Zn34)+EE/CW*Zh11
     &      *Zn14*dconjg(Zn33)-EE/CW*Zh21*Zn13*dconjg(Zn44
     &      )-EE/CW*Zh21*Zn14*dconjg(Zn43)-EE/SW*Zh11*Zn23
     &      *dconjg(Zn34)-EE/SW*Zh11*Zn24*dconjg(Zn33)+EE
     &      /SW*Zh21*Zn23*dconjg(Zn44)+EE/SW*Zh21*Zn24*dconjg(Zn43
     &      )+Sqrt2*Zh11*dconjg(Zn43)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh11*dconjg(Zn44)*dconjg(Zn53)*hL+Sqrt2*Zh21
     &      *dconjg(Zn33)*dconjg(Zn54)*hL+Sqrt2*Zh21*dconjg(Zn34
     &      )*dconjg(Zn53)*hL+Sqrt2*Zh31*dconjg(Zn33)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh31*dconjg(Zn34)*dconjg(Zn43)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn53)*dconjg(Zn54)*hK
      AAABC(534) = Sqrt2*Zh33*Zn53**2*hK-EE/CW*Zh13*Zn13
     &      *Zn33+EE/CW*Zh23*Zn13*Zn43+EE/SW*Zh13*Zn23*Zn33
     &      -EE/SW*Zh23*Zn23*Zn43-Sqrt2*Zh13*Zn43*Zn53*hL
     &      -Sqrt2*Zh23*Zn33*Zn53*hL-Sqrt2*Zh33*Zn33*Zn43*hL
      AAABC(535) = Sqrt2*Zh33*dconjg(Zn53)**2*hK-EE/CW*Zh13
     &      *Zn13*dconjg(Zn33)+EE/CW*Zh23*Zn13*dconjg(Zn43
     &      )+EE/SW*Zh13*Zn23*dconjg(Zn33)-EE/SW*Zh23*Zn23
     &      *dconjg(Zn43)-Sqrt2*Zh13*dconjg(Zn43)*dconjg(Zn53
     &      )*hL-Sqrt2*Zh23*dconjg(Zn33)*dconjg(Zn53)*hL
     &      -Sqrt2*Zh33*dconjg(Zn33)*dconjg(Zn43)*hL
      AAABC(536) = EE/CW*Zh12*Zn13*Zn34+EE/CW*Zh12*Zn14*Zn33
     &      -EE/CW*Zh22*Zn13*Zn44-EE/CW*Zh22*Zn14*Zn43-EE
     &      /SW*Zh12*Zn23*Zn34-EE/SW*Zh12*Zn24*Zn33+EE/SW
     &      *Zh22*Zn23*Zn44+EE/SW*Zh22*Zn24*Zn43+Sqrt2*Zh12
     &      *Zn43*Zn54*hL+Sqrt2*Zh12*Zn44*Zn53*hL+Sqrt2*Zh22
     &      *Zn33*Zn54*hL+Sqrt2*Zh22*Zn34*Zn53*hL+Sqrt2*Zh32
     &      *Zn33*Zn44*hL+Sqrt2*Zh32*Zn34*Zn43*hL-2D0*Sqrt2
     &      *Zh32*Zn53*Zn54*hK
      AAABC(537) = EE/CW*Zh12*Zn13*dconjg(Zn34)+EE/CW*Zh12
     &      *Zn14*dconjg(Zn33)-EE/CW*Zh22*Zn13*dconjg(Zn44
     &      )-EE/CW*Zh22*Zn14*dconjg(Zn43)-EE/SW*Zh12*Zn23
     &      *dconjg(Zn34)-EE/SW*Zh12*Zn24*dconjg(Zn33)+EE
     &      /SW*Zh22*Zn23*dconjg(Zn44)+EE/SW*Zh22*Zn24*dconjg(Zn43
     &      )+Sqrt2*Zh12*dconjg(Zn43)*dconjg(Zn54)*hL+Sqrt2
     &      *Zh12*dconjg(Zn44)*dconjg(Zn53)*hL+Sqrt2*Zh22
     &      *dconjg(Zn33)*dconjg(Zn54)*hL+Sqrt2*Zh22*dconjg(Zn34
     &      )*dconjg(Zn53)*hL+Sqrt2*Zh32*dconjg(Zn33)*dconjg(Zn44
     &      )*hL+Sqrt2*Zh32*dconjg(Zn34)*dconjg(Zn43)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn53)*dconjg(Zn54)*hK
      AAABC(538) = EE/CW*Za13*Zn13*Zn35+EE/CW*Za13*Zn15*Zn33
     &      -EE/CW*Za23*Zn13*Zn45-EE/CW*Za23*Zn15*Zn43-EE
     &      /SW*Za13*Zn23*Zn35-EE/SW*Za13*Zn25*Zn33+EE/SW
     &      *Za23*Zn23*Zn45+EE/SW*Za23*Zn25*Zn43-Sqrt2*Za13
     &      *Zn43*Zn55*hL-Sqrt2*Za13*Zn45*Zn53*hL-Sqrt2*Za23
     &      *Zn33*Zn55*hL-Sqrt2*Za23*Zn35*Zn53*hL-Sqrt2*Za33
     &      *Zn33*Zn45*hL-Sqrt2*Za33*Zn35*Zn43*hL+2D0*Sqrt2
     &      *Za33*Zn53*Zn55*hK
      AAABC(539) = EE/CW*Za13*Zn13*dconjg(Zn35)+EE/CW*Za13
     &      *Zn15*dconjg(Zn33)-EE/CW*Za23*Zn13*dconjg(Zn45
     &      )-EE/CW*Za23*Zn15*dconjg(Zn43)-EE/SW*Za13*Zn23
     &      *dconjg(Zn35)-EE/SW*Za13*Zn25*dconjg(Zn33)+EE
     &      /SW*Za23*Zn23*dconjg(Zn45)+EE/SW*Za23*Zn25*dconjg(Zn43
     &      )-Sqrt2*Za13*dconjg(Zn43)*dconjg(Zn55)*hL-Sqrt2
     &      *Za13*dconjg(Zn45)*dconjg(Zn53)*hL-Sqrt2*Za23
     &      *dconjg(Zn33)*dconjg(Zn55)*hL-Sqrt2*Za23*dconjg(Zn35
     &      )*dconjg(Zn53)*hL-Sqrt2*Za33*dconjg(Zn33)*dconjg(Zn45
     &      )*hL-Sqrt2*Za33*dconjg(Zn35)*dconjg(Zn43)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn53)*dconjg(Zn55)*hK
      AAABC(540) = EE/CW*Za11*Zn13*Zn35+EE/CW*Za11*Zn15*Zn33
     &      -EE/CW*Za21*Zn13*Zn45-EE/CW*Za21*Zn15*Zn43-EE
     &      /SW*Za11*Zn23*Zn35-EE/SW*Za11*Zn25*Zn33+EE/SW
     &      *Za21*Zn23*Zn45+EE/SW*Za21*Zn25*Zn43-Sqrt2*Za11
     &      *Zn43*Zn55*hL-Sqrt2*Za11*Zn45*Zn53*hL-Sqrt2*Za21
     &      *Zn33*Zn55*hL-Sqrt2*Za21*Zn35*Zn53*hL-Sqrt2*Za31
     &      *Zn33*Zn45*hL-Sqrt2*Za31*Zn35*Zn43*hL+2D0*Sqrt2
     &      *Za31*Zn53*Zn55*hK
      AAABC(541) = EE/CW*Za11*Zn13*dconjg(Zn35)+EE/CW*Za11
     &      *Zn15*dconjg(Zn33)-EE/CW*Za21*Zn13*dconjg(Zn45
     &      )-EE/CW*Za21*Zn15*dconjg(Zn43)-EE/SW*Za11*Zn23
     &      *dconjg(Zn35)-EE/SW*Za11*Zn25*dconjg(Zn33)+EE
     &      /SW*Za21*Zn23*dconjg(Zn45)+EE/SW*Za21*Zn25*dconjg(Zn43
     &      )-Sqrt2*Za11*dconjg(Zn43)*dconjg(Zn55)*hL-Sqrt2
     &      *Za11*dconjg(Zn45)*dconjg(Zn53)*hL-Sqrt2*Za21
     &      *dconjg(Zn33)*dconjg(Zn55)*hL-Sqrt2*Za21*dconjg(Zn35
     &      )*dconjg(Zn53)*hL-Sqrt2*Za31*dconjg(Zn33)*dconjg(Zn45
     &      )*hL-Sqrt2*Za31*dconjg(Zn35)*dconjg(Zn43)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn53)*dconjg(Zn55)*hK
      AAABC(542) = EE/CW*Za12*Zn13*Zn35+EE/CW*Za12*Zn15*Zn33
     &      -EE/CW*Za22*Zn13*Zn45-EE/CW*Za22*Zn15*Zn43-EE
     &      /SW*Za12*Zn23*Zn35-EE/SW*Za12*Zn25*Zn33+EE/SW
     &      *Za22*Zn23*Zn45+EE/SW*Za22*Zn25*Zn43-Sqrt2*Za12
     &      *Zn43*Zn55*hL-Sqrt2*Za12*Zn45*Zn53*hL-Sqrt2*Za22
     &      *Zn33*Zn55*hL-Sqrt2*Za22*Zn35*Zn53*hL-Sqrt2*Za32
     &      *Zn33*Zn45*hL-Sqrt2*Za32*Zn35*Zn43*hL+2D0*Sqrt2
     &      *Za32*Zn53*Zn55*hK
      AAABC(543) = EE/CW*Za12*Zn13*dconjg(Zn35)+EE/CW*Za12
     &      *Zn15*dconjg(Zn33)-EE/CW*Za22*Zn13*dconjg(Zn45
     &      )-EE/CW*Za22*Zn15*dconjg(Zn43)-EE/SW*Za12*Zn23
     &      *dconjg(Zn35)-EE/SW*Za12*Zn25*dconjg(Zn33)+EE
     &      /SW*Za22*Zn23*dconjg(Zn45)+EE/SW*Za22*Zn25*dconjg(Zn43
     &      )-Sqrt2*Za12*dconjg(Zn43)*dconjg(Zn55)*hL-Sqrt2
     &      *Za12*dconjg(Zn45)*dconjg(Zn53)*hL-Sqrt2*Za22
     &      *dconjg(Zn33)*dconjg(Zn55)*hL-Sqrt2*Za22*dconjg(Zn35
     &      )*dconjg(Zn53)*hL-Sqrt2*Za32*dconjg(Zn33)*dconjg(Zn45
     &      )*hL-Sqrt2*Za32*dconjg(Zn35)*dconjg(Zn43)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn53)*dconjg(Zn55)*hK
      AAABC(544) = EE/CW*Zh12*Zn13*Zn35+EE/CW*Zh12*Zn15*Zn33
     &      -EE/CW*Zh22*Zn13*Zn45-EE/CW*Zh22*Zn15*Zn43-EE
     &      /SW*Zh12*Zn23*Zn35-EE/SW*Zh12*Zn25*Zn33+EE/SW
     &      *Zh22*Zn23*Zn45+EE/SW*Zh22*Zn25*Zn43+Sqrt2*Zh12
     &      *Zn43*Zn55*hL+Sqrt2*Zh12*Zn45*Zn53*hL+Sqrt2*Zh22
     &      *Zn33*Zn55*hL+Sqrt2*Zh22*Zn35*Zn53*hL+Sqrt2*Zh32
     &      *Zn33*Zn45*hL+Sqrt2*Zh32*Zn35*Zn43*hL-2D0*Sqrt2
     &      *Zh32*Zn53*Zn55*hK
      AAABC(545) = EE/CW*Zh12*Zn13*dconjg(Zn35)+EE/CW*Zh12
     &      *Zn15*dconjg(Zn33)-EE/CW*Zh22*Zn13*dconjg(Zn45
     &      )-EE/CW*Zh22*Zn15*dconjg(Zn43)-EE/SW*Zh12*Zn23
     &      *dconjg(Zn35)-EE/SW*Zh12*Zn25*dconjg(Zn33)+EE
     &      /SW*Zh22*Zn23*dconjg(Zn45)+EE/SW*Zh22*Zn25*dconjg(Zn43
     &      )+Sqrt2*Zh12*dconjg(Zn43)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh12*dconjg(Zn45)*dconjg(Zn53)*hL+Sqrt2*Zh22
     &      *dconjg(Zn33)*dconjg(Zn55)*hL+Sqrt2*Zh22*dconjg(Zn35
     &      )*dconjg(Zn53)*hL+Sqrt2*Zh32*dconjg(Zn33)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh32*dconjg(Zn35)*dconjg(Zn43)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn53)*dconjg(Zn55)*hK
      AAABC(546) = EE/CW*Zh11*Zn13*Zn35+EE/CW*Zh11*Zn15*Zn33
     &      -EE/CW*Zh21*Zn13*Zn45-EE/CW*Zh21*Zn15*Zn43-EE
     &      /SW*Zh11*Zn23*Zn35-EE/SW*Zh11*Zn25*Zn33+EE/SW
     &      *Zh21*Zn23*Zn45+EE/SW*Zh21*Zn25*Zn43+Sqrt2*Zh11
     &      *Zn43*Zn55*hL+Sqrt2*Zh11*Zn45*Zn53*hL+Sqrt2*Zh21
     &      *Zn33*Zn55*hL+Sqrt2*Zh21*Zn35*Zn53*hL+Sqrt2*Zh31
     &      *Zn33*Zn45*hL+Sqrt2*Zh31*Zn35*Zn43*hL-2D0*Sqrt2
     &      *Zh31*Zn53*Zn55*hK
      AAABC(547) = EE/CW*Zh11*Zn13*dconjg(Zn35)+EE/CW*Zh11
     &      *Zn15*dconjg(Zn33)-EE/CW*Zh21*Zn13*dconjg(Zn45
     &      )-EE/CW*Zh21*Zn15*dconjg(Zn43)-EE/SW*Zh11*Zn23
     &      *dconjg(Zn35)-EE/SW*Zh11*Zn25*dconjg(Zn33)+EE
     &      /SW*Zh21*Zn23*dconjg(Zn45)+EE/SW*Zh21*Zn25*dconjg(Zn43
     &      )+Sqrt2*Zh11*dconjg(Zn43)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh11*dconjg(Zn45)*dconjg(Zn53)*hL+Sqrt2*Zh21
     &      *dconjg(Zn33)*dconjg(Zn55)*hL+Sqrt2*Zh21*dconjg(Zn35
     &      )*dconjg(Zn53)*hL+Sqrt2*Zh31*dconjg(Zn33)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh31*dconjg(Zn35)*dconjg(Zn43)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn53)*dconjg(Zn55)*hK
      AAABC(548) = EE/CW*Zh13*Zn13*Zn35+EE/CW*Zh13*Zn15*Zn33
     &      -EE/CW*Zh23*Zn13*Zn45-EE/CW*Zh23*Zn15*Zn43-EE
     &      /SW*Zh13*Zn23*Zn35-EE/SW*Zh13*Zn25*Zn33+EE/SW
     &      *Zh23*Zn23*Zn45+EE/SW*Zh23*Zn25*Zn43+Sqrt2*Zh13
     &      *Zn43*Zn55*hL+Sqrt2*Zh13*Zn45*Zn53*hL+Sqrt2*Zh23
     &      *Zn33*Zn55*hL+Sqrt2*Zh23*Zn35*Zn53*hL+Sqrt2*Zh33
     &      *Zn33*Zn45*hL+Sqrt2*Zh33*Zn35*Zn43*hL-2D0*Sqrt2
     &      *Zh33*Zn53*Zn55*hK
      AAABC(549) = EE/CW*Zh13*Zn13*dconjg(Zn35)+EE/CW*Zh13
     &      *Zn15*dconjg(Zn33)-EE/CW*Zh23*Zn13*dconjg(Zn45
     &      )-EE/CW*Zh23*Zn15*dconjg(Zn43)-EE/SW*Zh13*Zn23
     &      *dconjg(Zn35)-EE/SW*Zh13*Zn25*dconjg(Zn33)+EE
     &      /SW*Zh23*Zn23*dconjg(Zn45)+EE/SW*Zh23*Zn25*dconjg(Zn43
     &      )+Sqrt2*Zh13*dconjg(Zn43)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh13*dconjg(Zn45)*dconjg(Zn53)*hL+Sqrt2*Zh23
     &      *dconjg(Zn33)*dconjg(Zn55)*hL+Sqrt2*Zh23*dconjg(Zn35
     &      )*dconjg(Zn53)*hL+Sqrt2*Zh33*dconjg(Zn33)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh33*dconjg(Zn35)*dconjg(Zn43)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn53)*dconjg(Zn55)*hK
      AAABC(550) = Sqrt2*Za33*Zn54**2*hK+EE/CW*Za13*Zn14
     &      *Zn34-EE/CW*Za23*Zn14*Zn44-EE/SW*Za13*Zn24*Zn34
     &      +EE/SW*Za23*Zn24*Zn44-Sqrt2*Za13*Zn44*Zn54*hL
     &      -Sqrt2*Za23*Zn34*Zn54*hL-Sqrt2*Za33*Zn34*Zn44*hL
      AAABC(551) = Sqrt2*Za33*dconjg(Zn54)**2*hK+EE/CW*Za13
     &      *Zn14*dconjg(Zn34)-EE/CW*Za23*Zn14*dconjg(Zn44
     &      )-EE/SW*Za13*Zn24*dconjg(Zn34)+EE/SW*Za23*Zn24
     &      *dconjg(Zn44)-Sqrt2*Za13*dconjg(Zn44)*dconjg(Zn54
     &      )*hL-Sqrt2*Za23*dconjg(Zn34)*dconjg(Zn54)*hL
     &      -Sqrt2*Za33*dconjg(Zn34)*dconjg(Zn44)*hL
      AAABC(552) = Sqrt2*Za32*Zn54**2*hK+EE/CW*Za12*Zn14
     &      *Zn34-EE/CW*Za22*Zn14*Zn44-EE/SW*Za12*Zn24*Zn34
     &      +EE/SW*Za22*Zn24*Zn44-Sqrt2*Za12*Zn44*Zn54*hL
     &      -Sqrt2*Za22*Zn34*Zn54*hL-Sqrt2*Za32*Zn34*Zn44*hL
      AAABC(553) = Sqrt2*Za32*dconjg(Zn54)**2*hK+EE/CW*Za12
     &      *Zn14*dconjg(Zn34)-EE/CW*Za22*Zn14*dconjg(Zn44
     &      )-EE/SW*Za12*Zn24*dconjg(Zn34)+EE/SW*Za22*Zn24
     &      *dconjg(Zn44)-Sqrt2*Za12*dconjg(Zn44)*dconjg(Zn54
     &      )*hL-Sqrt2*Za22*dconjg(Zn34)*dconjg(Zn54)*hL
     &      -Sqrt2*Za32*dconjg(Zn34)*dconjg(Zn44)*hL
      AAABC(554) = Sqrt2*Za31*Zn54**2*hK+EE/CW*Za11*Zn14
     &      *Zn34-EE/CW*Za21*Zn14*Zn44-EE/SW*Za11*Zn24*Zn34
     &      +EE/SW*Za21*Zn24*Zn44-Sqrt2*Za11*Zn44*Zn54*hL
     &      -Sqrt2*Za21*Zn34*Zn54*hL-Sqrt2*Za31*Zn34*Zn44*hL
      AAABC(555) = Sqrt2*Za31*dconjg(Zn54)**2*hK+EE/CW*Za11
     &      *Zn14*dconjg(Zn34)-EE/CW*Za21*Zn14*dconjg(Zn44
     &      )-EE/SW*Za11*Zn24*dconjg(Zn34)+EE/SW*Za21*Zn24
     &      *dconjg(Zn44)-Sqrt2*Za11*dconjg(Zn44)*dconjg(Zn54
     &      )*hL-Sqrt2*Za21*dconjg(Zn34)*dconjg(Zn54)*hL
     &      -Sqrt2*Za31*dconjg(Zn34)*dconjg(Zn44)*hL
      AAABC(556) = Sqrt2*Zh32*Zn54**2*hK-EE/CW*Zh12*Zn14
     &      *Zn34+EE/CW*Zh22*Zn14*Zn44+EE/SW*Zh12*Zn24*Zn34
     &      -EE/SW*Zh22*Zn24*Zn44-Sqrt2*Zh12*Zn44*Zn54*hL
     &      -Sqrt2*Zh22*Zn34*Zn54*hL-Sqrt2*Zh32*Zn34*Zn44*hL
      AAABC(557) = Sqrt2*Zh32*dconjg(Zn54)**2*hK-EE/CW*Zh12
     &      *Zn14*dconjg(Zn34)+EE/CW*Zh22*Zn14*dconjg(Zn44
     &      )+EE/SW*Zh12*Zn24*dconjg(Zn34)-EE/SW*Zh22*Zn24
     &      *dconjg(Zn44)-Sqrt2*Zh12*dconjg(Zn44)*dconjg(Zn54
     &      )*hL-Sqrt2*Zh22*dconjg(Zn34)*dconjg(Zn54)*hL
     &      -Sqrt2*Zh32*dconjg(Zn34)*dconjg(Zn44)*hL
      AAABC(558) = Sqrt2*Zh31*Zn54**2*hK-EE/CW*Zh11*Zn14
     &      *Zn34+EE/CW*Zh21*Zn14*Zn44+EE/SW*Zh11*Zn24*Zn34
     &      -EE/SW*Zh21*Zn24*Zn44-Sqrt2*Zh11*Zn44*Zn54*hL
     &      -Sqrt2*Zh21*Zn34*Zn54*hL-Sqrt2*Zh31*Zn34*Zn44*hL
      AAABC(559) = Sqrt2*Zh31*dconjg(Zn54)**2*hK-EE/CW*Zh11
     &      *Zn14*dconjg(Zn34)+EE/CW*Zh21*Zn14*dconjg(Zn44
     &      )+EE/SW*Zh11*Zn24*dconjg(Zn34)-EE/SW*Zh21*Zn24
     &      *dconjg(Zn44)-Sqrt2*Zh11*dconjg(Zn44)*dconjg(Zn54
     &      )*hL-Sqrt2*Zh21*dconjg(Zn34)*dconjg(Zn54)*hL
     &      -Sqrt2*Zh31*dconjg(Zn34)*dconjg(Zn44)*hL
      AAABC(560) = Sqrt2*Zh33*Zn54**2*hK-EE/CW*Zh13*Zn14
     &      *Zn34+EE/CW*Zh23*Zn14*Zn44+EE/SW*Zh13*Zn24*Zn34
     &      -EE/SW*Zh23*Zn24*Zn44-Sqrt2*Zh13*Zn44*Zn54*hL
     &      -Sqrt2*Zh23*Zn34*Zn54*hL-Sqrt2*Zh33*Zn34*Zn44*hL
      AAABC(561) = Sqrt2*Zh33*dconjg(Zn54)**2*hK-EE/CW*Zh13
     &      *Zn14*dconjg(Zn34)+EE/CW*Zh23*Zn14*dconjg(Zn44
     &      )+EE/SW*Zh13*Zn24*dconjg(Zn34)-EE/SW*Zh23*Zn24
     &      *dconjg(Zn44)-Sqrt2*Zh13*dconjg(Zn44)*dconjg(Zn54
     &      )*hL-Sqrt2*Zh23*dconjg(Zn34)*dconjg(Zn54)*hL
     &      -Sqrt2*Zh33*dconjg(Zn34)*dconjg(Zn44)*hL
      AAABC(562) = EE/CW*Za13*Zn14*Zn35+EE/CW*Za13*Zn15*Zn34
     &      -EE/CW*Za23*Zn14*Zn45-EE/CW*Za23*Zn15*Zn44-EE
     &      /SW*Za13*Zn24*Zn35-EE/SW*Za13*Zn25*Zn34+EE/SW
     &      *Za23*Zn24*Zn45+EE/SW*Za23*Zn25*Zn44-Sqrt2*Za13
     &      *Zn44*Zn55*hL-Sqrt2*Za13*Zn45*Zn54*hL-Sqrt2*Za23
     &      *Zn34*Zn55*hL-Sqrt2*Za23*Zn35*Zn54*hL-Sqrt2*Za33
     &      *Zn34*Zn45*hL-Sqrt2*Za33*Zn35*Zn44*hL+2D0*Sqrt2
     &      *Za33*Zn54*Zn55*hK
      AAABC(563) = EE/CW*Za13*Zn14*dconjg(Zn35)+EE/CW*Za13
     &      *Zn15*dconjg(Zn34)-EE/CW*Za23*Zn14*dconjg(Zn45
     &      )-EE/CW*Za23*Zn15*dconjg(Zn44)-EE/SW*Za13*Zn24
     &      *dconjg(Zn35)-EE/SW*Za13*Zn25*dconjg(Zn34)+EE
     &      /SW*Za23*Zn24*dconjg(Zn45)+EE/SW*Za23*Zn25*dconjg(Zn44
     &      )-Sqrt2*Za13*dconjg(Zn44)*dconjg(Zn55)*hL-Sqrt2
     &      *Za13*dconjg(Zn45)*dconjg(Zn54)*hL-Sqrt2*Za23
     &      *dconjg(Zn34)*dconjg(Zn55)*hL-Sqrt2*Za23*dconjg(Zn35
     &      )*dconjg(Zn54)*hL-Sqrt2*Za33*dconjg(Zn34)*dconjg(Zn45
     &      )*hL-Sqrt2*Za33*dconjg(Zn35)*dconjg(Zn44)*hL
     &      +2D0*Sqrt2*Za33*dconjg(Zn54)*dconjg(Zn55)*hK
      AAABC(564) = EE/CW*Za12*Zn14*Zn35+EE/CW*Za12*Zn15*Zn34
     &      -EE/CW*Za22*Zn14*Zn45-EE/CW*Za22*Zn15*Zn44-EE
     &      /SW*Za12*Zn24*Zn35-EE/SW*Za12*Zn25*Zn34+EE/SW
     &      *Za22*Zn24*Zn45+EE/SW*Za22*Zn25*Zn44-Sqrt2*Za12
     &      *Zn44*Zn55*hL-Sqrt2*Za12*Zn45*Zn54*hL-Sqrt2*Za22
     &      *Zn34*Zn55*hL-Sqrt2*Za22*Zn35*Zn54*hL-Sqrt2*Za32
     &      *Zn34*Zn45*hL-Sqrt2*Za32*Zn35*Zn44*hL+2D0*Sqrt2
     &      *Za32*Zn54*Zn55*hK
      AAABC(565) = EE/CW*Za12*Zn14*dconjg(Zn35)+EE/CW*Za12
     &      *Zn15*dconjg(Zn34)-EE/CW*Za22*Zn14*dconjg(Zn45
     &      )-EE/CW*Za22*Zn15*dconjg(Zn44)-EE/SW*Za12*Zn24
     &      *dconjg(Zn35)-EE/SW*Za12*Zn25*dconjg(Zn34)+EE
     &      /SW*Za22*Zn24*dconjg(Zn45)+EE/SW*Za22*Zn25*dconjg(Zn44
     &      )-Sqrt2*Za12*dconjg(Zn44)*dconjg(Zn55)*hL-Sqrt2
     &      *Za12*dconjg(Zn45)*dconjg(Zn54)*hL-Sqrt2*Za22
     &      *dconjg(Zn34)*dconjg(Zn55)*hL-Sqrt2*Za22*dconjg(Zn35
     &      )*dconjg(Zn54)*hL-Sqrt2*Za32*dconjg(Zn34)*dconjg(Zn45
     &      )*hL-Sqrt2*Za32*dconjg(Zn35)*dconjg(Zn44)*hL
     &      +2D0*Sqrt2*Za32*dconjg(Zn54)*dconjg(Zn55)*hK
      AAABC(566) = EE/CW*Za11*Zn14*Zn35+EE/CW*Za11*Zn15*Zn34
     &      -EE/CW*Za21*Zn14*Zn45-EE/CW*Za21*Zn15*Zn44-EE
     &      /SW*Za11*Zn24*Zn35-EE/SW*Za11*Zn25*Zn34+EE/SW
     &      *Za21*Zn24*Zn45+EE/SW*Za21*Zn25*Zn44-Sqrt2*Za11
     &      *Zn44*Zn55*hL-Sqrt2*Za11*Zn45*Zn54*hL-Sqrt2*Za21
     &      *Zn34*Zn55*hL-Sqrt2*Za21*Zn35*Zn54*hL-Sqrt2*Za31
     &      *Zn34*Zn45*hL-Sqrt2*Za31*Zn35*Zn44*hL+2D0*Sqrt2
     &      *Za31*Zn54*Zn55*hK
      AAABC(567) = EE/CW*Za11*Zn14*dconjg(Zn35)+EE/CW*Za11
     &      *Zn15*dconjg(Zn34)-EE/CW*Za21*Zn14*dconjg(Zn45
     &      )-EE/CW*Za21*Zn15*dconjg(Zn44)-EE/SW*Za11*Zn24
     &      *dconjg(Zn35)-EE/SW*Za11*Zn25*dconjg(Zn34)+EE
     &      /SW*Za21*Zn24*dconjg(Zn45)+EE/SW*Za21*Zn25*dconjg(Zn44
     &      )-Sqrt2*Za11*dconjg(Zn44)*dconjg(Zn55)*hL-Sqrt2
     &      *Za11*dconjg(Zn45)*dconjg(Zn54)*hL-Sqrt2*Za21
     &      *dconjg(Zn34)*dconjg(Zn55)*hL-Sqrt2*Za21*dconjg(Zn35
     &      )*dconjg(Zn54)*hL-Sqrt2*Za31*dconjg(Zn34)*dconjg(Zn45
     &      )*hL-Sqrt2*Za31*dconjg(Zn35)*dconjg(Zn44)*hL
     &      +2D0*Sqrt2*Za31*dconjg(Zn54)*dconjg(Zn55)*hK
      AAABC(568) = EE/CW*Zh13*Zn14*Zn35+EE/CW*Zh13*Zn15*Zn34
     &      -EE/CW*Zh23*Zn14*Zn45-EE/CW*Zh23*Zn15*Zn44-EE
     &      /SW*Zh13*Zn24*Zn35-EE/SW*Zh13*Zn25*Zn34+EE/SW
     &      *Zh23*Zn24*Zn45+EE/SW*Zh23*Zn25*Zn44+Sqrt2*Zh13
     &      *Zn44*Zn55*hL+Sqrt2*Zh13*Zn45*Zn54*hL+Sqrt2*Zh23
     &      *Zn34*Zn55*hL+Sqrt2*Zh23*Zn35*Zn54*hL+Sqrt2*Zh33
     &      *Zn34*Zn45*hL+Sqrt2*Zh33*Zn35*Zn44*hL-2D0*Sqrt2
     &      *Zh33*Zn54*Zn55*hK
      AAABC(569) = EE/CW*Zh13*Zn14*dconjg(Zn35)+EE/CW*Zh13
     &      *Zn15*dconjg(Zn34)-EE/CW*Zh23*Zn14*dconjg(Zn45
     &      )-EE/CW*Zh23*Zn15*dconjg(Zn44)-EE/SW*Zh13*Zn24
     &      *dconjg(Zn35)-EE/SW*Zh13*Zn25*dconjg(Zn34)+EE
     &      /SW*Zh23*Zn24*dconjg(Zn45)+EE/SW*Zh23*Zn25*dconjg(Zn44
     &      )+Sqrt2*Zh13*dconjg(Zn44)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh13*dconjg(Zn45)*dconjg(Zn54)*hL+Sqrt2*Zh23
     &      *dconjg(Zn34)*dconjg(Zn55)*hL+Sqrt2*Zh23*dconjg(Zn35
     &      )*dconjg(Zn54)*hL+Sqrt2*Zh33*dconjg(Zn34)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh33*dconjg(Zn35)*dconjg(Zn44)*hL
     &      -2D0*Sqrt2*Zh33*dconjg(Zn54)*dconjg(Zn55)*hK
      AAABC(570) = EE/CW*Zh12*Zn14*Zn35+EE/CW*Zh12*Zn15*Zn34
     &      -EE/CW*Zh22*Zn14*Zn45-EE/CW*Zh22*Zn15*Zn44-EE
     &      /SW*Zh12*Zn24*Zn35-EE/SW*Zh12*Zn25*Zn34+EE/SW
     &      *Zh22*Zn24*Zn45+EE/SW*Zh22*Zn25*Zn44+Sqrt2*Zh12
     &      *Zn44*Zn55*hL+Sqrt2*Zh12*Zn45*Zn54*hL+Sqrt2*Zh22
     &      *Zn34*Zn55*hL+Sqrt2*Zh22*Zn35*Zn54*hL+Sqrt2*Zh32
     &      *Zn34*Zn45*hL+Sqrt2*Zh32*Zn35*Zn44*hL-2D0*Sqrt2
     &      *Zh32*Zn54*Zn55*hK
      AAABC(571) = EE/CW*Zh12*Zn14*dconjg(Zn35)+EE/CW*Zh12
     &      *Zn15*dconjg(Zn34)-EE/CW*Zh22*Zn14*dconjg(Zn45
     &      )-EE/CW*Zh22*Zn15*dconjg(Zn44)-EE/SW*Zh12*Zn24
     &      *dconjg(Zn35)-EE/SW*Zh12*Zn25*dconjg(Zn34)+EE
     &      /SW*Zh22*Zn24*dconjg(Zn45)+EE/SW*Zh22*Zn25*dconjg(Zn44
     &      )+Sqrt2*Zh12*dconjg(Zn44)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh12*dconjg(Zn45)*dconjg(Zn54)*hL+Sqrt2*Zh22
     &      *dconjg(Zn34)*dconjg(Zn55)*hL+Sqrt2*Zh22*dconjg(Zn35
     &      )*dconjg(Zn54)*hL+Sqrt2*Zh32*dconjg(Zn34)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh32*dconjg(Zn35)*dconjg(Zn44)*hL
     &      -2D0*Sqrt2*Zh32*dconjg(Zn54)*dconjg(Zn55)*hK
      AAABC(572) = EE/CW*Zh11*Zn14*Zn35+EE/CW*Zh11*Zn15*Zn34
     &      -EE/CW*Zh21*Zn14*Zn45-EE/CW*Zh21*Zn15*Zn44-EE
     &      /SW*Zh11*Zn24*Zn35-EE/SW*Zh11*Zn25*Zn34+EE/SW
     &      *Zh21*Zn24*Zn45+EE/SW*Zh21*Zn25*Zn44+Sqrt2*Zh11
     &      *Zn44*Zn55*hL+Sqrt2*Zh11*Zn45*Zn54*hL+Sqrt2*Zh21
     &      *Zn34*Zn55*hL+Sqrt2*Zh21*Zn35*Zn54*hL+Sqrt2*Zh31
     &      *Zn34*Zn45*hL+Sqrt2*Zh31*Zn35*Zn44*hL-2D0*Sqrt2
     &      *Zh31*Zn54*Zn55*hK
      AAABC(573) = EE/CW*Zh11*Zn14*dconjg(Zn35)+EE/CW*Zh11
     &      *Zn15*dconjg(Zn34)-EE/CW*Zh21*Zn14*dconjg(Zn45
     &      )-EE/CW*Zh21*Zn15*dconjg(Zn44)-EE/SW*Zh11*Zn24
     &      *dconjg(Zn35)-EE/SW*Zh11*Zn25*dconjg(Zn34)+EE
     &      /SW*Zh21*Zn24*dconjg(Zn45)+EE/SW*Zh21*Zn25*dconjg(Zn44
     &      )+Sqrt2*Zh11*dconjg(Zn44)*dconjg(Zn55)*hL+Sqrt2
     &      *Zh11*dconjg(Zn45)*dconjg(Zn54)*hL+Sqrt2*Zh21
     &      *dconjg(Zn34)*dconjg(Zn55)*hL+Sqrt2*Zh21*dconjg(Zn35
     &      )*dconjg(Zn54)*hL+Sqrt2*Zh31*dconjg(Zn34)*dconjg(Zn45
     &      )*hL+Sqrt2*Zh31*dconjg(Zn35)*dconjg(Zn44)*hL
     &      -2D0*Sqrt2*Zh31*dconjg(Zn54)*dconjg(Zn55)*hK
      AAABC(574) = Sqrt2*Za33*Zn55**2*hK+EE/CW*Za13*Zn15
     &      *Zn35-EE/CW*Za23*Zn15*Zn45-EE/SW*Za13*Zn25*Zn35
     &      +EE/SW*Za23*Zn25*Zn45-Sqrt2*Za13*Zn45*Zn55*hL
     &      -Sqrt2*Za23*Zn35*Zn55*hL-Sqrt2*Za33*Zn35*Zn45*hL
      AAABC(575) = Sqrt2*Za33*dconjg(Zn55)**2*hK+EE/CW*Za13
     &      *Zn15*dconjg(Zn35)-EE/CW*Za23*Zn15*dconjg(Zn45
     &      )-EE/SW*Za13*Zn25*dconjg(Zn35)+EE/SW*Za23*Zn25
     &      *dconjg(Zn45)-Sqrt2*Za13*dconjg(Zn45)*dconjg(Zn55
     &      )*hL-Sqrt2*Za23*dconjg(Zn35)*dconjg(Zn55)*hL
     &      -Sqrt2*Za33*dconjg(Zn35)*dconjg(Zn45)*hL
      AAABC(576) = Sqrt2*Za31*Zn55**2*hK+EE/CW*Za11*Zn15
     &      *Zn35-EE/CW*Za21*Zn15*Zn45-EE/SW*Za11*Zn25*Zn35
     &      +EE/SW*Za21*Zn25*Zn45-Sqrt2*Za11*Zn45*Zn55*hL
     &      -Sqrt2*Za21*Zn35*Zn55*hL-Sqrt2*Za31*Zn35*Zn45*hL
      AAABC(577) = Sqrt2*Za31*dconjg(Zn55)**2*hK+EE/CW*Za11
     &      *Zn15*dconjg(Zn35)-EE/CW*Za21*Zn15*dconjg(Zn45
     &      )-EE/SW*Za11*Zn25*dconjg(Zn35)+EE/SW*Za21*Zn25
     &      *dconjg(Zn45)-Sqrt2*Za11*dconjg(Zn45)*dconjg(Zn55
     &      )*hL-Sqrt2*Za21*dconjg(Zn35)*dconjg(Zn55)*hL
     &      -Sqrt2*Za31*dconjg(Zn35)*dconjg(Zn45)*hL
      AAABC(578) = Sqrt2*Za32*Zn55**2*hK+EE/CW*Za12*Zn15
     &      *Zn35-EE/CW*Za22*Zn15*Zn45-EE/SW*Za12*Zn25*Zn35
     &      +EE/SW*Za22*Zn25*Zn45-Sqrt2*Za12*Zn45*Zn55*hL
     &      -Sqrt2*Za22*Zn35*Zn55*hL-Sqrt2*Za32*Zn35*Zn45*hL
      AAABC(579) = Sqrt2*Za32*dconjg(Zn55)**2*hK+EE/CW*Za12
     &      *Zn15*dconjg(Zn35)-EE/CW*Za22*Zn15*dconjg(Zn45
     &      )-EE/SW*Za12*Zn25*dconjg(Zn35)+EE/SW*Za22*Zn25
     &      *dconjg(Zn45)-Sqrt2*Za12*dconjg(Zn45)*dconjg(Zn55
     &      )*hL-Sqrt2*Za22*dconjg(Zn35)*dconjg(Zn55)*hL
     &      -Sqrt2*Za32*dconjg(Zn35)*dconjg(Zn45)*hL
      AAABC(580) = Sqrt2*Zh33*Zn55**2*hK-EE/CW*Zh13*Zn15
     &      *Zn35+EE/CW*Zh23*Zn15*Zn45+EE/SW*Zh13*Zn25*Zn35
     &      -EE/SW*Zh23*Zn25*Zn45-Sqrt2*Zh13*Zn45*Zn55*hL
     &      -Sqrt2*Zh23*Zn35*Zn55*hL-Sqrt2*Zh33*Zn35*Zn45*hL
      AAABC(581) = Sqrt2*Zh33*dconjg(Zn55)**2*hK-EE/CW*Zh13
     &      *Zn15*dconjg(Zn35)+EE/CW*Zh23*Zn15*dconjg(Zn45
     &      )+EE/SW*Zh13*Zn25*dconjg(Zn35)-EE/SW*Zh23*Zn25
     &      *dconjg(Zn45)-Sqrt2*Zh13*dconjg(Zn45)*dconjg(Zn55
     &      )*hL-Sqrt2*Zh23*dconjg(Zn35)*dconjg(Zn55)*hL
     &      -Sqrt2*Zh33*dconjg(Zn35)*dconjg(Zn45)*hL
      AAABC(582) = Sqrt2*Zh31*Zn55**2*hK-EE/CW*Zh11*Zn15
     &      *Zn35+EE/CW*Zh21*Zn15*Zn45+EE/SW*Zh11*Zn25*Zn35
     &      -EE/SW*Zh21*Zn25*Zn45-Sqrt2*Zh11*Zn45*Zn55*hL
     &      -Sqrt2*Zh21*Zn35*Zn55*hL-Sqrt2*Zh31*Zn35*Zn45*hL
      AAABC(583) = Sqrt2*Zh31*dconjg(Zn55)**2*hK-EE/CW*Zh11
     &      *Zn15*dconjg(Zn35)+EE/CW*Zh21*Zn15*dconjg(Zn45
     &      )+EE/SW*Zh11*Zn25*dconjg(Zn35)-EE/SW*Zh21*Zn25
     &      *dconjg(Zn45)-Sqrt2*Zh11*dconjg(Zn45)*dconjg(Zn55
     &      )*hL-Sqrt2*Zh21*dconjg(Zn35)*dconjg(Zn55)*hL
     &      -Sqrt2*Zh31*dconjg(Zn35)*dconjg(Zn45)*hL
      AAABC(584) = Sqrt2*Zh32*Zn55**2*hK-EE/CW*Zh12*Zn15
     &      *Zn35+EE/CW*Zh22*Zn15*Zn45+EE/SW*Zh12*Zn25*Zn35
     &      -EE/SW*Zh22*Zn25*Zn45-Sqrt2*Zh12*Zn45*Zn55*hL
     &      -Sqrt2*Zh22*Zn35*Zn55*hL-Sqrt2*Zh32*Zn35*Zn45*hL
      AAABC(585) = Sqrt2*Zh32*dconjg(Zn55)**2*hK-EE/CW*Zh12
     &      *Zn15*dconjg(Zn35)+EE/CW*Zh22*Zn15*dconjg(Zn45
     &      )+EE/SW*Zh12*Zn25*dconjg(Zn35)-EE/SW*Zh22*Zn25
     &      *dconjg(Zn45)-Sqrt2*Zh12*dconjg(Zn45)*dconjg(Zn55
     &      )*hL-Sqrt2*Zh22*dconjg(Zn35)*dconjg(Zn55)*hL
     &      -Sqrt2*Zh32*dconjg(Zn35)*dconjg(Zn45)*hL
      AAABR(493) = EE*(SW/CW*Zm21*Zm22-2D0*CW/SW*Zm11*Zm12
     &      -CW/SW*Zm21*Zm22)
      AAABR(494) = EE*(SW/CW*Zp21*Zp22-2D0*CW/SW*Zp11*Zp12
     &      -CW/SW*Zp21*Zp22)
      AAABR(495) = EE*(SW/CW*Zm21**2-2D0*CW/SW*Zm11**2-CW
     &      /SW*Zm21**2)
      AAABR(496) = EE*(SW/CW*Zp21**2-2D0*CW/SW*Zp11**2-CW
     &      /SW*Zp21**2)
      AAABR(497) = EE*(SW/CW*Zm22**2-2D0*CW/SW*Zm12**2-CW
     &      /SW*Zm22**2)
      AAABR(498) = EE*(SW/CW*Zp22**2-2D0*CW/SW*Zp12**2-CW
     &      /SW*Zp22**2)
      AAABC(586) = EE/SW*(2D0*Zm12*Zn24+Sqrt2*Zm22*dconjg(Zn34
     &      ))
      AAABC(587) = EE/SW*(2D0*Zn24*Zp12-Sqrt2*Zn44*Zp22)
      AAABC(588) = EE/SW*(2D0*Zm12*Zn23+Sqrt2*Zm22*dconjg(Zn33
     &      ))
      AAABC(589) = EE/SW*(2D0*Zn23*Zp12-Sqrt2*Zn43*Zp22)
      AAABC(590) = EE/SW*(2D0*Zm12*Zn21+Sqrt2*Zm22*dconjg(Zn31
     &      ))
      AAABC(591) = EE/SW*(2D0*Zn21*Zp12-Sqrt2*Zn41*Zp22)
      AAABC(592) = EE/SW*(2D0*Zm11*Zn22+Sqrt2*Zm21*dconjg(Zn32
     &      ))
      AAABC(593) = EE/SW*(2D0*Zn22*Zp11-Sqrt2*Zn42*Zp21)
      AAABC(594) = EE/SW*(2D0*Zm12*Zn22+Sqrt2*Zm22*dconjg(Zn32
     &      ))
      AAABC(595) = EE/SW*(2D0*Zn22*Zp12-Sqrt2*Zn42*Zp22)
      AAABC(596) = EE/SW*(2D0*Zm11*Zn24+Sqrt2*Zm21*dconjg(Zn34
     &      ))
      AAABC(597) = EE/SW*(2D0*Zn24*Zp11-Sqrt2*Zn44*Zp21)
      AAABC(598) = EE/SW*(2D0*Zm11*Zn23+Sqrt2*Zm21*dconjg(Zn33
     &      ))
      AAABC(599) = EE/SW*(2D0*Zn23*Zp11-Sqrt2*Zn43*Zp21)
      AAABC(600) = EE/SW*(2D0*Zm11*Zn21+Sqrt2*Zm21*dconjg(Zn31
     &      ))
      AAABC(601) = EE/SW*(2D0*Zn21*Zp11-Sqrt2*Zn41*Zp21)
      AAABC(602) = EE/SW*(2D0*Zm11*Zn25+Sqrt2*Zm21*dconjg(Zn35
     &      ))
      AAABC(603) = EE/SW*(2D0*Zn25*Zp11-Sqrt2*Zn45*Zp21)
      AAABC(604) = EE/SW*(2D0*Zm12*Zn25+Sqrt2*Zm22*dconjg(Zn35
     &      ))
      AAABC(605) = EE/SW*(2D0*Zn25*Zp12-Sqrt2*Zn45*Zp22)
      AAABR(499) = EE/CW*SW
      AAABR(500) = EE/SW*Sqrt2
      AAABR(501) = EE*(SW/CW+3D0*CW/SW)
      AAABR(502) = EE/SW*Sqrt2*Vcd
      AAABR(503) = EE/SW*Sqrt2*Vcs
      AAABR(504) = EE/SW*Sqrt2*Vud
      AAABR(505) = EE/SW*Sqrt2*Vus
      AAABR(506) = EE*(SW/CW-3D0*CW/SW)
      AAABC(606) = EE/SW*(2D0*Zm11*Zn21+Sqrt2*Zm21*Zn31)
      AAABC(607) = EE/SW*(2D0*Zn21*Zp11-Sqrt2*dconjg(Zn41
     &      )*Zp21)
      AAABC(608) = EE/SW*(2D0*Zm11*Zn23+Sqrt2*Zm21*Zn33)
      AAABC(609) = EE/SW*(2D0*Zn23*Zp11-Sqrt2*dconjg(Zn43
     &      )*Zp21)
      AAABC(610) = EE/SW*(2D0*Zm11*Zn22+Sqrt2*Zm21*Zn32)
      AAABC(611) = EE/SW*(2D0*Zn22*Zp11-Sqrt2*dconjg(Zn42
     &      )*Zp21)
      AAABC(612) = EE/SW*(2D0*Zm11*Zn24+Sqrt2*Zm21*Zn34)
      AAABC(613) = EE/SW*(2D0*Zn24*Zp11-Sqrt2*dconjg(Zn44
     &      )*Zp21)
      AAABC(614) = EE/SW*(2D0*Zm11*Zn25+Sqrt2*Zm21*Zn35)
      AAABC(615) = EE/SW*(2D0*Zn25*Zp11-Sqrt2*dconjg(Zn45
     &      )*Zp21)
      AAABC(616) = EE/SW*(2D0*Zm12*Zn24+Sqrt2*Zm22*Zn34)
      AAABC(617) = EE/SW*(2D0*Zn24*Zp12-Sqrt2*dconjg(Zn44
     &      )*Zp22)
      AAABC(618) = EE/SW*(2D0*Zm12*Zn22+Sqrt2*Zm22*Zn32)
      AAABC(619) = EE/SW*(2D0*Zn22*Zp12-Sqrt2*dconjg(Zn42
     &      )*Zp22)
      AAABC(620) = EE/SW*(2D0*Zm12*Zn23+Sqrt2*Zm22*Zn33)
      AAABC(621) = EE/SW*(2D0*Zn23*Zp12-Sqrt2*dconjg(Zn43
     &      )*Zp22)
      AAABC(622) = EE/SW*(2D0*Zm12*Zn21+Sqrt2*Zm22*Zn31)
      AAABC(623) = EE/SW*(2D0*Zn21*Zp12-Sqrt2*dconjg(Zn41
     &      )*Zp22)
      AAABC(624) = EE/SW*(2D0*Zm12*Zn25+Sqrt2*Zm22*Zn35)
      AAABC(625) = EE/SW*(2D0*Zn25*Zp12-Sqrt2*dconjg(Zn45
     &      )*Zp22)
      AAABC(626) = EE*(SW/CW*dconjg(Zn31)*Zn32-SW/CW*dconjg(Zn41
     &      )*Zn42+CW/SW*dconjg(Zn31)*Zn32-CW/SW*dconjg(Zn41
     &      )*Zn42)
      AAABC(627) = EE*(SW/CW*dconjg(Zn31)*Zn33-SW/CW*dconjg(Zn41
     &      )*Zn43+CW/SW*dconjg(Zn31)*Zn33-CW/SW*dconjg(Zn41
     &      )*Zn43)
      AAABC(628) = EE*(SW/CW*dconjg(Zn31)*Zn34-SW/CW*dconjg(Zn41
     &      )*Zn44+CW/SW*dconjg(Zn31)*Zn34-CW/SW*dconjg(Zn41
     &      )*Zn44)
      AAABC(629) = EE*(SW/CW*Zn31*dconjg(Zn31)-SW/CW*Zn41
     &      *dconjg(Zn41)+CW/SW*Zn31*dconjg(Zn31)-CW/SW*Zn41
     &      *dconjg(Zn41))
      AAABC(630) = EE*(SW/CW*dconjg(Zn31)*Zn35-SW/CW*dconjg(Zn41
     &      )*Zn45+CW/SW*dconjg(Zn31)*Zn35-CW/SW*dconjg(Zn41
     &      )*Zn45)
      AAABC(631) = EE*(SW/CW*Zn32*dconjg(Zn32)-SW/CW*Zn42
     &      *dconjg(Zn42)+CW/SW*Zn32*dconjg(Zn32)-CW/SW*Zn42
     &      *dconjg(Zn42))
      AAABC(632) = EE*(SW/CW*dconjg(Zn32)*Zn34-SW/CW*dconjg(Zn42
     &      )*Zn44+CW/SW*dconjg(Zn32)*Zn34-CW/SW*dconjg(Zn42
     &      )*Zn44)
      AAABC(633) = EE*(SW/CW*dconjg(Zn32)*Zn33-SW/CW*dconjg(Zn42
     &      )*Zn43+CW/SW*dconjg(Zn32)*Zn33-CW/SW*dconjg(Zn42
     &      )*Zn43)
      AAABC(634) = EE*(SW/CW*dconjg(Zn32)*Zn35-SW/CW*dconjg(Zn42
     &      )*Zn45+CW/SW*dconjg(Zn32)*Zn35-CW/SW*dconjg(Zn42
     &      )*Zn45)
      AAABC(635) = EE*(SW/CW*dconjg(Zn33)*Zn34-SW/CW*dconjg(Zn43
     &      )*Zn44+CW/SW*dconjg(Zn33)*Zn34-CW/SW*dconjg(Zn43
     &      )*Zn44)
      AAABC(636) = EE*(SW/CW*Zn33*dconjg(Zn33)-SW/CW*Zn43
     &      *dconjg(Zn43)+CW/SW*Zn33*dconjg(Zn33)-CW/SW*Zn43
     &      *dconjg(Zn43))
      AAABC(637) = EE*(SW/CW*dconjg(Zn33)*Zn35-SW/CW*dconjg(Zn43
     &      )*Zn45+CW/SW*dconjg(Zn33)*Zn35-CW/SW*dconjg(Zn43
     &      )*Zn45)
      AAABC(638) = EE*(SW/CW*Zn34*dconjg(Zn34)-SW/CW*Zn44
     &      *dconjg(Zn44)+CW/SW*Zn34*dconjg(Zn34)-CW/SW*Zn44
     &      *dconjg(Zn44))
      AAABC(639) = EE*(SW/CW*dconjg(Zn34)*Zn35-SW/CW*dconjg(Zn44
     &      )*Zn45+CW/SW*dconjg(Zn34)*Zn35-CW/SW*dconjg(Zn44
     &      )*Zn45)
      AAABC(640) = EE*(SW/CW*Zn35*dconjg(Zn35)-SW/CW*Zn45
     &      *dconjg(Zn45)+CW/SW*Zn35*dconjg(Zn35)-CW/SW*Zn45
     &      *dconjg(Zn45))
      AAABR(507) = EE**2*(SW/CW-CW/SW)
      AAABR(508) = EE**2
      AAABR(509) = EE**2*(2D0-SW**2/CW**2-CW**2/SW**2)
      AAABR(510) = EE**2/SW**2
      AAABR(511) = EE**2/CW*(Za13*sb+Za23*cb)
      AAABR(512) = EE**2/SW*(Za13*sb+Za23*cb)
      AAABR(513) = EE**2/CW*(Za12*sb+Za22*cb)
      AAABR(514) = EE**2/SW*(Za12*sb+Za22*cb)
      AAABR(515) = EE**2/CW*(Za11*sb+Za21*cb)
      AAABR(516) = EE**2/SW*(Za11*sb+Za21*cb)
      AAABR(517) = EE**2/CW*(Zh12*sb-Zh22*cb)
      AAABR(518) = EE**2/SW*(Zh12*sb-Zh22*cb)
      AAABR(519) = EE**2/CW*(Zh11*sb-Zh21*cb)
      AAABR(520) = EE**2/SW*(Zh11*sb-Zh21*cb)
      AAABR(521) = EE**2/CW*(Zh13*sb-Zh23*cb)
      AAABR(522) = EE**2/SW*(Zh13*sb-Zh23*cb)
      AAABR(523) = EE**2*(SW/CW*Zl11**2+2D0*SW/CW*Zl41**2
     &      -CW/SW*Zl11**2)
      AAABR(524) = EE**2*(SW/CW*Zl33**2+2D0*SW/CW*Zl63**2
     &      -CW/SW*Zl33**2)
      AAABR(525) = EE**2*(SW/CW*Zl22**2+2D0*SW/CW*Zl52**2
     &      -CW/SW*Zl22**2)
      AAABR(526) = EE**2*(2D0*Zl22**2-SW**2/CW**2*Zl22**2
     &      -4D0*SW**2/CW**2*Zl52**2-CW**2/SW**2*Zl22**2)
      AAABR(527) = EE**2*(2D0*Zl33**2-SW**2/CW**2*Zl33**2
     &      -4D0*SW**2/CW**2*Zl63**2-CW**2/SW**2*Zl33**2)
      AAABR(528) = EE**2*(2D0*Zl11**2-SW**2/CW**2*Zl11**2
     &      -4D0*SW**2/CW**2*Zl41**2-CW**2/SW**2*Zl11**2)
      AAABR(529) = EE**2/SW**2*Zl33**2
      AAABR(530) = EE**2/SW**2*Zl22**2
      AAABR(531) = EE**2/SW**2*Zl11**2
      AAABR(532) = EE**2*(SW/CW*Zl11*Zl14+2D0*SW/CW*Zl41
     &      *Zl44-CW/SW*Zl11*Zl14)
      AAABR(533) = EE**2*(SW/CW*Zl33*Zl36+2D0*SW/CW*Zl63
     &      *Zl66-CW/SW*Zl33*Zl36)
      AAABR(534) = EE**2*(SW/CW*Zl22*Zl25+2D0*SW/CW*Zl52
     &      *Zl55-CW/SW*Zl22*Zl25)
      AAABR(535) = EE**2*(2D0*Zl22*Zl25-SW**2/CW**2*Zl22
     &      *Zl25-4D0*SW**2/CW**2*Zl52*Zl55-CW**2/SW**2*Zl22
     &      *Zl25)
      AAABR(536) = EE**2*(2D0*Zl33*Zl36-SW**2/CW**2*Zl33
     &      *Zl36-4D0*SW**2/CW**2*Zl63*Zl66-CW**2/SW**2*Zl33
     &      *Zl36)
      AAABR(537) = EE**2*(2D0*Zl11*Zl14-SW**2/CW**2*Zl11
     &      *Zl14-4D0*SW**2/CW**2*Zl41*Zl44-CW**2/SW**2*Zl11
     &      *Zl14)
      AAABR(538) = EE**2/SW**2*Zl11*Zl14
      AAABR(539) = EE**2/SW**2*Zl22*Zl25
      AAABR(540) = EE**2/SW**2*Zl33*Zl36
      AAABR(541) = EE**2/CW*Sqrt2*Zl22
      AAABR(542) = EE**2/SW*Sqrt2*Zl22
      AAABR(543) = EE**2/CW*Sqrt2*Zl33
      AAABR(544) = EE**2/SW*Sqrt2*Zl33
      AAABR(545) = EE**2/CW*Sqrt2*Zl11
      AAABR(546) = EE**2/SW*Sqrt2*Zl11
      AAABR(547) = EE**2*(SW/CW*Zl25**2+2D0*SW/CW*Zl55**2
     &      -CW/SW*Zl25**2)
      AAABR(548) = EE**2*(SW/CW*Zl36**2+2D0*SW/CW*Zl66**2
     &      -CW/SW*Zl36**2)
      AAABR(549) = EE**2*(SW/CW*Zl14**2+2D0*SW/CW*Zl44**2
     &      -CW/SW*Zl14**2)
      AAABR(550) = EE**2*(2D0*Zl36**2-SW**2/CW**2*Zl36**2
     &      -4D0*SW**2/CW**2*Zl66**2-CW**2/SW**2*Zl36**2)
      AAABR(551) = EE**2*(2D0*Zl14**2-SW**2/CW**2*Zl14**2
     &      -4D0*SW**2/CW**2*Zl44**2-CW**2/SW**2*Zl14**2)
      AAABR(552) = EE**2*(2D0*Zl25**2-SW**2/CW**2*Zl25**2
     &      -4D0*SW**2/CW**2*Zl55**2-CW**2/SW**2*Zl25**2)
      AAABR(553) = EE**2/SW**2*Zl36**2
      AAABR(554) = EE**2/SW**2*Zl14**2
      AAABR(555) = EE**2/SW**2*Zl25**2
      AAABR(556) = EE**2/CW*Sqrt2*Zl14
      AAABR(557) = EE**2/SW*Sqrt2*Zl14
      AAABR(558) = EE**2/CW*Sqrt2*Zl25
      AAABR(559) = EE**2/SW*Sqrt2*Zl25
      AAABR(560) = EE**2/CW*Sqrt2*Zl36
      AAABR(561) = EE**2/SW*Sqrt2*Zl36
      AAABR(562) = EE**2*(2D0+SW**2/CW**2+CW**2/SW**2)
      AAABR(563) = EE**2*(SW/CW*Zd11**2-2D0*SW/CW*Zd41**2
     &      +3D0*CW/SW*Zd11**2)
      AAABR(564) = EE**2*(SW/CW*Zd33**2-2D0*SW/CW*Zd63**2
     &      +3D0*CW/SW*Zd33**2)
      AAABR(565) = EE**2*(SW/CW*Zd22**2-2D0*SW/CW*Zd52**2
     &      +3D0*CW/SW*Zd22**2)
      AAABR(566) = EE*GG*(SW/CW*Zd33**2-2D0*SW/CW*Zd63**2
     &      +3D0*CW/SW*Zd33**2)
      AAABR(567) = EE*GG
      AAABR(568) = EE*GG*(SW/CW*Zd11**2-2D0*SW/CW*Zd41**2
     &      +3D0*CW/SW*Zd11**2)
      AAABR(569) = EE*GG*(SW/CW*Zd22**2-2D0*SW/CW*Zd52**2
     &      +3D0*CW/SW*Zd22**2)
      AAABR(570) = EE**2*(6D0*Zd22**2+SW**2/CW**2*Zd22**2
     &      +4D0*SW**2/CW**2*Zd52**2+9D0*CW**2/SW**2*Zd22
     &      **2)
      AAABR(571) = EE**2*(6D0*Zd11**2+SW**2/CW**2*Zd11**2
     &      +4D0*SW**2/CW**2*Zd41**2+9D0*CW**2/SW**2*Zd11
     &      **2)
      AAABR(572) = EE**2*(6D0*Zd33**2+SW**2/CW**2*Zd33**2
     &      +4D0*SW**2/CW**2*Zd63**2+9D0*CW**2/SW**2*Zd33
     &      **2)
      AAABR(573) = GG**2
      AAABR(574) = EE**2/SW**2*Zd33**2
      AAABR(575) = EE**2/SW**2*Zd22**2
      AAABR(576) = EE**2/SW**2*Zd11**2
      AAABR(577) = EE**2*(SW/CW*Zd33*Zd36-2D0*SW/CW*Zd63
     &      *Zd66+3D0*CW/SW*Zd33*Zd36)
      AAABR(578) = EE**2*(SW/CW*Zd11*Zd14-2D0*SW/CW*Zd41
     &      *Zd44+3D0*CW/SW*Zd11*Zd14)
      AAABR(579) = EE**2*(SW/CW*Zd22*Zd25-2D0*SW/CW*Zd52
     &      *Zd55+3D0*CW/SW*Zd22*Zd25)
      AAABR(580) = EE*GG*(SW/CW*Zd33*Zd36-2D0*SW/CW*Zd63
     &      *Zd66+3D0*CW/SW*Zd33*Zd36)
      AAABR(581) = EE*GG*(SW/CW*Zd22*Zd25-2D0*SW/CW*Zd52
     &      *Zd55+3D0*CW/SW*Zd22*Zd25)
      AAABR(582) = EE*GG*(SW/CW*Zd11*Zd14-2D0*SW/CW*Zd41
     &      *Zd44+3D0*CW/SW*Zd11*Zd14)
      AAABR(583) = EE**2*(6D0*Zd33*Zd36+SW**2/CW**2*Zd33
     &      *Zd36+4D0*SW**2/CW**2*Zd63*Zd66+9D0*CW**2/SW
     &      **2*Zd33*Zd36)
      AAABR(584) = EE**2*(6D0*Zd22*Zd25+SW**2/CW**2*Zd22
     &      *Zd25+4D0*SW**2/CW**2*Zd52*Zd55+9D0*CW**2/SW
     &      **2*Zd22*Zd25)
      AAABR(585) = EE**2*(6D0*Zd11*Zd14+SW**2/CW**2*Zd11
     &      *Zd14+4D0*SW**2/CW**2*Zd41*Zd44+9D0*CW**2/SW
     &      **2*Zd11*Zd14)
      AAABR(586) = EE**2/SW**2*Zd33*Zd36
      AAABR(587) = EE**2/SW**2*Zd11*Zd14
      AAABR(588) = EE**2/SW**2*Zd22*Zd25
      AAABR(589) = EE**2/CW*Sqrt2*Vud*Zd11*Zu11
      AAABR(590) = EE**2/SW*Sqrt2*Vud*Zd11*Zu11
      AAABR(591) = EE**2/CW*Sqrt2*Vus*Zd22*Zu11
      AAABR(592) = EE**2/SW*Sqrt2*Vus*Zd22*Zu11
      AAABR(593) = EE**2/CW*Sqrt2*Vcs*Zd22*Zu22
      AAABR(594) = EE**2/SW*Sqrt2*Vcs*Zd22*Zu22
      AAABR(595) = EE**2/CW*Sqrt2*Vcd*Zd11*Zu22
      AAABR(596) = EE**2/SW*Sqrt2*Vcd*Zd11*Zu22
      AAABR(597) = EE**2/CW*Sqrt2*Zd33*Zu33
      AAABR(598) = EE**2/SW*Sqrt2*Zd33*Zu33
      AAABR(599) = EE*GG/SW*Sqrt2*Vcs*Zd22*Zu22
      AAABR(600) = EE*GG/SW*Sqrt2*Vcd*Zd11*Zu22
      AAABR(601) = EE*GG/SW*Sqrt2*Vus*Zd22*Zu11
      AAABR(602) = EE*GG/SW*Sqrt2*Zd33*Zu33
      AAABR(603) = EE*GG/SW*Sqrt2*Vud*Zd11*Zu11
      AAABR(604) = EE**2/CW*Sqrt2*Vud*Zd11*Zu14
      AAABR(605) = EE**2/SW*Sqrt2*Vud*Zd11*Zu14
      AAABR(606) = EE**2/CW*Sqrt2*Vcd*Zd11*Zu25
      AAABR(607) = EE**2/SW*Sqrt2*Vcd*Zd11*Zu25
      AAABR(608) = EE**2/CW*Sqrt2*Zd33*Zu36
      AAABR(609) = EE**2/SW*Sqrt2*Zd33*Zu36
      AAABR(610) = EE**2/CW*Sqrt2*Vus*Zd22*Zu14
      AAABR(611) = EE**2/SW*Sqrt2*Vus*Zd22*Zu14
      AAABR(612) = EE**2/CW*Sqrt2*Vcs*Zd22*Zu25
      AAABR(613) = EE**2/SW*Sqrt2*Vcs*Zd22*Zu25
      AAABR(614) = EE*GG/SW*Sqrt2*Vcd*Zd11*Zu25
      AAABR(615) = EE*GG/SW*Sqrt2*Vus*Zd22*Zu14
      AAABR(616) = EE*GG/SW*Sqrt2*Vud*Zd11*Zu14
      AAABR(617) = EE*GG/SW*Sqrt2*Zd33*Zu36
      AAABR(618) = EE*GG/SW*Sqrt2*Vcs*Zd22*Zu25
      AAABR(619) = EE**2*(SW/CW*Zd14**2-2D0*SW/CW*Zd44**2
     &      +3D0*CW/SW*Zd14**2)
      AAABR(620) = EE**2*(SW/CW*Zd36**2-2D0*SW/CW*Zd66**2
     &      +3D0*CW/SW*Zd36**2)
      AAABR(621) = EE**2*(SW/CW*Zd25**2-2D0*SW/CW*Zd55**2
     &      +3D0*CW/SW*Zd25**2)
      AAABR(622) = EE*GG*(SW/CW*Zd14**2-2D0*SW/CW*Zd44**2
     &      +3D0*CW/SW*Zd14**2)
      AAABR(623) = EE*GG*(SW/CW*Zd25**2-2D0*SW/CW*Zd55**2
     &      +3D0*CW/SW*Zd25**2)
      AAABR(624) = EE*GG*(SW/CW*Zd36**2-2D0*SW/CW*Zd66**2
     &      +3D0*CW/SW*Zd36**2)
      AAABR(625) = EE**2*(6D0*Zd14**2+SW**2/CW**2*Zd14**2
     &      +4D0*SW**2/CW**2*Zd44**2+9D0*CW**2/SW**2*Zd14
     &      **2)
      AAABR(626) = EE**2*(6D0*Zd36**2+SW**2/CW**2*Zd36**2
     &      +4D0*SW**2/CW**2*Zd66**2+9D0*CW**2/SW**2*Zd36
     &      **2)
      AAABR(627) = EE**2*(6D0*Zd25**2+SW**2/CW**2*Zd25**2
     &      +4D0*SW**2/CW**2*Zd55**2+9D0*CW**2/SW**2*Zd25
     &      **2)
      AAABR(628) = EE**2/SW**2*Zd25**2
      AAABR(629) = EE**2/SW**2*Zd14**2
      AAABR(630) = EE**2/SW**2*Zd36**2
      AAABR(631) = EE**2/CW*Sqrt2*Zd36*Zu33
      AAABR(632) = EE**2/SW*Sqrt2*Zd36*Zu33
      AAABR(633) = EE**2/CW*Sqrt2*Vud*Zd14*Zu11
      AAABR(634) = EE**2/SW*Sqrt2*Vud*Zd14*Zu11
      AAABR(635) = EE**2/CW*Sqrt2*Vcs*Zd25*Zu22
      AAABR(636) = EE**2/SW*Sqrt2*Vcs*Zd25*Zu22
      AAABR(637) = EE**2/CW*Sqrt2*Vcd*Zd14*Zu22
      AAABR(638) = EE**2/SW*Sqrt2*Vcd*Zd14*Zu22
      AAABR(639) = EE**2/CW*Sqrt2*Vus*Zd25*Zu11
      AAABR(640) = EE**2/SW*Sqrt2*Vus*Zd25*Zu11
      AAABR(641) = EE*GG/SW*Sqrt2*Zd36*Zu33
      AAABR(642) = EE*GG/SW*Sqrt2*Vus*Zd25*Zu11
      AAABR(643) = EE*GG/SW*Sqrt2*Vcs*Zd25*Zu22
      AAABR(644) = EE*GG/SW*Sqrt2*Vcd*Zd14*Zu22
      AAABR(645) = EE*GG/SW*Sqrt2*Vud*Zd14*Zu11
      AAABR(646) = EE**2/CW*Sqrt2*Vcd*Zd14*Zu25
      AAABR(647) = EE**2/SW*Sqrt2*Vcd*Zd14*Zu25
      AAABR(648) = EE**2/CW*Sqrt2*Vus*Zd25*Zu14
      AAABR(649) = EE**2/SW*Sqrt2*Vus*Zd25*Zu14
      AAABR(650) = EE**2/CW*Sqrt2*Zd36*Zu36
      AAABR(651) = EE**2/SW*Sqrt2*Zd36*Zu36
      AAABR(652) = EE**2/CW*Sqrt2*Vcs*Zd25*Zu25
      AAABR(653) = EE**2/SW*Sqrt2*Vcs*Zd25*Zu25
      AAABR(654) = EE**2/CW*Sqrt2*Vud*Zd14*Zu14
      AAABR(655) = EE**2/SW*Sqrt2*Vud*Zd14*Zu14
      AAABR(656) = EE*GG/SW*Sqrt2*Vus*Zd25*Zu14
      AAABR(657) = EE*GG/SW*Sqrt2*Vcd*Zd14*Zu25
      AAABR(658) = EE*GG/SW*Sqrt2*Vud*Zd14*Zu14
      AAABR(659) = EE*GG/SW*Sqrt2*Vcs*Zd25*Zu25
      AAABR(660) = EE*GG/SW*Sqrt2*Zd36*Zu36
      AAABR(661) = EE**2*(SW/CW*Zu11**2+4D0*SW/CW*Zu41**2
     &      -3D0*CW/SW*Zu11**2)
      AAABR(662) = EE**2*(SW/CW*Zu33**2+4D0*SW/CW*Zu63**2
     &      -3D0*CW/SW*Zu33**2)
      AAABR(663) = EE**2*(SW/CW*Zu22**2+4D0*SW/CW*Zu52**2
     &      -3D0*CW/SW*Zu22**2)
      AAABR(664) = EE*GG*(SW/CW*Zu11**2+4D0*SW/CW*Zu41**2
     &      -3D0*CW/SW*Zu11**2)
      AAABR(665) = EE*GG*(SW/CW*Zu33**2+4D0*SW/CW*Zu63**2
     &      -3D0*CW/SW*Zu33**2)
      AAABR(666) = EE*GG*(SW/CW*Zu22**2+4D0*SW/CW*Zu52**2
     &      -3D0*CW/SW*Zu22**2)
      AAABR(667) = EE**2*(6D0*Zu22**2-SW**2/CW**2*Zu22**2
     &      -16D0*SW**2/CW**2*Zu52**2-9D0*CW**2/SW**2*Zu22
     &      **2)
      AAABR(668) = EE**2*(6D0*Zu11**2-SW**2/CW**2*Zu11**2
     &      -16D0*SW**2/CW**2*Zu41**2-9D0*CW**2/SW**2*Zu11
     &      **2)
      AAABR(669) = EE**2*(6D0*Zu33**2-SW**2/CW**2*Zu33**2
     &      -16D0*SW**2/CW**2*Zu63**2-9D0*CW**2/SW**2*Zu33
     &      **2)
      AAABR(670) = EE**2/SW**2*Zu33**2
      AAABR(671) = EE**2/SW**2*Zu11**2
      AAABR(672) = EE**2/SW**2*Zu22**2
      AAABR(673) = EE**2*(SW/CW*Zu33*Zu36+4D0*SW/CW*Zu63
     &      *Zu66-3D0*CW/SW*Zu33*Zu36)
      AAABR(674) = EE**2*(SW/CW*Zu11*Zu14+4D0*SW/CW*Zu41
     &      *Zu44-3D0*CW/SW*Zu11*Zu14)
      AAABR(675) = EE**2*(SW/CW*Zu22*Zu25+4D0*SW/CW*Zu52
     &      *Zu55-3D0*CW/SW*Zu22*Zu25)
      AAABR(676) = EE*GG*(SW/CW*Zu33*Zu36+4D0*SW/CW*Zu63
     &      *Zu66-3D0*CW/SW*Zu33*Zu36)
      AAABR(677) = EE*GG*(SW/CW*Zu22*Zu25+4D0*SW/CW*Zu52
     &      *Zu55-3D0*CW/SW*Zu22*Zu25)
      AAABR(678) = EE*GG*(SW/CW*Zu11*Zu14+4D0*SW/CW*Zu41
     &      *Zu44-3D0*CW/SW*Zu11*Zu14)
      AAABR(679) = EE**2*(6D0*Zu22*Zu25-SW**2/CW**2*Zu22
     &      *Zu25-16D0*SW**2/CW**2*Zu52*Zu55-9D0*CW**2/SW
     &      **2*Zu22*Zu25)
      AAABR(680) = EE**2*(6D0*Zu11*Zu14-SW**2/CW**2*Zu11
     &      *Zu14-16D0*SW**2/CW**2*Zu41*Zu44-9D0*CW**2/SW
     &      **2*Zu11*Zu14)
      AAABR(681) = EE**2*(6D0*Zu33*Zu36-SW**2/CW**2*Zu33
     &      *Zu36-16D0*SW**2/CW**2*Zu63*Zu66-9D0*CW**2/SW
     &      **2*Zu33*Zu36)
      AAABR(682) = EE**2/SW**2*Zu11*Zu14
      AAABR(683) = EE**2/SW**2*Zu22*Zu25
      AAABR(684) = EE**2/SW**2*Zu33*Zu36
      AAABR(685) = EE**2*(SW/CW*Zu25**2+4D0*SW/CW*Zu55**2
     &      -3D0*CW/SW*Zu25**2)
      AAABR(686) = EE**2*(SW/CW*Zu14**2+4D0*SW/CW*Zu44**2
     &      -3D0*CW/SW*Zu14**2)
      AAABR(687) = EE**2*(SW/CW*Zu36**2+4D0*SW/CW*Zu66**2
     &      -3D0*CW/SW*Zu36**2)
      AAABR(688) = EE*GG*(SW/CW*Zu25**2+4D0*SW/CW*Zu55**2
     &      -3D0*CW/SW*Zu25**2)
      AAABR(689) = EE*GG*(SW/CW*Zu36**2+4D0*SW/CW*Zu66**2
     &      -3D0*CW/SW*Zu36**2)
      AAABR(690) = EE*GG*(SW/CW*Zu14**2+4D0*SW/CW*Zu44**2
     &      -3D0*CW/SW*Zu14**2)
      AAABR(691) = EE**2*(6D0*Zu36**2-SW**2/CW**2*Zu36**2
     &      -16D0*SW**2/CW**2*Zu66**2-9D0*CW**2/SW**2*Zu36
     &      **2)
      AAABR(692) = EE**2*(6D0*Zu14**2-SW**2/CW**2*Zu14**2
     &      -16D0*SW**2/CW**2*Zu44**2-9D0*CW**2/SW**2*Zu14
     &      **2)
      AAABR(693) = EE**2*(6D0*Zu25**2-SW**2/CW**2*Zu25**2
     &      -16D0*SW**2/CW**2*Zu55**2-9D0*CW**2/SW**2*Zu25
     &      **2)
      AAABR(694) = EE**2/SW**2*Zu14**2
      AAABR(695) = EE**2/SW**2*Zu36**2
      AAABR(696) = EE**2/SW**2*Zu25**2
      AAABR(697) = EE**2/CW*(Za13*cb-Za23*sb)
      AAABR(698) = EE**2/SW*(Za13*cb-Za23*sb)
      AAABR(699) = EE**2/CW*(Za11*cb-Za21*sb)
      AAABR(700) = EE**2/SW*(Za11*cb-Za21*sb)
      AAABR(701) = EE**2/CW*(Za12*cb-Za22*sb)
      AAABR(702) = EE**2/SW*(Za12*cb-Za22*sb)
      AAABR(703) = EE**2/CW*(Zh13*cb+Zh23*sb)
      AAABR(704) = EE**2/SW*(Zh13*cb+Zh23*sb)
      AAABR(705) = EE**2/CW*(Zh11*cb+Zh21*sb)
      AAABR(706) = EE**2/SW*(Zh11*cb+Zh21*sb)
      AAABR(707) = EE**2/CW*(Zh12*cb+Zh22*sb)
      AAABR(708) = EE**2/SW*(Zh12*cb+Zh22*sb)
      AAABR(709) = EE**2*(2D0*Za13**2+2D0*Za23**2+SW**2/CW
     &      **2*Za13**2+SW**2/CW**2*Za23**2+CW**2/SW**2*Za13
     &      **2+CW**2/SW**2*Za23**2)
      AAABR(710) = EE**2/SW**2*(Za13**2+Za23**2)
      AAABR(711) = EE**2*(2D0*Za11*Za13+2D0*Za21*Za23+SW
     &      **2/CW**2*Za11*Za13+SW**2/CW**2*Za21*Za23+CW
     &      **2/SW**2*Za11*Za13+CW**2/SW**2*Za21*Za23)
      AAABR(712) = EE**2*(2D0*Za12*Za13+2D0*Za22*Za23+SW
     &      **2/CW**2*Za12*Za13+SW**2/CW**2*Za22*Za23+CW
     &      **2/SW**2*Za12*Za13+CW**2/SW**2*Za22*Za23)
      AAABR(713) = EE**2/SW**2*(Za12*Za13+Za22*Za23)
      AAABR(714) = EE**2/SW**2*(Za11*Za13+Za21*Za23)
      AAABR(715) = EE**2*(2D0*Za11**2+2D0*Za21**2+SW**2/CW
     &      **2*Za11**2+SW**2/CW**2*Za21**2+CW**2/SW**2*Za11
     &      **2+CW**2/SW**2*Za21**2)
      AAABR(716) = EE**2*(2D0*Za11*Za12+2D0*Za21*Za22+SW
     &      **2/CW**2*Za11*Za12+SW**2/CW**2*Za21*Za22+CW
     &      **2/SW**2*Za11*Za12+CW**2/SW**2*Za21*Za22)
      AAABR(717) = EE**2*(2D0*Za12**2+2D0*Za22**2+SW**2/CW
     &      **2*Za12**2+SW**2/CW**2*Za22**2+CW**2/SW**2*Za12
     &      **2+CW**2/SW**2*Za22**2)
      AAABR(718) = EE**2/SW**2*(Za11**2+Za21**2)
      AAABR(719) = EE**2/SW**2*(Za11*Za12+Za21*Za22)
      AAABR(720) = EE**2/SW**2*(Za12**2+Za22**2)
      AAABR(721) = EE**2*(2D0*Zh12*Zh13+2D0*Zh22*Zh23+SW
     &      **2/CW**2*Zh12*Zh13+SW**2/CW**2*Zh22*Zh23+CW
     &      **2/SW**2*Zh12*Zh13+CW**2/SW**2*Zh22*Zh23)
      AAABR(722) = EE**2*(2D0*Zh11*Zh12+2D0*Zh21*Zh22+SW
     &      **2/CW**2*Zh11*Zh12+SW**2/CW**2*Zh21*Zh22+CW
     &      **2/SW**2*Zh11*Zh12+CW**2/SW**2*Zh21*Zh22)
      AAABR(723) = EE**2*(2D0*Zh13**2+2D0*Zh23**2+SW**2/CW
     &      **2*Zh13**2+SW**2/CW**2*Zh23**2+CW**2/SW**2*Zh13
     &      **2+CW**2/SW**2*Zh23**2)
      AAABR(724) = EE**2*(2D0*Zh11*Zh13+2D0*Zh21*Zh23+SW
     &      **2/CW**2*Zh11*Zh13+SW**2/CW**2*Zh21*Zh23+CW
     &      **2/SW**2*Zh11*Zh13+CW**2/SW**2*Zh21*Zh23)
      AAABR(725) = EE**2*(2D0*Zh12**2+2D0*Zh22**2+SW**2/CW
     &      **2*Zh12**2+SW**2/CW**2*Zh22**2+CW**2/SW**2*Zh12
     &      **2+CW**2/SW**2*Zh22**2)
      AAABR(726) = EE**2*(2D0*Zh11**2+2D0*Zh21**2+SW**2/CW
     &      **2*Zh11**2+SW**2/CW**2*Zh21**2+CW**2/SW**2*Zh11
     &      **2+CW**2/SW**2*Zh21**2)
      AAABR(727) = EE**2/SW**2*(Zh12**2+Zh22**2)
      AAABR(728) = EE**2/SW**2*(Zh11*Zh13+Zh21*Zh23)
      AAABR(729) = EE**2/SW**2*(Zh11*Zh12+Zh21*Zh22)
      AAABR(730) = EE**2/SW**2*(Zh11**2+Zh21**2)
      AAABR(731) = EE**2/SW**2*(Zh13**2+Zh23**2)
      AAABR(732) = EE**2/SW**2*(Zh12*Zh13+Zh22*Zh23)
      AAABR(733) = CW*EE**2/SW
      AAABR(734) = CW**2*EE**2/SW**2
      end

      subroutine mtrini
      implicit none
#include "model.h"

      integer m1,m2,m3,m4

      do m1 = 1,3
      do m2 = 1,3
      MTR006(m1,m2)=0D0
      MTR007(m1,m2)=0D0
      MTR008(m1,m2)=0D0
      MTR009(m1,m2)=0D0
      MTR012(m1,m2)=0D0
      MTR013(m1,m2)=0D0
      MTR014(m1,m2)=0D0
      MTR015(m1,m2)=0D0
      MTR019(m1,m2)=0D0
      MTR020(m1,m2)=0D0
      MTR022(m1,m2)=0D0
      MTR023(m1,m2)=0D0
      MTR024(m1,m2)=0D0
      MTR025(m1,m2)=0D0
      MTR028(m1,m2)=0D0
      MTR029(m1,m2)=0D0
      MTR030(m1,m2)=0D0
      MTR031(m1,m2)=0D0
      MTR032(m1,m2)=0D0
      MTR033(m1,m2)=0D0
      MTR034(m1,m2)=0D0
      MTR037(m1,m2)=0D0
      MTR038(m1,m2)=0D0
      MTR039(m1,m2)=0D0
      MTR040(m1,m2)=0D0
      MTR041(m1,m2)=0D0
      MTR056(m1,m2)=0D0
      MTR057(m1,m2)=0D0
      MTR059(m1,m2)=0D0
      MTR060(m1,m2)=0D0
      MTR061(m1,m2)=0D0
      MTR062(m1,m2)=0D0
      MTR065(m1,m2)=0D0
      MTR066(m1,m2)=0D0
      MTR095(m1,m2)=0D0
      MTR115(m1,m2)=0D0
      MTR116(m1,m2)=0D0
      MTR117(m1,m2)=0D0
      MTR118(m1,m2)=0D0
      MTR119(m1,m2)=0D0
      MTR140(m1,m2)=0D0
      MTR141(m1,m2)=0D0
      MTR142(m1,m2)=0D0
      MTR143(m1,m2)=0D0
      MTR146(m1,m2)=0D0
      MTR209(m1,m2)=0D0
      MTR210(m1,m2)=0D0
      MTR244(m1,m2)=0D0
      MTR245(m1,m2)=0D0
      MTR246(m1,m2)=0D0
      MTR247(m1,m2)=0D0
      MTR248(m1,m2)=0D0
      MTR249(m1,m2)=0D0
      MTR254(m1,m2)=0D0
      MTR255(m1,m2)=0D0
      MTR256(m1,m2)=0D0
      MTR257(m1,m2)=0D0
      MTR258(m1,m2)=0D0
      MTR259(m1,m2)=0D0
      MTR260(m1,m2)=0D0
      MTR261(m1,m2)=0D0
      MTR262(m1,m2)=0D0
      MTR263(m1,m2)=0D0
      MTR264(m1,m2)=0D0
      MTR265(m1,m2)=0D0
      MTR274(m1,m2)=0D0
      MTR275(m1,m2)=0D0
      MTR276(m1,m2)=0D0
      MTR277(m1,m2)=0D0
      MTR278(m1,m2)=0D0
      MTR279(m1,m2)=0D0
      MTR292(m1,m2)=0D0
      MTR293(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,2
      do m2 = 1,3
      do m3 = 1,3
      MTR085(m1,m2,m3)=0D0
      MTR086(m1,m2,m3)=0D0
      MTR087(m1,m2,m3)=0D0
      MTR088(m1,m2,m3)=0D0
      enddo
      enddo
      enddo
      do m1 = 1,3
      MTR113(m1)=0D0
      MTR124(m1)=0D0
      MTR125(m1)=0D0
      MTR126(m1)=0D0
      MTR127(m1)=0D0
      MTR128(m1)=0D0
      MTR129(m1)=0D0
      MTR130(m1)=0D0
      MTR131(m1)=0D0
      MTR144(m1)=0D0
      MTR147(m1)=0D0
      MTR148(m1)=0D0
      MTR149(m1)=0D0
      MTR150(m1)=0D0
      MTR151(m1)=0D0
      MTR152(m1)=0D0
      MTR153(m1)=0D0
      MTR154(m1)=0D0
      MTR159(m1)=0D0
      MTR160(m1)=0D0
      MTR161(m1)=0D0
      MTR162(m1)=0D0
      MTR167(m1)=0D0
      MTR168(m1)=0D0
      MTR169(m1)=0D0
      MTR170(m1)=0D0
      MTR175(m1)=0D0
      MTR176(m1)=0D0
      MTR177(m1)=0D0
      MTR178(m1)=0D0
      MTR201(m1)=0D0
      MTR202(m1)=0D0
      MTR221(m1)=0D0
      MTR222(m1)=0D0
      MTR223(m1)=0D0
      MTR224(m1)=0D0
      MTR225(m1)=0D0
      MTR226(m1)=0D0
      MTR227(m1)=0D0
      MTR228(m1)=0D0
      MTR229(m1)=0D0
      MTR230(m1)=0D0
      MTR231(m1)=0D0
      MTR232(m1)=0D0
      MTR233(m1)=0D0
      MTR234(m1)=0D0
      MTR235(m1)=0D0
      MTR236(m1)=0D0
      MTR237(m1)=0D0
      MTR238(m1)=0D0
      MTR239(m1)=0D0
      MTR240(m1)=0D0
      MTR241(m1)=0D0
      MTR242(m1)=0D0
      MTR243(m1)=0D0
      MTR250(m1)=0D0
      MTR251(m1)=0D0
      MTR252(m1)=0D0
      MTR253(m1)=0D0
      MTR266(m1)=0D0
      MTR267(m1)=0D0
      MTR268(m1)=0D0
      MTR269(m1)=0D0
      MTR270(m1)=0D0
      MTR271(m1)=0D0
      MTR272(m1)=0D0
      MTR273(m1)=0D0
      MTR280(m1)=0D0
      MTR281(m1)=0D0
      MTR282(m1)=0D0
      MTR283(m1)=0D0
      MTR286(m1)=0D0
      MTR287(m1)=0D0
      enddo
      do m1 = 1,3
      do m2 = 1,2
      MTR114(m1,m2)=0D0
      MTR145(m1,m2)=0D0
      enddo
      enddo
      do m1 = 1,3
      do m2 = 1,2
      do m3 = 1,3
      MTR120(m1,m2,m3)=0D0
      MTR121(m1,m2,m3)=0D0
      MTR122(m1,m2,m3)=0D0
      MTR123(m1,m2,m3)=0D0
      MTR132(m1,m2,m3)=0D0
      MTR133(m1,m2,m3)=0D0
      MTR134(m1,m2,m3)=0D0
      MTR135(m1,m2,m3)=0D0
      enddo
      enddo
      enddo
      quuMass(1)=Mu
      quuMass(2)=Mc
      quuMass(3)=Mt
      qudMass(1)=Md
      qudMass(2)=Ms
      qudMass(3)=Mb
      lpdMass(1)=Me
      lpdMass(2)=Mm
      lpdMass(3)=Ml
      neuMass(1)=MNE1
      neuMass(2)=MNE2
      neuMass(3)=MNE3
      neuMass(4)=MNE4
      chaMass(1)=MC1
      chaMass(2)=MC2
      sluMass(1)=MSne
      sluMass(2)=MSnm
      sluMass(3)=MSnl
      sldMass(1)=MSe1
      sldMass(2)=MSm1
      sldMass(3)=MSl1
      sleMass(1)=MSe2
      sleMass(2)=MSm2
      sleMass(3)=MSl2
      squMass(1)=MSu1
      squMass(2)=MSc1
      squMass(3)=MSt1
      sqvMass(1)=MSu2
      sqvMass(2)=MSc2
      sqvMass(3)=MSt2
      sqdMass(1)=MSd1
      sqdMass(2)=MSs1
      sqdMass(3)=MSb1
      sqeMass(1)=MSd2
      sqeMass(2)=MSs2
      sqeMass(3)=MSb2
      hisMass(1)=Mh1
      hisMass(2)=Mh2
      hisMass(3)=Mh3
      hiaMass(1)=Mha
      hiaMass(2)=Mhb
      MTR001(1)=AAABR(11)
      MTR001(2)=AAABR(9)
      MTR001(3)=AAABR(10)
      MTR002(1)=AAABR(14)
      MTR002(2)=AAABR(16)
      MTR002(3)=AAABR(15)
      MTR003(1)=AAABR(20)
      MTR003(2)=AAABR(18)
      MTR003(3)=AAABR(19)
      MTR004(1)=AAABR(21)
      MTR004(2)=AAABR(22)
      MTR004(3)=AAABR(23)
      MTR005(1)=AAABR(26)
      MTR005(2)=AAABR(24)
      MTR005(3)=AAABR(25)
      MTR006(1,1)=2D0*AAABR(27)
      MTR006(1,2)=AAABR(29)
      MTR006(2,1)=AAABR(28)
      MTR006(2,2)=AAABR(31)
      MTR006(3,3)=AAABR(30)
      MTR007(1,1)=2D0*AAABR(35)
      MTR007(1,2)=AAABR(36)
      MTR007(2,1)=AAABR(33)
      MTR007(2,2)=AAABR(34)
      MTR007(3,3)=AAABR(32)
      MTR008(1,1)=2D0*AAABR(40)
      MTR008(1,2)=AAABR(41)
      MTR008(2,1)=AAABR(38)
      MTR008(2,2)=AAABR(39)
      MTR008(3,3)=AAABR(37)
      MTR009(1,1)=2D0*AAABR(43)
      MTR009(1,2)=AAABR(44)
      MTR009(2,1)=AAABR(45)
      MTR009(2,2)=AAABR(46)
      MTR009(3,3)=AAABR(42)
      MTR010(1)=AAABR(48)
      MTR010(2)=AAABR(49)
      MTR011(1)=AAABR(50)
      MTR011(2)=AAABR(51)
      MTR011(3)=AAABR(52)
      MTR012(1,1)=2D0*AAABR(27)
      MTR012(1,2)=AAABR(28)
      MTR012(2,1)=AAABR(29)
      MTR012(2,2)=AAABR(31)
      MTR012(3,3)=AAABR(30)
      MTR013(1,1)=2D0*AAABR(40)
      MTR013(1,2)=AAABR(38)
      MTR013(2,1)=AAABR(41)
      MTR013(2,2)=AAABR(39)
      MTR013(3,3)=AAABR(37)
      MTR014(1,1)=2D0*AAABR(35)
      MTR014(1,2)=AAABR(33)
      MTR014(2,1)=AAABR(36)
      MTR014(2,2)=AAABR(34)
      MTR014(3,3)=AAABR(32)
      MTR015(1,1)=2D0*AAABR(43)
      MTR015(1,2)=AAABR(45)
      MTR015(2,1)=AAABR(44)
      MTR015(2,2)=AAABR(46)
      MTR015(3,3)=AAABR(42)
      MTR016(1)=AAABR(54)
      MTR016(2)=AAABR(53)
      MTR016(3)=AAABR(55)
      MTR017(1)=AAABR(57)
      MTR017(2)=AAABR(56)
      MTR017(3)=AAABR(58)
      MTR018(1,1)=AAABR(62)
      MTR018(1,2)=AAABR(61)
      MTR018(1,3)=AAABR(59)
      MTR018(2,1)=AAABR(64)
      MTR018(2,2)=AAABR(60)
      MTR018(2,3)=AAABR(63)
      MTR019(1,1)=AAABR(65)
      MTR019(1,2)=AAABR(69)
      MTR019(1,3)=AAABR(66)
      MTR019(2,1)=AAABR(67)
      MTR019(2,2)=AAABR(72)
      MTR019(2,3)=AAABR(68)
      MTR019(3,1)=AAABR(71)
      MTR019(3,2)=AAABR(70)
      MTR019(3,3)=AAABR(73)
      MTR020(1,1)=AAABR(79)
      MTR020(1,2)=AAABR(74)
      MTR020(1,3)=AAABR(75)
      MTR020(2,1)=AAABR(81)
      MTR020(2,2)=AAABR(80)
      MTR020(2,3)=AAABR(82)
      MTR020(3,1)=AAABR(78)
      MTR020(3,2)=AAABR(77)
      MTR020(3,3)=AAABR(76)
      MTR021(1)=AAABR(83)
      MTR021(2)=AAABR(84)
      MTR021(3)=AAABR(85)
      MTR022(1,1)=AAABR(94)
      MTR022(1,2)=AAABR(87)
      MTR022(1,3)=AAABR(92)
      MTR022(2,1)=AAABR(91)
      MTR022(2,2)=AAABR(93)
      MTR022(2,3)=AAABR(86)
      MTR022(3,1)=AAABR(89)
      MTR022(3,2)=AAABR(90)
      MTR022(3,3)=AAABR(88)
      MTR023(1,1)=AAABR(96)
      MTR023(1,2)=AAABR(96)
      MTR023(1,3)=AAABR(96)
      MTR023(2,1)=AAABR(97)
      MTR023(2,2)=AAABR(97)
      MTR023(2,3)=AAABR(97)
      MTR023(3,1)=AAABR(95)
      MTR023(3,2)=AAABR(95)
      MTR023(3,3)=AAABR(95)
      MTR024(1,1)=AAABR(100)
      MTR024(1,2)=AAABR(101)
      MTR024(2,1)=AAABR(98)
      MTR024(2,2)=AAABR(99)
      MTR024(3,3)=AAABR(102)
      MTR025(1,1)=AAABR(104)
      MTR025(1,2)=AAABR(107)
      MTR025(2,1)=AAABR(105)
      MTR025(2,2)=AAABR(106)
      MTR025(3,3)=AAABR(103)
      MTR026(1)=AAABR(110)
      MTR026(2)=AAABR(108)
      MTR026(3)=AAABR(109)
      MTR027(1,1)=AAABR(115)
      MTR027(1,2)=AAABR(116)
      MTR027(1,3)=AAABR(113)
      MTR027(2,1)=AAABR(112)
      MTR027(2,2)=AAABR(114)
      MTR027(2,3)=AAABR(111)
      MTR028(1,1)=AAABR(125)
      MTR028(1,2)=AAABR(117)
      MTR028(1,3)=AAABR(120)
      MTR028(2,1)=AAABR(118)
      MTR028(2,2)=AAABR(124)
      MTR028(2,3)=AAABR(123)
      MTR028(3,1)=AAABR(121)
      MTR028(3,2)=AAABR(119)
      MTR028(3,3)=AAABR(122)
      MTR029(1,1)=AAABR(134)
      MTR029(1,2)=AAABR(127)
      MTR029(1,3)=AAABR(126)
      MTR029(2,1)=AAABR(133)
      MTR029(2,2)=AAABR(131)
      MTR029(2,3)=AAABR(132)
      MTR029(3,1)=AAABR(130)
      MTR029(3,2)=AAABR(128)
      MTR029(3,3)=AAABR(129)
      MTR030(1,1)=AAABR(136)
      MTR030(1,2)=AAABR(139)
      MTR030(2,1)=AAABR(135)
      MTR030(2,2)=AAABR(137)
      MTR030(3,3)=AAABR(138)
      MTR031(1,1)=AAABR(141)
      MTR031(1,2)=AAABR(144)
      MTR031(2,1)=AAABR(142)
      MTR031(2,2)=AAABR(143)
      MTR031(3,3)=AAABR(140)
      MTR032(1,1)=AAABR(149)
      MTR032(1,2)=AAABR(146)
      MTR032(1,3)=AAABR(153)
      MTR032(2,1)=AAABR(148)
      MTR032(2,2)=AAABR(150)
      MTR032(2,3)=AAABR(152)
      MTR032(3,1)=AAABR(147)
      MTR032(3,2)=AAABR(151)
      MTR032(3,3)=AAABR(145)
      MTR033(1,1)=AAABR(100)
      MTR033(1,2)=AAABR(98)
      MTR033(2,1)=AAABR(101)
      MTR033(2,2)=AAABR(99)
      MTR033(3,3)=AAABR(102)
      MTR034(1,1)=AAABR(136)
      MTR034(1,2)=AAABR(135)
      MTR034(2,1)=AAABR(139)
      MTR034(2,2)=AAABR(137)
      MTR034(3,3)=AAABR(138)
      MTR035(1)=AAABR(154)
      MTR035(2)=AAABR(155)
      MTR035(3)=AAABR(156)
      MTR036(1,1)=AAABR(158)
      MTR036(1,2)=AAABR(157)
      MTR036(1,3)=AAABR(161)
      MTR036(2,1)=AAABR(162)
      MTR036(2,2)=AAABR(160)
      MTR036(2,3)=AAABR(159)
      MTR037(1,1)=AAABR(163)
      MTR037(1,2)=AAABR(169)
      MTR037(1,3)=AAABR(164)
      MTR037(2,1)=AAABR(166)
      MTR037(2,2)=AAABR(171)
      MTR037(2,3)=AAABR(167)
      MTR037(3,1)=AAABR(170)
      MTR037(3,2)=AAABR(168)
      MTR037(3,3)=AAABR(165)
      MTR038(1,1)=AAABR(173)
      MTR038(1,2)=AAABR(178)
      MTR038(1,3)=AAABR(177)
      MTR038(2,1)=AAABR(176)
      MTR038(2,2)=AAABR(175)
      MTR038(2,3)=AAABR(179)
      MTR038(3,1)=AAABR(174)
      MTR038(3,2)=AAABR(180)
      MTR038(3,3)=AAABR(172)
      MTR039(1,1)=AAABR(104)
      MTR039(1,2)=AAABR(105)
      MTR039(2,1)=AAABR(107)
      MTR039(2,2)=AAABR(106)
      MTR039(3,3)=AAABR(103)
      MTR040(1,1)=AAABR(141)
      MTR040(1,2)=AAABR(142)
      MTR040(2,1)=AAABR(144)
      MTR040(2,2)=AAABR(143)
      MTR040(3,3)=AAABR(140)
      MTR041(1,1)=AAABR(183)
      MTR041(1,2)=AAABR(182)
      MTR041(1,3)=AAABR(187)
      MTR041(2,1)=AAABR(181)
      MTR041(2,2)=AAABR(188)
      MTR041(2,3)=AAABR(186)
      MTR041(3,1)=AAABR(185)
      MTR041(3,2)=AAABR(184)
      MTR041(3,3)=AAABR(189)
      MTR042(1)=AAABR(192)
      MTR042(2)=AAABR(190)
      MTR042(3)=AAABR(191)
      MTR043(1)=AAABR(195)
      MTR043(2)=AAABR(194)
      MTR043(3)=AAABR(193)
      MTR044(1,1)=AAABR(197)
      MTR044(1,2)=AAABR(199)
      MTR044(1,3)=AAABR(201)
      MTR044(2,1)=AAABR(200)
      MTR044(2,2)=AAABR(196)
      MTR044(2,3)=AAABR(198)
      MTR045(1,1,1)=AAABR(210)
      MTR045(1,1,2)=AAABR(202)
      MTR045(1,1,3)=AAABR(203)
      MTR045(1,2,1)=AAABR(209)
      MTR045(1,2,2)=AAABR(205)
      MTR045(1,2,3)=AAABR(207)
      MTR045(2,1,1)=MTR045(1,2,1)
      MTR045(2,1,2)=MTR045(1,2,2)
      MTR045(2,1,3)=MTR045(1,2,3)
      MTR045(2,2,1)=AAABR(204)
      MTR045(2,2,2)=AAABR(206)
      MTR045(2,2,3)=AAABR(208)
      MTR046(1,1,1)=AAABR(217)
      MTR046(1,1,2)=AAABR(211)
      MTR046(1,1,3)=AAABR(219)
      MTR046(1,2,1)=MTR046(1,1,2)
      MTR046(1,2,2)=AAABR(212)
      MTR046(1,2,3)=AAABR(213)
      MTR046(1,3,1)=MTR046(1,1,3)
      MTR046(1,3,2)=MTR046(1,2,3)
      MTR046(1,3,3)=AAABR(214)
      MTR046(2,1,1)=MTR046(1,1,2)
      MTR046(2,1,2)=MTR046(1,2,2)
      MTR046(2,1,3)=MTR046(1,2,3)
      MTR046(2,2,1)=MTR046(1,2,2)
      MTR046(2,2,2)=AAABR(220)
      MTR046(2,2,3)=AAABR(218)
      MTR046(2,3,1)=MTR046(1,2,3)
      MTR046(2,3,2)=MTR046(2,2,3)
      MTR046(2,3,3)=AAABR(215)
      MTR046(3,1,1)=MTR046(1,1,3)
      MTR046(3,1,2)=MTR046(1,2,3)
      MTR046(3,1,3)=MTR046(1,3,3)
      MTR046(3,2,1)=MTR046(1,2,3)
      MTR046(3,2,2)=MTR046(2,2,3)
      MTR046(3,2,3)=MTR046(2,3,3)
      MTR046(3,3,1)=MTR046(1,3,3)
      MTR046(3,3,2)=MTR046(2,3,3)
      MTR046(3,3,3)=AAABR(216)
      MTR047(1)=AAABR(224)
      MTR047(2)=AAABR(223)
      MTR048(1)=AAABR(227)
      MTR048(2)=AAABR(226)
      MTR048(3)=AAABR(225)
      MTR049(1)=AAABR(229)
      MTR049(2)=AAABR(230)
      MTR049(3)=AAABR(228)
      MTR050(1)=AAABR(232)
      MTR050(2)=AAABR(231)
      MTR050(3)=AAABR(233)
      MTR051(1)=AAABR(235)
      MTR051(2)=AAABR(236)
      MTR051(3)=AAABR(234)
      MTR052(1)=AAABR(237)
      MTR052(2)=AAABR(239)
      MTR052(3)=AAABR(238)
      MTR053(1)=AAABR(242)
      MTR053(2)=AAABR(241)
      MTR053(3)=AAABR(240)
      MTR054(1)=AAABR(246)
      MTR054(2)=AAABR(245)
      MTR054(3)=AAABR(244)
      MTR055(1)=AAABR(249)
      MTR055(2)=AAABR(247)
      MTR055(3)=AAABR(248)
      MTR056(1,1)=AAABR(254)
      MTR056(1,2)=AAABR(253)
      MTR056(2,1)=AAABR(250)
      MTR056(2,2)=AAABR(251)
      MTR056(3,3)=AAABR(252)
      MTR057(1,1)=AAABR(257)
      MTR057(1,2)=AAABR(255)
      MTR057(2,1)=AAABR(259)
      MTR057(2,2)=AAABR(256)
      MTR057(3,3)=AAABR(258)
      MTR058(1)=AAABR(262)
      MTR058(2)=AAABR(260)
      MTR058(3)=AAABR(261)
      MTR059(1,1)=AAABR(264)
      MTR059(1,2)=AAABR(265)
      MTR059(2,1)=AAABR(266)
      MTR059(2,2)=AAABR(263)
      MTR059(3,3)=AAABR(267)
      MTR060(1,1)=AAABR(272)
      MTR060(1,2)=AAABR(269)
      MTR060(2,1)=AAABR(271)
      MTR060(2,2)=AAABR(268)
      MTR060(3,3)=AAABR(270)
      MTR061(1,1)=AAABR(254)
      MTR061(1,2)=AAABR(250)
      MTR061(2,1)=AAABR(253)
      MTR061(2,2)=AAABR(251)
      MTR061(3,3)=AAABR(252)
      MTR062(1,1)=AAABR(264)
      MTR062(1,2)=AAABR(266)
      MTR062(2,1)=AAABR(265)
      MTR062(2,2)=AAABR(263)
      MTR062(3,3)=AAABR(267)
      MTR063(1)=AAABR(273)
      MTR063(2)=AAABR(275)
      MTR063(3)=AAABR(274)
      MTR064(1)=AAABR(276)
      MTR064(2)=AAABR(277)
      MTR064(3)=AAABR(278)
      MTR065(1,1)=AAABR(257)
      MTR065(1,2)=AAABR(259)
      MTR065(2,1)=AAABR(255)
      MTR065(2,2)=AAABR(256)
      MTR065(3,3)=AAABR(258)
      MTR066(1,1)=AAABR(272)
      MTR066(1,2)=AAABR(271)
      MTR066(2,1)=AAABR(269)
      MTR066(2,2)=AAABR(268)
      MTR066(3,3)=AAABR(270)
      MTR067(1)=AAABR(280)
      MTR067(2)=AAABR(281)
      MTR067(3)=AAABR(279)
      MTR068(1)=AAABR(284)
      MTR068(2)=AAABR(283)
      MTR069(1)=AAABR(285)
      MTR069(2)=AAABR(286)
      MTR069(3)=AAABR(287)
      MTR070(1)=AAABR(290)
      MTR070(2)=AAABR(288)
      MTR070(3)=AAABR(289)
      MTR071(1,1)=AAABR(291)
      MTR071(1,2)=AAABR(296)
      MTR071(1,3)=AAABR(292)
      MTR071(2,1)=AAABR(295)
      MTR071(2,2)=AAABR(293)
      MTR071(2,3)=AAABR(294)
      MTR072(1)=AAABR(298)
      MTR072(2)=AAABR(299)
      MTR072(3)=AAABR(300)
      MTR073(1,1)=AAABR(304)
      MTR073(1,2)=AAABR(303)
      MTR073(2,1)=AAABR(302)
      MTR073(2,2)=AAABR(301)
      MTR074(1,1)=AAABR(304)
      MTR074(1,2)=AAABR(302)
      MTR074(2,1)=AAABR(303)
      MTR074(2,2)=AAABR(301)
      MTR075(1,1,1)=AAABR(312)
      MTR075(1,1,2)=AAABR(305)
      MTR075(1,2,1)=AAABR(311)
      MTR075(1,2,2)=AAABR(306)
      MTR075(2,1,1)=AAABR(310)
      MTR075(2,1,2)=AAABR(307)
      MTR075(2,2,1)=AAABR(308)
      MTR075(2,2,2)=AAABR(309)
      MTR076(1,1,1)=AAABR(312)
      MTR076(1,1,2)=AAABR(305)
      MTR076(1,2,1)=AAABR(310)
      MTR076(1,2,2)=AAABR(307)
      MTR076(2,1,1)=AAABR(311)
      MTR076(2,1,2)=AAABR(306)
      MTR076(2,2,1)=AAABR(308)
      MTR076(2,2,2)=AAABR(309)
      MTR077(1,1,1)=AAABR(315)
      MTR077(1,1,2)=AAABR(319)
      MTR077(1,1,3)=AAABR(324)
      MTR077(1,2,1)=AAABR(320)
      MTR077(1,2,2)=AAABR(316)
      MTR077(1,2,3)=AAABR(322)
      MTR077(2,1,1)=AAABR(321)
      MTR077(2,1,2)=AAABR(317)
      MTR077(2,1,3)=AAABR(323)
      MTR077(2,2,1)=AAABR(313)
      MTR077(2,2,2)=AAABR(314)
      MTR077(2,2,3)=AAABR(318)
      MTR078(1,1,1)=AAABR(315)
      MTR078(1,1,2)=AAABR(319)
      MTR078(1,1,3)=AAABR(324)
      MTR078(1,2,1)=AAABR(321)
      MTR078(1,2,2)=AAABR(317)
      MTR078(1,2,3)=AAABR(323)
      MTR078(2,1,1)=AAABR(320)
      MTR078(2,1,2)=AAABR(316)
      MTR078(2,1,3)=AAABR(322)
      MTR078(2,2,1)=AAABR(313)
      MTR078(2,2,2)=AAABR(314)
      MTR078(2,2,3)=AAABR(318)
      MTR079(1,1)=AAABR(327)
      MTR079(1,2)=AAABR(328)
      MTR079(1,3)=AAABR(329)
      MTR079(2,1)=AAABR(326)
      MTR079(2,2)=AAABR(325)
      MTR079(2,3)=AAABR(330)
      MTR080(1,1)=AAABR(335)
      MTR080(1,2)=AAABR(332)
      MTR080(1,3)=AAABR(334)
      MTR080(2,1)=AAABR(333)
      MTR080(2,2)=AAABR(331)
      MTR080(2,3)=AAABR(336)
      MTR081(1,1)=AAABC(7)
      MTR081(1,2)=AAABC(1)
      MTR081(1,3)=AAABC(5)
      MTR081(1,4)=AAABC(3)
      MTR081(2,1)=AAABC(13)
      MTR081(2,2)=AAABC(11)
      MTR081(2,3)=AAABC(9)
      MTR081(2,4)=AAABC(15)
      MTR082(1,1)=AAABC(8)
      MTR082(1,2)=AAABC(2)
      MTR082(1,3)=AAABC(6)
      MTR082(1,4)=AAABC(4)
      MTR082(2,1)=AAABC(14)
      MTR082(2,2)=AAABC(12)
      MTR082(2,3)=AAABC(10)
      MTR082(2,4)=AAABC(16)
      MTR083(1,1)=AAABC(27)
      MTR083(1,2)=AAABC(29)
      MTR083(1,3)=AAABC(23)
      MTR083(1,4)=AAABC(21)
      MTR083(2,1)=AAABC(25)
      MTR083(2,2)=AAABC(19)
      MTR083(2,3)=AAABC(17)
      MTR083(2,4)=AAABC(31)
      MTR084(1,1)=AAABC(28)
      MTR084(1,2)=AAABC(30)
      MTR084(1,3)=AAABC(24)
      MTR084(1,4)=AAABC(22)
      MTR084(2,1)=AAABC(26)
      MTR084(2,2)=AAABC(20)
      MTR084(2,3)=AAABC(18)
      MTR084(2,4)=AAABC(32)
      MTR085(1,1,1)=2D0*AAABR(345)
      MTR085(1,1,2)=AAABR(348)
      MTR085(1,2,1)=2D0*AAABR(337)
      MTR085(1,2,2)=AAABR(343)
      MTR085(1,3,3)=AAABR(339)
      MTR085(2,1,1)=2D0*AAABR(342)
      MTR085(2,1,2)=AAABR(341)
      MTR085(2,2,1)=2D0*AAABR(351)
      MTR085(2,2,2)=AAABR(346)
      MTR085(2,3,3)=AAABR(349)
      MTR086(1,2,1)=AAABR(338)
      MTR086(1,2,2)=AAABR(344)
      MTR086(1,3,3)=AAABR(340)
      MTR086(2,2,1)=AAABR(352)
      MTR086(2,2,2)=AAABR(347)
      MTR086(2,3,3)=AAABR(350)
      MTR087(1,1,1)=2D0*AAABR(362)
      MTR087(1,1,2)=AAABR(360)
      MTR087(1,2,1)=2D0*AAABR(354)
      MTR087(1,2,2)=AAABR(363)
      MTR087(1,3,3)=AAABR(358)
      MTR087(2,1,1)=2D0*AAABR(353)
      MTR087(2,1,2)=AAABR(361)
      MTR087(2,2,1)=2D0*AAABR(367)
      MTR087(2,2,2)=AAABR(365)
      MTR087(2,3,3)=AAABR(356)
      MTR088(1,2,1)=AAABR(355)
      MTR088(1,2,2)=AAABR(364)
      MTR088(1,3,3)=AAABR(359)
      MTR088(2,2,1)=AAABR(368)
      MTR088(2,2,2)=AAABR(366)
      MTR088(2,3,3)=AAABR(357)
      MTR089(1)=AAABC(33)
      MTR089(2)=AAABC(35)
      MTR090(1)=AAABC(34)
      MTR090(2)=AAABC(36)
      MTR091(1)=AAABC(39)
      MTR091(2)=AAABC(37)
      MTR092(1)=AAABC(40)
      MTR092(2)=AAABC(38)
      MTR093(1)=AAABR(370)
      MTR093(2)=AAABR(369)
      MTR093(3)=AAABR(371)
      MTR094(1,1)=AAABR(374)
      MTR094(1,2)=AAABR(376)
      MTR094(2,1)=AAABR(375)
      MTR094(2,2)=AAABR(372)
      MTR094(3,1)=AAABR(373)
      MTR094(3,2)=AAABR(377)
      MTR095(1,1)=AAABR(383)
      MTR095(1,2)=AAABR(384)
      MTR095(1,3)=AAABR(381)
      MTR095(2,1)=AAABR(378)
      MTR095(2,2)=AAABR(382)
      MTR095(2,3)=AAABR(380)
      MTR095(3,1)=AAABR(385)
      MTR095(3,2)=AAABR(386)
      MTR095(3,3)=AAABR(379)
      MTR096(1)=AAABR(387)
      MTR096(2)=AAABR(388)
      MTR096(3)=AAABR(389)
      MTR097(1)=AAABR(391)
      MTR097(2)=AAABR(392)
      MTR097(3)=AAABR(390)
      MTR098(1,1)=AAABC(43)
      MTR098(1,2)=AAABC(63)
      MTR098(1,3)=AAABC(51)
      MTR098(2,1)=AAABC(41)
      MTR098(2,2)=AAABC(61)
      MTR098(2,3)=AAABC(55)
      MTR098(3,1)=AAABC(47)
      MTR098(3,2)=AAABC(53)
      MTR098(3,3)=AAABC(57)
      MTR098(4,1)=AAABC(59)
      MTR098(4,2)=AAABC(45)
      MTR098(4,3)=AAABC(49)
      MTR099(1,1)=AAABC(44)
      MTR099(1,2)=AAABC(64)
      MTR099(1,3)=AAABC(52)
      MTR099(2,1)=AAABC(42)
      MTR099(2,2)=AAABC(62)
      MTR099(2,3)=AAABC(56)
      MTR099(3,1)=AAABC(48)
      MTR099(3,2)=AAABC(54)
      MTR099(3,3)=AAABC(58)
      MTR099(4,1)=AAABC(60)
      MTR099(4,2)=AAABC(46)
      MTR099(4,3)=AAABC(50)
      MTR100(1,1)=AAABC(85)
      MTR100(1,2)=AAABC(73)
      MTR100(1,3)=AAABC(87)
      MTR100(2,1)=AAABC(81)
      MTR100(2,2)=AAABC(79)
      MTR100(2,3)=AAABC(83)
      MTR100(3,1)=AAABC(65)
      MTR100(3,2)=AAABC(77)
      MTR100(3,3)=AAABC(69)
      MTR100(4,1)=AAABC(71)
      MTR100(4,2)=AAABC(75)
      MTR100(4,3)=AAABC(67)
      MTR101(1,1)=AAABC(86)
      MTR101(1,2)=AAABC(74)
      MTR101(1,3)=AAABC(88)
      MTR101(2,1)=AAABC(82)
      MTR101(2,2)=AAABC(80)
      MTR101(2,3)=AAABC(84)
      MTR101(3,1)=AAABC(66)
      MTR101(3,2)=AAABC(78)
      MTR101(3,3)=AAABC(70)
      MTR101(4,1)=AAABC(72)
      MTR101(4,2)=AAABC(76)
      MTR101(4,3)=AAABC(68)
      MTR102(1,1)=AAABR(393)
      MTR102(1,2)=AAABR(393)
      MTR102(1,3)=AAABR(393)
      MTR102(2,1)=AAABR(397)
      MTR102(2,2)=AAABR(397)
      MTR102(2,3)=AAABR(397)
      MTR103(1,1)=AAABR(394)
      MTR103(1,2)=AAABR(396)
      MTR103(1,3)=AAABR(395)
      MTR103(2,1)=AAABR(399)
      MTR103(2,2)=AAABR(400)
      MTR103(2,3)=AAABR(398)
      MTR104(1)=AAABC(91)
      MTR104(2)=AAABC(89)
      MTR104(3)=AAABC(93)
      MTR105(1)=AAABC(92)
      MTR105(2)=AAABC(90)
      MTR105(3)=AAABC(94)
      MTR106(1)=AAABC(95)
      MTR106(2)=AAABC(97)
      MTR106(3)=AAABC(99)
      MTR107(1)=AAABC(96)
      MTR107(2)=AAABC(98)
      MTR107(3)=AAABC(100)
      MTR108(1,1)=AAABC(104)
      MTR108(1,2)=AAABC(104)
      MTR108(1,3)=AAABC(104)
      MTR108(2,1)=AAABC(103)
      MTR108(2,2)=AAABC(103)
      MTR108(2,3)=AAABC(103)
      MTR108(3,1)=AAABC(102)
      MTR108(3,2)=AAABC(102)
      MTR108(3,3)=AAABC(102)
      MTR108(4,1)=AAABC(101)
      MTR108(4,2)=AAABC(101)
      MTR108(4,3)=AAABC(101)
      MTR109(1,1)=AAABC(115)
      MTR109(1,2)=AAABC(123)
      MTR109(1,3)=AAABC(106)
      MTR109(2,1)=AAABC(109)
      MTR109(2,2)=AAABC(129)
      MTR109(2,3)=AAABC(124)
      MTR109(3,1)=AAABC(111)
      MTR109(3,2)=AAABC(119)
      MTR109(3,3)=AAABC(120)
      MTR109(4,1)=AAABC(117)
      MTR109(4,2)=AAABC(127)
      MTR109(4,3)=AAABC(112)
      MTR110(1,1)=2D0*AAABC(114)
      MTR110(1,2)=AAABC(122)
      MTR110(1,3)=AAABC(107)
      MTR110(2,1)=2D0*AAABC(108)
      MTR110(2,2)=AAABC(128)
      MTR110(2,3)=AAABC(125)
      MTR110(3,1)=2D0*AAABC(110)
      MTR110(3,2)=AAABC(118)
      MTR110(3,3)=AAABC(121)
      MTR110(4,1)=2D0*AAABC(116)
      MTR110(4,2)=AAABC(126)
      MTR110(4,3)=AAABC(113)
      MTR111(1,1)=2D0*AAABC(142)
      MTR111(1,2)=AAABC(146)
      MTR111(1,3)=AAABC(139)
      MTR111(2,1)=2D0*AAABC(144)
      MTR111(2,2)=AAABC(130)
      MTR111(2,3)=AAABC(141)
      MTR111(3,1)=2D0*AAABC(134)
      MTR111(3,2)=AAABC(148)
      MTR111(3,3)=AAABC(151)
      MTR111(4,1)=2D0*AAABC(152)
      MTR111(4,2)=AAABC(136)
      MTR111(4,3)=AAABC(133)
      MTR112(1,1)=AAABC(143)
      MTR112(1,2)=AAABC(147)
      MTR112(1,3)=AAABC(138)
      MTR112(2,1)=AAABC(145)
      MTR112(2,2)=AAABC(131)
      MTR112(2,3)=AAABC(140)
      MTR112(3,1)=AAABC(135)
      MTR112(3,2)=AAABC(149)
      MTR112(3,3)=AAABC(150)
      MTR112(4,1)=AAABC(153)
      MTR112(4,2)=AAABC(137)
      MTR112(4,3)=AAABC(132)
      MTR113(2)=AAABR(401)
      MTR113(3)=AAABR(402)
      MTR114(2,1)=AAABR(406)
      MTR114(2,2)=AAABR(403)
      MTR114(3,1)=AAABR(405)
      MTR114(3,2)=AAABR(404)
      MTR115(2,1)=AAABR(407)
      MTR115(2,2)=AAABR(409)
      MTR115(2,3)=AAABR(408)
      MTR115(3,1)=AAABR(410)
      MTR115(3,2)=AAABR(412)
      MTR115(3,3)=AAABR(411)
      MTR116(1,2)=AAABR(418)
      MTR116(2,2)=AAABR(413)
      MTR116(3,3)=AAABR(416)
      MTR117(2,1)=AAABR(415)
      MTR117(2,2)=AAABR(414)
      MTR117(3,3)=AAABR(417)
      MTR118(2,1)=AAABR(419)
      MTR118(2,2)=AAABR(421)
      MTR118(3,3)=AAABR(423)
      MTR119(1,2)=AAABR(424)
      MTR119(2,2)=AAABR(420)
      MTR119(3,3)=AAABR(422)
      MTR120(1,1,1)=2D0*AAABR(430)
      MTR120(1,1,2)=AAABR(429)
      MTR120(1,2,1)=2D0*AAABR(443)
      MTR120(1,2,2)=AAABR(446)
      MTR120(2,1,1)=2D0*AAABR(431)
      MTR120(2,1,2)=AAABR(425)
      MTR120(2,2,1)=2D0*AAABR(447)
      MTR120(2,2,2)=AAABR(444)
      MTR120(3,1,3)=AAABR(427)
      MTR120(3,2,3)=AAABR(441)
      MTR121(2,1,1)=AAABR(432)
      MTR121(2,1,2)=AAABR(426)
      MTR121(2,2,1)=AAABR(448)
      MTR121(2,2,2)=AAABR(445)
      MTR121(3,1,3)=AAABR(428)
      MTR121(3,2,3)=AAABR(442)
      MTR122(1,1,1)=2D0*AAABR(439)
      MTR122(1,1,2)=AAABR(440)
      MTR122(1,2,1)=2D0*AAABR(453)
      MTR122(1,2,2)=AAABR(456)
      MTR122(2,1,1)=2D0*AAABR(433)
      MTR122(2,1,2)=AAABR(437)
      MTR122(2,2,1)=2D0*AAABR(454)
      MTR122(2,2,2)=AAABR(451)
      MTR122(3,1,3)=AAABR(435)
      MTR122(3,2,3)=AAABR(449)
      MTR123(2,1,1)=AAABR(434)
      MTR123(2,1,2)=AAABR(438)
      MTR123(2,2,1)=AAABR(455)
      MTR123(2,2,2)=AAABR(452)
      MTR123(3,1,3)=AAABR(436)
      MTR123(3,2,3)=AAABR(450)
      MTR124(1)=AAABR(457)
      MTR124(2)=AAABR(459)
      MTR124(3)=AAABR(462)
      MTR125(1)=AAABR(458)
      MTR125(2)=AAABR(460)
      MTR125(3)=AAABR(461)
      MTR126(1)=AAABR(463)
      MTR126(2)=AAABR(467)
      MTR126(3)=AAABR(466)
      MTR127(1)=AAABR(464)
      MTR127(2)=AAABR(468)
      MTR127(3)=AAABR(465)
      MTR128(1)=2D0*AAABC(154)
      MTR128(2)=AAABC(156)
      MTR128(3)=AAABC(159)
      MTR129(1)=AAABC(155)
      MTR129(2)=AAABC(157)
      MTR129(3)=AAABC(158)
      MTR130(1)=AAABC(163)
      MTR130(2)=AAABC(165)
      MTR130(3)=AAABC(160)
      MTR131(1)=2D0*AAABC(162)
      MTR131(2)=AAABC(164)
      MTR131(3)=AAABC(161)
      MTR132(1,1,1)=2D0*AAABR(345)
      MTR132(1,1,2)=AAABR(348)
      MTR132(1,2,1)=2D0*AAABR(342)
      MTR132(1,2,2)=AAABR(341)
      MTR132(2,1,1)=2D0*AAABR(337)
      MTR132(2,1,2)=AAABR(343)
      MTR132(2,2,1)=2D0*AAABR(351)
      MTR132(2,2,2)=AAABR(346)
      MTR132(3,1,3)=AAABR(339)
      MTR132(3,2,3)=AAABR(349)
      MTR133(2,1,1)=AAABR(338)
      MTR133(2,1,2)=AAABR(344)
      MTR133(2,2,1)=AAABR(352)
      MTR133(2,2,2)=AAABR(347)
      MTR133(3,1,3)=AAABR(340)
      MTR133(3,2,3)=AAABR(350)
      MTR134(1,1,1)=2D0*AAABR(362)
      MTR134(1,1,2)=AAABR(360)
      MTR134(1,2,1)=2D0*AAABR(353)
      MTR134(1,2,2)=AAABR(361)
      MTR134(2,1,1)=2D0*AAABR(354)
      MTR134(2,1,2)=AAABR(363)
      MTR134(2,2,1)=2D0*AAABR(367)
      MTR134(2,2,2)=AAABR(365)
      MTR134(3,1,3)=AAABR(358)
      MTR134(3,2,3)=AAABR(356)
      MTR135(2,1,1)=AAABR(355)
      MTR135(2,1,2)=AAABR(364)
      MTR135(2,2,1)=AAABR(368)
      MTR135(2,2,2)=AAABR(366)
      MTR135(3,1,3)=AAABR(359)
      MTR135(3,2,3)=AAABR(357)
      MTR136(1,1)=AAABC(176)
      MTR136(1,2)=AAABC(168)
      MTR136(1,3)=AAABC(166)
      MTR136(2,1)=AAABC(180)
      MTR136(2,2)=AAABC(182)
      MTR136(2,3)=AAABC(170)
      MTR136(3,1)=AAABC(174)
      MTR136(3,2)=AAABC(188)
      MTR136(3,3)=AAABC(184)
      MTR136(4,1)=AAABC(172)
      MTR136(4,2)=AAABC(178)
      MTR136(4,3)=AAABC(186)
      MTR137(1,1)=4D0*AAABC(177)
      MTR137(1,2)=AAABC(169)
      MTR137(1,3)=AAABC(167)
      MTR137(2,1)=4D0*AAABC(181)
      MTR137(2,2)=AAABC(183)
      MTR137(2,3)=AAABC(171)
      MTR137(3,1)=4D0*AAABC(175)
      MTR137(3,2)=AAABC(189)
      MTR137(3,3)=AAABC(185)
      MTR137(4,1)=4D0*AAABC(173)
      MTR137(4,2)=AAABC(179)
      MTR137(4,3)=AAABC(187)
      MTR138(1,1)=AAABC(208)
      MTR138(1,2)=AAABC(206)
      MTR138(1,3)=AAABC(202)
      MTR138(2,1)=AAABC(194)
      MTR138(2,2)=AAABC(210)
      MTR138(2,3)=AAABC(204)
      MTR138(3,1)=AAABC(190)
      MTR138(3,2)=AAABC(192)
      MTR138(3,3)=AAABC(198)
      MTR138(4,1)=AAABC(200)
      MTR138(4,2)=AAABC(196)
      MTR138(4,3)=AAABC(212)
      MTR139(1,1)=4D0*AAABC(209)
      MTR139(1,2)=AAABC(207)
      MTR139(1,3)=AAABC(203)
      MTR139(2,1)=4D0*AAABC(195)
      MTR139(2,2)=AAABC(211)
      MTR139(2,3)=AAABC(205)
      MTR139(3,1)=4D0*AAABC(191)
      MTR139(3,2)=AAABC(193)
      MTR139(3,3)=AAABC(199)
      MTR139(4,1)=4D0*AAABC(201)
      MTR139(4,2)=AAABC(197)
      MTR139(4,3)=AAABC(213)
      MTR140(2,1)=AAABR(418)
      MTR140(2,2)=AAABR(413)
      MTR140(3,3)=AAABR(416)
      MTR141(1,2)=AAABR(415)
      MTR141(2,2)=AAABR(414)
      MTR141(3,3)=AAABR(417)
      MTR142(2,1)=AAABR(424)
      MTR142(2,2)=AAABR(420)
      MTR142(3,3)=AAABR(422)
      MTR143(1,2)=AAABR(419)
      MTR143(2,2)=AAABR(421)
      MTR143(3,3)=AAABR(423)
      MTR144(2)=AAABR(470)
      MTR144(3)=AAABR(469)
      MTR145(2,1)=AAABR(473)
      MTR145(2,2)=AAABR(474)
      MTR145(3,1)=AAABR(471)
      MTR145(3,2)=AAABR(472)
      MTR146(2,1)=AAABR(476)
      MTR146(2,2)=AAABR(477)
      MTR146(2,3)=AAABR(479)
      MTR146(3,1)=AAABR(480)
      MTR146(3,2)=AAABR(478)
      MTR146(3,3)=AAABR(475)
      MTR147(1)=AAABR(481)
      MTR147(2)=AAABR(485)
      MTR147(3)=AAABR(483)
      MTR148(1)=AAABR(482)
      MTR148(2)=AAABR(486)
      MTR148(3)=AAABR(484)
      MTR149(1)=AAABR(491)
      MTR149(2)=AAABR(489)
      MTR149(3)=AAABR(487)
      MTR150(1)=AAABR(492)
      MTR150(2)=AAABR(490)
      MTR150(3)=AAABR(488)
      MTR151(1)=AAABC(214)
      MTR151(2)=AAABC(218)
      MTR151(3)=AAABC(216)
      MTR152(1)=4D0*AAABC(215)
      MTR152(2)=AAABC(219)
      MTR152(3)=AAABC(217)
      MTR153(1)=AAABC(224)
      MTR153(2)=AAABC(222)
      MTR153(3)=AAABC(220)
      MTR154(1)=4D0*AAABC(225)
      MTR154(2)=AAABC(223)
      MTR154(3)=AAABC(221)
      MTR155(1,1)=AAABC(115)
      MTR155(1,2)=AAABC(333)
      MTR155(1,3)=AAABC(230)
      MTR155(2,1)=AAABC(109)
      MTR155(2,2)=AAABC(331)
      MTR155(2,3)=AAABC(232)
      MTR155(3,1)=AAABC(111)
      MTR155(3,2)=AAABC(327)
      MTR155(3,3)=AAABC(226)
      MTR155(4,1)=AAABC(117)
      MTR155(4,2)=AAABC(329)
      MTR155(4,3)=AAABC(228)
      MTR156(1,1)=2D0*AAABC(114)
      MTR156(1,2)=AAABC(332)
      MTR156(1,3)=AAABC(231)
      MTR156(2,1)=2D0*AAABC(108)
      MTR156(2,2)=AAABC(330)
      MTR156(2,3)=AAABC(233)
      MTR156(3,1)=2D0*AAABC(110)
      MTR156(3,2)=AAABC(326)
      MTR156(3,3)=AAABC(227)
      MTR156(4,1)=2D0*AAABC(116)
      MTR156(4,2)=AAABC(328)
      MTR156(4,3)=AAABC(229)
      MTR157(1,1)=AAABC(143)
      MTR157(1,2)=AAABC(341)
      MTR157(1,3)=AAABC(234)
      MTR157(2,1)=AAABC(145)
      MTR157(2,2)=AAABC(339)
      MTR157(2,3)=AAABC(236)
      MTR157(3,1)=AAABC(135)
      MTR157(3,2)=AAABC(337)
      MTR157(3,3)=AAABC(240)
      MTR157(4,1)=AAABC(153)
      MTR157(4,2)=AAABC(335)
      MTR157(4,3)=AAABC(238)
      MTR158(1,1)=2D0*AAABC(142)
      MTR158(1,2)=AAABC(340)
      MTR158(1,3)=AAABC(235)
      MTR158(2,1)=2D0*AAABC(144)
      MTR158(2,2)=AAABC(338)
      MTR158(2,3)=AAABC(237)
      MTR158(3,1)=2D0*AAABC(134)
      MTR158(3,2)=AAABC(336)
      MTR158(3,3)=AAABC(241)
      MTR158(4,1)=2D0*AAABC(152)
      MTR158(4,2)=AAABC(334)
      MTR158(4,3)=AAABC(239)
      MTR159(1)=AAABC(155)
      MTR159(2)=AAABC(343)
      MTR159(3)=AAABC(242)
      MTR160(1)=2D0*AAABC(154)
      MTR160(2)=AAABC(342)
      MTR160(3)=AAABC(243)
      MTR161(1)=AAABC(163)
      MTR161(2)=AAABC(345)
      MTR161(3)=AAABC(244)
      MTR162(1)=2D0*AAABC(162)
      MTR162(2)=AAABC(344)
      MTR162(3)=AAABC(245)
      MTR163(1,1)=AAABC(176)
      MTR163(1,2)=AAABC(246)
      MTR163(1,3)=AAABC(350)
      MTR163(2,1)=AAABC(180)
      MTR163(2,2)=AAABC(248)
      MTR163(2,3)=AAABC(352)
      MTR163(3,1)=AAABC(174)
      MTR163(3,2)=AAABC(250)
      MTR163(3,3)=AAABC(348)
      MTR163(4,1)=AAABC(172)
      MTR163(4,2)=AAABC(252)
      MTR163(4,3)=AAABC(346)
      MTR164(1,1)=4D0*AAABC(177)
      MTR164(1,2)=AAABC(247)
      MTR164(1,3)=AAABC(351)
      MTR164(2,1)=4D0*AAABC(181)
      MTR164(2,2)=AAABC(249)
      MTR164(2,3)=AAABC(353)
      MTR164(3,1)=4D0*AAABC(175)
      MTR164(3,2)=AAABC(251)
      MTR164(3,3)=AAABC(349)
      MTR164(4,1)=4D0*AAABC(173)
      MTR164(4,2)=AAABC(253)
      MTR164(4,3)=AAABC(347)
      MTR165(1,1)=AAABC(208)
      MTR165(1,2)=AAABC(254)
      MTR165(1,3)=AAABC(360)
      MTR165(2,1)=AAABC(194)
      MTR165(2,2)=AAABC(256)
      MTR165(2,3)=AAABC(358)
      MTR165(3,1)=AAABC(190)
      MTR165(3,2)=AAABC(258)
      MTR165(3,3)=AAABC(354)
      MTR165(4,1)=AAABC(200)
      MTR165(4,2)=AAABC(260)
      MTR165(4,3)=AAABC(356)
      MTR166(1,1)=4D0*AAABC(209)
      MTR166(1,2)=AAABC(255)
      MTR166(1,3)=AAABC(361)
      MTR166(2,1)=4D0*AAABC(195)
      MTR166(2,2)=AAABC(257)
      MTR166(2,3)=AAABC(359)
      MTR166(3,1)=4D0*AAABC(191)
      MTR166(3,2)=AAABC(259)
      MTR166(3,3)=AAABC(355)
      MTR166(4,1)=4D0*AAABC(201)
      MTR166(4,2)=AAABC(261)
      MTR166(4,3)=AAABC(357)
      MTR167(1)=AAABC(214)
      MTR167(2)=AAABC(262)
      MTR167(3)=AAABC(362)
      MTR168(1)=4D0*AAABC(215)
      MTR168(2)=AAABC(263)
      MTR168(3)=AAABC(363)
      MTR169(1)=AAABC(224)
      MTR169(2)=AAABC(264)
      MTR169(3)=AAABC(364)
      MTR170(1)=4D0*AAABC(225)
      MTR170(2)=AAABC(265)
      MTR170(3)=AAABC(365)
      MTR171(1,1)=AAABC(268)
      MTR171(1,2)=AAABC(310)
      MTR171(1,3)=AAABC(290)
      MTR171(2,1)=AAABC(272)
      MTR171(2,2)=AAABC(308)
      MTR171(2,3)=AAABC(288)
      MTR171(3,1)=AAABC(266)
      MTR171(3,2)=AAABC(306)
      MTR171(3,3)=AAABC(292)
      MTR171(4,1)=AAABC(270)
      MTR171(4,2)=AAABC(312)
      MTR171(4,3)=AAABC(286)
      MTR172(1,1)=AAABC(269)
      MTR172(1,2)=AAABC(311)
      MTR172(1,3)=AAABC(291)
      MTR172(2,1)=AAABC(273)
      MTR172(2,2)=AAABC(309)
      MTR172(2,3)=AAABC(289)
      MTR172(3,1)=AAABC(267)
      MTR172(3,2)=AAABC(307)
      MTR172(3,3)=AAABC(293)
      MTR172(4,1)=AAABC(271)
      MTR172(4,2)=AAABC(313)
      MTR172(4,3)=AAABC(287)
      MTR173(1,1)=AAABC(276)
      MTR173(1,2)=AAABC(320)
      MTR173(1,3)=AAABC(296)
      MTR173(2,1)=AAABC(278)
      MTR173(2,2)=AAABC(314)
      MTR173(2,3)=AAABC(300)
      MTR173(3,1)=AAABC(274)
      MTR173(3,2)=AAABC(316)
      MTR173(3,3)=AAABC(294)
      MTR173(4,1)=AAABC(280)
      MTR173(4,2)=AAABC(318)
      MTR173(4,3)=AAABC(298)
      MTR174(1,1)=AAABC(277)
      MTR174(1,2)=AAABC(321)
      MTR174(1,3)=AAABC(297)
      MTR174(2,1)=AAABC(279)
      MTR174(2,2)=AAABC(315)
      MTR174(2,3)=AAABC(301)
      MTR174(3,1)=AAABC(275)
      MTR174(3,2)=AAABC(317)
      MTR174(3,3)=AAABC(295)
      MTR174(4,1)=AAABC(281)
      MTR174(4,2)=AAABC(319)
      MTR174(4,3)=AAABC(299)
      MTR175(1)=AAABC(282)
      MTR175(2)=AAABC(322)
      MTR175(3)=AAABC(302)
      MTR176(1)=AAABC(283)
      MTR176(2)=AAABC(323)
      MTR176(3)=AAABC(303)
      MTR177(1)=AAABC(284)
      MTR177(2)=AAABC(324)
      MTR177(3)=AAABC(304)
      MTR178(1)=AAABC(285)
      MTR178(2)=AAABC(325)
      MTR178(3)=AAABC(305)
      MTR179(1,1)=AAABC(370)
      MTR179(1,2)=AAABC(368)
      MTR179(1,3)=AAABC(366)
      MTR179(1,4)=AAABC(372)
      MTR179(2,1)=AAABC(388)
      MTR179(2,2)=AAABC(390)
      MTR179(2,3)=AAABC(386)
      MTR179(2,4)=AAABC(392)
      MTR180(1,1)=AAABC(371)
      MTR180(1,2)=AAABC(369)
      MTR180(1,3)=AAABC(367)
      MTR180(1,4)=AAABC(373)
      MTR180(2,1)=AAABC(389)
      MTR180(2,2)=AAABC(391)
      MTR180(2,3)=AAABC(387)
      MTR180(2,4)=AAABC(393)
      MTR181(1,1)=AAABC(376)
      MTR181(1,2)=AAABC(374)
      MTR181(1,3)=AAABC(380)
      MTR181(1,4)=AAABC(378)
      MTR181(2,1)=AAABC(398)
      MTR181(2,2)=AAABC(400)
      MTR181(2,3)=AAABC(394)
      MTR181(2,4)=AAABC(396)
      MTR182(1,1)=AAABC(377)
      MTR182(1,2)=AAABC(375)
      MTR182(1,3)=AAABC(381)
      MTR182(1,4)=AAABC(379)
      MTR182(2,1)=AAABC(399)
      MTR182(2,2)=AAABC(401)
      MTR182(2,3)=AAABC(395)
      MTR182(2,4)=AAABC(397)
      MTR183(1)=AAABC(382)
      MTR183(2)=AAABC(402)
      MTR184(1)=AAABC(383)
      MTR184(2)=AAABC(403)
      MTR185(1)=AAABC(384)
      MTR185(2)=AAABC(404)
      MTR186(1)=AAABC(385)
      MTR186(2)=AAABC(405)
      MTR187(1,1)=2D0*AAABC(406)
      MTR187(1,2)=AAABC(410)
      MTR187(1,3)=AAABC(412)
      MTR187(1,4)=AAABC(408)
      MTR187(2,1)=AAABC(410)
      MTR187(2,2)=2D0*AAABC(470)
      MTR187(2,3)=AAABC(466)
      MTR187(2,4)=AAABC(468)
      MTR187(3,1)=AAABC(412)
      MTR187(3,2)=AAABC(466)
      MTR187(3,3)=2D0*AAABC(514)
      MTR187(3,4)=AAABC(516)
      MTR187(4,1)=AAABC(408)
      MTR187(4,2)=AAABC(468)
      MTR187(4,3)=AAABC(516)
      MTR187(4,4)=2D0*AAABC(550)
      MTR188(1,1)=2D0*AAABC(407)
      MTR188(1,2)=AAABC(411)
      MTR188(1,3)=AAABC(413)
      MTR188(1,4)=AAABC(409)
      MTR188(2,1)=AAABC(411)
      MTR188(2,2)=2D0*AAABC(471)
      MTR188(2,3)=AAABC(467)
      MTR188(2,4)=AAABC(469)
      MTR188(3,1)=AAABC(413)
      MTR188(3,2)=AAABC(467)
      MTR188(3,3)=2D0*AAABC(515)
      MTR188(3,4)=AAABC(517)
      MTR188(4,1)=AAABC(409)
      MTR188(4,2)=AAABC(469)
      MTR188(4,3)=AAABC(517)
      MTR188(4,4)=2D0*AAABC(551)
      MTR189(1,1,1)=2D0*AAABC(418)
      MTR189(1,1,2)=2D0*AAABC(422)
      MTR189(1,2,1)=AAABC(428)
      MTR189(1,2,2)=AAABC(414)
      MTR189(1,3,1)=AAABC(416)
      MTR189(1,3,2)=AAABC(424)
      MTR189(1,4,1)=AAABC(420)
      MTR189(1,4,2)=AAABC(426)
      MTR189(2,1,1)=AAABC(428)
      MTR189(2,1,2)=AAABC(414)
      MTR189(2,2,1)=2D0*AAABC(478)
      MTR189(2,2,2)=2D0*AAABC(476)
      MTR189(2,3,1)=AAABC(472)
      MTR189(2,3,2)=AAABC(480)
      MTR189(2,4,1)=AAABC(482)
      MTR189(2,4,2)=AAABC(474)
      MTR189(3,1,1)=AAABC(416)
      MTR189(3,1,2)=AAABC(424)
      MTR189(3,2,1)=AAABC(472)
      MTR189(3,2,2)=AAABC(480)
      MTR189(3,3,1)=2D0*AAABC(518)
      MTR189(3,3,2)=2D0*AAABC(522)
      MTR189(3,4,1)=AAABC(520)
      MTR189(3,4,2)=AAABC(524)
      MTR189(4,1,1)=AAABC(420)
      MTR189(4,1,2)=AAABC(426)
      MTR189(4,2,1)=AAABC(482)
      MTR189(4,2,2)=AAABC(474)
      MTR189(4,3,1)=AAABC(520)
      MTR189(4,3,2)=AAABC(524)
      MTR189(4,4,1)=2D0*AAABC(554)
      MTR189(4,4,2)=2D0*AAABC(552)
      MTR190(1,1,1)=2D0*AAABC(419)
      MTR190(1,1,2)=2D0*AAABC(423)
      MTR190(1,2,1)=AAABC(429)
      MTR190(1,2,2)=AAABC(415)
      MTR190(1,3,1)=AAABC(417)
      MTR190(1,3,2)=AAABC(425)
      MTR190(1,4,1)=AAABC(421)
      MTR190(1,4,2)=AAABC(427)
      MTR190(2,1,1)=AAABC(429)
      MTR190(2,1,2)=AAABC(415)
      MTR190(2,2,1)=2D0*AAABC(479)
      MTR190(2,2,2)=2D0*AAABC(477)
      MTR190(2,3,1)=AAABC(473)
      MTR190(2,3,2)=AAABC(481)
      MTR190(2,4,1)=AAABC(483)
      MTR190(2,4,2)=AAABC(475)
      MTR190(3,1,1)=AAABC(417)
      MTR190(3,1,2)=AAABC(425)
      MTR190(3,2,1)=AAABC(473)
      MTR190(3,2,2)=AAABC(481)
      MTR190(3,3,1)=2D0*AAABC(519)
      MTR190(3,3,2)=2D0*AAABC(523)
      MTR190(3,4,1)=AAABC(521)
      MTR190(3,4,2)=AAABC(525)
      MTR190(4,1,1)=AAABC(421)
      MTR190(4,1,2)=AAABC(427)
      MTR190(4,2,1)=AAABC(483)
      MTR190(4,2,2)=AAABC(475)
      MTR190(4,3,1)=AAABC(521)
      MTR190(4,3,2)=AAABC(525)
      MTR190(4,4,1)=2D0*AAABC(555)
      MTR190(4,4,2)=2D0*AAABC(553)
      MTR191(1,1,1)=2D0*AAABC(448)
      MTR191(1,1,2)=2D0*AAABC(452)
      MTR191(1,1,3)=2D0*AAABC(442)
      MTR191(1,2,1)=-AAABC(440)
      MTR191(1,2,2)=-AAABC(436)
      MTR191(1,2,3)=-AAABC(444)
      MTR191(1,3,1)=-AAABC(432)
      MTR191(1,3,2)=-AAABC(446)
      MTR191(1,3,3)=-AAABC(430)
      MTR191(1,4,1)=-AAABC(450)
      MTR191(1,4,2)=-AAABC(434)
      MTR191(1,4,3)=-AAABC(438)
      MTR191(2,1,1)=-AAABC(440)
      MTR191(2,1,2)=-AAABC(436)
      MTR191(2,1,3)=-AAABC(444)
      MTR191(2,2,1)=2D0*AAABC(498)
      MTR191(2,2,2)=2D0*AAABC(496)
      MTR191(2,2,3)=2D0*AAABC(494)
      MTR191(2,3,1)=-AAABC(492)
      MTR191(2,3,2)=-AAABC(486)
      MTR191(2,3,3)=-AAABC(500)
      MTR191(2,4,1)=-AAABC(488)
      MTR191(2,4,2)=-AAABC(490)
      MTR191(2,4,3)=-AAABC(484)
      MTR191(3,1,1)=-AAABC(432)
      MTR191(3,1,2)=-AAABC(446)
      MTR191(3,1,3)=-AAABC(430)
      MTR191(3,2,1)=-AAABC(492)
      MTR191(3,2,2)=-AAABC(486)
      MTR191(3,2,3)=-AAABC(500)
      MTR191(3,3,1)=2D0*AAABC(528)
      MTR191(3,3,2)=2D0*AAABC(526)
      MTR191(3,3,3)=2D0*AAABC(534)
      MTR191(3,4,1)=-AAABC(532)
      MTR191(3,4,2)=-AAABC(536)
      MTR191(3,4,3)=-AAABC(530)
      MTR191(4,1,1)=-AAABC(450)
      MTR191(4,1,2)=-AAABC(434)
      MTR191(4,1,3)=-AAABC(438)
      MTR191(4,2,1)=-AAABC(488)
      MTR191(4,2,2)=-AAABC(490)
      MTR191(4,2,3)=-AAABC(484)
      MTR191(4,3,1)=-AAABC(532)
      MTR191(4,3,2)=-AAABC(536)
      MTR191(4,3,3)=-AAABC(530)
      MTR191(4,4,1)=2D0*AAABC(558)
      MTR191(4,4,2)=2D0*AAABC(556)
      MTR191(4,4,3)=2D0*AAABC(560)
      MTR192(1,1,1)=2D0*AAABC(449)
      MTR192(1,1,2)=2D0*AAABC(453)
      MTR192(1,1,3)=2D0*AAABC(443)
      MTR192(1,2,1)=-AAABC(441)
      MTR192(1,2,2)=-AAABC(437)
      MTR192(1,2,3)=-AAABC(445)
      MTR192(1,3,1)=-AAABC(433)
      MTR192(1,3,2)=-AAABC(447)
      MTR192(1,3,3)=-AAABC(431)
      MTR192(1,4,1)=-AAABC(451)
      MTR192(1,4,2)=-AAABC(435)
      MTR192(1,4,3)=-AAABC(439)
      MTR192(2,1,1)=-AAABC(441)
      MTR192(2,1,2)=-AAABC(437)
      MTR192(2,1,3)=-AAABC(445)
      MTR192(2,2,1)=2D0*AAABC(499)
      MTR192(2,2,2)=2D0*AAABC(497)
      MTR192(2,2,3)=2D0*AAABC(495)
      MTR192(2,3,1)=-AAABC(493)
      MTR192(2,3,2)=-AAABC(487)
      MTR192(2,3,3)=-AAABC(501)
      MTR192(2,4,1)=-AAABC(489)
      MTR192(2,4,2)=-AAABC(491)
      MTR192(2,4,3)=-AAABC(485)
      MTR192(3,1,1)=-AAABC(433)
      MTR192(3,1,2)=-AAABC(447)
      MTR192(3,1,3)=-AAABC(431)
      MTR192(3,2,1)=-AAABC(493)
      MTR192(3,2,2)=-AAABC(487)
      MTR192(3,2,3)=-AAABC(501)
      MTR192(3,3,1)=2D0*AAABC(529)
      MTR192(3,3,2)=2D0*AAABC(527)
      MTR192(3,3,3)=2D0*AAABC(535)
      MTR192(3,4,1)=-AAABC(533)
      MTR192(3,4,2)=-AAABC(537)
      MTR192(3,4,3)=-AAABC(531)
      MTR192(4,1,1)=-AAABC(451)
      MTR192(4,1,2)=-AAABC(435)
      MTR192(4,1,3)=-AAABC(439)
      MTR192(4,2,1)=-AAABC(489)
      MTR192(4,2,2)=-AAABC(491)
      MTR192(4,2,3)=-AAABC(485)
      MTR192(4,3,1)=-AAABC(533)
      MTR192(4,3,2)=-AAABC(537)
      MTR192(4,3,3)=-AAABC(531)
      MTR192(4,4,1)=2D0*AAABC(559)
      MTR192(4,4,2)=2D0*AAABC(557)
      MTR192(4,4,3)=2D0*AAABC(561)
      MTR193(1)=AAABC(454)
      MTR193(2)=AAABC(502)
      MTR193(3)=AAABC(538)
      MTR193(4)=AAABC(562)
      MTR194(1)=AAABC(455)
      MTR194(2)=AAABC(503)
      MTR194(3)=AAABC(539)
      MTR194(4)=AAABC(563)
      MTR195(1,1)=AAABC(456)
      MTR195(1,2)=AAABC(458)
      MTR195(2,1)=AAABC(504)
      MTR195(2,2)=AAABC(506)
      MTR195(3,1)=AAABC(540)
      MTR195(3,2)=AAABC(542)
      MTR195(4,1)=AAABC(566)
      MTR195(4,2)=AAABC(564)
      MTR196(1,1)=AAABC(457)
      MTR196(1,2)=AAABC(459)
      MTR196(2,1)=AAABC(505)
      MTR196(2,2)=AAABC(507)
      MTR196(3,1)=AAABC(541)
      MTR196(3,2)=AAABC(543)
      MTR196(4,1)=AAABC(567)
      MTR196(4,2)=AAABC(565)
      MTR197(1,1)=AAABC(464)
      MTR197(1,2)=AAABC(462)
      MTR197(1,3)=AAABC(460)
      MTR197(2,1)=AAABC(512)
      MTR197(2,2)=AAABC(510)
      MTR197(2,3)=AAABC(508)
      MTR197(3,1)=AAABC(546)
      MTR197(3,2)=AAABC(544)
      MTR197(3,3)=AAABC(548)
      MTR197(4,1)=AAABC(572)
      MTR197(4,2)=AAABC(570)
      MTR197(4,3)=AAABC(568)
      MTR198(1,1)=AAABC(465)
      MTR198(1,2)=AAABC(463)
      MTR198(1,3)=AAABC(461)
      MTR198(2,1)=AAABC(513)
      MTR198(2,2)=AAABC(511)
      MTR198(2,3)=AAABC(509)
      MTR198(3,1)=AAABC(547)
      MTR198(3,2)=AAABC(545)
      MTR198(3,3)=AAABC(549)
      MTR198(4,1)=AAABC(573)
      MTR198(4,2)=AAABC(571)
      MTR198(4,3)=AAABC(569)
      MTR199(1)=AAABC(576)
      MTR199(2)=AAABC(578)
      MTR200(1)=AAABC(577)
      MTR200(2)=AAABC(579)
      MTR201(1)=AAABC(582)
      MTR201(2)=AAABC(584)
      MTR201(3)=AAABC(580)
      MTR202(1)=AAABC(583)
      MTR202(2)=AAABC(585)
      MTR202(3)=AAABC(581)
      MTR203(1,1)=AAABR(495)
      MTR203(1,2)=AAABR(493)
      MTR203(2,1)=AAABR(493)
      MTR203(2,2)=AAABR(497)
      MTR204(1,1)=AAABR(496)
      MTR204(1,2)=AAABR(494)
      MTR204(2,1)=AAABR(494)
      MTR204(2,2)=AAABR(498)
      MTR205(1,1)=AAABC(600)
      MTR205(1,2)=AAABC(592)
      MTR205(1,3)=AAABC(598)
      MTR205(1,4)=AAABC(596)
      MTR205(2,1)=AAABC(590)
      MTR205(2,2)=AAABC(594)
      MTR205(2,3)=AAABC(588)
      MTR205(2,4)=AAABC(586)
      MTR206(1,1)=AAABC(601)
      MTR206(1,2)=AAABC(593)
      MTR206(1,3)=AAABC(599)
      MTR206(1,4)=AAABC(597)
      MTR206(2,1)=AAABC(591)
      MTR206(2,2)=AAABC(595)
      MTR206(2,3)=AAABC(589)
      MTR206(2,4)=AAABC(587)
      MTR207(1)=AAABC(602)
      MTR207(2)=AAABC(604)
      MTR208(1)=AAABC(603)
      MTR208(2)=AAABC(605)
      MTR209(1,1)=AAABR(504)
      MTR209(1,2)=AAABR(502)
      MTR209(2,1)=AAABR(505)
      MTR209(2,2)=AAABR(503)
      MTR209(3,3)=AAABR(500)
      MTR210(1,1)=AAABR(504)
      MTR210(1,2)=AAABR(505)
      MTR210(2,1)=AAABR(502)
      MTR210(2,2)=AAABR(503)
      MTR210(3,3)=AAABR(500)
      MTR211(1,1)=AAABC(606)
      MTR211(1,2)=AAABC(610)
      MTR211(1,3)=AAABC(608)
      MTR211(1,4)=AAABC(612)
      MTR211(2,1)=AAABC(622)
      MTR211(2,2)=AAABC(618)
      MTR211(2,3)=AAABC(620)
      MTR211(2,4)=AAABC(616)
      MTR212(1,1)=AAABC(607)
      MTR212(1,2)=AAABC(611)
      MTR212(1,3)=AAABC(609)
      MTR212(1,4)=AAABC(613)
      MTR212(2,1)=AAABC(623)
      MTR212(2,2)=AAABC(619)
      MTR212(2,3)=AAABC(621)
      MTR212(2,4)=AAABC(617)
      MTR213(1)=AAABC(614)
      MTR213(2)=AAABC(624)
      MTR214(1)=AAABC(615)
      MTR214(2)=AAABC(625)
      MTR215(1,1)=AAABC(629)
      MTR215(1,2)=dconjg(AAABC(626))
      MTR215(1,3)=dconjg(AAABC(627))
      MTR215(1,4)=dconjg(AAABC(628))
      MTR215(2,1)=AAABC(626)
      MTR215(2,2)=AAABC(631)
      MTR215(2,3)=dconjg(AAABC(633))
      MTR215(2,4)=dconjg(AAABC(632))
      MTR215(3,1)=AAABC(627)
      MTR215(3,2)=AAABC(633)
      MTR215(3,3)=AAABC(636)
      MTR215(3,4)=dconjg(AAABC(635))
      MTR215(4,1)=AAABC(628)
      MTR215(4,2)=AAABC(632)
      MTR215(4,3)=AAABC(635)
      MTR215(4,4)=AAABC(638)
      MTR216(1,1)=AAABC(629)
      MTR216(1,2)=AAABC(626)
      MTR216(1,3)=AAABC(627)
      MTR216(1,4)=AAABC(628)
      MTR216(2,1)=dconjg(AAABC(626))
      MTR216(2,2)=AAABC(631)
      MTR216(2,3)=AAABC(633)
      MTR216(2,4)=AAABC(632)
      MTR216(3,1)=dconjg(AAABC(627))
      MTR216(3,2)=dconjg(AAABC(633))
      MTR216(3,3)=AAABC(636)
      MTR216(3,4)=AAABC(635)
      MTR216(4,1)=dconjg(AAABC(628))
      MTR216(4,2)=dconjg(AAABC(632))
      MTR216(4,3)=dconjg(AAABC(635))
      MTR216(4,4)=AAABC(638)
      MTR217(1)=AAABC(630)
      MTR217(2)=AAABC(634)
      MTR217(3)=AAABC(637)
      MTR217(4)=AAABC(639)
      MTR218(1)=dconjg(AAABC(630))
      MTR218(2)=dconjg(AAABC(634))
      MTR218(3)=dconjg(AAABC(637))
      MTR218(4)=dconjg(AAABC(639))
      MTR219(1)=AAABR(516)
      MTR219(2)=AAABR(514)
      MTR220(1)=AAABR(515)
      MTR220(2)=AAABR(513)
      MTR221(1)=AAABR(520)
      MTR221(2)=AAABR(518)
      MTR221(3)=AAABR(522)
      MTR222(1)=AAABR(519)
      MTR222(2)=AAABR(517)
      MTR222(3)=AAABR(521)
      MTR223(1)=AAABR(523)
      MTR223(2)=AAABR(525)
      MTR223(3)=AAABR(524)
      MTR224(1)=AAABR(531)
      MTR224(2)=AAABR(530)
      MTR224(3)=AAABR(529)
      MTR225(1)=AAABR(528)
      MTR225(2)=AAABR(526)
      MTR225(3)=AAABR(527)
      MTR226(1)=AAABR(532)
      MTR226(2)=AAABR(534)
      MTR226(3)=AAABR(533)
      MTR227(1)=AAABR(538)
      MTR227(2)=AAABR(539)
      MTR227(3)=AAABR(540)
      MTR228(1)=AAABR(537)
      MTR228(2)=AAABR(535)
      MTR228(3)=AAABR(536)
      MTR229(1)=AAABR(546)
      MTR229(2)=AAABR(542)
      MTR229(3)=AAABR(544)
      MTR230(1)=AAABR(545)
      MTR230(2)=AAABR(541)
      MTR230(3)=AAABR(543)
      MTR231(1)=AAABR(549)
      MTR231(2)=AAABR(547)
      MTR231(3)=AAABR(548)
      MTR232(1)=AAABR(554)
      MTR232(2)=AAABR(555)
      MTR232(3)=AAABR(553)
      MTR233(1)=AAABR(551)
      MTR233(2)=AAABR(552)
      MTR233(3)=AAABR(550)
      MTR234(1)=AAABR(557)
      MTR234(2)=AAABR(559)
      MTR234(3)=AAABR(561)
      MTR235(1)=AAABR(556)
      MTR235(2)=AAABR(558)
      MTR235(3)=AAABR(560)
      MTR236(1)=AAABR(563)
      MTR236(2)=AAABR(565)
      MTR236(3)=AAABR(564)
      MTR237(1)=AAABR(568)
      MTR237(2)=AAABR(569)
      MTR237(3)=AAABR(566)
      MTR238(1)=AAABR(576)
      MTR238(2)=AAABR(575)
      MTR238(3)=AAABR(574)
      MTR239(1)=AAABR(571)
      MTR239(2)=AAABR(570)
      MTR239(3)=AAABR(572)
      MTR240(1)=AAABR(578)
      MTR240(2)=AAABR(579)
      MTR240(3)=AAABR(577)
      MTR241(1)=AAABR(582)
      MTR241(2)=AAABR(581)
      MTR241(3)=AAABR(580)
      MTR242(1)=AAABR(587)
      MTR242(2)=AAABR(588)
      MTR242(3)=AAABR(586)
      MTR243(1)=AAABR(585)
      MTR243(2)=AAABR(584)
      MTR243(3)=AAABR(583)
      MTR244(1,1)=AAABR(590)
      MTR244(1,2)=AAABR(596)
      MTR244(2,1)=AAABR(592)
      MTR244(2,2)=AAABR(594)
      MTR244(3,3)=AAABR(598)
      MTR245(1,1)=AAABR(603)
      MTR245(1,2)=AAABR(600)
      MTR245(2,1)=AAABR(601)
      MTR245(2,2)=AAABR(599)
      MTR245(3,3)=AAABR(602)
      MTR246(1,1)=AAABR(589)
      MTR246(1,2)=AAABR(595)
      MTR246(2,1)=AAABR(591)
      MTR246(2,2)=AAABR(593)
      MTR246(3,3)=AAABR(597)
      MTR247(1,1)=AAABR(605)
      MTR247(1,2)=AAABR(607)
      MTR247(2,1)=AAABR(611)
      MTR247(2,2)=AAABR(613)
      MTR247(3,3)=AAABR(609)
      MTR248(1,1)=AAABR(616)
      MTR248(1,2)=AAABR(614)
      MTR248(2,1)=AAABR(615)
      MTR248(2,2)=AAABR(618)
      MTR248(3,3)=AAABR(617)
      MTR249(1,1)=AAABR(604)
      MTR249(1,2)=AAABR(606)
      MTR249(2,1)=AAABR(610)
      MTR249(2,2)=AAABR(612)
      MTR249(3,3)=AAABR(608)
      MTR250(1)=AAABR(619)
      MTR250(2)=AAABR(621)
      MTR250(3)=AAABR(620)
      MTR251(1)=AAABR(622)
      MTR251(2)=AAABR(623)
      MTR251(3)=AAABR(624)
      MTR252(1)=AAABR(629)
      MTR252(2)=AAABR(628)
      MTR252(3)=AAABR(630)
      MTR253(1)=AAABR(625)
      MTR253(2)=AAABR(627)
      MTR253(3)=AAABR(626)
      MTR254(1,1)=AAABR(634)
      MTR254(1,2)=AAABR(638)
      MTR254(2,1)=AAABR(640)
      MTR254(2,2)=AAABR(636)
      MTR254(3,3)=AAABR(632)
      MTR255(1,1)=AAABR(645)
      MTR255(1,2)=AAABR(644)
      MTR255(2,1)=AAABR(642)
      MTR255(2,2)=AAABR(643)
      MTR255(3,3)=AAABR(641)
      MTR256(1,1)=AAABR(633)
      MTR256(1,2)=AAABR(637)
      MTR256(2,1)=AAABR(639)
      MTR256(2,2)=AAABR(635)
      MTR256(3,3)=AAABR(631)
      MTR257(1,1)=AAABR(655)
      MTR257(1,2)=AAABR(647)
      MTR257(2,1)=AAABR(649)
      MTR257(2,2)=AAABR(653)
      MTR257(3,3)=AAABR(651)
      MTR258(1,1)=AAABR(658)
      MTR258(1,2)=AAABR(657)
      MTR258(2,1)=AAABR(656)
      MTR258(2,2)=AAABR(659)
      MTR258(3,3)=AAABR(660)
      MTR259(1,1)=AAABR(654)
      MTR259(1,2)=AAABR(646)
      MTR259(2,1)=AAABR(648)
      MTR259(2,2)=AAABR(652)
      MTR259(3,3)=AAABR(650)
      MTR260(1,1)=AAABR(590)
      MTR260(1,2)=AAABR(592)
      MTR260(2,1)=AAABR(596)
      MTR260(2,2)=AAABR(594)
      MTR260(3,3)=AAABR(598)
      MTR261(1,1)=AAABR(603)
      MTR261(1,2)=AAABR(601)
      MTR261(2,1)=AAABR(600)
      MTR261(2,2)=AAABR(599)
      MTR261(3,3)=AAABR(602)
      MTR262(1,1)=AAABR(589)
      MTR262(1,2)=AAABR(591)
      MTR262(2,1)=AAABR(595)
      MTR262(2,2)=AAABR(593)
      MTR262(3,3)=AAABR(597)
      MTR263(1,1)=AAABR(634)
      MTR263(1,2)=AAABR(640)
      MTR263(2,1)=AAABR(638)
      MTR263(2,2)=AAABR(636)
      MTR263(3,3)=AAABR(632)
      MTR264(1,1)=AAABR(645)
      MTR264(1,2)=AAABR(642)
      MTR264(2,1)=AAABR(644)
      MTR264(2,2)=AAABR(643)
      MTR264(3,3)=AAABR(641)
      MTR265(1,1)=AAABR(633)
      MTR265(1,2)=AAABR(639)
      MTR265(2,1)=AAABR(637)
      MTR265(2,2)=AAABR(635)
      MTR265(3,3)=AAABR(631)
      MTR266(1)=AAABR(661)
      MTR266(2)=AAABR(663)
      MTR266(3)=AAABR(662)
      MTR267(1)=AAABR(664)
      MTR267(2)=AAABR(666)
      MTR267(3)=AAABR(665)
      MTR268(1)=AAABR(671)
      MTR268(2)=AAABR(672)
      MTR268(3)=AAABR(670)
      MTR269(1)=AAABR(668)
      MTR269(2)=AAABR(667)
      MTR269(3)=AAABR(669)
      MTR270(1)=AAABR(674)
      MTR270(2)=AAABR(675)
      MTR270(3)=AAABR(673)
      MTR271(1)=AAABR(678)
      MTR271(2)=AAABR(677)
      MTR271(3)=AAABR(676)
      MTR272(1)=AAABR(682)
      MTR272(2)=AAABR(683)
      MTR272(3)=AAABR(684)
      MTR273(1)=AAABR(680)
      MTR273(2)=AAABR(679)
      MTR273(3)=AAABR(681)
      MTR274(1,1)=AAABR(605)
      MTR274(1,2)=AAABR(611)
      MTR274(2,1)=AAABR(607)
      MTR274(2,2)=AAABR(613)
      MTR274(3,3)=AAABR(609)
      MTR275(1,1)=AAABR(616)
      MTR275(1,2)=AAABR(615)
      MTR275(2,1)=AAABR(614)
      MTR275(2,2)=AAABR(618)
      MTR275(3,3)=AAABR(617)
      MTR276(1,1)=AAABR(604)
      MTR276(1,2)=AAABR(610)
      MTR276(2,1)=AAABR(606)
      MTR276(2,2)=AAABR(612)
      MTR276(3,3)=AAABR(608)
      MTR277(1,1)=AAABR(655)
      MTR277(1,2)=AAABR(649)
      MTR277(2,1)=AAABR(647)
      MTR277(2,2)=AAABR(653)
      MTR277(3,3)=AAABR(651)
      MTR278(1,1)=AAABR(658)
      MTR278(1,2)=AAABR(656)
      MTR278(2,1)=AAABR(657)
      MTR278(2,2)=AAABR(659)
      MTR278(3,3)=AAABR(660)
      MTR279(1,1)=AAABR(654)
      MTR279(1,2)=AAABR(648)
      MTR279(2,1)=AAABR(646)
      MTR279(2,2)=AAABR(652)
      MTR279(3,3)=AAABR(650)
      MTR280(1)=AAABR(686)
      MTR280(2)=AAABR(685)
      MTR280(3)=AAABR(687)
      MTR281(1)=AAABR(690)
      MTR281(2)=AAABR(688)
      MTR281(3)=AAABR(689)
      MTR282(1)=AAABR(694)
      MTR282(2)=AAABR(696)
      MTR282(3)=AAABR(695)
      MTR283(1)=AAABR(692)
      MTR283(2)=AAABR(693)
      MTR283(3)=AAABR(691)
      MTR284(1)=AAABR(700)
      MTR284(2)=AAABR(702)
      MTR285(1)=AAABR(699)
      MTR285(2)=AAABR(701)
      MTR286(1)=AAABR(706)
      MTR286(2)=AAABR(708)
      MTR286(3)=AAABR(704)
      MTR287(1)=AAABR(705)
      MTR287(2)=AAABR(707)
      MTR287(3)=AAABR(703)
      MTR288(1)=AAABR(714)
      MTR288(2)=AAABR(713)
      MTR289(1)=AAABR(711)
      MTR289(2)=AAABR(712)
      MTR290(1,1)=AAABR(718)
      MTR290(1,2)=AAABR(719)
      MTR290(2,1)=MTR290(1,2)
      MTR290(2,2)=AAABR(720)
      MTR291(1,1)=AAABR(715)
      MTR291(1,2)=AAABR(716)
      MTR291(2,1)=MTR291(1,2)
      MTR291(2,2)=AAABR(717)
      MTR292(1,1)=AAABR(730)
      MTR292(1,2)=AAABR(729)
      MTR292(1,3)=AAABR(728)
      MTR292(2,1)=MTR292(1,2)
      MTR292(2,2)=AAABR(727)
      MTR292(2,3)=AAABR(732)
      MTR292(3,1)=MTR292(1,3)
      MTR292(3,2)=MTR292(2,3)
      MTR292(3,3)=AAABR(731)
      MTR293(1,1)=AAABR(726)
      MTR293(1,2)=AAABR(722)
      MTR293(1,3)=AAABR(724)
      MTR293(2,1)=MTR293(1,2)
      MTR293(2,2)=AAABR(725)
      MTR293(2,3)=AAABR(721)
      MTR293(3,1)=MTR293(1,3)
      MTR293(3,2)=MTR293(2,3)
      MTR293(3,3)=AAABR(723)

      end

***********************************************

      subroutine ModelVarIni(sqrtS, *)
      implicit none
      double precision sqrtS
      double precision Alfas

#include "model.h"

c      double precision ALPHAS2
c      external ALPHAS2

c      Alfas = ALPHAS2(sqrtS)
c      GG = sqrt(4*pi*Alfas)
      end

************************************************

      subroutine ModelDigest
      implicit none

#include "model.h"


      end

#include "neutd5.F"

