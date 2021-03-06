* gzvert_his5_sld5.F
* form factors for {MNE1, MNE1} -> {0, MZ}
* generated by FormCalc 12 Apr 2013 17:41

	subroutine gzvert_his5_sld5
	implicit none

#include "gzvars.h"

	double complex dup1568, dup1569

        dup1568 = MTR026(his5,sld5)*MTR061(sld5) + 
     -   MTR022(his5,sld5)*MTR063(sld5)

        dup1569 = MTR026(his5,sld5)*MTR200(sld5) + 
     -   MTR022(his5,sld5)*MTR205(sld5)

        Cloop(1) = Cloop(1) + 
     -   (B0(S,sldMass(sld5)**2,sleMass(sld5)**2)*
     -       (1/(64.D0*Pi**2)*(Abb5*dup1569*MTR183(1,1,his5)) - 
     -         1/(64.D0*Pi**2)*(Abb7*dup1569*MTR184(1,1,his5))) + 
     -      (Cval(cc00,iint74(sld5)) + Cval(cc00,iint75(sld5)))*
     -       (-(1/(32.D0*Pi**2)*
     -            (Abb5*dup1568*EE*MTR183(1,1,his5))) + 
     -         1/(32.D0*Pi**2)*(Abb7*dup1568*EE*MTR184(1,1,his5))) + 
     -      (Cval(cc12,iint74(sld5)) + Cval(cc12,iint75(sld5)))*
     -       (1/(32.D0*Pi**2)*
     -          (AbbSum30*dup1568*EE*MTR183(1,1,his5)) - 
     -         1/(32.D0*Pi**2)*(AbbSum31*dup1568*EE*MTR184(1,1,his5))
     -         ) + Cval(cc00,iint73(sld5))*
     -       (-(1/(16.D0*Pi**2)*
     -            (Abb5*EE*MTR021(his5,sld5)*MTR060(sld5)*
     -              MTR183(1,1,his5))) + 
     -         1/(16.D0*Pi**2)*
     -          (Abb7*EE*MTR021(his5,sld5)*MTR060(sld5)*
     -            MTR184(1,1,his5))) + 
     -      Cval(cc12,iint73(sld5))*
     -       (1/(16.D0*Pi**2)*
     -          (AbbSum30*EE*MTR021(his5,sld5)*MTR060(sld5)*
     -            MTR183(1,1,his5)) - 
     -         1/(16.D0*Pi**2)*
     -          (AbbSum31*EE*MTR021(his5,sld5)*MTR060(sld5)*
     -            MTR184(1,1,his5))) + 
     -      B0(S,sldMass(sld5)**2,sldMass(sld5)**2)*
     -       (1/(64.D0*Pi**2)*
     -          (Abb5*MTR021(his5,sld5)*MTR183(1,1,his5)*
     -            MTR197(sld5)) - 
     -         1/(64.D0*Pi**2)*
     -          (Abb7*MTR021(his5,sld5)*MTR184(1,1,his5)*
     -            MTR197(sld5))))/(S - hisMass(his5)**2+(0D0,1D0)*hisMass(his5)*hisW(his5))

#ifdef DEBUG
	DEB 'vert_his5_sld5: Cloop(1) =', Cloop(1)
#endif

	end

