* se_HHG0.F
* this file is part of FeynHiggs
* generated 30 May 2012 12:16

#include "externals.h"
#include "types.h"
#include "debug.h"


	ComplexType function se_HHG0(k2, bX)
	implicit none
	RealType k2
	integer bX


#include "FH.h"

	ComplexType se_HHG0_mfv, se_HHG0_nmfv
	external se_HHG0_mfv, se_HHG0_nmfv

	if( ibits(fv, 3, 2) .eq. 0 ) then
	  se_HHG0 = se_HHG0_mfv(k2, bX)
	else
	  se_HHG0 = se_HHG0_nmfv(k2, bX)
	endif

#ifdef DEBUG
	DHIGGS 'k2 =', k2 ENDL
	DHIGGS 'se_HHG0 =', se_HHG0 ENDL
#endif
	end


************************************************************************


	ComplexType function se_HHG0_mfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup16, dup17, dup2, dup3, dup4, dup5, dup6, dup7
	ComplexType dup8, dup9, tmp1, tmp2

	se_HHG0_mfv = 0

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_HHG0_mfv = se_HHG0_mfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,tM,Gen3),MSf2(Sfe4,tM,Gen3))*
     &        (USf(Sfe4,1,tM,Gen3)*
     &           ((CAB*MW*MZ*SB*(3 - 4*SW2) + 
     &                6*CW*SA*Mf2(tM,Gen3))*USfC(Sfe3,1,tM,Gen3) + 
     &             3*CW*(SA*KfC(Gen3,Gen3,tM) - 
     &                CA*MUE*Mf(tM,Gen3))*USfC(Sfe3,2,tM,Gen3)) + 
     &          USf(Sfe4,2,tM,Gen3)*
     &           (3*CW*(SA*Kf(Gen3,Gen3,tM) - CA*MUEC*Mf(tM,Gen3))*
     &              USfC(Sfe3,1,tM,Gen3) + 
     &             2*(2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM,Gen3))*
     &              USfC(Sfe3,2,tM,Gen3)))*
     &        ((TB*Kf(Gen3,Gen3,tM) - MUEC*Mf(tM,Gen3))*
     &           USf(Sfe3,2,tM,Gen3)*USfC(Sfe4,1,tM,Gen3) + 
     &          (-(TB*KfC(Gen3,Gen3,tM)) + MUE*Mf(tM,Gen3))*
     &           USf(Sfe3,1,tM,Gen3)*USfC(Sfe4,2,tM,Gen3)))/
     &      (CW*MW2*SB*SW2*TB)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_mfv =", se_HHG0_mfv ENDL
#endif

	if( mssmpart .eq. 1 ) return

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        se_HHG0_mfv = se_HHG0_mfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3))*
     &        (USf(Sfe4,1,bX,Gen3)*
     &           ((CAB*CB*MW*MZ*(-3 + 2*SW2) + 
     &                6*CA*CW*Mf2(bX,Gen3))*USfC(Sfe3,1,bX,Gen3) + 
     &             3*CW*(CA*KfC(Gen3,Gen3,bX) - 
     &                MUE*SA*Mf(bX,Gen3))*USfC(Sfe3,2,bX,Gen3)) + 
     &          USf(Sfe4,2,bX,Gen3)*
     &           (3*CW*(CA*Kf(Gen3,Gen3,bX) - MUEC*SA*Mf(bX,Gen3))*
     &              USfC(Sfe3,1,bX,Gen3) - 
     &             2*(CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Gen3))*
     &              USfC(Sfe3,2,bX,Gen3)))*
     &        ((-Kf(Gen3,Gen3,bX) + MUEC*TB*Mf(bX,Gen3))*
     &           USf(Sfe3,2,bX,Gen3)*USfC(Sfe4,1,bX,Gen3) + 
     &          (KfC(Gen3,Gen3,bX) - MUE*TB*Mf(bX,Gen3))*
     &           USf(Sfe3,1,bX,Gen3)*USfC(Sfe4,2,bX,Gen3)))/
     &      (CB*CW*MW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_mfv =", se_HHG0_mfv ENDL
#endif

	if( mssmpart .eq. 2 ) return

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_HHG0_mfv = se_HHG0_mfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (USf(Sfe4,1,2,Gen3)*
     &           ((CAB*CB*MW*MZ*(-1 + 2*SW2) + 
     &                2*CA*CW*Mf2(2,Gen3))*USfC(Sfe3,1,2,Gen3) + 
     &             CW*(CA*KfC(Gen3,Gen3,2) - MUE*SA*Mf(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)) + 
     &          USf(Sfe4,2,2,Gen3)*
     &           (CW*(CA*Kf(Gen3,Gen3,2) - MUEC*SA*Mf(2,Gen3))*
     &              USfC(Sfe3,1,2,Gen3) + 
     &             2*(-(CAB*CB*MW*MZ*SW2) + CA*CW*Mf2(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)))*
     &        ((-Kf(Gen3,Gen3,2) + MUEC*TB*Mf(2,Gen3))*
     &           USf(Sfe3,2,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &          (KfC(Gen3,Gen3,2) - MUE*TB*Mf(2,Gen3))*
     &           USf(Sfe3,1,2,Gen3)*USfC(Sfe4,2,2,Gen3)))/
     &      (CB*CW*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_mfv =", se_HHG0_mfv ENDL
#endif

	if( mssmpart .eq. 3 ) return

	LOOP(Cha4, 1,2,1)

	tmp1 = A0(MCha2(Cha4))

	LOOP(Cha3, 1,2,1)

        dup1 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    SA*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = -(CB*UCha(Cha3,2)*VCha(Cha4,1)) + 
     &    SB*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = CB*UChaC(Cha4,2)*VChaC(Cha3,1) - 
     &    SB*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    SA*UChaC(Cha3,1)*VChaC(Cha4,2)

        se_HHG0_mfv = se_HHG0_mfv - 
     &    cI/(4.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (tmp1 + k2*B1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          B0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup3*dup4*MCha(Cha3)*MCha(Cha4) + 
     &             dup2*dup4*MCha2(Cha3) + 
     &             dup1*(dup2*MCha(Cha3)*MCha(Cha4) + 
     &                dup3*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)

	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)

	tmp2 = A0(MNeu2(Neu4))

	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = CB*ZNeu(Neu3,3) + SB*ZNeu(Neu3,4)

	dup7 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup8 = -(SW*ZNeu(Neu4,1)) + CW*ZNeu(Neu4,2)

	dup9 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup10 = CB*ZNeu(Neu4,3) + SB*ZNeu(Neu4,4)

	dup11 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup12 = CB*ZNeuC(Neu3,3) + SB*ZNeuC(Neu3,4)

	dup13 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup14 = -(SW*ZNeuC(Neu4,1)) + CW*ZNeuC(Neu4,2)

	dup15 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

	dup16 = CB*ZNeuC(Neu4,3) + SB*ZNeuC(Neu4,4)

        dup17 = (dup12*dup13 + dup11*dup16)*
     &     (dup5*dup9 + CA*dup7*ZNeu(Neu3,3) + 
     &       dup8*SA*ZNeu(Neu3,4)) - 
     &    (dup10*dup5 + dup6*dup7)*
     &     (dup11*dup15 + CA*dup13*ZNeuC(Neu3,3) + 
     &       dup14*SA*ZNeuC(Neu3,4))

        se_HHG0_mfv = se_HHG0_mfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (dup17*(tmp2 + k2*B1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          B0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup17*MNeu2(Neu3) - 
     &             MNeu(Neu3)*MNeu(Neu4)*
     &              ((dup10*dup5 + dup6*dup7)*
     &                 (dup5*dup9 + CA*dup7*ZNeu(Neu3,3) + 
     &                   dup8*SA*ZNeu(Neu3,4)) - 
     &                (dup12*dup13 + dup11*dup16)*
     &                 (dup11*dup15 + CA*dup13*ZNeuC(Neu3,3) + 
     &                   dup14*SA*ZNeuC(Neu3,4))))))/(CW2*SW2)

	ENDLOOP(Neu3)

	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_mfv =", se_HHG0_mfv ENDL
#endif

	end


************************************************************************


	ComplexType function se_HHG0_nmfv(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer All3, All4, Cha3, Cha4, Gen3, Ind1, Ind2, Ind3, Ind4
	integer Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup16, dup17, dup2, dup3, dup4, dup5, dup6, dup7
	ComplexType dup8, dup9, tmp1, tmp2, tmp3, tmp4

	se_HHG0_nmfv = 0

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp1 = B0(k2,MASf2(All3,tM),MASf2(All4,tM))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HHG0_nmfv = se_HHG0_nmfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*tmp1*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*MW*MZ*SB*(3 - 4*SW2) + 
     &                6*CW*SA*Mf2(tM,Ind3))*UASf(All4,Ind4,tM) + 
     &             3*CW*(SA*Kf(Ind3,Ind4,tM) - 
     &                CA*MUEC*Delta(Ind3,Ind4)*Mf(tM,Ind3))*
     &              UASf(All4,3 + Ind4,tM))*UASfC(All3,Ind3,tM) + 
     &          (3*CW*(SA*KfC(Ind4,Ind3,tM) - 
     &                CA*MUE*Delta(Ind3,Ind4)*Mf(tM,Ind3))*
     &              UASf(All4,Ind4,tM) + 
     &             2*Delta(Ind3,Ind4)*
     &              (2*CAB*MW*MZ*SB*SW2 + 3*CW*SA*Mf2(tM,Ind3))*
     &              UASf(All4,3 + Ind4,tM))*UASfC(All3,3 + Ind3,tM)
     &          )*((TB*Kf(Ind1,Ind2,tM) - 
     &             MUEC*Delta(Ind1,Ind2)*Mf(tM,Ind1))*
     &           UASf(All3,3 + Ind2,tM)*UASfC(All4,Ind1,tM) + 
     &          (-(TB*KfC(Ind2,Ind1,tM)) + 
     &             MUE*Delta(Ind1,Ind2)*Mf(tM,Ind1))*
     &           UASf(All3,Ind2,tM)*UASfC(All4,3 + Ind1,tM)))/
     &      (CW*MW2*SB*SW2*TB)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_nmfv =", se_HHG0_nmfv ENDL
#endif

	if( mssmpart .eq. 1 ) return

	LOOP(All4, 1,6,gM)
	LOOP(All3, 1,6,gM)

	tmp2 = B0(k2,MASf2(All3,bX),MASf2(All4,bX))

	LOOP(Ind4, 1,3,1)
	LOOP(Ind3, 1,3,1)
	LOOP(Ind2, 1,3,1)
	LOOP(Ind1, 1,3,1)

        se_HHG0_nmfv = se_HHG0_nmfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*tmp2*
     &        ((Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*(-3 + 2*SW2) + 
     &                6*CA*CW*Mf2(bX,Ind3))*UASf(All4,Ind4,bX) + 
     &             3*CW*(CA*Kf(Ind3,Ind4,bX) - 
     &                MUEC*SA*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,Ind3,bX) + 
     &          (3*CW*(CA*KfC(Ind4,Ind3,bX) - 
     &                MUE*SA*Delta(Ind3,Ind4)*Mf(bX,Ind3))*
     &              UASf(All4,Ind4,bX) - 
     &             2*Delta(Ind3,Ind4)*
     &              (CAB*CB*MW*MZ*SW2 - 3*CA*CW*Mf2(bX,Ind3))*
     &              UASf(All4,3 + Ind4,bX))*UASfC(All3,3 + Ind3,bX)
     &          )*((-Kf(Ind1,Ind2,bX) + 
     &             MUEC*TB*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &           UASf(All3,3 + Ind2,bX)*UASfC(All4,Ind1,bX) + 
     &          (KfC(Ind2,Ind1,bX) - 
     &             MUE*TB*Delta(Ind1,Ind2)*Mf(bX,Ind1))*
     &           UASf(All3,Ind2,bX)*UASfC(All4,3 + Ind1,bX)))/
     &      (CB*CW*MW2*SW2)

	ENDLOOP(Ind1)
	ENDLOOP(Ind2)
	ENDLOOP(Ind3)
	ENDLOOP(Ind4)

	ENDLOOP(All3)
	ENDLOOP(All4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_nmfv =", se_HHG0_nmfv ENDL
#endif

	if( mssmpart .eq. 2 ) return

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        se_HHG0_nmfv = se_HHG0_nmfv + 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        B0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (USf(Sfe4,1,2,Gen3)*
     &           ((CAB*CB*MW*MZ*(-1 + 2*SW2) + 
     &                2*CA*CW*Mf2(2,Gen3))*USfC(Sfe3,1,2,Gen3) + 
     &             CW*(CA*KfC(Gen3,Gen3,2) - MUE*SA*Mf(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)) + 
     &          USf(Sfe4,2,2,Gen3)*
     &           (CW*(CA*Kf(Gen3,Gen3,2) - MUEC*SA*Mf(2,Gen3))*
     &              USfC(Sfe3,1,2,Gen3) + 
     &             2*(-(CAB*CB*MW*MZ*SW2) + CA*CW*Mf2(2,Gen3))*
     &              USfC(Sfe3,2,2,Gen3)))*
     &        ((-Kf(Gen3,Gen3,2) + MUEC*TB*Mf(2,Gen3))*
     &           USf(Sfe3,2,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &          (KfC(Gen3,Gen3,2) - MUE*TB*Mf(2,Gen3))*
     &           USf(Sfe3,1,2,Gen3)*USfC(Sfe4,2,2,Gen3)))/
     &      (CB*CW*MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_nmfv =", se_HHG0_nmfv ENDL
#endif

	if( mssmpart .eq. 3 ) return

	LOOP(Cha4, 1,2,1)

	tmp3 = A0(MCha2(Cha4))

	LOOP(Cha3, 1,2,1)

        dup1 = CA*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    SA*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = -(CB*UCha(Cha3,2)*VCha(Cha4,1)) + 
     &    SB*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = CB*UChaC(Cha4,2)*VChaC(Cha3,1) - 
     &    SB*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = CA*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    SA*UChaC(Cha3,1)*VChaC(Cha4,2)

        se_HHG0_nmfv = se_HHG0_nmfv - 
     &    cI/(4.D0*Pi)*(Alfa1L*
     &        ((dup1*dup3 + dup2*dup4)*
     &           (tmp3 + k2*B1(k2,MCha2(Cha3),MCha2(Cha4))) + 
     &          B0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (dup3*dup4*MCha(Cha3)*MCha(Cha4) + 
     &             dup2*dup4*MCha2(Cha3) + 
     &             dup1*(dup2*MCha(Cha3)*MCha(Cha4) + 
     &                dup3*MCha2(Cha3)))))/SW2

	ENDLOOP(Cha3)

	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)

	tmp4 = A0(MNeu2(Neu4))

	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = CB*ZNeu(Neu3,3) + SB*ZNeu(Neu3,4)

	dup7 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup8 = -(SW*ZNeu(Neu4,1)) + CW*ZNeu(Neu4,2)

	dup9 = CA*ZNeu(Neu4,3) - SA*ZNeu(Neu4,4)

	dup10 = CB*ZNeu(Neu4,3) + SB*ZNeu(Neu4,4)

	dup11 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup12 = CB*ZNeuC(Neu3,3) + SB*ZNeuC(Neu3,4)

	dup13 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup14 = -(SW*ZNeuC(Neu4,1)) + CW*ZNeuC(Neu4,2)

	dup15 = CA*ZNeuC(Neu4,3) - SA*ZNeuC(Neu4,4)

	dup16 = CB*ZNeuC(Neu4,3) + SB*ZNeuC(Neu4,4)

        dup17 = (dup12*dup13 + dup11*dup16)*
     &     (dup5*dup9 + CA*dup7*ZNeu(Neu3,3) + 
     &       dup8*SA*ZNeu(Neu3,4)) - 
     &    (dup10*dup5 + dup6*dup7)*
     &     (dup11*dup15 + CA*dup13*ZNeuC(Neu3,3) + 
     &       dup14*SA*ZNeuC(Neu3,4))

        se_HHG0_nmfv = se_HHG0_nmfv - 
     &    cI/(16.D0*Pi)*(Alfa1L*
     &        (dup17*(tmp4 + k2*B1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          B0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup17*MNeu2(Neu3) - 
     &             MNeu(Neu3)*MNeu(Neu4)*
     &              ((dup10*dup5 + dup6*dup7)*
     &                 (dup5*dup9 + CA*dup7*ZNeu(Neu3,3) + 
     &                   dup8*SA*ZNeu(Neu3,4)) - 
     &                (dup12*dup13 + dup11*dup16)*
     &                 (dup11*dup15 + CA*dup13*ZNeuC(Neu3,3) + 
     &                   dup14*SA*ZNeuC(Neu3,4))))))/(CW2*SW2)

	ENDLOOP(Neu3)

	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "se_HHG0_nmfv =", se_HHG0_nmfv ENDL
#endif

	end

