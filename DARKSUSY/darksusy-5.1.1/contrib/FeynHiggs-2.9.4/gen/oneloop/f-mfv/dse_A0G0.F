* dse_A0G0.F
* generated 30-May-2012 11:58
* this file is part of FeynHiggs
* please do not edit directly

#include "externals.h"
#include "types.h"


	ComplexType function dse_A0G0(k2,bX)
	implicit none
	integer bX
	RealType k2

#include "FH.h"
#include "looptools.h"

	integer Cha3, Cha4, Gen3, Neu3, Neu4, Sfe3, Sfe4

	ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
	ComplexType dup16, dup17, dup2, dup3, dup4, dup5, dup6, dup7
	ComplexType dup8, dup9

	dse_A0G0 = 0

	LOOP(Gen3, gM,3,1)

        dse_A0G0 = dse_A0G0 - 
     &    3/(4.D0*Pi)*(Alfa1L*
     &        (B1(k2,Mf2(tM,Gen3),Mf2(tM,Gen3)) + 
     &          k2*DB1(k2,Mf2(tM,Gen3),Mf2(tM,Gen3)))*Mf2(tM,Gen3))
     &       /(MW2*SW2*TB)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        dse_A0G0 = dse_A0G0 + 
     &    3/(16.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,tM,Gen3),MSf2(Sfe4,tM,Gen3))*
     &        (-((Kf(Gen3,Gen3,tM) + MUEC*TB*Mf(tM,Gen3))*
     &             USf(Sfe4,2,tM,Gen3)*USfC(Sfe3,1,tM,Gen3)) + 
     &          (KfC(Gen3,Gen3,tM) + MUE*TB*Mf(tM,Gen3))*
     &           USf(Sfe4,1,tM,Gen3)*USfC(Sfe3,2,tM,Gen3))*
     &        ((TB*Kf(Gen3,Gen3,tM) - MUEC*Mf(tM,Gen3))*
     &           USf(Sfe3,2,tM,Gen3)*USfC(Sfe4,1,tM,Gen3) + 
     &          (-(TB*KfC(Gen3,Gen3,tM)) + MUE*Mf(tM,Gen3))*
     &           USf(Sfe3,1,tM,Gen3)*USfC(Sfe4,2,tM,Gen3)))/
     &      (MW2*SW2*TB2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_A0G0 =", dse_A0G0 ENDL
#endif

	if( mssmpart .eq. 1 ) return

	LOOP(Gen3, gM,3,1)

        dse_A0G0 = dse_A0G0 + 
     &    3/(4.D0*Pi)*(Alfa1L*TB*
     &        (B1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)) + 
     &          k2*DB1(k2,Mf2(bX,Gen3),Mf2(bX,Gen3)))*Mf2(bX,Gen3))
     &       /(MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Gen3, gM,3,1)
	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)

        dse_A0G0 = dse_A0G0 - 
     &    3/(16.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,bX,Gen3),MSf2(Sfe4,bX,Gen3))*
     &        (-((TB*Kf(Gen3,Gen3,bX) + MUEC*Mf(bX,Gen3))*
     &             USf(Sfe4,2,bX,Gen3)*USfC(Sfe3,1,bX,Gen3)) + 
     &          (TB*KfC(Gen3,Gen3,bX) + MUE*Mf(bX,Gen3))*
     &           USf(Sfe4,1,bX,Gen3)*USfC(Sfe3,2,bX,Gen3))*
     &        ((Kf(Gen3,Gen3,bX) - MUEC*TB*Mf(bX,Gen3))*
     &           USf(Sfe3,2,bX,Gen3)*USfC(Sfe4,1,bX,Gen3) + 
     &          (-KfC(Gen3,Gen3,bX) + MUE*TB*Mf(bX,Gen3))*
     &           USf(Sfe3,1,bX,Gen3)*USfC(Sfe4,2,bX,Gen3)))/
     &      (MW2*SW2)

	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)
	ENDLOOP(Gen3)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_A0G0 =", dse_A0G0 ENDL
#endif

	if( mssmpart .eq. 2 ) return

	LOOP(Gen3, 1,3,1)

        dse_A0G0 = dse_A0G0 + 
     &    1/(4.D0*Pi)*(Alfa1L*TB*
     &        (B1(k2,Mf2(2,Gen3),Mf2(2,Gen3)) + 
     &          k2*DB1(k2,Mf2(2,Gen3),Mf2(2,Gen3)))*Mf2(2,Gen3))/
     &      (MW2*SW2)

	ENDLOOP(Gen3)

	LOOP(Sfe4, 1,2,1)
	LOOP(Sfe3, 1,2,1)
	LOOP(Gen3, 1,3,1)

        dse_A0G0 = dse_A0G0 - 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        DB0(k2,MSf2(Sfe3,2,Gen3),MSf2(Sfe4,2,Gen3))*
     &        (-((TB*Kf(Gen3,Gen3,2) + MUEC*Mf(2,Gen3))*
     &             USf(Sfe4,2,2,Gen3)*USfC(Sfe3,1,2,Gen3)) + 
     &          (TB*KfC(Gen3,Gen3,2) + MUE*Mf(2,Gen3))*
     &           USf(Sfe4,1,2,Gen3)*USfC(Sfe3,2,2,Gen3))*
     &        ((Kf(Gen3,Gen3,2) - MUEC*TB*Mf(2,Gen3))*
     &           USf(Sfe3,2,2,Gen3)*USfC(Sfe4,1,2,Gen3) + 
     &          (-KfC(Gen3,Gen3,2) + MUE*TB*Mf(2,Gen3))*
     &           USf(Sfe3,1,2,Gen3)*USfC(Sfe4,2,2,Gen3)))/(MW2*SW2)

	ENDLOOP(Gen3)
	ENDLOOP(Sfe3)
	ENDLOOP(Sfe4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_A0G0 =", dse_A0G0 ENDL
#endif

	if( mssmpart .eq. 3 ) return

        dse_A0G0 = dse_A0G0 + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (C2B*MW2*S2B*
     &           (SAB**2*DB0(k2,MA02,Mh02) + 
     &             CAB**2*DB0(k2,MA02,MHH2)) - 
     &          (C2B*MW2*S2B*SAB**2 + 
     &             CBA*CW2*(3*k2*SBA + Mh02*SBA))*DB0(k2,Mh02,MZ2)+
     &            (-(C2B*CAB**2*MW2*S2B) + 
     &             CBA*CW2*(3*k2*SBA + MHH2*SBA))*DB0(k2,MHH2,MZ2)+
     &            CBA*CW2*SBA*
     &           (-3*B0(k2,Mh02,MZ2) + 3*B0(k2,MHH2,MZ2) + 
     &             2*(-B1(k2,MZ2,Mh02) + B1(k2,MZ2,MHH2) + 
     &                k2*(-DB1(k2,MZ2,Mh02) + DB1(k2,MZ2,MHH2))))))
     &       /(CW2**2*SW2)

	LOOP(Cha4, 1,2,1)
	LOOP(Cha3, 1,2,1)

        dup1 = SB*UCha(Cha4,2)*VCha(Cha3,1) + 
     &    CB*UCha(Cha4,1)*VCha(Cha3,2)

        dup2 = -(CB*UCha(Cha3,2)*VCha(Cha4,1)) + 
     &    SB*UCha(Cha3,1)*VCha(Cha4,2)

        dup3 = -(CB*UChaC(Cha4,2)*VChaC(Cha3,1)) + 
     &    SB*UChaC(Cha4,1)*VChaC(Cha3,2)

        dup4 = SB*UChaC(Cha3,2)*VChaC(Cha4,1) + 
     &    CB*UChaC(Cha3,1)*VChaC(Cha4,2)

        dse_A0G0 = dse_A0G0 + 
     &    1/(4.D0*Pi)*(Alfa1L*
     &        (-((dup1*dup3 + dup2*dup4)*
     &             (B1(k2,MCha2(Cha3),MCha2(Cha4)) + 
     &               k2*DB1(k2,MCha2(Cha3),MCha2(Cha4)))) + 
     &          DB0(k2,MCha2(Cha3),MCha2(Cha4))*
     &           (-((dup1*dup3 + dup2*dup4)*MCha2(Cha3)) + 
     &             MCha(Cha3)*MCha(Cha4)*
     &              (dup1*dup2 + 
     &                dup4*
     &                 (-(CB*UChaC(Cha4,2)*VChaC(Cha3,1)) + 
     &                   SB*UChaC(Cha4,1)*VChaC(Cha3,2))))))/SW2

	ENDLOOP(Cha3)
	ENDLOOP(Cha4)

	LOOP(Neu4, 1,4,1)
	LOOP(Neu3, 1,4,1)

	dup5 = SW*ZNeu(Neu3,1) - CW*ZNeu(Neu3,2)

	dup6 = CB*ZNeu(Neu3,3) + SB*ZNeu(Neu3,4)

	dup7 = SW*ZNeu(Neu4,1) - CW*ZNeu(Neu4,2)

	dup8 = -(SW*ZNeu(Neu4,1)) + CW*ZNeu(Neu4,2)

	dup9 = SB*ZNeu(Neu4,3) - CB*ZNeu(Neu4,4)

	dup10 = CB*ZNeu(Neu4,3) + SB*ZNeu(Neu4,4)

	dup11 = SW*ZNeuC(Neu3,1) - CW*ZNeuC(Neu3,2)

	dup12 = CB*ZNeuC(Neu3,3) + SB*ZNeuC(Neu3,4)

	dup13 = SW*ZNeuC(Neu4,1) - CW*ZNeuC(Neu4,2)

	dup14 = -(SW*ZNeuC(Neu4,1)) + CW*ZNeuC(Neu4,2)

	dup15 = SB*ZNeuC(Neu4,3) - CB*ZNeuC(Neu4,4)

	dup16 = CB*ZNeuC(Neu4,3) + SB*ZNeuC(Neu4,4)

        dup17 = (dup12*dup13 + dup11*dup16)*
     &     (dup5*dup9 + dup7*SB*ZNeu(Neu3,3) + 
     &       CB*dup8*ZNeu(Neu3,4)) + 
     &    (dup10*dup5 + dup6*dup7)*
     &     (dup11*dup15 + dup13*SB*ZNeuC(Neu3,3) + 
     &       CB*dup14*ZNeuC(Neu3,4))

        dse_A0G0 = dse_A0G0 + 
     &    1/(16.D0*Pi)*(Alfa1L*
     &        (dup17*(B1(k2,MNeu2(Neu3),MNeu2(Neu4)) + 
     &             k2*DB1(k2,MNeu2(Neu3),MNeu2(Neu4))) + 
     &          DB0(k2,MNeu2(Neu3),MNeu2(Neu4))*
     &           (dup17*MNeu2(Neu3) - 
     &             MNeu(Neu3)*MNeu(Neu4)*
     &              ((dup10*dup5 + dup6*dup7)*
     &                 (dup5*dup9 + dup7*SB*ZNeu(Neu3,3) + 
     &                   CB*dup8*ZNeu(Neu3,4)) + 
     &                (dup12*dup13 + dup11*dup16)*
     &                 (dup11*dup15 + dup13*SB*ZNeuC(Neu3,3) + 
     &                   CB*dup14*ZNeuC(Neu3,4))))))/(CW2*SW2)

	ENDLOOP(Neu3)
	ENDLOOP(Neu4)

#ifdef DETAILED_DEBUG
	DHIGGS "dse_A0G0 =", dse_A0G0 ENDL
#endif

	end

