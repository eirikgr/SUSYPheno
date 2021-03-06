#if 0
* huu-v1_All4_Cha4_Cha5_Ind1_Ind2.F
* form factors for {S[10, {hno}]} -> {F[3, {gno1, Index[Colour, 2]}], -F[3, {gno2, Index[Colour, 3]}]}
* generated by FormCalc 7.5 on 27-Sep-2012 9:19
* this file is part of FeynHiggs
#endif

#include "huu-vars.h"

	subroutine huu_v1_All4_Cha4_Cha5_Ind1Ind2
	implicit none

#include "huu-vars.h"

#include "inline.h"

        Cloop(1) = Cloop(1) + 
     &    1/(16.D0*Pi)*(-(16*Pi*
     &           (CB*MW2*SB*SW*SW2*
     &             (lint43(cc0,All4,Cha4,Cha5)*
     &                Sub265(All4,Cha4,Cha5,Ind1,Ind2) + 
     &               lint43(cc2,All4,Cha4,Cha5)*
     &                Sub266(All4,Cha4,Cha5,Ind1,Ind2) + 
     &               lint43(cc1,All4,Cha4,Cha5)*
     &                Sub267(All4,Cha4,Cha5,Ind1,Ind2)))) - 
     &        AlfaMZ*ELMZ*lint42(All4,Cha5)*
     &         Sub263(All4,Cha4,Cha5,Ind1,Ind2))/(CB*MW2*SB*SW*SW2)

#ifdef DETAILED_DEBUG
	DCOUP "v1_All4_Cha4_Cha5_Ind1_Ind2: Cloop(1) =", Cloop(1) ENDL
#endif

	end
