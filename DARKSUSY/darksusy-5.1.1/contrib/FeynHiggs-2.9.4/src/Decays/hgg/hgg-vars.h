#if 0
* hgg-vars.h
* variable declarations
* generated by FormCalc 7.5 on 27-Sep-2012 9:10
* this file is part of FeynHiggs
#endif

#ifndef VARS_H
#define VARS_H

#include "externals.h"
#include "types.h"
#include "debug.h"

#else

#include "Decay.h"

	ComplexType Sub14, Sub19, Sub24, Sub25, Sub8, Sub9
	ComplexType Sub10(6,3,3), Sub11(3), Sub12(3), Sub13(3)
	ComplexType Sub15(6,3), Sub16(3), Sub17(3), Sub18(3)
	ComplexType Sub20(6,3), Sub21(6,3), Sub22(6,3,3), Sub23(6)
	ComplexType Sub26(6,3,3), Sub27(3), Sub28(3), Sub29(3)
	ComplexType Sub30(6,3), Sub31(3), Sub32(3), Sub33(3)
	ComplexType Sub34(6,3), Sub35(6,3), Sub36(6,3,3), Sub37(6)
	common /hgg_abbrev1s/ Sub14, Sub19, Sub24, Sub25, Sub8, Sub9
	common /hgg_abbrev1s/ Sub10, Sub11, Sub12, Sub13, Sub15
	common /hgg_abbrev1s/ Sub16, Sub17, Sub18, Sub20, Sub21
	common /hgg_abbrev1s/ Sub22, Sub23, Sub26, Sub27, Sub28
	common /hgg_abbrev1s/ Sub29, Sub30, Sub31, Sub32, Sub33
	common /hgg_abbrev1s/ Sub34, Sub35, Sub36, Sub37

	ComplexType Abb1, Abb2, Abb3, AbbSum1, AbbSum2, Eps1, Pair1
	ComplexType Pair2, Pair3, Pair4, Pair5, Sub1, Sub2(3)
	ComplexType Sub3(3), Sub4(3), Sub5(3), Sub6(3), Sub7(3)
	common /hgg_abbrev1hel/ Abb1, Abb2, Abb3, AbbSum1, AbbSum2
	common /hgg_abbrev1hel/ Eps1, Pair1, Pair2, Pair3, Pair4
	common /hgg_abbrev1hel/ Pair5, Sub1, Sub2, Sub3, Sub4, Sub5
	common /hgg_abbrev1hel/ Sub6, Sub7

	ComplexType lint1(3), lint2(Ncc,3), lint3(Ncc,3), lint4(6)
	ComplexType lint5(Ncc,6), lint6(6), lint7(Ncc,6), lint8(3)
	common /hgg_loopint1s/ lint1, lint2, lint3, lint4, lint5
	common /hgg_loopint1s/ lint6, lint7, lint8

	integer seq(2), Hel(3)
	common /hgg_helvars/ seq, Hel

	integer All4, Gen4, Ind1, Ind2
	common /hgg_indices/ All4, Gen4, Ind1, Ind2

	ComplexType Cloop(1), MatSUN(1,1)
	common /hgg_formfactors/ Cloop, MatSUN

#endif
