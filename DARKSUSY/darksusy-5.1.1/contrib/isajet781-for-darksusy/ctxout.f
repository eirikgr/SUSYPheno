CDECK  ID>, CTXOUT.
      SUBROUTINE CTXOUT(NVC,VC,MXVC)
C-----------------------------------------------------------------------
C  Purpose:
C          Save the context for an ISAJET job:
C          Save in NVC words of VC all common blocks NOT associated only
C          with a single event. Call this and CTXIN to generate mixed
C          events.
C          PARAMETER (MXVC=20000)
C          REAL    VC(MXVC)
C          ...
C          CALL CTXIN(NVC,VC,MXVC)
C
C          Note that the MSSM common blocks are not saved, so different
C          SUSY runs cannot be mixed.
C
C          Ver. 7.02: Equivalenced dummy variables to avoid mixed 
C                     arguments in MOVLEV or multiple EQUIVALENCEd
C                     arguments to CTXIN/CTXOUT.
C
C  Author:
C          F.E. Paige, April 1992     
C-----------------------------------------------------------------------
      IMPLICIT NONE
C          LOOK must be dimensioned to the maximum value of INDEX.
      INTEGER   MXLOOK
      PARAMETER (MXLOOK=500)
      INTEGER   MXDKY
      PARAMETER (MXDKY=3000)
      COMMON/DKYTAB/LOOK(MXLOOK),CBR(MXDKY),MODE(5,MXDKY),MELEM(MXDKY)
      SAVE /DKYTAB/
      INTEGER   LOOK,MODE,MELEM
      REAL      CBR
      COMMON/DYLIM/QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     2  THWMAX,PHWMIN,PHWMAX
     3  ,SETLMQ(12)
      SAVE /DYLIM/
      LOGICAL SETLMQ
      EQUIVALENCE(BLIM1(1),QMIN)
      REAL      QMIN,QMAX,QTMIN,QTMAX,YWMIN,YWMAX,XWMIN,XWMAX,THWMIN,
     +          THWMAX,PHWMIN,PHWMAX,BLIM1(12)
      COMMON/DYPAR/FLW,RNU2(3),ANORM(3),QPOW(3),PTPOW(3)
      SAVE /DYPAR/
      LOGICAL FLW
      REAL      RNU2,ANORM,QPOW,PTPOW
      COMMON/EEPAR/SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,
     $UPSLON,SIGZ,IBREM,IBEAM,GAMGAM
      SAVE /EEPAR/
      REAL SGMXEE,PLEP,PLEM,RSHMIN,RSHMAX,UPSLON,SIGZ
      LOGICAL IBREM,IBEAM,GAMGAM
      COMMON/FINAL/NKINF,SIGF,ALUM,ACCEPT,NRECS
      SAVE /FINAL/
      INTEGER   NKINF,NRECS
      REAL      SIGF,ALUM,ACCEPT
      INTEGER   MXFORC
      PARAMETER (MXFORC=40)
      COMMON/FORCE/NFORCE,IFORCE(MXFORC),MFORCE(5,MXFORC)
     $,LOOK2(2,MXFORC),LOOKST(MXFORC),MEFORC(MXFORC)
      SAVE /FORCE/
      INTEGER   NFORCE,IFORCE,MFORCE,LOOK2,LOOKST,MEFORC
      COMMON/FRGPAR/PUD,PBARY,SIGQT,PEND,XGEN(8),PSPIN1(8),
     $PMIX1(3,2),PMIX2(3,2),XGENSS(9)
      SAVE /FRGPAR/
      EQUIVALENCE (PMIX1(1,1),PMIXX1(1))
      EQUIVALENCE (PMIX2(1,1),PMIXX2(1))
      EQUIVALENCE(FRPAR(1),PUD)
      REAL      PUD,PBARY,SIGQT,PEND,XGEN,PSPIN1,PMIX1,PMIX2,XGENSS,
     +          PMIXX1(6),PMIXX2(6),FRPAR(32)
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
      DOUBLE PRECISION ANWWWW,ADWWWW,AIWWWW
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS
      COMMON/IDRUN/IDVER,IDG(2),IEVT,IEVGEN
      SAVE /IDRUN/
      INTEGER   IDVER,IDG,IEVT,IEVGEN
      COMMON/ISLOOP/NEVOLV,NFRGMN,IEVOL,IFRG
      SAVE /ISLOOP/
      INTEGER NEVOLV,NFRGMN,IEVOL,IFRG
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
C          Jet limits
      INTEGER MXLIM
      PARAMETER (MXLIM=8)
      INTEGER MXLX12
      PARAMETER (MXLX12=12*MXLIM)
      COMMON/JETLIM/PMIN(MXLIM),PMAX(MXLIM),PTMIN(MXLIM),PTMAX(MXLIM),
     $YJMIN(MXLIM),YJMAX(MXLIM),PHIMIN(MXLIM),PHIMAX(MXLIM),
     $XJMIN(MXLIM),XJMAX(MXLIM),THMIN(MXLIM),THMAX(MXLIM),
     $SETLMJ(12*MXLIM)
      SAVE /JETLIM/
      COMMON/FIXPAR/FIXP(MXLIM),FIXPT(MXLIM),FIXYJ(MXLIM),
     $FIXPHI(MXLIM),FIXXJ(MXLIM),FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      SAVE /FIXPAR/
      COMMON/SGNPAR/CTHS(2,MXLIM),THS(2,MXLIM),YJS(2,MXLIM),XJS(2,MXLIM)  
      SAVE /SGNPAR/
      REAL      PMIN,PMAX,PTMIN,PTMAX,YJMIN,YJMAX,PHIMIN,PHIMAX,XJMIN,
     +          XJMAX,THMIN,THMAX,BLIMS(12*MXLIM),CTHS,THS,YJS,XJS
      LOGICAL SETLMJ
      LOGICAL FIXQM,FIXQT,FIXYW,FIXXW,FIXPHW
      LOGICAL FIXP,FIXPT,FIXYJ,FIXPHI,FIXXJ
      EQUIVALENCE(BLIMS(1),PMIN(1))
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
      COMMON /LIMEVL/ ETTHRS,CONCUT,USELIM
      SAVE /LIMEVL/
      REAL ETTHRS,CONCUT
      LOGICAL USELIM
      COMMON/LSTPRT/LSTPRT
      SAVE /LSTPRT/
      INTEGER   LSTPRT
      INTEGER   LIMPOM
      PARAMETER (LIMPOM=20)
      COMMON/MBGEN/POMWT(LIMPOM),POMGEN(LIMPOM),MNPOM,MXPOM,PDIFFR,
     $NPOM,XBARY(2),DXBARY(2),XPOM(LIMPOM,2)
      SAVE /MBGEN/
      INTEGER   MNPOM,MXPOM,NPOM
      REAL      POMWT,POMGEN,PDIFFR,XBARY,DXBARY,XPOM
      COMMON/MBPAR/PUD0,PJSPN,PISPN,SIGQT0,XGEN0(2),PMIX01(3,2)
     1,PMIX02(3,2),PBARY0
      SAVE /MBPAR/
      REAL      PUD0,PJSPN,PISPN,SIGQT0,XGEN0,PMIX01,PMIX02,PBARY0
      COMMON/NODCAY/NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      LOGICAL NODCAY,NOETA,NOPI0,NONUNU,NOEVOL,NOHADR,NOGRAV,
     $NOB,NOTAU
      SAVE /NODCAY/
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/PRTOUT/NEVPRT,NJUMP
      SAVE /PRTOUT/
      INTEGER   NEVPRT,NJUMP
      COMMON/PTPAR/PTFUN1,PTFUN2,PTGEN1,PTGEN2,PTGEN3,SIGMAX
      SAVE /PTPAR/
      REAL      PTFUN1,PTFUN2,PTGEN1,PTGEN2,PTGEN3,SIGMAX
      INTEGER MXGOQ,MXGOJ
      PARAMETER (MXGOQ=85,MXGOJ=8)
      COMMON/Q1Q2/GOQ(MXGOQ,MXGOJ),GOALL(MXGOJ),GODY(4),STDDY,
     $GOWW(25,2),ALLWW(2),GOWMOD(25,MXGOJ)
      SAVE /Q1Q2/
      LOGICAL GOQ,GOALL,GODY,STDDY,GOWW,ALLWW,GOWMOD
      COMMON/QCDPAR/ALAM,ALAM2,CUTJET,ISTRUC
      SAVE /QCDPAR/
      INTEGER   ISTRUC
      REAL      ALAM,ALAM2,CUTJET
      COMMON/QLMASS/AMLEP(100),NQLEP,NMES,NBARY
      SAVE /QLMASS/
      INTEGER   NQLEP,NMES,NBARY
      REAL      AMLEP
      COMMON/TCPAR/TCMRHO,TCGRHO
      SAVE /TCPAR/
      REAL TCMRHO,TCGRHO
      COMMON/TIMES/TIME1,TIME2
      SAVE /TIMES/
      REAL      TIME1,TIME2
      COMMON/TOTALS/NKINPT,NWGEN,NKEEP,SUMWT,WT
      SAVE /TOTALS/
      INTEGER   NKINPT,NWGEN,NKEEP
      REAL      SUMWT,WT
      INTEGER MXTYPE
      PARAMETER (MXTYPE=8)
      COMMON/TYPES/LOC(100),NTYP,NJTTYP(MXTYPE),NWWTYP(2),NWMODE(3)
      COMMON/XTYPES/PARTYP(40),TITLE(10),JETYP(30,MXTYPE),WWTYP(30,2)
     $,WMODES(30,3)
      SAVE /TYPES/,/XTYPES/
      CHARACTER*8 JETYP,WWTYP,TITLE,PARTYP,WMODES
      INTEGER   LOC,NTYP,NJTTYP,NWWTYP,NWMODE
      COMMON/WCON/SIN2W,WMASS(4),WGAM(4),AQ(12,4),BQ(12,4),COUT(4),
     1MATCH(25,4),WCBR(25,4),CUTOFF,CUTPOW,TBRWW(4,2),RBRWW(12,4,2),EZ,
     2AQDP(12,4),BQDP(12,4),EZDP,WFUDGE
      SAVE /WCON/
      DOUBLE PRECISION AQDP,BQDP,EZDP
      INTEGER   MATCH
      REAL      SIN2W,WMASS,WGAM,AQ,BQ,COUT,WCBR,CUTOFF,CUTPOW,TBRWW,
     +          RBRWW,EZ,WFUDGE
      COMMON/WCON2/CUMWBR(25,3)
      REAL CUMWBR
C
      INTEGER NVC,MXVC,NC,NN,I
      REAL VC(MXVC)
      CHARACTER*8 CLIST(290)
      EQUIVALENCE (CLIST(1),PARTYP(1))
C
C          Dummy real variables for integers
      REAL VLOOK(MXLOOK+6*MXDKY)
      EQUIVALENCE (VLOOK(1),LOOK(1))
      REAL VNKINF(5)
      EQUIVALENCE (VNKINF(1),NKINF)
      REAL VFORCE(9*MXFORC+1)
      EQUIVALENCE (VFORCE(1),NFORCE)
      REAL VIDVER(5)
      EQUIVALENCE (VIDVER(1),IDVER)
      REAL VEVOLV(4)
      EQUIVALENCE (VEVOLV(1),NEVOLV)
      REAL VITDKY(4)
      EQUIVALENCE (VITDKY(1),ITDKY)
      REAL VIKEYS(12)
      EQUIVALENCE (VIKEYS(1),IKEYS)
      REAL VSTPRT
      EQUIVALENCE (VSTPRT,LSTPRT)
      REAL VNJET(9)
      EQUIVALENCE (VNJET(1),NJET)
      REAL VEVPRT(2)
      EQUIVALENCE (VEVPRT(1),NEVPRT)
      REAL VKINPT(5)
      EQUIVALENCE (VKINPT(1),NKINPT)
      REAL VLOC(100)
      EQUIVALENCE (VLOC(1),LOC(1))
C          Dummy real variables for logicals
      REAL VFLW(13)
      EQUIVALENCE (VFLW(1),FLW)
      REAL VNODCY(6)
      EQUIVALENCE (VNODCY(1),NODCAY)
      REAL VGOQ(3*MXGOQ+135)
      EQUIVALENCE (VGOQ(1),GOQ(1,1))
C
      NC=0
C          DKYTAB
      NN=MXLOOK+6*MXDKY
      CALL MOVLEV(VLOOK(1),VC(NC+1),NN)
      NC=NC+NN
C          DYLIM
      CALL MOVLEV(QMIN,VC(NC+1),24)
      NC=NC+24
C          DYPAR
      CALL MOVLEV(VFLW(1),VC(NC+1),13)
      NC=NC+13
C          EEPAR
      CALL MOVLEV(SGMXEE,VC(NC+1),1)
      NC=NC+1
C          FINAL
      CALL MOVLEV(VNKINF(1),VC(NC+1),5)
      NC=NC+5
C          FORCE
      NN=9*MXFORC+1
      CALL MOVLEV(VFORCE(1),VC(NC+1),NN)
      NC=NC+NN
C          FRGPAR
      CALL MOVLEV(PUD,VC(NC+1),41)
      NC=NC+41
C          HCON
      CALL MOVLEV(HMASS,VC(NC+1),69)
      NC=NC+69
C          IDRUN
      CALL MOVLEV(VIDVER(1),VC(NC+1),5)
      NC=NC+5
C          ISLOOP
      CALL MOVLEV(VEVOLV(1),VC(NC+1),4)
      NC=NC+4
C          ITAPES
      CALL MOVLEV(VITDKY(1),VC(NC+1),4)
      NC=NC+4
C          JETLIM
      CALL MOVLEV(PMIN(1),VC(NC+1),72)
      NC=NC+72
C          KEYS
      CALL MOVLEV(VIKEYS(1),VC(NC+1),12)
      NC=NC+12
      CALL CTXC2I(REAC,VC(NC+1),8)
      NC=NC+8
C          LIMEVL
      CALL MOVLEV(ETTHRS,VC(NC+1),3)
      NC=NC+3
C          LSTPRT
      CALL MOVLEV(VSTPRT,VC(NC+1),1)
      NC=NC+1
C          MBGEN
      NN=4*LIMPOM+8
      CALL MOVLEV(POMWT(1),VC(NC+1),NN)
      NC=NC+NN
C          MBPAR
      CALL MOVLEV(PUD0,VC(NC+1),19)
      NC=NC+19
C          NODCAY
      CALL MOVLEV(VNODCY(1),VC(NC+1),6)
      NC=NC+6
C          PRIMAR
      CALL MOVLEV(VNJET(1),VC(NC+1),9)
      NC=NC+9
C          PRTOUT
      CALL MOVLEV(VEVPRT(1),VC(NC+1),2)
      NC=NC+2
C          PTPAR
      CALL MOVLEV(PTFUN1,VC(NC+1),6)
      NC=NC+6
C          Q1Q2
      CALL MOVLEV(VGOQ(1),VC(NC+1),3*MXGOQ+135)
      NC=NC+3*MXGOQ+135
C          QCDPAR
      CALL MOVLEV(ALAM,VC(NC+1),4)
      NC=NC+4
C          QLMASS
      CALL MOVLEV(AMLEP(1),VC(NC+1),55)
      NC=NC+55
C          TCPAR
      CALL MOVLEV(TCMRHO,VC(NC+1),2)
      NC=NC+2
C          TIMES
      CALL MOVLEV(TIME1,VC(NC+1),2)
      NC=NC+2
C          TOTALS
      CALL MOVLEV(VKINPT(1),VC(NC+1),5)
      NC=NC+5
C          TYPES
      CALL MOVLEV(VLOC(1),VC(NC+1),100)
      NC=NC+100
      DO 100 I=1,290
        CALL CTXC2I(CLIST(I),VC(NC+1),8)
        NC=NC+8
100   CONTINUE
C          WCON
      NN=514+97
      CALL MOVLEV(SIN2W,VC(NC+1),NN)
      NC=NC+NN
C
      IF(NC.LE.MXVC) THEN
        NVC=NC
        RETURN
      ELSE
        WRITE(ITLIS,9000) NC
9000    FORMAT(//' ERROR IN CTXOUT, NC = ',I5)
        STOP99
      ENDIF
      END
