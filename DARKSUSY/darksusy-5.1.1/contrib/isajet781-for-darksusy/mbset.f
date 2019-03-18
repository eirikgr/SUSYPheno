CDECK  ID>, MBSET.
      SUBROUTINE MBSET
C
C          SET PARAMETERS FOR GENERATING MINBIAS EVENTS OR BEAM JETS,
C          ALLOWING DIFFERENT PARAMETERS FOR TWO CASES.
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      COMMON/MBPAR/PUD0,PJSPN,PISPN,SIGQT0,XGEN0(2),PMIX01(3,2)
     1,PMIX02(3,2),PBARY0
      SAVE /MBPAR/
      REAL      PUD0,PJSPN,PISPN,SIGQT0,XGEN0,PMIX01,PMIX02,PBARY0
      INTEGER   LIMPOM
      PARAMETER (LIMPOM=20)
      COMMON/MBGEN/POMWT(LIMPOM),POMGEN(LIMPOM),MNPOM,MXPOM,PDIFFR,
     $NPOM,XBARY(2),DXBARY(2),XPOM(LIMPOM,2)
      SAVE /MBGEN/
      INTEGER   MNPOM,MXPOM,NPOM
      REAL      POMWT,POMGEN,PDIFFR,XBARY,DXBARY,XPOM
      COMMON/PRIMAR/NJET,SCM,HALFE,ECM,IDIN(2),NEVENT,NTRIES,NSIGMA,
     $WRTLHE
      SAVE /PRIMAR/
      INTEGER   NJET,IDIN,NEVENT,NTRIES,NSIGMA
      LOGICAL   WRTLHE
      REAL      SCM,HALFE,ECM
      COMMON/TOTALS/NKINPT,NWGEN,NKEEP,SUMWT,WT
      SAVE /TOTALS/
      INTEGER   NKINPT,NWGEN,NKEEP
      REAL      SUMWT,WT
      INTEGER MXKEYS
      PARAMETER (MXKEYS=20)
      COMMON/KEYS/IKEYS,KEYON,KEYS(MXKEYS)
      COMMON/XKEYS/REAC
      SAVE /KEYS/,/XKEYS/
      LOGICAL KEYS
      LOGICAL KEYON
      CHARACTER*8 REAC
      INTEGER   IKEYS
C
C
C          DN/DY INCREASES WITH LOG(S). INCLUDED IN SPLITTING FUNCTION
C          BECAUSE AVERAGE MULTIPLICITY COMES FROM SINGLE CHAIN GRAPH.
      XGEN0(1)=.9
      XGEN0(2)=1.+0.35*ALOG(ECM/60.)
C
C          POMWT ARE (RELATIVE) PROBABILITIES FOR N CUT POMERONS.
C          PDIFFR IS DIFFRACTIVE PROBABILITY.
C          SIGQT0 IS MEAN PT.
      IF(KEYS(4)) THEN
        PDIFFR=.15
        SIGQT0=.35
        PSUM=0.
        DO 100 I=1,LIMPOM
        POMWT(I)=(1.+4.*I**2)*EXP(-1.8*I)
        PSUM=PSUM+POMWT(I)
100     CONTINUE
      ELSE
        PDIFFR=0.
        SIGQT0=.45
        PSUM=0.
        DO 110 I=1,LIMPOM
        POMWT(I)=(1.+4.*I**2)*EXP(-1.8*I)
        PSUM=PSUM+POMWT(I)
110     CONTINUE
        POMWT(1)=.1*POMWT(1)
        POMWT(2)=.2*POMWT(2)
        POMWT(3)=.5*POMWT(3)
      ENDIF
C
C          RENORMALIZE POMWT.
      PSUM=1./PSUM
      DO 200 I=1,LIMPOM
      POMWT(I)=PSUM*POMWT(I)
200   CONTINUE
      PSUM=0.
      DO 210 I=MNPOM,MXPOM
      PSUM=PSUM+POMWT(I)
210   CONTINUE
C
C          POMGEN IS USED TO SELECT NUMBER OF POMERONS.
      PGEN=0.
      PSUM=1./PSUM
      DO 300 I=1,LIMPOM
      POMGEN(I)=0.
300   CONTINUE
      DO 310 I=MNPOM,MXPOM
      PGEN=PGEN+PSUM*POMWT(I)
      POMGEN(I)=PGEN
310   CONTINUE
      POMGEN(MXPOM)=1.
C
C          SET /TOTALS/ FOR MINBIAS EVENTS USING LOG**2(S) FIT TO
C          TOTAL CROSS SECTION.
      IF(KEYS(4)) THEN
        SIGTOT=25.65*(1.+.0102*ALOG(SCM/1.76)**2)
        SIGTOT=PSUM*SIGTOT
        NKINPT=NEVENT
        SUMWT=SIGTOT*NKINPT
      ENDIF
C
      RETURN
      END