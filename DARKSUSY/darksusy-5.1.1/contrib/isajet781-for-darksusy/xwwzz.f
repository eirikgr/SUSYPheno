CDECK  ID>, XWWZZ.
      SUBROUTINE XWWZZ
C
C          SET UP W+ W- -> Z0 Z0 AMPLITUDES AS RATIONAL FUNCTIONS OF Z
C
C          RE(F(Z,L)) = SUM(I,J)(ANWWWW(I+1,J,L)*Z**I
C                                  /(ADWWWW(1,J)+ADWWWW(2,J)*Z))
C          IM(F(Z,L)) = AIWWWW(L)   (INDEPENDENT OF Z)
C          J LABELS PIECES WITH SAME DENOMINATOR.
C          L=1 FOR 0,0; L=2 FOR 1,-1; L=3 FOR 1,1; L=4 FOR 0,1
C
C          *NOTE* A FACTOR OF SIN(THETA)/SQRT(2) IS REMOVED FROM F01
C
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
      COMMON/JETPAR/P(3),PT(3),YJ(3),PHI(3),XJ(3),TH(3),CTH(3),STH(3)
     1 ,JETTYP(3),SHAT,THAT,UHAT,QSQ,X1,X2,PBEAM(2)
     2 ,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,JWTYP
     3 ,ALFQSQ,CTHW,STHW,Q0W
     4 ,INITYP(2),ISIGS,PBEAMS(5)
      SAVE /JETPAR/
      INTEGER   JETTYP,JWTYP,INITYP,ISIGS
      REAL      P,PT,YJ,PHI,XJ,TH,CTH,STH,SHAT,THAT,UHAT,QSQ,X1,X2,
     +          PBEAM,QMW,QW,QTW,YW,XW,THW,QTMW,PHIW,SHAT1,THAT1,UHAT1,
     +          ALFQSQ,CTHW,STHW,Q0W,PBEAMS
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
      COMMON/HCON/ANWWWW(4,4,4),ADWWWW(2,4),AIWWWW(4)
     $,HMASS,HGAM,HGAMS(29),ETAHGG,MATCHH(29),ZSTARS(4,2)
     $,IHTYPE,HGAMSS(85,85)
      SAVE /HCON/
      DOUBLE PRECISION ANWWWW,ADWWWW,AIWWWW
      INTEGER   MATCHH,IHTYPE
      REAL      HMASS,HGAM,HGAMS,ETAHGG,ZSTARS,HGAMSS
      DOUBLE PRECISION WM,ZM,ZM2,ZM3,ZM4,ZM5,ZM6,HM,HM2,HM3,HM4,HG,HG2
     $,PROPH,RTS,S,S2,S3,SW,QQ0,QQI,QQF
C
C          USE UNITS OF WM TO AVOID LARGE NUMBERS - NOTE ANWWWW/ADWWWW
C          AND AIWWWW ARE DIMENSIONLESS
      WM=WMASS(2)
      ZM=WMASS(4)/WM
      ZM2=ZM**2
      ZM3=ZM**3
      ZM4=ZM**4
      ZM5=ZM**5
      ZM6=ZM**6
      HM=HMASS/WM
      HM2=HM**2
      HM3=HM**3
      HM4=HM**4
      HG=HGAM/WM
      HG2=HG**2
      RTS=QMW/WM
      S=RTS**2
      S2=S**2
      S3=S**3
      PROPH=(S-HM2)**2+(HM*HG)**2
C          CORRECT SIGN OF HIGGS AMPLITUDE.
      PROPH=-PROPH
C
      CW=1./ZM
      CW2=CW**2
      SW2=1.-CW2
      SW=SQRT(SW2)
      QQ0=.5*RTS
      QQI=.5*SQRT(S-4.)
      QQF=.5*SQRT(S-4.*ZM2)
      GSQ=4.*PI*ALFA/SW2
C
C          FROM WWZZ3.EX
      ANWWWW(1,1,1) = -1.60E+01 * ((HM2 * ZM3) / (CW * PROPH))
     $ + 1.60E+01 * ((S * ZM3) / (CW * PROPH) ) - 8.00E+00 * ((S2
     $ * ZM) / (CW * PROPH)) - 8.00E+00 * ((S2 * ZM3) / (CW * PROPH))
     $ + 4.00E+00 * ((S3 * ZM) / (CW * PROPH)) + 8.00E+00 * ((HM2 * S
     $ * ZM) / (CW * PROPH)) + 8.00E+00 * ((HM2 * S * ZM3) / (CW
     $ * PROPH)) - 4.00E+00 * ((HM2 * S2 * ZM) / (CW * PROPH))
     $ + 8.00E+00 * CW2 * S - 6.00E+00 * CW2 * S2 + 8.00E+00 * CW2
     $ * S * ZM2
      ANWWWW(1,1,2) = 4.00E+00 * CW2 * S * ZM2
      ANWWWW(1,1,3) = -1.60E+01 * ((HM2 * ZM3) / (CW * PROPH))
     $ + 1.60E+01 * ((S * ZM3) / (CW * PROPH) ) - 8.00E+00 * ((S2
     $ * ZM3) / (CW * PROPH)) + 8.00E+00 * ((HM2 * S * ZM3) / (CW
     $ * PROPH)) - 3.20E+01 * CW2 * ZM2 + 1.20E+01 * CW2 * S * ZM2
      ANWWWW(1,1,4) = 0.00E+00
      ANWWWW(1,2,1) = -4.00E+00 * CW2 * S2 + 3.00E+00 * CW2 * S3
     $ - 9.60E+01 * CW2 * ZM4 + 3.20E+01 * CW2 * ZM6 + 8.00E+00 * CW2
     $ * S * ZM2 + 1.60E+01 * CW2 * S * ZM4 - 8.00E+00 * CW2 * S
     $ * ZM6 - 1.00E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(1,2,2) = 6.40E+01 * CW2 * ZM2 + 1.20E+01 * CW2 * S * ZM2
     $ - 4.00E+00 * CW2 * S * ZM4 - 2.00E+00 * CW2 * S2 * ZM2
      ANWWWW(1,2,3) = -6.40E+01 * CW2 * ZM2 + 2.00E+01 * CW2 * S * ZM2
     $ + 4.00E+00 * CW2 * S * ZM4 - 6.00E+00 * CW2 * S2 * ZM2
      ANWWWW(1,2,4) = 1.92E+02 * CW2 * QQ0 * QQF * QQI * ZM
     $ - 3.20E+01 * CW2 * QQ0 * QQF * QQI * ZM3 + 3.20E+01 * CW2
     $ * QQ0 * QQF * QQI * ZM5 + 1.60E+01 * CW2 * QQ0 * QQF * QQI * S
     $ * ZM
      ANWWWW(1,3,1) = -4.00E+00 * CW2 * S2 + 3.00E+00 * CW2 * S3
     $ - 9.60E+01 * CW2 * ZM4 + 3.20E+01 * CW2 * ZM6 + 8.00E+00 * CW2
     $ * S * ZM2 + 1.60E+01 * CW2 * S * ZM4 - 8.00E+00 * CW2 * S
     $ * ZM6 - 1.00E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(1,3,2) = 6.40E+01 * CW2 * ZM2 + 1.20E+01 * CW2 * S * ZM2
     $ - 4.00E+00 * CW2 * S * ZM4 - 2.00E+00 * CW2 * S2 * ZM2
      ANWWWW(1,3,3) = -6.40E+01 * CW2 * ZM2 + 2.00E+01 * CW2 * S * ZM2
     $ + 4.00E+00 * CW2 * S * ZM4 - 6.00E+00 * CW2 * S2 * ZM2
      ANWWWW(1,3,4) = -1.92E+02 * CW2 * QQ0 * QQF * QQI * ZM
     $ + 3.20E+01 * CW2 * QQ0 * QQF * QQI * ZM3 - 3.20E+01 * CW2
     $ * QQ0 * QQF * QQI * ZM5 - 1.60E+01 * CW2 * QQ0 * QQF * QQI * S
     $ * ZM
      ANWWWW(1,4,1) = 0.00E+00
      ANWWWW(1,4,2) = 0.00E+00
      ANWWWW(1,4,3) = 0.00E+00
      ANWWWW(1,4,4) = 0.00E+00
      ANWWWW(2,1,1) = 0.00E+00
      ANWWWW(2,1,2) = 0.00E+00
      ANWWWW(2,1,3) = 0.00E+00
      ANWWWW(2,1,4) = 8.00E+00 * CW2 * QQ0 * S * ZM
      ANWWWW(2,2,1) = 4.80E+01 * CW2 * QQF * QQI * S - 2.00E+01 * CW2
     $ * QQF * QQI * S2 + 6.40E+01 * CW2 * QQF * QQI * ZM2 + 9.60E+01
     $ * CW2 * QQF * QQI * S * ZM2 - 1.60E+01 * CW2 * QQF * QQI * S
     $ * ZM4
      ANWWWW(2,2,2) = 8.00E+00 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(2,2,3) = 5.60E+01 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(2,2,4) = 1.28E+02 * CW2 * QQ0 * ZM + 1.92E+02 * CW2 * QQ0
     $ * ZM3 - 3.20E+01 * CW2 * QQ0 * S * ZM + 2.40E+01 * CW2 * QQ0
     $ * S * ZM3 - 8.00E+00 * CW2 * QQ0 * S * ZM5 - 8.00E+00 * CW2
     $ * QQ0 * S2 * ZM
      ANWWWW(2,3,1) = -4.80E+01 * CW2 * QQF * QQI * S + 2.00E+01 * CW2
     $ * QQF * QQI * S2 - 6.40E+01 * CW2 * QQF * QQI * ZM2
     $ - 9.60E+01 * CW2 * QQF * QQI * S * ZM2 + 1.60E+01 * CW2 * QQF
     $ * QQI * S * ZM4
      ANWWWW(2,3,2) = -8.00E+00 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(2,3,3) = -5.60E+01 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(2,3,4) = 1.28E+02 * CW2 * QQ0 * ZM + 1.92E+02 * CW2 * QQ0
     $ * ZM3 - 3.20E+01 * CW2 * QQ0 * S * ZM + 2.40E+01 * CW2 * QQ0
     $ * S * ZM3 - 8.00E+00 * CW2 * QQ0 * S * ZM5 - 8.00E+00 * CW2
     $ * QQ0 * S2 * ZM
      ANWWWW(2,4,1) = 0.00E+00
      ANWWWW(2,4,2) = 0.00E+00
      ANWWWW(2,4,3) = 0.00E+00
      ANWWWW(2,4,4) = 0.00E+00
      ANWWWW(3,1,1) = 2.00E+00 * CW2 * S2
      ANWWWW(3,1,2) = -4.00E+00 * CW2 * S * ZM2
      ANWWWW(3,1,3) = 4.00E+00 * CW2 * S * ZM2
      ANWWWW(3,1,4) = 0.00E+00
      ANWWWW(3,2,1) = 3.20E+01 * CW2 * S - 1.60E+01 * CW2 * S2 + CW2
     $ * S3 + 9.60E+01 * CW2 * S * ZM2 + 3.20E+01 * CW2 * S * ZM4
     $ - 2.20E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(3,2,2) = -6.40E+01 * CW2 * ZM2 - 1.20E+01 * CW2 * S
     $ * ZM2 + 4.00E+00 * CW2 * S * ZM4 + 2.00E+00 * CW2 * S2 * ZM2
      ANWWWW(3,2,3) = 6.40E+01 * CW2 * ZM2 + 1.20E+01 * CW2 * S * ZM2
     $ + 2.80E+01 * CW2 * S * ZM4 - 1.00E+01 * CW2 * S2 * ZM2
      ANWWWW(3,2,4) = 1.60E+01 * CW2 * QQ0 * QQF * QQI * S * ZM
      ANWWWW(3,3,1) = 3.20E+01 * CW2 * S - 1.60E+01 * CW2 * S2 + CW2
     $ * S3 + 9.60E+01 * CW2 * S * ZM2 + 3.20E+01 * CW2 * S * ZM4
     $ - 2.20E+01 * CW2 * S2 * ZM2 + 2.00E+00 * CW2 * S2 * ZM4
      ANWWWW(3,3,2) = -6.40E+01 * CW2 * ZM2 - 1.20E+01 * CW2 * S
     $ * ZM2 + 4.00E+00 * CW2 * S * ZM4 + 2.00E+00 * CW2 * S2 * ZM2
      ANWWWW(3,3,3) = 6.40E+01 * CW2 * ZM2 + 1.20E+01 * CW2 * S * ZM2
     $ + 2.80E+01 * CW2 * S * ZM4 - 1.00E+01 * CW2 * S2 * ZM2
      ANWWWW(3,3,4) = -1.60E+01 * CW2 * QQ0 * QQF * QQI * S * ZM
      ANWWWW(3,4,1) = 0.00E+00
      ANWWWW(3,4,2) = 0.00E+00
      ANWWWW(3,4,3) = 0.00E+00
      ANWWWW(3,4,4) = 0.00E+00
      ANWWWW(4,1,1) = 0.00E+00
      ANWWWW(4,1,2) = 0.00E+00
      ANWWWW(4,1,3) = 0.00E+00
      ANWWWW(4,1,4) = 0.00E+00
      ANWWWW(4,2,1) = 4.00E+00 * CW2 * QQF * QQI * S2
      ANWWWW(4,2,2) = -8.00E+00 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(4,2,3) = 8.00E+00 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(4,2,4) = 0.00E+00
      ANWWWW(4,3,1) = -4.00E+00 * CW2 * QQF * QQI * S2
      ANWWWW(4,3,2) = 8.00E+00 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(4,3,3) = -8.00E+00 * CW2 * QQF * QQI * S * ZM2
      ANWWWW(4,3,4) = 0.00E+00
      ANWWWW(4,4,1) = 0.00E+00
      ANWWWW(4,4,2) = 0.00E+00
      ANWWWW(4,4,3) = 0.00E+00
      ANWWWW(4,4,4) = 0.00E+00
C
      ADWWWW(1,1) = 1.00E+00
      ADWWWW(1,2) = S - 2.00E+00 * ZM2
      ADWWWW(1,3) = S - 2.00E+00 * ZM2
      ADWWWW(1,4) = 1.00E+00
      ADWWWW(2,1) = 0.00E+00
      ADWWWW(2,2) = -1.00E+00 * S + 4.00E+00 * (S / (S + 4.00E+00
     $ * QQF * QQI)) - 1.60E+01 * (ZM2 / (S + 4.00E+00 * QQF * QQI))
     $ + 4.00E+00 * ((S * ZM2) / (S + 4.00E+00 * QQF * QQI))
      ADWWWW(2,3) = S - 4.00E+00 * (S / (S + 4.00E+00 * QQF * QQI))
     $ + 1.60E+01 * (ZM2 / (S + 4.00E+00 * QQF * QQI)) - 4.00E+00
     $ * ((S * ZM2) / (S + 4.00E+00 * QQF * QQI))
      ADWWWW(2,4) = 0.00E+00
C
      AIWWWW(1) = 1.60E+01 * ((HG * HM * ZM3) / (CW * PROPH))
     $ - 8.00E+00 * ((HG * HM * S * ZM) / (CW * PROPH)) - 8.00E+00
     $ * ((HG * HM * S * ZM3) / (CW * PROPH)) + 4.00E+00 * ((HG * HM
     $ * S2 * ZM) / (CW * PROPH))
      AIWWWW(2) = 0.00E+00
      AIWWWW(3) = 1.60E+01 * ((HG * HM * ZM3) / (CW * PROPH))
     $ - 8.00E+00 * ((HG * HM * S * ZM3) / (CW * PROPH))
      AIWWWW(4) = 0.00E+00
C
C          RESTORE MISSING FACTORS
      DO 100 J=1,4
      AIWWWW(J)=AIWWWW(J)*GSQ/(16.*ZM2)
      DO 100 I=1,4
      DO 110 K=1,4
110   ANWWWW(K,I,J)=ANWWWW(K,I,J)*GSQ/(16.*ZM2)
100   CONTINUE
C
      RETURN
      END
