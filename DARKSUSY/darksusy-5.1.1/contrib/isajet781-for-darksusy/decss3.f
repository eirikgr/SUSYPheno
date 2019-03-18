CDECK  ID>, DECSS3.
      FUNCTION DECSS3(IP,MEA)
C
C          Compute matrix element for mode MEA of particle IP using
C          poles and couplings in /DKYSS3/.
C          Auxiliary routine for DECAY.
C
      IMPLICIT NONE
C
      COMMON/ITAPES/ITDKY,ITEVT,ITCOM,ITLIS
      SAVE /ITAPES/
      INTEGER   ITDKY,ITEVT,ITCOM,ITLIS
      INTEGER   MXPTCL,IPACK
      PARAMETER (MXPTCL=4000,IPACK=10000)
      COMMON/PARTCL/NPTCL,PPTCL(5,MXPTCL),IORIG(MXPTCL),IDENT(MXPTCL)
     1,IDCAY(MXPTCL)
      SAVE /PARTCL/
      INTEGER   NPTCL,IORIG,IDENT,IDCAY
      REAL      PPTCL
      COMMON/CONST/PI,SQRT2,ALFA,GF,UNITS
      SAVE /CONST/
      REAL      PI,SQRT2,ALFA,GF,UNITS
C
C          Data for SUSY 3-body matrix elements. There is a double 
C          pointer structure, first to modes, and then to poles that
C          make up the matrix element for that mode:
C          MELEM=-I in /DKYTAB/ points to the mode information:
C            J1SS3(I) = start of pole list for this mode
C            J2SS3(I) = end of pole list for this mode
C            WTSS3(I) = maximum weight for this mode
C          J1SS3<J<J2SS3 points to the corresponding poles:
C            KSS3(J)    = pole type
C            AMSS3(J)   = pole mass
C            ZISS3(2,J) = initial couplings
C            ZFSS3(2,J) = final couplings
C          For gaugino -> gaugino f fbar, the pole types are
C            KSS3=1: spin-1 pole in f-fbar channel
C            KSS3=2: spin-0 pole in gaugino-f channel
C            KSS3=3: spin-0 pole in gaugino-fbar channel
C            KSS3=4: spin-0 pole in f-fbar channel
C          The two couplings are the coefficients of 1,gamma_5 or of
C          gamma_mu,gamma_mu*gamma_5. 
C
      INTEGER MXMSS3,MXPSS3
      PARAMETER (MXMSS3=1000)
      PARAMETER (MXPSS3=2000)
      COMMON/DKYSS3/NMSS3,NPSS3,
     $J1SS3(MXMSS3),J2SS3(MXMSS3),WTSS3(MXMSS3),
     $KSS3(MXPSS3),AMSS3(MXPSS3),ZISS3(2,MXPSS3),ZFSS3(2,MXPSS3)
      INTEGER NMSS3,NPSS3,KSS3,J1SS3,J2SS3
      REAL WTSS3,AMSS3
      COMPLEX ZISS3,ZFSS3
C
      LOGICAL KIN(4),KINP(4)
      INTEGER IP,MEA,I,J,JP,II,PTYPE1,PTYPE2
      REAL DECSS3
      REAL AM0SQ,AM1SQ,AM2SQ,AM3SQ,S12,S13,S23
      REAL D12,D13,D23,D01,D02,D03,AS,BS,CS,DS,MSQ
      REAL DOT4
      COMPLEX A,B,C,D,AC,BC,CC,DC,AP,BP,CP,DP,APC,BPC,CPC,DPC,MMPD
C
      DOT4(I,J)=PPTCL(4,I)*PPTCL(4,J)-PPTCL(1,I)*PPTCL(1,J)-
     $PPTCL(2,I)*PPTCL(2,J)-PPTCL(3,I)*PPTCL(3,J)
C
C          Kinematics
C
      AM0SQ=PPTCL(5,IP)**2
      AM1SQ=PPTCL(5,NPTCL+1)**2
      AM2SQ=PPTCL(5,NPTCL+2)**2
      AM3SQ=PPTCL(5,NPTCL+3)**2
      D12=DOT4(NPTCL+1,NPTCL+2)
      D13=DOT4(NPTCL+1,NPTCL+3)
      D23=DOT4(NPTCL+2,NPTCL+3)
      D01=DOT4(IP,NPTCL+1)
      D02=DOT4(IP,NPTCL+2)
      D03=DOT4(IP,NPTCL+3)
      S12=2*D12+AM1SQ+AM2SQ
      S13=2*D13+AM1SQ+AM3SQ
      S23=2*D23+AM2SQ+AM3SQ
C
C          Generic matrix element
C
C          Loop over diagrams
      DECSS3=0.
      DO J=J1SS3(MEA),J2SS3(MEA)
       PTYPE1=KSS3(J)
       A=ZISS3(1,J)
       B=ZISS3(2,J)
       C=ZFSS3(1,J)
       D=ZFSS3(2,J)
       AC=CONJG(A)
       BC=CONJG(B)
       CC=CONJG(C)
       DC=CONJG(D)
       AS=A*AC
       BS=B*BC
       CS=C*CC
       DS=D*DC
       DO JP=J,J2SS3(MEA)
        MSQ=0.
        DO II=1,4
          KIN(II)=.FALSE.
          KINP(II)=.FALSE.
        END DO
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+1)).LT.AMSS3(J)) KIN(1)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+3)).LT.AMSS3(J)) KIN(2)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+2)).LT.AMSS3(J)) KIN(3)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+1)).LT.AMSS3(J)) KIN(4)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+1)).LT.AMSS3(JP)) KINP(1)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+3)).LT.AMSS3(JP)) KINP(2)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+2)).LT.AMSS3(JP)) KINP(3)=.TRUE.
        IF ((PPTCL(5,IP)-PPTCL(5,NPTCL+1)).LT.AMSS3(JP)) KINP(4)=.TRUE.
        IF (J.EQ.JP) THEN
         IF (PTYPE1.EQ.1.AND.KIN(1)) THEN
          MSQ=32*(((AS+BS)*(CS+DS)+4*REAL(A*BC)*REAL(C*DC))*D03*D12+
     $            ((AS+BS)*(CS+DS)-4*REAL(A*BC)*REAL(C*DC))*D02*D13+
     $             (BS-AS)*(CS+DS)*SQRT(AM0SQ*AM1SQ)*D23)/
     $            (S23-AMSS3(J)**2)**2
         ELSE IF (PTYPE1.EQ.2.AND.KIN(2)) THEN
          MSQ=16*(AS+BS)*(CS+DS)*D03*D12/(S12-AMSS3(J)**2)**2
         ELSE IF (PTYPE1.EQ.3.AND.KIN(3)) THEN
          MSQ=16*(AS+BS)*(CS+DS)*D02*D13/(S13-AMSS3(J)**2)**2
         ELSE IF (PTYPE1.EQ.4.AND.KIN(4)) THEN
          MSQ=16*((AS+BS)*(CS+DS)*D01*D23+(AS-BS)*(CS+DS)*D23*
     $        SQRT(AM0SQ*AM1SQ))/(S23-AMSS3(J)**2)**2
         END IF
        END IF          
        IF (J.NE.JP) THEN
        PTYPE2=KSS3(JP)
        AP=ZISS3(1,JP)
        BP=ZISS3(2,JP)
        CP=ZFSS3(1,JP)
        DP=ZFSS3(2,JP)
        APC=CONJG(AP)
        BPC=CONJG(BP)
        CPC=CONJG(CP)
        DPC=CONJG(DP)
         IF (PTYPE1.EQ.2.AND.PTYPE2.EQ.2.AND.KIN(2).AND.KINP(2)) THEN
          MMPD=16*D12*D03*(A*APC+B*BPC)*(C*CPC+D*DPC)/
     $        (S12-AMSS3(J)**2)/(S12-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.3.AND.PTYPE2.EQ.3.AND.KIN(3).AND.KINP(3)) THEN
          MMPD=16*D13*D02*(A*APC+B*BPC)*(C*CPC+D*DPC)/
     $        (S13-AMSS3(J)**2)/(S13-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.4.AND.PTYPE2.EQ.4.AND.KIN(4).AND.KINP(4)) THEN
          MMPD=16*D23*(D01*(A*APC+B*BPC)*(C*CPC+D*DPC)+
     $        SQRT(AM0SQ*AM1SQ)*(A*APC-B*BPC)*(C*CPC-D*DPC))/
     $        (S23-AMSS3(J)**2)/(S23-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.1.AND.PTYPE2.EQ.3.AND.KIN(1).AND.KINP(3)) THEN
          MMPD=(16*D13*D02*((A*C-B*D)*(-APC*CPC+BPC*DPC)+
     $         (A*D-B*C)*(APC*DPC-BPC*CPC))+
     $         8*D23*SQRT(AM0SQ*AM1SQ)*((A*C+B*D)*(APC*CPC-BPC*DPC)-
     $         (A*D+B*C)*(APC*DPC-BPC*CPC)))/
     $         (S23-AMSS3(J)**2)/(S13-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.1.AND.PTYPE2.EQ.2.AND.KIN(1).AND.KINP(2)) THEN
          MMPD=(16*D12*D03*((A*C+B*D)*(-APC*CPC+BPC*DPC)+
     $         (A*D+B*C)*(APC*DPC-BPC*CPC))+
     $         8*D23*SQRT(AM0SQ*AM1SQ)*((A*C-B*D)*(APC*CPC-BPC*DPC)+
     $         (-A*D+B*C)*(APC*DPC-BPC*CPC)))/
     $         (S23-AMSS3(J)**2)/(S12-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.3.AND.PTYPE2.EQ.4.AND.KIN(3).AND.KINP(4)) THEN
          MMPD=((8*D13*D23+4*D23*AM1SQ)*((A*C+B*D)*(APC*CPC+BPC*DPC)+
     $         (A*D+B*C)*(APC*DPC+BPC*CPC))+
     $         4*D23*SQRT(AM0SQ*AM1SQ)*((A*C+B*D)*(APC*CPC-BPC*DPC)+
     $         (A*D+B*C)*(APC*DPC-BPC*CPC)))/
     $         (S13-AMSS3(J)**2)/(S23-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.2.AND.PTYPE2.EQ.4.AND.KIN(2).AND.KINP(4)) THEN
          MMPD=-((8*D12*D23+4*D23*AM1SQ)*((A*C+B*D)*(APC*CPC+BPC*DPC)+
     $         (A*D+B*C)*(APC*DPC+BPC*CPC))+
     $         4*D23*SQRT(AM0SQ*AM1SQ)*((A*C+B*D)*(APC*CPC-BPC*DPC)+
     $         (A*D+B*C)*(APC*DPC-BPC*CPC)))/
     $         (S12-AMSS3(J)**2)/(S23-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
         IF (PTYPE1.EQ.2.AND.PTYPE2.EQ.3.AND.KIN(2).AND.KINP(3)) THEN
          MMPD=((8*D12*D13-4*D23*AM1SQ)*((A*C+B*D)*(APC*CPC+BPC*DPC)+
     $         (A*D+B*C)*(APC*DPC+BPC*CPC))-
     $         4*D23*SQRT(AM0SQ*AM1SQ)*((A*C-B*D)*(APC*CPC-BPC*DPC)+
     $         (A*D-B*C)*(APC*DPC-BPC*CPC)))/
     $         (S12-AMSS3(J)**2)/(S13-AMSS3(JP)**2)
          MSQ=2*REAL(MMPD)
         END IF
        END IF
        DECSS3=DECSS3+MSQ
       END DO
      END DO
C
      RETURN
      END
