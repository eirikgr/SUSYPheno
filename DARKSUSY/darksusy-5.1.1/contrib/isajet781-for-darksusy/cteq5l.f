CDECK  ID>, CTEQ5L.
      DOUBLE PRECISION FUNCTION CTEQ5L(IFL,X,Q)
C ----------------------------------------------------------------------
C          Parameterization of CTEQ5l parton distributions f(ifl,x,q)
C          IFL: 1=u,2=d,3=s,4=c,5=b
C               0=g
C              -1=ubar,-2=dbar,-3=sbar,-4=cbar,-5=bbar
C          Was faux5l by J. Pumplin, 9/99
C          Converted to strict Fortran 77 and Patchy by F. Paige
C ----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X,Q
      INTEGER IFL
      INTEGER NEX,NLF
      PARAMETER (NEX=8, NLF=2)
      DOUBLE PRECISION AM(0:NEX,0:NLF,-5:2)
      DOUBLE PRECISION ALFVEC(-5:2), QMAVEC(-5:2)
      DOUBLE PRECISION MEXVEC(-5:2)
      DOUBLE PRECISION UT1VEC(-5:2), UT2VEC(-5:2)
      DOUBLE PRECISION AF(0:NEX)
      DOUBLE PRECISION TMP,SB,SB1,SB2,SBX,Y,U,PART1,PART2,PART3,PART4
      INTEGER MLFVEC(-5:2)
      INTEGER I,K
C
      DATA MEXVEC( 2) / 8 /
      DATA MLFVEC( 2) / 2 /
      DATA UT1VEC( 2) /  0.4971265E+01 /
      DATA UT2VEC( 2) / -0.1105128E+01 /
      DATA ALFVEC( 2) /  0.2987216E+00 /
      DATA QMAVEC( 2) /  0.0000000E+00 /
      DATA (AM( 0,K, 2),K=0, 2)
     $ /  0.5292616E+01, -0.2751910E+01, -0.2488990E+01 /
      DATA (AM( 1,K, 2),K=0, 2)
     $ /  0.9714424E+00,  0.1011827E-01, -0.1023660E-01 /
      DATA (AM( 2,K, 2),K=0, 2)
     $ / -0.1651006E+02,  0.7959721E+01,  0.8810563E+01 /
      DATA (AM( 3,K, 2),K=0, 2)
     $ / -0.1643394E+02,  0.5892854E+01,  0.9348874E+01 /
      DATA (AM( 4,K, 2),K=0, 2)
     $ /  0.3067422E+02,  0.4235796E+01, -0.5112136E+00 /
      DATA (AM( 5,K, 2),K=0, 2)
     $ /  0.2352526E+02, -0.5305168E+01, -0.1169174E+02 /
      DATA (AM( 6,K, 2),K=0, 2)
     $ / -0.1095451E+02,  0.3006577E+01,  0.5638136E+01 /
      DATA (AM( 7,K, 2),K=0, 2)
     $ / -0.1172251E+02, -0.2183624E+01,  0.4955794E+01 /
      DATA (AM( 8,K, 2),K=0, 2)
     $ /  0.1662533E-01,  0.7622870E-02, -0.4895887E-03 /
C
      DATA MEXVEC( 1) / 8 /
      DATA MLFVEC( 1) / 2 /
      DATA UT1VEC( 1) /  0.2612618E+01 /
      DATA UT2VEC( 1) / -0.1258304E+06 /
      DATA ALFVEC( 1) /  0.3407552E+00 /
      DATA QMAVEC( 1) /  0.0000000E+00 /
      DATA (AM( 0,K, 1),K=0, 2)
     $ /  0.9905300E+00, -0.4502235E+00,  0.1624441E+00 /
      DATA (AM( 1,K, 1),K=0, 2)
     $ /  0.8867534E+00,  0.1630829E-01, -0.4049085E-01 /
      DATA (AM( 2,K, 1),K=0, 2)
     $ /  0.8547974E+00,  0.3336301E+00,  0.1371388E+00 /
      DATA (AM( 3,K, 1),K=0, 2)
     $ /  0.2941113E+00, -0.1527905E+01,  0.2331879E+00 /
      DATA (AM( 4,K, 1),K=0, 2)
     $ /  0.3384235E+02,  0.3715315E+01,  0.8276930E+00 /
      DATA (AM( 5,K, 1),K=0, 2)
     $ /  0.6230115E+01,  0.3134639E+01, -0.1729099E+01 /
      DATA (AM( 6,K, 1),K=0, 2)
     $ / -0.1186928E+01, -0.3282460E+00,  0.1052020E+00 /
      DATA (AM( 7,K, 1),K=0, 2)
     $ / -0.8545702E+01, -0.6247947E+01,  0.3692561E+01 /
      DATA (AM( 8,K, 1),K=0, 2)
     $ /  0.1724598E-01,  0.7120465E-02,  0.4003646E-04 /
C
      DATA MEXVEC( 0) / 8 /
      DATA MLFVEC( 0) / 2 /
      DATA UT1VEC( 0) / -0.4656819E+00 /
      DATA UT2VEC( 0) / -0.2742390E+03 /
      DATA ALFVEC( 0) /  0.4491863E+00 /
      DATA QMAVEC( 0) /  0.0000000E+00 /
      DATA (AM( 0,K, 0),K=0, 2)
     $ /  0.1193572E+03, -0.3886845E+01, -0.1133965E+01 /
      DATA (AM( 1,K, 0),K=0, 2)
     $ / -0.9421449E+02,  0.3995885E+01,  0.1607363E+01 /
      DATA (AM( 2,K, 0),K=0, 2)
     $ /  0.4206383E+01,  0.2485954E+00,  0.2497468E+00 /
      DATA (AM( 3,K, 0),K=0, 2)
     $ /  0.1210557E+03, -0.3015765E+01, -0.1423651E+01 /
      DATA (AM( 4,K, 0),K=0, 2)
     $ / -0.1013897E+03, -0.7113478E+00,  0.2621865E+00 /
      DATA (AM( 5,K, 0),K=0, 2)
     $ / -0.1312404E+01, -0.9297691E+00, -0.1562531E+00 /
      DATA (AM( 6,K, 0),K=0, 2)
     $ /  0.1627137E+01,  0.4954111E+00, -0.6387009E+00 /
      DATA (AM( 7,K, 0),K=0, 2)
     $ /  0.1537698E+00, -0.2487878E+00,  0.8305947E+00 /
      DATA (AM( 8,K, 0),K=0, 2)
     $ /  0.2496448E-01,  0.2457823E-02,  0.8234276E-03 /
C
      DATA MEXVEC(-1) / 8 /
      DATA MLFVEC(-1) / 2 /
      DATA UT1VEC(-1) /  0.3862583E+01 /
      DATA UT2VEC(-1) / -0.1265969E+01 /
      DATA ALFVEC(-1) /  0.2457668E+00 /
      DATA QMAVEC(-1) /  0.0000000E+00 /
      DATA (AM( 0,K,-1),K=0, 2)
     $ /  0.2647441E+02,  0.1059277E+02, -0.9176654E+00 /
      DATA (AM( 1,K,-1),K=0, 2)
     $ /  0.1990636E+01,  0.8558918E-01,  0.4248667E-01 /
      DATA (AM( 2,K,-1),K=0, 2)
     $ / -0.1476095E+02, -0.3276255E+02,  0.1558110E+01 /
      DATA (AM( 3,K,-1),K=0, 2)
     $ / -0.2966889E+01, -0.3649037E+02,  0.1195914E+01 /
      DATA (AM( 4,K,-1),K=0, 2)
     $ / -0.1000519E+03, -0.2464635E+01,  0.1964849E+00 /
      DATA (AM( 5,K,-1),K=0, 2)
     $ /  0.3718331E+02,  0.4700389E+02, -0.2772142E+01 /
      DATA (AM( 6,K,-1),K=0, 2)
     $ / -0.1872722E+02, -0.2291189E+02,  0.1089052E+01 /
      DATA (AM( 7,K,-1),K=0, 2)
     $ / -0.1628146E+02, -0.1823993E+02,  0.2537369E+01 /
      DATA (AM( 8,K,-1),K=0, 2)
     $ / -0.1156300E+01, -0.1280495E+00,  0.5153245E-01 /
C
      DATA MEXVEC(-2) / 7 /
      DATA MLFVEC(-2) / 2 /
      DATA UT1VEC(-2) /  0.1895615E+00 /
      DATA UT2VEC(-2) / -0.3069097E+01 /
      DATA ALFVEC(-2) /  0.5293999E+00 /
      DATA QMAVEC(-2) /  0.0000000E+00 /
      DATA (AM( 0,K,-2),K=0, 2)
     $ / -0.6556775E+00,  0.2490190E+00,  0.3966485E-01 /
      DATA (AM( 1,K,-2),K=0, 2)
     $ /  0.1305102E+01, -0.1188925E+00, -0.4600870E-02 /
      DATA (AM( 2,K,-2),K=0, 2)
     $ / -0.2371436E+01,  0.3566814E+00, -0.2834683E+00 /
      DATA (AM( 3,K,-2),K=0, 2)
     $ / -0.6152826E+01,  0.8339877E+00, -0.7233230E+00 /
      DATA (AM( 4,K,-2),K=0, 2)
     $ / -0.8346558E+01,  0.2892168E+01,  0.2137099E+00 /
      DATA (AM( 5,K,-2),K=0, 2)
     $ /  0.1279530E+02,  0.1021114E+00,  0.5787439E+00 /
      DATA (AM( 6,K,-2),K=0, 2)
     $ /  0.5858816E+00, -0.1940375E+01, -0.4029269E+00 /
      DATA (AM( 7,K,-2),K=0, 2)
     $ / -0.2795725E+02, -0.5263392E+00,  0.1290229E+01 /
C
      DATA MEXVEC(-3) / 7 /
      DATA MLFVEC(-3) / 2 /
      DATA UT1VEC(-3) /  0.3753257E+01 /
      DATA UT2VEC(-3) / -0.1113085E+01 /
      DATA ALFVEC(-3) /  0.3713141E+00 /
      DATA QMAVEC(-3) /  0.0000000E+00 /
      DATA (AM( 0,K,-3),K=0, 2)
     $ /  0.1580931E+01, -0.2273826E+01, -0.1822245E+01 /
      DATA (AM( 1,K,-3),K=0, 2)
     $ /  0.2702644E+01,  0.6763243E+00,  0.7231586E-02 /
      DATA (AM( 2,K,-3),K=0, 2)
     $ / -0.1857924E+02,  0.3907500E+01,  0.5850109E+01 /
      DATA (AM( 3,K,-3),K=0, 2)
     $ / -0.3044793E+02,  0.2639332E+01,  0.5566644E+01 /
      DATA (AM( 4,K,-3),K=0, 2)
     $ / -0.4258011E+01, -0.5429244E+01,  0.4418946E+00 /
      DATA (AM( 5,K,-3),K=0, 2)
     $ /  0.3465259E+02, -0.5532604E+01, -0.4904153E+01 /
      DATA (AM( 6,K,-3),K=0, 2)
     $ / -0.1658858E+02,  0.2923275E+01,  0.2266286E+01 /
      DATA (AM( 7,K,-3),K=0, 2)
     $ / -0.1149263E+02,  0.2877475E+01, -0.7999105E+00 /
C
      DATA MEXVEC(-4) / 7 /
      DATA MLFVEC(-4) / 2 /
      DATA UT1VEC(-4) /  0.4400772E+01 /
      DATA UT2VEC(-4) / -0.1356116E+01 /
      DATA ALFVEC(-4) /  0.3712017E-01 /
      DATA QMAVEC(-4) /  0.1300000E+01 /
      DATA (AM( 0,K,-4),K=0, 2)
     $ / -0.8293661E+00, -0.3982375E+01, -0.6494283E-01 /
      DATA (AM( 1,K,-4),K=0, 2)
     $ /  0.2754618E+01,  0.8338636E+00, -0.6885160E-01 /
      DATA (AM( 2,K,-4),K=0, 2)
     $ / -0.1657987E+02,  0.1439143E+02, -0.6887240E+00 /
      DATA (AM( 3,K,-4),K=0, 2)
     $ / -0.2800703E+02,  0.1535966E+02, -0.7377693E+00 /
      DATA (AM( 4,K,-4),K=0, 2)
     $ / -0.6460216E+01, -0.4783019E+01,  0.4913297E+00 /
      DATA (AM( 5,K,-4),K=0, 2)
     $ /  0.3141830E+02, -0.3178031E+02,  0.7136013E+01 /
      DATA (AM( 6,K,-4),K=0, 2)
     $ / -0.1802509E+02,  0.1862163E+02, -0.4632843E+01 /
      DATA (AM( 7,K,-4),K=0, 2)
     $ / -0.1240412E+02,  0.2565386E+02, -0.1066570E+02 /
C
      DATA MEXVEC(-5) / 6 /
      DATA MLFVEC(-5) / 2 /
      DATA UT1VEC(-5) /  0.5562568E+01 /
      DATA UT2VEC(-5) / -0.1801317E+01 /
      DATA ALFVEC(-5) /  0.4952010E-02 /
      DATA QMAVEC(-5) /  0.4500000E+01 /
      DATA (AM( 0,K,-5),K=0, 2)
     $ / -0.6031237E+01,  0.1992727E+01, -0.1076331E+01 /
      DATA (AM( 1,K,-5),K=0, 2)
     $ /  0.2933912E+01,  0.5839674E+00,  0.7509435E-01 /
      DATA (AM( 2,K,-5),K=0, 2)
     $ / -0.8284919E+01,  0.1488593E+01, -0.8251678E+00 /
      DATA (AM( 3,K,-5),K=0, 2)
     $ / -0.1925986E+02,  0.2805753E+01, -0.3015446E+01 /
      DATA (AM( 4,K,-5),K=0, 2)
     $ / -0.9480483E+01, -0.9767837E+00, -0.1165544E+01 /
      DATA (AM( 5,K,-5),K=0, 2)
     $ /  0.2193195E+02, -0.1788518E+02,  0.9460908E+01 /
      DATA (AM( 6,K,-5),K=0, 2)
     $ / -0.1327377E+02,  0.1201754E+02, -0.6277844E+01 /
C
      IF(Q.LE.QMAVEC(IFL).OR.X.GE.1D0) THEN
         CTEQ5L = 0.D0
         RETURN
      ENDIF
      TMP = LOG(Q/ALFVEC(IFL))
      IF(TMP .LE. 0.D0) THEN
         CTEQ5L = 0.D0
         RETURN
      ENDIF
      SB = LOG(TMP)
      SB1 = SB - 1.2D0
      SB2 = SB1*SB1
      DO 100 I = 0, NEX
         AF(I) = 0.D0
         SBX = 1.D0
         DO 110 K = 0, MLFVEC(IFL)
            AF(I) = AF(I) + SBX*AM(I,K,IFL)
            SBX = SB1*SBX
110      CONTINUE
100   CONTINUE
      Y = -LOG(X)
      U = LOG(X/0.00001D0)
      PART1 = AF(1)*Y**(1.D0+0.01D0*AF(4))*(1.D0+ AF(8)*U)
      PART2 = AF(0)*(1.D0 - X) + AF(3)*X 
      PART3 = X*(1.D0-X)*(AF(5)+AF(6)*(1.D0-X)+AF(7)*X*(1.D0-X))
      PART4 = UT1VEC(IFL)*LOG(1.D0-X) + 
     $            AF(2)*LOG(1.D0+EXP(UT2VEC(IFL))-X)
      CTEQ5L = EXP(LOG(X) + PART1 + PART2 + PART3 + PART4)
C          Include threshold factor...
      CTEQ5L = CTEQ5L * (1.D0 - QMAVEC(IFL)/Q)
C
      RETURN
      END
