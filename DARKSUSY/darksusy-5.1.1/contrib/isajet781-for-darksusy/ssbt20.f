CDECK  ID>, SSBT20.
      DOUBLE PRECISION FUNCTION SSBT20(M1,M2)
      IMPLICIT NONE
      DOUBLE PRECISION SSA0,SSB00,M1SQ,M2SQ
      COMPLEX*16 SSB0
      REAL M1,M2
      M1SQ=M1*M1
      M2SQ=M2*M2
c      SSBT20=((SSA0(M1)+SSA0(M2))/2.D0+(M1SQ+M2SQ)
c     $*SSB00(M1,M2)+M1SQ+M2SQ)/6.D0
c     $-(SSA0(M1)+SSA0(M2))/4.D0
      SSBT20=DBLE(((SSA0(M1)+SSA0(M2))/2.D0+(M1SQ+M2SQ)
     $*SSB0(0.,M1,M2)+M1SQ+M2SQ)/6.D0
     $-(SSA0(M1)+SSA0(M2))/4.D0)
      RETURN
      END
