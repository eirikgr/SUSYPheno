      
C-------------------------------------------------------
C#######################################################
C-------------------------------------------------------

      FUNCTION SQRLAM(x,y,z) 
      IMPLICIT  NONE
      REAL*8     SQRLAM,x,y,z  
      SQRLAM=sqrt(max(1.E-06,
     &       (x+(y+z))*(x-(y+z))*(x+(y-z))*(x-(y-z))))
      return
      END
