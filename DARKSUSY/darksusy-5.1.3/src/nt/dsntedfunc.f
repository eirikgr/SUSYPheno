
*************************
      real*8 function dsntedfunc(r)
      implicit none
      include 'dsmpconst.h'

      real*8 r,dsntearthdens

      dsntedfunc=dsntearthdens(r)*1000.0d0*4.0d0*pi*r**2

      return
      end
