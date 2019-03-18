*****************************************************************************
*** auxiliary routine called by dsIBffdxdy
*** author: Torsten Bringmann, 2007-07-05
*****************************************************************************

      real*8 function dsIBffdxdy_7(x,y,m0,ml,msl1)
      implicit none
      real*8 x,y,m0,ml,msl1

      dsIBffdxdy_7 = 
     - (-64*m0**3*ml*(-ml**8 + 
     -      4*m0**2*ml**4*(msl1**2*(2 - 2*x + x**2) - 
     -         2*ml**2*(1 + x**2 - 2*y)) + 
     -      64*m0**8*y*(2*x**4 + 2*y - 4*y**3 + x**2*(6 + y) - 
     -         x**3*(5 + 2*y) + x*(-2 - 6*y + 8*y**2)) + 
     -      4*m0**4*ml**2*(4*msl1**2*
     -          (x**3 - 4*y - 2*x**2*y + x*(2 + 4*y)) + 
     -         ml**2*(2 - 10*x**3 + 16*y - 24*y**2 + 2*x*(-7 + 4*y) + 
     -            x**2*(1 + 16*y))) - 
     -      16*m0**6*(4*msl1**2*(2 - 2*x + x**2)*(x - y)*y + 
     -         ml**2*(2*x**4 - x**3*(1 + 12*y) + 
     -            4*y*(1 + 2*y - 4*y**2) + 2*x**2*(2 + 5*y + 4*y**2) + 
     -            2*x*(-1 - 10*y + 8*y**2)))))/
     -  ((3*ml**2 - 2*msl1**2 + m0**2*(-2 + 4*x - 4*y))*
     -    (ml**2 + 4*m0**2*(x - y))**2*(ml**2 - 4*m0**2*y)**2*
     -    (ml**2 - 2*msl1**2 + m0**2*(-2 + 4*y)))

      return
      end   ! dsIBffdxdy_7
