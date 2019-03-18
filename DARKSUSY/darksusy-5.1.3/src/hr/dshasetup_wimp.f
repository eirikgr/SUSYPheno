*****************************************************************************
***   subroutine dshasetup_wimp prepares the common blocks for a basic 
***   generic WIMP that does not annihilate to Higgs. It sets:
***     - annihilation channel branching fractions,
***     - WIMP mass
***     - WIMP annihilation cross-section
***   Keys to entries in bfs:
***    Ch No  Particles                 
***    -----  ---------                 
***      1     S1 S1                  
***      2     S1 S2                  
***      3     S2 S2         
***      4     S3 S3              
***      5     S1 S3                  
***      6     S2 S3       
***      7     S- S+              
***      8     S1 Z              
***      9     S2 Z                
***     10     S3 Z	              
***     11     W- S+ and W+ S-     
***     12     Z0 Z0 	        
***     13     W+ W-            
***     14     nu_e nu_e-bar  
***     15     e+ e-             
***     16     nu_mu nu_mu-bar   
***     17     mu+ mu-             
***     18     nu_tau nu_tau-bar
***     19     tau+ tau-	     
***     20     u u-bar           
***     21     d d-bar         
***     22     c c-bar        
***     23     s s-bar         
***     24     t t-bar       
***     25     b b-bar       
***     26     gluon gluon    
***     27     q q gluon (not implemented yet, put to zero)
***     28     gamma gamma (1-loop)
***     29     Z0 gamma (1-loop)
***
***   Note: This routine needs to be called for each each WIMP model before
***   dshaloyield is called.
***
*** Author: Pat Scott
*** Date: 14-10-14
*****************************************************************************


      subroutine dshasetup_wimp(mchi, sigchi, bfs)
      implicit none
      include 'dshacom.h'
      include 'dsprep.h'
      
      double precision, intent(IN) :: mchi, sigchi, bfs(29)
      integer i

      haib='none'
      dshasetupcalled=.true.
      hasv=sigchi
      hamwimp=mchi

      if (abs(sum(bfs) - 1.d0) - 1.d2*epsilon(mchi) .gt. 0.d0 ) stop 'Error in dshasetup_wimp: BFs do not sum to 1.'

      do i=1,29
        habr(i) = bfs(i)
        sigv(i) = hasv * habr(i)
      enddo

      end



















