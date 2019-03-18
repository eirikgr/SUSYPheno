*         -*- mode: fortran -*-
*######################################################################*
*                       i n c l u d e     f i l e                      *
*######################################################################*

************************************************************************
***                        dsmpconst.h                               ***
***         this piece of code is needed as a separate file          ***
c----------------------------------------------------------------------c
c mathematical and physical constants, p. gondolo 2011-11-03, p. scott 2014-10-20

      real*8 m_p,m_n,m_d,n_avogadro,c_light,pi,gev2cm3s,fermiGeV,gev2cm2,atomicmassunit,m_p_amu,m_n_amu
      double precision secperyr, cmperpc, GeVperSolarMass, SolarMass, mperkpc 
      double precision CGRAV, omegamh2, zeq, Mhzeq, t_0, H_0, rho_cdm, rho_conh2, kappa
      double precision t_eq, dmfrac, omegacdmh2, omegabh2, keq, M_MW

      parameter(pi=4.0d0*datan(1.d0))            !pi
      parameter(m_p_amu=1.00727646688d0)         !proton mass [amu]
      parameter(m_n_amu=1.0086649156d0)          !neutron mass [amu]
      parameter(atomicmassunit=0.931494028d0)    !atomic mass unit [GeV/c^2]
      parameter(m_p=m_p_amu*atomicmassunit)      !m_p=0.938271998d0; proton mass [GeV/c^2]
      parameter(m_n=m_n_amu*atomicmassunit)      !m_n=0.9396d0; neutron mass [GeV/c^2]
      parameter(m_d=1.875612762d0)               !deuteron mass [GeV/c^2]
      parameter(n_avogadro=6.022d23)             !Avogradros number [number / mol ]
      parameter(c_light=299792.458d0)            !speed of light [km/s]
      parameter(gev2cm3s=0.38937966d-27*3.d10)   !conversion [GeV^2 cm^3/s]
      parameter(fermiGeV=1.d0/0.1973269602d0)    !conversion [GeV fm]
      parameter(gev2cm2=(197.327053d-16)**2)     !conversion [GeV^2 cm^2]

      parameter (secperyr = 365.25d0*24.d0*60.d0*60.d0) !seconds in a yr (average)
      parameter (cmperpc = 3.0857d18)            !cm per parsec
      parameter (mperkpc = 10.d0*cmperpc)        !meters per kpc
      parameter (GeVperSolarMass = 1.1157472757094956d57)! GeV per solar mass
      parameter (SolarMass = 1.98892d33)         !g per solar mass
      parameter (CGRAV = 6.67259d-8)             !Newton's constant in cgs units   
      parameter (M_MW = 0.94d12)                 !Milky Way Mass (in solar masses) 

      !From best fit Planck+WMAP (Planck 2013):
      parameter (t_0 = secperyr * 13.8242d9)     !Current age of the universe
      parameter (H_0 = 67.04d0)                  !Current Hubble constant 
      parameter (omegamh2 = 0.14305d0)           !Present matter content of the universe * h^2
      parameter (omegabh2 = 0.022032d0)          !Present baryon content of the universe * h^2
      parameter (omegacdmh2 = 0.12038d0)         !Present CDM content of the universe * h^2

      parameter (zeq = 2.32d4*omegamh2 - 1.d0)   !redshift of matter-radiation equality (Kolb & Turner)
      parameter (Mhzeq = 6.5d15/omegamh2**2)     !Horizon mass at zeq, in solar masses (Josan, Green & Malik 2009, PRD 79:103520)
      parameter (t_eq = secperyr * 59.073d3)     !Age at matter-radiation equality (Wright 2006, PASP 118:1711)
      parameter (keq = 0.07185*omegamh2)         !Mode entering horizon at equality (derived via Friedman Eq, in Mpc^-1)
      parameter (dmfrac = omegacdmh2/(omegacdmh2+omegabh2)) !Percentage of matter that is in DM
      parameter (rho_conh2 =3.d4/(8.d0*pi*CGRAV*mperkpc**2))!Critical density today / h^2 (g cm^-3)
      parameter (rho_cdm = rho_conh2*omegacdmh2) !Cosmological CDM density today (g cm^-3)
      parameter (kappa = 16.d0*pi*CGRAV*rho_cdm*(2.32d4*omegamh2)**3*t_eq**2/(3.d0*dmfrac)) !Dimensionless compound quantity

***                                                                 ***
*********************** end of dsmpconst.h ****************************
