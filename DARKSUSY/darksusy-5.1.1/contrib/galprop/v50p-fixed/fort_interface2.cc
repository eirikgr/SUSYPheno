
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * fort_interface2.cc *                          galprop package * 2001/05/11 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

// C++ wrapper for fortran routines 4/10/2000
// These routines are not to be renamed when changing platform 

using namespace std;//AWS20050624
#include"galprop_classes.h"
#include"galproph.h"

extern "C" void   yieldx_(int*,int*,int*,int*,float*,float*);  // TS code       IMOS20020502
extern "C" double wsigma_(int*,int*,int*,int*,double*);        // Webber's code IMOS20020502
extern "C" void sigtap2_(int*);                                          // IMOS20010511
extern "C" double sighad_(int*,double*,double*,double*,double*,double*); // IMOS20020502
extern "C" double
  cfactor_(int*,double*,double*,double*,double*,double*,double*,double*,double*);
extern "C" double antiproton_(int*,double*,double*,int*,int*,int*,int*); // IMOS20010511
extern "C" double synchrotron_(double*,double*,double*);
extern "C" double fjones_(double*,double*,double*);                      //IMOS20060420
extern "C" void                                                          //IMOS20060420
  aic_(int*,int*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Webber's isotopic production cross section  IMOS20020502
double wsigma_cc(int IZ, int IA, int IZF, int IAF, double E)
{
   return( wsigma_(&IZ,&IA,&IZF,&IAF,&E) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Silberberg & Tsao isotopic production cross section  IMOS20020502
double yieldx_cc(int IZ, int IA, int IZF, int IAF, float E)
{
   float CSmb;
   yieldx_(&IZ,&IA,&IZF,&IAF,&E,&CSmb);
   return( 1.*CSmb );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// Barashenkov & Polanski pA total cross section  IMOS20020502
double sighad_cc(int IS, double PA, double PZ, double TA, double TZ, double E)
{ 
   return( sighad_(&IS, &PA, &PZ, &TA, &TZ, &E) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

// initialization of the Barashenkov & Polanski cross section code
void sigtap_cc(int ISS)
{ 
   sigtap2_(&ISS);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double aic_cc(int key,int kbg,double E0,double E,double gamma,double RG,//IMOS20060420
	    double rho,double xi,double z,double RS,double DENS)
{
//   cout<<"aic_cc"<<endl;
  double SPEC;
  if(z==0.) z = 1.e-4;
  aic_(&key,&kbg,&E0,&E,&gamma,&RG,&rho,&xi,&z,&RS, &SPEC,&DENS);
  return(Pi*Rele*Rele/Mele *SPEC);
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double fjones_cc(double gam, double E1, double E4)   //IMOS20060420
{
//   cout<<"fjones_cc"<<endl;
  return ( Pi*Rele*Rele/Mele*fjones_(&gam,&E1,&E4) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double cfactor_cc(int kbg, double E0, double E, double PLindex,double RG,
		  double rho, double xi,double z,double RS)
{
//   cout<<"cfactor_cc"<<endl;
  if(z==0.) z = 1.e-4;   //IMOS20060420
  return ( cfactor_(&kbg,&E0,&E,&PLindex,&RG,&rho,&xi,&z,&RS) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int Galprop::test_cfactor()
{
   cout<<">>>>test_cfactor"<<endl;
   //   extern int isrf_energy_density_i_comp; AWS20050817 now in Galprop class
   int kbg=1; 
   double E0=1./.511e6;
   double E=1e3/.511; 
   double PLindex=3;
   double RG=15.;
   double rho=1;
   double xi=2;
   double z=3;
   double RS=8.5;
   
   for(xi= 0; xi<3; xi+=1)
      for(z= 0.1; z<4; z+=1)
         for (rho = 1; rho<15; rho+=2)
            for(int i_comp=0; i_comp<3; i_comp++)
            {
               isrf_energy_density_i_comp=i_comp;
               cout<<"i_comp="<<i_comp<<" rho="<<rho<<" xi="<<xi<<" z="<<z
                  <<" cfactor="<<cfactor_cc(kbg,E0,E,PLindex,RG,rho,xi,z,RS)
                  <<"IC_anisotropy_factor="
                  <<IC_anisotropy_factor(E,rho,xi,z,i_comp)<<endl;//IMOS20060420
            }
   cout<<"<<<<test_cfactor"<<endl;
   return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double antiproton_cc(int key,double Pap1,double Pp1,int NZ1,int NA1,int NZ2,int NA2) // IMOS20010511
{ 
//   cout<<"antiproton_cc"<<endl;
   return ( antiproton_(&key, &Pap1, &Pp1, &NZ1, &NA1, &NZ2, &NA2) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int test_antiproton_cc()
{
   cout<<">>>>test_antiproton_cc"<<endl;

   int NZ1=1,NA1=1,NZ2=1,NA2=1,key=0;         // IMOS20010511
   for(double Pap1=0.10; Pap1<1.e4; Pap1*=10) // antiproton momentum in GeV
   { 
      cout<<endl;
      for(double Pp1=1; Pp1<1.e7; Pp1*=2)
         cout<<"Pap1="<<Pap1<<" Pp1="<<Pp1<<" NZ1="<<NZ1<<" NA1="<<NA1
            <<" NZ2="<<NZ2<<" NA2="<<NA2
            <<" antiproton_cc="<<antiproton_cc(key,Pap1,Pp1,NZ1,NA1,NZ2,NA2)<<endl; // IMOS20010511
   }
   cout<<"<<<<antiproton_cc"<<endl;
   return 0;
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

double synchrotron_cc(double gamma, double nu, double B)
{ 
//   cout<<"synchrotron_cc_cc"<<endl;
   return ( synchrotron_(&gamma, &nu, &B) );
}

//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|

int test_synchrotron_cc()
{
   cout<<">>>>test_synchrotron_cc"<<endl;

   double B=1.0e-6;// Gauss 
   for( B=1.0e-6;B<1.0e-5;B*=2)
      for(double Eelectron=10;Eelectron<1.e6;Eelectron*=10) // (e-) energy in MeV
      {
         double gamma=Eelectron/m_electron;
         for(double nu=1.e6;nu<1.e12;nu*=5)
            cout<<"B="<<B<<" Eelectron="<<Eelectron<<" nu="<<nu
               <<" synchrotron_cc="<<synchrotron_cc(gamma,nu,B  )<<endl;
      }
   cout<<"<<<<test_synchrotron__cc"<<endl;
   return 0;
}
