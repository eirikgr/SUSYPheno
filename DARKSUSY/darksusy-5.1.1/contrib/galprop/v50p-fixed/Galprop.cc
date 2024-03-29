using namespace std;

#include"galprop_classes.h"
#include"galproph.h"

///////////////////////////////////////////////////////
int Galprop::run(int argc, char*argv[])  
{
   char version[]="50p"; //AWS 20060309

   cout<<" >>>> galprop version "<<version<<endl;

   if(argc!=2){cout<<" no galdef file specified!"<<endl;return -1;}
   cout<<argc<<" "<<argv[0]<<" "<<argv[1]<<endl;

   if(configure.init() !=0)return 1;
   cout<<configure.galdef_directory<<endl;// just a test

   if(galdef.read  (version,argv[1],configure.galdef_directory) !=0) return 1;

   if(galdef.test_suite>0) { test_suite(); return 0; }

//reading all isotopic cross sections, nuclear reaction network, nuclear fits
   read_nucdata();

//initialization of the Barashenkov & Polanski cross section code 
   sigtap_cc(-1);             // IMOS20010511 AWS20010620

//initialization of the Webber's routine
   set_sigma_cc();

   if(create_galaxy() !=0) return 1;
   if(create_gcr()    !=0) return 1;

//major routine
   if(propagate_particles() !=0) return 1;

//deleting all cross section arrays etc.
   cleanup_nucdata();

//   if(print_BC() !=0) return 1;
//   if(store_gcr() !=0) return 1; //IMOS20030129

   if(galdef.output_gcr_full!=0) if(store_gcr_full() !=0) return 1;

   if(galdef.synchrotron)
   {
      if(galdef.primary_electrons || galdef.secondary_electrons                       // IMOS20050912
	 || galdef.secondary_positrons || galdef.DM_positrons || galdef.DM_electrons)
      {
         if (  gen_synch_emiss()  !=0)  return 1;
         if (  gen_synch_skymap() !=0)  return 1;
         if (store_synch_skymap() !=0)  return 1;
      } //electrons
   }//galdef.synchrotron

   if(galdef.gamma_rays>=1)                          //AWS20050302
   {
     if(galdef.pi0_decay) // IMOS20060420
       { 
	 if (  gen_pi0_decay_emiss()     !=0)  return 1;
	 if (store_pi0_decay_emiss()     !=0)  return 1;
	 if (  gen_pi0_decay_skymap()    !=0)  return 1;
	 if (store_pi0_decay_skymap()    !=0)  return 1;
	 
	 if(galdef.gamma_rays==2)                       //AWS20050302
	   {
	     if (store_pi0_decay_HIR_skymap()!=0)  return 1;//AWS20041215
	     if (store_pi0_decay_H2R_skymap()!=0)  return 1;//AWS20041215
	   }
       }
      if(galdef.primary_electrons || galdef.secondary_electrons                       // IMOS20050912
	 || galdef.secondary_positrons || galdef.DM_positrons || galdef.DM_electrons)
      {
	if(galdef.bremss)                               // IMOS20060420
	  {
	    if (  gen_bremss_emiss()          !=0)  return 1;
	    if (store_bremss_emiss()          !=0)  return 1;
	    if (  gen_bremss_skymap        () !=0)  return 1;
	    if (store_bremss_skymap        () !=0)  return 1;
	    
	    if(galdef.gamma_rays==2)                    //AWS20050302
	      {
		if (store_bremss_HIR_skymap()     !=0)  return 1;//AWS20041215
		if (store_bremss_H2R_skymap()     !=0)  return 1;//AWS20041215
	      }
	    if (  gen_bremss_ionized_skymap() !=0)  return 1;
	    if (store_bremss_ionized_skymap() !=0)  return 1;

	  }

	if(galdef.IC_isotropic+galdef.IC_anisotropic) // IMOS20060420
	  {
	    if (               gen_IC_emiss() !=0)  return 1;
	    if (               gen_IC_skymap()!=0)  return 1;
	    if (store_IC_skymap("isotropic")!=0)  return 1;
	    if(galdef.IC_anisotropic) if(store_IC_skymap("anisotropic")!=0) return 1; // AWS20010206
            
	    if (store_IC_skymap_comp("isotropic")!=0)  return 1;
	    if(galdef.IC_anisotropic) if(store_IC_skymap_comp("anisotropic")!=0) return 1;
	  }
      } //electrons
   } //galdef.gamma_rays

   if(galdef.DM_gammas)                           // IMOS20050912
     {
       if (  gen_DM_emiss ()          !=0)  return 1;
       if (store_DM_emiss ()          !=0)  return 1;
       if (  gen_DM_skymap()          !=0)  return 1;
       if (store_DM_skymap()          !=0)  return 1;
     }
   
   if(galdef.ionization_rate) // IMOS20060420
     {
       if (  gen_ionization_rate()!=0)  return 1;
       if (store_ionization_rate()!=0)  return 1;
     }
   cout<<" completed processing of galdef_"<<version<<"_"<<argv[1]<<endl;
   cout<<" <<<< galprop"<<endl;
   cout << '\a';                            // IMOS20010511
   
   return 0;
};
