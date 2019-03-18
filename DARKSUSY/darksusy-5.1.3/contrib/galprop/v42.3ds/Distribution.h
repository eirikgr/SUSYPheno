
//**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
// * Distribution.h *                              galprop package * 4/14/2000 
//**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

#ifndef Distribution_h
#define Distribution_h

#include<iostream.h>
#include<math.h>
#include<stdlib.h> //IMOS20020112
#include<string.h> //IMOS20020112
#include"Spectrum.h"

class Distribution
{
 public:

   int n_pgrid;          // number of points in momentum
   int n_rgrid;          // number of points in radius (2D)      
   int n_zgrid;          // number of points in z (1D,2D)  
   int n_xgrid;          // number of points in x (3D)
   int n_ygrid;          // number of points in y (3D)  

   int n_spatial_dimensions;// 1,2,3D

   Spectrum   *d1; // 1D array
   Spectrum  **d2; // 2D array
   Spectrum ***d3; // 3D array

   Distribution(){};   // needed for initial construction
//   ~Distribution();   // needed for destruction
   Distribution(int n_rgrid_,int n_zgrid_,int n_pgrid_);
   Distribution(int n_xgrid_,int n_ygrid_,int n_zgrid_,int n_pgrid_);
   void delete_array();
   void init(int n_rgrid_,int n_zgrid_,int n_pgrid_);
   void init(int n_xgrid_,int n_ygrid_,int n_zgrid_,int n_pgrid_);
   void print();
   void print(double units);
   float max();
   Distribution operator=(double);
   Distribution operator+(double);
   Distribution operator+=(double);
   Distribution operator*(double);
   Distribution operator*=(double);
   Distribution operator/=(double);
   Distribution operator =(Distribution);
   Distribution operator+=(Distribution);
   Distribution operator*=(Distribution);
   Distribution operator/=(Distribution);
};

//  following cases have to be non-member functions since they have 2 arguments;
//  therefore their prototypes are defined outside the class

Distribution operator+(double,Distribution);
Distribution operator*(double,Distribution);
Distribution operator+(Distribution,Distribution);
Distribution operator*(Distribution,Distribution);
Distribution operator/(Distribution,Distribution);

#endif




