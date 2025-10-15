#include <iostream>
#include <vector>
#include <iomanip>

#include "stokeslets.h"

using namespace std;

//*********************************************************************
// evaluates stokes single layer potential (Surfc.f is density)

void eval_Stokes_SL_offSurf(double h, double DEL,
			    int N_quad, const vector<Surf_point>& Surfc,
			    int N_target, const vector<Target_point>& Target,
			    vector<double>& SL_vec) {

  initialize_vector(3*N_target, SL_vec);

  vector<double> f(3*N_quad,0);
  for (int i=0; i<N_quad; i++) {
    f[3*i  ] = Surfc[i].f[0];
    f[3*i+1] = Surfc[i].f[1];
    f[3*i+2] = Surfc[i].f[2];
  }

#pragma omp parallel 
  {
  vector<double> pt(3,0), SL(3,0), dx(3,0), f_temp(3,0);
  
#pragma omp for firstprivate(Surfc)
  
  for (int i=0; i<N_target; i++) {    

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;
    
    double b = Target[i].b;

    if ( abs(b/DEL) > 8.0 ) {

      SL[0] = 0.0;  SL[1] = 0.0;  SL[2] = 0.0;
      
      for (int j=0; j<N_quad; j++) {

	for (int k=0; k<3; k++)  dx[k] = pt[k] - Surfc[j].x[k];      

	double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	double r = sqrt(r2);
      
	f_temp[0] = Surfc[j].f[0];
	f_temp[1] = Surfc[j].f[1];
	f_temp[2] = Surfc[j].f[2];

	double f_dot_dx = dot_product(f_temp, dx);

	double area = Surfc[j].Area;
	
	double H1 = area / r;
	double H2 = f_dot_dx * area / (r2 * r);
	
	for (int k=0; k<3; k++) {
	  SL[k] += f_temp[k] * H1 + dx[k] * H2;
	}
      }// end for j
      for (int k=0; k<3; k++) {
	SL[k] *= denom_EightPI;
      }
    }
    else
    {
      double f0_dot_n0 = dot_product(Target[i].nrst_f, Target[i].nrst_Nrml);
    
      stokes_SL_7ord(DEL,
		     N_quad, Surfc,
		     f, f0_dot_n0,
		     pt[0], pt[1], pt[2], b,
		     SL);
    } // end if
    
    SL_vec[3*i  ] = SL[0];
    SL_vec[3*i+1] = SL[1];
    SL_vec[3*i+2] = SL[2];
  }
  }
}

//*********************************************************************
// evaluates Stokes single layer integral (Stokeslet) at any point
// f is density of SL
// subtraction is used in SL
// 7th order regularization

void stokes_SL_7ord(double DEL,
		    int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& f, double f0_dot_n0,
		    double x, double y, double z, double b,
		    vector<double>& SL) {
  
  double lam = b / DEL;
  double c1, c2, c3;
  get_reg_coeff_7ord(lam, c1, c2, c3);

  SL[0] = 0.0;  SL[1] = 0.0;  SL[2] = 0.0;

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];

    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double area = Surfc[j].Area;

    double f_x = f[3*j  ] - f0_dot_n0 * Surfc[j].Nrml[0];
    double f_y = f[3*j+1] - f0_dot_n0 * Surfc[j].Nrml[1];
    double f_z = f[3*j+2] - f0_dot_n0 * Surfc[j].Nrml[2];

    double f_dot_dx = f_x * dx + f_y * dy + f_z * dz;
    
    double H1, H2, H3;
    Stokes_SL_7ord(r, DEL, c1, c2, c3, H1, H2);
    H1 *= area;
    H2 *= area * f_dot_dx;
    
    SL[0] += f_x * H1 + dx * H2;
    SL[1] += f_y * H1 + dy * H2;
    SL[2] += f_z * H1 + dz * H2;
  }
  
  SL[0] *= denom_EightPI;
  SL[1] *= denom_EightPI;
  SL[2] *= denom_EightPI;  
}

//*********************************************************************
// evaluates stokes single layer potential on the surface
// Surfc.f is density

void eval_Stokes_SL_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   int N_target, const vector<Target_point>& Target,
			   vector<double>& SL_vec) {

  initialize_vector(3*N_target, SL_vec);

  vector<double> f(3*N_quad,0);
  for (int i=0; i<N_quad; i++) {
    f[3*i  ] = Surfc[i].f[0];
    f[3*i+1] = Surfc[i].f[1];
    f[3*i+2] = Surfc[i].f[2];
  }

#pragma omp parallel 
  {
  vector<double> pt(3,0), SL(3,0);
  
#pragma omp for firstprivate(Surfc)
  
  for (int i=0; i<N_target; i++) {    

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;
    
    double f0_dot_n0 = dot_product(Surfc[i].f, Surfc[i].Nrml);
    
    stokes_SL_7ord_onSurf(DEL,
			  N_quad, Surfc,
			  f, f0_dot_n0,
			  pt[0], pt[1], pt[2],
			  SL);

    SL_vec[3*i  ] = SL[0];
    SL_vec[3*i+1] = SL[1];
    SL_vec[3*i+2] = SL[2];
  }
  }
}

//*********************************************************************
// evaluates Stokes single layer integral (Stokeslet) at a surface point
// f is density of SL
// subtraction is used
// 7th order regularization

void stokes_SL_7ord_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   const vector<double>& f, double f0_dot_n0,
			   double x, double y, double z,
			   vector<double>& SL) {
  
  SL[0] = 0.0;  SL[1] = 0.0;  SL[2] = 0.0;

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];

    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double area = Surfc[j].Area;

    double f_x = f[3*j  ] - f0_dot_n0 * Surfc[j].Nrml[0];
    double f_y = f[3*j+1] - f0_dot_n0 * Surfc[j].Nrml[1];
    double f_z = f[3*j+2] - f0_dot_n0 * Surfc[j].Nrml[2];

    double f_dot_dx = f_x * dx + f_y * dy + f_z * dz;
    
    double H1, H2;
    Stokes_SL_7ord_onSurf(r, DEL, H1, H2);
    H1 *= area;
    H2 *= area * f_dot_dx;
    
    SL[0] += f_x * H1 + dx * H2;
    SL[1] += f_y * H1 + dy * H2;
    SL[2] += f_z * H1 + dz * H2;
  }
  
  SL[0] *= denom_EightPI;
  SL[1] *= denom_EightPI;
  SL[2] *= denom_EightPI;  
}

