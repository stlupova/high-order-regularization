#include <iostream>
#include <vector>
#include <iomanip>

#include "stresslets.h"

using namespace std;

//*********************************************************************

// evaluates stokes double layer potential (Surfc.g is density)

void eval_Stokes_DL_offSurf(double h, double DEL,
			    int N_quad, const vector<Surf_point>& Surfc,
			    int N_target, vector<Target_point>& Target,
			    vector<double>& DL_vec) {
  
  initialize_vector(3*N_target, DL_vec);
  
  vector<double> g(3*N_quad,0);
  for (int i=0; i<N_quad; i++) {
    g[3*i  ] = Surfc[i].g[0];
    g[3*i+1] = Surfc[i].g[1];
    g[3*i+2] = Surfc[i].g[2];
  }

#pragma omp parallel 
  {
  vector<double> pt(3,0), x0(3,0), n0(3,0), g0(3,0), DL(3,0), dx(3,0), g_temp(3,0);
  
#pragma omp for firstprivate(Surfc)

  for (int i=0; i<N_target; i++) {    

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;
        
    double b = Target[i].b;
    
    if ( abs(b)/DEL > 8.0 ) {
      
      DL[0] = 0.0;  DL[1] = 0.0;  DL[2] = 0.0;

      for (int j=0; j<N_quad; j++) {
	
	for (int k=0; k<3; k++) dx[k] = pt[k]-Surfc[j].x[k];      
	
	double r2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
	double r = sqrt(r2);

	g_temp[0] = Surfc[j].g[0]; 
	g_temp[1] = Surfc[j].g[1];  
	g_temp[2] = Surfc[j].g[2];  
	
	double g_dot_dx = dot_product(g_temp, dx);
	double n_dot_dx = dot_product(Surfc[j].Nrml, dx);
	
	double kernel = g_dot_dx * n_dot_dx /(r2*r2*r) * Surfc[j].Area;
	
	for (int k=0; k<3; k++) {
	  DL[k] += dx[k] * kernel;
	}
      }
      for (int k=0; k<3; k++) {
	DL[k] *= -6.0 * denom_EightPI;
      }
   }
    else {
      
      x0[0] = Target[i].nrst_x;
      x0[1] = Target[i].nrst_y;
      x0[2] = Target[i].nrst_z;
      
      n0[0] = Target[i].nrst_Nrml[0];
      n0[1] = Target[i].nrst_Nrml[1];
      n0[2] = Target[i].nrst_Nrml[2];
      
      g0[0] = Target[i].nrst_g[0];
      g0[1] = Target[i].nrst_g[1];
      g0[2] = Target[i].nrst_g[2];
        
      stokes_DL_7ord(DEL,
		     N_quad, Surfc,
		     g, g0,
		     x0, n0,
		     pt[0], pt[1], pt[2], b,
		     DL);
    }
    DL_vec[3*i  ] = DL[0];
    DL_vec[3*i+1] = DL[1];
    DL_vec[3*i+2] = DL[2];    
    
  } // end for i
  } // end pragma omp parallel
}

//*********************************************************************
// evaluates Stokes double layer integral (stresslet) at any point
// g is density
// subtraction is used
// 7th order regularization

void stokes_DL_7ord(double DEL,
		    int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& g, const vector<double>& g0,
		    const vector<double>& x0, const vector<double>& n0,
		    double x, double y, double z, double b,
		    vector<double>& DL) {
  
  double lam = b / DEL;
  double c1, c2, c3;
  get_reg_coeff_7ord(lam, c1, c2, c3);

  DL[0] = 0.0;  DL[1] = 0.0;  DL[2] = 0.0;

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];

    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double area = Surfc[j].Area;
    
    double H2, H3;
    Stokes_DL_7ord(r, DEL, c1, c2, c3, H2, H3);
    H2 *= area;
    H3 *= area;

    
    double x_hat0 = Surfc[j].x[0] - x0[0];
    double x_hat1 = Surfc[j].x[1] - x0[1];
    double x_hat2 = Surfc[j].x[2] - x0[2];
    
    double sig_sq = (x_hat0 * x_hat0 + x_hat1 * x_hat1 + x_hat2 * x_hat2)
      - 2.0 * b * (x_hat0 * n0[0] + x_hat1 * n0[1] + x_hat2 * n0[2]);

    double g_dl0 = g[3*j]   - g0[0];
    double g_dl1 = g[3*j+1] - g0[1];
    double g_dl2 = g[3*j+2] - g0[2];
    
    double g_dot_n0 = g_dl0 * n0[0] + g_dl1 * n0[1] + g_dl2 * n0[2];
    double g_dot_x_hat = g_dl0 * x_hat0 + g_dl1 * x_hat1 + g_dl2 * x_hat2;
    double n0_dot_n = n0[0] * Surfc[j].Nrml[0] + n0[1] * Surfc[j].Nrml[1]
      + n0[2] * Surfc[j].Nrml[2];
    double x_hat_dot_n = x_hat0 * Surfc[j].Nrml[0] + x_hat1 * Surfc[j].Nrml[1]
      + x_hat2 * Surfc[j].Nrml[2];
    
    double c_m1a = g_dot_n0 * n0_dot_n;
    double c_m1b = g_dot_x_hat * n0_dot_n  +  g_dot_n0 * x_hat_dot_n;
    
    double c_m2a = g_dot_x_hat * x_hat_dot_n;
    double c_m2b = b * (g_dot_x_hat * n0_dot_n  +  g_dot_n0 * x_hat_dot_n);

    double t1_x = dx * c_m1a - n0[0] * c_m1b;
    double t1_y = dy * c_m1a - n0[1] * c_m1b;
    double t1_z = dz * c_m1a - n0[2] * c_m1b;
    
    DL[0] += t1_x * H2;
    DL[1] += t1_y * H2;
    DL[2] += t1_z * H2;
    
    DL[0] += ( -sig_sq * t1_x  + dx * c_m2a + x_hat0 * c_m2b) * H3;
    DL[1] += ( -sig_sq * t1_y  + dy * c_m2a + x_hat1 * c_m2b) * H3;
    DL[2] += ( -sig_sq * t1_z  + dz * c_m2a + x_hat2 * c_m2b) * H3;
  }
  
  DL[0] *= -6.0 * denom_EightPI;
  DL[1] *= -6.0 * denom_EightPI;
  DL[2] *= -6.0 * denom_EightPI;

  if (( b < 0.0 ) || ( abs(b) < 1.0e-14 )) { // target point is inside or on the boundary 
    DL[0] += g0[0];
    DL[1] += g0[1];
    DL[2] += g0[2];
  }
}


//*********************************************************************
// evaluates stokes double layer potential on the surface
// Surfc.g is density

void eval_Stokes_DL_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   int N_target, vector<Target_point>& Target,
			   vector<double>& DL_vec) {
  
  initialize_vector(3*N_target, DL_vec);
  
  vector<double> g(3*N_quad,0);
  for (int i=0; i<N_quad; i++) {
    g[3*i  ] = Surfc[i].g[0];
    g[3*i+1] = Surfc[i].g[1];
    g[3*i+2] = Surfc[i].g[2];
  }

#pragma omp parallel 
  {
  vector<double> pt(3,0), g0(3,0), DL(3,0);
  
#pragma omp for firstprivate(Surfc)

  for (int i=0; i<N_target; i++) {    

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;
        
    g0[0] = Surfc[i].g[0];
    g0[1] = Surfc[i].g[1];
    g0[2] = Surfc[i].g[2];
        
    stokes_DL_7ord_onSurf(DEL,
			  N_quad, Surfc,
			  g, g0,
			  pt[0], pt[1], pt[2],
			  DL);
    
    DL_vec[3*i  ] = DL[0];
    DL_vec[3*i+1] = DL[1];
    DL_vec[3*i+2] = DL[2];    

  } // end for i
  } // end pragma omp parallel
}

//*********************************************************************
// evaluates Stokes double layer integral (stresslet) at any point
// g is density 
// subtraction is used
// 7th order regularization

void stokes_DL_7ord_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   const vector<double>& g, const vector<double>& g0,
			   double x, double y, double z,
			   vector<double>& DL) {
  
  DL[0] = 0.0;  DL[1] = 0.0;  DL[2] = 0.0;

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];

    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double area = Surfc[j].Area;
    
    double H3;
    Stokes_DL_7ord_onSurf(r, DEL, H3);
    
    double g_dl0 = g[3*j]   - g0[0];
    double g_dl1 = g[3*j+1] - g0[1];
    double g_dl2 = g[3*j+2] - g0[2];    
    
    double dx_dot_n = dx * Surfc[j].Nrml[0] + dy * Surfc[j].Nrml[1] + dz * Surfc[j].Nrml[2];
    double dx_dot_g = dx * g_dl0 + dy * g_dl1 + dz * g_dl2;

    H3 *= dx_dot_n * dx_dot_g * area;
    
    DL[0] += dx * H3;
    DL[1] += dy * H3;
    DL[2] += dz * H3;    
  }
  
  DL[0] *= -6.0 * denom_EightPI;
  DL[1] *= -6.0 * denom_EightPI;
  DL[2] *= -6.0 * denom_EightPI;

  // target point is on the boundary 
  DL[0] += 0.5 * g0[0];
  DL[1] += 0.5 * g0[1];
  DL[2] += 0.5 * g0[2];
}

