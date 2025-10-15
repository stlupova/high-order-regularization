#include <iostream>
#include <vector>
#include <iomanip>

#include "stokes.h" 

using namespace std;

//*********************************************************************
// evaluates Stokes single layer integral (Stokeslet) without special treatment
// f is density of SL

void stokes_SL_away(int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& f,
		    double x, double y, double z,
		    vector<double>& SL) {
  
  for (int i=0; i<3; i++) {
    SL[i] = 0.0;  
  }

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];

    double r2 = dx * dx + dy * dy + dz * dz;

    double area = Surfc[j].Area;

    double f_dot_dx = f[3*j] * dx + f[3*j+1] * dy + f[3*j+2] * dz;
    
    double H1 = 1.0 / sqrt(r2) * area;
    double H2 = H1 / r2 * f_dot_dx;
    
    SL[0] += f[3*j]   * H1 + dx * H2;
    SL[1] += f[3*j+1] * H1 + dy * H2;
    SL[2] += f[3*j+2] * H1 + dz * H2;
  }
  
  SL[0] *= denom_EightPI;
  SL[1] *= denom_EightPI;
  SL[2] *= denom_EightPI;  
}

//*********************************************************************
// evaluates Stokes single layer integral (Stokeslet) at any point
// f is density of SL
// subtraction is used
// 7th order regularization

void stokes_SL_7ord(double DEL,
		    int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& f, double f0_dot_n0,
		    double x, double y, double z, double b,
		    vector<double>& SL) {
  
  double lam = b / DEL;
  double c1, c2, c3;
  get_reg_coeff_7ord(lam, c1, c2, c3);

  for (int i=0; i<3; i++) {
    SL[i] = 0.0;  
  }

  vector<double> f_sl(3,0);

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];

    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double area = Surfc[j].Area;

    for (int k=0; k<3; k++) {
      f_sl[k] = f[3*j+k] - f0_dot_n0 * Surfc[j].Nrml[k];
    }
    double f_dot_dx = f_sl[0] * dx + f_sl[1] * dy + f_sl[2] * dz;
    
    double H1, H2, H3;
    Stokes_SL_7ord(r, DEL, c1, c2, c3, H1, H2);
    H1 *= area;
    H2 *= area * f_dot_dx;
    
    SL[0] += f_sl[0] * H1 + dx * H2;
    SL[1] += f_sl[1] * H1 + dy * H2;
    SL[2] += f_sl[2] * H1 + dz * H2;
  }
  
  SL[0] *= denom_EightPI;
  SL[1] *= denom_EightPI;
  SL[2] *= denom_EightPI;  
}

//*********************************************************************
// evaluates Stokes pressure at any point with any special treatment

double stokes_pressure_away(int N_quad, const vector<Surf_point>& Surfc,
			    const vector<double>& f,
			    double x, double y, double z) {
  
  double p = 0.0;  

  for (int j=0; j<N_quad; j++) {
    
    double dx = x - Surfc[j].x[0];
    double dy = y - Surfc[j].x[1];
    double dz = z - Surfc[j].x[2];
    
    double r2 = dx * dx + dy * dy + dz * dz;

    double H2 = 1.0 / (r2 * sqrt(r2));

    H2 *= f[3*j] * dx + f[3*j+1] * dy + f[3*j+2] * dz;
    
    p += H2 * Surfc[j].Area;
  }
    
  p *= denom_FourPI;

  return p;
}

//*********************************************************************
// evaluates Stokes pressure at any point
// f is surface force
// subtraction is used
// 7th order regularization

double stokes_pressure_7ord(double DEL,
			    int N_quad, const vector<Surf_point>& Surfc,
			    const vector<double>& f, const vector<double>& f0,
			    const vector<double>& n0,
			    double x, double y, double z, double b) {
  
  double lam = b / DEL;
  double c1, c2, c3;
  get_reg_coeff_7ord(lam, c1, c2, c3);
  
  double f0_dot_n0 = dot_product(f0, n0);
  vector<double> n0_cross_f0 = cross_product(n0, f0);
  
  vector<double> dx(3,0), n(3,0), ff(3,0), n_cross_dx(3,0), n_cross_f(3,0);

  double p = 0.0;  

  for (int j=0; j<N_quad; j++) {
    
    dx[0] = x - Surfc[j].x[0];
    dx[1] = y - Surfc[j].x[1];
    dx[2] = z - Surfc[j].x[2];
    
    double r = sqrt( dot_product(dx,dx) );

    double area = Surfc[j].Area;

    double H2;
    Lapl_DL_7ord(r, DEL, c1, c2, c3, H2);

    n[0] = Surfc[j].Nrml[0];
    n[1] = Surfc[j].Nrml[1];
    n[2] = Surfc[j].Nrml[2];
    
    ff[0] = f[3*j];
    ff[1] = f[3*j+1];
    ff[2] = f[3*j+2];

    double f_dot_n = dot_product(ff, Surfc[j].Nrml);
    //n_cross_f = cross_product(Surfc[j].Nrml, ff);
    
    n_cross_f[0] = n[1] * ff[2] - n[2] * ff[1];
    n_cross_f[1] = n[2] * ff[0] - n[0] * ff[2];
    n_cross_f[2] = n[0] * ff[1] - n[1] * ff[0];
    
    n_cross_f[0] -= n0_cross_f0[0];
    n_cross_f[1] -= n0_cross_f0[1];
    n_cross_f[2] -= n0_cross_f0[2];
    
    double n_dot_dx = dot_product(Surfc[j].Nrml, dx);
    //n_cross_dx = cross_product(Surfc[j].Nrml, dx);
    
    n_cross_dx[0] = n[1] * dx[2] - n[2] * dx[1];
    n_cross_dx[1] = n[2] * dx[0] - n[0] * dx[2];
    n_cross_dx[2] = n[0] * dx[1] - n[1] * dx[0];
    
    double normal = - (f_dot_n - f0_dot_n0) * n_dot_dx;
    double tangential = - dot_product(n_cross_f, n_cross_dx);
    
    p += ( normal + tangential ) * H2 * area;
  }
    
  p *= -denom_FourPI;

  if ( b < -1.e-12 ) { // target point is inside
    p -= f0_dot_n0;
  }

  return p;
}

