#include <iostream>
#include <vector>

#include "utilities.h"

using namespace std;

//*********************************************************************

double phi(const vector<double>& y) {

  double val;
  
  if ( surfc_elliptical ) {
    
    vector<double> x(3,0);
    x[0] = y[0] - ctr;
    x[1] = y[1] - ctr;
    x[2] = y[2] - ctr;
    
    vector<double> a2(3,0);
    a2[0] = ellipse_a * ellipse_a;
    a2[1] = ellipse_b * ellipse_b;
    a2[2] = ellipse_c * ellipse_c;
    
    val = x[0] * x[0] / a2[0]
      + x[1] * x[1] / a2[1]
      + x[2] * x[2] / a2[2] - 1.0;
  }

  else { // molecular surface
    
    vector<double> x(3,0);
    x[0] = y[0] - ctr;
    x[1] = y[1] - ctr;
    x[2] = y[2] - ctr;
    
    vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
    x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
    x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
    x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
    x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
    
    double r2 = 0.5*0.5;
    double c = 0.6;
    
    double a1, a2, a3, a4;
    a1 = (x[0]-x1[0])*(x[0]-x1[0]) + (x[1]-x1[1])*(x[1]-x1[1]) + (x[2]-x1[2])*(x[2]-x1[2]); 
    a2 = (x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]); 
    a3 = (x[0]-x3[0])*(x[0]-x3[0]) + (x[1]-x3[1])*(x[1]-x3[1]) + (x[2]-x3[2])*(x[2]-x3[2]); 
    a4 = (x[0]-x4[0])*(x[0]-x4[0]) + (x[1]-x4[1])*(x[1]-x4[1]) + (x[2]-x4[2])*(x[2]-x4[2]); 
    
    val = c - exp(-a1/r2) - exp(-a2/r2) - exp(-a3/r2) - exp(-a4/r2);
  }
  return val;
}

//*********************************************************************
// D_phi/D_x(i): i-th derivative of phi

double Dphi(int i, const vector<double>& y) {

  double val;
  
  if ( surfc_elliptical ) {

    vector<double> a2(3,0);
    a2[0] = ellipse_a * ellipse_a;
    a2[1] = ellipse_b * ellipse_b;
    a2[2] = ellipse_c * ellipse_c;
    
    val = 2.0 * (y[i]-ctr) / a2[i];
  }
  
  else { // molecular surface
    
    vector<double> x(3,0);
    x[0] = y[0] - ctr;
    x[1] = y[1] - ctr;
    x[2] = y[2] - ctr;
    
    vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
    x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
    x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
    x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
    x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
    
    double r2 = 0.5*0.5;
    
    double a1, a2, a3, a4;
    a1 = (x[0]-x1[0])*(x[0]-x1[0]) + (x[1]-x1[1])*(x[1]-x1[1]) + (x[2]-x1[2])*(x[2]-x1[2]); 
    a2 = (x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]); 
    a3 = (x[0]-x3[0])*(x[0]-x3[0]) + (x[1]-x3[1])*(x[1]-x3[1]) + (x[2]-x3[2])*(x[2]-x3[2]); 
    a4 = (x[0]-x4[0])*(x[0]-x4[0]) + (x[1]-x4[1])*(x[1]-x4[1]) + (x[2]-x4[2])*(x[2]-x4[2]); 
    
    double d1, d2, d3, d4;
    d1 = x[i] - x1[i];
    d2 = x[i] - x2[i];
    d3 = x[i] - x3[i];
    d4 = x[i] - x4[i];
    
    val = 2.0/r2* (exp(-a1/r2)*d1 + exp(-a2/r2)*d2 + exp(-a3/r2)*d3 + exp(-a4/r2)*d4);  
  }
  return val;
}

//*********************************************************************
// Second derivatives of phi

void D2phi(const vector<double>& y, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33) {
  
  if ( surfc_elliptical ) {

    vector<double> a2(3,0);
    a2[0] = ellipse_a * ellipse_a;
    a2[1] = ellipse_b * ellipse_b;
    a2[2] = ellipse_c * ellipse_c;
    
    phi11 = 2.0/a2[0]; phi12 = 0.0;       phi13 = 0.0;
    phi21 = 0.0;       phi22 = 2.0/a2[1]; phi23 = 0.0;
    phi31 = 0.0;       phi32 = 0.0;       phi33 = 2.0/a2[2];
  }
  
  else { // molecular surface
    
    vector<double> x(3,0);
    x[0] = y[0] - ctr;
    x[1] = y[1] - ctr;
    x[2] = y[2] - ctr;
    
    vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
    x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
    x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
    x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
    x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
    
    double r2 = 0.5*0.5;
    
    vector<double> d1 = x; d1[0]-=x1[0]; d1[1]-=x1[1]; d1[2]-=x1[2];
    vector<double> d2 = x; d2[0]-=x2[0]; d2[1]-=x2[1]; d2[2]-=x2[2];
    vector<double> d3 = x; d3[0]-=x3[0]; d3[1]-=x3[1]; d3[2]-=x3[2];
    vector<double> d4 = x; d4[0]-=x4[0]; d4[1]-=x4[1]; d4[2]-=x4[2];
    double a1, a2, a3, a4;
    a1 = exp(-dot_product(d1,d1)/r2);
    a2 = exp(-dot_product(d2,d2)/r2);
    a3 = exp(-dot_product(d3,d3)/r2); 
    a4 = exp(-dot_product(d4,d4)/r2); 
    
    phi11 = 2.0/r2*(a1*(1.0-2.0/r2*d1[0]*d1[0]) + a2*(1.0-2.0/r2*d2[0]*d2[0]) +
		    a3*(1.0-2.0/r2*d3[0]*d3[0]) + a4*(1.0-2.0/r2*d4[0]*d4[0]));
    phi22 = 2.0/r2*(a1*(1.0-2.0/r2*d1[1]*d1[1]) + a2*(1.0-2.0/r2*d2[1]*d2[1]) +
		    a3*(1.0-2.0/r2*d3[1]*d3[1]) + a4*(1.0-2.0/r2*d4[1]*d4[1]));
    phi33 = 2.0/r2*(a1*(1.0-2.0/r2*d1[2]*d1[2]) + a2*(1.0-2.0/r2*d2[2]*d2[2]) +
		    a3*(1.0-2.0/r2*d3[2]*d3[2]) + a4*(1.0-2.0/r2*d4[2]*d4[2]));
    
    phi12 = -4.0/(r2*r2)*(a1*d1[0]*d1[1] + a2*d2[0]*d2[1] +
			  a3*d3[0]*d3[1] + a4*d4[0]*d4[1]);
    phi13 = -4.0/(r2*r2)*(a1*d1[0]*d1[2] + a2*d2[0]*d2[2] +
			  a3*d3[0]*d3[2] + a4*d4[0]*d4[2]);
    phi23 = -4.0/(r2*r2)*(a1*d1[1]*d1[2] + a2*d2[1]*d2[2] +
			  a3*d3[1]*d3[2] + a4*d4[1]*d4[2]);
    
    phi21 = phi12;
    phi31 = phi13;
    phi32 = phi23;
  }
}

//*********************************************************************
//*********************************************************************

void get_reg_coeff_7ord(double lam, double& c1, double& c2, double& c3) {

  double lam2 = lam * lam;
  double e2 = exp(lam2);

  double I0, I2, I4;
  evaluate_In_024(lam, I0, I2, I4);
  
  c3 = e2 * (2.0 * I0 - 4.0 * I2 + I4) / 8.0;
  c2 = e2 * (I0-I2) + (4.0 * lam2 + 7.0) * c3;
  c1 = 2.0 * e2 * I0 + 2.0 * (lam2 + 1.0) * c2
    - (4.0*lam2*lam2 + 6.0*lam2 + 6.0) * c3;
}

//*********************************************************************

void Lapl_DL_7ord(double r, double d,
		  double c1, double c2, double c3,
		  double& s2) {
  if (r < 1e-14) {
    s2 = 256.0 / 15.0 / rootPI / (d * d * d);
  }
  else if ( r/d > 8.0) {
    s2 = 1.0 / (r * r * r);
  }
  else {
    double rd = r / d;
    double rd2 = rd * rd;
    double rd4 = rd2 * rd2;
    double errfn = erf(rd);
    double r_expfn = rd * exp(-rd2);
    
    double poly2 = -2.0/rootPI + 2.0 * (c1 + 2.0 * (c2 + c3)) * rd2
      - 4.0 * (c2 + 5.0 * c3) * rd4 + 8.0 * c3 * rd4 * rd2;
    
    s2 = ( errfn + poly2 * r_expfn ) / (r * r * r);
  }
}

//*********************************************************************

void Stokes_SL_7ord(double r, double d,
		    double c1, double c2, double c3,
		    double& s1, double& s2) {
  if (r < 1e-14) {
    s1 = 32.0 / 5.0 / rootPI / d;
    s2 = 256.0 / 15.0 / rootPI / (d * d * d);
  }
  else if ( r/d > 8.0) {
    double r2 = r * r;
    s1 = 1.0 / r;
    s2 = 1.0 / (r2 * r);
  }
  else {
    double rd = r / d;
    double rd2 = rd * rd;
    double rd4 = rd2 * rd2;
    double rd6 = rd4 * rd2;
    double errfn = erf(rd);
    double r_expfn = rd * exp(-rd2);
    
    double poly1 = c1 - 2.0*(c2+c3) * rd2 + 4.0*c3 * rd4;
    double poly2 = -2.0/rootPI + 2.0 * (c1 + 2.0 * (c2 + c3)) * rd2
      - 4.0 * (c2 + 5.0 * c3) * rd4 + 8.0 * c3 * rd6;
    
    double r2 = r * r;
    s1 = ( errfn + poly1 * r_expfn ) / r;
    s2 = ( errfn + poly2 * r_expfn ) / (r2 * r);
  }
}

//*********************************************************************

void evaluate_In_024(double lam, double& I0, double& I2, double& I4) {

  double lam2 = lam * lam;
  double em2 = exp(-lam2);
  double t = abs(lam) * erfc(abs(lam));

  I0 = em2/rootPI - t;
  I2 = 2.0/3.0 * ((0.5-lam2) * em2/rootPI + lam2 * t);
  I4 = 8.0/15.0 * ((0.75 - 0.5*lam2 + lam2 * lam2) * em2/rootPI
		   - lam2 * lam2 * t);
}

//*********************************************************************
//*********************************************************************
// Returns surface force at a given point
// Spheroid with semi-axes (ellipse_a, ellipse_b, ellipse_b)
// Assumes ellipse_a > ellipse_b
// Spheroid centered at (ctr,ctr,ctr)
// translates with velocity (1,0,0)

vector<double> stokeslet_density(const vector<double>& pt) {
  vector<double> f(3,0);

  double x = pt[0] - ctr;
  double y = pt[1] - ctr;
  double z = pt[2] - ctr;
  
  double a_sq = ellipse_a * ellipse_a;
  double b_sq = ellipse_b * ellipse_b;
  
  double x_sq = x * x;
  double e_sq = 1.0 - b_sq / a_sq;
  double e = sqrt(e_sq);
  double Le = log( (1.0 + e) / (1.0 - e) );
  f[0] = 4.0 * e_sq * e * ellipse_a / ellipse_b
    / ( (1.0 + e_sq) * Le - 2.0 * e) /sqrt(a_sq - e_sq * x_sq);
  
  return f;
}

//*********************************************************************
// Returns exact pressure value at a given point
// Spheroid with semi-axes (ellipse_a, ellipse_b, ellipse_b)
// Assumes ellipse_a > ellipse_b
// Spheroid centered at (ctr,ctr,ctr)
// translates with velocity (1,0,0)

double p_exact(double x_in, double y_in, double z_in) {

  double p_ex = 0.0;

  double x = x_in - ctr;
  double y = y_in - ctr;
  double z = z_in - ctr;
  
  double a_sq = ellipse_a * ellipse_a;
  double b_sq = ellipse_b * ellipse_b;

  double r_sq = y * y + z * z;
  double r = sqrt(r_sq);
  
  double c = sqrt(a_sq - b_sq);
  double R1 = sqrt((x+c)*(x+c)+r_sq);
  double R2 = sqrt((x-c)*(x-c)+r_sq);
  
  double dB1dx = 1.0/R1 - 1.0/R2;
    
  double e = c/ellipse_a;
  double Le = log((1.0+e)/(1.0-e));
  double a1 = e*e/((1.0+e*e)*Le-2.0*e);

  p_ex = - 2.0 * a1 * dB1dx; // (-) eq'n (29a) in Chwang & Wu'75

  return p_ex;
}
      
//*********************************************************************
// Returns the x-derivative of pressure given above

double dp_dx_exact(double x_in, double y_in, double z_in) {

  double dp_dx_ex = 0.0;

  double x = x_in - ctr;
  double y = y_in - ctr;
  double z = z_in - ctr;
  
  double a_sq = ellipse_a * ellipse_a;
  double b_sq = ellipse_b * ellipse_b;

  double r_sq = y * y + z * z;
  double r = sqrt(r_sq);
  
  double c = sqrt(a_sq - b_sq);
  double R1 = sqrt((x+c)*(x+c)+r_sq);
  double R2 = sqrt((x-c)*(x-c)+r_sq);
  
  double e = c/ellipse_a;
  double Le = log((1.0+e)/(1.0-e));
  double a1 = e*e/((1.0+e*e)*Le-2.0*e);

  dp_dx_ex = - 2.0 * a1 * (-(x+c)/pow(R1,3) + (x-c)/pow(R2,3));
  
  return dp_dx_ex;
}
      
//*********************************************************************
// Returns exact velocity (as Stokeslet) 
// Spheroid with semi-axes (ellipse_a, ellipse_b, ellipse_b)
// Assumes ellipse_a > ellipse_b
// Spheroid centered at (ctr,ctr,ctr)
// translates with velocity (1,0,0)

vector<double> SL_exact(double x_in, double y_in, double z_in) {

  vector<double> SL_ex(3,0);

  double x = x_in - ctr;
  double y = y_in - ctr;
  double z = z_in - ctr;
  
  double a_sq = ellipse_a * ellipse_a;
  double b_sq = ellipse_b * ellipse_b;

  double r_sq = y * y + z * z;
  double r = sqrt(r_sq);
  
  double c = sqrt(a_sq - b_sq);
  double R1 = sqrt((x+c)*(x+c)+r_sq);
  double R2 = sqrt((x-c)*(x-c)+r_sq);
  
  double R1c = R1-(x+c);
  double R2c = R2-(x-c);
  
  double B10;
  if (r_sq < 1.e-14) {
    if (x<0) B10 = log( (x-c) / (x+c) );
    else B10 = log( (x+c) / (x-c) );
  }
  else B10 = log(R2c/R1c);
  
  double B30 = (x+c)/R1 - (x-c)/R2;
  
  double dB1dx = 1.0/R1 - 1.0/R2;
  
  double dB1dyz = 0.0;
  if (r_sq > 1.e-14) {
    dB1dyz = 1.0/R2/R2c - 1.0/R1/R1c;
  }
  double dB3dyz = 1.0/R2 - 1.0/R1 + x * dB1dyz;
  
  double e = c/ellipse_a;
  double Le = log((1.0+e)/(1.0-e));
  double a1 = e*e/((1.0+e*e)*Le-2.0*e);
  double b1 = a1*(1.0-e*e)/(2.0*e*e);
  
  SL_ex[0] = 2.0*a1*B10 - a1*B30 + 2.0*b1
    *( (x-c)/R2 - (x+c)/R1 + B10 + x * dB1dx );
  
  double Er = a1 * (1.0/R2-1.0/R1) + 2.0 * b1 * dB3dyz;
  SL_ex[1] = Er * y;
  SL_ex[2] = Er * z;

  return SL_ex;
}
      
//*********************************************************************
// Signed distance between two points in 3D

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal) {
  
  return (x[0]-y[0])*normal[0] + (x[1]-y[1])*normal[1] + (x[2]-y[2])*normal[2];
}

//*********************************************************************

double dot_product(const vector<double>& x,
		   const vector<double>& y) {
  
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//*********************************************************************

vector<double> cross_product(const vector<double>& x,
			     const vector<double>& y) {
  
  vector<double> cross_prod(3,0);
  cross_prod[0] = x[1]*y[2] - x[2]*y[1];
  cross_prod[1] = x[2]*y[0] - x[0]*y[2];
  cross_prod[2] = x[0]*y[1] - x[1]*y[0];
  return cross_prod;
}

//*********************************************************************
// assumes neither vector is 0 length

void orthogonalize(vector<double>& x,
		   vector<double>& y) {

  double x_dot_x = dot_product(x,x);
  double x_dot_y = dot_product(x,y);

  double coeff = x_dot_y / x_dot_x;

  y[0] -= coeff * x[0];
  y[1] -= coeff * x[1];
  y[2] -= coeff * x[2];

  double y_dot_y = dot_product(y,y);

  double norm_x = sqrt(x_dot_x);
  double norm_y = sqrt(y_dot_y);

  for (int i=0; i<3; i++) {
    x[i] = x[i] / norm_x;
    y[i] = y[i] / norm_y;
  }      
}

//*****************************************************************************
// Solve linear system through LU decomposition with pivoting.                
//                                                                           
// Factorize PA = LU with pivoting:                                           
//   The lower and upper triangular matrices are still stored in the original 
//   matrix and the permutation matrix "P" is stored in the vector "int *p".  
//******************************************************************************
 
int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p) {
  int n = A.size();

  for (int j=0; j<n; j++) p[j] = j;

  for (int j=0; j<n; j++) {
    int k = j;
    double m = A[p[j]][j];
    // Search for maximum in this column
    for (int i=j+1; i<n; i++) {
      if (abs(A[p[i]][j]) > abs(m)) {
        m = A[p[i]][j];
        k = i;
      }
    }

    // "Swap" maximum row with current row (using the permutation vector)
    if (k != j) {
      int temp = p[j];
      p[j] = p[k];
      p[k] = temp;
    }

    double ajj = A[p[j]][j];
    if (abs(ajj) < 1.0e-15) {
      return (-1);
    }

    for (int i=j+1; i<n; i++) {
      double lij = A[p[i]][j] / ajj;
      A[p[i]][j] = lij; // lower triangular elements
      for (int k=j+1; k<n; k++) {
        A[p[i]][k] -= lij * A[p[j]][k]; // upper triangular elements
      }
    }
  }
  return 0;
}

//*****************************************************************************

vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b) {
  int n = A.size();
  vector<double> x(n,0);

  // Solve Ly=b by forward substitution
  x[0] = b[p[0]];
  for (int i=1; i<n; i++) {
    x[i] = b[p[i]];
    double rowsum = 0.0;
    for (int j=0; j<i; j++) {
      rowsum += A[p[i]][j] * x[j];
    }
    x[i] -= rowsum;
  }

  // Solve Ux=y by back substitution
  x[n-1] = x[n-1] / A[p[n-1]][n-1];
  for (int i=n-2; i>=0; i--) {
    double rowsum = 0.0;
    for (int j = n - 1; j > i; j--) {
      rowsum += A[p[i]][j] * x[j];
    }    
    x[i] = (x[i] - rowsum) / A[p[i]][i];
  }

  return x;
}

//*****************************************************************************

void initialize_vector(int n, vector<double>& vec) {
  for (int i=0; i<n; i++)  vec[i] = 0.0;
}
