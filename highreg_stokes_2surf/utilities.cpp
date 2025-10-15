#include <iostream>
#include <vector>

#include "utilities.h"

using namespace std;

//*********************************************************************

static const double alpha = 0.0;
static const double beta = 0.0;
static const double gam = 30.0 * PI/180.0;

// Generate the un-rotation matrix
// the inverse of rotation matrix is its transpose
vector< vector<double> > make_rot_mat() {

  double ca = cos(alpha);
  double sa = sin(alpha);
  double cb = cos(beta);
  double sb = sin(beta);
  double cg = cos(gam);
  double sg = sin(gam);
  
  vector< vector<double> > rot_mat(3, vector<double>(3,0));
  rot_mat[0][0] = cb*cg;
  rot_mat[0][1] = cb*sg;
  rot_mat[0][2] = -sb;

  rot_mat[1][0] = sa*sb*cg - ca*sg;
  rot_mat[1][1] = sa*sb*sg + ca*cg;
  rot_mat[1][2] = sa*cb;

  rot_mat[2][0] = ca*sb*cg + sa*sg;
  rot_mat[2][1] = ca*sb*sg - sa*cg;
  rot_mat[2][2] = ca*cb;

  return rot_mat;
}
static const vector< vector<double> > R = make_rot_mat();

//*********************************************************************

vector<double> unrotate(const vector<double>& x) {
  vector<double> x_unrt(3,0);

  // note: R is the "un-rotation" matrix"
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      x_unrt[i] += R[i][j] * x[j];
    }
  }
  return x_unrt;
}

//*********************************************************************

double phi(int surf_num, int rotated, const vector<double>& x) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;

  vector<double> xx = x;
  if (surf_num==1) {
    xx[0] = x[0]-ctr1_x;
    xx[1] = x[1]-ctr1_y;
    xx[2] = x[2]-ctr1_z;
  }
  if (surf_num==2) {
    xx[0] = x[0]-ctr2_x;
    xx[1] = x[1]-ctr2_y;
    xx[2] = x[2]-ctr2_z;
  }

  if ( rotated == 1) {
    vector<double> x_unrt = unrotate(xx); //un-rotated (x,y,z)
    
    return x_unrt[0] * x_unrt[0] / a2[0]
      + x_unrt[1] * x_unrt[1] / a2[1]
      + x_unrt[2] * x_unrt[2] / a2[2] - 1.0; 
  }
  else {
    return xx[0]*xx[0]/a2[0] + xx[1]*xx[1]/a2[1] + xx[2]*xx[2]/a2[2] - 1.0;
  }
}

//*********************************************************************
// D_phi/D_x(i): i-th derivative of phi

double Dphi(int surf_num, int rotated, int i, const vector<double>& x) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;
  
  vector<double> xx = x;
  if (surf_num==1) {
    xx[0] = x[0]-ctr1_x;
    xx[1] = x[1]-ctr1_y;
    xx[2] = x[2]-ctr1_z;
  }
  if (surf_num==2) {
    xx[0] = x[0]-ctr2_x;
    xx[1] = x[1]-ctr2_y;
    xx[2] = x[2]-ctr2_z;
  }

  if ( rotated == 1 ) {
    vector<double> x_unrt = unrotate(xx); //un-rotated (x,y,z)
    vector<double> Dphi_unrt(3,0);
    Dphi_unrt[0] = 2.0 * x_unrt[0] / a2[0]; //un-rotated grad_phi
    Dphi_unrt[1] = 2.0 * x_unrt[1] / a2[1];
    Dphi_unrt[2] = 2.0 * x_unrt[2] / a2[2];
    
    return Dphi_unrt[0] * R[0][i] + Dphi_unrt[1] * R[1][i] + Dphi_unrt[2] * R[2][i];    
  }
  else {
    return 2.0*xx[i]/a2[i];
  }
}

//*********************************************************************
// Second derivatives of phi

void D2phi(int surf_num, int rotated, const vector<double>& x,
	   double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23,
	   double& phi31, double& phi32, double& phi33) {
  
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;

  if ( rotated == 1 ) {
    double A0 = 2.0/a2[0];
    double A1 = 2.0/a2[1];
    double A2 = 2.0/a2[2];
    
    phi11 = A0 * R[0][0]*R[0][0] + A1 * R[1][0]*R[1][0] + A2 * R[2][0]*R[2][0]; 
    phi12 = A0 * R[0][0]*R[0][1] + A1 * R[1][0]*R[1][1] + A2 * R[2][0]*R[2][1]; 
    phi13 = A0 * R[0][0]*R[0][2] + A1 * R[1][0]*R[1][2] + A2 * R[2][0]*R[2][2];
    
    phi21 = phi12;
    phi22 = A0 * R[0][1]*R[0][1] + A1 * R[1][1]*R[1][1] + A2 * R[2][1]*R[2][1];
    phi23 = A0 * R[0][1]*R[0][2] + A1 * R[1][1]*R[1][2] + A2 * R[2][1]*R[2][2];
    
    phi31 = phi13;
    phi32 = phi23;
    phi33 = A0 * R[0][2]*R[0][2] + A1 * R[1][2]*R[1][2] + A2 * R[2][2]*R[2][2];  
  }
  else {
    phi11 = 2.0/a2[0]; phi12 = 0.0;       phi13 = 0.0;
    phi21 = 0.0;       phi22 = 2.0/a2[1]; phi23 = 0.0;
    phi31 = 0.0;       phi32 = 0.0;       phi33 = 2.0/a2[2];
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

void Stokes_SL_7ord_onSurf(double r, double d,
			   double& s1, double& s2) {
  if (r < 1e-14) {
    // 7th order
    
    s1 = 32.0 / 5.0 / rootPI / d;
    s2 = 256.0 / 15.0 / rootPI / (d * d * d);
    
    // 5th order
    /*
    s1 = 16.0 / 3.0 / rootPI / d;
    s2 = 32.0 / 3.0 / rootPI / (d * d * d);
    */
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
    double r_expfn = 2.0/rootPI * rd * exp(-rd2);
    
    // 7th order
    
    double poly1 = 11.0/5.0 - 26.0/15.0 * rd2 + 4.0/15.0 * rd4;
    double poly2 = -1.0 + (118.0 * rd2
			   - 68.0 * rd4 + 8.0 * rd6) / 15.0;
    
    // 5th order
    /*
    double poly1 = 5.0/3.0 - 2.0/3.0 * rd2;
    double poly2 = -1.0 + (14.0 * rd2 - 4.0 * rd4) / 3.0;
    */
    
    double r2 = r * r;
    s1 = ( errfn + poly1 * r_expfn ) / r;
    s2 = ( errfn + poly2 * r_expfn ) / (r2 * r);
  }
}

//*********************************************************************

void Stokes_DL_7ord(double r, double d,
		    double c1, double c2, double c3,
		    double& s2, double& s3) {
  if (r < 1e-14) {
    s2 = 256.0 / 15.0 / rootPI / (d * d * d);
    s3 = 0.0;
  }
  else if ( r/d > 8.0) {
    double r2 = r * r;
    s2 = 1.0 / (r2 * r);
    s3 = 1.0 / (r2 * r2 * r);
  }
  else {
    double rd = r / d;
    double rd2 = rd * rd;
    double rd4 = rd2 * rd2;
    double rd6 = rd4 * rd2;
    double errfn = erf(rd);
    double r_expfn = rd * exp(-rd2);
    
    double poly2 = -2.0/rootPI + 2.0 * (c1 + 2.0 * (c2 + c3)) * rd2
      - 4.0 * (c2 + 5.0 * c3) * rd4 + 8.0 * c3 * rd6;
    double poly3 = -2.0/rootPI * (1.0 + 2.0/3.0 * rd2) // original s3
      + 4.0/3.0 * ((c1 + 4.0 * c2 + 12.0 * c3) * rd4
		   - 2.0 * (c2 + 9.0 * c3) * rd6 + 4.0 * c3 * rd6 * rd2);
    
    double r2 = r * r;
    s2 = ( errfn + poly2 * r_expfn ) / (r2 * r);
    s3 = ( errfn + poly3 * r_expfn ) / (r2 * r2 * r);
  }
}

//*********************************************************************

void Stokes_DL_7ord_onSurf(double r, double d,
			   double& s3) {
  if (r < 1e-14) {
    // 7th order
    s3 = 128.0 / 45.0 / rootPI / (d * d * d * d * d);
    // 5th order
    //s3 = 128.0 / 15.0 / rootPI / (d * d * d * d * d);
  }
  else if ( r/d > 8.0) {
    double r2 = r * r;
    s3 = 1.0 / (r2 * r2 * r);
  }
  else {
    double rd = r / d;
    double rd2 = rd * rd;
    double rd4 = rd2 * rd2;
    double rd6 = rd4 * rd2;
    double errfn = erf(rd);
    double r_expfn = 2.0/rootPI * rd * exp(-rd2);
    
    // 7th order
    
    double poly3 = -1.0 - 2.0/3.0 * rd2
      + 52.0/45.0 * rd4 - 8.0/45.0 * rd6;
    
    // 5th order
    //double poly3 = -1.0 - 2.0/3.0 * rd2 + 4.0/9.0 * rd4;
    
    double r2 = r * r;
    s3 = ( errfn + poly3 * r_expfn ) / (r2 * r2 * r);
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
// Signed distance between two points in 3D

double distance(const vector<double>& x, const vector<double>& y,
		const vector<double>& normal) {
  return (x[0]-y[0])*normal[0] + (x[1]-y[1])*normal[1] + (x[2]-y[2])*normal[2];
}

//*********************************************************************

double dot_product(const vector<double>& x, const vector<double>& y) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//*********************************************************************
// Solve linear system through LU decomposition with pivoting.                
//                                                                           
// Factorize PA = LU with pivoting:                                           
//   The lower and upper triangular matrices are still stored in the original 
//   matrix and the permutation matrix "P" is stored in the vector "int *p".  
 
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

//*********************************************************************

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

//*********************************************************************

void initialize_vector(int n, vector<double>& vec) {
  for (int i=0; i<n; i++)  vec[i] = 0.0;
}

