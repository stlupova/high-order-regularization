#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <omp.h>

#include "stresslets.h"

using namespace std;

void test_Integral_Equation_gmres(double h,
				  int N_quad1, vector<Surf_point>& Surfc1,
				  int N_quad2, vector<Surf_point>& Surfc2,
				  int N_target1, vector<Target_point>& Target1,
				  int N_target2, vector<Target_point>& Target2);

void eval_SL(double h, double DEL_on, double DEL_off,
	     int N_quad1, vector<Surf_point>& Surfc1,
	     int N_quad2, vector<Surf_point>& Surfc2,
	     int N_target1, vector<Target_point>& Target1,
	     int N_target2, vector<Target_point>& Target2,
	     vector<double>& SL_1,
	     vector<double>& SL_2);

void solve_IE_1step(double h, double DEL_on, double DEL_off,
		    int N_quad1, vector<Surf_point>& Surfc1,
		    int N_quad2, vector<Surf_point>& Surfc2,
		    int N_target1, vector<Target_point>& Target1,
		    int N_target2, vector<Target_point>& Target2,
		    int n, const vector<double>& u_in,
		    vector<double>& u_out);


void interpolate_DL_density(int surf_num, double h, int N_quad,
			    const vector<Surf_point>& Surfc,
			    int N_target,
			    vector<Target_point>& Target);


void sixPoint_stencil(int j_start, int j_end,
		      int ind1, int ind2, int ind3,
		      double h, const vector<Surf_point>& Surfc,
		      double z0, double z_ctr, double tol,
		      double xLow, double yLow,
		      double xHigh, double yHigh,
		      vector< vector<int> >& stcl);
void sixPoint_coeff(double p, vector<double>& coeff);
double sixPointInterpolate(vector< vector<double> >& f, double p, double q);


void eightPoint_stencil(int j_start, int j_end,
			int ind1, int ind2, int ind3,
			double h, const vector<Surf_point>& Surfc,
			double z0, double z_ctr, double tol,
			double xLow, double yLow,
			double xHigh, double yHigh,
			vector< vector<int> >& stcl);
void eightPoint_coeff(double p, vector<double>& coeff);
double eightPointInterpolate(vector< vector<double> >& f, double p, double q);

void set_surface_force(int surf_num, int N_quad,
		       vector<Surf_point>& Surfc);

//*********************************************************************

int main(int argc, char** argv) {
  int N = 32;
  double L = 1.0; //bounding box: [-L,L]^3
  double h = L / N; //grid spacing

  double tm = omp_get_wtime();
  
  int surf_number = 0; // don't translate at creation
  int rotated = 0;     // don't rotate the first surface
  vector<Surf_point> Surfc1;
  Generate_Surface(surf_number, rotated, N, h, &Surfc1);
  int N_quad1 = Surfc1.size(); 
  cout << "Number of quadrature points surface ONE: " << N_quad1 << endl;  

  surf_number = 0;  // don't translate at creation
  rotated = 1;      // rotate the second surface
  vector<Surf_point> Surfc2; 
  Generate_Surface(surf_number, rotated, N, h, &Surfc2);
  int N_quad2 = Surfc2.size(); 
  cout << "Number of quadrature points surface TWO: " << N_quad2 << endl;  

  // translate both surfaces by ctr1, ctr2
  for (int i=0; i<N_quad1; i++) {
    Surfc1[i].x[0] += ctr1_x;
    Surfc1[i].x[1] += ctr1_y;
    Surfc1[i].x[2] += ctr1_z;
  }
  for (int i=0; i<N_quad2; i++) {
    Surfc2[i].x[0] += ctr2_x;
    Surfc2[i].x[1] += ctr2_y;
    Surfc2[i].x[2] += ctr2_z;
  }

  /*
  ofstream surf_file1("Spheroid1_N32_epsS");
  ofstream surf_file2("Spheroid2_N32_epsS");
  
  for (int i=0; i<N_quad1; i++) {
    surf_file1 << setprecision(16) << Surfc1[i].x[0] << " "
	      << Surfc1[i].x[1] << " " << Surfc1[i].x[2] << endl;
  }
  for (int i=0; i<N_quad2; i++) {
    surf_file2 << setprecision(16) << Surfc2[i].x[0] << " "
	      << Surfc2[i].x[1] << " " << Surfc2[i].x[2] << endl;
  }
  surf_file1.close();
  surf_file2.close(); 
  */
  
  //an array of target points with coordinates, and velocity values
  vector<Target_point> Target1, Target2; 
  Generate_Targets_OnSurface(N_quad1, Surfc1, &Target1);
  Generate_Targets_OnSurface(N_quad2, Surfc2, &Target2);

  int N_target1 = Target1.size();
  cout << "Number of points target ONE: " << N_target1 << endl;
  int N_target2 = Target2.size();
  cout << "Number of points target TWO: " << N_target2 << endl;

  test_Integral_Equation_gmres(h,
			       N_quad1, Surfc1,
			       N_quad2, Surfc2,
			       N_target1, Target1,
			       N_target2, Target2);
  
  tm = omp_get_wtime() - tm;
  cout << "CPU time elapsed: " << tm << " seconds" << endl << endl;
   
  return 0;
}

//*********************************************************************

void test_Integral_Equation_gmres(double h,
				  int N_quad1, vector<Surf_point>& Surfc1,
				  int N_quad2, vector<Surf_point>& Surfc2,
				  int N_target1, vector<Target_point>& Target1,
				  int N_target2, vector<Target_point>& Target2) {
  // 7th order
  
  double kappa0 = 4.0;
  double q = 5.0/7.0;
  
  // 5th order
  /*
  double kappa0 = 3.0;
  double q = 4.0/5.0;
  */
  double kappa = kappa0 * pow(1.0/64.0, 1.0-q);
  double DEL_off = kappa * pow(h, q);
  double DEL_on = DEL_off;

  cout << "kappa0 = " << kappa0 << "   DEL = "
		      << kappa << " * h^(" << q << ")" << endl;

  
  double tm = omp_get_wtime();
  
  Compute_nearest_data(1, 0, h, DEL_off, N_quad1, Surfc1, N_target2, Target2);
  Compute_nearest_data(2, 1, h, DEL_off, N_quad2, Surfc2, N_target1, Target1);
  
  tm = omp_get_wtime() - tm;      
  cout << "Nearest data computed in " << tm << " seconds" << endl;

  tm = omp_get_wtime();
  
  set_surface_force(1, N_quad1, Surfc1);
  set_surface_force(2, N_quad2, Surfc2);
  
  tm = omp_get_wtime() - tm;      
  cout << "Surface forces computed in " << tm << " seconds" << endl;

  
  double factor = 2.0/3.0;

  tm = omp_get_wtime();
  
  vector<double> SL_1(3*N_target1,0), SL_2(3*N_target2,0);
  
  eval_SL(h, DEL_on, DEL_off,
	  N_quad1, Surfc1,
	  N_quad2, Surfc2,
	  N_target1, Target1,
	  N_target2, Target2,
	  SL_1,
	  SL_2);  
 
  tm = omp_get_wtime() - tm;      
  cout << "SL evaluation CPU time = " << tm << " seconds" << endl << endl;
  


  int n = 3*N_target1 + 3*N_target2;
  vector<double> b(n,0);
  
  for (int i=0; i<3*N_target1; i++) {    
    b[              i]   = -SL_1[i] * factor;
  }
  for (int i=0; i<3*N_target2; i++) {    
    b[3*N_target1 + i]   = -SL_2[i] * factor;
  }
  
  
  int iter_max = 100;
  double tol = 1.e-10;

  int im = 15;

  double **v, *g, *y, **H, *c, *s, *resi;
  
  v = new double*[n];
  for (int k=0; k<n; k++)   v[k] = new double[im+1];
  
  g = new double[im+1];
  y = new double[im+1];
  H = new double*[im+1];
  for (int k=0; k<im+1; k++)   H[k] = new double[im];
  
  c = new double[im];
  s = new double[im];
  resi = new double[iter_max];

  int iter = 0;
  double error = 1.0;
  
  vector<double> x(n,0), Ax(n,0);
  initialize_vector(n, x);
  initialize_vector(n, Ax);
        
  while (iter <= iter_max) {
    
    //   STEP 1   
    if (error <= tol  ||  iter == iter_max) {
      cout << "GMRES converged in " << iter
	   << " iterations with residual = " << error << endl;

      for (int k=0; k<n; k++)   delete v[k];   delete v;
      delete g;   delete y;   delete c;   delete s;   delete resi;
      for (int k=0; k<im+1; k++)   delete H[k];   delete H;

      for (int i=0; i<N_target1; i++) {
	Target1[i].u = x[3*i];
	Target1[i].v = x[3*i+1];
	Target1[i].w = x[3*i+2];
      }
      for (int i=0; i<N_target2; i++) {
	Target2[i].u = x[3*N_target1 + 3*i];
	Target2[i].v = x[3*N_target1 + 3*i+1];
	Target2[i].w = x[3*N_target1 + 3*i+2];
      }
      
      ofstream out_file1("vel1_N32_highreg_57_k4");
      ofstream out_file2("vel2_N32_highreg_57_k4");
      
      for (int i=0; i<N_target1; i++) {
	out_file1 << setprecision(16) << Target1[i].u << " "
		  << Target1[i].v << " " << Target1[i].w << endl;
      }
      for (int i=0; i<N_target2; i++) {
	out_file2 << setprecision(16) << Target2[i].u << " "
		  << Target2[i].v << " " << Target2[i].w << endl;
      }
      out_file1.close();
      out_file2.close(); 
      
      return;
    }
    
    tm = omp_get_wtime();
    
    // compute A*x
    if (iter > 0) { //x=0 so Ax=0 at iteration 0, so skip the unnecessary calculation
      solve_IE_1step(h, DEL_on, DEL_off,
		     N_quad1, Surfc1,
		     N_quad2, Surfc2,
		     N_target1, Target1,
		     N_target2, Target2,
		     n, x, Ax);
    }
    for (int i=0; i<n; i++)  v[i][0] = b[i];// - Ax[i]; // r = b-Ax

    double r_norm = 0.0;    // r0 = norm(v(1:n,1),2);
    double b_norm = 0.0;
    for (int i=0; i<n; i++) {
      r_norm += v[i][0] * v[i][0];
      b_norm += b[i] * b[i];
    }
    r_norm = sqrt(r_norm);
    b_norm = sqrt(b_norm);
    error = r_norm;
    if (b_norm > 1.e-8)  error /= b_norm;


    tm = omp_get_wtime() - tm;      
    cout << "   step1    " << "  residual = " << error << "   tolerance = "
    	   << tol << "   CPU time = " << tm << " seconds" << endl;   

    double temp = 1.0 / r_norm;
    for (int i=0; i<n; i++)     v[i][0] = v[i][0] * temp;
    
    //   STEP 2
    
    g[0] = r_norm;
    
    int j = -1;
    
    while ((j < im-1) && (error > tol) && (iter < iter_max)) {
      
      j = j+1;
      iter = iter + 1;
      
      double tm = omp_get_wtime();
      
      // START Arnoldi's method
      //    Arnoldi's method uses the Gram-Schmidt method for computing
      //    an l2-orthonormal basis {v1,v2,...,vk} of the Krylov subspace
      //    K_k = span{v1,A*v1,...,A^(k-1)*v1}
      
      // perform Ax=A*v_j
      for (int i=0; i<n; i++)  x[i] = v[i][j];
      for (int i=0; i<n; i++)  Ax[i] = 0.0;

      solve_IE_1step(h, DEL_on, DEL_off,
		     N_quad1, Surfc1,
		     N_quad2, Surfc2,
		     N_target1, Target1,
		     N_target2, Target2,
		     n, x, Ax);
      // end of Ax
      
      for (int i=0; i<n; i++)  v[i][j+1] = Ax[i]; // Krylov vector v_j+1 = A*v_j
      
      // Modified Gram-Schmidt, keeping the Hessenberg matrix
      for (int i=0; i<=j; i++) {   
	double temp = 0.0;
	for (int k=0; k<n; k++)   temp += v[k][i] * v[k][j+1];
	H[i][j] = temp;  // (A*v_j,v_i)
	for (int k=0; k<n; k++)   v[k][j+1] -= temp * v[k][i];
      }
     
      double v_norm = 0.0;    // norm(v_j+1);
      for (int i=0; i<n; i++)   v_norm += v[i][j+1] * v[i][j+1];
      v_norm = sqrt(v_norm);
      
      H[j+1][j] = v_norm;
      if (v_norm != 0.0) {
	v_norm = 1.0 / v_norm;
	for (int k=0; k<n; k++)   v[k][j+1] *= v_norm;
      }
      // END Arnoldi method      
      
      // Update the QR factorization of bar(Hm) using plane rotations, i.e.,
      //   Apply Omega_i, i = 1,...,j-1, to the j-th column of bar(Hm)
      //   Perform previous transformations of j-th column of bar(Hm)
      
      if (j != 0) {
	for (int i=0; i<=j-1; i++) {
	  double temp = H[i][j];
	  H[i]  [j] =  c[i] * temp + s[i] * H[i+1][j];
	  H[i+1][j] = -s[i] * temp + c[i] * H[i+1][j];
	}
      }

      // Compute the plane rotation coefficients cm, sm
      
      double hh = sqrt(H[j][j] * H[j][j] + H[j+1][j] * H[j+1][j]);
      
      if (hh == 0.0)  hh = 1.e-16;
      
      c[j] = H[j][j] / hh;
      s[j] = H[j+1][j] / hh;
            
      // Premultiply the matrix by the new rotation matrix Omega_i
      H[j][j] = c[j] * H[j][j] + s[j] * H[j+1][j];

      // Eliminate the last row of matrix H(j+1,:)
      for (int i=0; i<=j; i++)  H[j+1][i] = 0.0;
      
      // Update rhs of the Least Squares problem, multiplying it by Omega_i
      g[j+1] = -s[j] * g[j];
      g[j] = c[j] * g[j];
    
      // Determine the residual norm from the last element of the rhs
      error = abs(g[j+1]);
      if (b_norm > 1.e-8)  error /= b_norm;
      resi[iter-1] = error;

      tm = omp_get_wtime() - tm;      

      cout << "Iteration " << j+1 << "   residual = " << error << "  "
         << tol << "   CPU time = " << tm << " seconds" << endl;      
            
      // Compute solution:
      //   solve upper triangular system by 'deleting' the last row of matrix
      //   Rm and right-hand side gm
      for (int i=0; i<=j; i++)   y[i] = g[i];
      
      y[j] = y[j] / H[j][j];
      for (int k=j-1; k>=0; k--) {
	double temp = y[k];
	for (int i=k+1; i<=j; i++)   temp -= H[k][i] * y[i];
	y[k] = temp / H[k][k];
      }
      
      for (int i=0; i<n; i++)   x[i] = 0.0;
      
      for (int i=0; i<=j; i++)
	for (int k=0; k<n; k++)   x[k] += y[i] * v[k][i];           
    }
  }
}

//*********************************************************************

void eval_SL(double h, double DEL_on, double DEL_off,
	     int N_quad1, vector<Surf_point>& Surfc1,
	     int N_quad2, vector<Surf_point>& Surfc2,
	     int N_target1, vector<Target_point>& Target1,
	     int N_target2, vector<Target_point>& Target2,
	     vector<double>& SL_1,
	     vector<double>& SL_2) {

  for (int i=0; i<3*N_target1; i++)  SL_1[i] = 0.0;    
  for (int i=0; i<3*N_target2; i++)  SL_2[i] = 0.0;

  double tm = omp_get_wtime();
  
  // Single layer using the first surface

  vector<double> SL_onSurf1(3*N_target1,0), SL_comp1(3*N_target2,0);

  if (use_tree == 1) { // use tree to speed up calculations
    eval_Stokes_SL_onSurf_tree(1, h, DEL_on,
			       N_quad1, Surfc1,
			       N_target1, Target1,
			       SL_onSurf1);
  }
  else {
    eval_Stokes_SL_onSurf(DEL_on,
			  N_quad1, Surfc1,
			  N_target1, Target1,
			  SL_onSurf1);
    
  }
  tm = omp_get_wtime() - tm;
  cout << "Single layer evaluated Surfc1->Target1  with CPU time = "
       << tm << " seconds" << endl;
  
  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    
    eval_Stokes_SL_offSurf_tree(1, h, DEL_off,
				N_quad1, Surfc1,
				N_target2, Target2,
				SL_comp1);  
    
  }
  else {
    
    eval_Stokes_SL_offSurf(h, DEL_off,
			   N_quad1, Surfc1,
			   N_target2, Target2,
			   SL_comp1);
    
  }

  tm = omp_get_wtime() - tm;
  cout << "Single layer evaluated Surfc1->Target2  with CPU time = "
       << tm << " seconds" << endl;
  
  
  for (int i=0; i<3*N_target1; i++) { 
    SL_1[i] = SL_onSurf1[i];    
  }
  for (int i=0; i<3*N_target2; i++) { 
    SL_2[i] = SL_comp1[i];
  }
  
  
  // Single layer using the second surface

  vector<double> SL_onSurf2(3*N_target2,0), SL_comp2(3*N_target1,0);

  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    
    eval_Stokes_SL_offSurf_tree(2, h, DEL_off,
				N_quad2, Surfc2,
				N_target1, Target1,
				SL_comp2);  
    
  }
  else {

    eval_Stokes_SL_offSurf(h, DEL_off,
			   N_quad2, Surfc2,
			   N_target1, Target1,
			   SL_comp2);
    
  }
  
  tm = omp_get_wtime() - tm;
  cout << "Single layer evaluated Surfc2->Target1  with CPU time = "
       << tm << " seconds" << endl;
  
  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    eval_Stokes_SL_onSurf_tree(2, h, DEL_on,
			       N_quad2, Surfc2,
			       N_target2, Target2,
			       SL_onSurf2);
  }
  else {
    eval_Stokes_SL_onSurf(DEL_on,
			  N_quad2, Surfc2,
			  N_target2, Target2,
			  SL_onSurf2);
    
  }
  tm = omp_get_wtime() - tm;
  cout << "Single layer evaluated Surfc2->Target2  with CPU time = "
       << tm << " seconds" << endl;
  
  
  for (int i=0; i<3*N_target1; i++) {
    SL_1[i] += SL_comp2[i];   
  }
  for (int i=0; i<3*N_target2; i++) {
    SL_2[i] += SL_onSurf2[i];
  }
}

//*********************************************************************

void solve_IE_1step(double h, double DEL_on, double DEL_off,
		    int N_quad1, vector<Surf_point>& Surfc1,
		    int N_quad2, vector<Surf_point>& Surfc2,
		    int N_target1, vector<Target_point>& Target1,
		    int N_target2, vector<Target_point>& Target2,
		    int n, const vector<double>& u_in,
		    vector<double>& u_out) {
  
  initialize_vector(n, u_out);
  double factor = 2.0/3.0;

  double tm;
  
  for (int i=0; i<N_target1; i++) {
    Surfc1[i].g[0] = u_in[3*i];
    Surfc1[i].g[1] = u_in[3*i+1];
    Surfc1[i].g[2] = u_in[3*i+2];
  }
  for (int i=0; i<N_target2; i++) {
    Surfc2[i].g[0] = u_in[3*N_target1 + 3*i];
    Surfc2[i].g[1] = u_in[3*N_target1 + 3*i+1];
    Surfc2[i].g[2] = u_in[3*N_target1 + 3*i+2];
  }
      
  interpolate_DL_density(1, h, N_quad1, Surfc1, N_target2, Target2);
  interpolate_DL_density(2, h, N_quad2, Surfc2, N_target1, Target1);


  // Double layer potential using the first surface

  vector<double> DL_onSurf1(3*N_target1,0), DL_comp1(3*N_target2,0);

  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    eval_Stokes_DL_onSurf_tree(1, h, DEL_on,
			       N_quad1, Surfc1,
			       N_target1, Target1,
			       DL_onSurf1);
  }
  else {
    eval_Stokes_DL_onSurf(DEL_on,
			  N_quad1, Surfc1,
			  N_target1, Target1,
			  DL_onSurf1);
    
  }
  tm = omp_get_wtime() - tm;
  cout << "Double layer evaluated Surfc1->Target1  with CPU time = "
       << tm << " seconds" << endl;
  
  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    
    eval_Stokes_DL_offSurf_tree(1, h, DEL_off,
				N_quad1, Surfc1,
				N_target2, Target2,
				DL_comp1);
    
  }
  else {
    eval_Stokes_DL_offSurf(h, DEL_off,
			   N_quad1, Surfc1,
			   N_target2, Target2,
			   DL_comp1);
    
  }
  
  tm = omp_get_wtime() - tm;
  cout << "Double layer evaluated Surfc1->Target2  with CPU time = "
       << tm << " seconds" << endl;
  
  for (int i=0; i<N_target1; i++) {    
    u_out[              3*i]   = u_in[              3*i]   - DL_onSurf1[3*i]   * factor;
    u_out[              3*i+1] = u_in[              3*i+1] - DL_onSurf1[3*i+1] * factor;
    u_out[              3*i+2] = u_in[              3*i+2] - DL_onSurf1[3*i+2] * factor;
  }
  for (int i=0; i<N_target2; i++) {    
    u_out[3*N_target1 + 3*i]   = u_in[3*N_target1 + 3*i]   - DL_comp1[3*i]   * factor;
    u_out[3*N_target1 + 3*i+1] = u_in[3*N_target1 + 3*i+1] - DL_comp1[3*i+1] * factor;
    u_out[3*N_target1 + 3*i+2] = u_in[3*N_target1 + 3*i+2] - DL_comp1[3*i+2] * factor;
  }
  
  
  // Double layer potential using the second surface
  
  vector<double> DL_onSurf2(3*N_target2,0), DL_comp2(3*N_target1,0);

  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    
    eval_Stokes_DL_offSurf_tree(2, h, DEL_off,
				N_quad2, Surfc2,
				N_target1, Target1,
				DL_comp2);  
    
  }
  else {
    eval_Stokes_DL_offSurf(h, DEL_off,
			   N_quad2, Surfc2,
			   N_target1, Target1,
			   DL_comp2);
    
  }
  
  tm = omp_get_wtime() - tm;
  cout << "Double layer evaluated Surfc2->Target1  with CPU time = "
       << tm << " seconds" << endl;
  
  tm = omp_get_wtime();
  if (use_tree == 1) { // use tree to speed up calculations
    eval_Stokes_DL_onSurf_tree(2, h, DEL_on,
			       N_quad2, Surfc2,
			       N_target2, Target2,
			       DL_onSurf2);
  }
  else {
    eval_Stokes_DL_onSurf(DEL_on,
			  N_quad2, Surfc2,
			  N_target2, Target2,
			  DL_onSurf2);
    
  }
  tm = omp_get_wtime() - tm;
  cout << "Double layer evaluated Surfc2->Target2  with CPU time = "
       << tm << " seconds" << endl;
  
  for (int i=0; i<N_target1; i++) {  
    u_out[              3*i]   -= DL_comp2[3*i]   * factor;
    u_out[              3*i+1] -= DL_comp2[3*i+1] * factor;
    u_out[              3*i+2] -= DL_comp2[3*i+2] * factor;
  }
  for (int i=0; i<N_target2; i++) {  
    u_out[3*N_target1 + 3*i]   -= DL_onSurf2[3*i]   * factor;
    u_out[3*N_target1 + 3*i+1] -= DL_onSurf2[3*i+1] * factor;
    u_out[3*N_target1 + 3*i+2] -= DL_onSurf2[3*i+2] * factor;
  }  
}

//*********************************************************************

// stencil for bivariate SIX point interpolation
void sixPoint_stencil(int j_start, int j_end,
		      int ind1, int ind2, int ind3,
		      double h, const vector<Surf_point>& Surfc,
		      double z0, double z_ctr, double tol,
		      double xLow, double yLow,
		      double xHigh, double yHigh,
		      vector< vector<int> >& stcl) {

    for (int j=j_start; j<j_end; j++) {
      
      double z = Surfc[j].x[ind3] - z_ctr;

      if ( ((z0 > 0) && (z > 0)) ||
	   ((z0 < 0) && (z < 0)) ) { // either "top" or "bottom"
	double X = Surfc[j].x[ind1];
	double Y = Surfc[j].x[ind2];
	if (abs(xLow-2.0*h - X)<tol) {         //first column
	  if (abs(yLow-2.0*h - Y)<tol)           stcl[0][0] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[0][1] = j;
	  else if (abs(yLow - Y)<tol)            stcl[0][2] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[0][3] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[0][4] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[0][5] = j;
	}
	else if (abs(xLow-h - X)<tol) {        //second column
	  if (abs(yLow-2.0*h - Y)<tol)           stcl[1][0] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[1][1] = j;
	  else if (abs(yLow - Y)<tol)            stcl[1][2] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[1][3] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[1][4] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[1][5] = j;
	}
	else if (abs(xLow - X)<1.e-14) {       //third column
	  if (abs(yLow-2.0*h - Y)<tol)           stcl[2][0] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[2][1] = j;
	  else if (abs(yLow - Y)<tol)            stcl[2][2] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[2][3] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[2][4] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[2][5] = j;
	}
	else if (abs(xHigh - X)<1.e-14) {      //fourth column
	  if (abs(yLow-2.0*h - Y)<tol)           stcl[3][0] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[3][1] = j;
	  else if (abs(yLow - Y)<tol)            stcl[3][2] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[3][3] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[3][4] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[3][5] = j;
	}
	else if (abs(xHigh+h - X)<1.e-14) {    //fifth column
	  if (abs(yLow-2.0*h - Y)<tol)           stcl[4][0] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[4][1] = j;
	  else if (abs(yLow - Y)<tol)            stcl[4][2] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[4][3] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[4][4] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[4][5] = j;
	}
	else if (abs(xHigh+2.0*h - X)<1.e-14) {//sixth column
	  if (abs(yLow-2.0*h - Y)<tol)           stcl[5][0] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[5][1] = j;
	  else if (abs(yLow - Y)<tol)            stcl[5][2] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[5][3] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[5][4] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[5][5] = j;
	}
      }
    }
}
    
void sixPoint_coeff(double p, vector<double>& coeff) {
  coeff[0] = -p*(p*p-1.0)*(p-2.0)*(p-3.0)/120.0;
  coeff[1] = p*(p-1.0)*(p*p-4.0)*(p-3.0)/24.0;
  coeff[2] = -(p*p-1.0)*(p*p-4.0)*(p-3.0)/12.0;
  coeff[3] = p*(p+1.0)*(p*p-4.0)*(p-3.0)/12.0;
  coeff[4] = -p*(p*p-1.0)*(p+2.0)*(p-3.0)/24.0;
  coeff[5] = p*(p*p-1.0)*(p*p-4.0)/120.0;
}

double sixPointInterpolate(vector< vector<double> >& f, double p, double q) {
  vector<double> a(6,0), b(6,0), f_temp(6,0);

  sixPoint_coeff(p, a);
  sixPoint_coeff(q, b);
  
  for (int j=0; j<6; j++) {
    f_temp[j] = 0.0;
    for (int i=0; i<6; i++) {
      f_temp[j] += a[i] * f[i][j];
    }
  }	
  
  double f_sixPoint = 0.0;
  for (int j=0; j<6; j++) {
    f_sixPoint += b[j] * f_temp[j];
  }
  
  return f_sixPoint;
}
//*********************************************************************

// stencil for bivariate EIGHT point interpolation
void eightPoint_stencil(int j_start, int j_end,
			int ind1, int ind2, int ind3,
			double h, const vector<Surf_point>& Surfc,
			double z0, double z_ctr, double tol,
			double xLow, double yLow,
			double xHigh, double yHigh,
			vector< vector<int> >& stcl) {

    for (int j=j_start; j<j_end; j++) {
      
      double z = Surfc[j].x[ind3] - z_ctr;

      if ( ((z0 > 0) && (z > 0)) ||
	   ((z0 < 0) && (z < 0)) ) { // either "top" or "bottom"
	double X = Surfc[j].x[ind1];
	double Y = Surfc[j].x[ind2];
	if (abs(xLow-3.0*h - X)<tol) {         //first column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[0][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[0][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[0][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[0][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[0][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[0][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[0][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[0][7] = j;
	}
	else if (abs(xLow-2.0*h - X)<tol) {    //second column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[1][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[1][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[1][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[1][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[1][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[1][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[1][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[1][7] = j;
	}
	else if (abs(xLow-h - X)<tol) {        //third column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[2][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[2][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[2][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[2][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[2][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[2][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[2][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[2][7] = j;
	}
	else if (abs(xLow - X)<1.e-14) {       //fourth column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[3][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[3][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[3][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[3][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[3][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[3][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[3][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[3][7] = j;
	}
	else if (abs(xHigh - X)<1.e-14) {      //fifth column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[4][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[4][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[4][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[4][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[4][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[4][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[4][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[4][7] = j;
	}
	else if (abs(xHigh+h - X)<1.e-14) {    //sixth column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[5][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[5][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[5][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[5][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[5][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[5][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[5][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[5][7] = j;
	}
	else if (abs(xHigh+2.0*h - X)<1.e-14) {//seventh column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[6][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[6][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[6][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[6][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[6][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[6][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[6][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[6][7] = j;
	}
	else if (abs(xHigh+3.0*h - X)<1.e-14) {//eighth column
	  if (abs(yLow-3.0*h - Y)<tol)           stcl[7][0] = j;
	  else if (abs(yLow-2.0*h - Y)<tol)      stcl[7][1] = j;
	  else if (abs(yLow-h - Y)<tol)          stcl[7][2] = j;
	  else if (abs(yLow - Y)<tol)            stcl[7][3] = j;
	  else if (abs(yHigh - Y)<tol)           stcl[7][4] = j;
	  else if (abs(yHigh+h - Y)<tol)         stcl[7][5] = j;
	  else if (abs(yHigh+2.0*h - Y)<tol)     stcl[7][6] = j;
	  else if (abs(yHigh+3.0*h - Y)<tol)     stcl[7][7] = j;
	}
      }
    }
}
    
void eightPoint_coeff(double p, vector<double>& coeff) {
  coeff[0] = -p*(p*p-1.0)*(p-3.0)*(p-4.0)*(p*p-4.0)/5040.0;
  coeff[1] = p*(p*p-1.0)*(p-2.0)*(p-4.0)*(p*p-9.0)/720.0;
  coeff[2] = -p*(p-1.0)*(p-4.0)*(p*p-4.0)*(p*p-9.0)/240.0;
  coeff[3] = (p*p-1.0)*(p-4.0)*(p*p-4.0)*(p*p-9.0)/144.0;
  coeff[4] = -p*(p+1.0)*(p-4.0)*(p*p-4.0)*(p*p-9.0)/144.0;
  coeff[5] = p*(p*p-1.0)*(p+2.0)*(p-4.0)*(p*p-9.0)/240.0;
  coeff[6] = -p*(p*p-1.0)*(p+3.0)*(p-4.0)*(p*p-4.0)/720.0;
  coeff[7] = p*(p*p-1.0)*(p*p-4.0)*(p*p-9.0)/5040.0;
}

double eightPointInterpolate(vector< vector<double> >& f, double p, double q) {
  int n=8;
  vector<double> a(n,0), b(n,0), f_temp(n,0);

  eightPoint_coeff(p, a);
  eightPoint_coeff(q, b);
  
  for (int j=0; j<n; j++) {
    f_temp[j] = 0.0;
    for (int i=0; i<n; i++) {
      f_temp[j] += a[i] * f[i][j];
    }
  }	
  
  double f_eightPoint = 0.0;
  for (int j=0; j<n; j++) {
    f_eightPoint += b[j] * f_temp[j];
  }
  
  return f_eightPoint;
}
//*********************************************************************

void interpolate_DL_density(int surf_num, double h, int N_quad,
			    const vector<Surf_point>& Surfc,
			    int N_target,
			    vector<Target_point>& Target) {
  
  vector<double> pt(3,0), x0(3,0), g0(3,0);

  double err = 0.0;
  double exact = 0.0;

  int num_interpolated = 0;
  
  for (int i=0; i<N_target; i++) {
    
    Target[i].nrst_g[0] = 0;
    Target[i].nrst_g[1] = 0;
    Target[i].nrst_g[2] = 0;

    if ( Target[i].nrst == 0 )  continue;

    num_interpolated++;
    
    int i0 = Target[i].nrst;

    vector<double> x0(3,0);
    x0[0] = Target[i].nrst_x;
    x0[1] = Target[i].nrst_y;
    x0[2] = Target[i].nrst_z;
    vector<double> norm0 = Target[i].nrst_Nrml;
    int ind1 = Target[i].i1;
    int ind2 = Target[i].i2;
    int ind3 = Target[i].i3;

    // determine which coordinate plane to interpolate in
    int j_start, j_end;

    
    if ( ( abs(norm0[2]) >= abs(norm0[0])) && (abs(norm0[2]) >= abs(norm0[1])) ) {
      //ind1 = 0; ind2 = 1; ind3 = 2; // interpolation variables are (x,y)
      if ( (ind1 != 0) || (ind2 != 1) || (ind3 != 2) ) cout << "wrong coordinate plane!!!" << endl; 
      j_start = 2*N_quad/3; j_end = N_quad; // will search in R3 portion of quadrature points
    }
    else if ( ( abs(norm0[1]) >= abs(norm0[0])) && (abs(norm0[1]) >= abs(norm0[2])) ) {
      //ind1 = 0; ind2 = 2; ind3 = 1; // interpolation variables are (x,z)
      if ( (ind1 != 0) || (ind2 != 2) || (ind3 != 1) ) cout << "wrong coordinate plane!!!" << endl; 
      j_start = N_quad/3; j_end = 2*N_quad/3; // will search in R2 portion of quadrature points
    }
    else {
      //ind1 = 1; ind2 = 2; ind3 = 0; // interpolation variables are (y,z)
      if ( (ind1 != 1) || (ind2 != 2) || (ind3 != 0) ) cout << "wrong coordinate plane!!!" << endl; 
      j_start = 0; j_end = N_quad/3; // will search in R1 portion of quadrature points
    }
    
    j_start = 0;
    j_end = N_quad;
    
    vector<double> ctr(3,0);
    if (surf_num == 1) {
      ctr[0] = ctr1_x;
      ctr[1] = ctr1_y;
      ctr[2] = ctr1_z;
    }
    if (surf_num == 2) {
      ctr[0] = ctr2_x;
      ctr[1] = ctr2_y;
      ctr[2] = ctr2_z;
    }
    
    // find and record indices of the 4 nearby points
    double xLow = floor((x0[ind1]-ctr[ind1])/h) * h + ctr[ind1];
    double yLow = floor((x0[ind2]-ctr[ind2])/h) * h + ctr[ind2];
    double xHigh = xLow + h;
    double yHigh = yLow + h;

    // generate interpolation points and np x np matrix
    int np = 6;
    vector<int> line_int(np,0);
    vector< vector<int> > stcl(np, line_int);

    double z0 = x0[ind3]-ctr[ind3];

    double tol = 1.e-14;

    if (np == 6) {
      sixPoint_stencil(j_start, j_end,
		       ind1, ind2, ind3,
		       h, Surfc,
		       z0, ctr[ind3], tol,
		       xLow, yLow,
		       xHigh, yHigh,
		       stcl);
    }
    else {
      eightPoint_stencil(j_start, j_end,
			 ind1, ind2, ind3,
			 h, Surfc,
			 z0, ctr[ind3], tol,
			 xLow, yLow,
			 xHigh, yHigh,
			 stcl);
    }
    

    vector<double> line_double(np,0);
    vector< vector<double> > g_x(np, line_double), g_y(np, line_double), g_z(np, line_double);
    
    for (int k1=0; k1<np; k1++) {
      for (int k2=0; k2<np; k2++) {
	if (stcl[k1][k2] == 0) cout << "stencil not set" << endl;
	g_x[k1][k2] = Surfc[stcl[k1][k2]].g[0];
	g_y[k1][k2] = Surfc[stcl[k1][k2]].g[1];
	g_z[k1][k2] = Surfc[stcl[k1][k2]].g[2];
      }
    }


    double p = (x0[ind1] - xLow)/h;
    double q = (x0[ind2] - yLow)/h;

    if (np == 6) {
      Target[i].nrst_g[0] = sixPointInterpolate(g_x, p, q);
      Target[i].nrst_g[1] = sixPointInterpolate(g_y, p, q);
      Target[i].nrst_g[2] = sixPointInterpolate(g_z, p, q);
    }
    else {
      Target[i].nrst_g[0] = eightPointInterpolate(g_x, p, q);
      Target[i].nrst_g[1] = eightPointInterpolate(g_y, p, q);
      Target[i].nrst_g[2] = eightPointInterpolate(g_z, p, q);
    }
  }
}

//*********************************************************************

void set_surface_force(int surf_num, int N_quad,
		       vector<Surf_point>& Surfc) {
  
  int rotated = 0;
  if (surf_num == 2) rotated = 1;
  
  for (int i=0; i<N_quad; i++) {
    Surfc[i].f = surfc_force(surf_num, rotated, Surfc[i].x, Surfc[i].Nrml);
  }   
}
