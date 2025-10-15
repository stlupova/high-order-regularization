#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <fftw3.h>
#include <omp.h>

#include "stokes.h"

using namespace std;

void test_Stokes_Poisson3d_FD(double L, int N, double h, double DEL);
void solve_Poisson(int N, double L, double h,
		   int*** nrst,
		   int N_target, const vector<Target_point>& Target,
		   double*** p,
		   double*** u,
		   int variable);

void eval_p_BI_cubeEdges(int N, double h,
		       int N_quad, const vector<Surf_point>& Surfc,
		       const vector<double>& f_ex,
		       double*** p_bi);
void eval_p_BI_nearSurfc(int N, double h, double DEL,
		       int N_quad, const vector<Surf_point>& Surfc,
		       int*** nrst,
		       int N_target, const vector<Target_point>& Target,
		       const vector<double>& f_ex,
		       double*** p_bi);

void eval_u_BI_cubeEdges(int N, double h,
		       int N_quad, const vector<Surf_point>& Surfc,
		       const vector<double>& f_ex,
		       double*** ux_bi, double*** uy_bi, double*** uz_bi);
void eval_u_BI_nearSurfc(int N, double h, double DEL,
		       int N_quad, const vector<Surf_point>& Surfc,
		       int*** nrst,
		       int N_target, const vector<Target_point>& Target,
		       const vector<double>& f_ex,
		       double*** ux_bi, double*** uy_bi, double*** uz_bi);

void fftw_DST_Poisson3d_FD(int N, double* f);

void extend_into_cube(int N, double L, double h,
		      double*** g, double*** G);

void lapl_4ord_inCube(int N, double h,
		      double*** G,
		      double* rhs);
void lapl_4ord_nearSurfc(int N, double h,
			 int*** nrst, const vector<Target_point>& Target,
			 double*** u, double* rhs);
void dp_dx_4ord(int N, double h,
		int*** nrst, const vector<Target_point>& Target,
		double*** p, double* rhs);
void dp_dy_4ord(int N, double h,
		int*** nrst, const vector<Target_point>& Target,
		double*** p, double* rhs);
void dp_dz_4ord(int N, double h,
		int*** nrst, const vector<Target_point>& Target,
		double*** p, double* rhs);

//*********************************************************************

int main(int argc, char** argv) {

  srand(time(NULL));
  
  double tm = omp_get_wtime(); 

  int N = 3 * 32;
  double L = 3.0;
  double h = L / N; // grid spacing in box [0,L]^3

  int order = 7;
  double kappa0 = 4.0;
  double q = 5.0/7.0;
  double kappa = kappa0 * pow(1.0/64.0, 1.0-q);
  double DEL = kappa * pow(h, q);

  cout << "h = 1/" << N/3 << "   order = " << order << "   kappa0 = " << kappa0
       << "   DEL = " << kappa << " * h^(" << q << ")" << endl;

  
  if ( surfc_is_sphere ) {
    cout << "Surface is a sphere;  ";
  }
  else if ( surfc_elliptical ) {
    cout << "Surface is an ellipsoid (" << ellipse_a << ", "
	 << ellipse_b << ", " << ellipse_c << ");  ";
  }
  else {
    cout << "Surface is the four-atom molecule;  ";
  }
  
  test_Stokes_Poisson3d_FD(L, N, h, DEL);

  tm = omp_get_wtime() - tm;
  cout << "Total CPU time = " << tm << " seconds" << endl;
  
  return 0;
}

//*********************************************************************

void test_Stokes_Poisson3d_FD(double L, int N, double h, double DEL) {

  int Nm1_cubed = pow(N-1,3);
  int Np1_cubed = pow(N+1,3);

  double BAND1 = 4.0; //must be greater than 2*sqrt(3) for stencil to cross b'ry
  double BAND2 = 2.0; //must be greater than sqrt(3) for discrete laplacian
  
  ofstream p_file("p_midY_h32");
  ofstream err_u_file("err_u_midY_h32");
  ofstream err_p_file("err_p_midY_h32");
  ofstream err_Du_file("err_gradu_midY_h32");
  
  // allocate memory
  int*** nrst  = new int**[N+1];
  double*** p_bi = new double**[N+1];
  for (int i=0; i<=N; i++) {
    nrst[i]  = new int*[N+1];
    p_bi[i] = new double*[N+1];
    for (int j=0; j<=N; j++) {
      nrst[i][j]  = new int[N+1];
      p_bi[i][j] = new double[N+1];
    }
  }

  clock_t tm = clock();

  vector<Surf_point> Surfc;
  Generate_Surface(N, h, &Surfc);
  int N_quad = Surfc.size();
  cout << "Number of quadrature points = " << N_quad << endl;

  tm = clock() - tm;
  cout << "   CPU time to create surface = "
       << ((float)tm)/CLOCKS_PER_SEC << " seconds" << endl;
   
  // set density function
  vector<double> f_ex(3*N_quad,0), fe(3,0);

  for (int i=0; i<N_quad; i++) {
    fe = stokeslet_density(Surfc[i].x);
    f_ex[3*i]   = fe[0];
    f_ex[3*i+1] = fe[1];
    f_ex[3*i+2] = fe[2];    
  }
  

  tm = clock();
  vector<Target_point> Target; // array of target points
  find_nrst_surfc(N, h,
		  BAND1, BAND2,
		  N_quad, Surfc,
		  nrst, &Target);
  int N_target = Target.size();
  cout << "Number of target points within " << BAND1 << "h of surface = " << N_target << endl;
  
  tm = clock() - tm;
  cout << "    CPU time to find nearest data = "
       << ((float)tm)/CLOCKS_PER_SEC << " seconds" << endl;
  

  // STAGE 1 : STOKES PRESSURE
  tm = omp_get_wtime();

  eval_p_BI_cubeEdges(N, h,
		      N_quad, Surfc,
		      f_ex,
		      p_bi);
  
  tm = omp_get_wtime() - tm;
  cout << "CPU time to compute P at edges of box = "
       << tm << " seconds" << endl;
  

  tm = omp_get_wtime();

  eval_p_BI_nearSurfc(N, h, DEL,
		      N_quad, Surfc,
		      nrst,
		      N_target, Target,
		      f_ex,
		      p_bi);

  tm = omp_get_wtime() - tm;
  cout << "CPU time to compute P in BAND1 = "
       << tm << " seconds" << endl << endl;

  
  solve_Poisson(N, L, h,
		nrst,
		N_target, Target,
		p_bi, // this won't be used
		p_bi,
		0);  // for p, use anything except 1,2,3
  
  double*** p_ex = new double**[N+1];
  for (int i=0; i<=N; i++) {
    p_ex[i] = new double*[N+1];
    for (int j=0; j<=N; j++) {
      p_ex[i][j] = new double[N+1];
    }
  }

  vector<double> x_vec(3,0);

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {
	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;

	double levl = phi(x_vec);
	if ( ( levl > 0 ) || ( abs(levl) < 1.e-12 ) ) { // outside or on surface
	  p_ex[i][j][k] = p_exact(i*h, j*h, k*h);
	}
	else {
	  p_ex[i][j][k] = 0.0;
	}
      }
    }
  }

  
  double err_max = 0.0;
  double err_l2 = 0.0;
  double val_max = 0.0;
  double val_l2 = 0.0;

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	// compute error max and l2 norms
	double err_temp = abs(p_bi[i][j][k] - p_ex[i][j][k]);
	if ( err_temp > err_max ) {
	  err_max = err_temp;
	}
	if ( abs(p_ex[i][j][k]) > val_max ) {
	  val_max = abs(p_ex[i][j][k]);
	}
	err_l2 += err_temp * err_temp;
	val_l2 += p_ex[i][j][k] * p_ex[i][j][k];

      } 
    }
  }
  err_l2 = sqrt(err_l2 / Np1_cubed);
  val_l2 = sqrt(val_l2 / Np1_cubed);
  cout << "entire grid  h=1/" << N/3 << "   max err (p) = " << err_max
       << "    l2 err (p) = " << err_l2
       << endl << "    max norm (p) = " << val_max
       << "    l2 norm (p) = " << val_l2 << endl << endl;

    
  err_max = 0.0;
  err_l2 = 0.0;
  val_max = 0.0;
  val_l2 = 0.0;
  int count = 0;
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	
	double levl = phi(x_vec);
	double err_temp;
	
	if ( ( levl > 0 )  || ( abs(levl) < 1.e-12 ) ) { // outside or on surface

	  // compute error max and l2 norms
	  err_temp = abs(p_bi[i][j][k] - p_ex[i][j][k]);
	  if ( err_temp > err_max ) {
	    err_max = err_temp;
	  }
	  if ( abs(p_ex[i][j][k]) > val_max ) {
	    val_max = abs(p_ex[i][j][k]);
	  }
	  err_l2 += err_temp * err_temp;
	  val_l2 += p_ex[i][j][k] * p_ex[i][j][k];

	  count++;
	}
	else {
	  err_temp = 0.0;
	}
	if ( abs(j*h - 1.5) < 1.e-8 ) {
	  p_file << setprecision(16) << p_ex[i][j][k] << "  ";
	  err_p_file << setprecision(16) << err_temp << "  ";
	}
      } 
    p_file << endl;
    err_p_file << endl;
    }
  }
  p_file.close();
  err_p_file.close();

  err_l2 = sqrt(err_l2 / count);
  val_l2 = sqrt(val_l2 / count);
  cout << "outside + on  h=1/" << N/3 << "   max err (p) = " << err_max
       << "    l2 err (p) = " << err_l2
       << endl << "    max norm (p) = " << val_max
       << "    l2 norm (p) = " << val_l2 << endl << endl;


  err_max = 0.0;
  err_l2 = 0.0;
  val_max = 0.0;
  val_l2 = 0.0;
  count = 0;
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	
	double levl = phi(x_vec);
	if ( ( levl > 0 ) ) { // outside surface only

	  // compute error max and l2 norms
	  double err_temp = abs(p_bi[i][j][k] - p_ex[i][j][k]);
	  if ( err_temp > err_max ) {
	    err_max = err_temp;
	  }
	  if ( abs(p_ex[i][j][k]) > val_max ) {
	    val_max = abs(p_ex[i][j][k]);
	  }
	  err_l2 += err_temp * err_temp;
	  val_l2 += p_ex[i][j][k] * p_ex[i][j][k];

	  count++;
	}
      } 
    }
  }
  err_l2 = sqrt(err_l2 / count);
  val_l2 = sqrt(val_l2 / count);
  cout << "outside only  h=1/" << N/3 << "   max err (p) = " << err_max
       << "    l2 err (p) = " << err_l2
       << endl << "    max norm (p) = " << val_max
       << "    l2 norm (p) = " << val_l2 << endl << endl;

  
  
  // STAGE 2 : STOKES VELOCITY

  // allocate memory
  double*** ux = new double**[N+1];
  double*** uy = new double**[N+1];
  double*** uz = new double**[N+1];

  for (int i=0; i<=N; i++) {
    ux[i] = new double*[N+1];
    uy[i] = new double*[N+1];
    uz[i] = new double*[N+1];
    for (int j=0; j<=N; j++) {
      ux[i][j] = new double[N+1];
      uy[i][j] = new double[N+1];
      uz[i][j] = new double[N+1];
    }
  }

  tm = omp_get_wtime();
  
  eval_u_BI_cubeEdges(N, h,
		      N_quad, Surfc,
		      f_ex,
		      ux, uy, uz);
  
  tm = omp_get_wtime() - tm;
  cout << "CPU time to compute U at edges of box = "
       << tm << " seconds" << endl;
  

  tm = omp_get_wtime();
  
  eval_u_BI_nearSurfc(N, h, DEL,
		      N_quad, Surfc,
		      nrst,
		      N_target, Target,
		      f_ex,
		      ux, uy, uz);
  
  tm = omp_get_wtime() - tm;
  cout << "CPU time to compute U in BAND1 = "
       << tm << " seconds" << endl << endl;


  // FIRST COMPONENT OF VELOCITY
  solve_Poisson(N, L, h,
		nrst,
		N_target, Target,
		p_bi,
		ux,
		1); 
  // SECOND COMPONENT OF VELOCITY
  solve_Poisson(N, L, h,
		nrst,
		N_target, Target,
		p_bi,
		uy,
		2); 

  // THIRD COMPONENT OF VELOCITY
  solve_Poisson(N, L, h,
		nrst,
		N_target, Target,
		p_bi,
		uz,
		3); 

  double*** ux_ex = new double**[N+1];
  double*** uy_ex = new double**[N+1];
  double*** uz_ex = new double**[N+1];
  for (int i=0; i<=N; i++) {
    ux_ex[i] = new double*[N+1];
    uy_ex[i] = new double*[N+1];
    uz_ex[i] = new double*[N+1];
    for (int j=0; j<=N; j++) {
      ux_ex[i][j] = new double[N+1];
      uy_ex[i][j] = new double[N+1];
      uz_ex[i][j] = new double[N+1];
    }
  }

  vector<double> u_ex(3,0);

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {
	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;

	double levl = phi(x_vec);
	if ( ( levl > 0 ) || ( abs(levl) < 1.e-12 ) ) { // outside or on surface
	  u_ex = SL_exact(i*h, j*h, k*h);
	  ux_ex[i][j][k] = u_ex[0];
	  uy_ex[i][j][k] = u_ex[1];
	  uz_ex[i][j][k] = u_ex[2];
	}
	else {
	  ux_ex[i][j][k] = 1.0;
	  uy_ex[i][j][k] = 0.0;
	  uz_ex[i][j][k] = 0.0;
	}
      }
    }
  }

  err_max = 0.0;
  err_l2 = 0.0;
  val_max = 0.0;
  val_l2 = 0.0;

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	// compute error max and l2 norms
	double err_temp = (ux[i][j][k]-ux_ex[i][j][k]) * (ux[i][j][k]-ux_ex[i][j][k])
	                + (uy[i][j][k]-uy_ex[i][j][k]) * (uy[i][j][k]-uy_ex[i][j][k])
	                + (uz[i][j][k]-uz_ex[i][j][k]) * (uz[i][j][k]-uz_ex[i][j][k]);
	double val_temp = ux[i][j][k] * ux[i][j][k]
	                + uy[i][j][k] * uy[i][j][k]
	                + uz[i][j][k] * uz[i][j][k];
	if ( sqrt(err_temp) > err_max ) {
	  err_max = sqrt(err_temp);
	}
	if ( sqrt(val_temp) > val_max ) {
	  val_max = sqrt(val_temp);
	}
	err_l2 += err_temp;
	val_l2 += val_temp;
      }
    }
  }
  err_l2 = sqrt(err_l2 / Np1_cubed);
  val_l2 = sqrt(val_l2 / Np1_cubed);
  cout << "entire grid  h=1/" << N/3 << "   max err (u) = " << err_max
       << "    l2 err (u) = " << err_l2
       << endl << "    with max norm (u) = " << val_max
       << "    l2 norm (u) = " << val_l2 << endl << endl;


  err_max = 0.0;
  err_l2 = 0.0;
  val_max = 0.0;
  val_l2 = 0.0;
  count = 0;

  int max_i, max_j, max_k;
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	
	double levl = phi(x_vec);
	double err_temp, val_temp;
	if ( ( levl > 0 )  || ( abs(levl) < 1.e-12 ) ) { // outside or on surface

	  // compute error max and l2 norms
	  err_temp = (ux[i][j][k]-ux_ex[i][j][k]) * (ux[i][j][k]-ux_ex[i][j][k])
	           + (uy[i][j][k]-uy_ex[i][j][k]) * (uy[i][j][k]-uy_ex[i][j][k])
	           + (uz[i][j][k]-uz_ex[i][j][k]) * (uz[i][j][k]-uz_ex[i][j][k]);
	  val_temp = ux[i][j][k] * ux[i][j][k]
	           + uy[i][j][k] * uy[i][j][k]
	           + uz[i][j][k] * uz[i][j][k];
	  if ( sqrt(err_temp) > err_max ) {
	    err_max = sqrt(err_temp);
	    max_i = i;
	    max_j = j;
	    max_k = k;
	  }
	  if ( sqrt(val_temp) > val_max ) {
	    val_max = sqrt(val_temp);
	  }
	  err_l2 += err_temp;
	  val_l2 += val_temp;

	  count++;
	}
	else {
	  err_temp = 0.0;
	}

	if ( abs(j*h - 1.5) < 1.e-8 ) {
	  err_u_file << setprecision(16) << sqrt(err_temp) << "  ";
	}
      }
    }
    err_u_file << endl;
  }
  err_u_file.close();

  err_l2 = sqrt(err_l2 / count);
  val_l2 = sqrt(val_l2 / count);
  cout << "outside + on  h=1/" << N/3 << "   max err (u) = " << err_max
       << "    l2 err (u) = " << err_l2
       << endl << "    with max norm (u) = " << val_max
       << "    l2 norm (u) = " << val_l2 << endl << endl;
  cout << "max i = " << max_i << "  max j = " << max_j << "  max k = " << max_k << endl;

  err_max = 0.0;
  err_l2 = 0.0;
  val_max = 0.0;
  val_l2 = 0.0;
  count = 0;
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	
	double levl = phi(x_vec);
	if ( ( levl > 0 ) ) { // outside surface only

	  // compute error max and l2 norms
	  double err_temp = (ux[i][j][k]-ux_ex[i][j][k]) * (ux[i][j][k]-ux_ex[i][j][k])
	                  + (uy[i][j][k]-uy_ex[i][j][k]) * (uy[i][j][k]-uy_ex[i][j][k])
	                  + (uz[i][j][k]-uz_ex[i][j][k]) * (uz[i][j][k]-uz_ex[i][j][k]);
	  double val_temp = ux[i][j][k] * ux[i][j][k]
	                  + uy[i][j][k] * uy[i][j][k]
	                  + uz[i][j][k] * uz[i][j][k];
	  if ( sqrt(err_temp) > err_max ) {
	    err_max = sqrt(err_temp);
	  }
	  if ( sqrt(val_temp) > val_max ) {
	    val_max = sqrt(val_temp);
	  }
	  err_l2 += err_temp;
	  val_l2 += val_temp;

	  count++;
	}
      }
    }
  }
  err_l2 = sqrt(err_l2 / count);
  val_l2 = sqrt(val_l2 / count);
  cout << "outside only  h=1/" << N/3 << "   max err (u) = " << err_max
       << "    l2 err (u) = " << err_l2
       << endl << "    with max norm (u) = " << val_max
       << "    l2 norm (u) = " << val_l2 << endl << endl;
  

  // deallocate memory
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      delete[] nrst[i][j];
      delete[] p_bi[i][j];
      delete[] p_ex[i][j];
    }
    delete[] nrst[i];
    delete[] p_bi[i];
    delete[] p_ex[i];
  }
  delete[] nrst;
  delete[] p_bi;
  delete[] p_ex;


  
  // Compute error in gradient of u to O(h^4) on entire grid

  err_max = 0.0;
  err_l2 = 0.0;
  count = 0;

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	// Derivatives in x
	double du_dx = 0.0;
	double du_dx_ex = 0.0;
	
	if ( i==0 || i==1 ) {
	  du_dx = - (25.0/12.0) * ux[i][j][k] + 4.0 * ux[i+1][j][k] - 3.0 * ux[i+2][j][k]
	    + (4.0/3.0) * ux[i+3][j][k] - 0.25 * ux[i+4][j][k];
	  du_dx /= h;
	  
	  du_dx_ex = - (25.0/12.0) * ux_ex[i][j][k] + 4.0 * ux_ex[i+1][j][k] - 3.0 * ux_ex[i+2][j][k]
	    + (4.0/3.0) * ux_ex[i+3][j][k] - 0.25 * ux_ex[i+4][j][k];
	  du_dx_ex /= h;
	}
	else if ( i==N || i==(N-1) ) {
	  du_dx = (25.0/12.0) * ux[i][j][k] - 4.0 * ux[i-1][j][k] + 3.0 * ux[i-2][j][k]
	    - (4.0/3.0) * ux[i-3][j][k] + 0.25 * ux[i-4][j][k];
	  du_dx /= h;
	  
	  du_dx_ex = (25.0/12.0) * ux_ex[i][j][k] - 4.0 * ux_ex[i-1][j][k] + 3.0 * ux_ex[i-2][j][k]
	    - (4.0/3.0) * ux_ex[i-3][j][k] + 0.25 * ux_ex[i-4][j][k];
	  du_dx_ex /= h;
	}
	else {
	  du_dx = - ux[i+2][j][k] + 8.0 * ux[i+1][j][k] - 8.0 * ux[i-1][j][k] + ux[i-2][j][k];
	  du_dx /= 12.0 * h;
	  
	  du_dx_ex = - ux_ex[i+2][j][k] + 8.0 * ux_ex[i+1][j][k] - 8.0 * ux_ex[i-1][j][k] + ux_ex[i-2][j][k];
	  du_dx_ex /= 12.0 * h;
	}
	
	double dv_dx = 0.0;
	double dv_dx_ex = 0.0;
	
	if ( i==0 || i==1 ) {
	  dv_dx = - (25.0/12.0) * uy[i][j][k] + 4.0 * uy[i+1][j][k] - 3.0 * uy[i+2][j][k]
	    + (4.0/3.0) * uy[i+3][j][k] - 0.25 * uy[i+4][j][k];
	  dv_dx /= h;
	  
	  dv_dx_ex = - (25.0/12.0) * uy_ex[i][j][k] + 4.0 * uy_ex[i+1][j][k] - 3.0 * uy_ex[i+2][j][k]
	    + (4.0/3.0) * uy_ex[i+3][j][k] - 0.25 * uy_ex[i+4][j][k];
	  dv_dx_ex /= h;
	}
	else if ( i==N || i==(N-1) ) {
	  dv_dx = (25.0/12.0) * uy[i][j][k] - 4.0 * uy[i-1][j][k] + 3.0 * uy[i-2][j][k]
	    - (4.0/3.0) * uy[i-3][j][k] + 0.25 * uy[i-4][j][k];
	  dv_dx /= h;
	  
	  dv_dx_ex = (25.0/12.0) * uy_ex[i][j][k] - 4.0 * uy_ex[i-1][j][k] + 3.0 * uy_ex[i-2][j][k]
	    - (4.0/3.0) * uy_ex[i-3][j][k] + 0.25 * uy_ex[i-4][j][k];
	  dv_dx_ex /= h;
	}
	else {
	  dv_dx = - uy[i+2][j][k] + 8.0 * uy[i+1][j][k] - 8.0 * uy[i-1][j][k] + uy[i-2][j][k];
	  dv_dx /= 12.0 * h;
	  
	  dv_dx_ex = - uy_ex[i+2][j][k] + 8.0 * uy_ex[i+1][j][k] - 8.0 * uy_ex[i-1][j][k] + uy_ex[i-2][j][k];
	  dv_dx_ex /= 12.0 * h;
	}
	
	double dw_dx = 0.0;
	double dw_dx_ex = 0.0;
	
	if ( i==0 || i==1 ) {
	  dw_dx = - (25.0/12.0) * uz[i][j][k] + 4.0 * uz[i+1][j][k] - 3.0 * uz[i+2][j][k]
	    + (4.0/3.0) * uz[i+3][j][k] - 0.25 * uz[i+4][j][k];
	  dw_dx /= h;
	  
	  dw_dx_ex = - (25.0/12.0) * uz_ex[i][j][k] + 4.0 * uz_ex[i+1][j][k] - 3.0 * uz_ex[i+2][j][k]
	    + (4.0/3.0) * uz_ex[i+3][j][k] - 0.25 * uz_ex[i+4][j][k];
	  dw_dx_ex /= h;
	}
	else if ( i==N || i==(N-1) ) {
	  dw_dx = (25.0/12.0) * uz[i][j][k] - 4.0 * uz[i-1][j][k] + 3.0 * uz[i-2][j][k]
	    - (4.0/3.0) * uz[i-3][j][k] + 0.25 * uz[i-4][j][k];
	  dw_dx /= h;
	  
	  dw_dx_ex = (25.0/12.0) * uz_ex[i][j][k] - 4.0 * uz_ex[i-1][j][k] + 3.0 * uz_ex[i-2][j][k]
	    - (4.0/3.0) * uz_ex[i-3][j][k] + 0.25 * uz_ex[i-4][j][k];
	  dw_dx_ex /= h;
	}
	else {
	  dw_dx = - uz[i+2][j][k] + 8.0 * uz[i+1][j][k] - 8.0 * uz[i-1][j][k] + uz[i-2][j][k];
	  dw_dx /= 12.0 * h;
	  
	  dw_dx_ex = - uz_ex[i+2][j][k] + 8.0 * uz_ex[i+1][j][k] - 8.0 * uz_ex[i-1][j][k] + uz_ex[i-2][j][k];
	  dw_dx_ex /= 12.0 * h;
	}
	

	// Derivatives in y
	double du_dy = 0.0;
	double du_dy_ex = 0.0;
	  
	if ( j==0 || j==1 ) {
	  du_dy = - (25.0/12.0) * ux[i][j][k] + 4.0 * ux[i][j+1][k] - 3.0 * ux[i][j+2][k]
	    + (4.0/3.0) * ux[i][j+3][k] - 0.25 * ux[i][j+4][k];
	  du_dy /= h;
	  
	  du_dy_ex = - (25.0/12.0) * ux_ex[i][j][k] + 4.0 * ux_ex[i][j+1][k] - 3.0 * ux_ex[i][j+2][k]
	    + (4.0/3.0) * ux_ex[i][j+3][k] - 0.25 * ux_ex[i][j+4][k];
	  du_dy_ex /= h;
	}
	else if ( j==N || j==(N-1) ) {
	  du_dy = (25.0/12.0) * ux[i][j][k] - 4.0 * ux[i][j-1][k] + 3.0 * ux[i][j-2][k]
	    - (4.0/3.0) * ux[i][j-3][k] + 0.25 * ux[i][j-4][k];
	  du_dy /= h;
	  
	  du_dy_ex = (25.0/12.0) * ux_ex[i][j][k] - 4.0 * ux_ex[i][j-1][k] + 3.0 * ux_ex[i][j-2][k]
	    - (4.0/3.0) * ux_ex[i][j-3][k] + 0.25 * ux_ex[i][j-4][k];
	  du_dy_ex /= h;
	}
	else {
	  du_dy = - ux[i][j+2][k] + 8.0 * ux[i][j+1][k] - 8.0 * ux[i][j-1][k] + ux[i][j-2][k];
	  du_dy /= 12.0 * h;
	  
	  du_dy_ex = - ux_ex[i][j+2][k] + 8.0 * ux_ex[i][j+1][k] - 8.0 * ux_ex[i][j-1][k] + ux_ex[i][j-2][k];
	  du_dy_ex /= 12.0 * h;
	}
	
	double dv_dy = 0.0;
	double dv_dy_ex = 0.0;
	  
	if ( j==0 || j==1 ) {
	  dv_dy = - (25.0/12.0) * uy[i][j][k] + 4.0 * uy[i][j+1][k] - 3.0 * uy[i][j+2][k]
	    + (4.0/3.0) * uy[i][j+3][k] - 0.25 * uy[i][j+4][k];
	  dv_dy /= h;
	  
	  dv_dy_ex = - (25.0/12.0) * uy_ex[i][j][k] + 4.0 * uy_ex[i][j+1][k] - 3.0 * uy_ex[i][j+2][k]
	    + (4.0/3.0) * uy_ex[i][j+3][k] - 0.25 * uy_ex[i][j+4][k];
	  dv_dy_ex /= h;
	}
	else if ( j==N || j==(N-1) ) {
	  dv_dy = (25.0/12.0) * uy[i][j][k] - 4.0 * uy[i][j-1][k] + 3.0 * uy[i][j-2][k]
	    - (4.0/3.0) * uy[i][j-3][k] + 0.25 * uy[i][j-4][k];
	  dv_dy /= h;
	  
	  dv_dy_ex = (25.0/12.0) * uy_ex[i][j][k] - 4.0 * uy_ex[i][j-1][k] + 3.0 * uy_ex[i][j-2][k]
	    - (4.0/3.0) * uy_ex[i][j-3][k] + 0.25 * uy_ex[i][j-4][k];
	  dv_dy_ex /= h;
	}
	else {
	  dv_dy = - uy[i][j+2][k] + 8.0 * uy[i][j+1][k] - 8.0 * uy[i][j-1][k] + uy[i][j-2][k];
	  dv_dy /= 12.0 * h;
	  
	  dv_dy_ex = - uy_ex[i][j+2][k] + 8.0 * uy_ex[i][j+1][k] - 8.0 * uy_ex[i][j-1][k] + uy_ex[i][j-2][k];
	  dv_dy_ex /= 12.0 * h;
	}
	
	double dw_dy = 0.0;
	double dw_dy_ex = 0.0;
	  
	if ( j==0 || j==1 ) {
	  dw_dy = - (25.0/12.0) * uz[i][j][k] + 4.0 * uz[i][j+1][k] - 3.0 * uz[i][j+2][k]
	    + (4.0/3.0) * uz[i][j+3][k] - 0.25 * uz[i][j+4][k];
	  dw_dy /= h;
	  
	  dw_dy_ex = - (25.0/12.0) * uz_ex[i][j][k] + 4.0 * uz_ex[i][j+1][k] - 3.0 * uz_ex[i][j+2][k]
	    + (4.0/3.0) * uz_ex[i][j+3][k] - 0.25 * uz_ex[i][j+4][k];
	  dw_dy_ex /= h;
	}
	else if ( j==N || j==(N-1) ) {
	  dw_dy = (25.0/12.0) * uz[i][j][k] - 4.0 * uz[i][j-1][k] + 3.0 * uz[i][j-2][k]
	    - (4.0/3.0) * uz[i][j-3][k] + 0.25 * uz[i][j-4][k];
	  dw_dy /= h;
	  
	  dw_dy_ex = (25.0/12.0) * uz_ex[i][j][k] - 4.0 * uz_ex[i][j-1][k] + 3.0 * uz_ex[i][j-2][k]
	    - (4.0/3.0) * uz_ex[i][j-3][k] + 0.25 * uz_ex[i][j-4][k];
	  dw_dy_ex /= h;
	}
	else {
	  dw_dy = - uz[i][j+2][k] + 8.0 * uz[i][j+1][k] - 8.0 * uz[i][j-1][k] + uz[i][j-2][k];
	  dw_dy /= 12.0 * h;
	  
	  dw_dy_ex = - uz_ex[i][j+2][k] + 8.0 * uz_ex[i][j+1][k] - 8.0 * uz_ex[i][j-1][k] + uz_ex[i][j-2][k];
	  dw_dy_ex /= 12.0 * h;
	}

	
	// Derivatives in z
	double du_dz = 0.0;
	double du_dz_ex = 0.0;
	  
	if ( k==0 || k==1 ) {
	  du_dz = - (25.0/12.0) * ux[i][j][k] + 4.0 * ux[i][j][k+1] - 3.0 * ux[i][j][k+2]
	    + (4.0/3.0) * ux[i][j][k+3] - 0.25 * ux[i][j][k+4];
	  du_dz /= h;
	  
	  du_dz_ex = - (25.0/12.0) * ux_ex[i][j][k] + 4.0 * ux_ex[i][j][k+1] - 3.0 * ux_ex[i][j][k+2]
	    + (4.0/3.0) * ux_ex[i][j][k+3] - 0.25 * ux_ex[i][j][k+4];
	  du_dz_ex /= h;
	}
	else if ( k==N || k==(N-1) ) {
	  du_dz = (25.0/12.0) * ux[i][j][k] - 4.0 * ux[i][j][k-1] + 3.0 * ux[i][j][k-2]
	    - (4.0/3.0) * ux[i][j][k-3] + 0.25 * ux[i][j][k-4];
	  du_dz /= h;
	  
	  du_dz_ex = (25.0/12.0) * ux_ex[i][j][k] - 4.0 * ux_ex[i][j][k-1] + 3.0 * ux_ex[i][j][k-2]
	    - (4.0/3.0) * ux_ex[i][j][k-3] + 0.25 * ux_ex[i][j][k-4];
	  du_dz_ex /= h;
	}
	else {
	  du_dz = - ux[i][j][k+2] + 8.0 * ux[i][j][k+1] - 8.0 * ux[i][j][k-1] + ux[i][j][k-2];
	  du_dz /= 12.0 * h;
	  
	  du_dz_ex = - ux_ex[i][j][k+2] + 8.0 * ux_ex[i][j][k+1] - 8.0 * ux_ex[i][j][k-1] + ux_ex[i][j][k-2];
	  du_dz_ex /= 12.0 * h;
	}
	
	double dv_dz = 0.0;
	double dv_dz_ex = 0.0;
	  
	if ( k==0 || k==1 ) {
	  dv_dz = - (25.0/12.0) * uy[i][j][k] + 4.0 * uy[i][j][k+1] - 3.0 * uy[i][j][k+2]
	    + (4.0/3.0) * uy[i][j][k+3] - 0.25 * uy[i][j][k+4];
	  dv_dz /= h;
	  
	  dv_dz_ex = - (25.0/12.0) * uy_ex[i][j][k] + 4.0 * uy_ex[i][j][k+1] - 3.0 * uy_ex[i][j][k+2]
	    + (4.0/3.0) * uy_ex[i][j][k+3] - 0.25 * uy_ex[i][j][k+4];
	  dv_dz_ex /= h;
	}
	else if ( k==N || k==(N-1) ) {
	  dv_dz = (25.0/12.0) * uy[i][j][k] - 4.0 * uy[i][j][k-1] + 3.0 * uy[i][j][k-2]
	    - (4.0/3.0) * uy[i][j][k-3] + 0.25 * uy[i][j][k-4];
	  dv_dz /= h;
	  
	  dv_dz_ex = (25.0/12.0) * uy_ex[i][j][k] - 4.0 * uy_ex[i][j][k-1] + 3.0 * uy_ex[i][j][k-2]
	    - (4.0/3.0) * uy_ex[i][j][k-3] + 0.25 * uy_ex[i][j][k-4];
	  dv_dz_ex /= h;
	}
	else {
	  dv_dz = - uy[i][j][k+2] + 8.0 * uy[i][j][k+1] - 8.0 * uy[i][j][k-1] + uy[i][j][k-2];
	  dv_dz /= 12.0 * h;
	  
	  dv_dz_ex = - uy_ex[i][j][k+2] + 8.0 * uy_ex[i][j][k+1] - 8.0 * uy_ex[i][j][k-1] + uy_ex[i][j][k-2];
	  dv_dz_ex /= 12.0 * h;
	}
	
	double dw_dz = 0.0;
	double dw_dz_ex = 0.0;
	  
	if ( k==0 || k==1 ) {
	  dw_dz = - (25.0/12.0) * uz[i][j][k] + 4.0 * uz[i][j][k+1] - 3.0 * uz[i][j][k+2]
	    + (4.0/3.0) * uz[i][j][k+3] - 0.25 * uz[i][j][k+4];
	  dw_dz /= h;
	  
	  dw_dz_ex = - (25.0/12.0) * uz_ex[i][j][k] + 4.0 * uz_ex[i][j][k+1] - 3.0 * uz_ex[i][j][k+2]
	    + (4.0/3.0) * uz_ex[i][j][k+3] - 0.25 * uz_ex[i][j][k+4];
	  dw_dz_ex /= h;
	}
	else if ( k==N || k==(N-1) ) {
	  dw_dz = (25.0/12.0) * uz[i][j][k] - 4.0 * uz[i][j][k-1] + 3.0 * uz[i][j][k-2]
	    - (4.0/3.0) * uz[i][j][k-3] + 0.25 * uz[i][j][k-4];
	  dw_dz /= h;
	  
	  dw_dz_ex = (25.0/12.0) * uz_ex[i][j][k] - 4.0 * uz_ex[i][j][k-1] + 3.0 * uz_ex[i][j][k-2]
	    - (4.0/3.0) * uz_ex[i][j][k-3] + 0.25 * uz_ex[i][j][k-4];
	  dw_dz_ex /= h;
	}
	else {
	  dw_dz = - uz[i][j][k+2] + 8.0 * uz[i][j][k+1] - 8.0 * uz[i][j][k-1] + uz[i][j][k-2];
	  dw_dz /= 12.0 * h;
	  
	  dw_dz_ex = - uz_ex[i][j][k+2] + 8.0 * uz_ex[i][j][k+1] - 8.0 * uz_ex[i][j][k-1] + uz_ex[i][j][k-2];
	  dw_dz_ex /= 12.0 * h;
	}	

	double err_temp = sqrt( (du_dx - du_dx_ex) * (du_dx - du_dx_ex) +
				(dv_dx - dv_dx_ex) * (dv_dx - dv_dx_ex) +
				(dw_dx - dw_dx_ex) * (dw_dx - dw_dx_ex) +
				(du_dy - du_dy_ex) * (du_dy - du_dy_ex) +
				(dv_dy - dv_dy_ex) * (dv_dy - dv_dy_ex) +
				(dw_dy - dw_dy_ex) * (dw_dy - dw_dy_ex)	+
				(du_dz - du_dz_ex) * (du_dz - du_dz_ex) +
				(dv_dz - dv_dz_ex) * (dv_dz - dv_dz_ex) +
				(dw_dz - dw_dz_ex) * (dw_dz - dw_dz_ex)	);

	  if ( err_temp > err_max ) {
	    err_max = err_temp;
	    max_i = i;
	    max_j = j;
	    max_k = k;
	  }
	  err_l2 += err_temp * err_temp;
	  count++;
	
	if ( abs(j*h - 1.5) < 1.e-8 ) {
	  err_Du_file << setprecision(16) << sqrt(err_temp) << "  ";
	}

      } // for k
    } // for j
    err_Du_file << endl;
  } // for i
  err_Du_file.close();

  err_l2 = sqrt(err_l2 / count);
  cout << "h=1/" << N/3 << "   max err (grad u) = " << err_max
       << "    l2 err (grad u) = " << err_l2 << endl << endl;
  cout << "max i = " << max_i << "  max j = " << max_j << "  max k = " << max_k << endl;
  

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      delete[] ux[i][j];
      delete[] uy[i][j];
      delete[] uz[i][j];
    }
    delete[] ux[i];
    delete[] uy[i];
    delete[] uz[i];
  }
  delete[] ux;
  delete[] uy;
  delete[] uz;
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      delete[] ux_ex[i][j];
      delete[] uy_ex[i][j];
      delete[] uz_ex[i][j];
    }
    delete[] ux_ex[i];
    delete[] uy_ex[i];
    delete[] uz_ex[i];
  }
  delete[] ux_ex;
  delete[] uy_ex;
  delete[] uz_ex;
  
}

//*********************************************************************
//*********************************************************************

void solve_Poisson(int N, double L, double h,
		   int*** nrst,
		   int N_target, const vector<Target_point>& Target,
		   double*** p,
		   double*** u,
		   int variable) {
  
  double*** w = new double**[N+1];
  for (int i=0; i<=N; i++) {
    w[i] = new double*[N+1];
    for (int j=0; j<=N; j++) {
      w[i][j] = new double[N+1];
    }
  }

  // compute G as extension of boundary values of ux
  extend_into_cube(N, L, h, u, w);
  
  int Nm1_cubed = pow(N-1,3);
  double* rhs = new double[Nm1_cubed];
  for (int i=0; i<Nm1_cubed; i++)  rhs[i] = 0.0; 

  // generate right hand side = -h^2 * Delta^2_h (w)
  lapl_4ord_inCube(N, h, w, rhs);		      

  
  // add to right hand side  h^2*f
  lapl_4ord_nearSurfc(N, h, nrst, Target, u, rhs);		      


  if ( variable == 1 ) { // first component of velocity
    // add to right hand side  h^2 * dp_dx
    dp_dx_4ord(N, h, nrst, Target, p, rhs);
  }
  else if ( variable == 2 ) {
    // add to right hand side  h^2 * dp_dy
    dp_dy_4ord(N, h, nrst, Target, p, rhs);
  }
  else if ( variable == 3 ) {
    // add to right hand side  h^2 * dp_dz
    dp_dz_4ord(N, h, nrst, Target, p, rhs);
  }

  
  fftw_DST_Poisson3d_FD(N-1, rhs);
    

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {
	u[i][j][k] = w[i][j][k];
      }
    }
  }
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	int ijk = ((i-1) * (N-1) + j-1) * (N-1) + k-1;
	u[i][j][k] += rhs[ijk];   // the solution 
      }
    }
  }

  // deallocate memory
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      delete[] w[i][j];
    }
    delete[] w[i];
  }
  delete[] w;

  delete[] rhs;
  
}

//*********************************************************************
//*********************************************************************

void eval_p_BI_cubeEdges(int N, double h,
			 int N_quad, const vector<Surf_point>& Surfc,
			 const vector<double>& f_ex,
			 double*** p_bi) {
  
#pragma omp parallel
  {
#pragma omp for
    
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	p_bi[i][j][k] = 0.0;

	if ( i==0 || i==N || j==0 || j==N || k==0 || k==N ) { 

	  double x = i * h;
	  double y = j * h;
	  double z = k * h;

	  // compute pressure
	  p_bi[i][j][k] = stokes_pressure_away(N_quad, Surfc,
					       f_ex,
					       x, y, z);
	  
	}
      }
    }
  }
  } // end pragma omp parallel
  
  double err_max_p = 0.0;
  double err_l2_p = 0.0;
  double val_max_p = 0.0;
  double val_l2_p = 0.0;
  int count_p = 0;

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	if ( i==0 || i==N || j==0 || j==N || k==0 || k==N ) { 

	  double x = i * h;
	  double y = j * h;
	  double z = k * h;

	  // compute error max and l2 norms
	  double p_ex = p_exact(x, y, z);
	  double err_temp = abs(p_bi[i][j][k] - p_ex);
	  if ( err_temp > err_max_p ) {
	    err_max_p = err_temp;
	    val_max_p = p_ex;
	  }
	  err_l2_p += err_temp * err_temp;
	  val_l2_p += p_ex * p_ex;
	  count_p++;	  
	}
      }
    }
  }
  cout << endl << "BI for pressure on box edges (" << count_p << " points)" << endl;
  cout << "--- max error = " << err_max_p << "   l2 error = " << sqrt(err_l2_p/count_p)
       << endl
       << "--- with exact = " << val_max_p
       << "   l2 (exact) = " << sqrt(val_l2_p/count_p) << endl << endl;  
}

//*********************************************************************

void eval_p_BI_nearSurfc(int N, double h, double DEL,
			 int N_quad, const vector<Surf_point>& Surfc,
			 int*** nrst,
			 int N_target, const vector<Target_point>& Target,
			 const vector<double>& f_ex,
			 double*** p_bi) {

#pragma omp parallel
  {
  vector<double> x0(3,0), n0(3,0), f0(3,0);

#pragma omp for
    
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	if ( nrst[i][j][k] != -1 ) {
	  double x = i * h;
	  double y = j * h;
	  double z = k * h;
	  int ind = nrst[i][j][k];
  
	  double b = Target[ind].b; // signed distance from interface

	  if ( (b > 0.0) || ( abs(b) < 1.e-12) ) { //outside or on surface only

	    x0[0] = Target[ind].nrst_x;
	    x0[1] = Target[ind].nrst_y;
	    x0[2] = Target[ind].nrst_z;
	    n0[0] = Target[ind].nrst_Nrml[0];
	    n0[1] = Target[ind].nrst_Nrml[1];
	    n0[2] = Target[ind].nrst_Nrml[2];
	    
	    f0 = stokeslet_density(x0);

	    // compute pressure
	    p_bi[i][j][k] = stokes_pressure_7ord(DEL,
						 N_quad, Surfc,
						 f_ex, f0, n0,
						 x, y, z, b);
	    
	  } // if outside or on surface
	  
	} // if ( nrst[i][j][k] != -1 )
      } // for k
    } // for j
  } // for i
  } // end pragma omp parallel

  
  double err_max_p = 0.0;
  double err_l2_p = 0.0;
  double val_max_p = 0.0;
  double val_l2_p = 0.0;
  int count_p = 0;

  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	if ( nrst[i][j][k] != -1 ) {
	  double x = i * h;
	  double y = j * h;
	  double z = k * h;
	  int ind = nrst[i][j][k];
  
	  double b = Target[ind].b; // signed distance from interface

	  if ( (b > 0.0) || ( abs(b) < 1.e-12) ) { //outside or on surface only

	    // compute error max and l2 norms
	    double p_ex = p_exact(x, y, z);
	    double err_temp = abs(p_bi[i][j][k] - p_ex);
	    if ( err_temp > err_max_p ) {
	      err_max_p = err_temp;
	      val_max_p = p_ex;
	    }
	    err_l2_p += err_temp * err_temp;
	    val_l2_p += p_ex * p_ex;
	    count_p++;

	  } // if outside or on surface
	  
	} // if ( nrst[i][j][k] != -1 )
      } // for k
    } // for j
  } // for i 	  
  cout << endl << "BI for pressure in BAND1 around interface (" << count_p << " points)" << endl;
  cout << "--- max error = " << err_max_p << "   l2 error = " << sqrt(err_l2_p/count_p)
       << endl
       << "--- with exact = " << val_max_p
       << "   l2 (exact) = " << sqrt(val_l2_p/count_p) << endl << endl;
}

//*********************************************************************
//*********************************************************************

void eval_u_BI_cubeEdges(int N, double h,
			 int N_quad, const vector<Surf_point>& Surfc,
			 const vector<double>& f_ex,
			 double*** ux_bi, double*** uy_bi, double*** uz_bi) {


#pragma omp parallel
  {
  vector<double> SL(3,0);
  
#pragma omp for
    
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	ux_bi[i][j][k] = 1.0; // SHOULD BE THE TRANSLATION VELOCITY
	uy_bi[i][j][k] = 0.0;
	uz_bi[i][j][k] = 0.0;

	if ( i==0 || i==N || j==0 || j==N || k==0 || k==N ) { 

	  double x = i * h;
	  double y = j * h;
	  double z = k * h;

	  // compute velocity
	  stokes_SL_away(N_quad, Surfc,
			 f_ex,
			 x, y, z,
			 SL);
	  ux_bi[i][j][k] = SL[0];
	  uy_bi[i][j][k] = SL[1];
	  uz_bi[i][j][k] = SL[2];	  	  
	}
      }
    }
  }
  } // end pragma omp parallel
  
  double err_max_u = 0.0;
  double err_l2_u = 0.0;
  double val_max_u = 0.0;
  double val_l2_u = 0.0;
  int count_u = 0;

  vector<double> SL_ex(3,0);
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {	

	if ( i==0 || i==N || j==0 || j==N || k==0 || k==N ) { 

	  double x = i * h;
	  double y = j * h;
	  double z = k * h;
	  
	  SL_ex = SL_exact(x, y, z);
	  double ex_norm = SL_ex[0]*SL_ex[0] + SL_ex[1]*SL_ex[1] + SL_ex[2]*SL_ex[2];
	  double err_temp = sqrt( (ux_bi[i][j][k] - SL_ex[0]) * (ux_bi[i][j][k] - SL_ex[0])
			        + (uy_bi[i][j][k] - SL_ex[1]) * (uy_bi[i][j][k] - SL_ex[1])
			        + (uz_bi[i][j][k] - SL_ex[2]) * (uz_bi[i][j][k] - SL_ex[2]) );
	  
	  if ( err_temp > err_max_u ) {
	    err_max_u = err_temp;
	    val_max_u = sqrt(ex_norm);
	  }
	  err_l2_u += err_temp * err_temp;
	  val_l2_u += ex_norm;
	  count_u++;
	  
	}
      }
    }
  }
  cout << endl << "BI for velocity on box edges (" << count_u << " points)" << endl;
  cout << "--- max error = " << err_max_u << "   l2 error = " << sqrt(err_l2_u/count_u)
       << endl
       << "--- with exact = " << val_max_u
       << "   l2 (exact) = " << sqrt(val_l2_u/count_u)<< endl << endl;
}

//*********************************************************************

void eval_u_BI_nearSurfc(int N, double h, double DEL,
			 int N_quad, const vector<Surf_point>& Surfc,
			 int*** nrst,
			 int N_target, const vector<Target_point>& Target,
			 const vector<double>& f_ex,
			 double*** ux_bi, double*** uy_bi, double*** uz_bi) {

#pragma omp parallel
  {
  vector<double> x0(3,0), n0(3,0), f0(3,0), SL(3,0);

#pragma omp for
    
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	if ( nrst[i][j][k] != -1 ) {
	  double x = i * h;
	  double y = j * h;
	  double z = k * h;
	  int ind = nrst[i][j][k];
  
	  double b = Target[ind].b; // signed distance from interface

	  if ( (b > 0.0) || ( abs(b) < 1.e-12) ) { //outside or on surface only

	    x0[0] = Target[ind].nrst_x;
	    x0[1] = Target[ind].nrst_y;
	    x0[2] = Target[ind].nrst_z;
	    n0[0] = Target[ind].nrst_Nrml[0];
	    n0[1] = Target[ind].nrst_Nrml[1];
	    n0[2] = Target[ind].nrst_Nrml[2];
	    
	    f0 = stokeslet_density(x0);
	    double f0_dot_n0 = dot_product(f0, n0);

	    // compute velocity
	    stokes_SL_7ord(DEL,
			   N_quad, Surfc,
			   f_ex, f0_dot_n0,
			   x, y, z, b,
			   SL);
	    ux_bi[i][j][k] = SL[0];
	    uy_bi[i][j][k] = SL[1];
	    uz_bi[i][j][k] = SL[2];
	    
	  } // if outside or on surface
	  
	} // if ( nrst[i][j][k] != -1 )
      } // for k
    } // for j
  } // for i
  } // end pragma omp parallel

  
  double err_max_u = 0.0;
  double err_l2_u = 0.0;
  double val_max_u = 0.0;
  double val_l2_u = 0.0;
  int count_u = 0;
  
  vector<double> SL_ex(3,0);

  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	if ( nrst[i][j][k] != -1 ) {
	  double x = i * h;
	  double y = j * h;
	  double z = k * h;
	  int ind = nrst[i][j][k];
  
	  double b = Target[ind].b; // signed distance from interface

	  if ( (b > 0.0) || ( abs(b) < 1.e-12) ) { //outside or on surface only
	    
	    SL_ex = SL_exact(x, y, z);
	    double ex_norm = SL_ex[0]*SL_ex[0] + SL_ex[1]*SL_ex[1] + SL_ex[2]*SL_ex[2];
	    double err_temp = sqrt( (ux_bi[i][j][k] - SL_ex[0]) * (ux_bi[i][j][k] - SL_ex[0])
		         	  + (uy_bi[i][j][k] - SL_ex[1]) * (uy_bi[i][j][k] - SL_ex[1])
			          + (uz_bi[i][j][k] - SL_ex[2]) * (uz_bi[i][j][k] - SL_ex[2]) );
    
	    if ( err_temp > err_max_u ) {
	      err_max_u = err_temp;
	      val_max_u = sqrt(ex_norm);
	    }
	    err_l2_u += err_temp * err_temp;
	    val_l2_u += ex_norm;
	    count_u++;

	  } // if outside or on surface
	  
	} // if ( nrst[i][j][k] != -1 )
      } // for k
    } // for j
  } // for i 
  cout << endl << "BI for velocity in BAND1 around interface (" << count_u << " points)" << endl;
  cout << "--- max error = " << err_max_u << "   l2 error = " << sqrt(err_l2_u/count_u)
       << endl
       << "--- with exact = " << val_max_u
       << "   l2 (exact) = " << sqrt(val_l2_u/count_u)<< endl << endl;  
}

//*********************************************************************
//*********************************************************************

// solves the Poisson equation nabla(u) = f (real f) on a 3D box [0,L]^3
// with Dirichlet B.C. u=0 on boundary of box
// uses the Discrete Sine Transform

// uses 4nd order finite differences from Iserles' book on Numerical Methods
// (p.169, exercise 8.11) to form a discrete Laplacian

// input: N = size of input/output vectors (N=n-1 where h=L/n)
//        f = h^2 * f_{i,j,k} : real-valued
// output: u = discrete solution on the 3D grid in a row-major format
//           = solution at interior points of cube

void fftw_DST_Poisson3d_FD(int N, double* f) {

  int N3 = N * N * N;
  
  double *in, *out;
  in  = (double*) fftw_malloc(sizeof(double) * N3);
  out = (double*) fftw_malloc(sizeof(double) * N3);

  fftw_plan pf, pb;
  pf = fftw_plan_r2r_3d(N, N, N, in, in,
			FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00,
			FFTW_ESTIMATE);
  pb = fftw_plan_r2r_3d(N, N, N, in, out,
			FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00,
			FFTW_ESTIMATE);


  for (int i=0; i<N3; i++) {
    in[i] = f[i];
  }
  fftw_execute(pf);

  double pi_over_n = PI / (N+1);
    
  for (int i=0; i<N; i++) {
    double c1 = cos((i+1) * pi_over_n);
    double s1 = sin((i+1) * pi_over_n / 2.0);
    s1 *= -2.0 * s1;

    for (int j=0; j<N; j++) {
      double c2 = cos((j+1) * pi_over_n);
      double s2 = sin((j+1) * pi_over_n / 2.0);
      s2 *= -2.0 * s2;

      for (int k=0; k<N; k++) {
	double c3 = cos((k+1) * pi_over_n);
	double s3 = sin((k+1) * pi_over_n / 2.0);
	s3 *= -2.0 * s3;

	int ind = (i * N + j) * N + k;

	// compute e-value factor for the laplacian
	// Laplacian 4th order finite difference from Iserles' book on Numerical Methods
	// (p.169, exercise 8.11)
	double lap7 = s1 + s2 + s3;
	double sum8 = c1 * c2 * c3;
	double factor = 2.0/3.0 * (2.0 * lap7 + sum8 - 1.0);
	  
	if ( abs(factor) > 1.e-14 ) { 
	  in[ind] /= factor;
	} else {
	  in[ind] = 0.0;  // DFT[0] = average of input values = 0
	}
	
      }
    }
  }
  fftw_execute(pb);

  double scale = 1./pow(2*(N+1),3); // the scaling for DST-I
  
  for (int i=0; i<N3; i++) {
    f[i] = out[i] * scale;  // FFTW doesn't scale the solution so scale it
  }
  
  fftw_destroy_plan(pf); fftw_destroy_plan(pb); fftw_cleanup();
}

//*********************************************************************
//*********************************************************************

void extend_into_cube(int N, double L, double h,
		      double*** g, double*** G) {

  for (int i=0; i<=N; i++) {
    double x = i * h / L;
    for (int j=0; j<=N; j++) {
      double y = j * h / L;
      for (int k=0; k<=N; k++) {	
	double z = k * h / L;

	double F, E1, E2, E3, C;
	
	F  = x * g[N][j][k] + (1.0-x) * g[0][j][k];
	F += y * g[i][N][k] + (1.0-y) * g[i][0][k];
	F += z * g[i][j][N] + (1.0-z) * g[i][j][0];

	E1 = (1.0-y) * (1.0-z) * g[i][0][0] + (1.0-y) * z * g[i][0][N]
	   +      y  * (1.0-z) * g[i][N][0] +      y  * z * g[i][N][N];
	
	E2 = (1.0-x) * (1.0-z) * g[0][j][0] + (1.0-x) * z * g[0][j][N]
	   +      x  * (1.0-z) * g[N][j][0] +      x  * z * g[N][j][N];
	
 	E3 = (1.0-x) * (1.0-y) * g[0][0][k] + (1.0-x) * y * g[0][N][k]
	   +      x  * (1.0-y) * g[N][0][k] +      x  * y * g[N][N][k];

	C = (1.0-x) * (1.0-y) * (1.0-z) * g[0][0][0] + (1.0-x) * (1.0-y) * z * g[0][0][N]
	  + (1.0-x) *      y  * (1.0-z) * g[0][N][0] + (1.0-x) *      y  * z * g[0][N][N]
	  +      x  * (1.0-y) * (1.0-z) * g[N][0][0] +      x  * (1.0-y) * z * g[N][0][N]
	  +      x  *      y  * (1.0-z) * g[N][N][0] +      x  *      y  * z * g[N][N][N];
	  
	G[i][j][k] = F - E1 - E2 - E3 + C;	
      }
    }
  }

  double err_corners = pow(G[0][0][0] - g[0][0][0],2) + pow(G[0][0][N] - g[0][0][N],2)
    + pow(G[0][N][0] - g[0][N][0],2) + pow(G[0][N][N] - g[0][N][N],2);

  //cout << "Extension error at corners = " << sqrt(err_corners/4) << endl;
  
  double err_edges = 0.0;
  for (int i=0; i<=N; i++) {
    err_edges += pow(G[i][0][0] - g[i][0][0],2) + pow(G[i][0][N] - g[i][0][N],2)
      + pow(G[i][N][0] - g[i][N][0],2) + pow(G[i][N][N] - g[i][N][N],2);

    err_edges += pow(G[0][i][0] - g[0][i][0],2) + pow(G[0][i][N] - g[0][i][N],2)
      + pow(G[N][i][0] - g[N][i][0],2) + pow(G[N][i][N] - g[N][i][N],2);
    
    err_edges += pow(G[0][0][i] - g[0][0][i],2) + pow(G[0][N][i] - g[0][N][i],2)
      + pow(G[N][0][i] - g[N][0][i],2) + pow(G[N][N][i] - g[N][N][i],2);
    
  }
  //cout << "Extension error on edges = " << sqrt(err_edges/12/N) << endl;

  double err_faces = 0.0;
  
  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      err_faces += pow(G[i][j][0] - g[i][j][0],2) + pow(G[i][j][N] - g[i][j][N],2)
	+ pow(G[0][i][j] - g[0][i][j],2) + pow(G[N][i][j] - g[N][i][j],2)
	+ pow(G[i][0][j] - g[i][0][j],2) + pow(G[i][N][j] - g[i][N][j],2);
    }
  }
  //cout << "Extension error on faces = " << sqrt(err_faces/6/N/N) << endl;
  
}

//*********************************************************************

// compute 4th order discrete laplacian with a factor
// return  -h^2 * Delta_h (G)

void lapl_4ord_inCube(int N, double h,
		      double*** G,
		      double* rhs) {
  
  double lap_err = 0.0;
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	int ijk = ((i-1) * (N-1) + j-1) * (N-1) + k-1;

	double lap7 = G[i+1][j][k] + G[i-1][j][k] + G[i][j+1][k] + G[i][j-1][k]
	  + G[i][j][k+1] + G[i][j][k-1] - 6.0 * G[i][j][k];

	double temp1 = G[i+1][j+1][k+1] + G[i+1][j-1][k+1]  
	  + G[i-1][j+1][k+1] + G[i-1][j-1][k+1] - 4.0 * G[i][j][k];

	double temp2 = G[i+1][j+1][k-1] + G[i+1][j-1][k-1]  
	  + G[i-1][j+1][k-1] + G[i-1][j-1][k-1] - 4.0 * G[i][j][k];

	rhs[ijk] = - 2.0/3.0 * (lap7 + (temp1 + temp2)/8.0 ); //scaled by (-h^2)
	
	lap_err += rhs[ijk] * rhs[ijk];
      }
    }
  }
  //cout << "Laplacian * h^2 in whole box --- error = "
  //     << sqrt(lap_err/pow(N-1,3)) << endl;
}

//*********************************************************************

void lapl_4ord_nearSurfc(int N, double h,
			 int*** nrst, const vector<Target_point>& Target,
			 double*** u, double* rhs) {

  double lap_err = 0.0;
  int count = 0;
  
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	int ijk = ((i-1) * (N-1) + j-1) * (N-1) + k-1;

	if ( nrst[i][j][k] != -1 ) {

	  int ind = nrst[i][j][k];
	  if ( Target[ind].band2 ) {
	    
	    double lap7 = u[i+1][j][k] + u[i-1][j][k]
	      + u[i][j+1][k] + u[i][j-1][k]
	      + u[i][j][k+1] + u[i][j][k-1] - 6.0 * u[i][j][k];
	    
	    double temp1 = u[i+1][j+1][k+1] + u[i+1][j-1][k+1]  
	      + u[i-1][j+1][k+1] + u[i-1][j-1][k+1] - 4.0 * u[i][j][k];
	    
	    double temp2 = u[i+1][j+1][k-1] + u[i+1][j-1][k-1]  
	      + u[i-1][j+1][k-1] + u[i-1][j-1][k-1] - 4.0 * u[i][j][k];
	    
	    rhs[ijk] += 2.0/3.0 * (lap7 + (temp1 + temp2)/8.0 ); //scaled by h^2

	    lap_err += rhs[ijk] * rhs[ijk];
	    count++;
	  }
	}
      }
    }
  }
  //cout << "Laplacian * h^2 in BAND2 --- error = "
  //     << sqrt(lap_err/count) << endl;
}

//*********************************************************************

void dp_dx_4ord(int N, double h,
		int*** nrst, const vector<Target_point>& Target,
		double*** p, double* rhs) {
  
  int count = 0;
  vector<double> x_vec(3,0);
  
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	int ijk = ((i-1) * (N-1) + j-1) * (N-1) + k-1;

	bool out_band2 = false;
	
	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	  
	double levl = phi(x_vec);
	
	if ( levl > 0 ) { // outside surface only
	  if ( nrst[i][j][k] == -1 ) { // outside band1
	    out_band2 = true;
	  }
	  else {
	    int ind = nrst[i][j][k];
	    if ( !(Target[ind].band2) ) { // inside band1 but outside band2
	      out_band2 = true;
	    }
	  }
	}

	if ( out_band2 ) {

	  double dp_dx;

	  if ( i==0 || i==1 ) {
	    dp_dx = - (25.0/12.0) * p[i][j][k] + 4.0 * p[i+1][j][k] - 3.0 * p[i+2][j][k]
	      + (4.0/3.0) * p[i+3][j][k] - 0.25 * p[i+4][j][k];
	    //dp_dx /= h;
	  }
	  else if ( i==N || i==(N-1) ) {
	    dp_dx = (25.0/12.0) * p[i][j][k] - 4.0 * p[i-1][j][k] + 3.0 * p[i-2][j][k]
	    - (4.0/3.0) * p[i-3][j][k] + 0.25 * p[i-4][j][k];
	    //dp_dx /= h;
	  }
	  else {
	    dp_dx = - p[i+2][j][k] + 8.0 * p[i+1][j][k] - 8.0 * p[i-1][j][k] + p[i-2][j][k];
	    dp_dx /= 12.0;
	    //dp_dx /= 12.0 * h;
	  }
	    
	  rhs[ijk] += dp_dx * h; //scaled by h^2 (h scale already in dp_dx)

	  count++;
	}
      }
    }
  }
}

//*********************************************************************

void dp_dy_4ord(int N, double h,
		int*** nrst, const vector<Target_point>& Target,
		double*** p, double* rhs) {
  
  int count = 0;
  vector<double> x_vec(3,0);
  
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	int ijk = ((i-1) * (N-1) + j-1) * (N-1) + k-1;

	bool out_band2 = false;
	
	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	  
	double levl = phi(x_vec);
	
	if ( levl > 0 ) { // outside surface only
	  if ( nrst[i][j][k] == -1 ) { // outside band1
	    out_band2 = true;
	  }
	  else {
	    int ind = nrst[i][j][k];
	    if ( !(Target[ind].band2) ) { // inside band1 but outside band2
	      out_band2 = true;
	    }
	  }
	}

	if ( out_band2 ) {

	  double dp_dy;

	  if ( j==0 || j==1 ) {
	    dp_dy = - (25.0/12.0) * p[i][j][k] + 4.0 * p[i][j+1][k] - 3.0 * p[i][j+2][k]
	      + (4.0/3.0) * p[i][j+3][k] - 0.25 * p[i][j+4][k];
	  }
	  else if ( j==N || j==(N-1) ) {
	    dp_dy = (25.0/12.0) * p[i][j][k] - 4.0 * p[i][j-1][k] + 3.0 * p[i][j-2][k]
	    - (4.0/3.0) * p[i][j-3][k] + 0.25 * p[i][j-4][k];
	  }
	  else {
	    dp_dy = - p[i][j+2][k] + 8.0 * p[i][j+1][k] - 8.0 * p[i][j-1][k] + p[i][j-2][k];
	    dp_dy /= 12.0;
	  }
	    
	  rhs[ijk] += dp_dy * h; //scaled by h^2 (h scale already in dp_dy)

	  count++;
	}
      }
    }
  }
}

//*********************************************************************

void dp_dz_4ord(int N, double h,
		int*** nrst, const vector<Target_point>& Target,
		double*** p, double* rhs) {
  
  int count = 0;
  vector<double> x_vec(3,0);
  
  for (int i=1; i<N; i++) {
    for (int j=1; j<N; j++) {
      for (int k=1; k<N; k++) {	

	int ijk = ((i-1) * (N-1) + j-1) * (N-1) + k-1;

	bool out_band2 = false;
	
	x_vec[0] = i*h;
	x_vec[1] = j*h;
	x_vec[2] = k*h;
	  
	double levl = phi(x_vec);
	
	if ( levl > 0 ) { // outside surface only
	  if ( nrst[i][j][k] == -1 ) { // outside band1
	    out_band2 = true;
	  }
	  else {
	    int ind = nrst[i][j][k];
	    if ( !(Target[ind].band2) ) { // inside band1 but outside band2
	      out_band2 = true;
	    }
	  }
	}

	if ( out_band2 ) {

	  double dp_dz;

	  if ( k==0 || k==1 ) {
	    dp_dz = - (25.0/12.0) * p[i][j][k] + 4.0 * p[i][j][k+1] - 3.0 * p[i][j][k+2]
	      + (4.0/3.0) * p[i][j][k+3] - 0.25 * p[i][j][k+4];
	  }
	  else if ( k==N || k==(N-1) ) {
	    dp_dz = (25.0/12.0) * p[i][j][k] - 4.0 * p[i][j][k-1] + 3.0 * p[i][j][k-2]
	    - (4.0/3.0) * p[i][j][k-3] + 0.25 * p[i][j][k-4];
	  }
	  else {
	    dp_dz = - p[i][j][k+2] + 8.0 * p[i][j][k+1] - 8.0 * p[i][j][k-1] + p[i][j][k-2];
	    dp_dz /= 12.0;
	  }
	    
	  rhs[ijk] += dp_dz * h; //scaled by h^2 (h scale already in dp_dz)

	  count++;
	}
      }
    }
  }
}
