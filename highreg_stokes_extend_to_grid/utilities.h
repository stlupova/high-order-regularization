#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double denom_FourPI = 1.0/(4.0*PI);
static const double denom_EightPI = 1.0/(8.0*PI);

static const double ellipse_a = 1.0;
static const double ellipse_b = 0.5;
static const double ellipse_c = 0.5;
static const bool surfc_is_sphere = false;

//set to true if sphere or ellipsoid; if false, a molecule is used
static const bool surfc_elliptical = true;
static const double ctr = 1.5; // shift in the center of the surface

struct Surf_point
{
  Surf_point() : x(3,0), Nrml(3,0), Area(0) {}
  vector<double> x;
  vector<double> Nrml;
  double Area;
};

struct Target_point
{
  Target_point() : x(0), y(0), z(0), nrst_x(0), nrst_y(0), nrst_z(0), b(0), nrst_Nrml(3,0), band2(false) {}
  double x, y, z;
  double nrst_x, nrst_y, nrst_z;
  double b;
  vector<double> nrst_Nrml;
  bool band2;
};

double phi(const vector<double>& x);
double Dphi(int i, const vector<double>& x);
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33);

void get_reg_coeff_7ord(double lam, double& c1, double& c2, double& c3);

void Lapl_DL_7ord(double r, double d,
		  double c1, double c2, double c3,
		  double& s2);

void Stokes_SL_7ord(double r, double d,
		    double c1, double c2, double c3,
		    double& s1, double& s2);

void evaluate_In_024(double lam, double& I0, double& I2, double& I4);


vector<double> stokeslet_density(const vector<double>& pt);

double p_exact(double x_in, double y_in, double z_in);
double dp_dx_exact(double x_in, double y_in, double z_in);

vector<double> SL_exact(double x_in, double y_in, double z_in);

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal);

double dot_product(const vector<double>& x,
		   const vector<double>& y);

vector<double> cross_product(const vector<double>& x,
			     const vector<double>& y);

void orthogonalize(vector<double>& x,
		   vector<double>& y);

int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p);

vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b);

void initialize_vector(int n, vector<double>& vec);
