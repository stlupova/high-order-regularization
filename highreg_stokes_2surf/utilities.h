#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

static const int use_tree = 1; //1: use treecode

static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double denom_EightPI = 1.0/(8.0*PI);

static const double theta = 70.0*PI/180.0; //angle in Beale et al (2.2)
static const double tol = 1e-14; //tolerance in the search of quadrature points (Newton's/bisection method)
static const double ellipse_a = 1.0;
static const double ellipse_b = 0.5;
static const double ellipse_c = 0.5;

static const double ctr1_x = 0.0;
static const double ctr1_y = 0.0;
static const double ctr1_z = -0.5 - 1.0/(16.0*16.0*16.0);
//static const double ctr1_z = -0.5 - 1.0/16.0;

static const double ctr2_x = 0.0;
static const double ctr2_y = 0.0;
static const double ctr2_z = 0.5;

using namespace std;

struct Surf_point
{
  Surf_point() : x(3,0), f(3,0), g(3,0), Nrml(3,0), Area(0), plane(0), left(0), right(0), up(0), down(0) {}
  vector<double> x;
  vector<double> f;
  vector<double> g;
  vector<double> Nrml;
  double Area;
  int plane, left, right, up, down;
};

struct Target_point
{
  Target_point() : x(0), y(0), z(0), nrst_x(0), nrst_y(0), nrst_z(0), nrst_Nrml(3,0), nrst_f(3,0), nrst_g(3,0), u(0), v(0), w(0), nrst(0), i1(0), i2(0), i3(0), b(0) {}
  double x, y, z;
  double nrst_x, nrst_y, nrst_z;
  vector<double> nrst_Nrml, nrst_f, nrst_g;
  double u, v, w;
  int nrst;
  int i1;
  int i2;
  int i3;
  double b;
};

vector<double> unrotate(const vector<double>& x);
double phi(int surf_num, int rotated, const vector<double>& x);
double Dphi(int surf_num, int rotated, int i, const vector<double>& x);
void D2phi(int surf_num, int rotated, const vector<double>& x,
	   double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23,
	   double& phi31, double& phi32, double& phi33);

void get_reg_coeff_7ord(double lam, double& c1, double& c2, double& c3);

void Stokes_SL_7ord(double r, double d,
		    double c1, double c2, double c3,
		    double& s1, double& s2);

void Stokes_SL_7ord_onSurf(double r, double d,
			   double& s1, double& s2);

void Stokes_DL_7ord(double r, double d,
		    double c1, double c2, double c3,
		    double& s2, double& s3);

void Stokes_DL_7ord_onSurf(double r, double d,
			   double& s3);

void evaluate_In_024(double lam, double& I0, double& I2, double& I4);


double distance(const vector<double>& x, const vector<double>& y,
		const vector<double>& normal);
double dot_product(const vector<double>& x, const vector<double>& y);

int LU_factorization(vector< vector<double> >& A,
		     vector<int>& p);
vector<double> solveByLU(const vector< vector<double> >& A,
			 const vector<int>& p,
			 const vector<double>& b);

void initialize_vector(int n, vector<double>& vec);
