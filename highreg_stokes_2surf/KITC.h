
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>

#include "surface.h"

using namespace std;

static const int P = 6; // degree of Barycentric polynomial
static const int Pflat = (P + 1)*(P + 1)*(P + 1);
static const int N0 = 1000; // leaf size
static const double sq_theta = 0.6*0.6; // theta squared

//**********************************************************//

struct xyz // particle coordinates (physical)
{
  double* x;
  double* y;
  double* z;
  double* nx;
  double* ny;
  double* nz;
  double* area;
  size_t* index;
  size_t* old_index;
  size_t size;
  xyz(size_t N_cube_in)
  {
    size = N_cube_in;
    x = new double[size];
    y = new double[size];
    z = new double[size];
    nx = new double[size];
    ny = new double[size];
    nz = new double[size];
    area = new double[size];
    index = new size_t[size];
    old_index = new size_t[size];
  }
  ~xyz()
  {
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] nx;
    delete[] ny;
    delete[] nz;
    delete[] area;
    delete[] index;
    delete[] old_index;
  }
};

//**********************************************************//

struct panel
{
  size_t members[2];
  double xinterval[2];
  double yinterval[2];
  double zinterval[2];
  double xc; // panel center x coordinate
  double yc; // panel center y coordinate
  double zc; // panel center z coordinate
  vector<size_t> children;
  double MAC; // r^2 / theta^2
  double moments[12][Pflat];
  int moment_flag;
  double t1[P + 1]; // interpolation points in x direction
  double t2[P + 1];
  double t3[P + 1];
  panel() // initialization
  {
    moment_flag = 0;
    members[0] = 0;
    members[1] = -1;
    for (size_t kk = 0; kk < Pflat + 1; kk++) {
      for (int i = 0; i < 12; i++)
	moments[i][kk] = 0;
    }
    
    for (int i = 0; i < P + 1; i++)
      {
	t1[i] = 0.0;
	t2[i] = 0.0;
	t3[i] = 0.0;
      }
  }
};

//********************************************************************//

void eval_Stokes_SL_onSurf_tree(int surf_num, double h, double DEL,
				int N_quad, const vector<Surf_point>& Surfc,
				int N_target, const vector<Target_point>& Target,
				vector<double>& SL);

void eval_Stokes_SL_offSurf_tree(int surf_num, double h, double DEL,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL);


void eval_Stokes_DL_onSurf_tree(int surf_num, double h, double DEL,
				int N_quad, const vector<Surf_point>& Surfc,
				int N_target, const vector<Target_point>& Target,
				vector<double>& DL);

void eval_Stokes_DL_offSurf_tree(int surf_num, double h, double DEL,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL);
