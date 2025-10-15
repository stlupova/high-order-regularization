#include <iostream>
#include <vector>
#include <cmath>

#include "utilities.h"

using namespace std;

void Generate_Surface(int surf_num, int rotated,
		      int N, double h, vector<Surf_point>* Surfc);
int sign_phi(int surf_num, int rotated, const vector<double>& x);
double Find_hypersurf_pt(int surf_num, int rotated,
			 const vector<double>& A, double H_seg, int i);

double Newton(int surf_num, int rotated,
	      const vector<double>& PT, double H_seg,
	      int i, double a, double b);
double bisection(int surf_num, int rotated,
		 const vector<double>& PT, double H_seg,
		 int i, double a, double b);

double Part_Unity(int i, const vector<double>& Nrml);
double bump(double r);

double mean_curvature(int surf_num, int rotated,
		      const vector<double>& x);

int find_Nearest_Start(int surf_num, int rotated, const vector<double>& pt,
		       int N_quad, const vector<Surf_point>& Surfc,
		       int nrst_i, Surf_point& nrst);
int find_Nearest_grid(const vector<double>& pt, int N_quad,
		      const vector<Surf_point>& Surfc, int i3);

void select_Monge_Patch(int surf_num, int rotated,
			const vector<double>& x0,
			int& i1, int& i2, int& i3);

void Generate_Targets_OnSurface(int N_quad, const vector<Surf_point>& Surfc,
				vector<Target_point>* Target);

void Compute_nearest_data(int surf_num, int rotated, double h, double DEL,
			  int N_quad, const vector<Surf_point>& Surfc,
			  int N_target, vector<Target_point>& Target);

vector<double> surfc_force(int surf_num, int rotated,
			   const vector<double>& x0,
			   const vector<double>& n0);
