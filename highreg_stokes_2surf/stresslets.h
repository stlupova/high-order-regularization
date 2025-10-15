#include <iostream>
#include <vector>
#include <cmath>

#include "stokeslets.h"

using namespace std;

void eval_Stokes_DL_offSurf(double h, double DEL,
			    int N_quad, const vector<Surf_point>& Surfc,
			    int N_target, vector<Target_point>& Target,
			    vector<double>& DL_vec);

void stokes_DL_7ord(double DEL,
		    int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& g, const vector<double>& g0,
		    const vector<double>& x0, const vector<double>& n0,
		    double x, double y, double z, double b,
		    vector<double>& DL);

void eval_Stokes_DL_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   int N_target, vector<Target_point>& Target,
			   vector<double>& DL_vec);

void stokes_DL_7ord_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   const vector<double>& g, const vector<double>& g0,
			   double x, double y, double z,
			   vector<double>& DL);

