#include <iostream>
#include <vector>
#include <cmath>

#include "KITC.h"

using namespace std;

void eval_Stokes_SL_offSurf(double h, double DEL,
			    int N_quad, const vector<Surf_point>& Surfc,
			    int N_target, const vector<Target_point>& Target,
			    vector<double>& SL);

void stokes_SL_7ord(double DEL,
		    int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& f, double f0_dot_n0,
		    double x, double y, double z, double b,
		    vector<double>& SL);
  
void eval_Stokes_SL_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   int N_target, const vector<Target_point>& Target,
			   vector<double>& SL_vec);

void stokes_SL_7ord_onSurf(double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   const vector<double>& f, double f0_dot_n0,
			   double x, double y, double z,
			   vector<double>& SL);

