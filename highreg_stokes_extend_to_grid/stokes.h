#include <iostream>
#include <vector>
#include <cmath>

#include "surface.h"

using namespace std;

void stokes_SL_away(int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& f,
		    double x, double y, double z,
		    vector<double>& SL);

void stokes_SL_7ord(double DEL,
		    int N_quad, const vector<Surf_point>& Surfc,
		    const vector<double>& f, double f0_dot_n0,
		    double x, double y, double z, double b,
		    vector<double>& SL);

double stokes_pressure_away(int N_quad, const vector<Surf_point>& Surfc,
			    const vector<double>& f,
			    double x, double y, double z);

double stokes_pressure_7ord(double DEL,
			    int N_quad, const vector<Surf_point>& Surfc,
			    const vector<double>& f, const vector<double>& f0,
			    const vector<double>& n0,
			    double x, double y, double z, double b);
