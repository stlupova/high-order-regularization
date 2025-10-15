#include <iostream>
#include <vector>

#include "surface.h"
#include <cmath>

using namespace std;

//*********************************************************************
// Find surface (quadrature) points using the method of Beale et al 2017
// Uses the simple bracketing algorithm (Wilson's thesis, sec. 2.2.1)

void Generate_Surface(int surf_num, int rotated,
		      int N, double h, vector<Surf_point>* Surfc) {
  double H_seg = h; // segment spacing
  vector<double> A(3,0), B(3,0); // second argument puts 0 as every value
  int count = 0;
  
  for (int i = 0; i < 3; i++) {
    for (int l1 = -N; l1 <= N; l1++) {
      
      for (int l2 = -N; l2 <= N; l2++) {
	if (i==0) {
	  A[1] = l1*h; A[2] = l2*h;
	  B[1] = l1*h; B[2] = l2*h;
	}
	if (i==1) {
	  A[0] = l1*h; A[2] = l2*h;
	  B[0] = l1*h; B[2] = l2*h;
	}
	if (i==2) {
	  A[0] = l1*h; A[1] = l2*h;
	  B[0] = l1*h; B[1] = l2*h;
	}
	for (int j = -N; j <= N; j++) { //j-th segment: x_i in [jH,(j+1)H]
	  A[i] = j*H_seg;      //end points of segment
	  B[i] = A[i] + H_seg;

	  //check that S(lambda,i,j) is a bracket
	  if ( sign_phi(surf_num,rotated,A) * sign_phi(surf_num,rotated,B) < 0 ) { 

	    //find a hypersurface point in S
	    double tau = Find_hypersurf_pt(surf_num, rotated, A, H_seg, i); 
	    Surf_point temp;
	    temp.x = A;
	    temp.x[i] += tau * H_seg; //the point then has i-th coordinate = (j+tau)*H

	    //compute the normal vector
	    for (int k=0; k<3; k++) temp.Nrml[k] = Dphi(surf_num,rotated,k,temp.x); 
	    double Norm_Dphi = sqrt(temp.Nrml[0]*temp.Nrml[0]
				  + temp.Nrml[1]*temp.Nrml[1]
				  + temp.Nrml[2]*temp.Nrml[2]);
	    for (int k=0; k<3; k++) temp.Nrml[k] /= Norm_Dphi; //make it a unit normal

	    if ( abs(temp.Nrml[i]) >= cos(theta) ) { //see Beale et al (2.4)
	      temp.Area = h*h*Part_Unity(i,temp.Nrml)/abs(temp.Nrml[i]); //"area element"
	      temp.left  = -1;
	      temp.right = -1;
	      temp.up    = -1;
	      temp.down  = -1;
	      temp.plane = i;
	      
	      Surfc->push_back(temp); //add the point to the quadrature list
	      count = count + 1;
	    }
	  }                
	}
      }
    }
  }
}

//*********************************************************************
// Sign function as defined in Wilson's thesis p.12

int sign_phi(int surf_num, int rotated, const vector<double>& x) {
  if ( phi(surf_num,rotated,x) < 0.0 )
      return -1;
  else
    return 1;
}

//*********************************************************************
// Find a quadrature point by a line search algorithm, essentially
// looking for a root of phi on a given interval using Newton's/bisection method

double Find_hypersurf_pt(int surf_num, int rotated,
			 const vector<double>& A, double H_seg, int i) {
  //return bisection(A,H_seg,i,0.0,1.0);
  return Newton(surf_num, rotated, A, H_seg, i, 0, 1);
}

//*********************************************************************

double Newton(int surf_num, int rotated,
	      const vector<double>& PT, double H_seg,
	      int i, double a, double b) {

    int n_max = 100;
    vector<double> x_left = PT;
    vector<double> x_right = PT;
    x_left[i] += a*H_seg;
    x_right[i] += b*H_seg;
    vector<double> x_n = PT;
    double phi_n, Dphi_n;
    int n = 0;
    double p = 0.0;

    if ( phi(surf_num, rotated, x_left) * phi(surf_num, rotated, x_right) > 0.0 ) 
      cout << "Root finding will fail: endpoints do not have opposite sign" << endl;
    else if ( abs(phi(surf_num, rotated, x_left)) <= tol )
        p = a;
    else if ( abs(phi(surf_num, rotated, x_right)) <= tol )
        p = b;
    else {
      // start Newton's iteration using the midpoint as initial guess
      p = (a + b)/2;
      x_n[i] = PT[i] + p*H_seg;
      phi_n = phi(surf_num, rotated, x_n);
      Dphi_n = Dphi(surf_num, rotated, i, x_n)*H_seg;
      while ((abs(phi_n) > tol) && (abs(Dphi_n) > tol) && (n < n_max)) {
	p -= phi_n/Dphi_n; 
	x_n[i] = PT[i] + p*H_seg;
	n = n+1;
	phi_n = phi(surf_num, rotated, x_n);
	Dphi_n = Dphi(surf_num, rotated, i, x_n)*H_seg;
      }
      // if Newton's failed to converge, restart with the bisection method
      if ((n >= n_max) || (abs(Dphi_n) <= tol)) {
	cout << "Newton's didn't converge... restart with bisection..." << endl;
	p = bisection(surf_num, rotated, PT,H_seg,i,a,b);
      }
    }
    return p;
}

//*********************************************************************

double bisection(int surf_num, int rotated,
		 const vector<double>& PT, double H_seg,
		 int i, double a, double b) {

    vector<double> x_left = PT;
    vector<double> x_right = PT;
    x_left[i] += a*H_seg;
    x_right[i] += b*H_seg;
    vector<double> x_mid = PT;
    double phi_mid;
    double p = 0.0;

    if ( phi(surf_num, rotated, x_left) * phi(surf_num, rotated, x_right) > 0.0 ) 
      cout << "Bisection will fail: endpoints do not have opposite sign" << endl;
    else if ( abs(phi(surf_num, rotated, x_left)) <= tol )
      p = a;
    else if ( abs(phi(surf_num, rotated, x_right)) <= tol )
      p = b;
    else {
      p = (a + b)/2;
      x_mid[i] = PT[i] + p*H_seg;
      phi_mid = phi(surf_num, rotated, x_mid);
      while ( abs(phi_mid) > tol ) {
	if ( phi(surf_num, rotated, x_left) * phi(surf_num, rotated, x_mid) < 0.0 ) 
	  b = p;
	else {
	  a = p;     
	  x_left = x_mid;
	}
	p = (a + b)/2; 
	x_mid[i] = PT[i] + p*H_seg;
	phi_mid = phi(surf_num, rotated, x_mid);
      }
    }
    return p;
}

//*********************************************************************
// Computes the sigma functions defined in Beale et al p.5.
//          These are universal partitions of unity on the unit sphere,
//          but applied to the unit normal of the given surface

double Part_Unity(int i, const vector<double>& Nrml) {
  double n1 = abs( Nrml[0] ); // abs(n*e1)
  double n2 = abs( Nrml[1] ); // abs(n*e2)
  double n3 = abs( Nrml[2] ); // abs(n*e3)

  vector<double> var(3);
  var[0] = bump( acos(n1) / theta );
  var[1] = bump( acos(n2) / theta );
  var[2] = bump( acos(n3) / theta );
        
  double varsum = var[0] + var[1] + var[2];
  return var[i]/varsum;
}

//*********************************************************************
// Smooth bump function (Beale et al p.4)

double bump(double r) {
  
  if ( abs(r) < 1.0 ) {
    double r2 = r*r;
    return exp( 2.0 * r2 / (r2-1.0) );
  }
  else {
    return 0.0;
  }
}

//*********************************************************************
//*********************************************************************
// Compute mean curvature at a given point on the surface
// assume mean curvature is positive when normal is outward

double mean_curvature(int surf_num, int rotated, const vector<double>& x) {

  // special formula if the surface is a spheroid
  double a = ellipse_b;
  double c = ellipse_a;

  double a_sq = a * a;
  
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

  double cos_sq;
  
  if ( rotated == 1) {
    vector<double> x_unrt = unrotate(xx); //un-rotated (x,y,z)
    cos_sq = (x_unrt[1] * x_unrt[1] + x_unrt[2] * x_unrt[2]) / a_sq;
  }
  else {
    cos_sq = (xx[1] * xx[1] + xx[2] * xx[2]) / a_sq;
  }

  double temp = a_sq + (c*c - a_sq) * cos_sq;

  double res = c * (a_sq + temp) / ( 2.0 * a * sqrt(temp * temp * temp));

  return res;

  // general formula if the surface is not a spheroid
  /*
  double phi1,phi2,phi3,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33;
  
  phi1 = Dphi(surf_num,rotated,0,x);
  phi2 = Dphi(surf_num,rotated,1,x);
  phi3 = Dphi(surf_num,rotated,2,x);
  D2phi(surf_num,rotated,x,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33);

  double phi1_sq = phi1 * phi1;
  double phi2_sq = phi2 * phi2;
  double phi3_sq = phi3 * phi3;
  double res;
  res  = phi1_sq * phi22 - 2.0 * phi1 * phi2 * phi12 + phi2_sq * phi11;
  res += phi1_sq * phi33 - 2.0 * phi1 * phi3 * phi13 + phi3_sq * phi11;
  res += phi2_sq * phi33 - 2.0 * phi2 * phi3 * phi23 + phi3_sq * phi22;

  double Dphi_sq = phi1_sq + phi2_sq + phi3_sq;
  res = res / (2.0 * Dphi_sq * sqrt(Dphi_sq));

  return res;
  */
}

//*********************************************************************
// Find nearest point on the surface given nearest quadrature point

int find_Nearest_Start(int surf_num, int rotated, const vector<double>& pt,
		       int N_quad, const vector<Surf_point>& Surfc,
		       int nrst_i, Surf_point& nrst) {

  int n_max = 100;
  
  vector<double> x = Surfc[nrst_i].x; //starting point in Newton's iteration
  
  double lam = 0.0;
  int n = 0;
  
  // Function in Newton's method
  vector<double> F(4,0);
  F[0] = phi(surf_num, rotated, x); 
  F[1] = x[0] - pt[0];
  F[2] = x[1] - pt[1];
  F[3] = x[2] - pt[2];
  
  double Fnorm2 = F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
  
  while ((Fnorm2 > tol) && (n <= n_max)) {
    
    vector<double> line(4,0), soln(4,0);
    vector< vector<double> > Jacobian(4,line);
    vector<int> permute(4,0);
    double phi1,phi2,phi3,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33;
    
    phi1 = Dphi(surf_num, rotated, 0, x);
    phi2 = Dphi(surf_num, rotated, 1, x);
    phi3 = Dphi(surf_num, rotated, 2, x);
    D2phi(surf_num,rotated,x,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33);
    
    // Function in Newton's method
    F[0] = phi(surf_num, rotated, x);
    F[1] = x[0]-pt[0] + lam * phi1;
    F[2] = x[1]-pt[1] + lam * phi2;
    F[3] = x[2]-pt[2] + lam * phi3;
    
    Fnorm2 = F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
    
    // Jacobian in Newton's method
    Jacobian[0][0] = 0.0;
    Jacobian[0][1] = phi1;
    Jacobian[0][2] = phi2;
    Jacobian[0][3] = phi3;
    
    Jacobian[1][0] = phi1;
    Jacobian[1][1] = 1.0 + lam * phi11;
    Jacobian[1][2] = lam * phi12;
    Jacobian[1][3] = lam * phi13;
    
    Jacobian[2][0] = phi2;
    Jacobian[2][1] = lam * phi21;
    Jacobian[2][2] = 1.0 + lam * phi22;
    Jacobian[2][3] = lam * phi23;
    
    Jacobian[3][0] = phi3;
    Jacobian[3][1] = lam * phi31;
    Jacobian[3][2] = lam * phi32;
    Jacobian[3][3] = 1.0 + lam * phi33;
    
    int result = LU_factorization(Jacobian, permute);
    if (result==-1) {
      return (-1);
    }
    soln = solveByLU(Jacobian,permute,F);			
    
    //update the solution (lambda,x,y,z)
    lam = lam - soln[0];
    for (int i=0; i<3; i++) x[i] = x[i] - soln[i+1];
    n = n+1;
  }
  
  if (n>n_max) {
    return (-1);
  }
  
  nrst.x = x;
  double n1 = Dphi(surf_num, rotated, 0, x);
  double n2 = Dphi(surf_num, rotated, 1, x);
  double n3 = Dphi(surf_num, rotated, 2, x);
  double nnorm = sqrt(n1 * n1 + n2 * n2 + n3 * n3);
  
  nrst.Nrml[0] = n1 / nnorm;
  nrst.Nrml[1] = n2 / nnorm;
  nrst.Nrml[2] = n3 / nnorm;
  
  return 0;
}

//*********************************************************************
// Find nearest quadrature point to a target point

int find_Nearest_grid(const vector<double>& pt, int N_quad,
		      const vector<Surf_point>& Surfc, int i3) {

  int i0 = -1; // -1 returned if search unsuccessful
  double min_dist = 100.0;
  
  for (int i=0; i<N_quad; i++) {
    
    if (Surfc[i].plane != i3) continue; // only want the same Monge patch

    double dist = (pt[0]-Surfc[i].x[0])*(pt[0]-Surfc[i].x[0])
                + (pt[1]-Surfc[i].x[1])*(pt[1]-Surfc[i].x[1])
                + (pt[2]-Surfc[i].x[2])*(pt[2]-Surfc[i].x[2]);
    
    if (dist < min_dist) {
      min_dist = dist;
      i0 = i;
    }
  }
  return i0;
}

//*********************************************************************
// decide which Monge parameterization to use based on normal vector

void select_Monge_Patch(int surf_num, int rotated,
			const vector<double>& x0,
			int& i1, int& i2, int& i3) {

  double nx = abs( Dphi(surf_num, rotated, 0, x0) );
  double ny = abs( Dphi(surf_num, rotated, 1, x0) );
  double nz = abs( Dphi(surf_num, rotated, 2, x0) );

  i1 = 0;
  i2 = 0;
  i3 = 0;
  
  if ( nx >= ny ) {
    if ( nx >= nz ) { //normal is mostly in the x-direction
      i1 = 1;
      i2 = 2;
      i3 = 0;
    }
    else { //normal is mostly in the z-direction
      i1 = 0;
      i2 = 1;
      i3 = 2;
    }
  }
  else {
    if ( ny >= nz ) { //normal is mostly in the y-direction
      i1 = 0;
      i2 = 2;
      i3 = 1;
    }
    else { //normal is mostly in the z-direction
      i1 = 0;
      i2 = 1;
      i3 = 2;
    }
  }
}

//*********************************************************************
//*********************************************************************

void Generate_Targets_OnSurface(int N_quad, const vector<Surf_point>& Surfc,
				vector<Target_point>* Target) { 

  Target_point temp;

  for (int i=0; i<N_quad; i++) {
    temp.x = Surfc[i].x[0];
    temp.y = Surfc[i].x[1];
    temp.z = Surfc[i].x[2];
    Target->push_back(temp); //add the point to the target list
  }
}
//*********************************************************************

void Compute_nearest_data(int surf_num, int rotated, double h, double DEL,
			  int N_quad, const vector<Surf_point>& Surfc,
			  int N_target, vector<Target_point>& Target) { 

  
#pragma omp parallel
  {
  vector<double> pt(3,0);
    
#pragma omp for firstprivate(Surfc)
  for (int i=0; i<N_target; i++) {

    Target[i].nrst_f[0] = 0;
    Target[i].nrst_f[1] = 0;
    Target[i].nrst_f[2] = 0;

    pt[0] = Target[i].x;
    pt[1] = Target[i].y;
    pt[2] = Target[i].z;

    // decide which Monge parameterization to use
    // based on components of normal vector
    int i1, i2, i3;      
    select_Monge_Patch(surf_num, rotated, pt, i1, i2, i3);
    
    int nrst_i = find_Nearest_grid(pt, N_quad, Surfc, i3);
    if (nrst_i == -1 ) cout << "Warning: Nearest grid point not found!" << endl;
    Target[i].nrst = nrst_i;
    
    // find nearest point on the surface
    Surf_point nrst;
    int near = find_Nearest_Start(surf_num, rotated, pt, N_quad, Surfc, nrst_i, nrst);       
    if ( near == -1 ) cout << "Nearest point was not found for subtraction at target "
			   << i << endl;
    if ( near == -1 ) continue;
    
    Target[i].nrst_x = nrst.x[0];
    Target[i].nrst_y = nrst.x[1];
    Target[i].nrst_z = nrst.x[2];
    
    Target[i].nrst_Nrml[0] = nrst.Nrml[0];
    Target[i].nrst_Nrml[1] = nrst.Nrml[1];
    Target[i].nrst_Nrml[2] = nrst.Nrml[2];
    
    Target[i].i1 = i1;
    Target[i].i2 = i2;
    Target[i].i3 = i3;
    
    Target[i].b = distance(pt, nrst.x, nrst.Nrml); // signed distance from surface
    
    if ( abs(Target[i].b)/DEL > 8.0 )  {
      Target[i].nrst = 0;
    }
    if ( abs(Target[i].b)/DEL > 8.0 )  continue;

    // find the value of surface force at the nearest point
    Target[i].nrst_f = surfc_force(surf_num, rotated, nrst.x, nrst.Nrml);

    
    if (isnan(Target[i].nrst_f[0])) {
      cout << "found nan 0 for target point number " << i
	   << "  away from surfc " << surf_num << endl;
      cout << Target[i].x << "  " << Target[i].y << "  " << Target[i].z << endl;
    }
    if (isnan(Target[i].nrst_f[1])) {
      cout << "found nan 1 for target point number " << i
	   << "  away from surfc " << surf_num << endl;
      cout << Target[i].x << "  " << Target[i].y << "  " << Target[i].z << endl;
    }
    if (isnan(Target[i].nrst_f[2])) {
      cout << "found nan 2 for target point number " << i
	   << "  away from surfc " << surf_num << endl;
      cout << Target[i].x << "  " << Target[i].y << "  " << Target[i].z << endl;
    }
  }
 }
}

//*********************************************************************

vector<double> surfc_force(int surf_num, int rotated,
			   const vector<double>& x0,
			   const vector<double>& n0) { 

  double mean_curv = mean_curvature(surf_num,rotated,x0);

  double zz = x0[2];
  if (surf_num==1)  zz -= ctr1_z;
  if (surf_num==2)  zz -= ctr2_z;

  double gamma = 1.0 + zz * zz;

  double Dgam_dx = 0.0;
  double Dgam_dy = 0.0;
  double Dgam_dz = 2.0 * zz;

  double Dgam_dn = Dgam_dx * n0[0] + Dgam_dy * n0[1] + Dgam_dz * n0[2];
  
  double surf_grad0 = Dgam_dx - Dgam_dn * n0[0];
  double surf_grad1 = Dgam_dy - Dgam_dn * n0[1];
  double surf_grad2 = Dgam_dz - Dgam_dn * n0[2];
  
  vector<double> f(3,0);
  gamma *= 2.0 * mean_curv;
  f[0] = gamma * n0[0] - surf_grad0;
  f[1] = gamma * n0[1] - surf_grad1;
  f[2] = gamma * n0[2] - surf_grad2;
  
  return f;
}
