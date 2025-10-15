#include <iostream>
#include <vector>

#include "surface.h" 

using namespace std;

//*********************************************************************
// Finds surface (quadrature) points using the method of Beale et al 2016
//      Uses the simple bracketing algorithm (Wilson's thesis, sec. 2.2.1)

void Generate_Surface(int N, double h, vector<Surf_point>* Surfc) {
  double theta = 70.0*PI/180.0; //angle in Beale et al (2.2)
  double H_seg = h; // segment spacing
  vector<double> A(3,0), B(3,0); // second argument puts 0 as every value
  int count = 0;
  
  for (int i = 0; i < 3; i++) {
    for (int l1 = -N; l1 <= N; l1++) {
      for (int l2 = -N; l2 <= N; l2++) {
          int ind1,ind2;
	if (i==0) {
        ind1=1; ind2=2;
	}
	if (i==1) {
        ind1=0; ind2=2;
	}
	if (i==2) {
        ind1=0; ind2=1;
	}
	
	A[ind1] = l1*h + ctr;   A[ind2] = l2*h + ctr;
	B[ind1] = l1*h + ctr;   B[ind2] = l2*h + ctr;

	for (int j = -N; j <= N; j++) { //j-th segment: x_i in [jH,(j+1)H]
	  A[i] = j*H_seg + ctr;      //end points of segment
	  B[i] = A[i] + H_seg;

	  if ( sign_phi(A)*sign_phi(B) < 0 ) { //check that S(lambda,i,j) is a bracket
	    
	    double tau = Find_hypersurf_pt(A,H_seg,i); //find a hypersurface point in S
	    Surf_point temp;
	    temp.x = A;
	    temp.x[i] += tau * H_seg; //the point then has i-th coordinate = (j+tau)*H
	    
	    for (int k=0; k<3; k++) temp.Nrml[k] = Dphi(k,temp.x); //compute the normal vector

	    double Norm_Dphi = sqrt(temp.Nrml[0]*temp.Nrml[0]
				  + temp.Nrml[1]*temp.Nrml[1]
				  + temp.Nrml[2]*temp.Nrml[2]);
	    
	    //make it a unit normal
	    for (int k=0; k<3; k++) temp.Nrml[k] /= Norm_Dphi; 

	    if ( abs(temp.Nrml[i]) >= cos(theta) ) { //see Beale et al (2.4)
	      temp.Area = h * h * Part_Unity(theta, i, temp.Nrml)
		                / abs(temp.Nrml[i]); //"area element"	      
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

int sign_phi(const vector<double>& x) {
  if ( phi(x) < 0.0 )
      return -1;
  else
    return 1;
}

//*********************************************************************
// Finds a quadrature point by a line search algorithm, essentially
// looking for a root of phi on a given interval using Newton's/bisection method

double Find_hypersurf_pt(const vector<double>& A, double H_seg, int i) {
  return Newton(A,H_seg,i,0,1);
}

//*********************************************************************

double Newton(const vector<double>& PT, double H_seg,
	      int i, double a, double b) {

  double tol = 1.e-14;
  int n_max = 100;
  vector<double> x_left = PT;
  vector<double> x_right = PT;
  x_left[i] += a*H_seg;
  x_right[i] += b*H_seg;
  vector<double> x_n = PT;
  double phi_n, Dphi_n;
  int n = 0;
  double p = 0.0;
  
  if ( phi(x_left)*phi(x_right) > 0.0 ) 
    cout << "Root finding will fail: endpoints do not have opposite sign" << endl;
  else if ( abs(phi(x_left)) <= tol )
    p = a;
  else if ( abs(phi(x_right)) <= tol )
    p = b;
  else {
    // start Newton's iteration using the midpoint as initial guess
    p = (a + b)/2;
    x_n[i] = PT[i] + p*H_seg;
    phi_n = phi(x_n);
    Dphi_n = Dphi(i,x_n)*H_seg;
    while ((abs(phi_n) > tol) && (abs(Dphi_n) > tol) && (n < n_max)) {
      p -= phi_n/Dphi_n; 
      x_n[i] = PT[i] + p*H_seg;
      n = n+1;
      phi_n = phi(x_n);
      Dphi_n = Dphi(i,x_n)*H_seg;
    }
    // if Newton's failed to converge, restart with the bisection method
    if ((n >= n_max) || (abs(Dphi_n) <= tol)) {
      cout << "Newton's didn't converge... restart with bisection..." << endl;
      p = bisection(PT,H_seg,i,a,b);
    }
  }
  return p;
}

//*********************************************************************

double bisection(const vector<double>& PT, double H_seg,
		 int i, double a, double b) {

  double tol = 1.e-14;
  vector<double> x_left = PT;
  vector<double> x_right = PT;
  x_left[i] += a*H_seg;
  x_right[i] += b*H_seg;
  vector<double> x_mid = PT;
  double phi_mid;
  double p = 0.0;
  
  if ( phi(x_left)*phi(x_right) > 0.0 ) 
    cout << "Bisection will fail: endpoints do not have opposite sign" << endl;
  else if ( abs(phi(x_left)) <= tol )
    p = a;
  else if ( abs(phi(x_right)) <= tol )
    p = b;
  else {
    p = (a + b)/2;
    x_mid[i] = PT[i] + p*H_seg;
    phi_mid = phi(x_mid);
    while ( abs(phi_mid) > tol ) {
      if ( phi(x_left)*phi(x_mid) < 0.0 ) 
	b = p;
      else {
	a = p;     
	x_left = x_mid;
      }
      p = (a + b)/2; 
      x_mid[i] = PT[i] + p*H_seg;
      phi_mid = phi(x_mid);
    }
  }
  return p;
}

//*********************************************************************
// Computes the sigma functions defined in Beale et al p.5.
//          These are universal partitions of unity on the unit sphere,
//          but applied to the unit normal of the given surface

double Part_Unity(double theta, int i, const vector<double>& Nrml) {
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
    
    double b = 0.0;
    if ( abs(r) < 1.0 ) {
      double r2 = r*r;
      b = exp( 2.0 * r2 / (r2 - 1.0) );
    }
    return b;
}

//*********************************************************************
// Finds nearest surface point to a target point pt

int find_Nearest(const vector<double>& pt,
		  int N_quad, const vector<Surf_point>& Surfc,
		  const vector<double>& x0,
		  Surf_point& nrst) {
  
  double tol = 1.e-14;
  int n_max = 100;

  vector<double> x = x0; //starting point in Newton's iteration
    
  // Function in Newton's method
  vector<double> F(4,0);
  F[0] = phi(x); 
  F[1] = x[0]-pt[0];
  F[2] = x[1]-pt[1];
  F[3] = x[2]-pt[2];

  double Fnorm2 = F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
  
  int n = 0;
  double lam = 0.0;

  while ((Fnorm2 > tol) && (n <= n_max)) {

    vector<double> line(4,0), soln(4,0);
    vector< vector<double> > Jacobian(4,line);
    vector<int> permute(4,0);
    double phi1,phi2,phi3,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33;
    
    phi1 = Dphi(0,x);
    phi2 = Dphi(1,x);
    phi3 = Dphi(2,x);
    D2phi(x,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33);
    
    // Function in Newton's method
    F[0] = phi(x);
    F[1] = x[0]-pt[0] + lam*phi1;
    F[2] = x[1]-pt[1] + lam*phi2;
    F[3] = x[2]-pt[2] + lam*phi3;
    
    Fnorm2 = F[0]*F[0] + F[1]*F[1] + F[2]*F[2] + F[3]*F[3];
    
    // Jacobian in Newton's method
    Jacobian[0][0] = 0.0;
    Jacobian[0][1] = phi1;
    Jacobian[0][2] = phi2;
    Jacobian[0][3] = phi3;
    
    Jacobian[1][0] = phi1;
    Jacobian[1][1] = 1.0 + lam*phi11;
    Jacobian[1][2] = lam*phi12;
    Jacobian[1][3] = lam*phi13;
    
    Jacobian[2][0] = phi2;
    Jacobian[2][1] = lam*phi21;
    Jacobian[2][2] = 1.0 + lam*phi22;
    Jacobian[2][3] = lam*phi23;
    
    Jacobian[3][0] = phi3;
    Jacobian[3][1] = lam*phi31;
    Jacobian[3][2] = lam*phi32;
    Jacobian[3][3] = 1.0 + lam*phi33;

    int result = LU_factorization(Jacobian,permute);
    if (result==-1) {     
      return (-1);
    }
    soln = solveByLU(Jacobian,permute,F);			
    
    //update the solution (lambda,x,y,z)
    lam = lam - soln[0];
    for (int i=0; i<3; i++) x[i] = x[i]-soln[i+1];
    n = n+1;

  }

  if (n>n_max) {
    return (-1);
  }

  nrst.x = x;
  double n1 = Dphi(0,x);
  double n2 = Dphi(1,x);
  double n3 = Dphi(2,x);
  double nnorm = sqrt(n1*n1 + n2*n2 + n3*n3);

  nrst.Nrml[0] = n1/nnorm;
  nrst.Nrml[1] = n2/nnorm;
  nrst.Nrml[2] = n3/nnorm;
  
  return 0;
}

//*********************************************************************
// find and tag grid points that are within BAND1*h of surface
// also tag them if they are within BAND2*h of surface (BAND2 < BAND1)
// for each of these points, find the nearest point on the surface
// find the signed distance to the surface b

void find_nrst_surfc(int N, double h,
		     double BAND1, double BAND2,
		     int N_quad, const vector<Surf_point>& Surfc,
		     int*** nrst, vector<Target_point>* Target) {

  double reject = 220.0 * h * h;

  for (int i=0; i<=N; i++) {
    for (int j=0; j<=N; j++) {
      for (int k=0; k<=N; k++) {
	nrst[i][j][k] = -1;
      }
    }
  }

  double tol = 1.e-12;
  vector<double> x_surf(3,0), pt(3,0);
  Target_point temp;
  int count = 0;
  
  int width = ceil(BAND1);
  
  for (int surf_ind=0; surf_ind<N_quad; surf_ind++) {

    x_surf[0] = Surfc[surf_ind].x[0];
    x_surf[1] = Surfc[surf_ind].x[1];
    x_surf[2] = Surfc[surf_ind].x[2];
    
    int i0 = floor(x_surf[0] / h - BAND1) - 1;
    int j0 = floor(x_surf[1] / h - BAND1) - 1;
    int k0 = floor(x_surf[2] / h - BAND1) - 1;

    for (int i=i0; i <= i0 + 2*width + 3; i++) {
      for (int j=j0; j <= j0 + 2*width + 3; j++) {
	for (int k=k0; k <= k0 + 2*width + 3; k++) {

	  pt[0] = i * h;
	  pt[1] = j * h;
	  pt[2] = k * h;

	  double d_sq = (pt[0] - x_surf[0]) * (pt[0] - x_surf[0])
	              + (pt[1] - x_surf[1]) * (pt[1] - x_surf[1])
	              + (pt[2] - x_surf[2]) * (pt[2] - x_surf[2]);
	  
	  if ( sqrt(d_sq) <= ((BAND1+1) * h + tol) ) {

	    // check that the point has not been tagged earlier
	    if ( nrst[i][j][k] == -1 ) {

	      // find the nearest point on surface
	      Surf_point pt_surf;
	      int near = find_Nearest(pt, 
				      N_quad, Surfc,
				      x_surf,
				      pt_surf);
	      if ( near == 0 ) {

		// b = signed distance from surface
		double b = distance(pt, pt_surf.x, pt_surf.Nrml); 

		if ( abs(b) <= (BAND1 * h + tol) ) {

		  temp.nrst_x = pt_surf.x[0]; // save nearest xyz
		  temp.nrst_y = pt_surf.x[1];
		  temp.nrst_z = pt_surf.x[2];
		  temp.nrst_Nrml[0] = pt_surf.Nrml[0];
		  temp.nrst_Nrml[1] = pt_surf.Nrml[1];
		  temp.nrst_Nrml[2] = pt_surf.Nrml[2];
		  temp.b = b;
		  
		  temp.band2 = false;
		  if ( abs(b) <= BAND2 * h + tol ) {
		    temp.band2 = true; // point is in band2
		  }
		  Target->push_back(temp); //add the point to the target list

		  nrst[i][j][k] = count; // save which target point it is for band1
		  count++;
		} // end tagging point
	      }
	      else {
		cout << "ERROR: Finding nearest point on surface failed for a point that needs to be tagged!" << endl;
	      }
	    } // end if point was not tagged earlier
	  } // end if grid point within band1 of a quadrature point
	} // end for k
      } // end for j
    } // end for i
  } // end for loop over surface points
}
