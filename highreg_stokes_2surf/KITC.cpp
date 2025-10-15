
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/times.h>

#include "KITC.h"

using namespace std;


//**********************************************************//
vector<panel> tree;
vector<size_t> leaf;
int max_level = 0;


static size_t node_count = 0;

//**********************************************************************//
long getTickCount() {
  tms tm;
  return times(&tm);
}

//**********************************************************************//
void build_tree_init(int surf_num, int N_cube) {
  panel temp_panel;
  
  // indices of particles belonging to panel
  temp_panel.members[0] = 0;
  temp_panel.members[1] = N_cube - 1;

  double eps = 1.e-8;
  temp_panel.xinterval[0] = -1.0-eps; // interval defining the panel
  temp_panel.xinterval[1] = 1.0+eps;
  temp_panel.yinterval[0] = -1.0-eps;
  temp_panel.yinterval[1] = 1.0+eps;
  temp_panel.zinterval[0] = -1.0-eps;
  temp_panel.zinterval[1] = 1.0+eps;
  
  if (surf_num == 1) {
    temp_panel.xinterval[0] += ctr1_x;
    temp_panel.xinterval[1] += ctr1_x;
    temp_panel.yinterval[0] += ctr1_y;
    temp_panel.yinterval[1] += ctr1_y;
    temp_panel.zinterval[0] += ctr1_z;
    temp_panel.zinterval[1] += ctr1_z;
  }
  else {
    temp_panel.xinterval[0] += ctr2_x;
    temp_panel.xinterval[1] += ctr2_x;
    temp_panel.yinterval[0] += ctr2_y;
    temp_panel.yinterval[1] += ctr2_y;
    temp_panel.zinterval[0] += ctr2_z;
    temp_panel.zinterval[1] += ctr2_z;
  }
  
  temp_panel.xc = 0.5 * (temp_panel.xinterval[0] + temp_panel.xinterval[1]);
  temp_panel.yc = 0.5 * (temp_panel.yinterval[0] + temp_panel.yinterval[1]);
  temp_panel.zc = 0.5 * (temp_panel.zinterval[0] + temp_panel.zinterval[1]);
  temp_panel.MAC = (3 * 10 * 10/ 4) / sq_theta; // MAC = r^2 / theta^2
  
  tree.push_back(temp_panel);
  node_count = 1;
}

//**********************************************************************//
void Swap(size_t i, size_t j, double *lambda[3], struct xyz &s) {
  if (i == j)
    return;
  
  double x = s.x[i];
  double y = s.y[i];
  double z = s.z[i];
  size_t index = s.index[i];
  size_t old_index = s.old_index[i];
  double nx = s.nx[i];
  double ny = s.ny[i];
  double nz = s.nz[i];
  double area = s.area[i];
  double lam0 = lambda[0][i];
  double lam1 = lambda[1][i];
  double lam2 = lambda[2][i];
  
  s.x[i] = s.x[j];
  s.y[i] = s.y[j];
  s.z[i] = s.z[j];
  s.index[i] = s.index[j];
  s.old_index[i] = s.old_index[j];
  s.nx[i] = s.nx[j];
  s.ny[i] = s.ny[j];
  s.nz[i] = s.nz[j];
  s.area[i] = s.area[j];
  lambda[0][i] = lambda[0][j];
  lambda[1][i] = lambda[1][j];
  lambda[2][i] = lambda[2][j];
  
  s.x[j] = x;
  s.y[j] = y;
  s.z[j] = z;
  s.index[j] = index;
  s.old_index[j] = old_index;
  s.nx[j] = nx;
  s.ny[j] = ny;
  s.nz[j] = nz;
  s.area[j] = area;
  lambda[0][j] = lam0;
  lambda[1][j] = lam1;
  lambda[2][j] = lam2;
}
//***********************************************************************//
void split_tree_node(size_t panel_index,
		     double *lambda[3],
		     struct xyz &particles) {

  panel child[8];
  
  double tp_x0 = tree[panel_index].xinterval[0];
  double tp_x1 = tree[panel_index].xinterval[1];
  double tp_y0 = tree[panel_index].yinterval[0];
  double tp_y1 = tree[panel_index].yinterval[1];
  double tp_z0 = tree[panel_index].zinterval[0];
  double tp_z1 = tree[panel_index].zinterval[1];
  
  double midpointx = (tp_x0 + tp_x1) / 2.0;
  double midpointy = (tp_y0 + tp_y1) / 2.0;
  double midpointz = (tp_z0 + tp_z1) / 2.0;
  
  double xc0 = (tp_x0 + midpointx) / 2.0;
  double xc1 = (tp_x1 + midpointx) / 2.0;
  double yc0 = (tp_y0 + midpointy) / 2.0;
  double yc1 = (tp_y1 + midpointy) / 2.0;
  double zc0 = (tp_z0 + midpointz) / 2.0;
  double zc1 = (tp_z1 + midpointz) / 2.0;
  
  child[0].xinterval[0] = tp_x0;
  child[0].xinterval[1] = midpointx;
  child[0].yinterval[0] = tp_y0;
  child[0].yinterval[1] = midpointy;
  child[0].zinterval[0] = tp_z0;
  child[0].zinterval[1] = midpointz;
  child[0].xc = xc0;
  child[0].yc = yc0;
  child[0].zc = zc0;
  child[0].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[1].xinterval[0] = midpointx;
  child[1].xinterval[1] = tp_x1;
  child[1].yinterval[0] = tp_y0;
  child[1].yinterval[1] = midpointy;
  child[1].zinterval[0] = tp_z0;
  child[1].zinterval[1] = midpointz;
  child[1].xc = xc1;
  child[1].yc = yc0;
  child[1].zc = zc0;
  child[1].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[2].xinterval[0] = tp_x0;
  child[2].xinterval[1] = midpointx;
  child[2].yinterval[0] = midpointy;
  child[2].yinterval[1] = tp_y1;
  child[2].zinterval[0] = tp_z0;
  child[2].zinterval[1] = midpointz;
  child[2].xc = xc0;
  child[2].yc = yc1;
  child[2].zc = zc0;
  child[2].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[3].xinterval[0] = midpointx;
  child[3].xinterval[1] = tp_x1;
  child[3].yinterval[0] = midpointy;
  child[3].yinterval[1] = tp_y1;
  child[3].zinterval[0] = tp_z0;
  child[3].zinterval[1] = midpointz;
  child[3].xc = xc1;
  child[3].yc = yc1;
  child[3].zc = zc0;
  child[3].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (midpointz - zc0) * (midpointz - zc0)) / sq_theta;
  
  child[4].xinterval[0] = tp_x0;
  child[4].xinterval[1] = midpointx;
  child[4].yinterval[0] = tp_y0;
  child[4].yinterval[1] = midpointy;
  child[4].zinterval[0] = midpointz;
  child[4].zinterval[1] = tp_z1;
  child[4].xc = xc0;
  child[4].yc = yc0;
  child[4].zc = zc1;
  child[4].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  child[5].xinterval[0] = midpointx;
  child[5].xinterval[1] = tp_x1;
  child[5].yinterval[0] = tp_y0;
  child[5].yinterval[1] = midpointy;
  child[5].zinterval[0] = midpointz;
  child[5].zinterval[1] = tp_z1;
  child[5].xc = xc1;
  child[5].yc = yc0;
  child[5].zc = zc1;
  child[5].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (midpointy - yc0) * (midpointy - yc0)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  child[6].xinterval[0] = tp_x0;
  child[6].xinterval[1] = midpointx;
  child[6].yinterval[0] = midpointy;
  child[6].yinterval[1] = tp_y1;
  child[6].zinterval[0] = midpointz;
  child[6].zinterval[1] = tp_z1;
  child[6].xc = xc0;
  child[6].yc = yc1;
  child[6].zc = zc1;
  child[6].MAC = ((midpointx - xc0) * (midpointx - xc0)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  child[7].xinterval[0] = midpointx;
  child[7].xinterval[1] = tp_x1;
  child[7].yinterval[0] = midpointy;
  child[7].yinterval[1] = tp_y1;
  child[7].zinterval[0] = midpointz;
  child[7].zinterval[1] = tp_z1;
  child[7].xc = xc1;
  child[7].yc = yc1;
  child[7].zc = zc1;
  child[7].MAC = ((tp_x1 - xc1) * (tp_x1 - xc1)
		  + (tp_y1 - yc1) * (tp_y1 - yc1)
		  + (tp_z1 - zc1) * (tp_z1 - zc1)) / sq_theta;
  
  vector<size_t> v[8];
  size_t start = tree[panel_index].members[0];
  size_t end = tree[panel_index].members[1];
  size_t* addr_table = new size_t[end - start + 1];
  
  size_t index;
  for (index = start; index <= end; index++)
    {
      particles.index[index] = index;
      addr_table[index - start] = index;
      
      if (particles.x[index] <= midpointx
	  && particles.y[index] <= midpointy &&
	  particles.z[index] <= midpointz)
	v[0].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] <= midpointy &&
	       particles.z[index] <= midpointz )
	v[1].push_back(index);
      else if (particles.x[index] <= midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index]<= midpointz)
	v[2].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index] <= midpointz)
	v[3].push_back(index);
      else if(particles.x[index] <= midpointx
	      && particles.y[index] <= midpointy &&
	      particles.z[index] > midpointz )
	v[4].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] <= midpointy &&
	       particles.z[index] > midpointz)
	v[5].push_back(index);
      else if (particles.x[index] <= midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index] > midpointz)
	v[6].push_back(index);
      else if (particles.x[index] > midpointx
	       && particles.y[index] > midpointy &&
	       particles.z[index] > midpointz)
	v[7].push_back(index);
    }
  
  size_t seq = start;
  for (size_t j = 0; j < 8; j++)
    {
      size_t size = v[j].size();
      
      if (size >= 1)
	{
	  for (size_t k = 0; k < size; k++)
	    {
	      if (k == 0)
		child[j].members[0] = seq;
	      if (k == size - 1)
		child[j].members[1] = seq;
	      
	      index = v[j][k];
	      // This uses an address table
	      size_t pos = addr_table[index - start];
	      size_t out = particles.index[seq];
	      Swap(pos, seq, lambda, particles);
	      addr_table[index - start] = seq;
	      addr_table[out - start] = pos;
	      
	      seq++;
	    }
	  
	  node_count++;
	  tree[panel_index].children.push_back(node_count - 1);
	  tree.push_back(child[j]);
	  v[j].clear();
	}
    }
  
  delete[] addr_table;
}

//***********************************************************************//
void build_tree_3D_Recursive(size_t panel_index,
			     double *lambda[3],
			     struct xyz &particles,
			     int level) {
  if (level > max_level)
    max_level = level;
  
  size_t n = tree[panel_index].members[1]
    - tree[panel_index].members[0] + 1;
  
  if (n >= (size_t)N0) {
    split_tree_node(panel_index, lambda, particles);
    
    for (size_t i = 0; i < tree[panel_index].children.size(); i++) {
      size_t panel_index_new = tree[panel_index].children[i];
      build_tree_3D_Recursive(panel_index_new,
			      lambda,
			      particles,
			      level + 1);
    }
  }
  else
    leaf.push_back(panel_index);
}

//***********************************************************************//
void Panel_Moments_SL(size_t panel_index,
		      double *lambda[3],
		      struct xyz &particles,
		      double m[][Pflat]) {
  double t1[P + 1];
  double t2[P + 1];
  double t3[P + 1];
  
  for (int i = 0; i < P + 1; i++) {
    t1[i] = tree[panel_index].t1[i];
    t2[i] = tree[panel_index].t2[i];
    t3[i] = tree[panel_index].t3[i];
  }
  
  double w1i[P + 1];
  double w2j[P + 1];
  double w3k[P + 1];
  double dj[P + 1];
  dj[0] = 0.5;
  dj[P] = 0.5;
  for (int j = 1; j<P; j++)
    dj[j] = 1;
  
  for (int j = 0; j < P + 1; j++)
    w3k[j] = w2j[j] = w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
  
  int tp0 = tree[panel_index].members[0];
  int tp1 = tree[panel_index].members[1];
  int cluster_size;
  cluster_size = tp1 - tp0 + 1;
  
  double **a1i;
  double **a2j;
  double **a3k;
  double *D;
  
  a1i = (double**)calloc(P + 1, sizeof(double*));
  a2j = (double**)calloc(P + 1, sizeof(double*));
  a3k = (double**)calloc(P + 1, sizeof(double*));
  D = (double*)calloc(cluster_size + 1, sizeof(double));
  for (int i = 0; i < (P + 1); i++) {
    a1i[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    a2j[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    a3k[i] = (double*)calloc(cluster_size + 1, sizeof(double));
  }
  

  for (size_t tp_j = tp0; tp_j <= tp1 ; tp_j++) {
    
    double x = particles.x[tp_j];
    double y = particles.y[tp_j];
    double z = particles.z[tp_j];

    int flag1 = -1;
    int flag2 = -1;
    int flag3 = -1;
    double SumA1 = 0.0;
    double SumA2 = 0.0;
    double SumA3 = 0.0;
    
    for (int j = 0; j < P + 1; j++) {
      if ( abs(x-t1[j])<1.e-07) { //removable singularity
	flag1 = j;
      }
      else {
	double temp1 = w1i[j] / (x - t1[j]);
	a1i[j][tp_j - tp0] = temp1;
	SumA1 += temp1;
      }
	
      if ( abs(y-t2[j])<1.e-07) { //removable singularity
	flag2 = j;
      }
      else {
	double temp2 = w2j[j] / (y - t2[j]);
	a2j[j][tp_j - tp0] = temp2;
	SumA2 += temp2;
      }
      
      if ( abs(z-t3[j])<1.e-07) { //removable singlarity
	flag3 = j;
      }
      else {
	double temp3 = w3k[j] / (z - t3[j]);
	a3k[j][tp_j - tp0] = temp3;
	SumA3 += temp3;
      }
      
    }

    // if a flag was set, adjust sum and a to handle removable singularity
    if (flag1 > -1) {
      SumA1 = 1.0;
      for (int j = 0; j < P + 1; j++) a1i[j][tp_j - tp0] = 0.0;
      a1i[flag1][tp_j - tp0] = 1.0;
    }
    if (flag2 > -1) {
      SumA2 = 1.0;
      for (int j = 0; j < P + 1; j++) a2j[j][tp_j - tp0] = 0.0;
      a2j[flag2][tp_j - tp0] = 1.0;
    }
    if (flag3 > -1) {
      SumA3 = 1.0;
      for (int j = 0; j < P + 1; j++) a3k[j][tp_j - tp0] = 0.0;
      a3k[flag3][tp_j - tp0] = 1.0;
    }
    
    D[tp_j - tp0] = 1.0 / (SumA1 * SumA2 * SumA3);
  }
  
  int kk = -1;
  for (int i = 0; i < P + 1; i++) {
    for (int j = 0; j < P + 1; j++) {
      for (int k = 0; k < P + 1; k++) {	
	kk = kk + 1;
	
	double sum[6] = {0};
	for (size_t tp_j = tp0; tp_j <= tp1; tp_j++) {
	  
	  double nx = particles.nx[tp_j];
	  double ny = particles.ny[tp_j];
	  double nz = particles.nz[tp_j];
	  
	  double area = particles.area[tp_j];
    
	  double s = a1i[i][tp_j - tp0] * a2j[j][tp_j - tp0]
	    * a3k[k][tp_j - tp0] * D[tp_j - tp0];
	  
	  s *= area;
	  
	  sum[0] += s * lambda[0][tp_j];
	  sum[1] += s * lambda[1][tp_j];
	  sum[2] += s * lambda[2][tp_j];
	  
	  sum[3] += s * nx;
	  sum[4] += s * ny;
	  sum[5] += s * nz;
	}
	for (int jj=0; jj<6; jj++)  m[jj][kk] = sum[jj];
      }
    }
  }
  
  for (int i = 0; i < P + 1; i++) {
    free(a1i[i]);
    free(a2j[i]);
    free(a3k[i]);
  }
  
  free(a1i);
  free(a2j);
  free(a3k);
  free(D); 
}

//***********************************************************************//
void Panel_Moments_DL(size_t panel_index,
		      double *lambda[3],
		      struct xyz &particles,
		      double m[][Pflat]) {
  double t1[P + 1];
  double t2[P + 1];
  double t3[P + 1];
  
  for (int i = 0; i < P + 1; i++) {
    t1[i] = tree[panel_index].t1[i];
    t2[i] = tree[panel_index].t2[i];
    t3[i] = tree[panel_index].t3[i];
  }
  
  double w1i[P + 1];
  double w2j[P + 1];
  double w3k[P + 1];
  double dj[P + 1];
  dj[0] = 0.5;
  dj[P] = 0.5;
  for (int j = 1; j<P; j++)
    dj[j] = 1;
  
  for (int j = 0; j < P + 1; j++)
    w3k[j] = w2j[j] = w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
  
  int tp0 = tree[panel_index].members[0];
  int tp1 = tree[panel_index].members[1];
  int cluster_size;
  cluster_size = tp1 - tp0 + 1;
  
  double **a1i;
  double **a2j;
  double **a3k;
  double *D;
  
  a1i = (double**)calloc(P + 1, sizeof(double*));
  a2j = (double**)calloc(P + 1, sizeof(double*));
  a3k = (double**)calloc(P + 1, sizeof(double*));
  D = (double*)calloc(cluster_size + 1, sizeof(double));
  for (int i = 0; i < (P + 1); i++) {
    a1i[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    a2j[i] = (double*)calloc(cluster_size + 1, sizeof(double));
    a3k[i] = (double*)calloc(cluster_size + 1, sizeof(double));
  }
  

  for (size_t tp_j = tp0; tp_j <= tp1 ; tp_j++) {
    
    double x = particles.x[tp_j];
    double y = particles.y[tp_j];
    double z = particles.z[tp_j];

    int flag1 = -1;
    int flag2 = -1;
    int flag3 = -1;
    double SumA1 = 0.0;
    double SumA2 = 0.0;
    double SumA3 = 0.0;
    
    for (int j = 0; j < P + 1; j++) {
      if ( abs(x-t1[j])<1.e-07) { //removable singularity
	flag1 = j;
      }
      else {
	double temp1 = w1i[j] / (x - t1[j]);
	a1i[j][tp_j - tp0] = temp1;
	SumA1 += temp1;
      }
	
      if ( abs(y-t2[j])<1.e-07) { //removable singularity
	flag2 = j;
      }
      else {
	double temp2 = w2j[j] / (y - t2[j]);
	a2j[j][tp_j - tp0] = temp2;
	SumA2 += temp2;
      }
      
      if ( abs(z-t3[j])<1.e-07) { //removable singlarity
	flag3 = j;
      }
      else {
	double temp3 = w3k[j] / (z - t3[j]);
	a3k[j][tp_j - tp0] = temp3;
	SumA3 += temp3;
      }
      
    }

    // if a flag was set, adjust sum and a to handle removable singularity
    if (flag1 > -1) {
      SumA1 = 1.0;
      for (int j = 0; j < P + 1; j++) a1i[j][tp_j - tp0] = 0.0;
      a1i[flag1][tp_j - tp0] = 1.0;
    }
    if (flag2 > -1) {
      SumA2 = 1.0;
      for (int j = 0; j < P + 1; j++) a2j[j][tp_j - tp0] = 0.0;
      a2j[flag2][tp_j - tp0] = 1.0;
    }
    if (flag3 > -1) {
      SumA3 = 1.0;
      for (int j = 0; j < P + 1; j++) a3k[j][tp_j - tp0] = 0.0;
      a3k[flag3][tp_j - tp0] = 1.0;
    }
    
    D[tp_j - tp0] = 1.0 / (SumA1 * SumA2 * SumA3);
  }
  
  int kk = -1;
  for (int i = 0; i < P + 1; i++) {
    for (int j = 0; j < P + 1; j++) {
      for (int k = 0; k < P + 1; k++) {	
	kk = kk + 1;
	
	double sum[12] = {0};
	for (size_t tp_j = tp0; tp_j <= tp1; tp_j++) {
	  
	  double nx = particles.nx[tp_j];
	  double ny = particles.ny[tp_j];
	  double nz = particles.nz[tp_j];
	  
	  double area = particles.area[tp_j];
    
	  double s = a1i[i][tp_j - tp0] * a2j[j][tp_j - tp0]
	    * a3k[k][tp_j - tp0] * D[tp_j - tp0];
	  
	  s *= area;
	  
	  sum[0] += s * lambda[0][tp_j] * nx;
	  sum[1] += s * lambda[0][tp_j] * ny;
	  sum[2] += s * lambda[0][tp_j] * nz;
	  
	  sum[3] += s * lambda[1][tp_j] * nx;
	  sum[4] += s * lambda[1][tp_j] * ny;
	  sum[5] += s * lambda[1][tp_j] * nz;
	  
	  sum[6] += s * lambda[2][tp_j] * nx;
	  sum[7] += s * lambda[2][tp_j] * ny;
	  sum[8] += s * lambda[2][tp_j] * nz;
	  
	  sum[9]  += s * nx;
	  sum[10] += s * ny;
	  sum[11] += s * nz;
	}
	for (int jj=0; jj<12; jj++)  m[jj][kk] = sum[jj];
      }
    }
  }
  
  for (int i = 0; i < P + 1; i++) {
    free(a1i[i]);
    free(a2j[i]);
    free(a3k[i]);
  }
  
  free(a1i);
  free(a2j);
  free(a3k);
  free(D); 
}

//**********************************************************************//
vector<double> Call_Treecode_SL(double x, double y, double z,
				double f0_dot_n0,
				int panel_index,
				double DEL) {
  
  vector<double> u(3,0);
    
  double dx[P + 1];
  double dy[P + 1];
  double dz[P + 1];
  
  for (int i = 0; i < P + 1; i++) {
    dx[i] = x - tree[panel_index].t1[i];
    dy[i] = y - tree[panel_index].t2[i];
    dz[i] = z - tree[panel_index].t3[i];
  }
  
  int kk = -1;
  for (int i = 0; i < P + 1; i++) {
    double dx0 = dx[i];
    double temp_i = dx0 * dx0;
    for (int j = 0; j < P + 1; j++) {
      double dx1 = dy[j];
      double temp_j =  dx1 * dx1;
      for (int k = 0; k < P + 1; k++) {
	double dx2 = dz[k];
	double temp_k = dx2 * dx2;
	kk = kk + 1;
	
	double temp_mom0 = tree[panel_index].moments[0][kk]
	  - f0_dot_n0 * tree[panel_index].moments[3][kk];
	double temp_mom1 = tree[panel_index].moments[1][kk]
	  - f0_dot_n0 * tree[panel_index].moments[4][kk];
	double temp_mom2 = tree[panel_index].moments[2][kk]
	  - f0_dot_n0 * tree[panel_index].moments[5][kk];

	double r2 = temp_i + temp_j + temp_k;
	double r = sqrt(r2);

	double temp = (temp_mom0 * dx0 + temp_mom1 * dx1 + temp_mom2 * dx2) / (r2*r);
	
	u[0] += temp_mom0/r + temp * dx0;
	u[1] += temp_mom1/r + temp * dx1;
	u[2] += temp_mom2/r + temp * dx2;

      }
    }
  }
  
  return u;
}

//**********************************************************************//
vector<double> Call_Treecode_DL(double x, double y, double z,
				const vector<double>& g0,
				int panel_index,
				double DEL) {

  vector<double> temp(3,0), u(3,0);
    
  double dx[P + 1];
  double dy[P + 1];
  double dz[P + 1];
  
  for (int i = 0; i < P + 1; i++) {
    dx[i] = x - tree[panel_index].t1[i];
    dy[i] = y - tree[panel_index].t2[i];
    dz[i] = z - tree[panel_index].t3[i];
  }
  
  int kk = -1;
  for (int i = 0; i < P + 1; i++) {
    double dx0 = dx[i];
    double temp_i = dx0 * dx0;
    for (int j = 0; j < P + 1; j++) {
      double dx1 = dy[j];
      double temp_j =  dx1 * dx1;
      for (int k = 0; k < P + 1; k++) {
	double dx2 = dz[k];
	double temp_k = dx2 * dx2;
	kk = kk + 1;
	
	for (int l = 0; l < 3; l++) {
	  double temp_mom1   = tree[panel_index].moments[3*l][kk]
	    - g0[l] * tree[panel_index].moments[9][kk];
	  double temp_mom2 = tree[panel_index].moments[3*l+1][kk]
	    - g0[l] * tree[panel_index].moments[10][kk];
	  double temp_mom3 = tree[panel_index].moments[3*l+2][kk]
	    - g0[l] * tree[panel_index].moments[11][kk];

	  temp[l] = temp_mom1 * dx0 + temp_mom2 * dx1 + temp_mom3 * dx2; 
	}

	double r2 = temp_i + temp_j + temp_k;
	double r = sqrt(r2);
	
	double kernel = (temp[0]*dx0 + temp[1]*dx1 + temp[2]*dx2) /(r2*r2*r);	

	u[0] += dx0 * kernel;
	u[1] += dx1 * kernel;
	u[2] += dx2 * kernel;
      }
    }
  }
  
  return u;
}

//**********************************************************************//
void Cluster_Chev_Points(size_t tree_size) {
  double h;
  h = 3.14159265358979323846/P;
  double t[P + 1] = {0.0};
  for (int i = 0; i < P + 1; i++)
    t[i] = cos(i * h);  //Chebyshev interpolation points [-1,1]
  
  double x1,x2,y1,y2,z1,z2;
  size_t tree_index;
  
  for (tree_index = 0; tree_index < tree_size ; tree_index++) {
    x1 = tree[tree_index].xinterval[0];
    x2 = tree[tree_index].xinterval[1];
    y1 = tree[tree_index].yinterval[0];
    y2 = tree[tree_index].yinterval[1];
    z1 = tree[tree_index].zinterval[0];
    z2 = tree[tree_index].zinterval[1];
    
    for (int i = 0; i < P + 1; i++) {// map to the cluster
      tree[tree_index].t1[i] =  x1 + (t[i] + 1)/2 * (x2 - x1);
      tree[tree_index].t2[i] =  y1 + (t[i] + 1)/2 * (y2 - y1);
      tree[tree_index].t3[i] =  z1 + (t[i] + 1)/2 * (z2 - z1);
    }
  }
}

//***********************************************************************//
vector<double> Call_Ds_SL_onSurf(int limit_1, int limit_2,
				 double p_x, double p_y, double p_z,
				 double f0_dot_n0,
				 struct xyz &particles,
				 double *lambda[3],
				 double DEL) {    

  vector<double> u(3,0);

  for (size_t j = limit_1; j <= limit_2; j++) {
    
    double dx = p_x - particles.x[j];
    double dy = p_y - particles.y[j];
    double dz = p_z - particles.z[j];
    
    double r = sqrt(dx * dx + dy * dy + dz * dz);
    
    double area = particles.area[j];

    double f_x = lambda[0][j] - f0_dot_n0 * particles.nx[j]; //use subtraction in Stokeslet
    double f_y = lambda[1][j] - f0_dot_n0 * particles.ny[j];
    double f_z = lambda[2][j] - f0_dot_n0 * particles.nz[j];
    
    double f_dot_dx = f_x * dx + f_y * dy + f_z * dz;
    
    double H1, H2;
    Stokes_SL_7ord_onSurf(r, DEL, H1, H2);
    
    H1 *= area;
    H2 *= area * f_dot_dx;

    u[0] += f_x * H1 + dx * H2;
    u[1] += f_y * H1 + dy * H2;
    u[2] += f_z * H1 + dz * H2;    
  }

  return u;
}

//***********************************************************************//
vector<double> Call_Ds_SL(int limit_1, int limit_2,
			  double p_x, double p_y, double p_z,
			  double f0_dot_n0, double b,
			  struct xyz &particles,
			  double *lambda[3],
			  double DEL) {    

  vector<double> u(3,0);

  if ( abs(b)/DEL > 8.0 ) {

    for (size_t j = limit_1; j <= limit_2; j++) {
    
      double dx = p_x - particles.x[j];
      double dy = p_y - particles.y[j];
      double dz = p_z - particles.z[j];
    
      double r2 = dx * dx + dy * dy + dz * dz;
      double r = sqrt(r2);
    
      double area = particles.area[j];

      double f_x = lambda[0][j];
      double f_y = lambda[1][j];
      double f_z = lambda[2][j];
    
      double f_dot_dx = f_x * dx + f_y * dy + f_z * dz;
      
      double H1 = area / r;
      double H2 = f_dot_dx * area / (r2 * r);

      u[0] += f_x * H1 + dx * H2;
      u[1] += f_y * H1 + dy * H2;
      u[2] += f_z * H1 + dz * H2;
    }
  }
  else {
       
    double lam = b / DEL;
    double c1, c2, c3;
    get_reg_coeff_7ord(lam, c1, c2, c3);

    for (size_t j = limit_1; j <= limit_2; j++) {
      
      double dx = p_x - particles.x[j];
      double dy = p_y - particles.y[j];
      double dz = p_z - particles.z[j];
      
      double r = sqrt(dx * dx + dy * dy + dz * dz);
      
      double area = particles.area[j];

      double f_x = lambda[0][j] - f0_dot_n0 * particles.nx[j]; //use subtraction in Stokeslet
      double f_y = lambda[1][j] - f0_dot_n0 * particles.ny[j];
      double f_z = lambda[2][j] - f0_dot_n0 * particles.nz[j];
            
      double f_dot_dx = f_x * dx + f_y * dy + f_z * dz;

      double H1, H2;
      Stokes_SL_7ord(r, DEL, c1, c2, c3, H1, H2);
      
      H1 *= area;
      H2 *= f_dot_dx * area;
      
      u[0] += f_x * H1 + dx * H2;
      u[1] += f_y * H1 + dy * H2;
      u[2] += f_z * H1 + dz * H2;
    }
  }
  
  return u;
}

//***********************************************************************//
vector<double> Call_Ds_DL_onSurf(int limit_1, int limit_2,
				 double p_x, double p_y, double p_z,
				 const vector<double>& g0,
				 struct xyz &particles,
				 double *lambda[3],
				 double DEL) {    
  
  vector<double> u(3,0);
  
  for (size_t j = limit_1; j <= limit_2; j++) {
    
    double dx = p_x - particles.x[j];
    double dy = p_y - particles.y[j];
    double dz = p_z - particles.z[j];
    
    double r = sqrt(dx * dx + dy * dy + dz * dz);

    double g_x = lambda[0][j] - g0[0]; //use subtraction in Stresslet
    double g_y = lambda[1][j] - g0[1];
    double g_z = lambda[2][j] - g0[2];
            
    double g_dot_dx = g_x * dx + g_y * dy + g_z * dz;
    double n_dot_dx = particles.nx[j] * dx + particles.ny[j] * dy + particles.nz[j] * dz;
    double area = particles.area[j];
      
    double H3;
    Stokes_DL_7ord_onSurf(r, DEL, H3);

    double kernel = n_dot_dx * g_dot_dx * H3 * area;

    u[0] += dx * kernel;
    u[1] += dy * kernel;
    u[2] += dz * kernel;
  }
  
  return u;
}

//***********************************************************************//
vector<double> Call_Ds_DL(int limit_1, int limit_2,
			  double p_x, double p_y, double p_z,
			  const vector<double>& g0, const vector<double>& x0,
			  const vector<double>& n0, double b,
			  struct xyz &particles,
			  double *lambda[3],
			  double DEL) {    
  
  vector<double> u(3,0);

  if ( abs(b)/DEL > 8.0 ) {

    for (size_t j = limit_1; j <= limit_2; j++) {

      double dx = p_x - particles.x[j];
      double dy = p_y - particles.y[j];
      double dz = p_z - particles.z[j];
      
      double r2 = dx * dx + dy * dy + dz * dz;
      double r = sqrt(r2);
      
      double g_x = lambda[0][j];
      double g_y = lambda[1][j];
      double g_z = lambda[2][j];
      
      double g_dot_dx = g_x * dx + g_y * dy + g_z * dz;
      double n_dot_dx = particles.nx[j] * dx + particles.ny[j] * dy + particles.nz[j] * dz;
      double area = particles.area[j];
      
      double kernel = g_dot_dx * n_dot_dx /(r2*r2*r) * area;

      u[0] += dx * kernel;
      u[1] += dy * kernel;
      u[2] += dz * kernel;
    }
  } 
  else {
    
    double lam = b / DEL;
    double c1, c2, c3;
    get_reg_coeff_7ord(lam, c1, c2, c3);

    for (size_t j = limit_1; j <= limit_2; j++) {
    
      double dx = p_x - particles.x[j];
      double dy = p_y - particles.y[j];
      double dz = p_z - particles.z[j];
      
      double r = sqrt(dx * dx + dy * dy + dz * dz);

      double x_hat0 = particles.x[j] - x0[0];
      double x_hat1 = particles.y[j] - x0[1];
      double x_hat2 = particles.z[j] - x0[2];
      
      double n_x = particles.nx[j];
      double n_y = particles.ny[j];
      double n_z = particles.nz[j];
      
      double g_x = lambda[0][j] - g0[0]; //use subtraction in Stresslet
      double g_y = lambda[1][j] - g0[1];
      double g_z = lambda[2][j] - g0[2];
      
      double sig_sq = ( x_hat0 * x_hat0 + x_hat1 * x_hat1 + x_hat2 * x_hat2 )
	- 2.0 * b * ( x_hat0 * n0[0] + x_hat1 * n0[1] + x_hat2 * n0[2] );
      
      double g_dot_n0 = g_x * n0[0] + g_y * n0[1] + g_z * n0[2];
      double g_dot_x_hat = g_x * x_hat0 + g_y * x_hat1 + g_z * x_hat2;
      double n0_dot_n = n0[0] * n_x + n0[1] * n_y + n0[2] * n_z;
      double x_hat_dot_n = x_hat0 * n_x + x_hat1 * n_y + x_hat2 * n_z;
      
      double c_m1a = g_dot_n0 * n0_dot_n;
      double c_m1b = g_dot_x_hat * n0_dot_n  +  g_dot_n0 * x_hat_dot_n;
      
      double c_m2a = g_dot_x_hat * x_hat_dot_n;
      double c_m2b = b * c_m1b;
      
      double H2, H3;
      Stokes_DL_7ord(r, DEL, c1, c2, c3, H2, H3);
      
      double area = particles.area[j];
      H2 *= area;
      H3 *= area;

      double m1, m2;
      m1 = dx * c_m1a - n0[0] * c_m1b;	
      m2 = -sig_sq * m1 + dx * c_m2a + x_hat0 * c_m2b;
      
      u[0] += m1 * H2 + m2 * H3;

      m1 = dy * c_m1a - n0[1] * c_m1b;
      m2 = -sig_sq * m1 + dy * c_m2a + x_hat1 * c_m2b;
      
      u[1] += m1 * H2 + m2 * H3;
      
      m1 = dz * c_m1a - n0[2] * c_m1b;
      m2 = -sig_sq * m1 + dz * c_m2a + x_hat2 * c_m2b;
      
      u[2] += m1 * H2 + m2 * H3;
    }
  }
  
  return u;
}

//**********************************************************************//

vector<double> eval_SL_Sum_cluster_onSurf(double *lambda[3],
					  struct xyz &particles,
					  double p_x, double p_y, double p_z,
					  double f0_dot_n0,
					  size_t panel_index,
					  double DEL) {	
  
  vector<double> Sum(3,0);
  
  size_t limit_1 = tree[panel_index].members[0];
  size_t limit_2 = tree[panel_index].members[1];  
  
  double xc = tree[panel_index].xc;
  double yc = tree[panel_index].yc;
  double zc = tree[panel_index].zc;
    
  double tpx = p_x - xc;
  double tpy = p_y - yc;
  double tpz = p_z - zc;
  
  double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

  
  if (tree[panel_index].MAC < R_sq) { //particle and cluster well-separated
    
    vector<double> far_result(3,0);
    
    far_result = Call_Treecode_SL(p_x, p_y, p_z,
				  f0_dot_n0,
				  panel_index,
				  DEL);

    for (int k=0; k<3; k++) { 
      Sum[k] += far_result[k];
    }
  }
  else {

    if ((limit_2 - limit_1) < N0) { //if cluster is a leaf, use direct sum  
   
      // this can be the "nearby" case, and the del values have to be done separately
      
      vector<double> DS_result(3,0);
      DS_result = Call_Ds_SL_onSurf(limit_1, limit_2,
				    p_x, p_y, p_z,
				    f0_dot_n0,
				    particles,
				    lambda,
				    DEL);
      
      for (int k=0; k<3; k++) {
	Sum[k] += DS_result[k];
      }
    }
    else {//if cluster is not a leaf, look at children recursively
      for (int k=0; k<3; k++) {
	Sum[k] = 0.0;
      }
      
      size_t length = tree[panel_index].children.size();
      for (size_t i = 0; i < length; i++) {
	
	size_t index = tree[panel_index].children[i];
	vector<double> temp_res(3,0);
	temp_res = eval_SL_Sum_cluster_onSurf(lambda,
					      particles,
					      p_x, p_y, p_z,
					      f0_dot_n0,
					      index,
					      DEL);
	for (int k=0; k<3; k++) {
	  Sum[k] += temp_res[k];
	}
      }// end for i
    }// end if (limit_2 - limit_1) < N0)
  }// end if (MAC < theta)
  
  return Sum;
}
//**********************************************************************//

void eval_SL_Sum_cluster(double *lambda[3],
			 struct xyz &particles,
			 double p_x, double p_y, double p_z,
			 double f0_dot_n0, double b,
			 size_t panel_index,
			 double DEL,
			 vector<double>& Sum) {	
  
  
  size_t limit_1 = tree[panel_index].members[0];
  size_t limit_2 = tree[panel_index].members[1];  
  
  double xc = tree[panel_index].xc;
  double yc = tree[panel_index].yc;
  double zc = tree[panel_index].zc;
    
  double tpx = p_x - xc;
  double tpy = p_y - yc;
  double tpz = p_z - zc;
  
  double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

  
  if (tree[panel_index].MAC < R_sq) { //particle and cluster well-separated
    
    vector<double> far_result(3,0);
    
    far_result = Call_Treecode_SL(p_x, p_y, p_z,
				  f0_dot_n0,
				  panel_index,
				  DEL);

    for (int k=0; k<3; k++) { 
      Sum[k] += far_result[k];
    }
  }
  else {
    if ((limit_2 - limit_1) < N0) { //if cluster is a leaf, use direct sum
      
      // this can be the "nearby" case, and the del values have to be done separately
      
      vector<double> DS_result(3,0);
      DS_result = Call_Ds_SL(limit_1, limit_2,
			     p_x, p_y, p_z,
			     f0_dot_n0, b,
			     particles,
			     lambda,
			     DEL);
      
      for (int k=0; k<3; k++) {
	Sum[k] += DS_result[k];
      }
    }
    else {//if cluster is not a leaf, look at children recursively
      for (int k=0; k<3; k++) {
	Sum[k] = 0.0;
      }
      
      size_t length = tree[panel_index].children.size();
      for (size_t i = 0; i < length; i++) {
	
	size_t index = tree[panel_index].children[i];
	vector<double> temp_res(3,0);
	eval_SL_Sum_cluster(lambda,
			    particles,
			    p_x, p_y, p_z,
			    f0_dot_n0, b,
			    index,
			    DEL,
			    temp_res);
	for (int k=0; k<3; k++) {
	  Sum[k] += temp_res[k];
	}
      }// end for i
    }// end if (limit_2 - limit_1) < N0)
  }// end if (MAC < theta)
}
//**********************************************************************//

vector<double> eval_DL_Sum_cluster_onSurf(double *lambda[3],
					  struct xyz &particles,
					  double p_x, double p_y, double p_z,
					  const vector<double>& g0,
					  size_t panel_index,
					  double DEL) {	
  
  vector<double> Sum(3,0);
  
  size_t limit_1 = tree[panel_index].members[0];
  size_t limit_2 = tree[panel_index].members[1];  
  
  double xc = tree[panel_index].xc;
  double yc = tree[panel_index].yc;
  double zc = tree[panel_index].zc;
    
  double tpx = p_x - xc;
  double tpy = p_y - yc;
  double tpz = p_z - zc;
  
  double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

  
  if (tree[panel_index].MAC < R_sq) { //particle and cluster well-separated
    
    vector<double> far_result(3,0);
    
    far_result = Call_Treecode_DL(p_x, p_y, p_z,
				  g0,
				  panel_index,
				  DEL);

    for (int k=0; k<3; k++) { 
      Sum[k] += far_result[k];
    }
  }
  else {

    if ((limit_2 - limit_1) < N0) { //if cluster is a leaf, use direct sum  
   
      // this can be the "nearby" case, and the del values have to be done separately
      
      vector<double> DS_result(3,0);
      DS_result = Call_Ds_DL_onSurf(limit_1, limit_2,
				    p_x, p_y, p_z,
				    g0,
				    particles,
				    lambda,
				    DEL);
      
      for (int k=0; k<3; k++) {
	Sum[k] += DS_result[k];
      }
    }
    else {//if cluster is not a leaf, look at children recursively
      for (int k=0; k<3; k++) {
	Sum[k] = 0.0;
      }
      
      size_t length = tree[panel_index].children.size();
      for (size_t i = 0; i < length; i++) {
	
	size_t index = tree[panel_index].children[i];
	vector<double> temp_res(3,0);
	temp_res = eval_DL_Sum_cluster_onSurf(lambda,
					      particles,
					      p_x, p_y, p_z,
					      g0,
					      index,
					      DEL);
	for (int k=0; k<3; k++) {
	  Sum[k] += temp_res[k];
	}
      }// end for i
    }// end if (limit_2 - limit_1) < N0)
  }// end if (MAC < theta)
  
  return Sum;
}
//**********************************************************************//

void eval_DL_Sum_cluster(double *lambda[3],
			 struct xyz &particles,
			 double p_x, double p_y, double p_z,
			 const vector<double>& g0, const vector<double>& x0,
			 const vector<double>& n0, double b,
			 size_t panel_index,
			 double DEL,
			 vector<double>& Sum) {	
  
  
  size_t limit_1 = tree[panel_index].members[0];
  size_t limit_2 = tree[panel_index].members[1];  
  
  double xc = tree[panel_index].xc;
  double yc = tree[panel_index].yc;
  double zc = tree[panel_index].zc;
    
  double tpx = p_x - xc;
  double tpy = p_y - yc;
  double tpz = p_z - zc;
  
  double R_sq = tpx * tpx + tpy * tpy + tpz * tpz;

  
  if (tree[panel_index].MAC < R_sq) { //particle and cluster well-separated

    vector<double> far_result(3,0);
    far_result = Call_Treecode_DL(p_x, p_y, p_z,
				    g0,
				    panel_index,
				    DEL);
    
    for (int k=0; k<3; k++) { 
      Sum[k] += far_result[k];
    }
  } // if (tree[panel_index].MAC < R_sq)
  else {
    
    if ((limit_2 - limit_1) < N0) { //if cluster is a leaf, use direct sum
      
      // this can be the "nearby" case, and the del values have to be done separately
      vector<double> DS_result(3,0);
      DS_result = Call_Ds_DL(limit_1, limit_2,
			     p_x, p_y, p_z,
			     g0, x0,
			     n0, b,
			     particles,
			     lambda,
			     DEL);
      
      for (int k=0; k<3; k++) {
	Sum[k] += DS_result[k];
      }
    }// end if-else (tree[panel_index].MAC < R_sq)
    
    else {//if cluster is not a leaf, look at children recursively
      for (int k=0; k<3; k++) {
	Sum[k] = 0.0;
      }
      
      size_t length = tree[panel_index].children.size();
      for (size_t i = 0; i < length; i++) {
	
	size_t index = tree[panel_index].children[i];
	vector<double> temp_res(3,0);
	eval_DL_Sum_cluster(lambda,
			    particles,
			    p_x, p_y, p_z,
			    g0, x0,
			    n0, b,
			    index,
			    DEL,
			    temp_res);
	for (int k=0; k<3; k++) {
	  Sum[k] += temp_res[k];
	}
      }// end for i
    }// end if (limit_2 - limit_1) < N0)
  }// end if (MAC < theta)
}

//*************************************************************************************//
void eval_Stokes_SL_onSurf_tree(int surf_num, double h, double DEL,
				int N_quad, const vector<Surf_point>& Surfc,
				int N_target, const vector<Target_point>& Target,
				vector<double>& SL) {
  
  double theta = sqrt(sq_theta);
  double EightPI = 1.0/(8.0*PI);

  int N_cube = N_quad;
  
  tree.reserve(5000);
  leaf.reserve(5000);
  struct xyz particles(N_cube);
  double *lambda[3];
  lambda[0] = new double[N_cube];
  lambda[1] = new double[N_cube];
  lambda[2] = new double[N_cube];

  //create a physical copy of the input data
  //so the treecode can modify the data
  for (int i=0; i<N_cube; i++) {

    particles.x[i] = Surfc[i].x[0];
    particles.y[i] = Surfc[i].x[1];
    particles.z[i] = Surfc[i].x[2];
    particles.index[i] = -1;
    particles.old_index[i] = i;

    particles.nx[i] = Surfc[i].Nrml[0];
    particles.ny[i] = Surfc[i].Nrml[1];
    particles.nz[i] = Surfc[i].Nrml[2];

    particles.area[i] = Surfc[i].Area;

    lambda[0][i] = Surfc[i].f[0];
    lambda[1][i] = Surfc[i].f[1];
    lambda[2][i] = Surfc[i].f[2];
  }
    
  //***************** Set up tree *******************************
  long Start_total, Start_btree;
  long End_total, End_btree;
  
  Start_total = getTickCount(); // Get currenct CPU time
  Start_btree = getTickCount();
  
  build_tree_init(surf_num, N_cube);
  build_tree_3D_Recursive(0, lambda, particles, 0);
 
  End_btree = getTickCount();
  
  // //***************** Compute moments for each panel ************
  size_t size = tree.size();
  Cluster_Chev_Points(size);


  long moments_start = getTickCount();
  for (size_t i = 1; i < size; i++) {// skip root
    Panel_Moments_SL(i, lambda, particles, tree[i].moments);
  }
  long moments_end = getTickCount();
  cout << "     CPU for moments " << (moments_end-moments_start)/100
       << " seconds" << endl;


#pragma omp parallel for
  for (int i=0; i<N_target; i++) {    

    double p_x = Target[i].x;
    double p_y = Target[i].y;
    double p_z = Target[i].z;

    double f0_dot_n0 = dot_product(Surfc[i].f, Surfc[i].Nrml);

    // tree creation re-orders the data
    // we want the output to have the original ordering
    // original order is in particles.old_index

    vector<double> u(3,0);
    u = eval_SL_Sum_cluster_onSurf(lambda,
				   particles,
				   p_x, p_y, p_z,
				   f0_dot_n0,
				   0,
				   DEL);
    
    for (int k=0; k<3; k++) {
      SL[3*i+k] = u[k] * EightPI;
    }
  }// end for i

  End_total = getTickCount(); 
  
  cout << "     build tree time = " << (End_btree - Start_btree)
       << "/100 seconds" << endl;  
  cout << "     total treecode time = " << (End_total - Start_total)/100
       << " seconds" << endl;
  
  delete [] lambda[0];
  delete [] lambda[1];
  delete [] lambda[2];
  tree.clear();
  leaf.clear();
}


//*************************************************************************************//
void eval_Stokes_SL_offSurf_tree(int surf_num, double h, double DEL,
				 int N_quad, const vector<Surf_point>& Surfc,
				 int N_target, const vector<Target_point>& Target,
				 vector<double>& SL) {
  
  double theta = sqrt(sq_theta);
  double EightPI = 1.0/(8.0*PI);

  int N_cube = N_quad;
  
  tree.reserve(5000);
  leaf.reserve(5000);
  struct xyz particles(N_cube);
  double *lambda[3];
  lambda[0] = new double[N_cube];
  lambda[1] = new double[N_cube];
  lambda[2] = new double[N_cube];

  //create a physical copy of the input data
  //so the treecode can modify the data
  for (int i=0; i<N_cube; i++) {

    particles.x[i] = Surfc[i].x[0];
    particles.y[i] = Surfc[i].x[1];
    particles.z[i] = Surfc[i].x[2];
    particles.index[i] = -1;
    particles.old_index[i] = i;

    particles.nx[i] = Surfc[i].Nrml[0];
    particles.ny[i] = Surfc[i].Nrml[1];
    particles.nz[i] = Surfc[i].Nrml[2];

    particles.area[i] = Surfc[i].Area;

    lambda[0][i] = Surfc[i].f[0];
    lambda[1][i] = Surfc[i].f[1];
    lambda[2][i] = Surfc[i].f[2];
  }
    
  //***************** Set up tree *******************************
  long Start_total, Start_btree;
  long End_total, End_btree;
  
  Start_total = getTickCount(); // Get currenct CPU time
  Start_btree = getTickCount();
  
  build_tree_init(surf_num, N_cube);
  build_tree_3D_Recursive(0, lambda, particles, 0);
 
  End_btree = getTickCount();
  
  // //***************** Compute moments for each panel ************
  size_t size = tree.size();
  Cluster_Chev_Points(size);


  long moments_start = getTickCount();
  for (size_t i = 1; i < size; i++) {// skip root
    Panel_Moments_SL(i, lambda, particles, tree[i].moments);
  }
  long moments_end = getTickCount();
  cout << "     CPU for moments " << (moments_end-moments_start)/100
       << " seconds" << endl;


#pragma omp parallel for
  for (int i=0; i<N_target; i++) {    

    double p_x = Target[i].x;
    double p_y = Target[i].y;
    double p_z = Target[i].z;

    double f0_dot_n0 = dot_product(Target[i].nrst_f, Target[i].nrst_Nrml);

    double b = Target[i].b;

    // tree creation re-orders the data
    // we want the output to have the original ordering
    // original order is in particles.old_index
    
    vector<double> u(3,0);
    eval_SL_Sum_cluster(lambda,
			particles,
			p_x, p_y, p_z,
			f0_dot_n0, b,
			0,
			DEL,
			u);
    
    for (int k=0; k<3; k++) {
      SL[3*i+k] = u[k] * EightPI;
    }
  }// end for i

  End_total = getTickCount(); 
  
  cout << "     build tree time = " << (End_btree - Start_btree)
       << "/100 seconds" << endl;  
  cout << "     total treecode time = " << (End_total - Start_total)/100
       << " seconds" << endl;
  
  delete [] lambda[0];
  delete [] lambda[1];
  delete [] lambda[2];
  tree.clear();
  leaf.clear();
}


//*************************************************************************************//
void eval_Stokes_DL_onSurf_tree(int surf_num, double h, double DEL,
				int N_quad, const vector<Surf_point>& Surfc,
				int N_target, const vector<Target_point>& Target,
				vector<double>& DL) {
  
  double theta = sqrt(sq_theta);
  double EightPI6 = -6.0/(8.0*PI);

  int N_cube = N_quad;
  
  tree.reserve(5000);
  leaf.reserve(5000);
  struct xyz particles(N_cube);
  double *lambda[3];
  lambda[0] = new double[N_cube];
  lambda[1] = new double[N_cube];
  lambda[2] = new double[N_cube];

  //create a physical copy of the input data
  //so the treecode can modify the data
  for (int i=0; i<N_cube; i++) {

    particles.x[i] = Surfc[i].x[0];
    particles.y[i] = Surfc[i].x[1];
    particles.z[i] = Surfc[i].x[2];
    particles.index[i] = -1;
    particles.old_index[i] = i;

    particles.nx[i] = Surfc[i].Nrml[0];
    particles.ny[i] = Surfc[i].Nrml[1];
    particles.nz[i] = Surfc[i].Nrml[2];

    particles.area[i] = Surfc[i].Area;

    lambda[0][i] = Surfc[i].g[0];
    lambda[1][i] = Surfc[i].g[1];
    lambda[2][i] = Surfc[i].g[2];
  }
    
  //***************** Set up tree *******************************
  long Start_total, Start_btree;
  long End_total, End_btree;
  
  Start_total = getTickCount(); // Get currenct CPU time
  Start_btree = getTickCount();
  
  build_tree_init(surf_num, N_cube);
  build_tree_3D_Recursive(0, lambda, particles, 0);
 
  End_btree = getTickCount();
  
  // //***************** Compute moments for each panel ************
  size_t size = tree.size();
  Cluster_Chev_Points(size);


  long moments_start = getTickCount();
  for (size_t i = 1; i < size; i++) {// skip root
    Panel_Moments_DL(i, lambda, particles, tree[i].moments);
  }
  long moments_end = getTickCount();
  cout << "     CPU for moments " << (moments_end-moments_start)/100
       << " seconds" << endl;


#pragma omp parallel for
  for (int i=0; i<N_target; i++) {    

    double p_x = Target[i].x;
    double p_y = Target[i].y;
    double p_z = Target[i].z;

    vector<double> g0 = Surfc[i].g;

    // tree creation re-orders the data
    // we want the output to have the original ordering
    // original order is in particles.old_index

    vector<double> u(3,0);
    u = eval_DL_Sum_cluster_onSurf(lambda,
				   particles,
				   p_x, p_y, p_z,
				   g0,
				   0,
				   DEL);

    for (int k=0; k<3; k++) {
      DL[3*i+k] = u[k] * EightPI6 + 0.5 * g0[k]; // target point is on the boundary
    }
  } // end for i

  End_total = getTickCount(); 
  
  cout << "     build tree time = " << (End_btree - Start_btree)
       << "/100 seconds" << endl;  
  cout << "     total treecode time = " << (End_total - Start_total)/100
       << " seconds" << endl;
  
  delete [] lambda[0];
  delete [] lambda[1];
  delete [] lambda[2];
  tree.clear();
  leaf.clear();
}

//*************************************************************************************//
void eval_Stokes_DL_offSurf_tree(int surf_num, double h, double DEL,
				      int N_quad, const vector<Surf_point>& Surfc,
				      int N_target, const vector<Target_point>& Target,
				      vector<double>& DL) {
  
  double theta = sqrt(sq_theta);
  double EightPI6 = -6.0/(8.0*PI);

  int N_cube = N_quad;
  
  tree.reserve(5000);
  leaf.reserve(5000);
  struct xyz particles(N_cube);
  double *lambda[3];
  lambda[0] = new double[N_cube];
  lambda[1] = new double[N_cube];
  lambda[2] = new double[N_cube];

  //create a physical copy of the input data
  //so the treecode can modify the data
  for (int i=0; i<N_cube; i++) {

    particles.x[i] = Surfc[i].x[0];
    particles.y[i] = Surfc[i].x[1];
    particles.z[i] = Surfc[i].x[2];
    particles.index[i] = -1;
    particles.old_index[i] = i;

    particles.nx[i] = Surfc[i].Nrml[0];
    particles.ny[i] = Surfc[i].Nrml[1];
    particles.nz[i] = Surfc[i].Nrml[2];

    particles.area[i] = Surfc[i].Area;

    lambda[0][i] = Surfc[i].g[0];
    lambda[1][i] = Surfc[i].g[1];
    lambda[2][i] = Surfc[i].g[2];
  }
    
  //***************** Set up tree *******************************
  long Start_total, Start_btree;
  long End_total, End_btree;
  
  Start_total = getTickCount(); // Get currenct CPU time
  Start_btree = getTickCount();
  
  build_tree_init(surf_num, N_cube);
  build_tree_3D_Recursive(0, lambda, particles, 0);
 
  End_btree = getTickCount();
  
  // //***************** Compute moments for each panel ************
  size_t size = tree.size();
  Cluster_Chev_Points(size);


  long moments_start = getTickCount();
  for (size_t i = 1; i < size; i++) {// skip root
    Panel_Moments_DL(i, lambda, particles, tree[i].moments);
  }
  long moments_end = getTickCount();
  cout << "     CPU for moments " << (moments_end-moments_start)/100
       << " seconds" << endl;


#pragma omp parallel for
  for (int i=0; i<N_target; i++) {    

    double p_x = Target[i].x;
    double p_y = Target[i].y;
    double p_z = Target[i].z;

    vector<double> x0(3,0);
    x0[0] = Target[i].nrst_x;
    x0[1] = Target[i].nrst_y;
    x0[2] = Target[i].nrst_z;

    double b = Target[i].b;

    // tree creation re-orders the data
    // we want the output to have the original ordering
    // original order is in particles.old_index

    vector<double> u(3,0);
    eval_DL_Sum_cluster(lambda,
			particles,
			p_x, p_y, p_z,
			Target[i].nrst_g, x0,
			Target[i].nrst_Nrml, b,
			0,
			DEL,
			u);

    double xi = 0.0; // target point is outside the boundary
    if ( abs(b) < 1.0e-12 ) { // target point is on the boundary
      xi = 0.5; 
    }
    else if ( b < 0.0 ) { // target point is inside the boundary
      xi = 1.0;
    }
    
    for (int k=0; k<3; k++) {
      DL[3*i+k] = u[k] * EightPI6 + xi * Target[i].nrst_g[k];
    }
  } // end for i

  End_total = getTickCount(); 
  
  cout << "     build tree time = " << (End_btree - Start_btree)
       << "/100 seconds" << endl;  
  cout << "     total treecode time = " << (End_total - Start_total)/100
       << " seconds" << endl;
  
  delete [] lambda[0];
  delete [] lambda[1];
  delete [] lambda[2];
  tree.clear();
  leaf.clear();
}




