/*
  n_body_lc.cpp

  Josh Carter, 2013

  Computes light curve from array of positions, radii, fluxes and quadratic limb darkening coefficients.
  Utilizes code by A. Pal (c) 2011, to compute overlap integrals of multiple dark bodies occulting
	a limb darkened star.  Please cite Pal (2012) whenever using this code for publication:

	MNRAS (2012) 420 (2): 1630-1635. doi: 10.1111/j.1365-2966.2011.20151.x

  This code will work for ANY number of limb-darkened spherical bodies in ANY configuration.

*/

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "icirc.h"
#include "n_body_lc.h"

#define PMAX	20 // Maximum number of bodies, for memory use.
#define ZPERS	1  // Observer perspective along z direction (-1 @ -z, 1 @ z)

double mttr_flux_general(circle *circles,int ncircle,double c0,double c1,double c2); //A. Pal's code

using namespace std;

void comb_sort(double ** pos, int input[], int size) {  // Quick comb sort to order bodies by z coordinate.
    int swap;
    int gap = size;
    bool swapped = false;
 
    while ((gap > 1) || swapped) {
        if (gap > 1) {
            gap = (int)((double)gap / 1.247330950103979);
        }
 
        swapped = false;
 
        for (int i = 0; gap + i < size; ++i) {
	  if ((pos[input[i]][2] - pos[input[i + gap]][2])*ZPERS > 0) {
                swap = input[i];
                input[i] = input[i + gap];
                input[i + gap] = swap;
                swapped = true;
            }
        }
    }
}

int uniq(int pairs[][2],int count, int uniq[]) {
  int n=2,l=2,k1,k2;

  uniq[0] = pairs[0][0];
  uniq[1] = pairs[0][1];

  for (int i = 1; i < count; i++) {
    k1 = 0; k2 = 0;
    for (int j = 0; j < n; j++) {
      if (pairs[i][0] != uniq[j]) k1++;	
      if (pairs[i][1] != uniq[j]) k2++;
    }
    if (k1 == n) uniq[n++] = pairs[i][0];
    if (k2 == n) uniq[n++] = pairs[i][1];
  }

  return n;
}


double occultn(double ** pos, double * radii, double * u1, double * u2, double * fluxes,int N) {

  int i,j;
  double tr = 0,rad,x0,y0,zij,rij;

  int pairs[PMAX][2];

  int count = 0; // tabulate overlapping pairs first - quicker than using A. Pal's generic code.
  for (j = 0; j < N; j++) { // N*(N-1)/2 computations - faster than doing arc intersects.
    for (i = 0; i < j; i++) {
      zij = ((pos[i][0]-pos[j][0])*(pos[i][0]-pos[j][0])+
	     (pos[i][1]-pos[j][1])*(pos[i][1]-pos[j][1]));
      rij = (radii[i]+radii[j]);
      if (zij < rij*rij) {
	pairs[count][0] = i;
	pairs[count++][1] = j;
      }
    }
  }

  if (count == 0) return 1; // no overlaps return 1
  
  int un[count*2];
  int cn = uniq(pairs,count,un); // Find unique pairs.
  
  circle circles[cn];

  int sorted[cn];
  for (i = 0; i < cn; i++) 
    sorted[i] = un[i];
  
  comb_sort(pos,sorted,cn); // Sort by z coordinate

  for (i = 0; i < cn; i++) { // record circles in list
    circles[i].x0 = pos[sorted[i]][0];
    circles[i].y0 = pos[sorted[i]][1];
    circles[i].r = radii[sorted[i]];
  }

  double mtt;
  for (i = 0; i< cn-1; i++) {
    if (fluxes[sorted[i]] != 0) { // compute overlap if occulted body has non-zero flux
      rad = (circles+i)[0].r;
      x0 = (circles+i)[0].x0;
      y0 = (circles+i)[0].y0;
      for (j = 0; j < cn-i; j++) { //compute hierarchy of overlaps.
	
	(circles+i+j)[0].x0 = ((circles+i+j)[0].x0-x0)/rad;
	(circles+i+j)[0].y0 = ((circles+i+j)[0].y0-y0)/rad;
	(circles+i+j)[0].r /= rad;
	
      }
      mtt = (1-mttr_flux_general((circles+i),cn-i,
			      1-u1[sorted[i]]-2*u2[sorted[i]],
				 u1[sorted[i]]+2*u2[sorted[i]],u2[sorted[i]]))*fluxes[sorted[i]];
	
      tr += mtt;
    }
  }

  return 1-tr;
}
