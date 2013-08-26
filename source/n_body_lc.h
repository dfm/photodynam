#ifndef _N_BODY_LC
#define _N_BODY_LC

#include "n_body_state.h"

/*
  n_body_lc.cpp

  Josh Carter, 2013

  Computes light curve from array of positions, radii, fluxes and quadratic limb darkening coefficients.
  Utilizes code by A. Pal (c) 2011, to compute overlap integrals of multiple dark bodies occulting
	a limb darkened star.  Please cite Pal (2012) whenever using this code for publication:

	MNRAS (2012) 420 (2): 1630-1635. doi: 10.1111/j.1365-2966.2011.20151.x

  This code will work for ANY number of limb-darkened spherical bodies in ANY configuration.

*/

// pos is NX3 matrix of coordinates for N bodies.  radii is length N array of object radii. 
// u1,u2 are length N arrays of quadratic limb-darkening coefficients. fl is length N array of
// relative flux contributions (fraction of total flux).  Is is assumed that sum(fl) = 1. 
//
double occultn(double ** pos, double * radii, double * u1, double * u2, double * fl, int N);

#endif
