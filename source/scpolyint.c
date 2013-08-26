/*****************************************************************************/
/* scpolyint.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* A library for integrating (a+r*cos(x))^p*(b+r*sin(x))^q*sin(x) like stuff */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2011; Pal, A. (apal@szofi.net)					     */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "scpolyint.h"

/*****************************************************************************/

typedef struct
 {	double	mpq1;
	double	mpq2;
	double	impq;
	double	impqc;
	double	impqs;
 } pqgrid;

double scpoly_integrate(double a,double b,double r,double x1,double x2,double **ccoeff,double **scoeff,int pmax,int qmax)
{
 int	p,q;
 pqgrid	**g,*glist;
 double	sx1,sx2,cx1,cx2,ret;

 glist=(pqgrid *)alloca((pmax+1)*(qmax+1)*sizeof(pqgrid));
 g=(pqgrid **)alloca((pmax+1)*sizeof(pqgrid *));
 for ( p=0 ; p<=pmax ; p++ )
  {	g[p]=glist+p*(qmax+1);		}

 sx1=sin(x1);
 sx2=sin(x2);
 cx1=cos(x1);
 cx2=cos(x2);

 g[0][0].mpq1=1;
 g[0][0].mpq2=1;
 g[0][0].impq=x2-x1;
 g[0][0].impqc=+sx2-sx1;
 g[0][0].impqs=-cx2+cx1;

 ret=g[0][0].impqc*ccoeff[0][0]+g[0][0].impqs*scoeff[0][0];

 for ( p=0 ; p<=pmax ; p++ )
  {	for ( q=(p?0:1) ; q<=qmax ; q++ )
	 {	double	wc,ws;

		if ( 0<p )
		 {	g[p][q].mpq1=g[p-1][q].mpq1*(a+r*cx1);
			g[p][q].mpq2=g[p-1][q].mpq1*(a+r*cx2);
			g[p][q].impq=a*g[p-1][q].impq+r*g[p-1][q].impqc;
		 }
		else
		 {	g[p][q].mpq1=g[p][q-1].mpq1*(b+r*sx1);
			g[p][q].mpq2=g[p][q-1].mpq1*(b+r*sx2);
			g[p][q].impq=b*g[p][q-1].impq+r*g[p][q-1].impqs;
		 }
		wc=+(g[p][q].mpq2*sx2-g[p][q].mpq1*sx1);
		if ( 0<p )	wc+=r*p*g[p-1][q].impq+a*p*g[p-1][q].impqc;
		if ( 0<q )	wc+=b*q*g[p][q-1].impqc;
		g[p][q].impqc=wc/(1+p+q);
		
		ws=-(g[p][q].mpq2*cx2-g[p][q].mpq1*cx1);
		if ( 0<p )	ws+=a*p*g[p-1][q].impqs;
		if ( 0<q )	ws+=r*q*g[p][q-1].impq+b*q*g[p][q-1].impqs;
		g[p][q].impqs=ws/(1+p+q);
		
		ret+=g[p][q].impqc*ccoeff[p][q]+g[p][q].impqs*scoeff[p][q];

	 }
  }

 return(ret);
}

/*****************************************************************************/
                                                           
