/*****************************************************************************/
/* mttr.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Mutual transit light curves						     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2011; Pal, A. (apal@szofi.net)					     */
/*****************************************************************************/

#include <stdio.h> 
#include <iostream>
#include <string.h> 
#include <memory.h> 
#include <stdlib.h> 
#include <math.h> 

//#ifndef MTTR_TEST_EXECUTABLE		/* we do not need lfit.h for testing */
//#include <lfit/lfit.h> 
//#endif

#ifdef  NBRV_TEST_EXECUTABLE 
#define STATIC 
#else 
#define STATIC  static 
#endif 

#include "icirc.h"
#include "scpolyint.h"
#include "elliptic.h"

using namespace std;

/*****************************************************************************/

double mttr_integral_primitive(double r,double c,double x)
{
 double	q2,s2,d2,sx,cx,w;
 double	rf,rd,rj;
 double	beta;
 double	iret;

 q2=r*r+c*c+2*r*c*cos(x);
 //d2=r*r+c*c-2*r*c;
 d2 = (r-c)*(r-c);
 s2=r*r+c*c+2*r*c;

 sx=sin(x/2);
 cx=cos(x/2);

 if ( 1.0<q2 )	q2=1.0;

 w=(1-q2)/(1-d2);
 if ( w<0.0 )	w=0.0;
 rf=carlson_elliptic_rf(w,sx*sx,1);
 rd=carlson_elliptic_rd(w,sx*sx,1);
 if ( r != c )	rj=carlson_elliptic_rj(w,sx*sx,1,q2/d2);
 else		rj=0.0;

 beta=atan2((c-r)*sx,(c+r)*cx);
 iret=-beta/3;
 iret+=x/6;

 w=cx/sqrt(1-d2);
 

 iret+=
    +2.0/ 9.0*c*r*sin(x)*sqrt(1-q2)
    +1.0/ 3.0*(1+2*r*r*r*r-4*r*r)*w*rf
    +2.0/ 9.0*r*c*(4-7*r*r-c*c+5*r*c)*w*rf
    -4.0/27.0*r*c*(4-7*r*r-c*c)*w*cx*cx*rd;
 if ( r != c )
	iret += 1.0/3.0*w*(r+c)/(r-c)*(rf-(q2-d2)/(3*d2)*rj);
 else
	iret -= 1.0/3.0*w*(r+c)*(q2-d2)*M_PI/(2*q2*sqrt(q2));

 if (rj != rj) {
   cout << "NAN in MTTR" << endl;
 }

 return(iret);
}

double mttr_integral_definite(double r,double c,double x0,double dx)
{
 double	dc,nx;
 double	ret;

 if ( c<=0.0 )
  {	if ( r<1.0 )
		ret=(1-(1-r*r)*sqrt(1-r*r))*dx/3.0;
	else /* this case implies r==1: */
		ret=dx/3.0;

	return(ret);
  }

 if ( dx<0.0 )
  {	x0+=dx;
	dx=-dx;
  }
 while ( x0<0.0 )	x0+=2*M_PI;
 while ( 2*M_PI<=x0 )	x0-=2*M_PI;

 ret=0.0;
 while ( 0.0<dx )
  {	dc=2*M_PI-x0;
	if ( dx<dc )	dc=dx,nx=x0+dx;
	else		nx=0.0;

	ret+=mttr_integral_primitive(r,c,x0+dc)-mttr_integral_primitive(r,c,x0);

	x0=nx;
	dx-=dc;
  }
 
 return(ret); 
}

/*****************************************************************************/

double mttr_flux_general(circle *circles,int ncircle,double c0,double c1,double c2)
{
 arc	*arcs,*a;
 int	i,narc;
 double	fc,f0;

 icirc_arclist_intersections(circles,ncircle,&arcs,&narc); 

 

 fc=0.0;

 for ( i=0 ; i<narc ; i++ )
  {	double	sign,x0,y0,r,p0,dp,p1,df0,df1,df2;
	double	x2,y2,r2;

	a=&arcs[i];
	if ( a->cidx==0 && a->noidx<=0 )
		sign=+1;
	else if ( a->cidx != 0 && a->noidx==1 && a->oidxs[0]==0 )
		sign=-1;
	else
		continue;

	x0=circles[a->cidx].x0;
	y0=circles[a->cidx].y0;
	r =circles[a->cidx].r;
	p0=a->phi0;
	dp=a->dphi;
	p1=p0+dp;

	x2=x0*x0;
	y2=y0*y0;
	r2=r*r;

	df0=0.5*r*(x0*(sin(p1)-sin(p0))+y0*(-cos(p1)+cos(p0)))+0.5*r*r*dp;

	if ( c1 != 0.0 )
	 {	double	delta,rho;
		delta=atan2(y0,x0);
		rho=sqrt(x2+y2);
		
		df1=mttr_integral_definite(r,rho,p0-delta,dp);
		/* fprintf(stderr,"p0-delta=%g dp=%g\n",p0-delta,dp); */
		/* fprintf(stderr,"r=%g rho=%g df1=%g\n",r,rho,df1);  */
	 }
	else
		df1=0.0;

	if ( c2 != 0.0 )
	 {	df2=(r/48.0)*(	+(24*(x2+y2)+12*r2)*r*dp
				-4*y0*(6*x2+2*y2+9*r2)*(cos(p1)-cos(p0))
				-24*x0*y0*r*(cos(2*p1)-cos(2*p0))
				-4*y0*r2*(cos(3*p1)-cos(3*p0))
				+4*x0*(2*x2+6*y2+9*r2)*(sin(p1)-sin(p0))
				-4*x0*r2*(sin(3*p1)-sin(3*p0))
				-r2*r*(sin(4*p1)-sin(4*p0)) );
	 }
	else
		df2=0.0;

	fc += sign*(c0*df0+c1*df1+c2*df2);


  }

 f0=2.0*M_PI*(c0/2.0+c1/3.0+c2/4.0);

 icirc_arclist_free(arcs,narc);

 if ( 0.0<f0 )
	return(fc/f0);
 else
	return(0.0);
}

double mttr_flux(double x1,double y1,double r1,double x2,double y2,double r2,double c)
{
 circle	circles[3];
 arc	*arcs,*a;
 int	i,narc;
 double	fc,f0;

 circles[0].x0=0.0;
 circles[0].y0=0.0;
 circles[0].r =1.0;
 circles[1].x0=x1;
 circles[1].y0=y1;
 circles[1].r =r1;
 circles[2].x0=x2;
 circles[2].y0=y2;
 circles[2].r =r2;

 icirc_arclist_intersections(circles,3,&arcs,&narc); 

 fc=0.0;

 for ( i=0 ; i<narc ; i++ )
  {	double	sign,x0,y0,r,p0,dp,p1,df0,df2;
	double	x2,y2,r2;

	a=&arcs[i];
	if ( a->cidx==0 && a->noidx<=0 )
		sign=+1;
	else if ( a->cidx != 0 && a->noidx==1 && a->oidxs[0]==0 )
		sign=-1;
	else
		continue;

	x0=circles[a->cidx].x0;
	y0=circles[a->cidx].y0;
	r =circles[a->cidx].r;
	p0=a->phi0;
	dp=a->dphi;
	p1=p0+dp;

	x2=x0*x0;
	y2=y0*y0;
	r2=r*r;

	df0=0.5*r*(x0*(sin(p1)-sin(p0))+y0*(-cos(p1)+cos(p0)))+0.5*r*r*dp;

	df2=(r/48.0)*(	+(24*(x2+y2)+12*r2)*r*dp
			-4*y0*(6*x2+2*y2+9*r2)*(cos(p1)-cos(p0))
			-24*x0*y0*r*(cos(2*p1)-cos(2*p0))
			-4*y0*r2*(cos(3*p1)-cos(3*p0))
			+4*x0*(2*x2+6*y2+9*r2)*(sin(p1)-sin(p0))
			-4*x0*r2*(sin(3*p1)-sin(3*p0))
			-r2*r*(sin(4*p1)-sin(4*p0)) );

	fc += sign*(df0-c*df2);


  }

 f0=M_PI*(1-0.5*c);

 icirc_arclist_free(arcs,narc);

 if ( 0.0<f0 )
	return(fc/f0);
 else
	return(0.0);
}

#ifndef MTTR_TEST_EXECUTABLE

STATIC double mf1u(double x1,double y1,double r1)
{
 circle	circles[2];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 return(mttr_flux_general(circles,2,1.0,0.0,0.0));
}

STATIC double mf1l(double x1,double y1,double r1,double cl)
{
 circle	circles[2];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 return(mttr_flux_general(circles,2,1-cl,cl,0.0));
}

STATIC double mf1q(double x1,double y1,double r1,double g1,double g2)
{
 circle	circles[2];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 return(mttr_flux_general(circles,2,1-g1-2*g2,g1+2*g2,g2));
}

STATIC double mf1s(double x1,double y1,double r1,double s)
{
 circle	circles[2];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 return(mttr_flux_general(circles,2,1,0.0,-s));
}

STATIC double mf2u(double x1,double y1,double r1,double x2,double y2,double r2)
{
 circle	circles[3];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 circles[2].x0=x2 ,circles[2].y0=y2 ,circles[2].r =r2;
 return(mttr_flux_general(circles,3,1.0,0.0,0.0));
}

STATIC double mf2l(double x1,double y1,double r1,double x2,double y2,double r2,double cl)
{
 circle	circles[3];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 circles[2].x0=x2 ,circles[2].y0=y2 ,circles[2].r =r2;
 return(mttr_flux_general(circles,2,1-cl,cl,0.0));
}

STATIC double mf2q(double x1,double y1,double r1,double x2,double y2,double r2,double g1,double g2)
{
 circle	circles[3];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 circles[2].x0=x2 ,circles[2].y0=y2 ,circles[2].r =r2;
 return(mttr_flux_general(circles,2,1-g1-2*g2,g1+2*g2,g2));
}

STATIC double mf2s(double x1,double y1,double r1,double x2,double y2,double r2,double s)
{
 circle	circles[3];
 circles[0].x0=0.0,circles[0].y0=0.0,circles[0].r =1.0;
 circles[1].x0=x1 ,circles[1].y0=y1 ,circles[1].r =r1 ;
 circles[2].x0=x2 ,circles[2].y0=y2 ,circles[2].r =r2;
 return(mttr_flux_general(circles,2,1,0.0,-s));
}

STATIC int mttr_funct_flux(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mttr_flux(vars[0],vars[1],vars[2],vars[3],vars[4],vars[5],idvs[0]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2,v3,v4,v5;
	v0=vars[0],v1=vars[1],v2=vars[2];
	v3=vars[3],v4=vars[4],v5=vars[5];
	d=0.001;
	diff[0]=(mttr_flux(v0+d/2,v1,v2,v3,v4,v5,idvs[0])-mttr_flux(v0-d/2,v1,v2,v3,v4,v5,idvs[0]))/d;
	diff[1]=(mttr_flux(v0,v1+d/2,v2,v3,v4,v5,idvs[0])-mttr_flux(v0,v1-d/2,v2,v3,v4,v5,idvs[0]))/d;
	diff[2]=(mttr_flux(v0,v1,v2+d/2,v3,v4,v5,idvs[0])-mttr_flux(v0,v1,v2-d/2,v3,v4,v5,idvs[0]))/d;
	diff[3]=(mttr_flux(v0,v1,v2,v3+d/2,v4,v5,idvs[0])-mttr_flux(v0,v1,v2,v3-d/2,v4,v5,idvs[0]))/d;
	diff[4]=(mttr_flux(v0,v1,v2,v3,v4+d/2,v5,idvs[0])-mttr_flux(v0,v1,v2,v3,v4-d/2,v5,idvs[0]))/d;
	diff[5]=(mttr_flux(v0,v1,v2,v3,v4,v5+d/2,idvs[0])-mttr_flux(v0,v1,v2,v3,v4,v5-d/2,idvs[0]))/d;
  }

 return(0);
}

STATIC int lmf1u(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf1u(vars[0],vars[1],vars[2]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2;
	v0=vars[0],v1=vars[1],v2=vars[2];
	d=0.001;
	diff[0]=(mf1u(v0+d/2,v1,v2)-mf1u(v0-d/2,v1,v2))/d;
	diff[1]=(mf1u(v0,v1+d/2,v2)-mf1u(v0,v1-d/2,v2))/d;
	diff[2]=(mf1u(v0,v1,v2+d/2)-mf1u(v0,v1,v2-d/2))/d;
  }

 return(0);
}

STATIC int lmf1l(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf1l(vars[0],vars[1],vars[2],idvs[0]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2;
	v0=vars[0],v1=vars[1],v2=vars[2];
	d=0.001;
	diff[0]=(mf1l(v0+d/2,v1,v2,idvs[0])-mf1l(v0-d/2,v1,v2,idvs[0]))/d;
	diff[1]=(mf1l(v0,v1+d/2,v2,idvs[0])-mf1l(v0,v1-d/2,v2,idvs[0]))/d;
	diff[2]=(mf1l(v0,v1,v2+d/2,idvs[0])-mf1l(v0,v1,v2-d/2,idvs[0]))/d;
  }

 return(0);
}

STATIC int lmf1s(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf1s(vars[0],vars[1],vars[2],idvs[0]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2;
	v0=vars[0],v1=vars[1],v2=vars[2];
	d=0.001;
	diff[0]=(mf1s(v0+d/2,v1,v2,idvs[0])-mf1s(v0-d/2,v1,v2,idvs[0]))/d;
	diff[1]=(mf1s(v0,v1+d/2,v2,idvs[0])-mf1s(v0,v1-d/2,v2,idvs[0]))/d;
	diff[2]=(mf1s(v0,v1,v2+d/2,idvs[0])-mf1s(v0,v1,v2-d/2,idvs[0]))/d;
  }

 return(0);
}

STATIC int lmf1q(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf1q(vars[0],vars[1],vars[2],idvs[0],idvs[1]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2;
	v0=vars[0],v1=vars[1],v2=vars[2];
	d=0.001;
	diff[0]=(mf1q(v0+d/2,v1,v2,idvs[0],idvs[1])-mf1q(v0-d/2,v1,v2,idvs[0],idvs[1]))/d;
	diff[1]=(mf1q(v0,v1+d/2,v2,idvs[0],idvs[1])-mf1q(v0,v1-d/2,v2,idvs[0],idvs[1]))/d;
	diff[2]=(mf1q(v0,v1,v2+d/2,idvs[0],idvs[1])-mf1q(v0,v1,v2-d/2,idvs[0],idvs[1]))/d;
  }

 return(0);
}

STATIC int lmf2u(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf2u(vars[0],vars[1],vars[2],vars[3],vars[4],vars[5]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2,v3,v4,v5;
	v0=vars[0],v1=vars[1],v2=vars[2];
	v3=vars[3],v4=vars[4],v5=vars[5];
	d=0.001;
	diff[0]=(mf2u(v0+d/2,v1,v2,v3,v4,v5)-mf2u(v0-d/2,v1,v2,v3,v4,v5))/d;
	diff[1]=(mf2u(v0,v1+d/2,v2,v3,v4,v5)-mf2u(v0,v1-d/2,v2,v3,v4,v5))/d;
	diff[2]=(mf2u(v0,v1,v2+d/2,v3,v4,v5)-mf2u(v0,v1,v2-d/2,v3,v4,v5))/d;
	diff[3]=(mf2u(v0,v1,v2,v3+d/2,v4,v5)-mf2u(v0,v1,v2,v3-d/2,v4,v5))/d;
	diff[4]=(mf2u(v0,v1,v2,v3,v4+d/2,v5)-mf2u(v0,v1,v2,v3,v4-d/2,v5))/d;
	diff[5]=(mf2u(v0,v1,v2,v3,v4,v5+d/2)-mf2u(v0,v1,v2,v3,v4,v5-d/2))/d;
  }

 return(0);
}

STATIC int lmf2l(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf2l(vars[0],vars[1],vars[2],vars[3],vars[4],vars[5],idvs[0]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2,v3,v4,v5;
	v0=vars[0],v1=vars[1],v2=vars[2];
	v3=vars[3],v4=vars[4],v5=vars[5];
	d=0.001;
	diff[0]=(mf2l(v0+d/2,v1,v2,v3,v4,v5,idvs[0])-mf2l(v0-d/2,v1,v2,v3,v4,v5,idvs[0]))/d;
	diff[1]=(mf2l(v0,v1+d/2,v2,v3,v4,v5,idvs[0])-mf2l(v0,v1-d/2,v2,v3,v4,v5,idvs[0]))/d;
	diff[2]=(mf2l(v0,v1,v2+d/2,v3,v4,v5,idvs[0])-mf2l(v0,v1,v2-d/2,v3,v4,v5,idvs[0]))/d;
	diff[3]=(mf2l(v0,v1,v2,v3+d/2,v4,v5,idvs[0])-mf2l(v0,v1,v2,v3-d/2,v4,v5,idvs[0]))/d;
	diff[4]=(mf2l(v0,v1,v2,v3,v4+d/2,v5,idvs[0])-mf2l(v0,v1,v2,v3,v4-d/2,v5,idvs[0]))/d;
	diff[5]=(mf2l(v0,v1,v2,v3,v4,v5+d/2,idvs[0])-mf2l(v0,v1,v2,v3,v4,v5-d/2,idvs[0]))/d;
  }

 return(0);
}

STATIC int lmf2q(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf2q(vars[0],vars[1],vars[2],vars[3],vars[4],vars[5],idvs[0],idvs[1]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2,v3,v4,v5;
	v0=vars[0],v1=vars[1],v2=vars[2];
	v3=vars[3],v4=vars[4],v5=vars[5];
	d=0.001;
	diff[0]=(mf2q(v0+d/2,v1,v2,v3,v4,v5,idvs[0],idvs[1])-mf2q(v0-d/2,v1,v2,v3,v4,v5,idvs[0],idvs[1]))/d;
	diff[1]=(mf2q(v0,v1+d/2,v2,v3,v4,v5,idvs[0],idvs[1])-mf2q(v0,v1-d/2,v2,v3,v4,v5,idvs[0],idvs[1]))/d;
	diff[2]=(mf2q(v0,v1,v2+d/2,v3,v4,v5,idvs[0],idvs[1])-mf2q(v0,v1,v2-d/2,v3,v4,v5,idvs[0],idvs[1]))/d;
	diff[3]=(mf2q(v0,v1,v2,v3+d/2,v4,v5,idvs[0],idvs[1])-mf2q(v0,v1,v2,v3-d/2,v4,v5,idvs[0],idvs[1]))/d;
	diff[4]=(mf2q(v0,v1,v2,v3,v4+d/2,v5,idvs[0],idvs[1])-mf2q(v0,v1,v2,v3,v4-d/2,v5,idvs[0],idvs[1]))/d;
	diff[5]=(mf2q(v0,v1,v2,v3,v4,v5+d/2,idvs[0],idvs[1])-mf2q(v0,v1,v2,v3,v4,v5-d/2,idvs[0],idvs[1]))/d;
  }

 return(0);
}

STATIC int lmf2s(double *vars,double *idvs,double *ret,double *diff)
{
 double	r;

 r=mf2s(vars[0],vars[1],vars[2],vars[3],vars[4],vars[5],idvs[0]);

 if ( ret != NULL )
	*ret=r;
 if ( diff != NULL )
  {	double	d,v0,v1,v2,v3,v4,v5;
	v0=vars[0],v1=vars[1],v2=vars[2];
	v3=vars[3],v4=vars[4],v5=vars[5];
	d=0.001;
	diff[0]=(mf2s(v0+d/2,v1,v2,v3,v4,v5,idvs[0])-mf2s(v0-d/2,v1,v2,v3,v4,v5,idvs[0]))/d;
	diff[1]=(mf2s(v0,v1+d/2,v2,v3,v4,v5,idvs[0])-mf2s(v0,v1-d/2,v2,v3,v4,v5,idvs[0]))/d;
	diff[2]=(mf2s(v0,v1,v2+d/2,v3,v4,v5,idvs[0])-mf2s(v0,v1,v2-d/2,v3,v4,v5,idvs[0]))/d;
	diff[3]=(mf2s(v0,v1,v2,v3+d/2,v4,v5,idvs[0])-mf2s(v0,v1,v2,v3-d/2,v4,v5,idvs[0]))/d;
	diff[4]=(mf2s(v0,v1,v2,v3,v4+d/2,v5,idvs[0])-mf2s(v0,v1,v2,v3,v4-d/2,v5,idvs[0]))/d;
	diff[5]=(mf2s(v0,v1,v2,v3,v4,v5+d/2,idvs[0])-mf2s(v0,v1,v2,v3,v4,v5-d/2,idvs[0]))/d;
  }

 return(0);
}


/* lfitfunction mttr [] = */
/*  { */
/*         { "mttr_flux", LFITFUNCTION_DIFFERENTIABLE,  6, 1, mttr_funct_flux 	}, */
/*         { "mttr1u",    LFITFUNCTION_DIFFERENTIABLE,  3, 0, lmf1u 	}, */
/*         { "mttr1l",    LFITFUNCTION_DIFFERENTIABLE,  3, 1, lmf1l 	}, */
/*         { "mttr1s",    LFITFUNCTION_DIFFERENTIABLE,  3, 1, lmf1s 	}, */
/*         { "mttr1q",    LFITFUNCTION_DIFFERENTIABLE,  3, 2, lmf1q 	}, */
/*         { "mttr2u",    LFITFUNCTION_DIFFERENTIABLE,  6, 0, lmf2u 	}, */
/*         { "mttr2l",    LFITFUNCTION_DIFFERENTIABLE,  6, 1, lmf2l 	}, */
/*         { "mttr2s",    LFITFUNCTION_DIFFERENTIABLE,  6, 1, lmf2s 	}, */
/*         { "mttr2q",    LFITFUNCTION_DIFFERENTIABLE,  6, 2, lmf2q 	}, */
/* 	{ NULL, 0, 0, 0, NULL } */
/*  }; */

#else
int main(int argc,char *argv[])
{
 return(0);
}
#endif

/*****************************************************************************/
             
