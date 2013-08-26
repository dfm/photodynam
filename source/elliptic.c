/*****************************************************************************/
/* elliptic.c								     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Functions related to the calculation of elliptic integrals		     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* (c) 2007; Pal, A. (apal@szofi.elte.hu)				     */
/* Based on the implementation of Numerical Recipes (Press et al. 1992)	     */
/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "elliptic.h"

/*****************************************************************************/

#define		FMIN(a,b)	((a)<(b)?(a):(b))
#define		FMAX(a,b)	((a)<(b)?(a):(b))
#define		SQR(a)		((a)*(a))

#define C1 0.3
#define C2 (1.0/7.0)
#define C3 0.375
#define C4 (9.0/22.0)

double carlson_elliptic_rc(double x,double y)
{
 double alamb,ave,s,w,xt,yt,ans;

 if ( y > 0.0 )
  {	xt=x;
	yt=y;
	w=1.0;
  }
 else
  {	xt=x-y;
	yt = -y;
	w=sqrt(x)/sqrt(xt);
  }
 do 
  {	alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
	xt=0.25*(xt+alamb);
	yt=0.25*(yt+alamb);
	ave=(1.0/3.0)*(xt+yt+yt);
	s=(yt-ave)/ave;
  } while ( fabs(s) > 0.0012 );

 ans=w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);

 return(ans);
}

#undef	C4
#undef	C3
#undef	C2
#undef	C1

double carlson_elliptic_rf(double x,double y,double z)
{
 double	alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;
 xt=x;
 yt=y;
 zt=z;
 do 
  {	sqrtx=sqrt(xt);
	sqrty=sqrt(yt);
	sqrtz=sqrt(zt);
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
	xt=0.25*(xt+alamb);
	yt=0.25*(yt+alamb);
	zt=0.25*(zt+alamb);
	ave=(1.0/3.0)*(xt+yt+zt);
	delx=(ave-xt)/ave;
	dely=(ave-yt)/ave;
	delz=(ave-zt)/ave;
  } while ( fabs(delx) > 0.0025 || fabs(dely) > 0.0025 || fabs(delz) > 0.0025 );
 e2=delx*dely-delz*delz;
 e3=delx*dely*delz;
 return((1.0+((1.0/24.0)*e2-(0.1)-(3.0/44.0)*e3)*e2+(1.0/14.0)*e3)/sqrt(ave));
}

#define C1 (3.0/14.0)
#define C2 (1.0/6.0)
#define C3 (9.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.25*C3)
#define C6 (1.5*C4)

double carlson_elliptic_rd(double x,double y,double z)
{
 double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,
	sqrtx,sqrty,sqrtz,sum,xt,yt,zt,ans;

 xt=x;
 yt=y;
 zt=z;
 sum=0.0;
 fac=1.0;
 do
  {	sqrtx=sqrt(xt);
	sqrty=sqrt(yt);
	sqrtz=sqrt(zt);
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
	sum+=fac/(sqrtz*(zt+alamb));
	fac=0.25*fac;
	xt=0.25*(xt+alamb);
	yt=0.25*(yt+alamb);
	zt=0.25*(zt+alamb);
	ave=0.2*(xt+yt+3.0*zt);
	delx=(ave-xt)/ave;
	dely=(ave-yt)/ave;
	delz=(ave-zt)/ave;
 } while ( fabs(delx) > 0.0015 || fabs(dely) > 0.0015 || fabs(delz) > 0.0015 );
 ea=delx*dely;
 eb=delz*delz;
 ec=ea-eb;
 ed=ea-6.0*eb;
 ee=ed+ec+ec;
 ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
 return(ans);
}

#undef	C6
#undef	C5
#undef	C4
#undef	C3
#undef	C2
#undef	C1

#define C1 (3.0/14.0)
#define C2 (1.0/3.0)
#define C3 (3.0/22.0)
#define C4 (3.0/26.0)
#define C5 (0.75*C3)
#define C6 (1.5*C4)
#define C7 (0.5*C2)
#define C8 (C3+C3)

double carlson_elliptic_rj(double x,double y,double z,double p)
{
 double	a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,
	ed,ee,fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;

 sum=0.0;
 fac=1.0;
 if ( p > 0.0 )
  {	xt=x;
	yt=y;
	zt=z;
	pt=p;
	a=b=rcx=0.0;
  }
 else
  {	xt=FMIN(FMIN(x,y),z);
	zt=FMAX(FMAX(x,y),z);
	yt=x+y+z-xt-zt;
	a=1.0/(yt-p);
	b=a*(zt-yt)*(yt-xt);
	pt=yt+b;
	rho=xt*zt/yt;
	tau=p*pt/yt;
	rcx=carlson_elliptic_rc(rho,tau);
  }
 do
  {	sqrtx=sqrt(xt);
	sqrty=sqrt(yt);
	sqrtz=sqrt(zt);
	alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
	alpha=SQR(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
	beta=pt*SQR(pt+alamb);
	sum += fac*carlson_elliptic_rc(alpha,beta);
	fac=0.25*fac;
	xt=0.25*(xt+alamb);
	yt=0.25*(yt+alamb);
	zt=0.25*(zt+alamb);
	pt=0.25*(pt+alamb);
	ave=0.2*(xt+yt+zt+pt+pt);
	delx=(ave-xt)/ave;
	dely=(ave-yt)/ave;
	delz=(ave-zt)/ave;
	delp=(ave-pt)/ave;
 } while ( fabs(delx)>0.0015 || fabs(dely)>0.0015 || fabs(delz)>0.0015 || fabs(delp)>0.0015 );
 ea=delx*(dely+delz)+dely*delz;
 eb=delx*dely*delz;
 ec=delp*delp;
 ed=ea-3.0*ec;
 ee=eb+2.0*delp*(ea-ec);

 ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
	+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));

 if ( p <= 0.0 ) ans=a*(b*ans+3.0*(rcx-carlson_elliptic_rf(xt,yt,zt)));

 return(ans);
}

#undef	C6
#undef	C5
#undef	C4
#undef	C3
#undef	C2
#undef	C1
#undef	C8
#undef	C7

#undef			SQR
#undef			FMAX
#undef			FMIN

/*****************************************************************************/

double elliptic_complete_first(double ak)
{
 double	q;
 q=(1.0-ak)*(1.0+ak);
 if ( q>0.0 )
	return(carlson_elliptic_rf(0.0,q,1.0));
 else
	return(-1.0);
}

double elliptic_complete_second(double ak)
{
 double	q;
 q=(1.0-ak)*(1.0+ak);
 if ( q>0.0 )
	return(carlson_elliptic_rf(0.0,q,1.0)-(ak*ak)*carlson_elliptic_rd(0.0,q,1.0)/3.0);
 else if ( q<0.0 )
	return(-1.0);
 else
	return(1.0);
}

double elliptic_complete_third(double en,double ak)
{
 double	q;
 q=(1.0-ak)*(1.0+ak);
 if ( q>0.0 )
	return(carlson_elliptic_rf(0.0,q,1.0)+en*carlson_elliptic_rj(0.0,q,1.0,1.0-en)/3.0);
 else
	return(-1.0);
}

/*****************************************************************************/

/* this is the primitive function of 1/\sqrt{1-k^2\sin^2x}, sx:=\sin x, ak:=k*/
double elliptic_incomplete_first(double sx,double ak)
{
 double	s2x,c2x,k2,p2;
 k2=ak*ak;
 s2x=sx*sx;
 c2x=1.0-s2x;
 p2=k2*s2x;
 return(sx*carlson_elliptic_rf(c2x,1.0-p2,1.0));
}

/* this is the primitive function of \sqrt{1-k^2\sin^2x}, sx:=\sin x, ak:=k*/
double elliptic_incomplete_second(double sx,double ak)
{
 double	s2x,c2x,k2,p2;

 k2=ak*ak;
 s2x=sx*sx;
 c2x=1.0-s2x;
 p2=k2*s2x;
 return(sx*carlson_elliptic_rf(c2x,1.0-p2,1.0)-k2*sx*s2x*carlson_elliptic_rd(c2x,1.0-p2,1.0)/3.0);
}

/* this is the primitive function of 1/((1-n\sin^2x)\sqrt{1-k^2\sin^2x}), 
  sx:=\sin x, en:=n, ak:=k */
double elliptic_incomplete_third(double sx,double en,double ak)
{
 double	s2x,c2x,k2,p2;

 k2=ak*ak;
 s2x=sx*sx;
 c2x=1.0-s2x;
 p2=k2*s2x;
 return(sx*carlson_elliptic_rf(c2x,1.0-p2,1.0)+en*sx*s2x*carlson_elliptic_rj(c2x,1.0-p2,1.0,1.0-en*s2x)/3.0);
}

/*****************************************************************************/
                                         
