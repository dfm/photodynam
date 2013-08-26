/*****************************************************************************/
/* elliptic.h								     */
/*****************************************************************************/

#ifndef	__ELLIPTIC_H_INCLUDED
#define	__ELLIPTIC_H_INCLUDED	1

/*****************************************************************************/

double	carlson_elliptic_rc(double x,double y);
double	carlson_elliptic_rf(double x,double y,double z);
double	carlson_elliptic_rd(double x,double y,double z);
double	carlson_elliptic_rj(double x,double y,double z,double p);

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

double	elliptic_complete_first(double k);
double	elliptic_complete_second(double k);
double	elliptic_complete_third(double n,double k);

double	elliptic_incomplete_first (double sx,double k);
double	elliptic_incomplete_second(double sx,double k);
double	elliptic_incomplete_third (double sx,double n,double k);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                                               
              
