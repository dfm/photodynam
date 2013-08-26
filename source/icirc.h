/*****************************************************************************/
/* icirc.h								     */
/*****************************************************************************/

#ifndef	__ICIRC_H_INCLUDED
#define	__ICIRC_H_INCLUDED	1

/*****************************************************************************/

typedef struct
 {	double	x0,y0;
	double	r;
 } circle;

typedef struct
 {	int	cidx;
	double	phi0,dphi;
	int	noidx;
	int	*oidxs;
 } arc;

int	icirc_arclist_intersections(circle *circles,int ncircle,arc **rarcs,int *rnarc);
int	icirc_arclist_free(arc *arcs,int narc);

/*****************************************************************************/

#endif

/*****************************************************************************/
                                       
