/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

#include "decs.h"

/***************************************************************************/
/***************************************************************************
    coord():
    -------
       -- given the indices i,j and location in the cell, return with 
          the values of X1,X2 there;  
       -- the locations are defined by : 
           -----------------------
           |                     |
           |                     |
           |FACE1   CENT         |
           |                     |
           |CORN    FACE2        |
           ----------------------
***************************************************************************/

/**************************************************************
 * The indices (i,j,k) by themselves refer to the corner (CORN) of a cell. 
 * All other positions (FACE1, FACE2, FACE3, CENT, EDGE1, EDGE2, EDGE3) are defined relative to this corner position.
 * (i,j)-------FACE1--------(i+1,j)
   |           |            |
   |           |            |
FACE2    CENT(i+½,j+½)   FACE2
   |           |            |
   |           |            |
(i,j+1)----FACE1------(i+1,j+1)

CORN = (i,j)         : Cell corner
FACE1 = (i+1,j+½)   : x1-face center
FACE2 = (i+½,j)     : x2-face center
CENT = (i+½,j+½)    : Cell center
EDGE1 = (i+½,j)     : x1-edge
 

 3D Cell Indexing (i,j,k):
FACE1 (x1-face):          FACE2 (x2-face):          FACE3 (x3-face):
X[1]: i      (corner)     X[1]: i+½   (center)      X[1]: i+½   (center)
X[2]: j+½    (center)     X[2]: j     (corner)      X[2]: j+½   (center)
X[3]: k+½    (center)     X[3]: k+½   (center)      X[3]: k     (corner)
    
    ●---------------●
   /|              /|
  / |             / |
 /  |   FACE1    /  |
●---+------------●  |
|   |            |  |
|   ●------------+--● 
|  /    FACE2    |  /
| /              | /
|/     FACE3     |/
●----------------●

EDGE1 (x1-edge):          EDGE2 (x2-edge):          EDGE3 (x3-edge):
X[1]: i+½    (center)     X[1]: i     (corner)      X[1]: i     (corner)
X[2]: j      (corner)     X[2]: j+½   (center)      X[2]: j     (corner)
X[3]: k      (corner)     X[3]: k     (corner)      X[3]: k+½   (center)

CENT (cell center):        CORN (cell corner):
X[1]: i+½    (center)     X[1]: i     (corner)
X[2]: j+½    (center)     X[2]: j     (corner)
X[3]: k+½    (center)     X[3]: k     (corner)
 */
void coord(int i, int j, int k, int loc, double *X)
{
  //Staggered grid is best explained here.
  //The following if conditions tells us the location of FACES, EDGES, CENT and CORN wrt to (i,j,k) indices
  /*
  Key points:
+0.5 means centered in that direction
No +0.5 means on the corner in that direction
startx[] is the starting coordinate
dx[] is the grid spacing
mpi_startn[] handles MPI domain decomposition offsets
  */
  if(loc == FACE1) {
    X[1] = startx[1] + (i+mpi_startn[1])*dx[1] ; //No +0.5: on corner in x1
    X[2] = startx[2] + (j+mpi_startn[2] + 0.5)*dx[2] ; //+0.5: centered in x2
    X[3] = startx[3] + (k+mpi_startn[3] + 0.5)*dx[3] ; //+0.5: centered in x3
  }
  else if(loc == FACE2) {
    X[1] = startx[1] + (i+mpi_startn[1] + 0.5)*dx[1] ;
    X[2] = startx[2] + (j+mpi_startn[2])*dx[2] ;
    X[3] = startx[3] + (k+mpi_startn[3] + 0.5)*dx[3] ;
  }
  else if(loc == FACE3) {
    X[1] = startx[1] + (i+mpi_startn[1] + 0.5)*dx[1] ;
    X[2] = startx[2] + (j+mpi_startn[2] + 0.5)*dx[2] ;
    X[3] = startx[3] + (k+mpi_startn[3])*dx[3] ;
  }
  else if(loc == EDGE1) {
    X[1] = startx[1] + (i+mpi_startn[1] + 0.5)*dx[1] ; //+0.5: centered in x1
    X[2] = startx[2] + (j+mpi_startn[2])*dx[2] ; //No +0.5: on corner in x2
    X[3] = startx[3] + (k+mpi_startn[3])*dx[3] ; //No +0.5: on corner in x3
  }
  else if(loc == EDGE2) {
    X[1] = startx[1] + (i+mpi_startn[1])*dx[1] ;
    X[2] = startx[2] + (j+mpi_startn[2] + 0.5)*dx[2] ;
    X[3] = startx[3] + (k+mpi_startn[3])*dx[3] ;
  }
  else if(loc == EDGE3) {
    X[1] = startx[1] + (i+mpi_startn[1])*dx[1] ;
    X[2] = startx[2] + (j+mpi_startn[2])*dx[2] ;
    X[3] = startx[3] + (k+mpi_startn[3] + 0.5)*dx[3] ;
  }
  else if(loc == CENT) {
    X[1] = startx[1] + (i+mpi_startn[1] + 0.5)*dx[1] ;
    X[2] = startx[2] + (j+mpi_startn[2] + 0.5)*dx[2] ;
    X[3] = startx[3] + (k+mpi_startn[3] + 0.5)*dx[3] ;
  }
  else {
    X[1] = startx[1] + (i+mpi_startn[1])*dx[1] ; //No +0.5: on corners in all directions
    X[2] = startx[2] + (j+mpi_startn[2])*dx[2] ;
    X[3] = startx[3] + (k+mpi_startn[3])*dx[3] ;
  }

  return ;
}

/* assumes gcov has been set first; returns determinant */
double gdet_func(double gcov[][NDIM]) 
{
  int i,j,k;
  int permute[NDIM]; 
  double gcovtmp[NDIM][NDIM];
  double detg;

  for( i = 0 ; i < NDIM*NDIM ; i++ ) {  gcovtmp[0][i] = gcov[0][i]; }
  if( LU_decompose( gcovtmp,  permute ) != 0  ) { 
    fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
    fail(FAIL_METRIC);
  }
  detg = 1.;
  DLOOPA detg *= gcovtmp[j][j];
  return( sqrt(fabs(detg)) );

}

/* invert gcov to get gcon */
void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
  invert_matrix( gcov, gcon );
}

/***************************************************************************/
/***************************************************************************
  conn_func():
  -----------

   -- this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   --  where i = {1,2,3,4} corresponds to {t,r,theta,phi}

***************************************************************************/

/* Sets the spatial discretization in numerical derivatives : */
#define DELTA 1.e-5

/* NOTE: parameter hides global variable */
void conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM])
{
	int i,j,k,l ;
	double tmp[NDIM][NDIM][NDIM] ;
	double Xh[NDIM],Xl[NDIM] ;
	double gh[NDIM][NDIM] ;
	double gl[NDIM][NDIM] ;

	for(k=0;k<NDIM;k++) {
		for(l=0;l<NDIM;l++) Xh[l] = X[l] ;
		for(l=0;l<NDIM;l++) Xl[l] = X[l] ;
		Xh[k] += DELTA ;
		Xl[k] -= DELTA ;
		gcov_func(Xh,gh) ;
		gcov_func(Xl,gl) ;

		for(i=0;i<NDIM;i++)
		for(j=0;j<NDIM;j++) 
			conn[i][j][k] = (gh[i][j] - gl[i][j])/(Xh[k] - Xl[k]) ;
	}

	/* now rearrange to find \Gamma_{ijk} */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++) 
		tmp[i][j][k] = 0.5*(conn[j][i][k] + conn[k][i][j] - conn[k][j][i]) ;

	/* finally, raise index */
	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++)
	for(k=0;k<NDIM;k++)  {
		conn[i][j][k] = 0. ;
		for(l=0;l<NDIM;l++) conn[i][j][k] += geom->gcon[i][l]*tmp[l][j][k] ;
	}

	/* done! */
}

/* NOTE: parameter hides global variable */
void dxdxp_func(double *X, double dxdxp[][NDIM])
{
  int i,j,k,l ;
  double Xh[NDIM],Xl[NDIM] ;
  double Vh[NDIM],Vl[NDIM] ;

  if(BL){
    for(k=0;k<NDIM;k++) {
      for(l=0;l<NDIM;l++) Xh[l] = X[l] ;
      for(l=0;l<NDIM;l++) Xl[l] = X[l] ;
      Xh[k] += DELTA ;
      Xl[k] -= DELTA ;

      bl_coord_vec(Xh,Vh) ;
      bl_coord_vec(Xl,Vl) ;
      
      for(j=0;j<NDIM;j++)
        dxdxp[j][k] = (Vh[j]-Vl[j])/(Xh[k] - Xl[k]) ;
    }
  }
  else{
    for(i=0;i<NDIM;i++) {
      for(j=0;j<NDIM;j++) {
        dxdxp[i][j] = delta(i,j);
      }
    }
  }
}

/* Lowers a contravariant rank-1 tensor to a covariant one */
void lower(double *ucon, struct of_geom *geom, double *ucov)
{

	ucov[0] = geom->gcov[0][0]*ucon[0] 
		+ geom->gcov[0][1]*ucon[1] 
		+ geom->gcov[0][2]*ucon[2] 
		+ geom->gcov[0][3]*ucon[3] ;
	ucov[1] = geom->gcov[1][0]*ucon[0] 
		+ geom->gcov[1][1]*ucon[1] 
		+ geom->gcov[1][2]*ucon[2] 
		+ geom->gcov[1][3]*ucon[3] ;
	ucov[2] = geom->gcov[2][0]*ucon[0] 
		+ geom->gcov[2][1]*ucon[1] 
		+ geom->gcov[2][2]*ucon[2] 
		+ geom->gcov[2][3]*ucon[3] ;
	ucov[3] = geom->gcov[3][0]*ucon[0] 
		+ geom->gcov[3][1]*ucon[1] 
		+ geom->gcov[3][2]*ucon[2] 
		+ geom->gcov[3][3]*ucon[3] ;

        return ;
}

/* Raises a covariant rank-1 tensor to a contravariant one */
void raise(double *ucov, struct of_geom *geom, double *ucon)
{

	ucon[0] = geom->gcon[0][0]*ucov[0] 
		+ geom->gcon[0][1]*ucov[1] 
		+ geom->gcon[0][2]*ucov[2] 
		+ geom->gcon[0][3]*ucov[3] ;
	ucon[1] = geom->gcon[1][0]*ucov[0] 
		+ geom->gcon[1][1]*ucov[1] 
		+ geom->gcon[1][2]*ucov[2] 
		+ geom->gcon[1][3]*ucov[3] ;
	ucon[2] = geom->gcon[2][0]*ucov[0] 
		+ geom->gcon[2][1]*ucov[1] 
		+ geom->gcon[2][2]*ucov[2] 
		+ geom->gcon[2][3]*ucov[3] ;
	ucon[3] = geom->gcon[3][0]*ucov[0] 
		+ geom->gcon[3][1]*ucov[1] 
		+ geom->gcon[3][2]*ucov[2] 
		+ geom->gcon[3][3]*ucov[3] ;

        return ;
}

/* load local geometry into structure geom */
void get_geometry(int ii, int jj, int kk, int loc, struct of_geom *geom)
{
	int j,k ;

	//-new DLOOP geom->gcov[j][k] = gcov[ii][jj][kk][j][k] ;
	//-new DLOOP geom->gcon[j][k] = gcon[ii][jj][kk][j][k] ;
	for(j=0;j<=NDIM*NDIM-1;j++){
	  geom->gcon[0][j] = gcon[ii][jj][kk][loc][0][j];
	  geom->gcov[0][j] = gcov[ii][jj][kk][loc][0][j];
	}
	geom->g = gdet[ii][jj][kk][loc] ;
	icurr = ii ;
	jcurr = jj ;
        kcurr = kk ;
	pcurr = loc ;
}

/* load coordinates into V array */
void get_phys_coord_vec(int ii, int jj, int kk, double *V)
{
  int j ;
  
  SLOOPA V[j] = phys_coords[j][ii][jj][kk];
}

/* load r-coordinate value into *r */
void get_phys_coord_r(int ii, int jj, int kk, double *r)
{
  *r = phys_coords[1][ii][jj][kk];
}

/* load coordinates value into r, theta, phi */
void get_phys_coord(int ii, int jj, int kk, double *r, double *theta, double *phi)
{
  *r     = phys_coords[1][ii][jj][kk];
  *theta = phys_coords[2][ii][jj][kk];
  *phi   = phys_coords[3][ii][jj][kk];
}

#undef DELTA

/* Minkowski metric; signature +2 */
double mink(int i, int j)
{
	if(i == j) {
		if(i == 0) return(-1.) ;
		else return(1.) ;
	}
	else return(0.) ;
}

/* Boyer-Lindquist ("bl") metric functions */
void blgset(int i, int j, int k, struct of_geom *geom)
{
  // Metric terms computed at cell centers
  // Used for both regular variables and magnetic fields 
  // Suggests geometric factors handle the effective staggering 
	double r,th,phi,X[NDIM] ;

	coord(i,j,k,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;

	if(th < 0) th *= -1. ;
	if(th > M_PI) th = 2.*M_PI - th ;

	geom->g = bl_gdet_func(r,th,phi) ;
	bl_gcov_func(r,th,phi,geom->gcov) ;
	bl_gcon_func(r,th,phi,geom->gcon) ;
}

double bl_gdet_func(double r, double th, double phi)
{
	double a2,r2 ;

	a2 = a*a ;
	r2 = r*r ;
	return( 
		r*r*fabs(sin(th))*(1. + 0.5*(a2/r2)*(1. + cos(2.*th)))
	) ;
}

void bl_gcov_func(double r, double th, double phi, double gcov[][NDIM])
{
	int j,k ;
	double sth,cth,s2,a2,r2,DD,mu ;

	DLOOP gcov[j][k] = 0. ;

	sth = fabs(sin(th)) ;
	s2 = sth*sth ;
	cth = cos(th) ;
	a2 = a*a ;
	r2 = r*r ;
	DD = 1. - 2./r + a2/r2 ;
	mu = 1. + a2*cth*cth/r2 ;
	
	gcov[TT][TT] = -(1. - 2./(r*mu)) ;
	gcov[TT][3] = -2.*a*s2/(r*mu) ;
	gcov[3][TT] = gcov[TT][3] ;
	gcov[1][1] = mu/DD ;
	gcov[2][2] = r2*mu ;
	gcov[3][3] = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu)) ;

}

void bl_gcon_func(double r, double th, double phi, double gcon[][NDIM])
{
	int j,k ;
	double sth,cth,a2,r2,r3,DD,mu ;

	DLOOP gcon[j][k] = 0. ;

	sth = sin(th) ;
	cth = cos(th) ;

#if(COORDSINGFIX && BL)
	if (fabs(sth) < SINGSMALL) {
	  if(sth>=0) sth=SINGSMALL;
	  if(sth<0) sth=-SINGSMALL;
	}
#endif

	a2 = a*a ;
	r2 = r*r ;
	r3 = r2*r ;
	DD = 1. - 2./r + a2/r2 ;
	mu = 1. + a2*cth*cth/r2 ;

	gcon[TT][TT] = -1. - 2.*(1. + a2/r2)/(r*DD*mu) ;
	gcon[TT][3] = -2.*a/(r3*DD*mu) ;
	gcon[3][TT] = gcon[TT][3] ;
	gcon[1][1] = DD/mu ;
	gcon[2][2] = 1./(r2*mu) ;
	gcon[3][3] = (1. - 2./(r*mu))/(r2*sth*sth*DD) ;


}

/**********************************************************************
 * For transmissive boundary implementation, this staggering is crucial because:
  1. Face-centered quantities:
    Magnetic fields (B1, B2, B3) are naturally face-centered
    For transmissive boundaries, you need to properly interpolate these face-centered values

    void bound_x2dn_transmissive(double prim[][N2M][N3M][NPR]) {
    // For face-centered B-fields:
    // Need to properly interpolate at faces
    for (i=-N1G; i<N1+N1G; i++) {
        for (k=-N3G; k<N3+N3G; k++) {
            // B2 is face-centered in theta direction
            // Use appropriate staggering for interpolation
            double B2_face = 0.5*(prim[i][0][k][B2] + prim[i][-1][k][B2]);
            
            // Extrapolate into ghost zones maintaining face-centering
            for (j=-N2G; j<0; j++) {
                prim[i][j][k][B2] = B2_face;
            }
        }
    }
}

2. Cell-centered quantities:
Density, pressure, velocities are cell-centered (CENT)
Transmissive boundaries need different interpolation for these
// For cell-centered quantities
for (m=0; m<NPR; m++) {
    if (m != B1 && m != B2 && m != B3) {  // Skip B-fields
        // Use cell-centered interpolation
        double val_cent = prim[i][0][k][m];
        
        // Extrapolate maintaining cell-centering
        for (j=-N2G; j<0; j++) {
            prim[i][j][k][m] = val_cent;
        }
    }
}

3. Corner values:
Important for computing gradients and fluxes
Your transmissive boundary needs to maintain consistency between corners, faces, and centers
 * 
 */