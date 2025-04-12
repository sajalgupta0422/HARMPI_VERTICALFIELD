//
//  madfield.c
//  HARM2D
//
//  Created by Alexander Tchekhovskoy on 7/6/15.
//  Copyright (c) 2015 Alexander Tchekhovskoy. All rights reserved.
//

#include "decs.h"
#include <float.h>

void getmax_densities(double (*prim)[N2M][N3M][NPR], double *rhomax, double *umax)
{
  int i,j,k;

  *rhomax = 0;
  *umax = 0;
  
  ZLOOP {
    if (prim[i][j][k][RHO] > *rhomax)   *rhomax = prim[i][j][k][RHO];
    if (prim[i][j][k][UU] > *umax )    *umax = prim[i][j][k][UU];
  }

#ifdef MPI
  MPI_Allreduce(MPI_IN_PLACE,rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,umax,  1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
}

double get_maxprimvalrpow(double (*prim)[N2M][N3M][NPR], double rpow, int m )
{
  int i,j,k;
  double X[NDIM], V[NDIM];
  double r;
  
  double val;
  double maxval = -DBL_MAX;
  
  ZLOOP {
    coord(i,j,k,CENT,X);
    bl_coord_vec(X,V);
    
    r = V[1];
    val = pow(r,rpow)*prim[i][j][k][m];
    if ( val > maxval ) {
      maxval = val;
    }
  }
  
#ifdef MPI
  MPI_Allreduce(MPI_IN_PLACE,&maxval,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  return(maxval);
  
}

int normalize_field_local_nodivb(double targbeta, double rhomax, double amax,
                                 double (*prim)[N2M][N3M][NPR],
                                 double (*A)[N2+D2][N3+D3], int dir)
{
  //is called only two times: in init_torus() and init_torus_grb(), only if WHICHFIELD = MADFIELD
   // Normalizes field components individually
  // Note the different array dimensions for A vs prim
  // This suggests vector potential A is truly staggered in storage as it is corner centered array
  // while B-fields in prim are logically staggered and are face centered by interpretation but are stored as cell centered
  // This is why the indexing is different for A and prim
  int i,j,k;
  double ratc_ij;
  double compute_rat(double (*prim)[N2M][N3M][NPR], double (*A)[N2+D2][N3+D3],
                     double rhomax, double amax, double targbeta, int loc, int i, int j, int k);
  
  bound_prim(prim);
  
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) { 
    //The loop goes from 0 to Nq, in every direction for q = 1,2,3 because
    //The vector potential A is defined on cell corners
  //For N cells, we need N+1 corners
    //So if we have N1 cells in r-direction, we need N1+1 points for A.
    // This extra point is crucial for properly computing B-fields at the last cell using curl of A.
    //cell centered ratio in this cell
    /*In more detail:
    The A parameter is passed because this is part of an iterative process:
    Normalize B fields in primitive array
    Recompute vector potential A from the normalized B fields
    Recompute B fields from new A
    Repeat until desired configuration is achieved
    So while  normalize_field_local_nodivb doesn't directly use A, it's part of a sequence where A and B are kept consistent with each other. 
    The loop bounds might be set to match A's dimensions for consistency in this iterative process. */
    ratc_ij   = compute_rat(prim, A, rhomax, amax, targbeta, CENT, i, j, k);
    
    // normalize staggered field primitive
    if(dir == 1) prim[i][j][k][B1] *= ratc_ij;
    if(dir == 2) prim[i][j][k][B2] *= ratc_ij;
    if(dir == 3) prim[i][j][k][B3] *= ratc_ij;
  }
  
  //bound_allprim(STAGEM1,t,prim,pstag,ucons, 1, USEMPI);
  
  return(0);
}

#define MYSMALL (1.e-300)

//Returns: factor to multiply field components by to get the desired
//value of beta: targbeta = p_g/p_mag
double compute_rat(double (*prim)[N2M][N3M][NPR], double (*A)[N2+D2][N3+D3],
                  double rhomax, double amax, double targbeta, int loc, int i, int j, int k)
{
  double bsq_ij,pg_ij,beta_ij,rat_ij;
  struct of_geom geom;
  double X[NDIM], V[NDIM];
  double rat, ratc;
  double profile;
  double r;
  double rho, u;
  double compute_profile( double (*prim)[N2M][N3M][NPR], double amax, double aphipow, int loc, int i, int j, int k );
  
  get_geometry(i, j, k, loc, &geom);
  coord(i,j,k,loc,X);
  bl_coord_vec(X,V);
  r = V[1];
  
  bsq_ij = bsq_calc(prim[i][j][k], &geom);
  
  rho = prim[i][j][k][RHO];
  //use the following instead of MACP0A1(prim,i,j,k,UU) because the latter
  //can be perturbed by random noise, which we want to avoid
  u = global_kappa * pow(rho, gam) / (gam - 1.);

  //EOSMARK
  pg_ij= (gam-1)*u;
  beta_ij=2*pg_ij/(bsq_ij+MYSMALL);
  rat_ij = sqrt(beta_ij / targbeta); //ratio at CENT
  
  //rescale rat_ij so that:
  // rat_ij = 1 inside the main body of torus
  // rat_ij = 0 outside of main body of torus
  // rat_ij ~ rho in between
  //ASSUMING DENSITY HAS ALREADY BEEN NORMALIZED -- SASMARK
  profile = compute_profile(prim,amax,aphipow,loc,i,j,k);
  rat_ij *= profile;
  
  return(rat_ij);
}

double compute_profile( double (*prim)[N2M][N3M][NPR], double amax, double aphipow, int loc, int i, int j, int k )
{
  /*
  The  compute_profile function creates a tapering profile for magnetic field strength that:
Returns 1 in the main torus body
Returns 0 outside the torus
Smoothly transitions between these values:*/
  double X[NDIM], V[NDIM], r;
  struct of_geom geom;
  double profile;
  
  get_geometry(i, j, k, loc, &geom);
  coord(i,j,k,loc,X);
  bl_coord_vec(X,V);
  r = V[1];
  profile = ( log10(pow(r,aphipow)*prim[i][j][k][RHO]/amax+SMALL) + 3. ) / 1.0;
  if(profile<0.) profile = 0.;
  if(profile>1.) profile = 1.;
  profile=1.;
  return(profile);
}

//compute vector potential assuming B_\phi = 0 and zero flux at poles
//(not tested in non-axisymmetric field distribution but in principle should work)
//is called only two times: in init_torus() and init_torus_grb(), only if WHICHFIELD = MADFIELD
int compute_vpot_from_gdetB1( double (*prim)[N2M][N3M][NPR],  double (*A)[N2+D2][N3+D3] )
{
  /*
  Physics:
    In axisymmetry, B^r = (1/√g)∂(√g A_φ)/∂θ
    Integrating this: A_φ = ∫B^r √g dθ
    This preserves ∇·B = 0 by construction
  Mathematics:
    Integrates B1 (radial field) along θ-direction to get A_φ
    Uses two integrations:
    From pole (θ=0) to equator
    From opposite pole (θ=π) to equator
Ensures A_φ = 0 at both poles
  */
  int i, j, k;
  int jj;
  int cj;
  int dj, js, je, jsb, jeb;
  int dointegration;
  struct of_geom geom;
  double gdet;
  int finalstep;
  
  //first, bound to ensure consistency of magnetic fields across tiles
  bound_prim(prim);
  
  if( mpi_dims[2] == 1 ) {
    //1-cpu version
    for (i=0; i<N1+D1; i++) {
      for (k=0; k<N3+D3; k++) {
        //zero out starting element of vpot
        A[i][0][k] = 0.0;
        //integrate vpot along the theta line
        for (j=0; j<N2/2; j++) {
          get_geometry(i, j, k, CENT, &geom);
          gdet = geom.g;
          //take a loop along j-line at a fixed i,k and integrate up vpot
          A[i][j+1][k] = A[i][j][k] + prim[i][j][k][B1]*gdet*dx[2];
        }
        A[i][N2][k] = 0.0;
        //integrate vpot along the theta line
        for (j=N2; j>N2/2; j--) {
          get_geometry(i, j-1, k, CENT, &geom);
          gdet = geom.g;
          //take a loop along j-line at a fixed i,k and integrate up vpot
          A[i][j-1][k] = A[i][j][k] - prim[i][j-1][k][B1]*gdet*dx[2];
        }
      }
    }
  }
  else {
#ifdef MPI
    for( cj = 0; cj < mpi_dims[2]/2; cj++ ) {
      if( mpi_coords[2] == cj ){
        dj = +1;
        js = 0;
        jsb = 0;
        je = N2;
        jeb = N2 - 1;
        dointegration = 1;
      }
      else if( mpi_dims[2] - mpi_coords[2] - 1 == cj ){
        dj = -1;
        js = N2;
        jsb = N2-1;
        je = 0;
        jeb = 0;
        dointegration = 1;
      }
      else {
        dointegration = 0;  //skip directly to bounding
      }
      
      if( 1 == dointegration ) {
        //then it's the turn of the current row of CPUs to pick up where the previous row has left it off
        //since pstag is bounded unlike A, use pstag[B3] as temporary space to trasnfer values of A[3] between CPUs
        //initialize lowest row of A[3]
        for (i=0; i<N1+D1; i++) {
          for (k=0; k<N3+D3; k++) {
            //zero out or copy starting element of vpot
            if( 0 == cj ) {
              //if CPU is at physical boundary, initialize (zero out) A[3]
              A[i][js][k] = 0.0;
            }
            else {
              //else copy B[3] (which was bounded below) -> A[3]
              A[i][js][k] = prim[i][jsb-dj][k][B3];
            }
            //integrate vpot along the theta line
            for (j=js; j!=je; j+=dj) {
              get_geometry(i, j-js+jsb, k, CENT, &geom);
              gdet = geom.g;
              //take a loop along j-line at a fixed i,k and integrate up vpot
              A[i][j+dj][k] = A[i][j][k] + dj * prim[i][j-js+jsb][k][B1]*gdet*dx[2];
            }
            //copy A[3] -> B[3] before bounding
            prim[i][jeb][k][B3] = A[i][je][k];
          }
        }
      }
      //just in case, wait until all CPUs get here
      MPI_Barrier(MPI_COMM_WORLD);
      //bound here
      bound_prim(prim);
    }
    //ensure consistency of vpot across the midplane
    if( mpi_coords[2] == mpi_dims[2]/2 ) {
      for (i=0; i<N1+D1; i++) {
        for (k=0; k<N3+D3; k++) {
          A[i][0][k] = prim[i][-1][k][B3];
        }
      }
    }
#endif
  }  
  //need to zero out prim[B3] everywhere
  ZSLOOP(-N1G,N1-1+N1G,-N2G,N2-1+N2G,-N3G,N3-1+N3G){
    prim[i][j][k][B3]=0.0;
  }
  return(0);
}


