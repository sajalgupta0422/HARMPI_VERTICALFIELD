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

/* all diagnostics subroutine */

void append_rank(char *name)
{
  char suff[MAXLEN];
  //file suffixes that will be different by MPI rank
  sprintf(suff, "_%04d", mpi_rank);
  strcat(name,suff);
}

/*
diag() is called with different  call_code values at different points in the simulation:
 INIT_OUT: At initialization if not restarting
 DIVB_OUT: After restart
*/
void diag(call_code)
int call_code ;
{
  int di, dj, dk;
  char efnam[MAXLEN], dfnam[MAXLEN],ifnam[MAXLEN] ;
  int i,j,k,m ;
  FILE *dump_file;
  double U[NPR],pp,e,rmed,divb,divbmax,divbmax_global,e_fin,m_fin,gamma ;
  struct of_geom geom ;
  struct of_state q ;
  int imax,jmax,kmax ;
  int istart,istop,jstart,jstop,kstart,kstop;
  
  static double e_init,m_init ;
  static FILE *ener_file = NULL;

#if(DOENER)
  if(call_code==INIT_OUT || NULL == ener_file) {
    /* set things up */
    strcpy(efnam, "ener.out");
    append_rank(efnam);
    ener_file = fopen(efnam,"a") ;
    if(ener_file==NULL) {
      fprintf(stderr,"error opening energy output file, %s\n", efnam) ;
      exit(1) ;
    }
  }
#endif
  
  /* calculate conserved quantities */
  if(call_code==INIT_OUT ||
     call_code==LOG_OUT ||
     call_code==DIVB_OUT ||
     (call_code==FINAL_OUT &&
     !failed)) {
    pp = 0. ;
    e = 0. ;
    rmed = 0. ;
    divbmax = 0. ;
    imax = 0 ;
    jmax = 0 ;
    kmax = 0;
    di = (N1>1); //boolean loops ; di = 1 if N1 > 1 otherwise 0
    dj = (N2>1); //boolean loops ; dj = 1 if N2 > 1 otherwise 0
    dk = (N3>1); //boolean loops ; dk = 1 if N3 > 1 otherwise 0
    //all cells not immediately adjacent to physical boundaries
    istart = is_physical_bc(1, 0);
    istop  = N1 - 1 - is_physical_bc(1, 1);
    jstart = 0*is_physical_bc(2, 0);
    jstop  = N2 - 0*is_physical_bc(2, 1);
    kstart = 0*is_physical_bc(3, 0);
    kstop  = N3 - 0*is_physical_bc(3, 1);
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) {
      get_geometry(i,j,k,CENT,&geom) ;
      get_state(p[i][j][k],&geom,&q) ;
      primtoU(p[i][j][k],&q,&geom,U) ;
      
      rmed += U[RHO]*dV ;
      pp += U[U3]*dV ;
      e += U[UU]*dV ;

      /* flux-ct defn */
      /* divergence of the magnetic field, divB, approximated by finite differences applied to the cell‐centered primitive magnetic fields multiplied by the metric determinant (gdet). 
      In each coordinate direction, the derivative is approximated using a centered difference‐stencil that averages over four neighboring points in the transverse directions. 
       In each coordinate direction, the stencil uses values from the current cell and the immediately previous cell (i–1 in x1, j–1 in x2, and k–1 in x3) combined with an average over four transverse “corner” points to yield a robust, second-order approximation of the derivative of B−g
      The chosen stencil computes a centered difference. The derivative is approximated as:
              (⟨B1−g⟩_i−⟨B1−g⟩_{i−1})/dx[1]
where the averages ⟨⋅⟩⟨⋅⟩ over the four adjacent cells (at j, j–1, k, k–1) yield a second-order accurate derivative.
    Definition Location:
      The computed divB is a cell-centered quantity since it uses cell-centered values (denoted with [CENT]) and is intended as a diagnostic at the cell center.
    Choice of Stencil (Why Not i+1, etc.):
The stencil uses i and i–1 (rather than i+1) because this backward difference (averaged over transverse directions) provides a centered approximation consistent with the finite-volume methodology. The scheme is set up so that the difference over the interface between cell i and i–1 is representative of the divergence at that cell center. Inserting i+1 instead would change the location and could disrupt the symmetry and second-order accuracy of the derivative.
*/
      divb = fabs(
#if(N1>1)
                  0.25*(
                        + p[i  ][j   ][k   ][B1]*gdet[i  ][j   ][k   ][CENT]
                        + p[i  ][j   ][k-dk][B1]*gdet[i  ][j   ][k-dk][CENT]
                        + p[i  ][j-dj][k   ][B1]*gdet[i  ][j-dj][k   ][CENT]
                        + p[i  ][j-dj][k-dk][B1]*gdet[i  ][j-dj][k-dk][CENT]
                        - p[i-1][j   ][k   ][B1]*gdet[i-1][j   ][k   ][CENT]
                        - p[i-1][j   ][k-dk][B1]*gdet[i-1][j   ][k-dk][CENT]
                        - p[i-1][j-dj][k   ][B1]*gdet[i-1][j-dj][k   ][CENT]
                        - p[i-1][j-dj][k-dk][B1]*gdet[i-1][j-dj][k-dk][CENT]
                        )/dx[1] //This term approximates the radial derivative ∂x1(B1\sqrt{−g}) at the cell interface between cells i-1 and i
#endif
#if(N2>1)
                 +0.25*(
                        + p[i   ][j  ][k   ][B2]*gdet[i   ][j  ][k   ][CENT]
                        + p[i   ][j  ][k-dk][B2]*gdet[i   ][j  ][k-dk][CENT]
                        + p[i-di][j  ][k   ][B2]*gdet[i-di][j  ][k   ][CENT]
                        + p[i-di][j  ][k-dk][B2]*gdet[i-di][j  ][k-dk][CENT]
                        - p[i   ][j-1][k   ][B2]*gdet[i   ][j-1][k   ][CENT]
                        - p[i   ][j-1][k-dk][B2]*gdet[i   ][j-1][k-dk][CENT]
                        - p[i-di][j-1][k   ][B2]*gdet[i-di][j-1][k   ][CENT]
                        - p[i-di][j-1][k-dk][B2]*gdet[i-di][j-1][k-dk][CENT]
                       )/dx[2]
#endif
#if(N3>1)
                 +0.25*(
                        + p[i   ][j   ][k  ][B3]*gdet[i   ][j   ][k  ][CENT]
                        + p[i-di][j   ][k  ][B3]*gdet[i-di][j   ][k  ][CENT]
                        + p[i   ][j-dj][k  ][B3]*gdet[i   ][j-dj][k  ][CENT]
                        + p[i-di][j-dj][k  ][B3]*gdet[i-di][j-dj][k  ][CENT]
                        - p[i   ][j   ][k-1][B3]*gdet[i   ][j   ][k-1][CENT]
                        - p[i-di][j   ][k-1][B3]*gdet[i-di][j   ][k-1][CENT]
                        - p[i   ][j-dj][k-1][B3]*gdet[i   ][j-dj][k-1][CENT]
                        - p[i-di][j-dj][k-1][B3]*gdet[i-di][j-dj][k-1][CENT]
                   )/dx[3]
#endif
                  );
      if(divb > divbmax) {
        imax = i+mpi_startn[1] ;
        jmax = j+mpi_startn[2] ;
        kmax = k+mpi_startn[3] ;
        divbmax = divb ;
      }
    }
  }
  
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(&divbmax,&divbmax_global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&e,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&rmed,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&pp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
  divbmax_global = divbmax;
#endif
  
  if(call_code == INIT_OUT) {
    e_init = e ;
    m_init = rmed ;
  }
  
  if(MASTER == mpi_rank && call_code == FINAL_OUT) {
    e_fin = e ;
    m_fin = rmed ;
    fprintf(stderr,"\n\nEnergy: ini,fin,del: %g %g %g\n",
            e_init,e_fin,(e_fin-e_init)/e_init) ;
    fprintf(stderr,"mass: ini,fin,del: %g %g %g\n",
            m_init,m_fin,(m_fin-m_init)/m_init) ;
  }
  
  if(call_code == INIT_OUT ||
     call_code == LOG_OUT ||
     call_code == FINAL_OUT ||
     call_code == DIVB_OUT) {
    if (divbmax==divbmax_global) {
      fprintf(stderr,"LOG      t=%g \t divbmax: %d %d %d %g\n",
            t,imax,jmax,kmax,divbmax) ;
    }
  }

#if(DOENER)
  if(call_code == INIT_OUT ||
     call_code == LOG_OUT ||
     call_code == FINAL_OUT) {
    fprintf(ener_file,"%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
            t,rmed,pp,e,p[N1/2][N2/2][N3/2][UU]*pow(p[N1/2][N2/2][N3/2][RHO],-gam),
            p[N1/2][N2/2][N3/2][UU]) ;
    fprintf(ener_file,"%15.8g %15.8g %15.8g ",mdot,edot,ldot) ;
    fprintf(ener_file,"\n") ;
    fflush(ener_file) ;
  }
#endif
  
  /* gdump only at code start */
  if(call_code == INIT_OUT) {
    /* make grid dump file */
    
    gdump(0) ;
    gdump2(0) ;
  }
  
  /* dump at regular intervals */
  if(call_code == INIT_OUT ||
     call_code == DUMP_OUT ||
     call_code == FINAL_OUT) {
    /* make regular dump file */
    dump(dump_cnt,1,0) ;
#if(DOPARTICLES)
    pdump(dump_cnt,1);
#endif
    dump_cnt++ ;
  }

  /* output full dump right after restart to check divb being zero */
  if(call_code == DIVB_OUT) {
    /* make regular dump file */
    dump(-dump_cnt,0,0) ;
    //don't do particles for now
    //pdump(-pdump_cnt,0);
    //dump_cnt++ ;
  }
  
  /* rdump at regular intervals */
  if(call_code == INIT_OUT ||
     call_code == RDUMP_OUT ||
     call_code == FINAL_OUT) {
    /* make rdump file */
    //increment rdump count *before* writing the restart dump file
    //so that get the correct counter value if restart from the file
    rdump_cnt++ ;
    restart_write(rdump_cnt-1);
    dump(rdump_cnt-1,0,0);
#if(DOPARTICLES)
      pdump(rdump_cnt-1,0);
#endif
   
  }

  /* image dump at regular intervals */
  if(call_code == IMAGE_OUT ||
     call_code == INIT_OUT ||
     call_code == FINAL_OUT) {
    
    fdump( image_cnt );
    
    image_cnt++ ;
  }
}


/** some diagnostic routines **/


void fail(int fail_type)
{

	failed = 1 ;

	fprintf(stderr,"\n\nfail: %d %d %d %d\n",icurr,jcurr,kcurr,fail_type) ;

	area_map(icurr,jcurr,kcurr,p) ;
	
	fprintf(stderr,"fail\n") ;

	diag(FINAL_OUT) ;

#ifdef MPI
        MPI_Finalize();
#endif
	/* for diagnostic purposes */
	exit(0) ;
}



/* map out region around failure point */
void area_map(int i, int j, int k, double prim[][N2M][N3M][NPR])
{
	int m ;

	fprintf(stderr,"area map\n") ;

	PLOOP {
		fprintf(stderr,"variable %d \n",m) ;
		fprintf(stderr,"i = \t %12d %12d %12d\n",i-1,i,i+1) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j+1,
				prim[i-1][j+1][k][m],
				prim[i][j+1][k][m],
				prim[i+1][j+1][k][m]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j,
				prim[i-1][j][k][m],
				prim[i][j][k][m],
				prim[i+1][j][k][m]) ;
		fprintf(stderr,"j = %d \t %12.5g %12.5g %12.5g\n",
				j-1,
				prim[i-1][j-1][k][m],
				prim[i][j-1][k][m],
				prim[i+1][j-1][k][m]) ;
	}

	/* print out other diagnostics here */

}

/* evaluate fluxed based diagnostics: mdot, edot, jdot only; put results in
 * global variables */
void diag_flux(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR])
{
	int j,k ;

        mdot = edot = ldot = 0. ;
//#pragma omp parallel for schedule(static,MY_MAX(N2*N3/nthreads,1)) collapse(2) \
//    default(none) \
//    shared(F1,F2,dx,nthreads) \
//    private(j,k) \
//    reduction(+:mdot) reduction(-:edot) reduction(+:ldot)
        for(j=0;j<N2;j++) {
          for(k=0;k<N3;k++) {
                mdot += F1[0][j][k][RHO]*dx[2]*dx[3] ;
                edot -= (F1[0][j][k][UU] - F1[0][j][k][RHO])*dx[2]*dx[3] ;
                ldot += F1[0][j][k][U3] *dx[2]*dx[3] ;
          }
        }
}

/*divb calculation in detail:
Even though all values p[i][j][k] are cell‐centered, when we compute the divergence we are really calculating the net flux through each cell’s boundaries (its “faces”). 
In a finite–volume method, the divergence in a cell is found by summing the fluxes through its faces and then dividing by the cell volume. 
Because our variables are stored at cell centers, we must “reconstruct” or “approximate” what the value would be at the face.
What “Corner Values” Means in This Context

    The Cell Face Is a 2D Surface:
    The face between cells i−1 and i (in the x₁ direction) extends in the other two directions (x₂ and x₃). In other words, the face is not a single point but a surface.

    Averaging over Transverse Cells:
    To approximate the flux at that face using cell–centered data, the scheme averages the values from several cell centers that are “adjacent” to the face. In the code, they take four cell–centered values from the “i” side:
    p[i][j][k],p[i][j][k−dk],p[i][j−dj][k],p[i][j−dj][k−dk]

    The factor 0.25 multiplies these four contributions to give their average. This average is intended to represent the value at the face on the “i” side.
  Why four corner points instead of eight?
  In a three-dimensional structured grid, each cell face is two‐dimensional and—geometrically speaking—has four corners. 
  In other words, when you’re trying to approximate the value at the face (which is an extended 2D surface), 
  you need to sample (or “average”) the four cell‐centered values that are adjacent to that face.

  Using four values is the natural choice because:
    A cell has six faces, and each face is a flat (or approximately flat) quadrilateral.
    Each quadrilateral has four vertices (or corners).
    By averaging the four cell‐centered values that are “around” a particular face, you obtain a good second–order accurate approximation of the flux across that face.

    Similarly for the “i–1” Side:
    They also average the corresponding four cell–centered values from the cell immediately to the left (i.e. cell i−1i−1).

    Computing the Difference:
    By subtracting the average from the i–1 side from the average on the i side and then dividing by dx[1], 
    the code approximates the derivative at the face (which is then associated with the center of the cell). 
    This is a standard centered finite difference—except that here you average over four “corner” (transverse) cell centers to approximate the value at the extended two–dimensional face.

In Summary

    Even Though Data Are Cell Centered:
    The divergence is defined by flux differences over the cell faces. The “corner values” are simply the four cell–centered values used to approximate the average state at the face. 
    They do not mean that divB is defined at the geometric corners of the cell—instead, they are used to recover a properly averaged face value while keeping the divergence itself as a cell–centered diagnostic.

    Why Not Use Directly One Cell Center?
    Using only one cell center (say, j = 0 on the “i” side) would be less accurate because the cell face spans an area in the transverse directions. Averaging four cell–centered values 
    (which come from the four adjacent cells that share the face) gives a higher–order, more robust approximation.

So, although the variables p[i][j][k] are cell centered, the finite-volume method reconstructs “face values” by averaging values from the cells that border that face. The difference between the face value as approximated from the ii-side and that from the i−1i−1-side, divided by the spacing, yields a centered derivative. This derivative is then associated with the cell center.


*/