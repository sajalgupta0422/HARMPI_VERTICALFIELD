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

/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables
 *
 **/


#include "decs.h"
#include "electron.h"
#include "nuclear.h"

/** algorithmic choices **/

/* use local lax-friedrichs or HLL flux:  these are relative weights on each numerical flux */
#define HLLF  (0.0)
#define LAXF  (1.0)

/** end algorithmic choices **/


double advance( double pi[][N2M][N3M][NPR], double pb[][N2M][N3M][NPR],
               double Dt, double pf[][N2M][N3M][NPR], double *ndt1, double *ndt2, double *ndt3) ;
void advance_vpot( double vpot[][N1+D1][N2+D2][N3+D3],
                  double emf[][N1+D1][N2+D2][N3+D3],
                  double Dt);

double fluxcalc( double pr[][N2M][N3M][NPR], double F[][N2M][N3M][NPR],
                int dir ) ;
void   flux_cd(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR]) ;
void   Kick(void);
void   DELVFIXUP(double Ein);
int    ShockDetect(int i, int j, int k);
double del2calc(int i, int j, int k, int m, int dir, int itemp);
double etacalc(int i, int j, int k, int m, int dir, int itemp);



/***********************************************************************************************/
/***********************************************************************************************
 step_ch():
 ---------
 -- handles the sequence of making the time step, the fixup of unphysical values,
 and the setting of boundary conditions;
 
 -- also sets the dynamically changing time step size;

This is a standard second-order accurate predictor-corrector scheme where:
    Variables are first advanced to t+dt/2 (predictor)
    These intermediate values are used to compute accurate fluxes using advance() function
    Variables are then advanced to t+dt (corrector)
    Various fixes ensure physical validity of the solution
 
 ***********************************************************************************************/
void step_ch(double *ndt1, double *ndt2, double *ndt3)
{
  double ndt;
    int i,j,k,m;
#ifdef MPI
  double mpi_buf[3]; // Buffer for MPI communication
#endif
  
    
  if(MASTER==mpi_rank) fprintf(stderr,"h") ; // Print 'h' to indicate half-step
  /*Following line advances primitive variables (p) to intermediate state ( ph) by half timestep (0.5*dt)*/
  advance(p, p, 0.5*dt, ph, ndt1, ndt2, ndt3) ;   /* time step primitive variables to the half step; Advance from p to ph (half-step) */
  
    
	fixup(ph) ;         /* Set updated densities to floor, set limit for gamma */
  fixupuel(ph);       // Fix electron internal energy
	bound_prim(ph) ;    /* Set boundary conditions for primitive variables, flag bad ghost zones */
	fixup_utoprim(ph);  /* Fix the failure points using interpolation and updated ghost zone values, i.e. Fix failures in conservative-to-primitive conversion */
    
        //ZLOOP eVol(p,p,ph,0.5*dt,i,j,k,0,1); /* Account for floor heat */
	bound_prim(ph) ;    /* Reset boundary conditions with fixed up points */
    
#if(DOPARTICLES)
    /* step Lagrangian tracer particles forward */
    if(myNp > 0)
        advance_particles(ph, dt);
#endif
    
  /* Repeat and rinse for the full time (aka corrector) step:  */
  
  if(MASTER==mpi_rank) fprintf(stderr,"f") ; // Print 'f' to indicate full-step
  /* Save current state */
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) shared(psave,p,nthreads) private(i,j,k,m)
  ZLOOP PLOOP psave[i][j][k][m] = p[i][j][k][m] ;
  /* Advance from ph to p (full-step) */
  advance(p, ph, dt,    p, ndt1, ndt2, ndt3) ;

#if(EVOLVEVPOT)
  //in-place update of vector potential
  //do so only on full step since no need to know vpot on half-steps
  advance_vpot(vpot, emf, dt);
#endif
  
  //Post Full-Step Corrections:
  fixup(p) ;
  fixupuel(p);
  bound_prim(p) ;
  fixup_utoprim(p);
  //ZLOOP eVol(psave,ph,p,dt,i,j,k,0,1);  /* Account for floor heat */

  bound_prim(p) ;
    
  //Turbulence Driving (if applicable)
  if(WHICHPROBLEM == TURB){
    Kick(); // Add turbulent forcing
    bound_prim(p) ;
  }
  
  
  /* Determine next time increment based on current characteristic speeds: */
  if(dt < 1.e-9 && MASTER==mpi_rank) {
    fprintf(stderr,"timestep too small: dt = %g\n", dt) ;
    exit(11) ;
  }
  
  
  /* increment time */
  
  t += dt ; 
  
  /* set next timestep */
#ifdef MPI
  //find the minimum ndt across all MPI processes
  mpi_buf[0] = *ndt1;
  mpi_buf[1] = *ndt2;
  mpi_buf[2] = *ndt3;
  MPI_Allreduce(MPI_IN_PLACE,mpi_buf,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
  *ndt1 = mpi_buf[0];
  *ndt2 = mpi_buf[1];
  *ndt3 = mpi_buf[2];
#endif
// Harmonic mean of directional timesteps with safety factor
  ndt = defcon * 1./(1./(*ndt1) + 1./(*ndt2) + 1./(*ndt3)) ;
  if(ndt > SAFE*dt) ndt = SAFE*dt ;
  dt = ndt ; //dt changes
  if(t + dt > tf) dt = tf - t ;  /* but don't step beyond end of run */
  
  //Final MPI Synchronization:
#ifdef MPI
  mpi_buf[0] = t;
  mpi_buf[1] = dt;
  MPI_Bcast(mpi_buf,2,MPI_DOUBLE,MASTER,MPI_COMM_WORLD);
  t = mpi_buf[0];
  dt = mpi_buf[1];
#endif
  
  /* done! */
}



/***********************************************************************************************/
/***********************************************************************************************
 advance():
 ---------
 -- responsible for what happens during a time step update, including the 
      flux calculation,
      the constrained transport calculation (aka flux_ct(),i.e div(B) = 0), 
      the finite difference form of the time integral, 
      and the calculation of the primitive variables from the update conserved variables;

 -- also handles the "fix_flux()" call that sets the boundary condition on the fluxes;
 
 ***********************************************************************************************/
double advance(
               double pi[][N2M][N3M][NPR], // initial state
               double pb[][N2M][N3M][NPR], // base state for flux calculation
               double Dt, // timestep     
               double pf[][N2M][N3M][NPR], // final state
               /*timestep constraints*/
               double *ndt1,
               double *ndt2,
               double *ndt3
               )
{ 
  /*
  Declares necessary variables including arrays for conserved variables (U), which are {\sqrt{-g}* (ρu^t, T^t_t , T^t_i , B^i)}, 
  their changes ( dU), and geometric structures.
  */
  int i,j,k,m, dummy;
  double ndt,U[NPR],dU[NPR], usave;
  struct of_geom geom ;
  struct of_state q ; //for four-velocities and magnetic fields (both covariant and contravariant)
  double kappa;
  
  double Uo[NPR],po[NPR] ;
  
  int was_floor_activated = 0;
  
  //Copies initial state to final state array in parallel.
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) shared(pf,pi,nthreads) private(i,j,k,m)
  ZLOOP PLOOP pf[i][j][k][m] = pi[i][j][k][m] ;        /* needed for Utoprim */
  
  if(MASTER==mpi_rank) fprintf(stderr,"0") ;

  //Calculates fluxes in all three directions using the base state pb. Returns timestep constraints
  *ndt1 = fluxcalc(pb, F1, 1) ;
  *ndt2 = fluxcalc(pb, F2, 2) ;
  *ndt3 = fluxcalc(pb, F3, 3) ;
  
  fix_flux(F1,F2,F3) ; //defined in coord.c; Applies boundary conditions to fluxes
  
  flux_ct(F1,F2,F3) ; //performs constrained transport to maintain div(B) = 0.
  
  /* evaluate diagnostics based on fluxes */
  diag_flux(F1,F2) ; //in diag.c
  
  if(MASTER==mpi_rank) fprintf(stderr,"1") ;
  /** now update pi to pf **/
  
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) \
    shared(pb,pi,pf,pflag,F1,F2,F3,dx,nthreads,failimage,Katm,Dt) \
    private(i,j,k,m,q,dU,geom,U,usave,was_floor_activated,kappa)
  ZLOOP {
    /*
    For each cell:
      Gets geometric factors
      Calculates source terms
      Gets state variables
    */
    get_geometry(i,j,k,CENT,&geom) ;

    source(pb[i][j][k],&geom,i,j,k,dU,Dt) ; //in phys.c
    get_state(pb[i][j][k],&geom,&q) ;
    //misc_source(pb[i][j][k],pb[i+1][j][k],pb[i-1][j][k], pb[i][j+1][k],pb[i][j-1][k],i, j, k, &geom,&q,dU, Dt);

    get_state(pi[i][j][k],&geom,&q) ;
    primtoU(pi[i][j][k],&q,&geom,U) ;

      
    PLOOP {
      //Updates conserved variables using flux differences and source terms.
      U[m] += Dt*(
#if( N1 != 1 )
		  - (F1[i+1][j][k][m] - F1[i][j][k][m])/dx[1] //\partial F^r / \partial x1
#endif
#if( N2 != 1 )
		  - (F2[i][j+1][k][m] - F2[i][j][k][m])/dx[2] //\partial F^θ / \partial x2
#endif
#if( N3 != 1 )
		  - (F3[i][j][k+1][m] - F3[i][j][k][m])/dx[3] //\partial F^φ / \partial x3
#endif
		  + dU[m]

		  ) ;

    }

    
    //Converts updated conserved variables back to primitive variables
    //and stores them in pf (the states at timesetp Dt).
    /*
    This Utoprim_2d() function converts the updated conserved variables back to primitive variables and stores them in  pf. So the sequence is:
      Copy initial state: pf = pi (done earlier)
      Convert primitives to conserved: pi -> U (done earlier)
      Update conserved variables using fluxes and sources (done earlier)
      Convert updated conserved back to primitives: U -> pf (doing below)
    */
    pflag[i][j][k] = Utoprim_2d(U, geom.gcov, geom.gcon, geom.g, pf[i][j][k]);
    
    
    if( pflag[i][j][k] ) {
      failimage[i][j][k][0]++ ;
      was_floor_activated = 1;
    }
    else {
      was_floor_activated = 0;
    }

    
#if( DO_FONT_FIX || eHEAT || DOKTOT )
    if (eHEAT || DOKTOT) {
      kappa = pf[i][j][k][KTOT];
    }
    else {
      kappa = Katm[i];
    }
    if( pflag[i][j][k] ) {
      
      pflag[i][j][k] = Utoprim_1dvsq2fix1(U, geom.gcov, geom.gcon, geom.g, pf[i][j][k], kappa );
      if( pflag[i][j][k] ) {
        was_floor_activated = 2;
	failimage[i][j][k][1]++ ;
	pflag[i][j][k] = Utoprim_1dfix1(U, geom.gcov, geom.gcon, geom.g, pf[i][j][k], kappa);
        if( pflag[i][j][k] ) {
          was_floor_activated = 3;
          failimage[i][j][k][2]++ ;
        }
      }
    }
#if(DOFLR)
    if( was_floor_activated ) {
      pf[i][j][k][FLR] = 1.;
    }
#endif
#endif
    
#if(eHEAT || eCOND)
    eVol(pi,pb,pf,Dt,i,j,k,was_floor_activated,0); //ELECTRONS
#endif
    
#if(DONUCLEAR)
    nuc_evol(pi,pb,pf,Dt,i,j,k,was_floor_activated); //nuclear physics
#endif
    
  }
  
    
  
  
  if(MASTER==mpi_rank)  fprintf(stderr,"2") ;
  
  return(0) ;
}

//in-place update of vector potential, vpot, according to the electric fields, emf
//both vpot and emf are located at cell corners
//update ONLY on the full time step since that's what interested in
void advance_vpot(
               double vpot[][N1+D1][N2+D2][N3+D3],
               double emf[][N1+D1][N2+D2][N3+D3],
               double Dt
               )
{
    int i, j, k;
    #pragma omp parallel for schedule(static,(N1+D1)*(N2+D2)*(N3+D3)/nthreads) collapse(3) default(none) shared(vpot,emf,nthreads,Dt) private(i,j,k)
    //Loop for cell-cornered variables
    //E⃗ = ∂A⃗/∂t −∇V; but V = 0 //Faraday's law
    ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
        vpot[1][i][j][k] += emf[1][i][j][k] * Dt;
        vpot[2][i][j][k] += emf[2][i][j][k] * Dt;
        vpot[3][i][j][k] += emf[3][i][j][k] * Dt;
    }
}

/***********************************************************************************************/
/***********************************************************************************************
 fluxcalc():
 ---------
 -- sets the numerical fluxes, avaluated at the cell boundaries using the slope limiter
 slope_lim();
 
 -- only has HLL and Lax-Friedrichs  approximate Riemann solvers implemented;
 
 The fluxcalc() function calculates numerical fluxes at cell boundaries using either HLL (Harten-Lax-van Leer) or
  Lax-Friedrichs approximate Riemann solvers
The function also returns a double (the local timestep, ndt) and computes the fluxes F in the chosen direction (dir).
 ***********************************************************************************************/
double fluxcalc(
                double pr[][N2M][N3M][NPR], //Array of primitive variables
                double F[][N2M][N3M][NPR], //Array to store computed fluxes
                int dir //Direction (1, 2, or 3) for flux calculation
                )
{
  int i,j,k,m,idel,jdel,kdel,face ; //Various local loop counters (i, j, k, m) and helper variables are declared.
  double p_l[NPR],p_r[NPR],F_l[NPR],F_r[NPR],U_l[NPR],U_r[NPR] ; //Arrays (like p_l[], p_r[], etc.) store reconstructed left and right states
  double cmax_l,cmax_r,cmin_l,cmin_r,cmax,cmin,ndt,dtij ;
  double ctop, s_l, s_r ;
  double cmax_lw,cmax_rw,cmin_lw,cmin_rw,cmaxw,cminw ;
  double ctopw, dummy;
  double  tmp;
  
  struct of_geom geom ;
  struct of_state state_l,state_r ; //State structs (state_l, state_r) store the result of converting primitives to conserved or flux–related state.
  struct of_state state_lw,state_rw ;
  
  void rescale(double *pr, int which, int dir, int ii, int jj, int kk, int face, struct of_geom *geom) ;
  double bsq ;
  
  if( (N1==1 && dir == 1) || (N2==1 && dir == 2) || (N3==1 && dir == 3)) {
    //if current direction is trivial, reset fluxes in the current direction to zero so no evolution due to this direction
#pragma omp parallel for schedule(static,(N1+2*D1)*(N2+2*D2)*(N3+2*D3)/nthreads) collapse(3) default(none) shared(dq,F,nthreads) private(i,j,k,m)
    ZSLOOP(-D1,N1-1+D1,-D2,N2-1+D2,-D3,N3-1+D3) PLOOP {
      dq[i][j][k][m] = 0.; //slopes are resetting to 0 if N1 =N2 = N3 =1
      F[i][j][k][m] = 0.; //fluxes are setting to 0
    }
    ndt = 1.e9; //returning very large timestep so that no restrive CFL condition comes from a trivial dimension
    return(ndt) ;
  }

  //idel, jdel, kdel to indicate the shift vector used for reconstructing values from neighboring cells
  //face indicates which face (of the cell) the flux will be evaluated at (e.g., FACE1 is the radial face).
  if     (dir == 1) {idel = 1; jdel = 0; kdel = 0; face = FACE1;} //idel = 1 means “look one cell ahead in the x1 direction.”;jdel and kdel are zero because the reconstruction is only along the x1 axis.
  else if(dir == 2) {idel = 0; jdel = 1; kdel = 0; face = FACE2;}
  else if(dir == 3) {idel = 0; jdel = 0; kdel = 1; face = FACE3;}
  else { exit(10); }
  
#if(RESCALE)
  /** evaluate slopes of primitive variables **/
  /* first rescale */
#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) default(none) collapse(3) shared(pr,dir,nthreads) private(i,j,k,m,geom)
  ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G-1) {
    get_geometry(i,j,k,CENT,&geom) ;
    rescale(pr[i][j][k],FORWARD, dir, i,j,k,CENT,&geom) ;
  }
#endif

  /* then evaluate slopes on active grid plus 1-cell-wide layer of ghost cells */
#pragma omp parallel for schedule(static,(N1+2*D1)*(N2+2*D2)*(N3+2*D3)/nthreads) collapse(3) default(none) shared(dq,pr,idel,jdel,kdel,nthreads) private(i,j,k,m)
 ZSLOOP(-D1,N1-1+D1,-D2,N2-1+D2,-D3,N3-1+D3){
    PLOOP {
        //Purpose:This loop computes “slopes” for each primitive variable using a chosen slope limiter (slope_lim).
        //Mechanism: For each cell, it uses the values from the neighboring cells only for the passed direction (shifted by ± idel, etc.) 
        //to estimate the gradient.
        //Role: These slopes (stored in array dq) are used later for reconstruction of left and right states at cell faces.
      //-new get_geometry(i,j,CENT,&geom) ;
      //-new bsq = bsq_calc(pr[i][j],&geom) ;
      //-new if(bsq/pr[i][j][RHO] > 10. ||
      //-new    bsq/pr[i][j][UU]  > 1.e3) lim = MINM ;
      //-new else lim = MC ;
      
      dq[i][j][k][m] = slope_lim(
                              pr[i-idel][j-jdel][k-kdel][m],
                              pr[i][j][k][m],
                              pr[i+idel][j+jdel][k+kdel][m]
                              ) ;
    }
  
//!!!ATCH: took the following out of PLOOP since does not depend on the primitive
//    dsource[i][j][k] = slope_lim(
//                              qrat[i-idel][j-jdel][k-kdel],
//                              qrat[i][j][k],
//                              qrat[i+idel][j+jdel][k+kdel]
//                              ) ;
  }
  
  
  ndt = 1.e9 ;
#pragma omp parallel for schedule(static,(D1*(N1+!idel)+1)*(D2*(N2+!jdel)+1)*(D3*(N3+!kdel)+1)/nthreads) \
    collapse(3) \
    shared(F,dq,idel,jdel,kdel,nthreads,dir,pr,face,dx,cour,mpi_startn,mpi_ntot) \
    private(i,j,k,m,geom,p_l,p_r,state_l,state_r,F_l,F_r,U_l,U_r,cmin_l,cmin_r,cmin,cmax,ctop,dtij,tmp,cmax_l,cmax_r) \
    default(none) \
    reduction(min:ndt)
  ZSLOOP(-D1*!idel,D1*N1,-D2*!jdel,D2*N2,-D3*!kdel,D3*N3){ //-jdel,N1,-idel,N2) {
    
    //if dir=1, idel=1,jdel = 0,kdel = 0, then ZSLOOP(0,N1,-1,N2,-1,N3)
    
    /* this avoids problems on the pole 
    Special Handling for the θ Direction Near the Pole 
    When computing fluxes in the θ (x2) direction on a face that lies at the polar axis (either j = 0 or j = mpi_ntot[2] in global coordinates) 
    and if Boyer-Lindquist coordinates (BL) are in use, the fluxes are set to zero to avoid numerical difficulties at the pole.
    
    The rationale for setting the flux to zero near the pole (θ = 0 or θ = π) in reflective schemes is to avoid numerical 
    singularities due to the vanishing of sin⁡θ and degeneracies in the grid.
    This “flux zeroing” forces the intercell flux to vanish, effectively decoupling the polar ghost zones from the interior 
    to prevent unphysical inflow or spurious currents due to the coordinate singularity.
    */
    if(dir == 2 && (j+mpi_startn[2] == 0 || j+mpi_startn[2] == mpi_ntot[2]) && BL) {
      //FIXFLUXMARK: fix_flux-like zero out of fluxes in x2-dir , putting flux = 0 near poles -- change it for transmissive BC
      PLOOP F[i][j][k][m] = 0. ;
    }
    else {
      
#if(RECONSTRUCT == LIN) //Linear reconstruction
      PLOOP {
        //p_l​ (the left state at the cell face) is obtained by taking the value in the cell just “upwind” (at i–idel, etc.) 
        //and then adding half of the slope (from dq).
        p_l[m] = pr[i-idel][j-jdel][k-kdel][m]
        + 0.5*dq[i-idel][j-jdel][k-kdel][m] ; 
        //Similarly, p_r​ (the right state) is computed from the cell at (i, j, k) by subtracting half the slope.
        p_r[m] = pr[i][j][k][m]
        - 0.5*dq[i][j][k][m] ;
      }
#elif (RECONSTRUCT ==PARA) //Parabolic reconstruction using  a five-point stencil

          PLOOP
          {
            double x1 = pr[i-3*idel][j-3*jdel][k-3*kdel][m];
            double x2 = pr[i-2*idel][j-2*jdel][k-2*kdel][m];
            double x3 = pr[i-1*idel][j-1*jdel][k-1*kdel][m];
            double x4 = pr[i][j][k][m];
            double x5 = pr[i+1*idel][j+1*jdel][k+1*kdel][m];
            
            para(x1, x2, x3, x4, x5, &tmp, &p_l[m]);
            
            x1 = pr[i-2*idel][j-2*jdel][k-2*kdel][m];
            x2 = pr[i-1*idel][j-1*jdel][k-1*kdel][m];
            x3 = pr[i][j][k][m];
            x4 = pr[i+1*idel][j+1*jdel][k+1*kdel][m];
            x5 = pr[i+2*idel][j+2*jdel][k+2*kdel][m];
            
            para(x1, x2, x3, x4, x5, &p_r[m], &tmp);
          }
#endif
        
      
      
      
      
      
      
      get_geometry(i,j,k,face,&geom) ;//The geometry is obtained at the cell face (rather than the center) using the “face” index determined earlier.
      
#if(RESCALE)
      rescale(p_l,REVERSE,dir,i,j,k,face,&geom) ;
      rescale(p_r,REVERSE,dir,i,j,k,face,&geom) ;
#endif
      // Converts primitive variables (p_l and p_r) into a state structure that typically contains conserved variables or other derived quantities.
      get_state(p_l,&geom,&state_l) ;
      get_state(p_r,&geom,&state_r) ;
      
      //calculate fluxes in direction dir and store in F_l and F_r, The fluxes are defined at zone faces.
      primtoflux(p_l,&state_l,dir,&geom,F_l) ;
      primtoflux(p_r,&state_r,dir,&geom,F_r) ;
      
      //calculate fluxes in time-direction and store in U_l and U_r, The fluxes are defined at zone faces.
      primtoflux(p_l,&state_l,TT, &geom,U_l) ;
      primtoflux(p_r,&state_r,TT, &geom,U_r) ;
      
      //calculating characteristic wave speed for courant or CFL condition 
      //(It returns maximum and minimum speeds (cmax_l, cmin_l for the left state and similarly for the right).
      vchar(p_l,&state_l,&geom,dir,&cmax_l,&cmin_l) ;
      vchar(p_r,&state_r,&geom,dir,&cmax_r,&cmin_r) ;
      
      cmax = fabs(MY_MAX(MY_MAX(0., cmax_l),  cmax_r)) ;
      cmin = fabs(MY_MAX(MY_MAX(0.,-cmin_l), -cmin_r)) ;
      ctop = MY_MAX(cmax,cmin) ;
      
      PLOOP {
        //A hybrid Riemann solver is then applied:
        //The HLLF flux is computed using an HLL type formula. (weightage = 0); a purely HLL-type solver might be too diffusive
        //A Lax–Friedrichs flux is also computed. (weightage = 1), and can improve stability

        /*Mathematical rationale:
         1.First, Averages the physical fluxes from left and right states (0.5*(F_l[m] + F_r[m]))
         2.Them Subtracts numerical dissipation (-ctop*(U_r[m] - U_l[m])) proportional to: 
              The maximum wave speed, ctop
              The jump in conserved variables (U_r[m] - U_l[m])
        The subtraction of ctop(Ur−Ul)) acts as an artificial dissipation term, ensuring that waves traveling through the interface do so stably.
        The dissipation term stabilizes the scheme by damping high-frequency oscillations, 
        but can make the solution more diffusive 
        
        Physical Rationale:
          The formula represents an upwinded flux in which information from the faster-moving waves dominates the interaction at the interface.
          The “dissipation” added by the Lax–Friedrichs flux is consistent with physical diffusion that might arise in highly nonlinear regimes; 
          however, its main purpose here is numerical stability.
          The scheme approximates the true intercell flux while ensuring that the numerical method remains robust even in the presence of steep gradients or shocks.
        */
          F[i][j][k][m] =
          HLLF*(
                (cmax*F_l[m] + cmin*F_r[m]
                 - cmax*cmin*(U_r[m] - U_l[m]))/
                (cmax + cmin + SMALL)
                ) +
          LAXF*(
                0.5*(F_l[m] + F_r[m]
                     - ctop*(U_r[m] - U_l[m]))
                
                ) ;
        }
        
//                        }
//                    }
      
      /* evaluate restriction on timestep */
      cmax = MY_MAX(cmax,cmin) ;
      dtij = cour*dx[dir]/cmax ; //For each interface, the code computes a local timestep restriction dtij=cour⋅dx[dir]/cmax according to the Courant condition.
      if(dtij < ndt) ndt = dtij ;
        
    }
  }
  
#if(RESCALE)
#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) default(none) collapse(3) shared(pr,dir,nthreads) private(i,j,k,m,geom)
  ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G) {
    get_geometry(i,j,k,CENT,&geom) ;
    rescale(pr[i][j][k],REVERSE,dir,i,j,k,CENT,&geom) ;
  }
#endif
  
  
  return(ndt) ;
  
}

/***********************************************************************************************/
/***********************************************************************************************
 flux_ct():
 ---------
 -- performs the flux-averaging used to preserve the del.B = 0 constraint (see Toth 2000); by doing:
(1) computing the electromotive forces (EMFs) at cell edges from face‐fluxes and
 (2) “rewriting” these EMFs into updated fluxes on cell faces for each magnetic field component.
 ***********************************************************************************************/
void flux_ct(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR], double F3[][N2M][N3M][NPR])
{
  int i,j,k ;

  /* calculate EMFs */
#define DOE1 (N2>1 && N3>1) //boolen: 1 if N2>1 & N3>1
#define DOE2 (N1>1 && N3>1) //boolen: 1 if N1>1 & N3>1
#define DOE3 (N1>1 && N2>1) //boolen: 1 if N1>1 & N2>1
  
  /* Toth approach: just average */
#pragma omp parallel for schedule(static,(N1+D1)*(N2+D2)*(N3+D3)/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {

//Rationale for: E_i = e_{ijk} F_j(B_k) 
/*
In ideal MHD the electric field EE (or EMF) is typically given by a cross‐product, E=−v×B ,
but in constrained transport (CT) schemes it is often more useful to work directly with the electromotive force (EMF) 
defined at the cell edges using the antisymmetric permutation tensor ϵijk​. In this context a generic relation is written as
Ei=ϵijk Fj(Bk) ,
*/
#if(DOE1)
    //E1
    /*
    The fluxes in the x₂ direction (F2) contribute the terms involving B3.
    The fluxes in the x₃ direction (F3) contribute the terms involving B2.
    The sum:
      1/4[(F2 at (i,j,k)+F2 at (i,j,k−1))−(F3 at (i,j,k)+F3 at (i,j−1,k))]
      represents an average over four adjacent faces meeting at the cell edge where E₁ is defined.

    The factor 0.25 averages over the four contributions 
    (the face values are not points but averages over the cell face’s “corners” in the transverse directions).
    */
    emf[1][i][j][k] = 0.25*((F2[i][j][k][B3] + F2[i][j][k-1][B3])
                          - (F3[i][j][k][B2] + F3[i][j-1][k][B2])
                            ) ;
#endif

#if(DOE2)
    //E2; Again, a four-point average is taken over the values at adjacent locations in the i and k directions.
    emf[2][i][j][k] = 0.25*((F3[i][j][k][B1] + F3[i-1][j][k][B1])
                          - (F1[i][j][k][B3] + F1[i][j][k-1][B3])
                            ) ;
#endif

#if(DOE3)
    //E3; The contributions come from the x₁ flux (F1 with B2) and the x₂ flux (F2 with B1).
    emf[3][i][j][k] = 0.25*((F1[i][j][k][B2] + F1[i][j-1][k][B2])
                          - (F2[i][j][k][B1] + F2[i-1][j][k][B1])
                           ) ;
#endif
  }

  /* rewrite EMFs as fluxes, after Toth */
#pragma omp parallel for schedule(static,(N1+D1)*N2*N3/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1+D1,0,N2-1,0,N3-1) {

/*
Why Averages?
These averages are taken because the E-field is defined at the cell edge, but the fluxes need to be defined at the cell face. 
The averaging over two adjacent edge EMFs yields a face-centered flux that is consistent with the discrete curl representation of Faraday’s law.
*/
#if(N1>1)
    F1[i][j][k][B1] = 0. ; //B1 flux is set to 0 (there’s no update in that component from the CT step at an x₁ face).
#endif
#if(DOE3)
    F1[i][j][k][B2] =  0.5*(emf[3][i][j][k] + emf[3][i][j+1][k]) ; //B2 flux is assigned as the average of two E₃ values (from j and j+1).
#endif
#if(DOE2)
    F1[i][j][k][B3] = -0.5*(emf[2][i][j][k] + emf[2][i][j][k+1]) ; //B3 flux is set as negative one-half the average of two E₂ values (from k and k+1).
#endif
  }

#pragma omp parallel for schedule(static,N1*(N2+D2)*N3/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1,0,N2-1+D2,0,N3-1) {
#if(DOE3)
    F2[i][j][k][B1] = -0.5*(emf[3][i][j][k] + emf[3][i+1][j][k]) ; //B1 flux now becomes the negative half–average of two adjacent E₃ values in the x₁ direction.
#endif
#if(DOE1)
    F2[i][j][k][B3] =  0.5*(emf[1][i][j][k] + emf[1][i][j][k+1]) ; //B3 flux is the half–average of two adjacent E₁ values in the x₃ direction
#endif
#if(N2>1)
    F2[i][j][k][B2] = 0. ; //B2 is set to zero on x₂ faces.
#endif
  }

#pragma omp parallel for schedule(static,N1*N2*(N3+D3)/nthreads) collapse(3) default(none) shared(emf,F1,F2,F3,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1+D3) {
#if(DOE2)
    F3[i][j][k][B1] =  0.5*(emf[2][i][j][k] + emf[2][i+1][j][k]) ; //B1 flux is given by the half–average of two E₂ values along the x₁ direction.
#endif
#if(DOE1)
    F3[i][j][k][B2] = -0.5*(emf[1][i][j][k] + emf[1][i][j+1][k]) ; //B2 flux is the negative half–average of two E₁ values along the x₂ direction.
#endif
#if(N3>1)
    F3[i][j][k][B3] = 0. ; //B3 flux is set to 0 on x₃ faces.
#endif
  }
}
/***********************************************************************************************/
/***********************************************************************************************
 Kick():
 ---------
 -- Gives velocity kicks at each time step ;
 --implement turbulent velocity perturbations in Fourier space:
 ***********************************************************************************************/
void Kick()
{
  //Declares indices and variables for wave numbers, velocities, and energies
  int i,j,k ;
  double  delvmag;
    double cs, kx, ky, dk, kmag, kpeak,sigv;
  double u1r, u2r, zr, eps;
  double vxavg, vyavg;
    double Ekick, Ein, norm;
  
  cs =sqrt(gam * (gam-1) * (gam-1.)/pow(1000.,2));
  
  //Set up velocity purturbations in k-space (Fourier space)//
  dk = 2.*M_PI/1.; // Wave number spacing
    kpeak = dk*2.; // Peak wave number
    eps = .5 ; //.1;  //fraction of energy below sound speed
    Ein = eps*cs*cs*cs*dt; // Input energy scaling
  
    //SAS3D: need to make 3D (so far set k = 0)
    ZSLOOP(0,N1-1,0,N2-1,0,0) {
        kx = (1.*i)*dk; // x-component of wave vector
        ky = (1.*j)*dk; // y-component of wave vector
    kmag = sqrt(kx*kx+ky*ky); // Wave vector magnitude
    
    if (kmag<M_PI*N1){ // Check if within Nyquist frequency
            sigv = sqrt(pow(kmag,6.)*exp(-8.*kmag/kpeak)); // Velocity spectrum
      
            
            //Implements Box-Muller transform to generate Gaussian random numbers
      u1r = ranc(0); // Random number 1
      u2r = ranc(0); // Random number 2
      zr = sin(2.*M_PI*u1r)*sqrt(-2.*log(u2r)); //// Box-Muller transform
      delvmag = sigv*zr; // Random velocity magnitude
      
            
            //delv dot k = 0
            if(kmag>=dk){
      delvx[i][j] = -ky/kmag * delvmag; // x-component
      delvy[i][j] = kx/kmag * delvmag; // y-component
                
    }
    else{
                delvx[i][j] = 0.+I*0.;
                delvy[i][j] = 0.+I*0.;
    }
        }
        else{
            delvx[i][j] = 0.+I*0.; // Zero for k=0 mode
            delvy[i][j] = 0.+I*0.;
            
        }
        
       
        

    
  }

  //Fourier transform to real space;  Converts velocity perturbations from k-space to real space
  i = FFT2D(delvx,N1,N2,-1); // Transform x-velocities
  j = FFT2D(delvy,N1,N2,-1); // Transform y-velocities

    DELVFIXUP(Ein);  //Normalize velocity kicks, subtracts off average, and makes real
  
    //SAS3D: note, that here simple copy in the 3rd direction
    ZSLOOP(0,N1-1,0,N2-1,0,N3-1){
        p[i][j][k][U1] += creal(delvx[i][j]);
        p[i][j][k][U2] += creal(delvy[i][j]);
    }
  

    
  
  }

/***********************************************************************************************/
/***********************************************************************************************
 delvfixup():
 ---------
 -- Gives velocity kicks at each time step ;
-- The function serves three main purposes:
      Removes mean flow (ensures zero net momentum addition)
      Makes velocities real (removes imaginary components)
      Normalizes velocities to achieve desired energy input
 ***********************************************************************************************/

void DELVFIXUP(double Ein){
    
    double vxavg,vyavg, Ekick;
    double norm;
    double ptmp[NPR];
    double aco, bco, cco, delE;
    struct of_geom geom ;
    struct of_state q ;
    double mhdf[NDIM][NDIM];
    double mhdi[NDIM][NDIM];
    int i,j,k,m;
    
    vxavg = 0;
    vyavg = 0;
    
    // 1. Calculate and remove mean velocities
    ZSLOOP(0,N1-1,0,N2-1,0,0)
    {
        delvx[i][j] = creal(delvx[i][j]);  //Keep only real part
        delvy[i][j] = creal(delvy[i][j]);
        vxavg += delvx[i][j]/(1.*N1*N2);   //Compute average
        vyavg += delvy[i][j]/(1.*N1*N2);
     }
    
    // 2. Energy normalization setup
    bco = 0.;
    aco = 0.;
    Ekick = 0;
    delE = 0;

    // 3. Remove mean flow and compute energy coefficients
    ZSLOOP(0,N1-1,0,N2-1,0,0)
    {
        
        delvx[i][j] = delvx[i][j]-vxavg;  //Keep average velocity zero
        delvy[i][j] = delvy[i][j]-vyavg;

        
//        PLOOP ptmp[m] = p[i][j][k][m];
//        
//        ptmp[i][j][k][U1] += norm*delvx[i][j];
//        ptmp[i][j][k][U2] += norm*delvy[i][j];
//        get_geometry(i,j,k,CENT,&geom) ;
//        get_state(ptmp, &geom, &q) ;
//        mhd_calc(ptmp, 0, &q, mhdf[0]) ;
//        get_geometry(i,j,k,CENT,&geom) ;
//        get_state(p[i][j][k], &geom, &q) ;
//        mhd_calc(p[i][j][k], 0, &q, mhdi[0]) ;
//        
//        delE += mhdf[0,0]-mhdi[0,0];

        // Calculate coefficients for quadratic equation
        // aco: coefficient of norm²
        aco  += .5*p[i][j][k][RHO]*(pow(creal(delvx[i][j]),2.)+pow(creal(delvy[i][j]),2.));
        // bco: coefficient of norm
        bco  += p[i][j][k][RHO]*(creal(delvx[i][j])*p[i][j][k][U1]+creal(delvy[i][j])*p[i][j][k][U2]);
  }
  
     // 4. Set up energy constraint
    cco = -Ein*N1*N2;
  
    // 5. Solve quadratic equation for normalization factor
    if(aco>0.){
    norm = (-bco + sqrt(bco*bco-4.*aco*cco))/(2.*aco);
    }
    else{
        norm = 0.;
    }
  
    // 6. Apply normalization
    ZSLOOP(0,N1-1,0,N2-1,0,0)
    {
        delvx[i][j] = delvx[i][j]*norm;
        delvy[i][j] = delvy[i][j]*norm;
}


}

/*
For transmissive polar BC, possible solutions from AI in fluxcalc() function:
if(dir == 2 && (j+mpi_startn[2] == 0 || j+mpi_startn[2] == mpi_ntot[2]) && BL) {
    // Instead of zeroing fluxes, compute them with extra care
    PLOOP {
        // Use a more diffusive but stable flux near poles
        double pole_factor = sin(get_geometry(i,j,k,CENT).gcov[2][2]); // Get sin(θ)
        
        // Gradually blend regular flux with a more stable version
        F[i][j][k][m] = pole_factor * F[i][j][k][m] + 
            (1.0 - pole_factor) * 0.5*(F_l[m] + F_r[m]); 
        
        // Add extra numerical viscosity near poles if needed
        if(m == U2) { // Special treatment for θ-velocity component
            F[i][j][k][m] *= exp(-1.0/pole_factor);
        }
    }
}

ChatGPT:
For example, if the cell-centered primitives near the pole are smoothly varying, you might compute the flux as:
Fpole≈1/2(F(p_{active}(θ=0+))+ F(p_{active}(θ=0+)))−1/2 cmax[U(p_{active}(θ=0+))−U(p_{active}(θ=0+))],

which essentially reduces to using the active zone flux rather than zero
 */

/* 
   Schematic diagram for computing E1 (emf[1]) using flux CT

   We want to compute:
     emf[1][i][j][k] = 0.25 * [ (F2[i][j][k][B3] + F2[i][j][k-1][B3])
                                 - (F3[i][j][k][B2] + F3[i][j-1][k][B2]) ]
   
   Diagram overview:
   ---------------------------
   Consider a cell (with center at index (i,j,k)) and its neighboring cell centers.
   In our scheme, the electric field E1 is defined at the edge where the x2 and x3 faces meet.
   The following diagram shows the transverse (j, k) plane at fixed i (the x1 direction):

                          k (cell-center)
                             ●  (F2[i][j][k] with B3)
                             |
                             |     ← (EMF is computed on this edge)
                             |
                ---------------------------
                |          |            |
         j (cell)●----------+------------●
                |          |            |
                |          |            |
                |          |            |
                ---------------------------
                             | 
                             |  (F2[i][j][k-1] with B3 is from the cell shifted in k)
                             ●  (k-1, cell-center)
                             |
                (F3 contributions come from the j shift)
         F3[i][j][k][B2] comes from the top cell,
         F3[i][j-1][k][B2] from the bottom cell.
   
   More explicitly:
   -------------------------------------------------------------------
                   <-- x2 direction (j) -->
              j (active)      j-1 (active)
                ●---------------●  
                |               |  
       k (active) |    E1       |  (The red double-arrow here indicates
      cell-center |  (edge)     |   that the E1 is conceptually located along 
                |   where the face values are combined)
                ●---------------●  
              j (active)      j-1 (active)
                (F3 fluxes are taken from these cells in the x3 direction)
   -------------------------------------------------------------------
   
   In words:
   - The F2 fluxes (in the x2 direction) provide contributions from two adjacent cell faces 
     in the x3 direction:
       • F2[i][j][k][B3] from the cell at (i,j,k)
       • F2[i][j][k-1][B3] from the cell adjacent in the k-direction.
   - The F3 fluxes (in the x3 direction) provide contributions from two adjacent cell faces 
     in the x2 direction:
       • F3[i][j][k][B2] from the cell at (i,j,k)
       • F3[i][j-1][k][B2] from the cell adjacent in the j-direction.
   - These are averaged (with weight 0.25) to represent the flux at the edge.
   - The subtraction (F2 sum minus F3 sum) is dictated by the antisymmetry of the cross–product,
     corresponding to E1 = (F2 B3 - F3 B2).
     
   Note: Even though pr, F2, F3 are cell-centered quantities, the CT scheme reconstructs the 
   “face” or “edge” values by averaging over the four nearby cell centers that border the physical face.

  /*
   3D Schematic for Computing E₁ (emf[1]) via Constrained Transport (CAN BE WRONG!)
   --------------------------------------------------------------
             (x₃ direction)
                ↑
                |
                |        [F₂ flux contributions]
                |       (B₃ component)
                |      ↖       ↗ 
                |    +-----+-----+
                |    |     |     |   <-- Cell centers in the (x₂,x₃) plane at fixed x₁
                |    |  A  |  B  |   where the fluxes F₂ are defined.
                |    +-----+-----+
                |    |     |     |
           (x₂)← |----+-----+-----+----→ (x₂ direction)
                |    |  C  |  D  |   The adjacent cells have centers at C and D.
                |    +-----+-----+
                |      
                |     (F₃ flux contributions)
                |      (B₂ component)
                |      ↙       ↘ 
                |
                | 
         --------------------------------------------------------------
                ↑
                |   <-- This vertical line (in x₁) indicates the 
                |       “cell face” separating cells along x₁.
                |
                *  <-- The edge (at the intersection of the two faces)
                     is where E₁ is defined.
                     
   Explanation:
   - The diagram shows a slice of the grid in the plane transverse to x₁.
   - The four active cell centers (labeled A, B, C, D) are used to approximate the
     flux through the x₁-face.
   - F₂ (flux in the x₂ direction) is associated with the B₃ component; its values 
     are taken from, for example, cells A and B (or with one shift in the k direction).
   - F₃ (flux in the x₃ direction) is associated with the B₂ component; its values 
     are taken from cells A and C (or with one shift in the j direction).
   - At the edge (marked by the asterisk “*”), the constrained transport (CT)
     algorithm computes:
         E₁ ≈ 0.25 * [ (F₂ from A + F₂ from the cell behind in x₃)
                        - (F₃ from A + F₃ from the cell behind in x₂) ]
   - This averaging over four contributions (one from each “corner” of the face) 
     provides a second–order accurate estimate for the edge–centered electric field.
   - Although all the primitive and flux values are originally cell-centered, the 
     CT update reconstructs an effective face value (by averaging the adjacent cell 
     centers) and then computes its difference across the face; the result (E₁) is 
     then assigned to the cell edge.
*/
