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

/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"
#include <float.h>


void coord_transform(double *pr,int i, int j,int k) ;
double compute_Amax( double (*A)[N2+D2][N3+D3] );
double compute_B_from_A( double (*A)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] );
double compute_B_from_AAth( double (*A)[N2+D2][N3+D3], double (*A2)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] );
double normalize_B_by_maxima_ratio(double beta_target, double (*p)[N2M][N3M][NPR], double *norm_value);
double normalize_B_by_beta(double beta_target, double (*p)[N2M][N3M][NPR], double rmax, double *norm_value);

/////////////////////
//magnetic field geometry and normalization
#define NORMALFIELD (0)
#define MADFIELD (1)
#define SEMIMAD (2)

#define WHICHFIELD SEMIMAD

#define NORMALIZE_FIELD_BY_MAX_RATIO (1)
#define NORMALIZE_FIELD_BY_BETAMIN (2)
#define WHICH_FIELD_NORMALIZATION NORMALIZE_FIELD_BY_MAX_RATIO
//end magnetic field
//////////////////////

//////////////////////
//torus density normalization
#define THINTORUS_NORMALIZE_DENSITY (1)
#define DOAUTOCOMPUTEENK0 (1)

#define NORMALIZE_BY_TORUS_MASS (1)
#define NORMALIZE_BY_DENSITY_MAX (2)

#if(DONUCLEAR)
#define DENSITY_NORMALIZATION NORMALIZE_BY_TORUS_MASS
#else
#define DENSITY_NORMALIZATION NORMALIZE_BY_DENSITY_MAX
#endif
//torus density normalization
//////////////////////

double rmax = 0.;
double rhomax = 1.;

void init()
{
  void init_bondi(void);
  void init_torus(void);
  void init_torus_grb(void);
  void init_Shock(void);
  void init_sndwave(void);
  void init_entwave(void);
  void init_OT(void) ;
  void init_monopole(double Rout_val);
  void init_SET(void);
  void init_modes(void);
  void init_vtran(void);
  void init_turb(void);
  void init_grad(void);
  void init_statcond(void);
  void init_statcond1D(void);
  void init_atm(void);
  void init_bondicon(void);
  void init_Hubble(void);

  switch( WHICHPROBLEM ) {
  case MONOPOLE_PROBLEM_1D:
  case MONOPOLE_PROBLEM_2D:
    init_monopole(1e3);
    break;
  case BZ_MONOPOLE_2D:
    init_monopole(100.);
    break;
  case TORUS_PROBLEM:
    init_torus();
    break;
  case TORUS_GRB:
    init_torus_grb();
    break;
  case SHOCK_TEST :
    init_Shock() ;
    break ;
  case SNDWAVE_TEST :
    init_sndwave() ;
    break ;
  case ENTWAVE_TEST :
    init_entwave() ;
    break ;
  case OT_TEST :
    init_OT() ;
    break ;
  case SHOCK_ENTROPY_TEST:
    init_SET();
    break;
  case LIN_MODES :
    init_modes();
    break;
  case VTRAN_TEST:
    init_vtran();
    break;
  case TURB:
    init_turb();
    break;
  case GRAD:
    init_grad();
    break;
  case STATCOND:
    init_statcond();
    break;
case STATCOND1D:
    init_statcond1D();
    break;
  case ATM_TEST:
    init_atm();
    break;
  case BONDI_CON:
    init_bondicon();
    break;
      case HUBBLE:
          init_Hubble();
          break;
          

          
          
  case BONDI_PROBLEM_1D:
  case BONDI_PROBLEM_2D:
    init_bondi();
    break;
  }

}

void init_torus()
{
  int i,j,k ; //Loop indices for iterating over the grid for each sub-domain depending upon the MPI Procs, i.e. N1, N2 and N3 defined in decs.h 
  double r,th,phi,sth,cth ; //Spherical polar coordinates (radius, theta, and azimuthal angle). sin(theta) and cos(theta) for efficiency in trigonometric calculations.
  double ur,uh,up,u,rho ; // Radial, polar, and azimuthal velocity components. Internal energy and density, respectively.
  double X[NDIM] ; //Array to store coordinates (likely 3D for space-time in GR due 3+1 decomposition).
  struct of_geom geom ; //A structure to store geometric quantities (likely related to spacetime or grid).

  /* tilted coordinates */
  double tilt, ctilt, stilt ; //Angle of tilt (in radians). ctilt, stilt: cos(tilt) and sin(tilt) for efficiency. 
  double sph,cph, thp, php, xp, yp, zp ; //Transformed (tilted) theta and phi coordinates. Transformed Cartesian coordinates (related to the tilt).
  double sth0, cth0, sph0, cph0 ; //What is the difference between sth0 and sth? and others... 
  double dth_dthp, dth_dphp, dph_dthp, dph_dphp ; // Derivatives of transformed coordinates for mapping between systems.\partial(θ)/\partial(θ_tilted), etc.
  double dphp_dth, dphp_dph ; // Derivatives of transformed coordinates for mapping between systems.\partial(φ_tilted)/\partial(φ), etc.

  /* for disk interior */
  /*l: Specific angular momentum, an important quantity for determining disk stability.
  rin: Inner radius of the torus.
  lnh: Logarithm of the density or enthalpy, often used in GRMHD.
  DD, AA, SS: Geometric factors related to the metric in GR.
  thin, sthin, cthin: Spherical coordinates (theta) specific to the torus structure.
  kappa: Entropy parameter.
  hm1: Enthalpy minus one, used in relativistic hydrodynamics.*/
  double l,rin,lnh,expm2chi,up1 ; 
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;

  /* for magnetic field */
  double Aphip[N1+D1][N2+D2][N3+D3] ; //Vector potential phi component in tilted coordinates
  double Ath[N1+D1][N2+D2][N3+D3] ; //Vector potential Theta component in UN-tilted coordinates
  double Aphi[N1+D1][N2+D2][N3+D3] ; //Vector potential Phicomponent in UN-tilted coordinates
  double rho_av,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ; //Average density, (norm) Normalization factor for ensuring consistent field strength.
  double lfish_calc(double rmax) ;
  
  int iglob, jglob, kglob;
  double rancval;
  
  double amax, aphipow; //aphipow to chose the radial scaling/power for vector potential
  
  double Amin, Amax, cutoff_frac = 0.01;
  
  /* some physics parameters */
  gam = 5./3. ;

  /* disk parameters (use fishbone.m to select new solutions) */
  a = 0.9375 ;
  rin = 12. ;
  rmax = 25. ;
  l = lfish_calc(rmax) ; //Function to compute the specific angular momentum (l) at rmax, based on Fishbone-Moncrief solutions for a relativistic torus.

  /* tilt angle parameter in radians */
  tilt = 0.0 ;

  kappa =1.e-3;
  beta = 100. ;

  // Define new parameters R0_field and w_field at disk initialization
  double Rcyl; 
  double R0_field = rmax; // in units of rg Central radius of the magnetic field distribution.
  double w_field = rmax - 5;  // in units of rg Width of the magnetic field profile.
  double C_norm_beta = 1./beta ;    // Normalization constant for vector potential  Normalization constant to ensure a specific magnetization (ββ) is achieved.
  double Rout_field = 100; 

  /* some numerical parameters */
  lim = MC ; //Chooses a flux limiter for numerical schemes.
  failed = 0 ;	/* start slow */ //Initializes failure flag for diagnostics.
  cour = .8 ; //Courant number, controls time step size based on grid resolution.
  dt = 1.e-5 ; // Initial time step size.
  R0 = 0.0 ;
  Rin = 0.87*(1. + sqrt(1. - a*a)) ;  //.98
  Rout = 1e3; //changed from 1e5 to 1e3
  rbr = 400.;
  npow2=4.0; //power exponent
  cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)

  t = 0. ;
  hslope = 0.3 ;

  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.0;

#if(BL==2 || BL==3)
  /////////////////////
  //ANGULAR GRID SETUP (so far irrelevant for WHICHPROBLEM==NSTAR)
  /////////////////////
  
  //transverse resolution fraction devoted to different components
  //(sum should be <1)
  global_fracdisk = 0.25; //fraction of resolution going into the disk
  global_fracjet = 0.40; //fraction of resolution going into the jets
  
  global_disknu1 = -2;
  global_disknu2 = 0.75;
  
  global_jetnu1 = -2;  //the nu-parameter that determines jet shape at small radii (r < min(r0jet,r0disk))
  global_jetnu2 = 0.75;  //the nu-parameter that determines jet shape at large radii (r > r0jet)
  
  //not used any more
  global_rsjet = 0.0;
  
  //distance at which theta-resolution is *exactly* uniform in the jet grid -- want to have this at Rin;
  //the larger r0grid, the larger the thickness of the jet
  global_r0grid = Rin;
  
  //distance out to which jet decollimates
  //the larger it is, the wider is the jet grid
  global_r0jet = 2*Rin;
  
  //distance beyond which the jet grid stops collimating and becomes radial
  global_rjetend = 1e3;
  
  //distance out to which the disk decollimates
  //the larger r0disk, the thinner is the disk grid
  global_r0disk = 2*Rin;
  
  //not used any more
  global_rdiskend = 5*Rin;
#endif
  //cylindrification parameters
  global_x10 = 5.;  //radial distance in MCOORD until which the innermost angular cell is cylinrdical
  global_x20 = -1. + 1./256.; //1./mpi_ntot[2];     //This restricts grid cylindrification to the one
  //single grid cell closest to the pole (other cells virtually unaffeced, so there evolution is accurate).
  //This trick minimizes the resulting pole deresolution and relaxes the time step.
  //The innermost grid cell is evolved inaccurately whether you resolve it or not, and it will be fixed
  //by POLEFIX (see bounds.c).

  set_arrays() ;
  set_grid() ;

  get_phys_coord(5,0,0,&r,&th,&phi) ;
  if(MASTER==mpi_rank) {
    fprintf(stderr,"r[5]: %g\n",r) ;
    fprintf(stderr,"r[5]/rhor: %g",r/(1. + sqrt(1. - a*a))) ;
    if( r > 1. + sqrt(1. - a*a) ) {
      fprintf(stderr, ": INSUFFICIENT RESOLUTION, ADD MORE CELLS INSIDE THE HORIZON\n" );
    }
    else {
      fprintf(stderr, "\n");
    }
  }

  /* output choices */
  tf = 20000.0 ;

  DTd = 10.; /* dumping frequency, in units of M */
  //increase values to match those chosen upon restarting:
  DTl = 100. ;	/* logfile frequency, in units of M */
  DTi = 100. ; 	/* image file frequ., in units of M */
  DTr = 100. ; /* restart file frequ., in units of M */
  DTr01 = 2500. ; /* restart file frequ., in timesteps */

  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;

  rhomax = 0. ; //Track the maximum density and internal energy in the domain for normalization or diagnostic purposes.
  umax = 0. ;
  //ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
  //Loops through the global grid using a decomposition structure (mpi_ntot for global dimensions and mpi_startn for offsets)
  //iglob, jglob, kglob: Iterate over the global indices of the simulation grid. hese loops cover the entire domain defined by mpi_ntot
  for(iglob=0;iglob<mpi_ntot[1];iglob++) {
    for(jglob=0;jglob<mpi_ntot[2];jglob++) {
      for(kglob=0;kglob<mpi_ntot[3];kglob++) {
        
        rancval = ranc(0); //Random value (seeded externally) used for small perturbations later. 
        /*Converts global indices (iglob, etc.) to local indices (i, etc.) for this processor using the offset mpi_startn*/
        i = iglob-mpi_startn[1];
        j = jglob-mpi_startn[2];
        k = kglob-mpi_startn[3];
        //Ensures that indices fall within the local grid on this processor. If not, skips the iteration.
        if(i<0 ||
           j<0 ||
           k<0 ||
           i>=N1 ||
           j>=N2 ||
           k>=N3){
          continue;
        }
        get_phys_coord(i,j,k,&r,&th,&phi) ; //Retrieves physical coordinates (r,θ,ϕr,θ,ϕ) corresponding to the grid cell (i, j, k).

	// primed coordinates 
  //Precomputes trigonometric values of the tilt angle and spherical coordinates to simplify calculations later.
	ctilt = cos(tilt) ;
	stilt = sin(tilt) ;

  sth0 = sin(th) ;
  cth0 = cos(th) ;
	sph0 = sin(phi) ;
	cph0 = cos(phi) ;
        
  /* Transforms spherical coordinates (r,θ,ϕr,θ,ϕ) into the tilted coordinate system (r,θp,ϕpr,θp​,ϕp​).
    xp,yp,zp​: Cartesian-like intermediate variables in the tilted frame.
    θp,ϕp​: Tilted spherical angles calculated using inverse trigonometric functions.*/
	xp = ctilt*sth0*cph0 - stilt*cth0 ;
	yp = sth0*sph0 ;
	zp = stilt*sth0*cph0 + ctilt*cth0 ;

	thp = acos(zp) ;
	php = atan2(yp,xp) ;

  //Precomputes trigonometric functions of the tilted coordinates for later use.
	sth = sin(thp) ;
	cth = cos(thp) ;
	sph = sin(php) ;
	cph = cos(php) ;

        /* calculate lnh */
        /*DD,AA,SS: Standard Kerr metric quantities in Boyer-Lindquist coordinates.
          Depend on rr, θθ, and black hole spin aa.*/
        DD = r*r - 2.*r + a*a ;
        AA = (r*r + a*a)*(r*r + a*a) - DD*a*a*sth*sth ;
        SS = r*r + a*a*cth*cth ;

        /*Defines similar metric quantities at the inner boundary radius rinrin​.
          Assumes θin=π/2 (midplane of the torus).*/
        thin = M_PI/2. ; //This is π/2
        sthin = sin(thin) ;
        cthin = cos(thin) ;
        DDin = rin*rin - 2.*rin + a*a ;
        AAin = (rin*rin + a*a)*(rin*rin + a*a) 
                - DDin*a*a*sthin*sthin ;
        SSin = rin*rin + a*a*cthin*cthin ;

        /*Computes ln⁡h (logarithmic specific enthalpy) inside the torus (r≥rin​).
          Outside the torus (r<rin​), assigns ln⁡h=1*/
        if(r >= rin) {
          lnh = 0.5*log((1. + sqrt(1. + 4.*(l*l*SS*SS)*DD/
                  (AA*sth*AA*sth)))/(SS*DD/AA)) 
                  - 0.5*sqrt(1. + 4.*(l*l*SS*SS)*DD/(AA*AA*sth*sth))
                  - 2.*a*r*l/AA 
                  - (0.5*log((1. + sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
                  (AAin*AAin*sthin*sthin)))/(SSin*DDin/AAin)) 
                  - 0.5*sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
                          (AAin*AAin*sthin*sthin)) 
                  - 2.*a*rin*l/AAin ) ;
        }
        else
          lnh = 1. ;

        /* regions outside torus */
        /*Sets density (ρ) and internal energy (u) to zero outside the torus.
          Velocities (ur,uθ,uϕ) are set to zero to avoid unphysical flows. */
        if(lnh < 0. || r < rin) {
          //reset density and internal energy to zero outside torus
          rho = 0.; //1.e-7*RHOMIN ;
          u = 0.; //1.e-7*UUMIN ;

          /* these values are demonstrably physical
             for all values of a and r */
          /*
          ur = -1./(r*r) ;
          uh = 0. ;
          up = 0. ;
          */

          ur = 0. ;
          uh = 0. ;
          up = 0. ;

          /*
          get_geometry(i,j,CENT,&geom) ;
          ur = geom.gcon[0][1]/geom.gcon[0][0] ;
          uh = geom.gcon[0][2]/geom.gcon[0][0] ;
          up = geom.gcon[0][3]/geom.gcon[0][0] ;
          */

          p[i][j][k][RHO] = rho ;
          p[i][j][k][UU] = u ;
          p[i][j][k][U1] = ur ;
          p[i][j][k][U2] = uh ;
          p[i][j][k][U3] = up ;
        }
        /* region inside magnetized torus; u^i is calculated in
         * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
         * so it needs to be transformed at the end */
        //Uses ln⁡h to compute density and internal energy based on a polytropic equation of state (P=κ ρ^γ)
        else { 
          hm1 = exp(lnh) - 1. ;
          rho = pow(hm1*(gam - 1.)/(kappa*gam),
                                  1./(gam - 1.)) ; 
          u = kappa*pow(rho,gam)/(gam - 1.) ;
          ur = 0. ;
          uh = 0. ;

          /* calculate u^phi => Computes the azimuthal velocity (uϕ) based on specific angular momentum l. */
          expm2chi = SS*SS*DD/(AA*AA*sth*sth) ;
          up1 = sqrt((-1. + sqrt(1. + 4.*l*l*expm2chi))/2.) ;
          up = 2.*a*r*sqrt(1. + up1*up1)/sqrt(AA*SS*DD) +
                  sqrt(SS/AA)*up1/sth ;


          p[i][j][k][RHO] = rho ;
          if(rho > rhomax) rhomax = rho ;
          p[i][j][k][UU] = u*(1. + 4.e-2*(rancval-0.5)) ;
          if(u > umax && r > rin) umax = u ;
          p[i][j][k][U1] = ur ;

    // Transforms the θ- and ϕ-components of velocity (uθ,uϕ) to account for the tilt.
	  /* transform uh = dthp/dtau as dth/dtau = dthp/dtau dth/dthp + dphp/dtau dth/dphp */
          dth_dthp = (ctilt*sth + stilt*cth*cph) / sth0;
          dth_dphp = -stilt*sth*sph / sth0;

          /* transform dph/dtau = dthp/dtau dthp/dph + dphp/dtau dph/dphp */
          dph_dthp = stilt*sph/(sth0*sth0);
          dph_dphp = sth/(sth0*sth0)*(ctilt*sth + stilt*cth*cph);

          p[i][j][k][U2] = uh*dth_dthp + up*dth_dphp;
          p[i][j][k][U3] = uh*dph_dthp + up*dph_dphp ;

          //p[i][j][k][U2] = uh ;

	  //          p[i][j][k][U3] = up ;

          /* convert from 4-vel in BL coords to relative 4-vel in code coords */
          coord_transform(p[i][j][k],i,j,k) ;
        }

        //Initializes magnetic field components to zero for now. 
        p[i][j][k][B1] = 0. ;
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
      }
    }
  }

/*MPI_Allreduce synchronizes values across all MPI processes.
The operation ensures all processes have the global maximum density (rhomax) and internal energy (umax).
The MPI_IN_PLACE directive tells MPI to overwrite rhomax and umax with the reduced (maximum) values.*/
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  /* Normalize the densities so that max(rho) = 1 */
  if(MASTER==mpi_rank) fprintf(stderr,"rhomax: %g\n",rhomax) ; //fprintf: The master MPI process outputs the global maximum density (rhomax).
  /* Loops over all grid cells within the physical domain (N1, N2, N3).

    The macro iterates through i,j,ki,j,k indices within their specified ranges.
    Inside the loop: Density (RHO) and internal energy (UU) are normalized to make the maximum density equal to 1.*/
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
          p[i][j][k][RHO] /= rhomax ;
          p[i][j][k][UU]  /= rhomax ;
  }
  umax /= rhomax ; //(maximum internal energy) is also normalized.
  kappa *= pow(rhomax,gam-1);
  global_kappa = kappa; //kappa is adjusted using the equation of state, and its value is stored globally.
  rhomax = 1. ; //rhomax is set to 1 as a new reference

  bound_prim(p) ; //Applies boundary conditions to the primitive variables p => Ensures consistency of density, velocity, and magnetic fields at the domain boundaries.

  /*Determines the configuration of the magnetic field (aphipow) based on the problem setup.
    NORMALFIELD: Sets aphipow to 0 (no radial dependency).
    MADFIELD and SEMIMAD: Set aphipow to 2 (stronger radial dependency, suitable for MAD initial conditions).
    If WHICHFIELD is undefined, an error is raised.*/
  if (WHICHFIELD == NORMALFIELD) {
    aphipow = 0.;
  }
  // with a big torus this initial condition will be a MAD
  else if(WHICHFIELD == MADFIELD || WHICHFIELD == SEMIMAD) {
    aphipow = 2;
  }
  else {
    fprintf(stderr, "Unknown field type: %d\n", (int)WHICHFIELD);
    exit(321);
  }

  /* first find corner-centered vector potential */
  /*The ZSLOOP iterates over the entire grid, including ghost zones (as extended by D1, D2, D3); NG = 3 (Ghost Cells) as defined in decs.h.
    Initializes all components of the vector potential (Aphip, Aphi, Ath) to zero. This ensures a clean slate for further computations.*/
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) Aphip[i][j][k] = 0. ;
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) Aphi[i][j][k] = 0. ;
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) Ath[i][j][k] = 0. ;


  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
          /* radial field version */
          /*
          coord(i,j,k,CORN,X) ;
          bl_coord(X,&r,&th,&phi) ;
         
          A[i][j][k] = (1-cos(th)) ;
          */

    
          /* vertical field version */
          /*
          coord(i,j,k,CORN,X) ;
          bl_coord(X,&r,&th,&phi) ;

          A[i][j][k] = r*r*sin(th)*sin(th) ;
          */
    
    

          /* field-in-disk version */
          /* flux_ct */
    
      //cannot use get_phys_coords() here because it can only provide coords at CENT
      coord(i,j,k,CORN,X) ; //coord() computes Cartesian coordinates (X) at cell corners.
      bl_coord(X,&r,&th,&phi) ; //bl_coord() converts X to Boyer-Lindquist coordinates (r,θ,ϕr,θ,ϕ)
      ctilt = cos(tilt) ;
      stilt = sin(tilt) ;

      sth0 = sin(th) ;
      cth0 = cos(th) ;
      sph0 = sin(phi) ;
      cph0 = cos(phi) ;
        
      // rotated cartesian coordinates based on tilt angle.
      xp = ctilt*sth0*cph0 - stilt*cth0 ;
      yp = sth0*sph0 ;
      zp = stilt*sth0*cph0 + ctilt*cth0 ;

      // spherical coordinates in tilted frame.
      thp = acos(zp) ;
      php = atan2(yp,xp) ;

      sth = sin(thp) ;
      cth = cos(thp) ;
      sph = sin(php) ;
      cph = cos(php) ;

      rho_av = 0.25*(
		     p[i][j][k][RHO] +
		     p[i-1][j][k][RHO] +
		     p[i][j-1][k][RHO] +
		     p[i-1][j-1][k][RHO]) ; //Average density is computed from surrounding grid cells, hower 'k' is not changed, so its like a 2D square, not 3D cube for average

      q = pow(r,aphipow)*rho_av/rhomax ; //The magnitude of Aϕ_tilted​ depends on: q=r^aphipow (ρ_av/rhomax)
      if (WHICHFIELD == NORMALFIELD) {
	q -= 0.2;
      }
      if (WHICHFIELD == MADFIELD || WHICHFIELD == SEMIMAD) {
	q *= q;
      }
      
      //Aphip[i][j][k] = -1*q; // Aphip polarity is switched to change direction of magnetic field
      
  //     if (sth*sth > 0.) {

	// dphp_dth = -stilt * sph0 / (sth*sth);
	// dphp_dph = sth0/(sth*sth)*(ctilt * sth0 - stilt * cth0 * cph0);
	
  // //Aphi and Ath components are vector pot components in untilted coords derived from Aphip using partial derivatives of the transformed coordinates.
	// Aphi[i][j][k] = Aphip[i][j][k]*dphp_dph ;
	// Ath[i][j][k] = Aphip[i][j][k]*dphp_dth ;
      // TEST with axisymmetry but Ath non-zero: divb = 0 still
      //Aphi[i][j][k] = Aphip[i][j][k] ;
      //Ath[i][j][k] = 0.1*Aphip[i][j][k] ;
      // TEST with non-axisymmetry but Ath = 0
      //Aphi[i][j][k] = Aphip[i][j][k]*dphp_dph ;
      //}
      
      Rcyl = r * sin(th);  // Cylindrical radius
      // Apply the piecewise conditions for A_phi
      if (Rcyl < rin) {
          Aphi[i][j][k] = C_norm_beta * tanh((rin - R0_field) / w_field);
      }
      else if (Rcyl >= rin && Rcyl <= Rout_field) {
          Aphi[i][j][k] = C_norm_beta * tanh((Rcyl - R0_field) / w_field);
      }
      else {
          Aphi[i][j][k] = C_norm_beta * tanh((Rout_field - R0_field) / w_field);
      }
    
      
  }

  fixup(p) ; //in fixup.c; Ensures primitive variables p are physically valid and satisfy constraints (e.g., positive density, proper magnetic field divergence).

  /* now differentiate to find cell-centered B,
     and begin normalization */
  
  //bsq_max = compute_B_from_A(Aphi,p);
  /*This line calls the function compute_B_from_AAth to compute the maximum magnetic field strength, bsq_max. 
  It uses the magnetic vector potential components Aphi and Ath (with Ath = 0 for axisymmetric or vertical field case) and the primitive variables p */
  bsq_max = compute_B_from_AAth(Aphi,Ath,p); 
  //what goes wrong here?
  //if(MASTER==mpi_rank)
  //  fprintf(stderr,"A B1 B2 B3: %g %g %g %g\n",Aphi[30][18][10],p[30][18][10][B1],p[30][18][10][B2],p[30][18][10][B3]) ;
  
  //This checks the field type (NORMALFIELD or SEMIMAD)
  if(WHICHFIELD == NORMALFIELD || WHICHFIELD == SEMIMAD) {
    if(MASTER==mpi_rank)
      fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ; //Prints the initial value of bsq_max if the current processor is the master in an MPI setup.

    /* finally, normalize to set field strength */
    beta_act =(gam - 1.)*umax/(0.5*bsq_max) ; //actual beta is calculated, where umax is the maximum internal energy density and bsq_max is the maximum magnetic field squared.

    if(MASTER==mpi_rank)
      fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ; //This prints the calculated beta_act and compares it with the target beta value.
    
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(beta, p, 10*rmax, &norm);
    }
    else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      /*the code is calling the normalize_B_by_maxima_ratio function to adjust the field strength to achieve the target beta*/
      beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        MPI_Finalize();
        exit(2345);
      }
    }

    if(MASTER==mpi_rank)
      fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ; //    This prints the final value of beta_act after the normalization step.
  }
  
  // The code will not go ino else if(WHICHFIELD == MADFIELD) conditional block as WHICHFIELD is set to SEMIMAD, and thus, I don't need to change Amax for polarity changed in Aphip
  else if(WHICHFIELD == MADFIELD) {
    getmax_densities(p, &rhomax, &umax);
    amax = get_maxprimvalrpow( p, aphipow, RHO );
    if(MASTER==mpi_rank)
      fprintf(stderr, "amax = %g\n", amax);
    
    //by now have the fields (centered and stag) computed from vector potential
    
    //here need to
    //1) compute bsq
    //2) rescale field components such that beta = p_g/p_mag is what I want (constant in the main disk body and tapered off to zero near torus edges)
    normalize_field_local_nodivb( beta, rhomax, amax, p, Aphi, 1 ); //defined in madfield.c
    
    //3) re-compute vector potential by integrating up ucons (plays with MPI and zeros out B[3])
    compute_vpot_from_gdetB1( p, Aphi );
    
    Amax = compute_Amax( Aphi );
    Amin = cutoff_frac * Amax;
    
    //chop off magnetic field close to the torus boundaries
    ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
      if (Aphi[i][j][k] < Amin) {
        Aphi[i][j][k] = Amin;
      }
    }
  
    //4) call user1_init_vpot2field_user() again to recompute the fields
    //convert A to staggered pstag, centered prim and ucons, unsure about Bhat
    compute_B_from_A( Aphi, p );
    
    //5) normalize the field
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(beta, p, 10*rmax, &norm);
    }
    else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        MPI_Finalize();
        exit(2345);
      }
    }
  }

    
  /* enforce boundary conditions */
  fixup(p) ;
  bound_prim(p) ; // These functions apply boundary conditions to the variables p (which hold the primitive quantities).

//Initializing the Global Vector Potential Array:
/*     If the EVOLVEVPOT flag is set (It is set for torus problem), 
this loop initializes the global potential array vpot with the components of the magnetic vector potential (Ath and Aphi), 
scaled by a normalization constant norm. */
#if(EVOLVEVPOT)
  //initialize the global vpot array with the magnetic vector potential
  //currently, only the phi-component is supported initially
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
    vpot[0][i][j][k] = 0.;
    vpot[1][i][j][k] = 0.;
    vpot[2][i][j][k] = Ath[i][j][k]*norm;
    vpot[3][i][j][k] = Aphi[i][j][k]*norm;
  }
#endif



#if( DO_FONT_FIX )
  set_Katm();
#endif 

#if(DOPARTICLES)
  init_particles();
#endif
  


}

//note that only axisymmetric A is supported
double compute_Amax( double (*A)[N2+D2][N3+D3] )
{
  double Amax = 0.;
  int i, j, k;
  struct of_geom geom;
  
  ZLOOP {
    if(A[i][j][k] > Amax) Amax = A[i][j][k];
  }
  
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&Amax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  return(Amax);
}


//note that only axisymmetric A is supported
double compute_B_from_A( double (*A)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] )
{
  double bsq_max = 0., bsq_ij ;
  int i, j, k;
  struct of_geom geom;
  
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;
    
    /* flux-ct */
    p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
                       + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
    p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
                      - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;
    
    p[i][j][k][B3] = 0. ;
    
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  return(bsq_max);
}

// This function calculates the magnetic field from the vector potential.
  double compute_B_from_AAth( double (*A)[N2+D2][N3+D3], double (*A2)[N2+D2][N3+D3], double (*p)[N2M][N3M][NPR] ) //We called bsq_max = compute_B_from_AAth(Aphi,Ath,p); 
{
  {
    /*This initializes bsq_max to track the maximum magnetic field and sets up a loop (ZLOOP) to iterate over the grid. 
    The get_geometry function is used to retrieve the geometry at each point.*/
    double bsq_max = 0., bsq_ij ;
    int i, j, k;
    struct of_geom geom;

    ZLOOP {
      get_geometry(i,j,k,CENT,&geom) ;

      /* flux-ct */
      /* non-axisymmetric, Ar = 0, compared with Athena++ and axisymmetric harmpi, see notes 2/3/2021 */
      
      /*
      In a finite‐volume or staggered‐grid scheme, the “curl” is naturally computed as an average circulation around the cell, 
      so its result is attributed to the cell center. The discrete curl operator is essentially a “circulation” integral around the cell’s boundary divided by the cell’s area. 
      Although the raw data A are at corners, the average of the appropriate differences gives you the cell-centered B.
      The averaging of 8 corner values (4 for each term) directly gives us the cell-centered B field:
Corner(A)   Corner(A)  Corner(A)
    *----------*----------*
    |          |          |
    |     Cell(B)        |
    |          |          |
    *----------*----------*
Corner(A)   Corner(A)  Corner(A)
      In a finite-difference scheme, the derivative ∂Aϕ/∂θ can be approximated by taking the difference between the 
      values of Aϕ​ on the “top” and “bottom” faces of the cell. But since these faces do not coincide with the cell center,
       the algorithm reconstructs an estimate at the center by averaging the four cell-corner values that lie on that face.
      
      */
      // Bi=ϵ^ijk ∂A^k/∂x^j​​
      // p[i][j][k][B1]=1./(\sqrt{g}​) (∂Aϕ/∂θ​​−∂Aθ/∂ϕ​​) =>  Include the metric determinant \sqrt{g}to account for coordinate transformations and ensure the physical interpretation of the field.
      // The first term corresponds to the derivative of Aϕ​ with respect to θ. The second term corresponds to the derivative of Aθ​ with respect to ϕ
      p[i][j][k][B1] = ((A[i][j+1][k] + A[i+1][j+1][k] + A[i][j+1][k+1] + A[i+1][j+1][k+1]
			 - A[i][j][k] - A[i+1][j][k] - A[i][j][k+1] - A[i+1][j][k+1])/dx[2]
			- (A2[i][j][k+1] + A2[i+1][j][k+1] + A2[i][j+1][k+1] + A2[i+1][j+1][k+1]
			   - A2[i][j][k] - A2[i+1][j][k] - A2[i][j+1][k] - A2[i+1][j+1][k])/dx[3])/(4.*geom.g) ;

      // p[i][j][k][B2]=1./(\sqrt{g}​) (∂Ar/∂φ​​−∂Aφ/∂r​​) 
      // The first term corresponds to the derivative of Ar with respect to φ, which is = 0 as Ar = 0. The second term corresponds to the derivative of Aφ​ with respect to r
      p[i][j][k][B2] = -(A[i+1][j][k] + A[i+1][j+1][k] + A[i+1][j][k+1] + A[i+1][j+1][k+1]
			- A[i][j][k] - A[i][j+1][k] - A[i][j][k+1] - A[i][j+1][k+1])/(4.*dx[1]*geom.g) ;

      // p[i][j][k][B3]=1./(\sqrt{g}​) (∂Aθ/∂r​​−∂Ar/∂θ​​) 
      // The first term corresponds to the derivative of Aθ with respect to r. The second term corresponds to the derivative of Ar​ with respect to θ,which is = 0 as Ar = 0
      p[i][j][k][B3] = (A2[i+1][j][k] + A2[i+1][j+1][k] + A2[i+1][j][k+1] + A2[i+1][j+1][k+1]
      		       - A2[i][j][k] - A2[i][j+1][k] - A2[i][j][k+1] - A2[i][j+1][k+1])/(4.*dx[1]*geom.g) ;
      
      /*
      
      //TESTING 
      p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
			 + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
      p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
			- A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;
      
      p[i][j][k][B3] = 0. ;
      //END TESTING
      */

     /*
     returns b^2 (i.e., twice magnetic pressure)
      double bsq_calc(double *pr, struct of_geom *geom)
      {
        struct of_state q ;
        
        get_state(pr,geom,&q) ;
        return( dot(q.bcon,q.bcov) ) ;
      }*/
      bsq_ij = bsq_calc(p[i][j][k],&geom) ;
      if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    }

//If using MPI, this line ensures that bsq_max is updated across all processors. The function then returns the maximum magnetic field squared value.
#ifdef MPI
    //exchange the info between the MPI processes to get the true max
    MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

    return(bsq_max);
  }
}


//This function normalizes the magnetic field based on the target beta value.
double normalize_B_by_maxima_ratio(double beta_target, double (*p)[N2M][N3M][NPR], double *norm_value) //We called it as: beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
{
  double beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  int i, j, k;
  struct of_geom geom;
  
  // This loop computes the maximum values of bsq_ij and umax across the grid
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    u_ij = p[i][j][k][UU];
    if(u_ij > umax) umax = u_ij;
  }

// These MPI calls ensure that umax and bsq_max are the maximum values across all processors.
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  /* finally, normalize to set field strength */
  //This calculates the actual beta_act value and computes the normalization factor norm
  beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;
  norm = sqrt(beta_act/beta_target) ;
  bsq_max = 0. ;

  //This applies the normalization to the magnetic field components (B1, B2, and B3) and recalculates the maximum bsq_max
  ZLOOP {
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
    p[i][j][k][B3] *= norm ;
    
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }

//This MPI call updates the maximum bsq_max across all processors
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  //This recalculates the beta_act after normalization.
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;

  if(norm_value) {
    *norm_value = norm; //If the norm_value pointer is provided, it stores the normalization factor
  }
  return(beta_act); //The function returns the final beta_act value after normalization
}

//normalize the magnetic field using the values inside r < rmax
double normalize_B_by_beta(double beta_target, double (*p)[N2M][N3M][NPR], double rmax, double *norm_value)
{
  double beta_min = 1e100, beta_ij, beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  int i, j, k;
  struct of_geom geom;
  double X[NDIM], r, th, ph;
  
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th, &ph);
    if (r>rmax) {
      continue;
    }
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    u_ij = p[i][j][k][UU];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  /* finally, normalize to set field strength */
  beta_act = beta_min;
  
  norm = sqrt(beta_act/beta_target) ;
  beta_min = 1e100;
  ZLOOP {
    p[i][j][k][B1] *= norm ;
    p[i][j][k][B2] *= norm ;
    p[i][j][k][B3] *= norm ;
    get_geometry(i,j,k,CENT,&geom) ;
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    u_ij = p[i][j][k][UU];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  beta_act = beta_min;

  if(norm_value) {
    *norm_value = norm;
  }

  return(beta_act);
}

/////////////////////////////////////////////////////////////////
//
// init_torus_grb() preliminaries
//
/////////////////////////////////////////////////////////////////
#define JMAX 100

double rtbis(double (*func)(double,double*), double *parms, double x1, double x2, double xacc)
//Using bisection, find the root of a function func known to lie between x1 and x2. The root,
//returned as rtbis, will be refined until its accuracy is \pm xacc.
//taken from http://gpiserver.dcom.upv.es/Numerical_Recipes/bookcpdf/c9-1.pdf
{
  int j;
  double dx,f,fmid,xmid,rtb;
  f=(*func)(x1, parms);
  fmid=(*func)(x2, parms);
  if (f*fmid >= 0.0) {
    fprintf( stderr, "f(%g)=%g f(%g)=%g\n", x1, f, x2, fmid );
    fprintf(stderr, "Root must be bracketed for bisection in rtbis\n");
    exit(434);
  }
  rtb = (f < 0.0) ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so that f>0 lies at x+dx.
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5),parms); //Bisection loop.
    if (fmid <= 0.0) {
      rtb=xmid;
    }
    if (fabs(dx) < xacc || fmid == 0.0) {
      return rtb;
    }
  }
  fprintf(stderr, "Too many bisections in rtbis");
  return 0.0; //Never get here.
}

double lfunc( double lin, double *parms )
{
  double gutt, gutp, gupp, al, c;
  double ans;
  
  gutt = parms[0];
  gutp = parms[1];
  gupp = parms[2];
  al = parms[3];
  c = parms[4];
  
  ans = (gutp - lin * gupp)/( gutt - lin * gutp) - c *pow(lin/c,al); // (lin/c) form avoids catastrophic cancellation due to al = 2/n - 1 >> 1 for 2-n << 1
  
  return(ans);
}

void compute_gu( double r, double th, double a, double *gutt, double *gutp, double *gupp )
{
  //metric (expressions taken from eqtorus_c.nb):
  *gutt = -1 - 4*r*(pow(a,2) + pow(r,2))*
  pow((-2 + r)*r + pow(a,2),-1)*
  pow(pow(a,2) + cos(2*th)*pow(a,2) + 2*pow(r,2),-1);
  
  *gutp = -4*a*r*pow((-2 + r)*r + pow(a,2),-1)*
  pow(pow(a,2) + cos(2*th)*pow(a,2) + 2*pow(r,2),-1);
  
  *gupp = 2*((-2 + r)*r + pow(a,2)*pow(cos(th),2))*
  pow(sin(th),-2)*pow((-2 + r)*r + pow(a,2),-1)*
  pow(pow(a,2) + cos(2*th)*pow(a,2) + 2*pow(r,2),-1);
}

double thintorus_findl( double r, double th, double a, double c, double al )
{
  double gutt, gutp, gupp;
  double parms[5];
  double l;
  
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  //store params in an array before function call
  parms[0] = gutt;
  parms[1] = gutp;
  parms[2] = gupp;
  parms[3] = al;
  parms[4] = c;
  
  //solve for lin using bisection, specify large enough root search range, (1e-3, 1e3)
  //demand accuracy 5x machine prec.
  //in non-rel limit l_K = sqrt(r), use 10x that as the upper limit:
  l = rtbis( &lfunc, parms, 1, 10*sqrt(r), 5.*DBL_EPSILON );
  
  return( l );
}

double compute_udt( double r, double th, double a, double l )
{
  double gutt, gutp, gupp;
  double udt;
  
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  udt = -sqrt(- 1/(gutt - 2 * l * gutp + l * l * gupp));
  
  return( udt );
}


double compute_omega( double r, double th, double a, double l )
{
  double gutt, gutp, gupp;
  double omega;
  
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  omega = (gutp - gupp*l)*pow(gutt - gutp*l,-1);
  
  return( omega );
}

double compute_uuphi( double r, double th, double a, double l )
{
  double gutt, gutp, gupp;
  double udt, udphi, uuphi;
  
  //u_t
  udt = compute_udt(r, th, a, l);

  //u_phi
  udphi = -udt * l;
  
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  //u^phi
  uuphi = gutp * udt + gupp * udphi;
  return( uuphi );
}

double compute_l_from_omega( double r, double th, double a, double omega )
{
  double gutt, gutp, gupp;
  double l;
  
  compute_gu( r, th, a, &gutt, &gutp, &gupp );
  
  l = (gutp - omega * gutt)/(gupp - omega * gutp);
  
  return( l );
}

void init_torus_grb()
{
  double compute_udt( double r, double th, double a, double l );
  double compute_uuphi( double r, double th, double a, double l );
  double compute_omega( double r, double th, double a, double l );

  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  
  /* for disk interior */
  double l,rin,lnh,expm2chi,up1 ;
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;
  
  /* for magnetic field */
  double A[N1+D1][N2+D2][N3+D3] ;
  double rho_av,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double rmax, lfish_calc(double rmax) ;
  
  int iglob, jglob, kglob;
  double rancval;
  double omk, lk, kk, c, n, al, lin, utin, udt, flin, rho_scale_factor, f, hh, eps, rhofloor, ufloor, rho_at_pmax;
  double torus_mass;
  
  double amax, aphipow;
  
  double Amin, Amax, cutoff_frac = 0.001;
  
  const double frac_pert = 4.e-2;

  /* some physics parameters */
#if(DONUCLEAR)
  gam = 4./3. ;
  /* radial distribution of angular momentum */
  n = 0.0;  // = 0 constant ang. mom. torus; 0.25 standard setting for large MAD torii
  
  /* disk parameters (use fishbone.m to select new solutions) */
  a = 0.8 ;
  rmax = 11.; //11.209750769799999 ;
  rin = 0.616 * rmax ;
#else
  gam = 5./3. ;
  /* radial distribution of angular momentum */
  n = 0.25;  // = 0 constant ang. mom. torus; 0.25 standard setting for large MAD torii
  
  /* disk parameters (use fishbone.m to select new solutions) */
  a = 0.5 ;
  rin = 15. ;
  rmax = 34.475 ;
#endif
  
  kappa =1.e-3;
  beta = 1.e2 ;
  
  /* some numerical parameters */
  lim = MC ;
  failed = 0 ;	/* start slow */
  cour = .8 ;
  dt = 1.e-5 ;
  R0 = 0.0 ;
  Rin = 0.87*(1. + sqrt(1. - a*a)) ;  //.98
  Rout = 1e5;
  rbr = 200.;
  npow2=4.0; //power exponent
  cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
  
  
  t = 0. ;
  hslope = 0.3 ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;
  
#if(BL==2 || BL==3)
  ////////////////////
  //RADIAL GRID SETUP
  ////////////////////
  rbr = 400.;  //make it smaller than the default 1000. while trying to make thin disks work
  
  /////////////////////
  //ANGULAR GRID SETUP (so far irrelevant for WHICHPROBLEM==NSTAR)
  /////////////////////
  
  //transverse resolution fraction devoted to different components
  //(sum should be <1)
  global_fracdisk = 0.25; //fraction of resolution going into the disk
  global_fracjet = 0.40; //fraction of resolution going into the jets
  
  global_disknu1 = -2.;
  global_disknu2 = -2.;
  
  global_jetnu1 = -2.;  //the nu-parameter that determines jet shape at small radii (r < min(r0jet,r0disk))
  global_jetnu2 = 0.75;  //the nu-parameter that determines jet shape at large radii (r > r0jet)
  
  //not used any more
  global_rsjet = 0.0;
  
  //distance at which theta-resolution is *exactly* uniform in the jet grid -- want to have this at Rin;
  //the larger r0grid, the larger the thickness of the jet
  global_r0grid = Rin;
  
  //distance out to which jet decollimates
  //the larger it is, the wider is the jet grid
  global_r0jet = 2*Rin;
  
  //distance beyond which the jet grid stops collimating and becomes radial
  global_rjetend = 1e3;
  
  //distance out to which the disk decollimates
  //the larger r0disk, the thinner is the disk grid
  global_r0disk = Rin;
  
  //not used any more
  global_rdiskend = 1.e7;
#endif
  //cylindrification parameters
  global_x10 = 3.0;  //radial distance in MCOORD until which the innermost angular cell is cylinrdical
  global_x20 = -1. + 1./mpi_ntot[2];     //This restricts grid cylindrification to the one
  //single grid cell closest to the pole (other cells virtually unaffeced, so there evolution is accurate).
  //This trick minimizes the resulting pole deresolution and relaxes the time step.
  //The innermost grid cell is evolved inaccurately whether you resolve it or not, and it will be fixed
  //by POLEFIX (see bounds.c).
  
  set_arrays() ;
  set_grid() ;
  
  get_phys_coord(5,0,0,&r,&th,&phi) ;
  if(MASTER==mpi_rank) {
    fprintf(stderr,"r[5]: %g\n",r) ;
    fprintf(stderr,"r[5]/rhor: %g",r/(1. + sqrt(1. - a*a))) ;
    if( r > 1. + sqrt(1. - a*a) ) {
      fprintf(stderr, ": INSUFFICIENT RESOLUTION, ADD MORE CELLS INSIDE THE HORIZON\n" );
    }
    else {
      fprintf(stderr, "\n");
    }
  }
  
  /* output choices */
  tf = 20000.0 ;
  
  DTd = 10.; /* dumping frequency, in units of M */
  DTl = 10. ;	/* logfile frequency, in units of M */
  DTi = 10. ; 	/* image file frequ., in units of M */
  DTr = 10. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  
  rhomax = 0. ;
  umax = 0. ;
  torus_mass = 0.;
  //ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
  for(iglob=0;iglob<mpi_ntot[1];iglob++) {
    for(jglob=0;jglob<mpi_ntot[2];jglob++) {
      for(kglob=0;kglob<mpi_ntot[3];kglob++) {
        
        rancval = ranc(0);
        i = iglob-mpi_startn[1];
        j = jglob-mpi_startn[2];
        k = kglob-mpi_startn[3];
        if(i<0 ||
           j<0 ||
           k<0 ||
           i>=N1 ||
           j>=N2 ||
           k>=N3){
          continue;
        }
        
        ///
        /// Computations at pressure max
        ///
        
        r = rmax;
        th = M_PI_2;
        omk = 1./(pow(rmax,1.5) + a);
        lk = compute_l_from_omega( r, th, a, omk );
        //log(omk) == (2/n) log(c) + (1 - 2./n) log(lk) <-- solve for c:
        c = pow( lk, 1 - n/2. ) * pow( omk, n/2. );
        
        if (0 != n) {
          //variable l case
          kk = pow(c, 2./n);
          al = (n - 2.)/n;
          
          
          ///
          /// Computations at torus inner edge
          ///
          
          //l = lin at inner edge, r = rin
          r = rin;
          th = M_PI_2;
          lin = thintorus_findl( r, th, a, c, al );
          
          //finding DHK03 lin, utin, f (lin)
          utin = compute_udt( r, th, a, lin );
          flin = pow(fabs(1 - kk*pow(lin,1 + al)),pow(1 + al,-1));
          
          
          ///
          /// EXTRA CALCS AT PRESSURE MAX TO NORMALIZE DENSITY
          ///
          //l at pr. max r, th
  #if( THINTORUS_NORMALIZE_DENSITY )
          r = rmax;
          th = M_PI_2;
          l = lk;
          udt = compute_udt( r, th, a, l );
          f = pow(fabs(1 - kk*pow(l,1 + al)),pow(1 + al,-1));
          hh = flin*utin*pow(f,-1)*pow(udt,-1);
          eps = (-1 + hh)*pow(gam,-1);
          rho_at_pmax = pow((-1 + gam)*eps*pow(kappa,-1),pow(-1 + gam,-1));
          rho_scale_factor = 1.0/rho_at_pmax;
  #if( DOAUTOCOMPUTEENK0 )
          //will be recomputed for every (ti,tj,tk), but is same for all of them, so ok
          global_kappa = kappa * pow( rho_scale_factor, 1 - gam );
  #endif
          
  #else
          rho_scale_factor = 1.0;
  #endif
          
          ///
          /// Computations at current point: r, th
          ///
          
          coord(i, j, k, CENT, X);
          bl_coord(X, &r,&th,&phi);
          
          //l at current r, th
          if (r>=rin) {
            l = thintorus_findl( r, th, a, c, al );
            udt = compute_udt( r, th, a, l );
            f = pow(fabs(1 - kk*pow(l,1 + al)),pow(1 + al,-1));
            hh = flin*utin*pow(f,-1)*pow(udt,-1);
          }
          else {
            l = udt = f = hh = 0.;
          }
        }
        else{
          //l = constant case
          l = c;
          ///
          /// Computations at torus inner edge
          ///
          //l = lin at inner edge, r = rin
          r = rin;
          th = M_PI_2;
          lin = l; //constant ang. mom.
          
          //finding DHK03 lin, utin, f (lin)
          utin = compute_udt( r, th, a, lin );
          flin = 1.;  //f(l) is unity everywhere according to Chakrabarti

          ///
          /// EXTRA CALCS AT PRESSURE MAX TO NORMALIZE DENSITY
          ///
          //l at pr. max r, th
#if( THINTORUS_NORMALIZE_DENSITY )
          r = rmax;
          th = M_PI_2;
          l = lk;
          udt = compute_udt( r, th, a, l );
          f = 1.;
          hh = flin*utin*pow(f,-1)*pow(udt,-1);
          eps = (-1 + hh)*pow(gam,-1);
          rho_at_pmax = pow((-1 + gam)*eps*pow(kappa,-1),pow(-1 + gam,-1));
          rho_scale_factor = 1.0/rho_at_pmax;
#if( DOAUTOCOMPUTEENK0 )
          //will be recomputed for every (ti,tj,tk), but is same for all of them, so ok
          global_kappa = kappa * pow( rho_scale_factor, 1 - gam );
#endif
          
#else
          rho_scale_factor = 1.0;
#endif

          ///
          /// Computations at current point: r, th
          ///
          
          coord(i, j, k, CENT, X);
          bl_coord(X, &r,&th,&phi);

          if (r>=rin) {
            udt = compute_udt( r, th, a, l );
            f = 1.;
            hh = utin/udt;
          }
          else {
            hh = 0;
          }
        }
        eps = (-1 + hh)*pow(gam,-1);
        rho = pow((-1 + gam)*eps*pow(kappa,-1),pow(-1 + gam,-1));
        
        //compute atmospheric values
        get_geometry(i, j, k, CENT, &geom); // true coordinate system
        get_phys_coord(i,j,k,&r,&th,&phi) ;
        get_rho_u_floor(r, th, phi, &rhofloor, &ufloor);
#       if DONUCLEAR
        //floor values according to Rodrigo
          rhofloor *= 1.5;
          ufloor *= 1.5;
#       endif
        /* regions outside torus */
        if( r < rin || isnan(eps) || eps < 0 || rho*rho_scale_factor < rhofloor) {
          /* these values are demonstrably physical
           for all values of a and r */
          rho = 0;
          u = 0;
          
          /*
           ur = -1./(r*r) ;
           uh = 0. ;
           up = 0. ;
           */
          
          ur = 0. ;
          uh = 0. ;
          up = 0. ;
          
          /*
           get_geometry(i,j,CENT,&geom) ;
           ur = geom.gcon[0][1]/geom.gcon[0][0] ;
           uh = geom.gcon[0][2]/geom.gcon[0][0] ;
           up = geom.gcon[0][3]/geom.gcon[0][0] ;
           */
          
          p[i][j][k][RHO] = rho ;
          p[i][j][k][UU] = u ;
          p[i][j][k][U1] = ur ;
          p[i][j][k][U2] = uh ;
          p[i][j][k][U3] = up ;
        }
        /* region inside magnetized torus; u^i is calculated in
         * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
         * so it needs to be transformed at the end */
        else {
          u = kappa * pow(rho, gam) / (gam - 1.);
          
          rho *= rho_scale_factor;
          u *= rho_scale_factor;
          
          ur = 0. ;
          uh = 0. ;
          up = compute_uuphi( r, th, a, l );
          
          p[i][j][k][RHO] = rho ;
          if(rho > rhomax) rhomax = rho ;
          p[i][j][k][UU] = u*(1. + frac_pert*(rancval-0.5)) ;
          if(u > umax && r > rin) umax = u ;
          p[i][j][k][U1] = ur ;
          p[i][j][k][U2] = uh ;
          
          p[i][j][k][U3] = up ;
          
          /* convert from 4-vel in BL coords to relative 4-vel in code coords */
          coord_transform(p[i][j][k],i,j,k) ;
          
          //add up mass to compute total torus mass
          torus_mass += gdet[i][j][k][CENT] * rho * dV;
        }
        
        p[i][j][k][B1] = 0. ;
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
      }
    }
  }
  torus_mass /= fracphi;  //account for the size of phi wedge when computing torus mass
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&torus_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  
  /* Normalize the densities so that max(rho) = 1 */
  if(MASTER==mpi_rank) fprintf(stderr,"Before normalization: rhomax: %g, torus_mass: %g\n",rhomax, torus_mass) ;
  if (DENSITY_NORMALIZATION == NORMALIZE_BY_DENSITY_MAX) {
    rho_scale_factor = 1./rhomax;
    if(MASTER==mpi_rank) fprintf(stderr,"Normalizing by rhomax = 1:\n");
  }
  else if (DENSITY_NORMALIZATION == NORMALIZE_BY_TORUS_MASS) {
    //a factor of fracphi accounts for missing mass outside the wedge
    rho_scale_factor = 0.01*fracphi/torus_mass;
    if(MASTER==mpi_rank) fprintf(stderr,"Normalizing by torus_mass = 0.01:\n");
  }
  torus_mass = 0.;
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    p[i][j][k][RHO] *= rho_scale_factor ;
    p[i][j][k][UU]  *= rho_scale_factor ;
    //add up mass to compute total torus mass
    torus_mass += gdet[i][j][k][CENT] * p[i][j][k][RHO] * dV;
  }
  torus_mass /= fracphi;  //account for the size of phi wedge when computing torus mass
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&torus_mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
  umax *= rho_scale_factor ;
  rhomax *= rho_scale_factor ;
  if(MASTER==mpi_rank) fprintf(stderr,"After normalization: rhomax: %g, torus_mass: %g\n",rhomax, torus_mass) ;

  if (WHICHFIELD == NORMALFIELD) {
    aphipow = 0.;
  }
  else if(WHICHFIELD == MADFIELD || WHICHFIELD == SEMIMAD) {
    //dependence on polytropic index such that r^aphipow*rho ~ constant along disk midplane
    aphipow = 2.5 / (3.*(gam-1.));
  }
  else {
    fprintf(stderr, "Unknown field type: %d\n", (int)WHICHFIELD);
    exit(321);
  }

  //need to bound density before computing vector potential
  bound_prim(p) ;

  /* first find corner-centered vector potential */
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) A[i][j][k] = 0. ;
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
    /* radial field version */
    /*
     coord(i,j,k,CORN,X) ;
     bl_coord(X,&r,&th,&phi) ;
     
     A[i][j][k] = (1-cos(th)) ;
     */
    
    
    /* vertical field version */
    /*
     coord(i,j,k,CORN,X) ;
     bl_coord(X,&r,&th,&phi) ;
     
     A[i][j][k] = r*r*sin(th)*sin(th) ;
     */
    
    
    
    /* field-in-disk version */
    /* flux_ct */
    
    //cannot use get_phys_coords() here because it can only provide coords at CENT
    coord(i,j,k,CORN,X) ;
    bl_coord(X,&r,&th,&phi) ;
    
    
    rho_av = 0.25*(
                   p[i][j][k][RHO] +
                   p[i-1][j][k][RHO] +
                   p[i][j-1][k][RHO] +
                   p[i-1][j-1][k][RHO]) ;
    
    q = pow(r,aphipow) * rho_av/rhomax ;
    if (WHICHFIELD == NORMALFIELD) {
      q -= 0.2;
    }
    if (WHICHFIELD == MADFIELD || WHICHFIELD == SEMIMAD) {
      //this gives Aphi ~ (r^aphipow*rho)^2
      q = q*q;
    }
    if(q > 0.) A[i][j][k] = q ;  //cos(th) gives two loops.
    
  }

  //need to apply the floor on density to avoid beta ~ ug/(bsq+SMALL) = 0 outside torus when normalizing B
  fixup(p) ;

  /* now differentiate to find cell-centered B,
   and begin normalization */
  bsq_max = compute_B_from_A(A,p);
  
  if(WHICHFIELD == NORMALFIELD || WHICHFIELD == SEMIMAD) {
    if(MASTER==mpi_rank)
      fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;
    
    /* finally, normalize to set field strength */
    beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;
    
    if(MASTER==mpi_rank)
      fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
    
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(beta, p, rmax, &norm);
    }
    else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        MPI_Finalize();
        exit(2345);
      }
    }
  
    if(MASTER==mpi_rank)
      fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;
  }
  else if(WHICHFIELD == MADFIELD) {
    getmax_densities(p, &rhomax, &umax);
    amax = get_maxprimvalrpow( p, aphipow, RHO );
    if(MASTER==mpi_rank)
      fprintf(stderr, "amax = %g\n", amax);
    
    //by now have the fields computed from vector potential
    
    //here:
    //1) computing bsq
    //2) rescaling field components such that beta = p_g/p_mag is what I want
    //   (constant in the main disk body and tapered off to zero near torus edges)
    normalize_field_local_nodivb( beta, rhomax, amax, p, A, 1 );

    //3) re-compute vector potential by integrating up \int B^r dA in theta
    //   (this uses MPI and zeros out B[3] because B[3] is used to communicate the integration results)
    compute_vpot_from_gdetB1( p, A );
    
    Amax = compute_Amax( A );
    Amin = cutoff_frac * Amax;

    //chop off magnetic field close to the torus boundaries
    ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
      if (A[i][j][k] < Amin) {
        A[i][j][k] = Amin;
      }
    }

    //4) recompute the fields by converting the new A to B
    compute_B_from_A( A, p );

    //5) normalize the field
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(beta, p, rmax, &norm);
    }
    else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      beta_act = normalize_B_by_maxima_ratio(beta, p, &norm);
    }
    else {
      if(i_am_the_master) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        MPI_Finalize();
        exit(2345);
      }
    }
  }
  
#if DONUCLEAR
  //enforce the floors
  ZLOOP {
    p[i][j][k][YE] = 0.1; //initially start with Ye = 0.1
    //compute atmospheric values
    get_geometry(i, j, k, CENT, &geom); // true coordinate system
    get_phys_coord(i,j,k,&r,&th,&phi) ;
    get_rho_u_floor(r, th, phi, &rhofloor, &ufloor);
    //floor values slightly elevated for Rodrigo
    rhofloor *= 1.5;
    ufloor *= 1.5;
    if (p[i][j][k][RHO]<rhofloor) {
      p[i][j][k][RHO]=rhofloor;
      p[i][j][k][RHOFLOOR]=rhofloor;
    }
    if (p[i][j][k][UU]<ufloor) {
      p[i][j][k][UU]=ufloor;
    }
  }
#endif
  /* enforce boundary conditions */
  fixup(p) ;
  bound_prim(p) ;
  
#if(EVOLVEVPOT)
  //initialize the global vpot array with the magnetic vector potential
  //currently, only the phi-component is supported initially
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
    vpot[0][i][j][k] = 0.;
    vpot[1][i][j][k] = 0.;
    vpot[2][i][j][k] = 0.;
    vpot[3][i][j][k] = A[i][j][k] * norm;
  }
#endif

#if( DO_FONT_FIX )
  set_Katm();
#endif 
  
#if(DOPARTICLES)
  init_particles();
#endif
  
}

void init_bondi()
{
	int i,j,k ;
	double r,th,phi,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1][N3+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax, lfish_calc(double rmax) ;

	/* some physics parameters */
	gam = 4./3. ;

	/* black hole parameters */
        a = 0.9375 ;

	kappa = 1.e-3 ;

	/* radius of the inner edge of the initial density distribution */
	rin = 10.;

        /* some numerical parameters */
        lim = MC ;
        failed = 0 ;	/* start slow */
        cour = 0.9 ;
        dt = 1.e-5 ;
	rhor = (1. + sqrt(1. - a*a)) ;
	R0 = -2*rhor ;
        Rin = 0.5*rhor ;
        Rout = 1e3 ;
        rbr = Rout*10.;
        npow2=4.0; //power exponent
        cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)


        t = 0. ;
        hslope = 1.0 ; //uniform angular grid

	if(N2!=1) {
	  //2D problem, use full pi-wedge in theta
	  fractheta = 1.;
	}
	else{
	  //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
	  fractheta = 1.e-2;
	}
        fracphi = 1.;

        set_arrays() ;
        set_grid() ;

	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;
	fprintf(stderr,"rmin: %g\n",r) ;
	fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        /* output choices */
	tf = Rout ;

	DTd = 10. ;	/* dumping frequency, in units of M */
	DTl = 10. ;	/* logfile frequency, in units of M */
	DTi = 10. ; 	/* image file frequ., in units of M */
        DTr = 10 ; /* restart file frequ., in units of M */
        DTr01 = 1000 ; /* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
        rdump01_cnt = 0 ;
	defcon = 1. ;

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
		coord(i,j,k,CENT,X) ;
		bl_coord(X,&r,&th,&phi) ;

		sth = sin(th) ;
		cth = cos(th) ;

		/* regions outside uniform density distribution */
		if(r < rin) {
			rho = 1.e-7*RHOMIN ;
                        u = 1.e-7*UUMIN ;

			/* these values are demonstrably physical
			   for all values of a and r */
			/*
                        ur = -1./(r*r) ;
                        uh = 0. ;
			up = 0. ;
			*/

			ur = 0. ;
			uh = 0. ;
			up = 0. ;

			/*
			get_geometry(i,j,CENT,&geom) ;
                        ur = geom.gcon[0][1]/geom.gcon[0][0] ;
                        uh = geom.gcon[0][2]/geom.gcon[0][0] ;
                        up = geom.gcon[0][3]/geom.gcon[0][0] ;
			*/

			p[i][j][k][RHO] = rho ;
			p[i][j][k][UU] = u ;
			p[i][j][k][U1] = ur ;
			p[i][j][k][U2] = uh ;
			p[i][j][k][U3] = up ;
		}
		/* region inside initial uniform density */
		else { 
		  rho = 1.;
		  u = kappa*pow(rho,gam)/(gam - 1.) ;
		  ur = 0. ;
		  uh = 0. ;


		  p[i][j][k][RHO] = rho ;
		  if(rho > rhomax) rhomax = rho ;
		  p[i][j][k][UU] = u;
		  if(u > umax && r > rin) umax = u ;
		  p[i][j][k][U1] = ur ;
		  p[i][j][k][U2] = uh ;
		  p[i][j][k][U3] = up ;
		  
		  /* convert from 4-vel to 3-vel */
		  coord_transform(p[i][j][k],i,j,k) ;
		}

		p[i][j][k][B1] = 0. ;
		p[i][j][k][B2] = 0. ;
		p[i][j][k][B3] = 0. ;

	}

	fixup(p) ;
	bound_prim(p) ;
    

#if(0) //disable for now
	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2,0,N3) A[i][j][k] = 0. ;
        ZSLOOP(0,N1,0,N2,0,N3) {
                /* vertical field version */
                /*
                coord(i,j,l,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;

                A[i][j][k] = 0.5*r*sin(th) ;
                */

                /* field-in-disk version */
		/* flux_ct */
                rho_av = 0.25*(
                        p[i][j][RHO] +
                        p[i-1][j][RHO] +
                        p[i][j-1][RHO] +
                        p[i-1][j-1][RHO]) ;

                q = rho_av/rhomax - 0.2 ;
                if(q > 0.) A[i][j][k] = q ;

        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,k,CENT,&geom) ;

		/* flux-ct */
		p[i][j][B1] = -(A[i][j][k] - A[i][j+1][k]
				+ A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
		p[i][j][B2] = (A[i][j][k] + A[i][j+1][k]
				- A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;

		p[i][j][B3] = 0. ;

		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
	norm = sqrt(beta_act/beta) ;
	bsq_max = 0. ;
	ZLOOP {
		p[i][j][k][B1] *= norm ;
		p[i][j][k][B2] *= norm ;

		get_geometry(i,j,k,CENT,&geom) ;
		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;
	fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;

	/* enforce boundary conditions */
	fixup(p) ;
	bound_prim(p) ;
    
#endif

    

    
#if( DO_FONT_FIX ) 
	set_Katm();
#endif 

#if(DOPARTICLES)
  init_particles();
#endif


}

void init_monopole(double Rout_val)
{
	int i,j,k ;
	double r,th,phi,sth,cth ;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor;

	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1 ;

	/* for magnetic field */
	double A[N1+1][N2+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax, lfish_calc(double rmax) ;

	/* some physics parameters */
	gam = 4./3. ;

	/* disk parameters (use fishbone.m to select new solutions) */
        a = 0.9375 ;
        rin = 6. ;
        rmax = 12. ;
        l = lfish_calc(rmax) ;

	kappa = 1.e-3 ;
	beta = 1.e2 ;

        /* some numerical parameters */
        lim = MC ;
        failed = 0 ;	/* start slow */
        cour = 0.9 ;
        dt = 1.e-5 ;
	rhor = (1. + sqrt(1. - a*a)) ;
	R0 = -4*rhor ;
        Rin = 0.7*rhor ;
        Rout = Rout_val ;
        rbr = Rout*10.;
    npow2=4.0; //power exponent
    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)

        t = 0. ;
        hslope = 1. ;

	if(N2!=1) {
	  //2D problem, use full pi-wedge in theta
	  fractheta = 1.;
	}
	else{
	  //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
	  fractheta = 1.e-2;
	}

        fracphi = 1.;

        set_arrays() ;
        set_grid() ;

	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&r,&th,&phi) ;
	fprintf(stderr,"rmin: %g\n",r) ;
	fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;

        /* output choices */
	tf = 2*Rout ;

	DTd = 1. ;	/* dumping frequency, in units of M */
	DTl = 50. ;	/* logfile frequency, in units of M */
	DTi = 50. ; 	/* image file frequ., in units of M */
        DTr = 1. ; /* restart file frequ., in units of M */
	DTr01 = 1000 ; 	/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0 ;
	image_cnt = 0 ;
	rdump_cnt = 0 ;
        rdump01_cnt = 0 ;
	defcon = 1. ;

	rhomax = 0. ;
	umax = 0. ;
	ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
	  coord(i,j,k,CENT,X) ;
	  bl_coord(X,&r,&th,&phi) ;

	  sth = sin(th) ;
	  cth = cos(th) ;

	  /* rho = 1.e-7*RHOMIN ; */
	  /* u = 1.e-7*UUMIN ; */

	  /* rho = pow(r,-4.)/BSQORHOMAX; */
	  /* u = pow(r,-4.*gam)/BSQOUMAX; */

	  rho = RHOMINLIMIT+(r/10./rhor)/pow(r,4)/BSQORHOMAX;
	  u = UUMINLIMIT+(r/10./rhor)/pow(r,4)/BSQORHOMAX;

	  /* these values are demonstrably physical
	     for all values of a and r */
	  /*
	    ur = -1./(r*r) ;
	    uh = 0. ;
	    up = 0. ;
	  */

	  ur = 0. ;
	  uh = 0. ;
	  up = 0. ;

	  /*
	    get_geometry(i,j,CENT,&geom) ;
	    ur = geom.gcon[0][1]/geom.gcon[0][0] ;
	    uh = geom.gcon[0][2]/geom.gcon[0][0] ;
	    up = geom.gcon[0][3]/geom.gcon[0][0] ;
	  */

	  p[i][j][k][RHO] = rho ;
	  p[i][j][k][UU] = u ;
	  p[i][j][k][U1] = ur ;
	  p[i][j][k][U2] = uh ;
	  p[i][j][k][U3] = up ;
	  p[i][j][k][B1] = 0. ;
	  p[i][j][k][B2] = 0. ;
	  p[i][j][k][B3] = 0. ;
	}

	rhomax = 1. ;
	fixup(p) ;
	bound_prim(p) ;

        //leave A[][] a 2D array for now, which means that magnetic field will be axisymmetric
	/* first find corner-centered vector potential */
	ZSLOOP(0,N1,0,N2,0,0) A[i][j] = 0. ;
        ZSLOOP(0,N1,0,N2,0,0) {
#if(0)
                /* vertical field version */
                coord(i,j,k,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;
                A[i][j] = 0.5*pow(r*sin(th),2);
#elif(1)
                /* radial (monopolar) field version */
                coord(i,j,k,CORN,X) ;
                bl_coord(X,&r,&th,&phi) ;
                A[i][j] = (1-cos(th)) ;
#endif

        }

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0. ;
	ZLOOP {
		get_geometry(i,j,k,CENT,&geom) ;

		/* flux-ct */
		p[i][j][k][B1] = -(A[i][j] - A[i][j+1]
				+ A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*geom.g) ;
		p[i][j][k][B2] = (A[i][j] + A[i][j+1]
				- A[i+1][j] - A[i+1][j+1])/(2.*dx[1]*geom.g) ;

		p[i][j][k][B3] = 0. ;

		bsq_ij = bsq_calc(p[i][j][k],&geom) ;
		if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
	}
	fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

	/* enforce boundary conditions */
	fixup(p) ;
	bound_prim(p) ;
    
    



#if( DO_FONT_FIX )
	set_Katm();
#endif 


}

void init_Shock()
{
	int i,j,k,m ;
	double x,y,z,r,sth,cth ;
    double slope, xdis;
	double ur,uh,up,u,rho ;
	double X[NDIM] ;
	struct of_geom geom ;
	double rhor, inc;
    
	/* for disk interior */
	double l,rin,lnh,expm2chi,up1 ;
	double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
	double kappa,hm1,rarefac ;
    
	/* for magnetic field */
	double A[N1+1][N2+1] ;
	double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
	double rmax;
    
	/* some physics parameters */
	gam = 5./3. ;
    rarefac = 1.;  //-1 for diverging, +1 for converging flows
    inc = 0.;
    slope = tan(inc);
    
    
    /* some numerical parameters */
    lim = VANL ;
    failed = 0 ;	/* start slow */
    cour = .9;
    dt = 1.e-5 ;
    
    t = 0. ;
    hslope = 1. ;
    
    if(N2!=1) {
      //2D problem, use full pi-wedge in theta
      fractheta = 1.;
    }
    else{
      //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
      fractheta = 1.e-2;
    }
  
    fracphi = 1.;

    set_arrays() ;
    set_grid() ;
    
	coord(-2,0,0,CENT,X) ;
	bl_coord(X,&x,&y,&z) ;
	fprintf(stderr,"xmin: %g\n",x) ;
  
    /* output choices */
    tf = 600.;

    DTd = tf/1. ;	/* dumping frequency, in units of M */
    DTl = 50000. ;	/* logfile frequency, in units of M */
    DTi = 50000. ; 	/* image file frequ., in units of M */
    DTr = tf/1. ; /* restart file frequ., in units of M */
    DTr01 = 10000 ; /* restart file frequ., in timesteps */

    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;

    rhomax = 0. ;
    umax = 0. ;
    ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
      coord(i,j,k,CENT,X) ;
      bl_coord(X,&x,&y,&z) ;

        xdis = slope*(y-.5)+.5;
      
        if(x<xdis){
        ur = 1./1000. *rarefac ;
      }
      else {
        ur = -1./1000.*rarefac ;
      }
      
      rho = 1. ;
      
        u = (gam-1.)*pow(1./1000.,3.);

      uh = 0. ;
      up = 0. ;
      
        p[i][j][k][RHO] = rho ;
        p[i][j][k][UU] = u ;
        p[i][j][k][U1] = ur*cos(inc) ;
        p[i][j][k][U2] = -ur*sin(inc) ;
        p[i][j][k][U3] = up ;
        p[i][j][k][B1] = 0. ;
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
    }

    rhomax = 1. ;
    fixup(p) ;
    bound_prim(p) ;
  
    /* first find corner-centered vector potential */
    ZSLOOP(0,N1,0,N2,0,0) A[i][j] = 0. ;
  
    
    /* now differentiate to find cell-centered B,
     and begin normalization */
    bsq_max = 0. ;
    ZLOOP {
      get_geometry(i,j,k,CENT,&geom) ;

      /* flux-ct */
      p[i][j][k][B1] = -(A[i][j] - A[i][j+1]
              + A[i+1][j] - A[i+1][j+1])/(2.*dx[2]*geom.g) ;
      p[i][j][k][B2] = (A[i][j] + A[i][j+1]
             - A[i+1][j] - A[i+1][j+1])/(2.*dx[1]*geom.g) ;

      p[i][j][k][B3] = 0. ;

      bsq_ij = bsq_calc(p[i][j][k],&geom) ;
      if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    }
    fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

    /* enforce boundary conditions */


    fixup(p) ;
    bound_prim(p) ;

    
  
#if( DO_FONT_FIX ) 
    set_Katm();
#endif 
    
    
}

void init_entwave()
{
  int i,j,k ;
  double x,y,z,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double rhor;
  
  double myrho, myu, mycs, myv;
  double delta_rho;
  double cosa, sina;
  double delta_ampl = 1.e-1; //amplitude of the wave
  double k_vec_x = 2 * M_PI;  //wavevector
  double k_vec_y = 2 * M_PI;
  double k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );
  double tfac = 1.e3; //factor by which to reduce velocity
  

  /* some physics parameters */
  gam = 5./3. ;
  
  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;

  set_arrays() ;
  set_grid() ;
  
  
  myrho = 1.;
  myu = 4.*myrho / (gam * (gam-1));  //so that mycs is unity
  
  mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed
  myv = 1.;  //velocity with a magnitude of 1
  
  /* output choices */

    tf = tfac;///mycs;
  
  DTd = tf/10. ;	/* dumping frequency, in units of M */
  DTl = tf/10. ;	/* logfile frequency, in units of M */
  DTi = tf/10. ; 	/* image file frequ., in units of M */
  DTr = tf/10. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; 	/* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  

  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&x,&y,&z) ;
    
    //applying the perturbations
    delta_rho = delta_ampl * cos( k_vec_x * x + k_vec_y * y );
    
  //  p[i][j][k][RHO] = myrho + delta_rho;
      if(x<.2){
          p[i][j][k][RHO] = 1.;
          p[i][j][k][U2] = 0.;

      }
      else if (x>.8){
          p[i][j][k][RHO] = 1.;
          p[i][j][k][U2] = 0.;
      }
      else{
          p[i][j][k][RHO] = 1e4;
          p[i][j][k][U2] = 0.;

      }
    p[i][j][k][UU] = myu/(tfac*tfac);
    p[i][j][k][U1] = myv/tfac;
    //p[i][j][k][U2] = myv/tfac;
    //p[i][j][k][U2] = myv/tfac;//cos(k_vec_x*x);
    p[i][j][k][U3] = 0 ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
  }
  
  /* enforce boundary conditions */
  
  
  fixup(p) ;
  bound_prim(p) ;
  
  
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
  
  
}

void init_sndwave()
{
  int i,j,k ;
  double x,y,z,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  
  double myrho, myu, mycs, myv;
  double delta_rho;
  double cosa, sina;
  double delta_ampl = 1e-5; //amplitude of the wave
  double k_vec_x = 2 * M_PI;  //wavevector
  double k_vec_y = 0;
  double k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );
  double tfac = 1e4; //factor by which to reduce velocity
  
  
  /* some physics parameters */
  gam = 5./3. ;
  
  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;
  
  set_arrays() ;
  set_grid() ;
  
  
  myrho = 1.;
  myu = myrho / (gam * (gam-1));  //so that mycs is unity
  
  mycs = sqrt(gam * (gam-1) * myu / myrho);  //background sound speed
  
  /* output choices */
  
  tf = tfac/mycs;
  
  DTd = tf/10. ;	/* dumping frequency, in units of M */
  DTl = tf/10. ;	/* logfile frequency, in units of M */
  DTi = tf/10. ; 	/* image file frequ., in units of M */
  DTr = tf/10. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&x,&y,&z) ;
    
    //applying the perturbations
    delta_rho = delta_ampl * cos( k_vec_x * x + k_vec_y * y );
    
    p[i][j][k][RHO] = myrho + delta_rho;
    p[i][j][k][UU] = (myu + gam * myu * delta_rho / myrho)/(tfac*tfac);
    p[i][j][k][U1] = (delta_rho/myrho * mycs * k_vec_x / k_vec_len)/tfac;
    p[i][j][k][U2] = (delta_rho/myrho * mycs * k_vec_y / k_vec_len)/tfac;
    p[i][j][k][U3] = 0 ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
  }
  
  /* enforce boundary conditions */
  
  
  fixup(p) ;
  bound_prim(p) ;
  
  
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
  
  
}

void init_OT()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double rhor;


  /* for magnetic field */
  double A[N1+1][N2+1][N3+1] ;
  double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double rmax;

  /* some physics parameters */
  gam = 5./3. ;

  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
 
  
  set_arrays() ;
  set_grid() ;
  
  coord(-2,0,0,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  fprintf(stderr,"rmin: %g\n",r) ;

  /* output choices */
  tf = 2000.;

  DTd = 10. ;	/* dumping frequency, in units of M */
  DTl = 50. ;	/* logfile frequency, in units of M */
  DTi = 50000. ; 	/* image file frequ., in units of M */
  DTr = 10 ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */

  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;

  rhomax = 0. ;
  umax = 0. ;
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&r,&th,&phi) ;
    
    rho = 25./36./M_PI ;
    
    u = (gam-1.)*5./12./M_PI/pow(1000.,2);
    
    ur = -sin(2.*M_PI*th)/1000.   ;
    uh = sin(2.*M_PI*r)/1000. ;
    up = 0. ;
    
    p[i][j][k][RHO] = rho ;
    p[i][j][k][UU] = u ;
    p[i][j][k][U1] = ur ;
    p[i][j][k][U2] = uh ;
    p[i][j][k][U3] = up ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
  }

  rhomax = 1. ;
  fixup(p) ;
  bound_prim(p) ;

  /* first find corner-centered vector potential */
  ZSLOOP(0,N1,0,N2,0,N3) A[i][j][k] = 0. ;

  ZSLOOP(0,N1,0,N2,0,N3) {
     coord(i,j,k,CORN,X) ;
     bl_coord(X,&r,&th,&phi) ;
     
     A[i][j][k] = 1./sqrt(4.*M_PI)*(cos(4.*M_PI*r)/(4*M_PI)+cos(2.*M_PI*th)/(2.*M_PI))/1000.;
  }

  
  
  /* now differentiate to find cell-centered B,
   and begin normalization */
  bsq_max = 0. ;
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;

    /* flux-ct */
    p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
            + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
    p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
           - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;

    p[i][j][k][B3] = 0. ;

    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
  fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

  /* enforce boundary conditions */
  fixup(p) ;
  bound_prim(p) ;

  

  
#if( DO_FONT_FIX )
  set_Katm();
#endif
    
    
}

void init_SET()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double rhor;

  /* for disk interior */
  double l,rin,lnh,expm2chi,up1 ;
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;

  /* for magnetic field */
  double A[N1+1][N2+1][N3+1] ;
  double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double rmax;

  /* some physics parameters */
  gam = 5./3. ;

  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }

  fracphi = 1.;

  set_arrays() ;
  set_grid() ;
  
  coord(-2,0,0,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  fprintf(stderr,"rmin: %g\n",r) ;

  /* output choices */
  tf = 600.;

  DTd = 1. ;	/* dumping frequency, in units of M */
  DTl = 50. ;	/* logfile frequency, in units of M */
  DTi = 50. ; 	/* image file frequ., in units of M */
  DTr = 1. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */

  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;

  rhomax = 0. ;
  umax = 0. ;
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&r,&th,&phi) ;
    
    
    if(r<.05){
      rho = 3.857143;
      ur = 2.629369/1000.  ;
       u = (gam-1.)*pow(1./1000.,2.)*10.3333333;
    }
    else {
      ur = 0. ;
      rho = 1.+.2*sin(50.*r);
       u = (gam-1.)*pow(1./1000.,2.);
    }
    

    uh = 0. ;
    up = 0. ;
    
    p[i][j][k][RHO] = rho ;
    p[i][j][k][UU] = u ;
    p[i][j][k][U1] = ur ;
    p[i][j][k][U2] = uh ;
    p[i][j][k][U3] = up ;
    p[i][j][k][B1] = 0. ;
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
  }

  rhomax = 1. ;
  fixup(p) ;
  bound_prim(p) ;

  /* first find corner-centered vector potential */
  ZSLOOP(0,N1,0,N2,0,N3) A[i][j][k] = 0. ;


  /* now differentiate to find cell-centered B,
and begin normalization */
  bsq_max = 0. ;
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;

    /* flux-ct */
    p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
            + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
    p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
           - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;

    p[i][j][k][B3] = 0. ;

    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
  fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

  /* enforce boundary conditions */


  fixup(p) ;
  bound_prim(p) ;
  
  

  
#if( DO_FONT_FIX )
  set_Katm();
#endif
    
    
}
void init_vtran()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double rhor;
  double k_vec_x = 2. * M_PI;  //wavevector
  double k_vec_y = 0.;
  double k_vec_len = sqrt( k_vec_x * k_vec_x + k_vec_y * k_vec_y );
  double incl = M_PI/3;
  double upar, uperp;
  double ADVEL = 1./1000.;
  
#define TWOD (0)
  
  /* for disk interior */
  double l,rin,lnh,expm2chi,up1 ;
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;

  /* for magnetic field */
  double A[N1+1][N2+1][N3+1] ;
  double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double rmax;

  /* some physics parameters */
  gam = 5./3. ;

  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;
  
  set_arrays() ;
  set_grid() ;
  
  coord(-2,0,0,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  fprintf(stderr,"rmin: %g\n",r) ;

  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;

  rhomax = 0. ;
  umax = 0. ;

#if(TWOD)
  double elx,ely,nx,ny ;
  elx = 1. ;
  ely = 0.5 ;
  nx = 2. ;
  ny = 1. ;
  
  
  /* output choices */
  tf = 1./(ADVEL*ely) ;
  DTd = tf/10.;		/* dumping frequency, in units of M */
  DTl = tf/100.;		/* logfile frequency, in units of M */
  DTi = tf/100.;		/* image file frequ., in units of M */
  DTr = tf/10. ; /* restart file frequ., in units of M */
  DTr01 = 512 ; /* restart file frequ., in timesteps */
  
  
  ZLOOP {
    coord(i, j, CENT, X);
    bl_coord(X,&r,&th) ;
    
    p[i][j][k][RHO] = 1. ;
    p[i][j][k][UU] = 1.*ADVEL*ADVEL/(gam*(gam - 1.)) ;
    p[i][j][k][U1] = ADVEL*elx +
                ADVEL*(-ely)*cos(2.*M_PI*(nx*r + ny*th)) ;
    p[i][j][k][U2] = ADVEL*ely +
                ADVEL*(elx)*cos(2.*M_PI*(nx*r + ny*th)) ;
    p[i][j][k][U3] = 0.;
    p[i][j][k][B1] = 0.;
    p[i][j][k][B2] = 0.;
    p[i][j][k][B3] = 0.;
  }

#else
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    
    coord(i, j, k, CENT, X);
    bl_coord(X,&r,&th,&phi) ;
    
    p[i][j][k][RHO] = 1. ;
    p[i][j][k][UU] = 1.*ADVEL*ADVEL/(gam*(gam - 1.)) ;
    p[i][j][k][U1] = ADVEL ;
    p[i][j][k][U2] =ADVEL*cos(2.*M_PI*r) ;
    p[i][j][k][U3] = 0.;
    p[i][j][k][B1] = 0.;
    p[i][j][k][B2] = 0.;
    p[i][j][k][B3] = 0.;
    
    /* output choices */
    tf = 1./(ADVEL) ;
    DTd = tf/10.;		/* dumping frequency, in units of M */
    DTl = tf/100.;		/* logfile frequency, in units of M */
    DTi = tf/100.;		/* image file frequ., in units of M */
    DTr = tf/10. ; /* restart file frequ., in units of M */
    DTr01 = 512 ; /* restart file frequ., in timesteps */
  }
  
#endif

  rhomax = 1. ;
  fixup(p) ;
  bound_prim(p) ;

  /* first find corner-centered vector potential */
  ZSLOOP(0,N1,0,N2,0,N3) A[i][j][k] = 0. ;


  /* now differentiate to find cell-centered B,
and begin normalization */
  bsq_max = 0. ;
  ZLOOP {
    get_geometry(i,j,k,CENT,&geom) ;

    /* flux-ct */
    p[i][j][k][B1] = -(A[i][j][k] - A[i][j+1][k]
            + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*geom.g) ;
    p[i][j][k][B2] = (A[i][j][k] + A[i][j+1][k]
           - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*geom.g) ;

    p[i][j][k][B3] = 0. ;

    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
  fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

  /* enforce boundary conditions */


  fixup(p) ;
  bound_prim(p) ;
  
  
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
  
  
}



void init_modes()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double rhor, delphi,delkel, Tel,kpert, kdisc,phase;
  
  /* for disk interior */
  double l,rin,lnh,expm2chi,up1 ;
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;
  
  /* for magnetic field */
  double A[N1+1][N2+1][N3+1] ;
  double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double rmax;
  
  /* some physics parameters */
  gam = 5./3. ;
    
    
    
  
  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = .5 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  if(N2!=1) {
    //2D problem, use full pi-wedge in theta
    fractheta = 1.;
  }
  else{
    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
    fractheta = 1.e-2;
  }
  
  fracphi = 1.;

  set_arrays() ;
  set_grid() ;
  
  coord(-2,0,0,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  fprintf(stderr,"rmin: %g\n",r) ;
  
  /* output choices */
  tf = 2.*M_PI/1.0871304099;
  
  DTd = tf/100. ;	/* dumping frequency, in units of M */
  DTl = 50. ;	/* logfile frequency, in units of M */
  DTi = 50. ; 	/* image file frequ., in units of M */
  DTr = tf/100. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  
  rhomax = 0. ;
  umax = 0. ;
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&r,&th,&phi) ;
    
    ur = 0.;
    rho = 1. ;
    
    u = (gam-1.);
      
      delphi =0.361157559257*u*(game4-1.)/pow(rho,-game4)*.01 ;
      delkel = 0.93250480824*u*(game4-1.)/pow(rho,-game4)*.01;
      phase = -3.0956324000290;
      
      kdisc = 2.*M_PI;
      kpert =kdisc;///dx[1];
    
    uh = 0. ;
    up = 0. ;
    
    p[i][j][k][RHO] = rho ;
    p[i][j][k][UU] = u ;
    p[i][j][k][U1] = ur ;
    p[i][j][k][U2] = uh ;
    p[i][j][k][U3] = up ;
    p[i][j][k][B1] = 1. ;
    p[i][j][k][B2] = 2. ;
    p[i][j][k][B3] = 0. ;
    p[i][j][k][PHI] = delphi*sin(kpert*r+kpert*th+phase);
    p[i][j][k][KEL4] = u*(game4-1.)/pow(rho,-game4) *(1.+delkel*sin(kpert*r+kpert*th));
  }
  
  rhomax = 1. ;
  fixup(p) ;
  bound_prim(p) ;
  
  
  /* enforce boundary conditions */
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
    
  
}

void init_turb()
{
  int i,j,k ;
  double r,th,phi,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;
  struct of_geom geom ;
  double  delvmag;
  double cs, kx, ky, dk, kmag, kpeak, Anorm, sigv, kang;
  double u1r, u2r, zr, eps;
  double vxavg, vyavg;
  
  
  /* for magnetic field */
  double A[N1+1][N2+1][N3+1] ;
  double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
  double rmax;
  
  /* some physics parameters */
  gam = 5./3. ;
  
  /* some numerical parameters */
  lim = VANL ;
  failed = 0 ;	/* start slow */
  cour = 0.9 ;
  dt = 1.e-5 ;
  
  t = 0. ;
  hslope = 1. ;
  
  
  set_arrays() ;
  set_grid() ;
  
  coord(-2,0,0,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  fprintf(stderr,"rmin: %g\n",r) ;
  
  /* output choices */
    tf = 30000.;
  
	DTd = tf/20. ;	/* dumping frequency, in units of M */
  DTl = 50. ;	/* logfile frequency, in units of M */
	DTi = 500000. ; 	/* image file frequ., in units of M */
  DTr = tf/20. ; /* restart file frequ., in units of M */
  DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
  /* start diagnostic counters */
  dump_cnt = 0 ;
  image_cnt = 0 ;
  rdump_cnt = 0 ;
  rdump01_cnt = 0 ;
  defcon = 1. ;
  
  rhomax = 0. ;
  umax = 0. ;
  bsq_max = 0. ;
  ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    coord(i,j,k,CENT,X) ;
    bl_coord(X,&r,&th,&phi) ;
    
    
    rho = 1. ;
    
    u = (gam-1.)*rho/pow(1000.,2);
    cs =sqrt(gam * (gam-1) * u / rho);
    beta = 10.;
    
    ur = 0.   ;
    uh = 0.;
    up = 0. ;
    
    p[i][j][k][RHO] = rho ;
    p[i][j][k][UU] = u ;
    p[i][j][k][U1] = ur ;
    p[i][j][k][U2] = uh ;
    p[i][j][k][U3] = up ;
    p[i][j][k][B1] = sqrt(2.*cs*cs*rho/beta);
    p[i][j][k][B2] = 0. ;
    p[i][j][k][B3] = 0. ;
    
    bsq_ij = bsq_calc(p[i][j][k],&geom) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
  
  
  
  rhomax = 1. ;
  fixup(p) ;
  bound_prim(p) ;
  
  
  fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;
  
  /* enforce boundary conditions */
  fixup(p) ;
  bound_prim(p) ;
  
  
  
  
  
#if( DO_FONT_FIX )
  set_Katm();
#endif
  
  
}
void init_grad()
{
    int i,j,k ;
    double x,y,z,sth,cth ;
    double ur,uh,up,u,rho ;
    double X[NDIM] ;
    struct of_geom geom ;
    double rhor;
    double tfac = 1.e3; //factor by which to reduce velocity
    
    
    /* some physics parameters */
    gam = 5./3. ;
        
    /* some numerical parameters */
    lim = VANL ;
    failed = 0 ;	/* start slow */
    cour = 0.9 ;
    dt = 1.e-5 ;
    
    t = 0. ;
    hslope = 1. ;
    
    if(N2!=1) {
        //2D problem, use full pi-wedge in theta
        fractheta = 1.;
    }
    else{
        //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
        fractheta = 1.e-2;
    }
    
    set_arrays() ;
    set_grid() ;
    
    
    /* output choices */
    
    tf = 15.; //tfac*2./100.*5.;
    
    DTd = tf/100. ;	/* dumping frequency, in units of M */
    DTl = tf/10. ;	/* logfile frequency, in units of M */
    DTi = tf/10. ; 	/* image file frequ., in units of M */
    DTr = tf/100. ; /* restart file frequ., in units of M */
    DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;
    
    
    ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&x,&y,&z) ;
        
        
        p[i][j][k][RHO] = 1. ;//1./(1.+exp(-(sqrt(y*y+x*x)-.5)/.0001));
        p[i][j][k][UU] = pow(p[i][j][k][RHO],gam)/1000./1000.;
        p[i][j][k][U1] = 0.;
        p[i][j][k][U2] = 0.;
        p[i][j][k][U3] = 0 ;
        p[i][j][k][B1] = 0.;//1./(1.+exp(-(x-.5)/.0001))*p[i][j][k][UU]*1000. ;
        p[i][j][k][B2] = 1./(1.+exp(-(x-.5)/.0001))*p[i][j][k][UU]*10000.  ;
        p[i][j][k][B3] = 0. ;
    }
    
    /* enforce boundary conditions */
    

    
    
    
    fixup(p) ;
    bound_prim(p) ;
    
#if( DO_FONT_FIX )
    set_Katm();
#endif
    
    
}

void init_statcond()
{
    int i,j,k ;
    double x,y,z,sth,cth ;
    double ur,uh,up,u,rho ;
    double X[NDIM] ;
    struct of_geom geom ;
    double rhor, r;
    
    
    /* some physics parameters */
    gam = 4./3. ;
    
    /* some numerical parameters */
    lim = VANL ;
    failed = 0 ;	/* start slow */
    cour = 0.9 ;
    dt = 1.e-5 ;
    
    t = 0. ;
    hslope = 1. ;
    
    if(N2!=1) {
        //2D problem, use full pi-wedge in theta
        fractheta = 1.;
    }
    else{
        //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
        fractheta = 1.e-2;
    }
    
    set_arrays() ;
    set_grid() ;
    
    
    /* output choices */
    
    tf = 100.;
    
    DTd = tf/100. ;	/* dumping frequency, in units of M */
    DTl = tf/10. ;	/* logfile frequency, in units of M */
    DTi = tf/10. ; 	/* image file frequ., in units of M */
    DTr = tf/100. ; /* restart file frequ., in units of M */
    DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;
    
    
    ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&x,&y,&z) ;
        
        r = sqrt((x-.5)*(x-.5)+(y-.5)*(y-.5));
        
        p[i][j][k][RHO] = 1.-.9*exp(-(r*r)/.005);
        p[i][j][k][UU] = 10.;
        p[i][j][k][U1] = 0.;
        p[i][j][k][U2] = 0.;
        p[i][j][k][U3] = 0 ;
        p[i][j][k][B1] = .00000000001;
        p[i][j][k][B2] = .00000000001*sin(8.*M_PI*(x));
        p[i][j][k][B3] = 0. ;
    }
    
    /* enforce boundary conditions */
    
    
    
    
    
    fixup(p) ;
    bound_prim(p) ;
    
#if( DO_FONT_FIX )
    set_Katm();
#endif
}
void init_statcond1D()
{
    int i,j,k ;
    double x,y,z,sth,cth ;
    double ur,uh,up,u,rho ;
    double X[NDIM] ;
    struct of_geom geom ;
    double rhor, r;
    
    
    /* some physics parameters */
    gam = 4./3. ;
    
    /* some numerical parameters */
    lim = VANL ;
    failed = 0 ;	/* start slow */
    cour = 0.9 ;
    dt = 1.e-5 ;
    
    t = 0. ;
    hslope = 1. ;
    
    if(N2!=1) {
        //2D problem, use full pi-wedge in theta
        fractheta = 1.;
    }
    else{
        //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
        fractheta = 1.e-2;
    }
    
    set_arrays() ;
    set_grid() ;
    
    
    /* output choices */
    
    tf = 1000.;
    
    DTd = tf/100. ;	/* dumping frequency, in units of M */
    DTl = tf/10. ;	/* logfile frequency, in units of M */
    DTi = tf/10. ; 	/* image file frequ., in units of M */
    DTr = tf/100. ; /* restart file frequ., in units of M */
    DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;
    
    
    ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&x,&y,&z) ;
        
        r = sqrt((x-.5)*(x-.5)+(y-.5)*(y-.5));
        
        p[i][j][k][RHO] = 1.-.9999999*exp(-(r*r)/.005);
        p[i][j][k][UU] = 10./1000./1000.;
        p[i][j][k][U1] = 0.;
        p[i][j][k][U2] = 0.;
        p[i][j][k][U3] = 0 ;
        p[i][j][k][B1] = .00000000001;
        p[i][j][k][B2] = 0;
        p[i][j][k][B3] = 0. ;
    }
    
    /* enforce boundary conditions */
    
    
    
    
    
    fixup(p) ;
    bound_prim(p) ;
    
#if( DO_FONT_FIX )
    set_Katm();
#endif
}
void init_Hubble()
{
    int i,j,k,m ;
    double x,y,z,sth,cth ;
    double ur,uh,up,u,rho ;
    double X[NDIM] ;
    struct of_geom geom ;
    double rhor, r;
    double v0;
    
    
    /* some physics parameters */
    gam = 5./3. ;
    
    /* some numerical parameters */
    lim = VANL ;
    failed = 0 ;	/* start slow */
    cour = 0.9 ;
    dt = 1.e-5 ;
    
    t = 0. ;
    hslope = 1. ;
    
    if(N2!=1) {
        //2D problem, use full pi-wedge in theta
        fractheta = 1.;
    }
    else{
        //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
        fractheta = 1.e-2;
    }
    
    set_arrays() ;
    set_grid() ;
    
    
    /* output choices */
    
    tf = 1000.;
    
    DTd = tf/100. ;	/* dumping frequency, in units of M */
    DTl = tf/10. ;	/* logfile frequency, in units of M */
    DTi = tf/10. ; 	/* image file frequ., in units of M */
    DTr = tf/100. ; /* restart file frequ., in units of M */
    DTr01 = 1000 ; /* restart file frequ., in timesteps */
  
    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;
    
    
    v0 =1./1000.;
    
    ZSLOOP(-N1G,N1-1+N1G,-N2G,N2-1+N2G,0,N3-1) {
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&x,&y,&z) ;
        
        p[i][j][k][RHO] = 1.;
        p[i][j][k][UU] = 1./1000./1000.;
        p[i][j][k][U1] = v0*(x-.5);
        p[i][j][k][U2] = 0.;
        p[i][j][k][U3] = 0 ;
        p[i][j][k][B1] = 0.;
        p[i][j][k][B2] = 0;
        p[i][j][k][B3] = 0. ;
    }
    
    /* enforce boundary conditions */
    
    
    
    PLOOP{
        
        ZSLOOP(-N1G,N1+N1G-1,-N2G,N2-1+N2G,0,N3-1){
            
            pbound[i][j][k][m]=p[i][j][k][m];
            
        }
    }
    
    
    fixup(p) ;
    bound_prim(p) ;
    
#if( DO_FONT_FIX )
    set_Katm();
#endif
}



void init_atm()
{
    int i,j,k,m;
    double r,th,phi,sth,cth ;
    double ur,uh,up,u ;
    double X[NDIM] ;
    struct of_geom geom ;
    
    /* for disk interior */
    double l,rin,lnh,expm2chi,up1 ;
    double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
    double kappa,hm1 ;
    
    /* for magnetic field */
    double A[N1+2*NG][N2+2*NG][N3+2*NG] ;
    double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
    double rmax, lfish_calc(double rmax) ;
    
    /* some physics parameters */
    gam = 4./3. ;
    
    /* disk parameters (use fishbone.m to select new solutions) */
    a = 0. ;
    
    kappa =1.e-3;
    beta = 1.e2 ;
    
    /* some numerical parameters */
    lim = MC ;
    failed = 0 ;	/* start slow */
    cour = 0.5 ;
    dt = 1.e-5 ;
    R0 = 0.0 ;
    Rin = exp(1.)/2.*(1. + sqrt(1. - a*a)) ;  //.98
    Rout = exp(4.5) ;
    rbr = 5000.;
    npow2=4.0; //power exponent
    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    
    t = 0. ;
    hslope = 0.3 ;
    
    if(N2!=1) {
        //2D problem, use full pi-wedge in theta
        fractheta = 1.;
    }
    else{
        //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
        fractheta = 1.e-2; //.998;
    }
    
    fracphi = 1.;
    
    set_arrays() ;
    set_grid() ;
    
    coord(-2,0,0,CENT,X) ;
    bl_coord(X,&r,&th,&phi) ;
    fprintf(stderr,"rmin: %g\n",r) ;
    fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;
    
    /* output choices */
    tf = 100.0 ;
    
    DTd = tf/1.; //10.; /* dumping frequency, in units of M */
    DTl = 2. ;	/* logfile frequency, in units of M */
    DTi = 20000. ; 	/* image file frequ., in units of M */
    DTr = tf/1.; /* restart file frequ., in units of M */
    DTr01 = 100 ; /* restart file frequ., in timesteps */
  
    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;
    
    rhomax = 0. ;
    umax = 0. ;
    
    FILE *atmosphereSolnRHO, *atmosphereSolnUU, *atmosphereSolnPHI;
    FILE *atmosphereSolnRCoords, *atmosphereSolnUe;
    atmosphereSolnRHO = fopen("atmosphere_soln_rho.txt", "r");
    atmosphereSolnUU = fopen("atmosphere_soln_u.txt", "r");
    atmosphereSolnPHI = fopen("atmosphere_soln_phi.txt", "r");
    atmosphereSolnRCoords = fopen("atmosphere_soln_rCoords.txt", "r");
    
    char *rhoLine = NULL, *uLine = NULL, *rLine = NULL, *phiLine = NULL;
    size_t rhoLen=0; ssize_t rhoRead;
    size_t uLen=0; ssize_t uRead;
    size_t rLen=0; ssize_t rRead;
    size_t phiLen=0; ssize_t phiRead;
    
    double rho[N1+2*NG], uu[N1+2*NG], phicond[N1+2*NG], rCoords[N1+2*NG];
    
    for ( i=-NG; i<N1+NG; i++) {
        rhoRead = getline(&rhoLine, &rhoLen, atmosphereSolnRHO);
        uRead = getline(&uLine, &uLen, atmosphereSolnUU);
        phiRead = getline(&phiLine, &phiLen, atmosphereSolnPHI);
        rRead = getline(&rLine, &rLen, atmosphereSolnRCoords);

        
        rho[i+NG] = atof(rhoLine);
        uu[i+NG] = atof(uLine);
        phicond[i+NG] = atof(phiLine);
        rCoords[i+NG] = atof(rLine);
    }
    
    free(rhoLine); free(uLine); free(rLine); free(phiLine);
    fclose(atmosphereSolnRHO);
    fclose(atmosphereSolnUU);
    fclose(atmosphereSolnPHI);
    fclose(atmosphereSolnRCoords);
    
    ZSLOOP(-N1G,N1+N1G-1,-N2G,N2-1+N2G,0,N3-1) {
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&r,&th,&phi) ;
        get_geometry(i,j,k,CENT,&geom);
        double a1 = geom.gcov[1][1];
        double b1 = geom.gcon[0][1];
        double c1 = geom.gcon[0][0];
        double v1 = (c1*0./r - sqrt(-a1*b1*b1*b1*b1 -a1*b1*b1*c1*0./(r*r) - b1*b1*c1))/(a1*b1*b1 + c1);

        
        if (abs(r - rCoords[i+NG])>1e-15)
        {
            fprintf(stderr,"r = %f, rCoords = %f, i= %d, DX1 = %f\n", r,
                        rCoords[i+NG],i, dx[1]);
            fprintf(stderr, "Mismatch in rCoords! Check r coords in python script\n");
            exit(1);
        }
        
        p[i][j][k][RHO] = rho[i+NG] ;
        p[i][j][k][UU] = uu[i+NG]; ;
        p[i][j][k][U1] = 0 ;
        p[i][j][k][U2] = 0. ;
        p[i][j][k][U3] = 0. ;
        
        p[i][j][k][B1] = 0. ;
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
        p[i][j][k][PHI] = phicond[i+NG];
        p[i][j][k][KEL4] = uu[i+NG]*(game4-1.)/pow(rho[i+NG],game4);
        /* convert from 4-vel to 3-vel */
        coord_transform(p[i][j][k],i,j,k);
        
    
    }
    
    PLOOP{
        
        ZSLOOP(-N1G,N1+N1G-1,-N2G,N2-1+N2G,0,N3-1){
    
            pbound[i][j][k][m]=p[i][j][k][m];
            
        }
    }
    
    
    fixup(p) ;
    bound_prim(p) ;
    
        
    ZSLOOP(-N1G,N1+NG-1,-NG,N2-1+NG,0,N3-1) {
         coord(i,j,k,CORN,X) ;
         bl_coord(X,&r,&th,&phi) ;
        A[i+NG][j+NG][k] = .001*(1-cos(th));
        
    }
    /* now differentiate to find cell-centered B,
     and begin normalization */
    ZSLOOP(-N1G+1,N1+N1G-2,0,N2-1,0,N3-1) {
        get_geometry(i,j,k,CENT,&geom) ;
        
        /* flux-ct */
        p[i][j][k][B1] = -(A[i+NG][j+NG][k] - A[i+NG][j+NG+1][k]
                           + A[i+NG+1][j+NG][k] - A[i+1+NG][j+1+NG][k])/(2.*dx[2]*geom.g) ;
        p[i][j][k][B2] = (A[i+NG][j+NG][k] + A[i+N1G][j+1+NG][k]
                          - A[i+1+NG][j+NG][k] - A[i+1+N1G][j+1+NG][k])/(2.*dx[1]*geom.g) ;
        
        p[i][j][k][B3] = 0. ;
        
            }
    
    PLOOP{
        
        ZSLOOP(-N1G,N1+N1G-1,-N2G,N2-1+N2G,0,N3-1){
            
            pbound[i][j][k][m]=p[i][j][k][m];
            
        }
    }

    
    
    /* enforce boundary conditions */
    fixup(p) ;
    bound_prim(p) ;
    
    
        
    
    
    
#if( DO_FONT_FIX )
    set_Katm();
#endif 
    
    
}

void init_bondicon()
{
    int i,j,k,m;
    double r,th,phi,sth,cth ;
    double uh,up,u ;
    double X[NDIM] ;
    struct of_geom geom ;
    
    /* for disk interior */
    double l,rin,lnh,expm2chi,up1 ;
    double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
    double kappa,hm1 ;
    
    /* for magnetic field */
    double A[N1+1][N2+1][N3+1] ;
    double rho_av,rhomax,umax,beta,bsq_ij,bsq_max,norm,q,beta_act ;
    double rmax, lfish_calc(double rmax) ;
    
    /* some physics parameters */
    gam = 4./3. ;
    
    /* disk parameters (use fishbone.m to select new solutions) */
    a = 0.;
    
    kappa =1.e-3;
    beta = 1.e2 ;
    
    /* some numerical parameters */
    lim = MC ;
    failed = 0 ;	/* start slow */
    cour = 0.9 ;
    dt = 1.e-5 ;
    R0 = 0. ;
    Rin = .8*(1. + sqrt(1. - a*a)) ;
    Rout = 40. ;
    rbr = 5000000.;
    npow2=4.0; //power exponent
    cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)
    
    t = 0. ;
    hslope = 0.3 ;
    
    if(N2!=1) {
        //2D problem, use full pi-wedge in theta
        fractheta = 1.;
    }
    
    else{
        //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
        fractheta = 1.e-2;
    }
    
    fracphi = 1.;
    
    set_arrays() ;
    set_grid() ;
    
    coord(-2,0,0,CENT,X) ;
    bl_coord(X,&r,&th,&phi) ;
    fprintf(stderr,"rmin: %g\n",r) ;
    fprintf(stderr,"rmin/rm: %g\n",r/(1. + sqrt(1. - a*a))) ;
    
    /* output choices */
    tf = 200.0 ;
    
    DTd = tf/1.; //10.; /* dumping frequency, in units of M */
    DTl = 2. ;	/* logfile frequency, in units of M */
    DTi = 20000. ; 	/* image file frequ., in units of M */
    DTr = tf/1.; /* restart file frequ., in units of M */
    DTr01 = 100 ; /* restart file frequ., in timesteps */
  
    /* start diagnostic counters */
    dump_cnt = 0 ;
    image_cnt = 0 ;
    rdump_cnt = 0 ;
    rdump01_cnt = 0 ;
    defcon = 1. ;
    
    rhomax = 0. ;
    umax = 0. ;
    
    FILE *bondiSolnRHO, *bondiSolnUU, *bondiSolnU1, *bondiSolnRCoords;
    FILE *bondiSolnPHI;
    bondiSolnRHO = fopen("bondi_soln_rho.txt", "r");
    bondiSolnUU = fopen("bondi_soln_u.txt", "r");
    bondiSolnU1 = fopen("bondi_soln_ur.txt", "r");
    bondiSolnPHI = fopen("bondi_soln_phi_no_backreaction.txt", "r");
    bondiSolnRCoords = fopen("bondi_soln_rCoords.txt", "r");
    
    char *rhoLine = NULL, *uLine = NULL, *rLine = NULL, *urLine = NULL;
    char *phiLine = NULL;
    size_t rhoLen=0; ssize_t rhoRead;
    size_t uLen=0; ssize_t uRead;
    size_t urLen=0; ssize_t urRead;
    size_t phiLen=0; ssize_t phiRead;
    size_t rLen=0; ssize_t rRead;
    
    double rho[N1+2*NG], uu[N1+2*NG], ur[N1+2*NG], rCoords[N1+2*NG], phicond[N1+2*NG];
    
    for ( i=-NG; i<N1+NG; i++) {
        rhoRead = getline(&rhoLine, &rhoLen, bondiSolnRHO);
        uRead = getline(&uLine, &uLen, bondiSolnUU);
        urRead = getline(&urLine, &urLen, bondiSolnU1);
        phiRead = getline(&phiLine, &phiLen, bondiSolnPHI);
        rRead = getline(&rLine, &rLen, bondiSolnRCoords);
        
        rho[i+NG] = atof(rhoLine);
        uu[i+NG] = atof(uLine);
        ur[i+NG] = atof(urLine);
        phicond[i+NG] = atof(phiLine);
        rCoords[i+NG] = atof(rLine);
    }
    
    free(rhoLine); free(uLine); free(urLine); free(phiLine);
    fclose(bondiSolnRHO);
    fclose(bondiSolnUU);
    fclose(bondiSolnU1);
    fclose(bondiSolnPHI);
    fclose(bondiSolnRCoords);
    
    
    
    
    ZSLOOP(-NG,N1+NG-1,0,N2-1,0,N3-1) {
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&r,&th,&phi) ;
        get_geometry(i,j,k,CENT,&geom);
        double a = geom.gcov[1][1];
        double b = geom.gcon[0][1];
        double c = geom.gcon[0][0];
        double v1 = (c*ur[i+NG]/r - sqrt(-a*b*b*b*b -a*b*b*c*ur[i+NG]*ur[i+NG]/(r*r) - b*b*c))/(a*b*b + c);

        
        if (abs(r - rCoords[i+NG])>1e-15)
        {
            fprintf(stderr,"r = %f, rCoords = %f, i= %d, DX1 = %f\n", r,
                    rCoords[i+NG],i, dx[1]);
            fprintf(stderr, "Mismatch in rCoords! Check r coords in python script\n");
            exit(1);
        }
        
        p[i][j][k][RHO] = rho[i+NG] ;
        p[i][j][k][UU] = uu[i+NG]; ;
        p[i][j][k][U1] = v1;
        p[i][j][k][U2] = 0. ;
        p[i][j][k][U3] = 0. ;
        
        double qB = 0.0000000001;
        p[j][i][k][B1] = qB/(r*r*r);
        p[i][j][k][B2] = 0. ;
        p[i][j][k][B3] = 0. ;
        p[i][j][k][PHI] = phicond[i+NG];
        p[i][j][k][KEL4] = uu[i+NG]*(game4-1.)/pow(rho[i+NG],game4);

        
        /* convert from 4-vel to 3-vel */
        //coord_transform(p[i][j][k],i,j,k);
    }
    
    PLOOP{
        
        ZSLOOP(-NG,N1+NG-1,0,N2-1,0,N3-1){
            
            pbound[i][j][k][m]=p[i][j][k][m];
            
        }
    }

    /* enforce boundary conditions */
    fixup(p) ;
    bound_prim(p) ;
    
    
    
    PLOOP{
        
        ZSLOOP(-NG,N1+NG-1,0,N2-1,0,N3-1){
            
            pbound[i][j][k][m]=p[i][j][k][m];
            
        }
    }
    
    
    
#if( DO_FONT_FIX )
    set_Katm();
#endif
    
    
}



/* this version starts w/ BL 4-velocity and
 * converts to relative 4-velocities in modified
 * Kerr-Schild coordinates */

void coord_transform(double *pr,int ii, int jj, int kk)
{
  double X[NDIM],r,th,phi,ucon[NDIM],uconp[NDIM],trans[NDIM][NDIM],tmp[NDIM] ;
  double AA,BB,CC,discr ;
  double utconp[NDIM], dxdxp[NDIM][NDIM], dxpdx[NDIM][NDIM] ;
  struct of_geom geom ;
  struct of_state q ;
  int i,j,k,m ;

  coord(ii,jj,kk,CENT,X) ;
  bl_coord(X,&r,&th,&phi) ;
  blgset(ii,jj,kk,&geom) ;

  ucon[1] = pr[U1] ;
  ucon[2] = pr[U2] ;
  ucon[3] = pr[U3] ;

  AA =     geom.gcov[TT][TT] ;
  BB = 2.*(geom.gcov[TT][1]*ucon[1] +
           geom.gcov[TT][2]*ucon[2] +
           geom.gcov[TT][3]*ucon[3]) ;
  CC = 1. +
          geom.gcov[1][1]*ucon[1]*ucon[1] +
          geom.gcov[2][2]*ucon[2]*ucon[2] +
          geom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(geom.gcov[1][2]*ucon[1]*ucon[2] +
          geom.gcov[1][3]*ucon[1]*ucon[3] +
          geom.gcov[2][3]*ucon[2]*ucon[3]) ;

  discr = BB*BB - 4.*AA*CC ;
  ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;
  /* now we've got ucon in BL coords */

  /* transform to Kerr-Schild */
  /* make transform matrix */
  DLOOP trans[j][k] = 0. ;
  DLOOPA trans[j][j] = 1. ;
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a) ;
  trans[3][1] = a/(r*r - 2.*r + a*a) ;

  /* transform ucon */
  DLOOPA tmp[j] = 0. ;
  DLOOP tmp[j] += trans[j][k]*ucon[k] ;
  DLOOPA ucon[j] = tmp[j] ;
  /* now we've got ucon in KS coords */

  /* transform to KS' coords */
  /* dr^\mu/dx^\nu jacobian, where x^\nu are internal coords */
  dxdxp_func(X, dxdxp);
  /* dx^\mu/dr^\nu jacobian */
  invert_matrix(dxdxp, dxpdx);
  
  for(i=0;i<NDIM;i++) {
    uconp[i] = 0;
    for(j=0;j<NDIM;j++){
      uconp[i] += dxpdx[i][j]*ucon[j];
    }
  }
  //old way of doing things for Gammie coords
  //ucon[1] *= (1./(r - R0)) ;
  //ucon[2] *= (1./(M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]))) ;
  //ucon[3] *= 1.; //!!!ATCH: no need to transform since will use phi = X[3]

  get_geometry(ii, jj, kk, CENT, &geom);
  
  /* now solve for relative 4-velocity that is used internally in the code:
   * we can use the same u^t because it didn't change under KS -> KS' */
  ucon_to_utcon(uconp,&geom,utconp);
  
  pr[U1] = utconp[1] ;
  pr[U2] = utconp[2] ;
  pr[U3] = utconp[3] ;

  /* done! */
}


double lfish_calc(double r)
{
	return(
   ((pow(a,2) - 2.*a*sqrt(r) + pow(r,2))*
      ((-2.*a*r*(pow(a,2) - 2.*a*sqrt(r) + pow(r,2)))/
         sqrt(2.*a*sqrt(r) + (-3. + r)*r) +
        ((a + (-2. + r)*sqrt(r))*(pow(r,3) + pow(a,2)*(2. + r)))/
         sqrt(1 + (2.*a)/pow(r,1.5) - 3./r)))/
    (pow(r,3)*sqrt(2.*a*sqrt(r) + (-3. + r)*r)*(pow(a,2) + (-2. + r)*r))
	) ;
}


