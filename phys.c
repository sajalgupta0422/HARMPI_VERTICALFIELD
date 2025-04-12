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

/***********************************************************************************************/
/***********************************************************************************************
 primtoflux():
 ---------
 --  calculate fluxes in direction dir, The fluxes are defined at zone faces.
 
 ***********************************************************************************************/

void primtoflux(double *pr, struct of_state *q, int dir,
		struct of_geom *geom, double *flux)
{
	int j,k, m ;
	double mhd[NDIM], EE4[NDIM], EE5[NDIM], EUAD[NDIM], EEVAR[NDIM];

	/* particle number flux */ 
	flux[RHO] = pr[RHO]*q->ucon[dir] ; //q is passed as a pointer to struct of_state which contains various state variables including the four-velocities ( ucon,  ucov), magnetic field four-vectors ( bcon,  bcov)

	mhd_calc(pr, dir, q, mhd) ; //calculating mhd stress tensor = T^μ_ν, and stores in mhd variable

	/* MHD stress-energy tensor w/ first index up, 
	 * second index down. 
     mhd[0] = Tᵈ₀: Energy flux component in dir = d
    mhd[1] = Tᵈ₁: Radial momentum flux component in dir = d
    mhd[2] = Tᵈ₂: Theta momentum flux component in dir = d 
    mhd[3] = Tᵈ₃: Phi momentum flux component in dir = d
   */
	flux[UU] = mhd[0] + flux[RHO] ;
	flux[U1] = mhd[1] ;
	flux[U2] = mhd[2] ;
	flux[U3] = mhd[3] ;

	/* dual of Maxwell tensor (Faraday tensor)
  F*ᵃᵈ = bᵃuᵈ - bᵈuᵃ

  where:
- F*ᵃᵈ is the dual of Faraday tensor
- bᵃ is the magnetic four-vector
- uᵈ is the four-velocity
- a = 1,2,3 (spatial components)
- d = dir (direction of flux) 

  B^a = F*^at = b^a u^t - b^t u^a; 
  so, following terms are not B^a, but are actually fluxes of B^a defined as F*^{ad}

  Fluxes are basically terms that comes induction eq.
  \partial B/\partial t = Curl(v x B) (for η=0)
  In GR: Induction equation: ∂ₜ(*F^{it}) = -∂ᵢ(*F^{ij}) (due to anti-symmetric nature of tensor)
  ∂ₜ(√(-g)B^i) = -∂ⱼ[√(-g)(b^j u^i - b^i u^j)]

  so for i =r, ∂ₜ(√(-g)B^r) = ∂ⱼ[√(-g)(b^r u^j - b^j u^r)]
  In terms of the vector potential A_μ (when using constrained transport):
  ∂ₜA_i = ϵᵢⱼₖ E^k
  */
	flux[B1]  = q->bcon[1]*q->ucon[dir] - q->bcon[dir]*q->ucon[1] ; // F*¹ᵈ; Radial component of magnetic field flux
	flux[B2]  = q->bcon[2]*q->ucon[dir] - q->bcon[dir]*q->ucon[2] ; // F*²ᵈ; Theta component of magnetic field flux
	flux[B3]  = q->bcon[3]*q->ucon[dir] - q->bcon[dir]*q->ucon[3] ; // F*³ᵈ; Phi component of magnetic field flux

#if(DOKTOT || eHEAT || eCOND)
        flux[KTOT] = flux[RHO]*pr[KTOT] ; //heat flux
#endif
	
#if( eHEAT || eCOND)
	/* Flux of Entropy*/
    flux[KELDIS] = flux[RHO]*pr[KELDIS];
    flux[KEL4] = flux[RHO]*pr[KEL4];
    flux[KEL4A] = flux[RHO]*pr[KEL4A];
    flux[KEL4B] = flux[RHO]*pr[KEL4B];
    flux[KEL4C] = flux[RHO]*pr[KEL4C];
    flux[KEL4D] = flux[RHO]*pr[KEL4D];
    flux[KEL4E] = flux[RHO]*pr[KEL4E];
    flux[KEL5] = flux[RHO]*pr[KEL5];
    flux[PHI] = pr[PHI]*flux[RHO];
#if(DOFLR)
    flux[FLR] = flux[RHO]*pr[FLR];
#endif
#endif

#if(DONUCLEAR)
  flux[RHONP] = flux[RHO]*pr[RHONP];
  flux[RHOALPHA] = flux[RHO]*pr[RHOALPHA];
  flux[RHOFLOOR] = flux[RHO]*pr[RHOFLOOR];
  flux[YE] = flux[RHO]*pr[YE];
#endif
  
  PLOOP flux[m] *= geom->g ; //multiplying by metric factor
  
  
}

/* calculate "conserved" quantities; provided strictly for
 * historical reasons 
Conserved quantities are U = {\sqrt{-g} (ρu^t, T^t_t , T^t_i , B^i)}
Conserved variables are updates using fluxes
 */
void primtoU(double *pr, struct of_state *q, struct of_geom *geom, double *U)
{
  
  primtoflux(pr,q,0,geom, U) ;
  return ;
}

/* calculate magnetic field four-vector 
b^t = B^i u^μ g_{i μ} 
b^i =(B^i + b^t u^i) / u^t
*/
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon)
{
  int j ;
  
  bcon[TT] = pr[B1]*ucov[1] + pr[B2]*ucov[2] + pr[B3]*ucov[3] ;
  for(j=1;j<4;j++)
    bcon[j] = (pr[B1-1+j] + bcon[TT]*ucon[j])/ucon[TT] ;
  
  return ;
}

/* MHD stress tensor, with first index up, second index down */
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)
{
  int j ;
  double r,u,P,w,bsq,eta,ptot ;
  
  r = pr[RHO] ;
  u = pr[UU] ; //internal energy
  P = (gam - 1.)*u ;
  w = P + r + u ;
  bsq = dot(q->bcon,q->bcov) ;
  eta = w + bsq ;
  ptot = P + 0.5*bsq;
  
  /*
  mhd[j] = Tᵈⱼ (where d = dir) is:

  Tᵈⱼ = (ρh + b²)uᵈuⱼ + (P + ½b²)δᵈⱼ - bᵈbⱼ

  where:
  - ρ = rest-mass density (r)
  - h = specific enthalpy (w/ρ)
  - b² = magnetic field squared (bsq)
  - uᵈ = four-velocity (q->ucon[dir])
  - uⱼ = four-velocity covariant (q->ucov[j])
  - P = pressure
  - δᵈⱼ = Kronecker delta
  - bᵈ = magnetic four-vector (q->bcon[dir])
  - bⱼ = magnetic four-vector covariant (q->bcov[j])
  */
  /* single row of mhd stress tensor,
   * first index up, second index down */
  DLOOPA mhd[j] = eta*q->ucon[dir]*q->ucov[j]
  + ptot*delta(dir,j) - q->bcon[dir]*q->bcov[j] ;
  
}

/* add in source terms to equations of motion */
void source(double *ph, struct of_geom *geom, int ii, int jj, int kk, double *dU,
            double Dt)
{
  double mhd[NDIM][NDIM], EE4[NDIM][NDIM],EE5[NDIM][NDIM], EUAD[NDIM][NDIM], EEVAR[NDIM][NDIM];
  int j,k,m ;
  struct of_state q ;
  
  get_state(ph, geom, &q) ;
  mhd_calc(ph, 0, &q, mhd[0]) ;
  mhd_calc(ph, 1, &q, mhd[1]) ;
  mhd_calc(ph, 2, &q, mhd[2]) ;
  mhd_calc(ph, 3, &q, mhd[3]) ;
  
  
  
  
  /* contract mhd stress tensor with connection 
  T^μ_ν ζ^ν_{κ μ}*/
  PLOOP dU[m] = 0. ;
  DLOOP {
    dU[UU] += mhd[j][k]*conn[ii][jj][kk][k][0][j] ;
    dU[U1] += mhd[j][k]*conn[ii][jj][kk][k][1][j] ;
    dU[U2] += mhd[j][k]*conn[ii][jj][kk][k][2][j] ;
    dU[U3] += mhd[j][k]*conn[ii][jj][kk][k][3][j] ;

  }


  
  
  
  
  //misc_source(ph, ii, jj, geom, &q, dU, Dt) ;
  
  PLOOP dU[m] *= geom->g ; //multiplying by metric factor
  
  /* done! */
}

/* returns b^2 (i.e., twice magnetic pressure) */
double bsq_calc(double *pr, struct of_geom *geom)
{
  // Calculates B^2 using centered values
  // Shows magnetic energy density computed at centers
  // despite fields being face-centered
  struct of_state q ;
  
  get_state(pr,geom,&q) ;
  return( dot(q.bcon,q.bcov) ) ;
}

/* find ucon, ucov, bcon, bcov from primitive variables */
void get_state(double *pr, struct of_geom *geom, struct of_state *q)
{
  // Converts primitives to four-vectors
  // Handles both centered and face-centered quantities
  // in a unified way
  
  /* get ucon */
  ucon_calc(pr, geom, q->ucon) ;
  lower(q->ucon, geom, q->ucov) ;
  bcon_calc(pr, q->ucon, q->ucov, q->bcon) ;
  lower(q->bcon, geom, q->bcov) ;
  
  return ;
}

/* find relative 4-velocity from 4-velocity (both in code coords) 
The 4-velocity relative to this ZAMO is ũ^μ = u^μ −γ η^μ ; 
η^μ = (1/α,-β^i/α) is 4-contra velocity of ZAMO observer
In KS spacetime, β^r = g^{tr}/α^2; but β^θ = 0; β^φ = 0 
It is also comes with a spatial drift that is unique by always being related 
to the physical observer for any space-time*/
void ucon_to_utcon(double *ucon,struct of_geom *geom, double *utcon)
{
  double alpha, beta[NDIM], gamma;
  int j;
  
  /* now solve for v-- we can use the same u^t because
   * it didn't change under KS -> KS' 
   SLOOPA for(j=1;j<NDIM;j++) -- looping over spatial dimensions*/
  alpha = 1./sqrt(-geom->gcon[TT][TT]) ; //time lapse factor; α = 1/sqrt(-g^tt)
  SLOOPA beta[j] = geom->gcon[TT][j]*alpha*alpha ; //β^μ = g^tμ α^2
  gamma = alpha*ucon[TT] ; // γ = -u^α η_α; where η_μ = (-α,0,0,0) is 4-velocity of ZAMO observer 
  
 
  utcon[0] = 0;
  /*
  //(gamma*beta[j]/alpha) is the additional term represents the spatial drift
 of the ZAMO frame four velocity: η_μ = (-α,0,0,0)
  */
  SLOOPA utcon[j] = ucon[j] + gamma*beta[j]/alpha ; 

}

/* find contravariant four-velocity */
void ucon_calc(double *pr, struct of_geom *geom, double *ucon)
{
  double alpha,gamma ;
  double beta[NDIM] ;
  int j ;
  
  alpha = 1./sqrt(-geom->gcon[TT][TT]) ; //time lapse factor; α = 1/sqrt(-g^tt)
  SLOOPA beta[j] = geom->gcon[TT][j]*alpha*alpha ; //β^μ = g^tμ α^2
  
  if( gamma_calc(pr,geom,&gamma) ) { //γ = \sqrt{1 + q^2}; q^2 = g_{ij} \tild{u}^i \tild{u}^j
    fflush(stderr);
    fprintf(stderr,"\nucon_calc(): gamma failure \n");
    fflush(stderr);
    fail(FAIL_GAMMA);
  }
  
  ucon[TT] = gamma/alpha ; //u^t = γ/α
  SLOOPA ucon[j] = pr[U1+j-1] - gamma*beta[j]/alpha ; //getting u^i = \tild{u}^i + γη^i, where η^μ = (1/α,-β^i/α) is 4-contra velocity of ZAMO observer
  
  return ;
}

void ut_calc_3vel(double *vcon, struct of_geom *geom, double *ut)
{
  double AA, BB, CC, DD, one_over_alpha_sq;
  //compute the Lorentz factor based on contravariant 3-velocity
  AA =     geom->gcov[TT][TT] ;
  BB = 2.*(geom->gcov[TT][1]*vcon[1] +
           geom->gcov[TT][2]*vcon[2] +
           geom->gcov[TT][3]*vcon[3]) ;
  CC = geom->gcov[1][1]*vcon[1]*vcon[1] +
  geom->gcov[2][2]*vcon[2]*vcon[2] +
  geom->gcov[3][3]*vcon[3]*vcon[3] +
  2.*(geom->gcov[1][2]*vcon[1]*vcon[2] +
      geom->gcov[1][3]*vcon[1]*vcon[3] +
      geom->gcov[2][3]*vcon[2]*vcon[3]) ;
  
  DD = -1./(AA+BB+CC);
  
  one_over_alpha_sq = -geom->gcon[TT][TT];

  if(DD<one_over_alpha_sq) {
    DD = one_over_alpha_sq;
  }
  
  *ut = sqrt(DD);

}

/* find gamma-factor wrt normal observer 
γ = \sqrt{1 + q^2}; q^2 = g_{ij} \tilde{u}^i v\tilde{u}^j, 
where \tilde{u}^i = Spatial drift velocity relative to ZAMO
*/
int gamma_calc(double *pr, struct of_geom *geom, double *gamma)
{
  double qsq ;
  
  qsq =     geom->gcov[1][1]*pr[U1]*pr[U1]
  + geom->gcov[2][2]*pr[U2]*pr[U2]
  + geom->gcov[3][3]*pr[U3]*pr[U3]
  + 2.*(geom->gcov[1][2]*pr[U1]*pr[U2]
        + geom->gcov[1][3]*pr[U1]*pr[U3]
        + geom->gcov[2][3]*pr[U2]*pr[U3]) ;
  
  if( qsq < 0. ){
    if( fabs(qsq) > 1.E-10 ){ // then assume not just machine precision
      fprintf(stderr,"gamma_calc():  failed: i,j,k,qsq = %d %d %d %28.18e \n", icurr,jcurr,kcurr,qsq);
      fprintf(stderr,"v[1-3] = %28.18e %28.18e %28.18e  \n",pr[U1],pr[U2],pr[U3]);
      *gamma = 1.;
      return (1);
    }
    else qsq=1.E-10; // set floor
  }
  
  *gamma = sqrt(1. + qsq) ;
  
  return(0) ;
}


/*
 * VCHAR():
 *
 * calculate components of magnetosonic velocity (characteristic wave speeds in relativistic MHD.)
 * corresponding to primitive variables p
 *
 * cfg 7-10-01
 *
 */

void vchar(double *pr, struct of_state *q, struct of_geom *geom, int js,
           double *vmax, double *vmin)
{
  double discr,vp,vm,bsq,EE,EF,va2,cs2,cms2,rho,u ;
  double Acov[NDIM],Bcov[NDIM],Acon[NDIM],Bcon[NDIM] ;
  double Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,C ;
  int j ;
  
  //Setting  up basis vectors for the wave propagation:
  DLOOPA Acov[j] = 0. ;
  Acov[js] = 1. ;
  raise(Acov,geom,Acon) ;
  
  DLOOPA Bcov[j] = 0. ;
  Bcov[TT] = 1. ;
  raise(Bcov,geom,Bcon) ;
  
  /* find fast magnetosonic speed */
  bsq = dot(q->bcon,q->bcov) ;
  rho = pr[RHO] ;
  u = pr[UU] ; // internal energy
  EF = rho + gam*u ; // fluid energy density
  EE = bsq + EF ;  // total energy density

  va2 = bsq/EE ; // Alfvén speed squared 
  cs2 = gam*(gam - 1.)*u/EF ; // Sound speed squared 
  
  //	if(cs2 < 0.) cs2 = SMALL ;
  //	if(cs2 > 1.) cs2 = 1. ;
  //	if(va2 < 0.) va2 = SMALL ;
  //	if(va2 > 1.) va2 = 1. ;
  
  cms2 = cs2 + va2 - cs2*va2 ;	/* and there it is... */
  
  //cms2 *= 1.1 ;
  
  /* check on it! */
  if(cms2 < 0.) {
    fail(FAIL_COEFF_NEG) ;
    cms2 = SMALL ;
  }
  if(cms2 > 1.) {
    fail(FAIL_COEFF_SUP) ;
    cms2 = 1. ;
  }
  
  /* now require that speed of wave measured by observer
   q->ucon is cms2  
   Computing various dot products needed for the characteristic equation, which is
   [-cms²(g^μν + u^μ u^ν) + b^μ b^ν/ρh]k_μ k_ν = 0 
   where:
k_μ is the wave vector
g^μν is the metric
u^μ is the four-velocity
b^μ is the magnetic four-vector
ρh is the enthalpy density
*/
  Asq = dot(Acon,Acov) ;
  Bsq = dot(Bcon,Bcov) ;
  Au =  dot(Acov,q->ucon) ;
  Bu =  dot(Bcov,q->ucon) ;
  AB =  dot(Acon,Bcov) ;
  Au2 = Au*Au ;
  Bu2 = Bu*Bu ;
  AuBu = Au*Bu ;
  
  /* Coefficients of the quadratic equation for wave speeds */
  A =      Bu2  - (Bsq + Bu2)*cms2 ;
  B = 2.*( AuBu - (AB + AuBu)*cms2 ) ;
  C =      Au2  - (Asq + Au2)*cms2 ;
  
  /* Solving the quadratic equation for wave speeds */
  discr = B*B - 4.*A*C ;
  if((discr<0.0)&&(discr>-1.e-10)) discr=0.0;
  else if(discr < -1.e-10) {
    fprintf(stderr,"\n\t %g %g %g %g %g\n",A,B,C,discr,cms2) ;
    fprintf(stderr,"\n\t q->ucon: %g %g %g %g\n",q->ucon[0],q->ucon[1],
            q->ucon[2],q->ucon[3]) ;
    fprintf(stderr,"\n\t q->bcon: %g %g %g %g\n",q->bcon[0],q->bcon[1],
            q->bcon[2],q->bcon[3]) ;
    fprintf(stderr,"\n\t Acon: %g %g %g %g\n",Acon[0],Acon[1],
            Acon[2],Acon[3]) ;
    fprintf(stderr,"\n\t Bcon: %g %g %g %g\n",Bcon[0],Bcon[1],
            Bcon[2],Bcon[3]) ;
    fail(FAIL_VCHAR_DISCR) ;
    discr = 0. ;
  }
  
  discr = sqrt(discr) ;
  vp = -(-B + discr)/(2.*A) ; // Fast wave speed
  vm = -(-B - discr)/(2.*A) ; // Slow wave speed
  
  /*
  The physical meaning:
    vmax and vmin represent the maximum and minimum characteristic speeds
    These speeds are crucial for:
    Determining the CFL (courant) condition for numerical stability
    Computing the numerical fluxes in Godunov-type schemes
    Resolving wave propagation in the MHD system
  */
  if(vp > vm) {
    *vmax = vp ;
    *vmin = vm ;
  }
  else {
    *vmax = vm ;
    *vmin = vp ;
  }
  
  return ;
}


/* Add any additional source terms (e.g. cooling functions) */

void misc_source(double *ph, double *phxp1, double *phxm1, double *phyp1, double *phym1,int ii, int jj, int kk, struct of_geom *geom,
                 struct of_state *q, double *dU, double Dt)
{
  
  double dudx1,dudx2;
  
  
  /* This is merely an example and does not represent any physical source term that I can think of */
  /* Place your calculation for the extra source terms here */
  //  dU[RHO]  += ph[RHO] ;
  //  dU[UU ]  += ph[UU ] ;
  //  dU[U1 ]  += ph[U1 ] ;
  //  dU[U2 ]  += ph[U2 ] ;
  //  dU[U3 ]  += ph[U3 ] ;
  
  
  
  
  
  
}
