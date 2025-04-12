//
//  nuclear.c
//  HARM2D
//
//  Created by Alexander Tchekhovskoy on 9/8/15.
//  Copyright (c) 2015 Home. All rights reserved.
//

#include "decs.h"
#include "nuclear.h"

#if(DONUCLEAR)

void nuc_evol(double pi[][N2M][N3M][NPR],double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], double Dt, int i, int j, int k, int was_floor_activated)
{
  //all sorts of nuclear physics
  double Xalpha, Xn, Xp, Xnp, rho, rhofloor, rhotot;
  double Xalphanew, dXalpha;
  double fac;
  double T, Ye;
  const double rho_unit = compute_rhounit();
  const double T_unit = compute_Tunit();
  const double dq_unit = compute_dq_unit();
  struct of_geom geom;
  FTYPE ucon[NDIM];
  
  /////////
  //
  // Compute mass fractions: Xnp, Xalpha
  //
  
  //re-normalize densities
  if( pr[i][j][k][RHONP] < 0. ) pr[i][j][k][RHONP] = 0.;
  if( pr[i][j][k][RHOALPHA] < 0. ) pr[i][j][k][RHOALPHA] = 0.;
  if( pr[i][j][k][RHOFLOOR] < 0. ) pr[i][j][k][RHOFLOOR] = 0.;
  
  //total density
  rho = pr[i][j][k][RHO];
  //floor+ambient density
  rhofloor = pr[i][j][k][RHOFLOOR];

  //total physical density
  if(rho < 0. ) rho = 0.;
  if(rhofloor < 0. ) rho = 0.;
  
  fac = rho/(rhofloor+pr[i][j][k][RHONP]+pr[i][j][k][RHOALPHA]+SMALL);
  
  //rescale rho_alkpha and rho_np to give the total density
  pr[i][j][k][RHOALPHA] *= fac;
  pr[i][j][k][RHONP]    *= fac;
  pr[i][j][k][RHOFLOOR]    *= fac;
  
  Xalpha = pr[i][j][k][RHOALPHA] / (rho+SMALL);
  Xnp = pr[i][j][k][RHONP] / (rho+SMALL);
  Xfloor = pr[i][j][k][RHOFLOOR] / (rho+SMALL);
  
  //
  // End compute mass fractions
  //
  /////////
  
  //advected Ye
  Ye = pr[i][j][k][YE];

  //compute temperature accounting for a mixture of gas and neutrinos
  T = compute_temperature(rhotot, (gam-1)*ug, Ye);
  
  if( T < SMALL ) T = SMALL; //to avoid division by zero later on

  //compute degeneracy
  etae = compute_degeneracy(rho, T, Ye);
  
  //Xalpha + dXalpha = min(2Ye,2-2Ye) * (1-min(1,Xwb))
  k_cgs = 1.380658e-16; //erg/K
  eV_cgs = 1.6021772e-12; //erg/eV
  MeV_cgs = 1.e-6 * k_cgs / eV_cgs;
  T_MeV = T*T_unit*MeV_cgs; //T in units of MeV
  rho_10 = rho * r_unit * 1.e-10; //rho in units of 1e10 g/cm^3
  Xwb = 15.58*pow(T_MeV,1.125)*pow(rho_10,-0.75)*exp(-7.074/T_MeV);

  Xalphanew = MY_MIN(2.*Ye,2.*(1.-Ye)) * (1.-MY_MIN(1.,Xwb)) - Xfloor;
  Xalphanew = MY_MAX(Xalphanew, 1.e-10);

  Xp = Ye - 0.5*Xalpha - Xfloor;
  Xp = MY_MAX(Xp, 1.e-10);
  
  Xn = 1. - Xp - Xalpha - Xfloor;
  Xn = MY_MAX(Xn, 1.e-10);
  
  
  dXalpha = Xalphanew - Xalpha;
  
  //no need for computing u^t because dq is already in the fluid frame
  //and don't care about the time it took to take a step
  //get_geometry(i,j,k,CENT,&geom) ;
  //ucon_calc(pr[i][j][k], &geom, ucon) ;

  //update Xalpha, Xn, and Xp
  Xalpha = Xalphanew;

  //update the mass fractions
  pr[i][j][k][RHONP] = (Xn+Xp)*rho;
  pr[i][j][k][RHOALPHA] = Xalpha*rho;
  
  //heat per unit mass converted to code units from cgs
  dqalpha = 6.8e18*dXalpha/dq_unit;  //heating in a time Dt
  
  //apply nuclear heating directly to internal energy
  pr[i][j][k][UU] += dqalpha * rho;
  
  pr[i][j][k][YE] += dG;
  
}

double compute_rhounit()
{
  const double Mbh_cgs = 3.*1.99e33;
  const double mneutron_cgs = 1.6749286e-24;
  const double mneutron = mneutron_cgs/Mbh_cgs;
  const double G = 6.67259e-8;
  const double c = 2.99792458e10;
  const double hbar = 1.05457266e-27;
  const double rg_cgs = G*Mbh_cgs/(c*c);
  const double rho_unit = Mbh_cgs / (rg_cgs*rg_cgs*rg_cgs);
  return( rho_unit );
}

double compute_dq_unit()
{
  const double Mbh_cgs = 3.*1.99e33;
  const double mneutron_cgs = 1.6749286e-24;
  const double mneutron = mneutron_cgs/Mbh_cgs;
  const double G = 6.67259e-8;
  const double c = 2.99792458e10;
  const double hbar = 1.05457266e-27;
  const double rg_cgs = G*Mbh_cgs/(c*c);
  const double rho_unit = Mbh_cgs / (rg_cgs*rg_cgs*rg_cgs);
  const double mass_unit = Mbh_cgs;  //[g], BH mass
  //const double time_unit = rg_cgs/c; //[s], light crossing time of r_g
  const double energy_unit = mass_unit*c*c; //[erg], energy unit
  const double dq_unit = energy_unit / (mass_unit);
  return( dq_unit );
}


double compute_Tunit()
{
  const double mneutron_cgs = 1.6749286e-24;
  const double T_unit = mneutron_cgs*c*c/k_cgs;
  return(T_unit);
}

//equation of state stuff: rho, p -> T
double compute_temperature(double rho, double p, double Ye)
{
  //COMPUTE THE TEMPERATURE (no nuclear recombination yet)
  
  const double Mbh_cgs = 3 * 1.99e33;
  const double mneutron_cgs = 1.6749286e-24;
  const double c = 2.99792458e10;
  const double G = 6.67259e-8;
  const double rg_cgs = G*Mbh_cgs/(c*c);
  const double hbar = 1.05457266e-27;
  const double lambdac_cgs = hbar / (mneutron_cgs*c);
  const double mneutron = mneutron_cgs/Mbh_cgs;
  const double lambdac = lambdac_cgs/rg_cgs;
  double atemp = (M_PI*M_PI/45.) * (1./rho) * mneutron * pow(lambdac,-3.);
  double btemp = (1. + Ye)*(rho/p);
  double Fzero, dFzero;
  
  //initial guess assuming radiation-pressure-dominance
  double x = pow(atemp,-0.25);
  
  double prec_newton = 1e-6;
  int ind;
  const int max_iter = 50;
  
  for( ind=1; ind < max_iter; ind++ ) {
    Fzero  =    atemp*pow(x,4.) + btemp*x - 1.;
    dFzero = 4.*atemp*pow(x,3.) + btemp;
    
    x = x - Fzero/dFzero;
    if (fabs(Fzero/dFzero) < prec_newton*x) break;
  }
  if (max_iter == ind) {
    fprintf(stderr, "Cool: N-R for temperature did not converge\n" );
    exit(2345);
  }
  return(x);
}

//accepts temperature in code units T = k Tcgs / mneutron_cgs c^2
double compute_degeneracy(double rho, double T, double Ye)
{
  //Compute degeneracy parameter
  const double Mbh_cgs = 3.*1.99e33;
  const double mneutron_cgs = 1.6749286e-24;
  const double mneutron = mneutron_cgs/Mbh_cgs;
  const double G = 6.67259e-8;
  const double c = 2.99792458e10;
  const double hbar = 1.05457266e-27;
  const double rg_cgs = G*Mbh_cgs/(c*c);
  const double rho_unit = Mbh_cgs / (rg_cgs*rg_cgs*rg_cgs);
  const double k_cgs = 1.380658e-16;
  double ANN, BOB;
  const double CHARLIE = (M_PI*M_PI*M_PI)*(M_PI*M_PI*M_PI)/27.;
  double pf_cgs = hbar * pow(3.*M_PI*M_PI*Ye*rho*rho_unit/mneutron_cgs,1./3.);
  double T_cgs = T*mneutron_cgs*c*c/k_cgs;
  double pfcokT = pf_cgs * c / (k_cgs * T_cgs);
  double pfcokT3 = (pfcokT*pfcokT)*pfcokT;
  double pfcokT6 = pfcokT3*pfcokT3;
  double etae;
  
  if (T_cgs > 1.1604519308e+10) {
    ANN = 0.5*pfcokT3;
    BOB = sqrt( 0.25*pfcokT6 + CHARLIE );
    etae = pow(ANN+BOB,1./3.) - pow(BOB-ANN,1./3.);
  }
  else {
    etae = 0.;
  }
  return(etae);
}

#endif