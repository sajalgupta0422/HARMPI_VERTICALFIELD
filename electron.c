//
//  electron.c
//  HARM2D
//
//  Created by Sean Ressler on 1/13/15.
//  Copyright (c) 2015 Home. All rights reserved.
//

#include "decs.h"
#include "electron.h"

/***********************************************************************************************/
/***********************************************************************************************
 eVol():
 ---------
 -- Evolves the electron entropy;
 
 ***********************************************************************************************/
void eVol(double pi[][N2M][N3M][NPR], double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], double Dt, int i , int j, int k, int was_floor_activated,int is_after_fixup)
{
    

    
#if(eHEAT)
        Heating(pi,prh,pr,i,j,k, Dt, was_floor_activated,is_after_fixup);
#endif
    
    
    if (!is_after_fixup){
#if(eCOND)
        Conduction(pi,prh,pr,i,j,k, Dt, was_floor_activated);
#endif
    }
    
  
}
/***********************************************************************************************/
/***********************************************************************************************
 felcalc():
 ---------
 -- Calculates the fraction of heat given to electrons as a function of beta;
 ***********************************************************************************************/
        
double felcalc(double pi[][N2M][N3M][NPR],int i, int j, int k, double Tel){
    struct of_geom geom ;
    double beta, fel, c1, Trat, Tpr,mbeta, qrat;
    double ppres, bsq;
    double c2,c3,c22,c32;

        
    //Calculations for fel
    c1 = .92 ;//heating constant
        
    Tpr = (gam-1.)*pi[i][j][k][UU]/pi[i][j][k][RHO];
        
    if(Tel<=0.) Tel = SMALL;
    Trat = fabs(Tpr/Tel) ;
    
    ppres = pi[i][j][k][RHO]*Tpr ;  //Proton pressure
    
    
    get_geometry(i,j,k,CENT,&geom) ;
    bsq = bsq_calc(pi[i][j][k],&geom) ; //Magnetic pressure times 2


    /* magnetization parameter */

    double sig = bsq/pi[i][j][k][RHO];
    
    beta = ppres/bsq * 2. ;
    if(beta>1.e20)beta = 1.e20;
    mbeta = 2.-.2*log10(Trat);
    
    if(Trat<=1.){
        
        c2 = 1.6 / Trat ;
        c3 = 18. + 5.*log10(Trat);
        
    }
    else{
        c2 = 1.2/ Trat ;
        c3 = 18. ;
    }
    
    c22 = pow(c2,2.);
    c32 = pow(c3,2.);
    
    qrat = c1 * (c22+pow(beta,mbeta))/(c32 + pow(beta,mbeta)) * exp(-1./beta)*pow(mrat*Trat,.5) ;
        
    fel = 1./(1.+qrat);

    // if (sig>1.){
    //     fel = 0.;
    // }
    
    
    return fel;
}
/***********************************************************************************************/
/***********************************************************************************************
 felcalc2():
 ---------
 -- Calculates the fraction of heat given to electrons as a function of beta;
 -- This version is more accurate than felcalc(): it accounts for diff btw Tgas and Tp
 ***********************************************************************************************/

double felcalc2(double pi[][N2M][N3M][NPR],int i, int j, int k, double Tel){
    struct of_geom geom ;
    double beta, fel, c1, Trat, Tg, Tpr, mbeta, qrat;
    double ppres, bsq;
    double c2,c3,c22,c32;


    //Calculations for fel
    c1 = .92 ;//heating constant

    //total gas and electron temperatures
    Tg = (gam-1.)*pi[i][j][k][UU]/pi[i][j][k][RHO];

    //limit electron temperature to be no more than gas temperature
    if(Tel>=Tg) Tel = Tg;
    if(Tel<=SMALL) Tel = SMALL;

    //proton temperature
    Tpr = Tg - Tel;
    if(Tpr<=SMALL) Tpr = SMALL;

    Trat = fabs(Tpr/Tel) ;

    ppres = pi[i][j][k][RHO]*Tpr ;  //Proton pressure


    get_geometry(i,j,k,CENT,&geom) ;
    bsq = bsq_calc(pi[i][j][k],&geom) ; //Magnetic pressure times 2


    /* magnetization parameter */

    double sig = bsq/pi[i][j][k][RHO];

    beta = 2.*ppres/(bsq+SMALL);
    if(beta>1.e20)beta = 1.e20;
    mbeta = 2.-.2*log10(Trat);

    if(Trat<=1.){

        c2 = 1.6 / Trat ;
        c3 = 18. + 5.*log10(Trat);

    }
    else{
        c2 = 1.2/ Trat ;
        c3 = 18. ;
    }

    c22 = pow(c2,2.);
    c32 = pow(c3,2.);

    qrat = c1 * (c22+pow(beta,mbeta))/(c32 + pow(beta,mbeta)) * exp(-1./beta)*pow(mrat*Trat,.5) ;




    fel = 1./(1.+qrat);

    // if (sig>1.){
    //     fel = 0.;
    // }


    return fel;
}

/* JAD adding Kawazura fitting formula equation 2 using felcalc2 as template
   which is supposed to be more accurate than felcalc */

double felkawa(double pi[][N2M][N3M][NPR],int i, int j, int k, double Tel){
    struct of_geom geom ;
    double beta, fel, c1, Trat, Tg, Tpr, mbeta, qrat;
    double ppres, bsq;
    double c2,c3,c22,c32;

    //total gas and electron temperatures
    Tg = (gam-1.)*pi[i][j][k][UU]/pi[i][j][k][RHO];

    //limit electron temperature to be no more than gas temperature
    if(Tel>=Tg) Tel = Tg;
    if(Tel<=SMALL) Tel = SMALL;

    //proton temperature
    Tpr = Tg - Tel;
    if(Tpr<=SMALL) Tpr = SMALL;

    Trat = fabs(Tpr/Tel) ;

    ppres = pi[i][j][k][RHO]*Tpr ;  //Proton pressure

    get_geometry(i,j,k,CENT,&geom) ;
    bsq = bsq_calc(pi[i][j][k],&geom) ; //Magnetic pressure times 2

    /* magnetization parameter */

    //    double sig = bsq/pi[i][j][k][RHO];

    beta = 2.*ppres/(bsq+SMALL);
    if(beta>1.e20)beta = 1.e20;
    //    mbeta = 2.-.2*log10(Trat);

    /*    if(Trat<=1.){

        c2 = 1.6 / Trat ;
        c3 = 18. + 5.*log10(Trat);

    }
    else{
        c2 = 1.2/ Trat ;
        c3 = 18. ;
	}

    c22 = pow(c2,2.);
    c32 = pow(c3,2.);*/

    qrat = 35./(1.+pow(beta/15.,-1.4)*exp(-0.1/Trat));

    fel = 1./(1.+qrat);

    return fel;
}

/*felwerner:
 ---------
 -- approximate heating fraction from reconnection following Werner+2018 fitting formula
 -- simply varies between fel = 1/4 and 1/2 as a function of sigma_ion
 ***********************************************************************************************/

double felwerner(double pi[][N2M][N3M][NPR],int i, int j, int k){
    struct of_geom geom ;
    double sig, fel;
    double bsq;

    //Calculations
    get_geometry(i,j,k,CENT,&geom) ;
    bsq = bsq_calc(pi[i][j][k],&geom) ; //Magnetic pressure times 2
    /* magnetization parameter */

    sig = bsq/pi[i][j][k][RHO];
    fel = 1./4.+1./4.*sqrt(sig/5./(2.+sig/5.));

    return fel;
}
/****************************************************************
felchael():
 ---------
 -- Calculates the fraction of heat given to electrons as a function of beta, sigma following Chael+ based on Rowan+ reconnection simulations
 ***********************************************************************************************/

double felchael(double pi[][N2M][N3M][NPR],int i, int j, int k){
    struct of_geom geom ;
    double beta, fel, betamax;
    double bsq, sigmaw;

    get_geometry(i,j,k,CENT,&geom) ;
    bsq = bsq_calc(pi[i][j][k],&geom) ; //Magnetic pressure times 2

    /* magnetization parameter */

    //double sig = bsq/pi[i][j][k][RHO];

    //beta = 2.*ppres/(bsq+SMALL);
    //if(beta>1.e20)beta = 1.e20;

    // Chael+2018 equation 9, this version assumes gam=game=gamp OR 
    // up >> ue and gam=gamp
    sigmaw = bsq/(pi[i][j][k][RHO]+gam*pi[i][j][k][UU]);
   
    betamax=1./4./sigmaw;
    // Chael+2018 equation 13
    if (beta < betamax) {
      fel=1./2.*exp(-(1.-beta/betamax)/(0.8+sigmaw));
    }
    else {
      fel=1./2.;
    }

    return fel;
}


/***********************************************************************************************/
/***********************************************************************************************
 Heating():
 ---------
 -- performs the update to the elctron entropy given by the heating terms and resets the total entropy to its actual value;
 
 ***********************************************************************************************/
void Heating(double pi[][N2M][N3M][NPR],double prh[][N2M][N3M][NPR], double pr[][N2M][N3M][NPR], int i, int j, int k, double Dt, int was_floor_activated, int is_after_fixup)
{
#if(eHEAT)
    struct of_geom geom ;
    struct of_state q;
    double ucon[NDIM];
    double kad;
    double kappatot;
    double ughat, ugharm;
    double ueldis, uel4, uel4a, uel4b, uel4c, uel4d, uel4e, uel5;
    double Teldis,ueldisi;
    double Ttot, Tel4, Tel4a, Tel4b, Tel4c, Tel4d, Tel4e, Tel5, Telvar;
    double uel4i, uel4ai, uel4bi, uel4ci, uel4di, uel4ei, uel5i,du;
    double fel5,fel4,fel4a,fel4b,fel4c,fel4d,fel4e,feldis;
    double fac, fac9;
    double rho, irho, rhofac;
    double dq;
    double X[NDIM], r, th, phi;
    double felcalc2(double pi[][N2M][N3M][NPR],int i, int j, int k, double Tel);
    double felwerner(double pi[][N2M][N3M][NPR],int i, int j, int k);
    double felchael(double pi[][N2M][N3M][NPR],int i, int j, int k);
        
    fel4 = 0.5;
    fel5 = 0.5 ;
    feldis =  0.5;
    
    kappatot =(gam-1.)*pr[i][j][k][UU]*pow(pr[i][j][k][RHO],-gam);
    
    ueldis = pr[i][j][k][KELDIS]*pow(pr[i][j][k][RHO],game)/(game-1.);
    rhofac = pow(pr[i][j][k][RHO],game4)/(game4-1.);
    uel4 = pr[i][j][k][KEL4]*rhofac;
    uel4a = pr[i][j][k][KEL4A]*rhofac;
    uel4b = pr[i][j][k][KEL4B]*rhofac;
    uel4c = pr[i][j][k][KEL4C]*rhofac;
    uel4d = pr[i][j][k][KEL4D]*rhofac;
    uel4e = pr[i][j][k][KEL4E]*rhofac;
    uel5 = pr[i][j][k][KEL5]*pow(pr[i][j][k][RHO],game5)/(game5-1.);
        
    //Temperatures
    Ttot =(gam-1.)*pr[i][j][k][UU]/pr[i][j][k][RHO];
    rhofac = pow(prh[i][j][k][RHO],game4);
    uel4i = 1./(game4-1.)*prh[i][j][k][KEL4]*rhofac;
    uel4ai = 1./(game4-1.)*prh[i][j][k][KEL4A]*rhofac;
    uel4bi = 1./(game4-1.)*prh[i][j][k][KEL4B]*rhofac;
    uel4ci = 1./(game4-1.)*prh[i][j][k][KEL4C]*rhofac;
    uel4di = 1./(game4-1.)*prh[i][j][k][KEL4D]*rhofac;
    uel4ei = 1./(game4-1.)*prh[i][j][k][KEL4E]*rhofac;
    uel5i = 1./(game5-1.)*prh[i][j][k][KEL5]*pow(prh[i][j][k][RHO],game5);
    ueldisi =1./(game-1.)*prh[i][j][k][KELDIS]*pow(prh[i][j][k][RHO],game);
    irho = 1./prh[i][j][k][RHO];
    Tel4 = (game4-1.)*uel4i*irho;
    Tel4a = (game4-1.)*uel4ai*irho;
    Tel4b = (game4-1.)*uel4bi*irho;
    Tel4c = (game4-1.)*uel4ci*irho;
    Tel4d = (game4-1.)*uel4di*irho;
    Tel4e = (game4-1.)*uel4ei*irho;
    Tel5 = (game5-1.)*uel5i*irho;;
    Teldis =(game-1.)*ueldisi*irho;;
        
    //Calculations for fel
    
#if (BETAHEAT)
    fel4  = felcalc(prh,i,j,k,Tel4);
    fel4a = felcalc(prh,i,j,k,Tel4a); fel4a = MY_MAX(fel4a,0.1);
    //    fel4b = felcalc2(prh,i,j,k,Tel4b);
    //    fel4c = felcalc(prh,i,j,k,Tel4c);
    //    fel4d = felcalc(prh,i,j,k,Tel4d); fel4d = MY_MAX(fel4d,0.1);
    // REPLACING THESE WITH NEW VERSIONS OF CONSTANT/WERNER/CHAEL heating
    // REPLACING AGAIN CONSTANT FEL4B WITH KAWAZURA+2018 AND FLOOR OF 1/10
    fel4b = felkawa(prh,i,j,k,Tel4b); fel4b = MY_MAX(fel4b,0.1);
    fel4c = felwerner(prh,i,j,k);
    fel4d = felchael(prh,i,j,k);
    fel4e = felcalc(prh,i,j,k,Tel4e);
    fel5 = felcalc(prh,i,j,k,Tel5);
    feldis = felcalc(prh,i,j,k,Teldis);
    
#endif
    
    
    
    
    kad = pr[i][j][k][KTOT];
    
    ughat = kad*pow(pr[i][j][k][RHO],gam)/(gam-1.);
    ugharm =pr[i][j][k][UU];
    
    
    
    
#if(0&&DOFLR) //disabled the floor suppression
    //smoothly reset Theta_e to 0.1 if the floor was activated
    fac = (0.6-pr[i][j][k][FLR])*10.;
    fac = MY_MAX(0.,fac);
    fac = MY_MIN(1.,fac);
    fac9 = (0.3-pr[i][j][k][FLR])*10.;
    fac9 = MY_MAX(0.,fac9);
    fac9 = MY_MIN(1.,fac9);

#else
    fac = 1.;
    fac9 = 1.;
#endif
    
    //update electron internal energy
    
    du = (ugharm-ughat)*pow(prh[i][j][k][RHO]/pr[i][j][k][RHO],gam-game4);
    ueldis+= fac*feldis*du;
    uel4  += fac*fel4*du;
    uel4a += fac*fel4a*du;
    uel4b += fel4b*du;
    uel4c += fel4c*du;
    uel4d += fel4d*du;
    uel4e += fel4e*du;
    uel5  += fac*fel5*(ugharm-ughat)*pow(prh[i][j][k][RHO]/pr[i][j][k][RHO],gam-game5);
    
    //convert back to electron entropy
    pr[i][j][k][KELDIS] =ueldis*(game-1.)*pow(pr[i][j][k][RHO],-game);
    rhofac = pow(pr[i][j][k][RHO],-game4);
    pr[i][j][k][KEL4] = (game4-1.)*uel4*rhofac;
    pr[i][j][k][KEL4A] = (game4-1.)*uel4a*rhofac;
    pr[i][j][k][KEL4B] = (game4-1.)*uel4b*rhofac;
    pr[i][j][k][KEL4C] = (game4-1.)*uel4c*rhofac;
    pr[i][j][k][KEL4D] = (game4-1.)*uel4d*rhofac;
    pr[i][j][k][KEL4E] = (game4-1.)*uel4e*rhofac;
    pr[i][j][k][KEL5] =(game5-1.)*uel5*pow(pr[i][j][k][RHO],-game5);    
    
    get_geometry(i,j,k,CENT,&geom);
    get_state(prh[i][j][k],&geom,&q);
    
    dq = q.ucon[0]*q.ucov[0]*(ugharm-ughat)*pow(prh[i][j][k][RHO]/pr[i][j][k][RHO],gam);
    qdot[i][j][k] = dq/Dt;
    //only account for cumulative heating on a full timestep
#if(DOQDOTAVG)
    if(Dt == dt) {
        qdotavg[i][j][k] += dq;
#       if(!BETAHEAT)
            fel4e = felcalc(prh,i,j,k,Tel4e);
#       endif
        feqdotavg[i][j][k] += fel4e * dq;
    }
#endif
    
    //reset to actual internal energy/Entropy
    pr[i][j][k][KTOT] += (kappatot-pr[i][j][k][KTOT]);
    
#endif //eHEAT

}

/***********************************************************************************************/
/***********************************************************************************************
 Conduction():
 ---------
 -- performs the update to the elctron entropy given by electron conduction;
 
 ***********************************************************************************************/
void Conduction( double pi[][N2M][N3M][NPR], double prh[][N2M][N3M][NPR],double pr[][N2M][N3M][NPR], int i, int j, int k, double Dt, int was_floor_activated)
{
    
    int n,idel,kdel,jdel, beta,mu,alpha,nprime;
    struct of_state q, ql, qr,qc,qf,qi;
    struct of_geom geom, geoml,geomr,geomc;
    double bsq, grhofac, Teli;
    double phih, phif, kelf, rhof, rhoh;
    double bhat[NDIM],bhatl[NDIM],bhatr[NDIM], bhatf[NDIM], bhati[NDIM];
    double qcond[NDIM], qcondl[NDIM], qcondr[NDIM],qcondf[NDIM],qcondi[NDIM];
    double uconl[NDIM], ucon[NDIM], uconr[NDIM], uconf[NDIM], uconi[NDIM];
    double dcommau_dx, dcommaq_dx, dTe_dx[NDIM], dTdis_dx[NDIM];
    double dTdis_dxcon[NDIM];
    double bdotdTe, bdotdTenolim;
    double ugradu[NDIM];
    double acov[NDIM];
    double acon[NDIM];
    double Telr, Tell, Telc, Tdisl, Tdisr, Tdis;
    double bsqr,bsql, bsqf, bsqi, ugf;
    double qdota, bdota;
    double a11, a12, a21, a22,b1,b2;
    double kappacond, safety, tau;
    double vtherm, cse,theta, ugel4, kappamax, Txscale;
    double det;
    double adotu, bdotu, ucov[NDIM], qiso[NDIM],qisocov[NDIM];
    double X[NDIM], r,th,phi;
    double keff;
    double dr, dth;
    double phimax, fphi;
    double betae;
    double oldkel4;
    double xx1, xx2, upperLimit1, myrhof, myugf;
    
    
    
     Telc = prh[i][j][k][KEL4]*pow(prh[i][j][k][RHO],game4-1.);
     Tdis = prh[i][j][k][KELDIS]*pow(prh[i][j][k][RHO],game -1.);

    
    
    //initialize variables calculated cumulatively to zero
    
    dcommau_dx = 0;
    dcommaq_dx = 0;
    bdotdTe = 0;
    bdotdTenolim = 0;
    
    for(n=0;n<NDIM;n++){
        ugradu[n] =0;
        acon[n] = 0;
        acov[n] = 0;
    }
    //half time
    get_geometry(i,j,k,CENT,&geom);
    get_state(prh[i][j][k],&geom,&q);
    bsq = dot(q.bcon,q.bcov);

    
    for(n=0;n<NDIM;n++){
        bhat[n] =q.bcon[n]/(sqrt(bsq)+SMALL);
        qcond[n] = bhat[n]*prh[i][j][k][PHI];
        ucon[n] = q.ucon[n];
        
    }
    
    //initial time
    get_state(pi[i][j][k],&geom,&qi);
    bsqi = dot(qi.bcon,qi.bcov);
    
    for(n=0;n<NDIM;n++){
        bhati[n] =qi.bcon[n]/(sqrt(bsqi)+SMALL);
        qcondi[n] = bhati[n]*pi[i][j][k][PHI];
        uconi[n] = qi.ucon[n];
        
    }
    
    //final time
    get_state(pr[i][j][k],&geom,&qf);
    bsqf = dot(qf.bcon,qf.bcov);
    ugf = pr[i][j][k][UU];
  
    
    for(n=0;n<NDIM;n++){
        bhatf[n] =qf.bcon[n]/(sqrt(bsqf)+SMALL);
        qcondf[n] = bhatf[n]*pr[i][j][k][PHI];
        uconf[n] = qf.ucon[n];
        
    }
    
    
    
    //Now spatial derivatives
    
    for(n=1;n<NDIM;n++){
        
        idel = 0;
        jdel = 0;
        kdel = 0;
        if(n==1) idel = 1;
        if(n==2) jdel = 1;
        if(n==3) kdel = 1;
        
        if(N2==1) jdel = 0;
        if(N3==1) kdel = 0;
        
        get_geometry(i+idel,j+jdel,k+kdel,CENT,&geomr);
        get_state(prh[i+idel][j+jdel][k+kdel],&geomr,&qr);
        
        get_geometry(i-idel,j-jdel,k-kdel,CENT,&geoml);
        get_state(prh[i-idel][j-jdel][k-kdel],&geoml,&ql);
        
        bsqr = dot(qr.bcon,qr.bcov);
        bhatr[n] =qr.bcon[n]/(sqrt(bsqr)+SMALL);
        qcondr[n] = bhatr[n]*prh[i+idel][j+jdel][k+kdel][PHI];
        uconr[n] = qr.ucon[n];
        Telr = prh[i+idel][j+jdel][k+kdel][KEL4]*pow(prh[i+idel][j+jdel][k+kdel][RHO],game4-1.);
        Tdisr =prh[i+idel][j+jdel][k+kdel][KELDIS]*pow(prh[i+idel][j+jdel][k+kdel][RHO],game-1.);
        
        bsql = dot(ql.bcon,ql.bcov);
        bhatl[n] =ql.bcon[n]/(sqrt(bsql)+SMALL);
        qcondl[n] =bhatl[n]*prh[i-idel][j-jdel][k-kdel][PHI];
        uconl[n] = ql.ucon[n];
        Tell = prh[i-idel][j-jdel][k-kdel][KEL4]*pow(prh[i-idel][j-jdel][k-kdel][RHO],game4-1.);
        Tdisl = prh[i-idel][j-jdel][k-kdel][KELDIS]*pow(prh[i-idel][j-jdel][k-kdel][RHO],game-1.);
        
        
        dcommau_dx +=  (slope_lim(geoml.g*uconl[n], geom.g*ucon[n],geomr.g*uconr[n]))/dx[n];
        dcommaq_dx +=  (slope_lim(geoml.g*qcondl[n], geom.g*qcond[n],geomr.g*qcondr[n]))/dx[n];
        
        dTe_dx[n] = slope_lim(Tell,Telc,Telr)/dx[n];
        dTdis_dx[n] = slope_lim(Tdisl,Tdis,Tdisr)/dx[n];
        bdotdTe += dTe_dx[n]*bhat[n];
        
        for ( nprime=1;nprime<NDIM;nprime++){  //ugradu[i]
            
            idel = 0;
            jdel = 0;
            kdel = 0;
            if(nprime==1) idel = 1;
            if(nprime==2) jdel = 1;
            if(nprime==3) kdel = 1;
            
            if(N2==1) jdel = 0;
            if(N3==1) kdel = 0;
            
            
            get_geometry(i+idel,j+jdel,k+kdel,CENT,&geomr);
            get_state(prh[i+idel][j+jdel][k+kdel],&geomr,&qr);
            
            uconr[n] = qr.ucon[n];
            
            
            get_geometry(i-idel,j-jdel,k-kdel,CENT,&geoml);
            get_state(prh[i-idel][j-jdel][k-kdel],&geoml,&ql);
            
            uconl[n] = ql.ucon[n];
            
            ugradu[n] += ucon[nprime]*slope_lim(uconl[n],ucon[n],uconr[n])/dx[nprime];
        }
    
        ugradu[n] += ucon[0]*(uconf[n]-uconi[n])/Dt;
    }
    
    
    //add temporal components if necessary
    
    dcommau_dx += geom.g*(uconf[0]-uconi[0])/Dt;
    
    for ( nprime=1;nprime<NDIM;nprime++){   //ugradu[0]
        
        idel = 0;
        jdel = 0;
        kdel = 0;
        if(nprime==1) idel = 1;
        if(nprime==2) jdel = 1;
        if(nprime==3) kdel = 1;
        
        if(N2==1) jdel = 0;
        if(N3==1) kdel = 0;
        
        
        get_geometry(i+idel,j+jdel,k+kdel,CENT,&geomr);
        get_state(prh[i+idel][j+jdel][k+kdel],&geomr,&qr);
        
        uconr[0] = qr.ucon[0];
        
        
        get_geometry(i-idel,j-jdel,k-kdel,CENT,&geoml);
        get_state(prh[i-idel][j-jdel][k-kdel],&geoml,&ql);
        
        uconl[0] = ql.ucon[0];
        
        ugradu[0] += ucon[nprime]* slope_lim(uconl[0],ucon[0],uconr[0])/dx[nprime];
        
        
        
    }
    
    ugradu[0] += ucon[0]*(uconf[0]-uconi[0])/Dt;  //u[0]d/dt u[0]
    
    
    //add connection coefficients, if necessary
    
    for ( mu = 0;mu<NDIM;mu++){
        for ( alpha=0; alpha<NDIM; alpha++){
            
            for ( beta=0; beta<NDIM; beta++){
                acon[mu] += conn[i][j][k][mu][alpha][beta]*ucon[alpha]*ucon[beta];
            }
        }
        acon[mu]+=ugradu[mu];
        
    }
    
    
    lower(acon,&geom,acov);
    
    
    qdota = dot(qcond,acov);
    bdota = dot(bhat,acov);
    adotu = dot(acov,ucon); //should be zero!
    lower(ucon,&geom,ucov);
    bdotu = dot(bhat,ucov);
    
    double Telf, Tdisf, Tdisi;
    grhofac =(game4-1.)*pow(prh[i][j][k][RHO],1.-game4);
    Teli =pi[i][j][k][KEL4]*pow(pi[i][j][k][RHO],game4-1.);
    Telf =pr[i][j][k][KEL4]*pow(pr[i][j][k][RHO],game4-1.);
    Tdisf =pr[i][j][k][KELDIS]*pow(pr[i][j][k][RHO],game-1.);
    Tdisi =pi[i][j][k][KELDIS]*pow(pi[i][j][k][RHO],game-1.);
    phih = prh[i][j][k][PHI];
    rhof = pr[i][j][k][RHO];
    rhoh = prh[i][j][k][RHO];
    kelf = pr[i][j][k][KEL4];
    phif = pr[i][j][k][PHI];
    ugel4 = 1./(game4-1.)*prh[i][j][k][KEL4]*pow(rhoh,game4);

    betae = 2./(bsq+SMALL) * ugel4*(game4-1.);
    
    //kappacond = .5*rhoh;
    kappacond = .5;
#if(BL)
    coord(i,j,k,CENT,X) ;
    get_phys_coord(i,j,k,&r,&th,&phi);
    
    kappacond = 10.*r*rhoh;
    kappamax =rhoh*sqrt(r)*r * pow(4.*betae+SMALL,-.8);
    
    if (kappacond>kappamax) kappacond = kappamax;
    
#endif
    
    
    safety = 2.;
    
    theta = Telc*mrat;
    vtherm = sqrt(theta)*sqrt(theta+2)/(theta+1);
    cse = sqrt(game4*(game4-1)*ugel4/(rhoh/mrat+game4*ugel4));
    Txscale = Telc/(fabs(bdotdTe));
    phimax = (ugel4+rhoh/mrat)*cse;
    
    
    double upperLimit, fermiDirac,lambda;
    lambda = 0.1;
  
    //smoothly switch off conduction in floored regions
    /* myrhof = MY_MAX(rhof, SMALL); */
    /* myugf = MY_MAX(ugf, SMALL); */
    /* xx1 = bsqf/(20.*myrhof); */
    /* xx2 = bsqf/(200.*myugf); */
    /* upperLimit1 = MY_MAX(xx1, xx2); */
    upperLimit      = fabs(pr[i][j][k][PHI])/(phimax);
#if(DOFLR&&0) //do not suppress conduction in floored regions
    upperLimit1 = 2.*pr[i][j][k][FLR];
    upperLimit = MY_MAX( upperLimit, upperLimit1 );
#endif
    /* upperLimit      = MY_MAX(upperLimit,upperLimit1); */
    upperLimit      = (upperLimit-1.)/lambda;
    fermiDirac      = exp(-upperLimit)/(exp(-upperLimit) + 1.) + 1e-10;
    //fphi = (1.- tanh(abs(phih)/(phimax+SMALL)-1.))/2.;
    
    //fermiDirac = 1.;
    fphi = 1.;
    //if(abs(phih)>phimax) fphi =0.;
    
    keff = kappacond*fermiDirac;
    
    
    
    

    
    //dr = dx[1]*r;
    //dth = r*(M_PI+(1.-hslope)/2.*cos(2.*M_PI*X[2]))*dx[2];

    
    //phimax = MY_MIN(rhoh*.01, .1*vtherm*rhoh * Telc);
    // if(kappacond>kappamax) keff = kappamax;
    //keff = (kappacond-kappamax)*tanh(1.*kappamax/kappacond)+kappamax;
    

    
    tau = keff/rhoh/vtherm/vtherm + 0.01;
    
   
    if (tau < safety*dt/(1.-cour)) tau = safety*dt/(1.-cour);
    
    
    if(WHICHPROBLEM == ATM_TEST){
        tau = 10.;
        kappacond = .5;
        keff = kappacond;
        phimax = 1000000000000000000.;
    }
    
    if(WHICHPROBLEM == BONDI_CON){
        coord(i,j,k,CENT,X) ;
        get_phys_coord(i,j,k,&r,&th,&phi);
        kappacond = 10.*rhoh*sqrt(r);
        keff = kappacond;
        tau = 2.*keff*prh[i][j][k][UU]*(gam-1.)/rhoh/(rhoh+gam*prh[i][j][k][UU]) + .25;
        phimax = 1000000000000000000.;
    }


    
    a11 = geom.g*rhof*uconf[0]/Dt;
    a12 = grhofac*geom.g/Dt*bhatf[0];
    a21 = geom.g*rhoh*keff/Dt/tau*bhat[0]*pow(pr[i][j][k][RHO],game4-1.);
    a22 = geom.g*rhof*uconf[0]/Dt;
    
    det = a11*a22-a21*a12;
    
    b1 = geom.g*kelf*rhof*uconf[0]/Dt  + (grhofac)*(geom.g*qcondi[0]/Dt  - dcommaq_dx - qdota*geom.g);
    
    b2 = geom.g*phif*rhof*uconf[0]/Dt  - geom.g*rhoh*phih/tau + geom.g*rhoh*keff/tau * (bhat[0] * Teli/Dt - bdotdTe - bdota*Telc);
    
    pr[i][j][k][PHI] = (-b1*a21 + b2*a11)/det;
    pr[i][j][k][KEL4] = (b1*a22-b2*a12)/det;
  
    
//    if(pr[i][j][k][PHI]>phimax){
//        pr[i][j][k][PHI]= phimax;
//        pr[i][j][k][KEL4]= (b1-a12*pr[i][j][k][PHI])/a11; //Be consistent with new flux
//    }
//    if(pr[i][j][k][PHI]<-phimax){
//        pr[i][j][k][PHI]= -phimax;
//        pr[i][j][k][KEL4]= (b1-a12*pr[i][j][k][PHI])/a11;
//    }
//
//    
    //pr[i][j][k][KEL4]= (b1-a12*pr[i][j][k][PHI])/a11;
    
    //pr[i][j][k][PHI] = (b2-a21*kelf)/a22;
    
    
    
    //pr[i][j][k][PHI] = -keff*(bdota*Telc + bdotdTe + bhat[0]*(Telf-Teli)/Dt);
    
    
   // pr[i][j][k][PHI] = -keff*(bhat[1]*dTe_dx[1]+bhat[1]*Telc*acov[1]);
    
    
    dTdis_dx[0] = (Tdisf-Tdisi)/Dt;

    raise(dTdis_dx,&geom,dTdis_dxcon);
    
    for ( n=0;n<NDIM;n++){
        qiso[n] = dot(dTdis_dx,ucon)*ucon[n]+ dTdis_dxcon[n] + Tdis*acon[n];
    }
    
    lower(qiso,&geom,qisocov);
    
    
    
    qisodotb[i][j][k] = fabs(dot(bhat,qisocov));
    //compute q^2 avoiding nans if dot product is slightly negative
    qisosq[i][j][k] = dot(qiso,qisocov);
    if (qisosq[i][j][k] < 0.) {
      qisosq[i][j][k] = 0.;
    }
    qisosq[i][j][k] = sqrt(qisosq[i][j][k]);
    
}

void fixupuel(double (* pv)[N2M][N3M][NPR])
{
#if( eCOND || eHEAT)
    int i,j, k ;
#pragma omp parallel for schedule(static,N1*N2*N3/nthreads) collapse(3) default(none) shared(pv,nthreads) private(i,j,k)
    ZLOOP {
        fixupuel1zone( i, j, k, pv[i][j][k] );
    }
#endif
}

void fixupuel1zone( int i, int j, int k, double pv[NPR])
{
#if( eCOND || eHEAT)
    void get_rho_u_floor( double r, double th, double phi, double *rho_floor, double *u_floor );
    double uel4, uel4a, uel4b, uel4c, uel4d, uel4e, uel5, uel, uels;
    double fmin;
    double r,th,phi,X[NDIM],uuscal,rhoscal, rhoflr,uuflr ;
    double bsq, thex;
    double thold,aold;
    struct of_geom geom ;
    int dofloor = 0;
  
    fmin = felfloor;
    
    coord(i,j,k,CENT,X) ;
    get_phys_coord(i,j,k,&r,&th,&phi) ;
    
    get_rho_u_floor(r,th,phi,&rhoflr,&uuflr);
  
    //SEAN: I added the below 3 lines, which were apparently missing:
    //compute the square of fluid frame magnetic field (twice magnetic pressure)
    get_geometry(i,j,k,CENT,&geom) ;
    bsq = bsq_calc(pv,&geom) ;
    
    //tie floors to the local values of magnetic field and internal energy density
#if(1)
    if( rhoflr < bsq / BSQORHOMAX ) rhoflr = bsq / BSQORHOMAX;
    if( uuflr < bsq / BSQOUMAX ) uuflr = bsq / BSQOUMAX;
    if( rhoflr < pv[UU] / UORHOMAX ) rhoflr = pv[UU] / UORHOMAX;
#endif
    
    
    if (BL){  //SEAN: this line should be here for consistency with fixup1zone()
        
        
        if( rhoflr < RHOMINLIMIT ) rhoflr = RHOMINLIMIT;
        if( uuflr  < UUMINLIMIT  ) uuflr  = UUMINLIMIT;
        
    } //SEAN: this line should be here for consistency with fixup1zone()
    uel4 = 1./(game4-1.)*pv[KEL4]*pow(pv[RHO],game4);
    uel4a = 1./(game4-1.)*pv[KEL4A]*pow(pv[RHO],game4);
    uel4b = 1./(game4-1.)*pv[KEL4B]*pow(pv[RHO],game4);
    uel4c = 1./(game4-1.)*pv[KEL4C]*pow(pv[RHO],game4);
    uel4d = 1./(game4-1.)*pv[KEL4D]*pow(pv[RHO],game4);
    uel4e = 1./(game4-1.)*pv[KEL4E]*pow(pv[RHO],game4);
    uel5 = 1./(game5-1.)*pv[KEL5]*pow(pv[RHO],game5);
    uel = pv[KELDIS]*1./(game-1.)*pow(pv[RHO],game);
    
    if(uel4  < uuflr*fmin || isnan(uel4))   {
        uel4  = uuflr*fmin;
        pv[KEL4] = (game4-1.)*uel4*pow(pv[RHO],-game4);
        dofloor = 1;
      
    }
    if(uel4a < uuflr*fmin || isnan(uel4a))   {
        uel4a  = uuflr*fmin;
        pv[KEL4A] = (game4-1.)*uel4a*pow(pv[RHO],-game4);
        //dofloor = 1;
        
    }
    if(uel4b < uuflr*fmin || isnan(uel4b))   {
        uel4b  = uuflr*fmin;
        pv[KEL4B] = (game4-1.)*uel4b*pow(pv[RHO],-game4);
        //dofloor = 1;
        
    }
    if(uel4c < uuflr*fmin || isnan(uel4c))   {
        uel4c  = uuflr*fmin;
        pv[KEL4C] = (game4-1.)*uel4c*pow(pv[RHO],-game4);
        //dofloor = 1;
        
    }
    // CHANGING THIS ONE to treat floors the same as uel4abc
    // hope the new e- heating don't have floor problems but who knows
    if(uel4d < uuflr*fmin || isnan(uel4d))   {
        uel4d  = uuflr*fmin;
        pv[KEL4D] = (game4-1.)*uel4d*pow(pv[RHO],-game4);
        //dofloor = 1;
        
    }
    if(isnan(uel4e))   {
        uel4e  = uuflr*fmin;
        pv[KEL4E] = (game4-1.)*uel4*pow(pv[RHO],-game4);
        //dofloor = 1;
        
    }
    if(uel5  < uuflr*fmin || isnan(uel5) )   {
        uel5  = uuflr*fmin;
        pv[KEL5] =(game5-1.)*uel5*pow(pv[RHO],-game5);
        dofloor = 1;
      
    }
    
    if(uel < uuflr*fmin || isnan(uel))   {
        uel  = uuflr*fmin;
        pv[KELDIS] =1./(game-1.)*uel*pow(pv[RHO],-game);
        dofloor = 1;
      
    }
  
#   if(eCOND)
      //if the heat flux is a nan, reset it to zero
      if(isnan(pv[PHI])) {
        pv[PHI] = 0.0;
        dofloor = 1;
      }
#   endif
#if(0 && DOFLR) //"0 &&" disables floor activation if electron int. en. goes below floor
    if(dofloor) {
       pv[FLR] = 1.;
    }
#endif
#endif
}

/* this function sets the entropy
 * based on the initial conditions */
void init_entropy(){
    
    double thex, gex, kappatot, agame, kappa, kappae;
    double kappae4, kappae5, kappaes;
    double uel4, uel5, uels;
    double Tel4, Tel5, Ttot;
    double r,th,phi, X[NDIM];
    int i,j,k ;
    int iglob, jglob, kglob;
    double rancval;
  
    fel0 =0.1;  // Initial ugel/ugtot
    felfloor = .01;
    
    kelmin =1.e40;
    game = 4./3.;
    
    
    ZSLOOP(-N1G,N1+N1G-1,-N2G,N2G+N2-1,-N3G,N3G+N3-1){
        kappatot = (gam-1.)*p[i][j][k][UU]*pow(p[i][j][k][RHO],-gam);
        
        if(fel0*kappatot<kelmin) kelmin = fel0*kappatot;
    }
    
     kelmin = 1.e-5*kelmin;
     
  
    //ZSLOOP(-N1G,N1+N1G-1,-N2G,N2G+N2-1,-N3G,N3G+N3-1){
    for(iglob=-N1G;iglob<mpi_ntot[1]+N1G;iglob++) {
      for(jglob=-N2G;jglob<mpi_ntot[2]+N2G;jglob++) {
        for(kglob=-N3G;kglob<mpi_ntot[3]+N3G;kglob++) {
      
          rancval = ranc(0);

          i = iglob-mpi_startn[1];
          j = jglob-mpi_startn[2];
          k = kglob-mpi_startn[3];
          
          if(i<-N1G ||
             j<-N2G ||
             k<-N3G ||
             i>=N1+N1G ||
             j>=N2+N2G ||
             k>=N3+N3G){
            continue;
          }
          //Initialize Temperatures
          Ttot   =  (gam-1.)*p[i][j][k][UU]/p[i][j][k][RHO];
          
          Tel4 = fel0*Ttot;
          Tel5 = fel0*Ttot;
          
          coord(i,j,k,CENT,X);
          get_phys_coord(i,j,k,&r,&th,&phi);

          uel4 = fel0*p[i][j][k][UU]; //*(1.+.2*sin(2.*M_PI*r));
          uel5 = fel0*p[i][j][k][UU]*(1. + 4.e-2*(rancval-0.5)); //*(1.+.2*sin(2.*M_PI*r));
          uels = fel0*p[i][j][k][UU];
       
         
          //Initialize Entropies
          
          kappa =(gam-1.)*fel0*p[i][j][k][UU]*pow(p[i][j][k][RHO],-gam);
          kappae = (game-1.)*fel0*p[i][j][k][UU]*pow(p[i][j][k][RHO],-game);
          kappae4 = (game4-1.)*uel4*pow(p[i][j][k][RHO],-game4);
          kappae5 = (game5-1.)*uel5*pow(p[i][j][k][RHO],-game5);
                    
#if( eCOND || eHEAT || DOKTOT )
          kappatot = (gam-1.)*p[i][j][k][UU]*pow(p[i][j][k][RHO],-gam);
          p[i][j][k][KTOT] = kappatot;
#endif
          
#if( eCOND || eHEAT )
          p[i][j][k][KELDIS] = kappae;
          p[i][j][k][KEL4] = kappae4;
          p[i][j][k][KEL4A] = kappae4;
          p[i][j][k][KEL4B] = kappae4;
          p[i][j][k][KEL4C] = kappae4;
          p[i][j][k][KEL4D] = kappae4;
          p[i][j][k][KEL4E] = kappae4;
          p[i][j][k][KEL5] = kappae5;
#endif
          
#if( DOFLR )
          //initially nothing is floored
          p[i][j][k][FLR] = 0;
#endif
        }
      }
    }
  
     bound_prim(p);
     
  
    
    
}


/*
For changing to transmissive polar boundaries:
The most important considerations in  electron.c are:
init_entropy() initializes electron properties including at boundaries:
You should ensure that:
Electron properties remain well-behaved near poles with transmissive boundaries
Floor values (fixupuel()) are appropriate for transmissive conditions
Heat conduction (Conduction()) handles the pole regions properly with the new boundary conditions
However, you won't need to modify this file directly for implementing transmissive boundaries - 
it will just need careful testing to ensure electron physics remains stable with the new boundary conditions.

Here are the specific tests you should run to verify electron behavior with transmissive boundaries:
Conservation Tests:
// Monitor these quantities near poles:
kappatot = (gam-1.)*p[i][j][k][UU]*pow(p[i][j][k][RHO],-gam);  // Total entropy
kappae = (game-1.)*fel0*p[i][j][k][UU]*pow(p[i][j][k][RHO],-game);  // Electron entropy
Check:
Total electron energy conservation across polar boundaries
No artificial accumulation/loss of electron energy near poles
Entropy conservation for electrons (in fixup.c)
if(pv[UU] < uuflr) {
    pv[UU] = uuflr;
    dofloor = 1;
}
Floor Activation Monitoring:
Monitor:
Frequency of floor activation near poles
Compare with previous reflective boundary behavior
Check if floors are being hit more often with transmissive conditions
Heat Conduction Tests:
#if(eCOND)
    //if the heat flux is a nan, reset it to zero
    if(isnan(pv[PHI])) {
        pv[PHI] = 0.0;
        dofloor = 1;
    }
#endif
Verify:
Heat flux remains well-behaved across polar boundaries
No unphysical temperature gradients develop
Monitor  PHI (heat flux) for NaN values near poles
Temperature Profile Tests:
- Plot electron temperature profiles across θ = 0 and θ = π
- Check for discontinuities or unphysical jumps
- Compare with analytical solutions where available

Output Diagnostics:
- Monitor KTOT (total entropy)
- Track KELDIS (electron entropy)
- Check KEL4/KEL5 (different electron components)
- Watch FLR (floor activation indicator)


*/
