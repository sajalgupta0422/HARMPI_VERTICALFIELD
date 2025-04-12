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

/*
MPI Boundary Communication (for parallel runs)
Function	                                                                              Purpose	                                                            Comments
bound_mpi_dim(int dim, int ispack, double prim[][N2M][N3M][NPR], int pflag[][N2M][N3M])	Boundary communication in direction dim (1=x1, 2=x2, 3=x3). 
                                                                                        Does packing and sending ghost zone data using MPI.	            ispack probably tells whether it is packing (send) or unpacking (receive).
pack_prim(int ispack, int dim, int isup, double prim[][N2M][N3M][NPR], double *mpi_buf)	Packs or unpacks primitive variable data for MPI exchange 
                                                                                        at boundary in a given dimension and side.	                    isup=0 for lower boundary, isup=1 for upper boundary.
pack_particles(...)	                                                                    Same idea, but for particles data (if particles are evolved). 
                                                                                                              Not our focus now.	
pack_pflag(...)	                                                                        Same, but for packing/unpacking a flag array pflag 
                                                                                        (likely marking zones needing special treatment).
*/
void bound_mpi_dim(int dim, int ispack, double prim[][N2M][N3M][NPR], int pflag[][N2M][N3M]); //Handle MPI ghost zone exchange in direction dim	Needed for multi-processor runs
int pack_prim(int ispack, int dim, int isup, double prim[][N2M][N3M][NPR], double *mpi_buf ); //Packs prim array for MPI send/recv	Memory layout stuff
int pack_particles(int ispack, int dim, int isup, double xp[][NDIM], double *mpi_buf, int *count_particles ); //Packs particles	If particles exist in the simulation
int pack_pflag(int ispack, int dim, int isup, int pflag[][N2M][N3M], int *mpi_buf ); //Packs pflag array for MPI	Same idea for flags
void bound_x1dn(double prim[][N2M][N3M][NPR] ); //Apply boundary condition at lower x1 boundary	: Inner r
void bound_x2dn(double prim[][N2M][N3M][NPR] ); //Apply BC at lower x2 boundary:	Lower θ
void bound_x2dn_polefix(double prim[][N2M][N3M][NPR] ); //Special treatment near pole at θ=0:	Lower θ pole fix
void bound_x3dn(double prim[][N2M][N3M][NPR] ); //Apply BC at lower x3 boundary:	Lower φ
void bound_x1up(double prim[][N2M][N3M][NPR] ); //Apply BC at upper x1 boundary: Outer r
void bound_x2up(double prim[][N2M][N3M][NPR] ); //Apply BC at upper x2 boundary: Upper θ
void bound_x2up_polefix(double prim[][N2M][N3M][NPR] );  //Special treatment near pole at θ=π:	Upper θ pole fix
void bound_x3up(double prim[][N2M][N3M][NPR] ); //Apply BC at upper x3 boundary: Upper φ


/* bound array containing entire set of primitive variables */

void bound_x1dn(double prim[][N2M][N3M][NPR] ) //Purpose: Apply the physical boundary condition at the lower x1 boundary → usually inner radial boundary. 
/*Argument:
    prim[][N2M][N3M][NPR] → This is the full 3D array of primitive variables,or you can say 4D array holding primitive variables
    Dimensions: N1M = N1 + 2 * NG (ghost zones on x1); N2M = N2 + 2 * NG (ghost zones on x2); N3M = N3 + 2 * NG (ghost zones on x3)
    NPR = number of primitive variables at each zone. 
    This function updates prim values at ghost zones on the inner x1 boundary (lower radial side)*/
{
  int i,j,k,m,jref ; //i,j,k = x1,x2 and x3 index; m = variable index (runs over NPR); jref = to reflect accross pole (used sometimes)
  void inflow_check(double *pr, int ii, int jj, int kk, int type ); //This is a key utility function. Its role is:
//Prevent unphysical inflow into the domain at inner or outer boundaries. Typically kills radial inflow if not allowed (by setting v^r ≥ 0 at inner boundary, or ≤ 0 at outer boundary).
  struct of_geom geom ; //Declaring a local geometry structure to store metric coefficients (will be filled later, at a given (i,j,k) cell): Metric components, Coordinate mappings, Volume factors, Likely comes from set_geometry(...)
  
  int iNg, jNg, kNg; //Just a local alias for NG (number of ghost zones) = 3 
  
#if(N1!=1) // True
  if (!is_physical_bc(1, 0)) return;  //True; int is_physical_bc( int dim, int isup ) dim → which dimension: 1 → x1 (radial), 2 → x2 (theta), 3 → x3 (phi)
                                                                                  //isup → is it upper boundary? (1) or lower boundary? (0)

  /* inner r boundary condition */
  for(j=0;j<N2;j++)
  {
    for(k=0; k<N3;k++)  //!!!ATCH: make sure don't need to expand to transverse ghost cells
    {
#if( RESCALE ) //This rescales variables at i=0 (physical boundary) → if code wants rescaling, but False as RESCALE is 0 in decs.h
      get_geometry(0,j,k,CENT,&geom) ;
      rescale(prim[0][j][k],FORWARD, 1, 0,j,k,CENT,&geom) ;
#endif
      
      for (iNg=-N1G; iNg<0; iNg++) //iNg runs from -N1G to -1 → Ghost cells at lower radial boundary
      {
#if (PERIODIC==1)
        PLOOP prim[iNg][j][k][m] = prim[N1+iNg][j][k][m]; //Periodic boundary: Wrap around to i=N1+iNg
        pflag[iNg][j][k] = pflag[N1+iNg][j][k];
        //isdis[iNg][j][k] = isdis[N1+iNg][j][k];
#elif(DIRICHLET==1)
        PLOOP prim[iNg][j][k][m] = pbound[iNg][j][k][m]; // Dirichlet boundary: Use pre-specified pbound[iNg][j][k][m]
        pflag[iNg][j][k] = pflag[0][j][k] ;
        //isdis[iNg][j][k] = isdis[0][j][k];
#else //outflow Outflow boundary allows matter to leave (flow out of) the domain freely. But it tries to prevent inflow (stuff coming back in from the ghost zones into the domain). Copy boundary-adjacent value at i=0
        PLOOP prim[iNg][j][k][m] = prim[0][j][k][m]; //PLOOP → Loop over m=0 to NPR-1 → all primitive variables.
        pflag[iNg][j][k] = pflag[0][j][k] ; //Also copying pflag for each ghost zone: pflag[iNg][j][k] = pflag[0][j][k]; pflag likely marks things like: floor-flag, shocked, atmosphere ,etc.
        //isdis[iNg][j][k] = isdis[0][j][k];
#endif
      }
      
#if( RESCALE ) //False as RESCALE is 0 in decs.h
      for (iNg = -N1G; iNg<=0; iNg++)
      {
        get_geometry(iNg,j,k,CENT,&geom) ;
        rescale(prim[iNg][j][k],REVERSE, 1, iNg,j,k,CENT,&geom) ;
      }
#endif
    }
  }

  /* make sure there is no inflow at the inner boundary  
   Outflow boundary condition should prevent inflow into the domain. 
   In MHD codes: Outflow boundary ≠ setting velocities to zero. Instead: Let fluid exit freely. But if any fluid tries to enter from ghost zones → kill the inflow.*/
  if(1!=INFLOW) { //INFLOW IS 0 in decs.h, so TRUE
    for(i=-N1G;i<=-1;i++)  for(j=-N2G;j<N2+N2G;j++) for(k=-N3G;k<N3+N3G;k++) // This loop ensures: No inflow through any ghost cell at lower x1.
    {
      //!!!ATCH: eliminated one loop that seemed excessive. verify.
      inflow_check(prim[i][j][k],i,j,k,0) ; //0 stands for -x1 boundary , so do not allow ucon[1] > 0, as it Checks v^1 (radial velocity) in ghost zones. If v^1 points into the domain (i.e., v^1 > 0 at -x1 boundary) → set it to zero or flip it negative. (see function defn below)
    }
  }
#endif
  
}

/*1. What exactly is "Outflow" boundary condition in this context?

→ Outflow BC does not set velocity to zero.
→ It allows matter to leave (flow out of) the domain freely.
→ But it tries to prevent inflow (stuff coming back in from the ghost zones into the domain).
How does it achieve that?

    Copy the primitive variables in the first active zone to all ghost zones near the boundary.

    If after copying, the radial velocity points inward (into domain), we flip it to zero or outward.
    → This part may be elsewhere in the code (sometimes via fixup functions).

Why is prim[0][j][k][m] the first active (physical) zone in bound_x1dn(..) function?
i-index	              Meaning
-3,-2,-1	        Ghost zones (inner)
0 → N1-1	        Active zones (physical)
N1 → N1+2	        Ghost zones (outer)
→ So if N1=16, then valid physical i runs from 0 to 15.
Thus:
    prim[-3][j][k][m] → ghost zone
    prim[0][j][k][m] → first physical zone (inner boundary)
    prim[N1-1][j][k][m] → last physical zone (outer boundary)
    prim[N1][j][k][m] to prim[N1+2][j][k][m] → outer ghost zones
*/
void bound_x1up(double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  void inflow_check(double *pr, int ii, int jj, int kk, int type );
  struct of_geom geom ;
  
  int iNg, jNg, kNg;
#if(N1!=1)
  if (!is_physical_bc(1, 1)) return;

  /* Outer r boundary condition */
  for(j=0;j<N2;j++)
  {
    for(k=0; k<N3;k++)  //!!!ATCH: make sure don't need to expand to transverse ghost cells
    {
#if( RESCALE )
      get_geometry(N1-1,j,k,CENT,&geom) ;
      rescale(prim[N1-1][j][k],FORWARD, 1, N1-1,j,k,CENT,&geom) ;
#endif
      
      for (iNg=0; iNg<N1G; iNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[N1+iNg][j][k][m] = prim[iNg][j][k][m];
        pflag[N1+iNg][j][k] = pflag[iNg][j][k] ;
        //isdis[N1+iNg][j][k] = isdis[iNg][j][k];
#elif(DIRICHLET==1)
        PLOOP prim[N1+iNg][j][k][m] = pbound[N1+iNg][j][k][m];
        pflag[N1+iNg][j][k] = pflag[N1-1][j][k] ;
        //isdis[N1+iNg][j][k] = isdis[N1-1][j][k];
        
#else //outflow
        PLOOP prim[N1+iNg][j][k][m] = prim[N1-1][j][k][m];
        pflag[N1+iNg][j][k] = pflag[N1-1][j][k] ;
        //isdis[N1+iNg][j][k] = isdis[N1-1][j][k];
#endif
      }
#if( RESCALE )
      for (iNg= 0; iNg<N1G+1; iNg++) //!!!ATCH: added +1 to N1G to ensure that all ghost cells are looped over
      {
        get_geometry(N1-1+iNg,j,k,CENT,&geom) ;
        rescale(prim[N1-1+iNg][j][k],REVERSE, 1, N1-1+iNG,j,k,CENT,&geom) ;
      }
#endif
    }
  }
  /* make sure there is no inflow at the outer boundary */
  if(1!=INFLOW) {
    for(i=N1;i<=N1+N1G-1;i++)  for(j=-N2G;j<N2+N2G;j++) for(k=-N3G;k<N3+N3G;k++)
    {
      //!!!ATCH: eliminated one loop that seemed excessive. verify.
      inflow_check(prim[i][j][k],i,j,k,1) ; //1 stands for +x1 boundary
    }
  }

#endif
}

void bound_x2dn_polefix( double prim[][N2M][N3M][NPR] ) // It provides a special fix near the inner θ boundary (θ → 0), i.e., the pole. (e.g., terms like 1/sinθ in metric components become divergent)
{
  /*
  Physical Motivation
    Avoid Metric Singularities: Near θ=0, terms like sinθ ≈ 0 appear in the metric. If variables like U2 or density have large gradients here, calculations of fluxes or source terms (e.g., ΓνλμΓνλμ​) can become unstable.
    Example: A large U2 (θ-velocity) in j=0 could lead to division by sinθ ≈ 0 in the connection coefficients
  Why Skip Magnetic Fields?
    Divergence Constraint: Magnetic fields must satisfy ∇⋅B=0∇⋅B=0. Arbitrary overwriting of B1/B2/B3 could violate this.
    Specialized Treatment: Magnetic fields near poles are likely handled via constrained transport (CT) or staggered grid algorithms elsewhere.    
  However, this function introduces localized violations of the GRMHD equations (continuity, energy-momentum conservation) in the polar regions to prioritize numerical stability.
  For instance, Continuity Equation: ∇_μ(ρu^μ)=0
    Issue: By overwriting ρ and U2 independently, the divergence term ∇_μ(ρu^μ) is no longer guaranteed to vanish in cells j=0,1.
    Example:
        Suppose ρ is copied from j=2, while U2 is interpolated.
        The original j=2 cell might satisfy ∇_μ(ρu^μ)=0, but the modified j=0,1 cells likely do not. 
    But, 
        Numerical Dissipation: The code’s shock-capturing schemes (e.g., HLL fluxes) implicitly add dissipation, damping these errors.
        Time Integration: Over multiple timesteps, the errors are integrated into the solution but remain confined near the pole.
    */
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;
  
#if(POLEFIX && POLEFIX < N2/2 && BL) //TRUE as polefix = 2 < N2/2 and BL = 3
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(2, 0)) return; //    Only MPI processes containing the physical θ=0 boundary execute this function. Avoids redundant operations on processes not adjacent to the pole.

  //copy all densities and B^phi in; interpolate linearly transverse velocity
  jref = POLEFIX;
  for(i=-N1G;i<N1+N1G;i++) {
    for(k=-N3G;k<N3+N3G;k++) {
      for(j=0;j<jref;j++) { //Looping over only j = 0.. jref = 2, so only 0 and 1. So basically first two cells near θ ~ 0
        PLOOP {
          if(m==B1 || m==B2 || (N3>1 && m==B3)) //Checking if primitive variables are magnetic field
            //don't touch magnetic fields
            continue;
          else if(m==U2) {
            //linear interpolation of transverse velocity (both poles) Creates a smooth transition from j = 2 to j = 0 or 1
            /* Example: For j = 0:
              U2(j=0)=(0.5/2.5)*U2(j=2) = 0.2*U2(j=2)
              Prevents abrupt velocity jumps near the pole.*/
            prim[i][j][k][m] = (j+0.5)/(jref+0.5) * prim[i][jref][k][m];
          }
          else {
            //everything else copy (both poles)     Directly copy values from j = 2 to j = 0, 1. Effect: Dampens fluctuations near the pole by enforcing uniformity.
            prim[i][j][k][m] = prim[i][jref][k][m];
          }
        }
      }
    }
  }
#endif
}

void bound_x2up_polefix( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ; 
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(POLEFIX && POLEFIX < N2/2 && BL)
  //only do anything if physically in an MPI process that touches the outer pole
  if (!is_physical_bc(2, 1)) return;

  //copy all densities and B^phi in; interpolate linearly transverse velocity
  jref = POLEFIX;
  for(i=-N1G;i<N1+N1G;i++) {
    for(k=-N3G;k<N3+N3G;k++) {
      for(j=0;j<jref;j++) {
        PLOOP {
          if(m==B1 || m==B2 || (N3>1 && m==B3))
            //don't touch magnetic fields
            continue;
          else if(m==U2) {
            //linear interpolation of transverse velocity (both poles)
            prim[i][N2-1-j][k][m] = (j+0.5)/(jref+0.5) * prim[i][N2-1-jref][k][m];
          }
          else {
            //everything else copy (both poles)
            prim[i][N2-1-j][k][m] = prim[i][N2-1-jref][k][m];
          }
        }
      }
    }
  }
#endif
  
}

/* polar BCs */
//inner theta boundary 
void bound_x2dn( double prim[][N2M][N3M][NPR] ) //Purpose: Apply boundary condition at lower θ boundary (near θ = 0, i.e., north pole)
{
  int i,j,k,m,jref ; //i,j,k = x1,x2 and x3 index; m = variable index (runs over NPR); jref = to reflect accross pole (used sometimes)
  struct of_geom geom ; //Declaring a local geometry structure to store metric coefficients (will be filled later, at a given (i,j,k) cell): Metric components, Coordinate mappings, Volume factors, Likely comes from set_geometry(...)
  int iNg, jNg, kNg; //Just a local alias for NG (number of ghost zones) = 3 

#if(N2!=1) //TRUE; Only apply this if there is more than 1 θ-zone.
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(2, 0)) return; /* Very important. This checks: Are we really sitting on the lower θ boundary in this MPI process? 
  dim=2 → means θ direction. isup=0 → means lower boundary (θ=0). If not → exit the function. This is MPI-safe design*/

  bound_x2dn_polefix(prim);
  
  for (i=-N1G; i<N1+N1G; i++) //All radial zones (i)
  {
    for (k=-N3G; k<N3+N3G; k++) //All φ zones (k)
    {
      for (jNg=-N2G; jNg<0; jNg++) //All θ ghost zones near θ=0 (jNg=-3,-2,-1 if N2G=3)
      {
#if (PERIODIC==1) //Copy data from opposite θ boundary.
        PLOOP prim[i][jNg][k][m] = prim[i][N2+jNg][k][m];
        pflag[i][jNg][k] = pflag[i][N2+jNg][k]; //N2+jNg is mapped to upper θ side.
        //isdis[i][jNg][k] = isdis[i][N2+jNg][k];
#elif(DIRICHLET==1) //Fixed boundary values from pbound array.
        PLOOP prim[i][jNg][k][m] = pbound[i][jNg][k][m];
        pflag[i][jNg][k] = pflag[i][-jNg-1][k];
        //isdis[i][jNg][k] = isdis[i][-jNg-1][k];
        
#elif (OUTFLOW==1) //Let the stuff flow out of θ=0 if it wants to, but don't artificially enforce inflow.
        PLOOP prim[i][jNg][k][m] = prim[i][0][k][m]; //Just copy whatever is in the first active θ zone (j=0) into the ghost zones at jNg=-3,-2,-1.
        pflag[i][jNg][k] = pflag[i][0][k]; //Copy the pflag (physical flag) from active zone.
        //isdis[i][jNg][k] = isdis[i][0][k];
#else //symmetric/asymmetric ; Reflective BCs enforce symmetry across the polar boundaries (θ = 0 and θ = π)
/*This is a clean and typical reflecting boundary at the pole:
    Even symmetry for scalars
    Odd symmetry for θ-velocity and magnetic field 
     Result: Variables like density (RHO), φ-velocity (U3), and radial velocity (U1) are mirrored */
        PLOOP prim[i][jNg][k][m] = prim[i][-jNg-1][k][m]; //For a ghost cell at jNg (negative index near θ=0), copy values from the mirrored physical cell at -jNg-1. Example: Ghost cell j = -1 (θ < 0) mirrors physical cell j = 0 (θ > 0).
        pflag[i][jNg][k] = pflag[i][-jNg-1][k];
        //isdis[i][jNg][k] = isdis[i][-jNg-1][k];
#endif
        
      }
    }
  }
  
  /* polar BCs */
  /* make sure b and u are antisymmetric at the poles */
  if(BL){ //TRUE as BL = 3
    for(i=-N1G;i<N1+N1G;i++) {
      for(k=-N3G;k<N3+N3G;k++) {
        for(j=-N2G;j<0;j++) {
          prim[i][j][k][U2] *= -1. ; //Flip sign of θ-velocity. Why? Because θ-direction reverses across the pole. Velocity U2 is a contravariant vector; reflection across θ=0 flips its direction
          prim[i][j][k][B2] *= -1. ; //Flip sign of θ-magnetic field to maintain consistency with Maxwell’s equations.
          /*This mimics a "mirror" at the pole, ensuring no net flux of mass, momentum, or magnetic field through θ = 0/π. 
          The fluid and fields "bounce" symmetrically off the pole. 
          When the solver computes fluxes or gradients using ghost cells, the mirrored/antisymmetric values enforce zero net flux through the pole.
          Example: If fluid flows toward θ=0 (negative U2 in physical cells), ghost cells have positive U2, canceling the flux at the boundary*/
        }
      }
    }
  }

#endif

}

/*
Staggered Grid Positioning
In staggered grids, variables are stored at different locations (cell centers, faces, edges). This affects BC implementation:
  Cell-Centered Variables (e.g., density RHO, internal energy UU):
        Stored at cell centers
        Symmetric reflection (no sign flip) unless explicitly modified.
    Face-Centered Variables (e.g., magnetic fields B2 on θ-faces):
        Staggered at cell faces (edges in θ-direction).
        Reflection must account for their position.
        Example: B2 at θ-face j+1/2 is mirrored to j-1/2 with a sign flip.

Why Staggered Grids Matter:
    Fluxes (e.g., F2) are computed at faces. If B2 is face-centered, flipping its sign in ghost cells ensures proper antisymmetry when calculating divergence terms like ∇·B.*/
/* polar BCs 

4. Impact Near θ ≈ 0
Near the pole (θ ≈ 0), reflective BCs ensure:
    Velocity:
        Fluid approaching θ=0 (negative U2) is reflected back (positive U2 in ghost cells).
        This prevents unphysical inflow/outflow through the pole.

    Magnetic Fields:
        B2 (θ-component) reverses sign, preserving the solenoidal condition (∇·B = 0).

    Stability:
        Without reflection, singularities at θ=0/π (due to coordinate terms like sinθ) could cause numerical instabilities.

5. Code-Specific Insights
    Metric Handling (BL=3):
    The code uses Boyer-Lindquist (BL) coordinates for black holes. Reflective BCs here must respect the curved spacetime geometry:
        The √-g term (metric determinant) is implicitly handled by the grid structure (phys_coords and a_gdet arrays).
        Sign flips for U2/B2 ensure vectors transform correctly under reflection in BL coordinates.

    Pole Stabilization (POLEFIX=2):
    The code modifies variables in the 2 cells nearest the pole to avoid numerical issues (e.g., division by sinθ ≈ 0).*/
//outer theta boundary
void bound_x2up( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;
  
#if(N2!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(2, 1)) return;

  bound_x2up_polefix(prim);

  for (i=-N1G; i<N1+N1G; i++)
  {
    for (k=-N3G; k<N3+N3G; k++)
    {
      //outer theta boundary
      for (jNg=0; jNg<N2G; jNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[i][N2+jNg][k][m] = prim[i][jNg][k][m];
        pflag[i][N2+jNg][k] = pflag[i][jNg][k];
        //isdis[i][N2+jNg][k] = isdis[i][jNg][k];
#elif(DIRICHLET==1)
        PLOOP prim[i][N2+jNg][k][m] = pbound[i][N2+jNg][k][m];
        pflag[i][N2+jNg][k] = pflag[i][N2-jNg-1][k];
        //isdis[i][N2+jNg][k] = isdis[i][N2-jNg-1][k];
#elif (OUTFLOW==1)
        PLOOP prim[i][N2+jNg][k][m] = prim[i][N2-1][k][m];
        pflag[i][N2+jNg][k] = pflag[i][N2-1][k];
        //isdis[i][N2+jNg][k] = isdis[i][N2-1][k];
#else //symmetric/asymmetric
        PLOOP prim[i][N2+jNg][k][m] = prim[i][N2-jNg-1][k][m];
        pflag[i][N2+jNg][k] = pflag[i][N2-jNg-1][k];
        //isdis[i][N2+jNg][k] = isdis[i][N2-jNg-1][k];
        
#endif
      }
    }
  }
  
  /* polar BCs */
  /* make sure b and u are antisymmetric at the poles */
  //SEANMARK: do we need the "!= ATM_TEST" condition here? This should not be activated for N2==1 as in ATM_TEST
  if(BL && WHICHPROBLEM != ATM_TEST){
    for(i=-N1G;i<N1+N1G;i++) {
      for(k=-N3G;k<N3+N3G;k++) {
        for(j=N2;j<N2+N2G;j++) {
          prim[i][j][k][U2] *= -1. ;
          prim[i][j][k][B2] *= -1. ;
        }
      }
    }
  }

#endif
}

void bound_x3dn( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(N3!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(3, 0)) return;

  /* phi BCs */
  //inner phi-boundary
  for (i=-N1G; i<N1+N1G; i++)
  {
    for (j=-N2G; j<N2+N2G; j++)
    {
      for (kNg=-N3G; kNg<0; kNg++)
      {
#if (PERIODIC==1)
        PLOOP prim[i][j][kNg][m] = prim[i][j][N3+kNg][m];
        pflag[i][j][kNg] = pflag[i][j][N3+kNg];
        //isdis[i][j][kNg] = isdis[i][j][N3+kNg];
#elif(DIRICHLET==1)
        PLOOP prim[i][j][kNg][m] = pbound[i][j][kNg][m];
        pflag[i][j][kNg] = pflag[i][j][0];
        //isdis[i][j][kNg] = isdis[i][j][0];
        
#elif (OUTFLOW==1)
        PLOOP prim[i][j][kNg][m] = prim[i][j][0][m];
        pflag[i][j][kNg] = pflag[i][j][0];
        //isdis[i][j][kNg] = isdis[i][j][0];
#else
        //periodic by default
        PLOOP prim[i][j][kNg][m] = prim[i][j][N3+kNg][m];
        pflag[i][j][kNg] = pflag[i][j][N3+kNg];
        //isdis[i][j][kNg] = isdis[i][j][N3+kNg];
#endif
      }
    }
  }
#endif

}

void bound_x3up( double prim[][N2M][N3M][NPR] )
{
  int i,j,k,m,jref ;
  struct of_geom geom ;
  int iNg, jNg, kNg;

#if(N3!=1)
  //only do anything if physically in an MPI process that touches the inner pole
  if (!is_physical_bc(3, 1)) return;

  /* phi BCs */
  //outer phi-boundary
  for (i=-N1G; i<N1+N1G; i++)
  {
    for (j=-N2G; j<N2+N2G; j++)
    {
      for (kNg=-N3G; kNg<0; kNg++)
      {
        for (kNg=0; kNg<N3G; kNg++)
        {
#if (PERIODIC==1)
          PLOOP prim[i][j][N3+kNg][m] = prim[i][j][kNg][m];
          pflag[i][j][N3+kNg] = pflag[i][j][kNg];
          //isdis[i][j][N3+kNg] = isdis[i][j][kNg];
#elif(DIRICHLET==1)
          PLOOP prim[i][j][N3+kNg][m] = pbound[i][j][N3+kNg][m];
          pflag[i][j][N3+kNg] = pflag[i][j][N3-1];
          //isdis[i][j][N3+kNg] = isdis[i][j][N3-1];
          
#elif (OUTFLOW==1)
          PLOOP prim[i][j][N3+kNg][m] = prim[i][j][N3-1][m];
          pflag[i][j][N3+kNg] = pflag[i][j][N3-1];
          //isdis[i][j][N3+kNg] = isdis[i][j][N3-1];
#else
          //periodic by default
          PLOOP prim[i][j][N3+kNg][m] = prim[i][j][kNg][m];
          pflag[i][j][N3+kNg] = pflag[i][j][kNg];
          //isdis[i][j][N3+kNg] = isdis[i][j][kNg];
#endif
        }
      }
    }
  }
#endif
}

void bound_prim( double prim[][N2M][N3M][NPR] )
{
  /*
  It is a master boundary function that applies boundary conditions to the primitive variables (prim) in all three spatial dimensions (radial, polar, azimuthal). It combines:
    MPI Communication for inter-process boundary updates (for parallel simulations).
    Physical BCs (e.g., reflective, outflow) at domain edges.
    Overlap of Computation and Communication to improve efficiency.
    Its job?: Make sure all boundary zones (ghost zones) have correct values.
  */
  int ispack, dim; //dim → which coordinate direction: x1 (radial), x2 (theta), x3 (phi); ispack → flag if we are "packing" or "unpacking" data (MPI comm)
  
  //MPIMARK: could be optimized by individually passing corner zones:
  //         then, speed up by not doing comm dimension by dimension

  /*
  1. Dimension Loop Structure
      The function processes each spatial dimension sequentially (x1, x2, x3). For each dimension:
      Step 1: Initiate MPI communication to exchange ghost cell data with neighboring MPI processes.
      Step 2: Apply physical BCs to the local domain’s boundaries (non-MPI edges).
      Step 3: Finalize MPI communication and unpack received data into ghost cells.
  */
  //x1-dim
  //packing, putting send and receive requests
  dim = 1;
  ispack = 1; //ispack = 1 triggers packing/sending ghost cell data to neighboring MPI processes. Its like, Hey, pack the data & send it to neighbor MPI ranks
  bound_mpi_dim(dim, ispack, prim, pflag); //bound_mpi_dim handles the MPI send/receive for the x1 ghost zones.
  /*while waiting for MPI comm to complete, do physical boundaries
  Apply physical BCs to the inner (x1dn) and outer (x1up) radial boundaries of the local domain (not handled by MPI)
  */
  bound_x1dn(prim);
  bound_x1up(prim);
  //waiting for comm to complete and unpacking. Now MPI should be done.
  ispack = 0; //ispack = 0 triggers unpacking received ghost cell data from neighboring processes. That is, Unpack received ghost zone data back into prim.
  bound_mpi_dim(dim, ispack, prim, pflag);

  //x2-dim
  //packing, putting send and receive requests
  dim = 2;
  ispack = 1;
  bound_mpi_dim(dim, ispack, prim, pflag);
  //while waiting for MPI comm to complete, do physical boundaries
  bound_x2dn(prim);
  bound_x2up(prim);
  //waiting for comm to complete and unpacking. 
  ispack = 0;
  bound_mpi_dim(dim, ispack, prim, pflag);
  

  //x3-dim
  //waiting for comm to complete and unpacking
  dim = 3;
  ispack = 1;
  bound_mpi_dim(dim, ispack, prim, pflag);
  //while waiting for MPI comm to complete, do physical boundaries
  bound_x3dn(prim);
  bound_x3up(prim);
  //waiting for comm to complete and unpacking
  ispack = 0;
  bound_mpi_dim(dim, ispack, prim, pflag);


  /*Key Concepts
1. Overlapping Computation and Communication
    MPI Non-Blocking Calls: While waiting for MPI communication to complete (e.g., between ispack=1 and ispack=0), the code applies physical BCs to local boundaries. This hides communication latency, improving performance.
2. Ghost Cell Population
    MPI Ghost Cells: Updated via bound_mpi_dim (data from neighboring MPI processes).
    Physical Ghost Cells: Updated via bound_x*dn/up (e.g., reflective/outflow BCs).
3. Dimension-Wise Processing
    Sequential Handling: Each dimension (x1, x2, x3) is processed separately to simplify logic, though this may not optimize corner zones.
4. Role of ispack
    ispack=1: Initiate MPI communication (pack/send ghost cells).
    ispack=0: Finalize MPI communication (receive/unpack ghost cells)*/
}

/*
Explaing bound_prim(..) function with an example
Example: Polar Boundary Handling with MPI
Imagine a 2D grid in (r, θ) split across 2 MPI processes along the θ-direction:
    Process A: Handles θ = [0, π/2] (lower hemisphere).
    Process B: Handles θ = [π/2, π] (upper hemisphere).

Each process has 3 ghost cells (NG=3) at their θ-boundaries. Here’s how bound_prim works for the θ-dimension:
Step 1: MPI Communication (θ-Dimension)
    Goal: Exchange ghost cell data between Process A and B at their shared boundary (θ = π/2).
    Action:
        Process A sends its upper ghost cells (θ = π/2 - Δθ, π/2 - 2Δθ, π/2 - 3Δθ) to Process B.
        Process B sends its lower ghost cells (θ = π/2 + Δθ, π/2 + 2Δθ, π/2 + 3Δθ) to Process A.

    Code Execution:
    dim = 2; ispack = 1;
    bound_mpi_dim(dim, ispack, prim, pflag); // Start MPI communication

Step 2: Apply Physical BCs (Local Poles)
While MPI data is in transit:

    Process A applies lower θ BC (θ=0):
    bound_x2dn(prim); // Calls reflective BCs + polefix for θ=0; Mirrors cells near θ=0 into ghost zones, flips U2/B2, and stabilizes with bound_x2dn_polefix.
    Process B applies upper θ BC (θ=π):
    bound_x2up(prim); // Calls reflective BCs + polefix for θ=π; Mirrors cells near θ=π into ghost zones, flips U2/B2, and stabilizes.

Step 3: Finalize MPI Communication
    Process A receives data from Process B and populates its upper ghost cells (θ > π/2).
    Process B receives data from Process A and populates its lower ghost cells (θ < π/2).

    Code Execution:
    ispack = 0;
    bound_mpi_dim(dim, ispack, prim, pflag); // Unpack received MPI data

Visualization
After MPI Communication:
    Process A’s Ghost Cells (θ > π/2): Filled with data from Process B’s physical cells (θ = π/2 + Δθ, etc.).
    Process B’s Ghost Cells (θ < π/2): Filled with data from Process A’s physical cells (θ = π/2 - Δθ, etc.).

After Physical BCs:
    Process A’s Ghost Cells (θ < 0): Mirrored from θ=0 with U2/B2 sign flips and stabilization.
    Process B’s Ghost Cells (θ > π): Mirrored from θ=π with U2/B2 sign flips and stabilization.

Key Takeaways
    MPI Handles Internal Boundaries:
        Exchanges data between processes at θ = π/2 (shared edge). Ensures smooth transitions across subdomain boundaries.

    Physical BCs Handle Domain Edges:
        θ=0 (Process A) and θ=π (Process B) use reflective BCs with stabilization.

    Overlap Improves Efficiency:
        MPI communication (latency-heavy) overlaps with physical BC computation (CPU-heavy).

Why This Works
    Ghost Cells as Buffers: MPI-populated ghost cells allow each process to compute fluxes as if it had the full grid.

    Pole Stabilization: polefix functions prevent numerical issues at θ=0/π without affecting MPI communication.
This approach balances parallel scalability (via MPI) with physical accuracy (via BCs). Let’s examine bound_mpi_dim next to see how data is packed/unpacked!
THE END  */

//packs (ispack=1) or unpacks (ispack=0) the cells to be communicated along
//dimension (dim), either upper (isup=1) or lower (isup=0) boundary
//returns the number of items packed (count)
int pack_prim(int ispack, int dim, int isup, double prim[][N2M][N3M][NPR], double *mpi_buf )
{
  /*
  The pack_prim function is designed to handle packing and unpacking primitive variables for MPI communication
  Arguments:
    dim = which direction (x1, x2, x3)
    isup = 0 (lower neighbor) or 1 (upper neighbor)
    ispack = 1 → packing for sending;  0 → unpacking for receiving (but this routine only packs!)
  */
  int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count; //i,j,k → spatial indices, m → which primitive variable (e.g., density, velocity components), count → counter for 1D mpi buffer position

  //if true, ensure it has value of unity for the below to work
  if(ispack) ispack=1; //Redundant unless ispack comes in as something weird like 3 — ensures it’s 0 or 1 only.
  
  //do not do empty dimensions
  if(mpi_ntot[dim] == 1) return(0); // If this dimension has only 1 cell → nothing to pack.
  
  /* visual map for x1-boundary
  isup	              ispack	          istart calculation	                          Result
  0 (lower)	          1 (packing)	        -N1G + N1G = 0	            Packing lower boundary cells at i=0 to i=N1G-1
  0 (lower)	          0 (unpacking)	      -N1G + 0 = -N1G	            Unpacking into ghost zones at i=-N1G to i=-1
  1 (upper)	          1 (packing)	         N1 - N1G	                  Packing upper boundary cells at i=N1-N1G to i=N1-1
  1 (upper)	          0 (unpacking)	     N1 + 0 = N1	                Unpacking into ghost zones at i=N1 to i=N1+N1G-1

So, basically,
When packing (ispack=1), it's pulling data from active zones near the boundary.
When unpacking (ispack=0), it's writing data into ghost zones outside the active domain.
*/
  //figure out the range of indices to transfer
  //x1: in x1-dim transfer only ghost cells immediately adjacent to active cells
  if(1==dim){
    jstart=0; jstop=N2-1;
    kstart=0; kstop=N3-1;
    if(isup) istart=N1-ispack*N1G; else istart=-N1G+ispack*N1G;
    istop=istart+N1G-1;
  }
  //x2: in x2-dim, additionally trasfer the ghost cells that have been communicated in x1-dim
  else if(2==dim){
    istart=-N1G; istop=N1+N1G-1;
    kstart=0; kstop=N3-1;
    if(isup) jstart=N2-ispack*N2G; else jstart=-N2G+ispack*N2G;
    jstop=jstart+N2G-1;
  }
  //x3, in x3-dim, additionally trasfer the ghost cells that have been communicated in x1-dim and x2-dim
  else if(3==dim){
    istart=-N1G; istop=N1+N1G-1;
    jstart=-N2G; jstop=N2+N2G-1;
    if(isup) kstart=N3-ispack*N3G; else kstart=-N3G+ispack*N3G;
    kstop=kstart+N3G-1;
  }
  
  //initialize the counter of the number of doubles (un)packed
  count = 0;
  
  //packing; The classic triple loop to pack data into mpi_buf
  if(ispack){
    //OMPMARK
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) PLOOP {
      mpi_buf[count++] = prim[i][j][k][m];
    }
  }
  ///unpacking
  else {
    //OMPMARK
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) PLOOP {
      prim[i][j][k][m] = mpi_buf[count++];
    }
  }
  return(count); //Return number of packed elements:
}

#if(DOPARTICLES)
//packs (ispack=1) or unpacks (ispack=0) the particles to be communicated along
//dimension (dim), either upper (isup=1) or lower (isup=0) boundary
//returns the number of items packed, count = NDIM*[number of particles]. Eliminates the packed particles from xp[][].
int pack_particles(int ispack, int dim, int isup, double xp[][NDIM], double *mpi_buf, int *count_particles )
{
  int j,count,l;
  
  //do not do empty dimensions
  if(mpi_ntot[dim] == 1) {
    *count_particles = 0;
    return(0);
  }
  
  if (ispack) {
    //loop over all particles and ship away those that have left this tile
    for (count = 0, l = 0; l < myNp; ) {
      if( (xp[l][dim] >= mpi_stopx[dim] && isup) ||
          (xp[l][dim] <  mpi_startx[dim] && !isup) ) {
        //place the particle into the outgoing MPI buffer
        DLOOPA mpi_buf[count++] = xp[l][j];
        //eliminate the particle that we are shipping off
        //first, overwrite the shipped particle with the last one
        DLOOPA xp[l][j] = xp[myNp-1][j];
        //second, decrement the total number of particles
        myNp--;
        //since the current particle was replaced with a different one, process it again
        continue;
      }
      //current particle has not left the current tile, so move on to the next one
      l++;
    }
  }
  else {
    for (count = 0; count < *count_particles; count+=NDIM) {
      //place each newly acquired particle at the end of the list
      DLOOPA xp[myNp][j] = mpi_buf[count+j];
      //increment particle count
      myNp++;
    }
  }
  *count_particles = count;
  return(count);
}
#endif

//packs (ispack=1) or unpacks (ispack=0) the cells to be communicated along
//dimension (dim), either upper (isup=1) or lower (isup=0) boundary
//returns the number of items packed (count)
int pack_pflag(int ispack, int dim, int isup, int pflag[][N2M][N3M], int *mpi_buf )
{
  int istart,istop,jstart,jstop,kstart,kstop,i,j,k,m,count;
  //number of ghost cells to copy: need only layer of thickness one for pflag
  int n1g = (N1G>0), n2g = (N2G>0), n3g = (N3G>0);
  
  //if true, ensure it has value of unity for the below to work
  if(ispack) ispack=1;
  
  //do not do empty dimensions
  if(mpi_ntot[dim] == 1) return(0);
  
  //figure out the range of indices to transfer
  //x1: in x1-dim transfer only ghost cells immediately adjacent to active cells
  if(1==dim){
    jstart=0; jstop=N2-1;
    kstart=0; kstop=N3-1;
    if(isup) istart=N1-ispack*n1g; else istart=-n1g+ispack*n1g;
    istop=istart+n1g-1;
  }
  //x2: in x2-dim, additionally trasfer the ghost cells that have been communicated in x1-dim
  else if(2==dim){
    istart=-n1g; istop=N1+n1g-1;
    kstart=0; kstop=N3-1;
    if(isup) jstart=N2-ispack*n2g; else jstart=-n2g+ispack*n2g;
    jstop=jstart+n2g-1;
  }
  //x3, in x3-dim, additionally trasfer the ghost cells that have been communicated in x1-dim and x2-dim
  else if(3==dim){
    istart=-n1g; istop=N1+n1g-1;
    jstart=-n2g; jstop=N2+n2g-1;
    if(isup) kstart=N3-ispack*n3g; else kstart=-n3g+ispack*n3g;
    kstop=kstart+n3g-1;
  }
  
  //initialize the counter of the number of doubles (un)packed
  count = 0;
  
  //packing
  if(ispack){
    //OMPMARK
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) {
      mpi_buf[count++] = pflag[i][j][k];
    }
  }
  ///unpacking
  else {
    //OMPMARK
    ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) {
      pflag[i][j][k] = mpi_buf[count++];
    }
  }
  return(count);
}


int is_physical_bc( int dim, int isup )  //This is a helper function that answers: "Should I even bother applying a physical boundary condition in this direction?"

{
  /*Explanation:
    dim → which dimension: 1 → x1 (radial), 2 → x2 (theta), 3 → x3 (phi)
    isup → is it upper boundary? (1) or lower boundary? (0)
→ If mpi_ntot[dim] == 1 → There is only 1 MPI domain in this direction → No need for physical boundary condition (since the domain is periodic or fully internal).
Else → Return 1 → Apply physical boundary condition.
Example:
Situation	                                    is_physical_bc result
Single MPI domain in x1	                        return 0 (skip BC)
Multiple MPI domains in x1	                    return 1 (apply BC)

What is mpi_ntot[dim]?
    mpi_ntot[dim] stores how many MPI domains exist along direction dim.
    For example, if mpi_ntot[1] = 1 → There is only one domain in x1. → No physical boundary is needed because we must have MPI exchange with the neighbor domain.
    But if mpi_ntot[1] > 1, the domain touches a physical boundary → We should apply boundary conditions.*/
  
  //dimension is trivial => boundary is not physical
  if( 1 == mpi_ntot[dim] )
    return(0);
  
  //lower boundary is physical
  if( 0 == isup && 0 == mpi_coords[dim] && (1 == mpi_dims[dim] || 1 != mpi_periods[dim]) )
    return(1);
  
  //upper boundary is physical
  if( 1 == isup && mpi_dims[dim]-1 == mpi_coords[dim] && (1 == mpi_dims[dim] || 1 != mpi_periods[dim]) )
    return(1);
  
  //MPI boundary
  return(0);
}

//initiates send (isrecv = 0) or receive (isrecv = 1) operation in dimension dim (=0,1,2)
void bound_mpi_dim(int dim, int ispack, double prim[][N2M][N3M][NPR], int pflag[][N2M][N3M])
{
  /*
  Arguments:
        Variable	        Purpose
          dim	        Which dimension? (1 → x1, 2 → x2, 3 → x3)
          ispack	    1 → pack & send; 0 → receive & unpack
          prim	      Primitive variables array
          pflag	      Flags for zones (troubled cells etc.)
          */
#ifdef MPI
#define DOMPIPFLAG (0)
    
  int count, count_pflag, count_prims, count_particles, count_total, tagsend, tagrecv, isup; //Declare counters and flags: count_prims = # of primitives packed ; isup = 0 (lower neighbor) or 1 (upper neighbor)
  
  //packing and sending
  if(1 == ispack) { //We're in "packing & sending" mode
    for(isup=0;isup<=1;isup++) { //Loop twice: isup=0 → send/recv to lower neighbor; isup=1 → send/recv to upper neighbor
      //skip if on the physical boundary and the BC is not periodic
      if(is_physical_bc(dim,isup)) continue; // Skip if at Physical Boundary (no neighbor). Example: outer r boundary in a non-periodic run → skip send.

      //pack the data from prim[] into mpi send buffers; Pack prim values for ghost zones along dim and store in buffer → mpi_buf_send[..]. Returns number of elements packed → stored in count
      count_prims = count = pack_prim(ispack, dim, isup, prim, mpi_buf_send[dim][isup]); //

#if(DOMPIPFLAG) // Here it's disabled (DOMPIPFLAG=0), so skipped.
      //pack the data from pflag[] into mpi send buffers
      count_pflag = pack_pflag(ispack, dim, isup, pflag, mpi_buf_send_pflag[dim][isup]);
#endif
        
#if(DOPARTICLES) //Disabled
      //pack the data from xp[][] into mpi send buffers
      pack_particles(ispack, dim, isup, xp, mpi_buf_send[dim][isup]+count, &count_particles);
      count += count_particles;
#endif

      //prims
      /*MPI tags to differentiate lower vs upper communication.
          isup=0 → tagsend=0 (lower), tagrecv=2 (upper)
          isup=1 → tagsend=2 (upper), tagrecv=0 (lower)*/
      tagsend = 0+2*isup;
      tagrecv = 0+2*!isup;
      MPI_Isend(mpi_buf_send[dim][isup], //buffer that's being sent
                count,              //number of items sent
                MPI_DOUBLE,         //data type
                mpi_nbrs[dim][isup],//the rank of destination process
                tagsend,                //tag
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_send[dim][isup] //error
                );

      MPI_Irecv(mpi_buf_recv[dim][isup], //buffer that's being received
                count_prims+NPTOT*4,//max possible number of items to be received (>= those sent)
                MPI_DOUBLE,         //data type
                mpi_nbrs[dim][isup],//the rank of source process
                tagrecv,                //tag (should be same as in the send process)
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_recv[dim][isup] //error
                );
#if(DOMPIPFLAG) //Disabled
      //pflags
      tagsend = 1+2*isup;
      tagrecv = 1+2*!isup;
      MPI_Isend(mpi_buf_send_pflag[dim][isup], //buffer that's being sent
                count_pflag,              //number of items sent
                MPI_INT,            //data type
                mpi_nbrs[dim][isup],//the rank of destination process
                tagsend,                //tag
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_send_pflag[dim][isup] //error
                );
      
      MPI_Irecv(mpi_buf_recv_pflag[dim][isup], //buffer that's being received
                count_pflag,              //number of items received (same as those sent)
                MPI_INT,            //data type
                mpi_nbrs[dim][isup],//the rank of source process
                tagrecv,                //tag (should be same as in the send process)
                MPI_COMM_WORLD,     //communicator
                &mpi_reqs_recv_pflag[dim][isup] //error
                );
#endif
    }
  }
  //waiting, unpacking, and putting results back
  else {
    for(isup=0;isup<=1;isup++) {
      //skip if on the physical boundary and the BC is not periodic
      if(is_physical_bc(dim,isup)) continue;
      
      //wait for comminication to complete
      //prims
      MPI_Wait(&mpi_reqs_recv[dim][isup],&mpi_stat_recv[dim][isup]);
      MPI_Wait(&mpi_reqs_send[dim][isup],&mpi_stat_send[dim][isup]);
#if(DOMPIPFLAG)
      //pflags
      MPI_Wait(&mpi_reqs_recv_pflag[dim][isup],&mpi_stat_recv_pflag[dim][isup]);
      MPI_Wait(&mpi_reqs_send_pflag[dim][isup],&mpi_stat_send_pflag[dim][isup]);
#endif
      //unpack the data from mpi recv buffers into prim[]
      count_prims = pack_prim(ispack, dim, isup, prim, mpi_buf_recv[dim][isup]);

#if(DOPARTICLES)
      //need to find the size of the buffer, including particles
      MPI_Get_count(&mpi_stat_recv[dim][isup], MPI_DOUBLE, &count_total);
      count_particles = count_total - count_prims;
      //unpack the data from mpi receive buffers into xp[][]
      pack_particles(ispack, dim, isup, xp, mpi_buf_recv[dim][isup]+count_prims, &count_particles);
#endif
      
#if(DOMPIPFLAG)
      //pack the data from mpi send buffers into  pflag[]
      pack_pflag(ispack, dim, isup, pflag, mpi_buf_recv_pflag[dim][isup]);
#endif
    }
  }
#endif
}

//do not allow "sucking" on the boundary:
//type = 0: do not allow ucon[1] > 0
//type = 1: do not allow ucon[1] < 0
//if a disallowed value detected, reset vcon[1] to zero
void inflow_check(double *pr, int ii, int jj, int kk, int type )
{
        struct of_geom geom ;
        double ucon[NDIM] ;
        int j,k ;
        double alpha,beta1,gamma,vsq ;

        get_geometry(ii,jj,kk,CENT,&geom) ;
        ucon_calc(pr, &geom, ucon) ;

        if( ((ucon[1] > 0.) && (type==0)) || ((ucon[1] < 0.) && (type==1)) ) { 
                /* find gamma and remove it from primitives */
	  if( gamma_calc(pr,&geom,&gamma) ) { 
	    fflush(stderr);
	    fprintf(stderr,"\ninflow_check(): gamma failure, (%d,%d,%d) \n",
                    ii+mpi_startn[1], jj+mpi_startn[2], kk+mpi_startn[3]);
	    fflush(stderr);
	    fail(FAIL_GAMMA);
	  }
	    pr[U1] /= gamma ;
	    pr[U2] /= gamma ;
	    pr[U3] /= gamma ;
	    alpha = 1./sqrt(-geom.gcon[0][0]) ;
	    beta1 = geom.gcon[0][1]*alpha*alpha ;

	    /* reset radial velocity so radial 4-velocity
	     * is zero */
	    pr[U1] = beta1/alpha ;

	    /* now find new gamma and put it back in */
	    vsq = 0. ;
	    SLOOP vsq += geom.gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;
	    if( fabs(vsq) < 1.e-13 )  vsq = 1.e-13;
	    if( vsq >= 1. ) { 
	      vsq = 1. - 1./(GAMMAMAX*GAMMAMAX) ;
	    }
	    gamma = 1./sqrt(1. - vsq) ;
	    pr[U1] *= gamma ;
	    pr[U2] *= gamma ;
	    pr[U3] *= gamma ;

	    /* done */
	  }
	  else
	    return ;

}


/*
Boundary function Uses same indices for all variables because:
  // 1. The actual staggering is handled in flux calculations for fields
  // 2. Storage is cell-centered but interpretation is face-centered for fields

Storage vs. Interpretation:
  Magnetic fields are stored in the same prim array with cell-centered indexing
  The staggering is handled through interpretation during calculations, not in storage
  This is a common technique called "logical centering" where variables are stored cell-centered but interpreted as face-centered when needed
Staggering Implementation:
  Staggering is primarily implemented in the constrained transport algorithm
  The flux calculations handle the actual face-centered nature of magnetic fields
  This approach simplifies memory management while maintaining magnetic field evolution properties
Boundary Conditions:
  Use same indexing because the storage is unified
  The physical interpretation of face-centering happens during:
    Flux calculations
    Field evolution
    Constrained transport steps
  This design choice simplifies boundary condition implementation while maintaining consistency
This explains why the boundary conditions can use the same indexing - it's a design choice that separates:
  Storage (unified, cell-centered arrays)
  Physical interpretation (face-centered for magnetic fields)
  Numerical implementation (staggering in flux and evolution calculations)
*/