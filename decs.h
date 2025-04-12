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

#ifdef MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdarg.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/*************************************************************************
      COMPILE-TIME PARAMETERS : 
*************************************************************************/

/*
Grid Dimensions:

    N1, N2, N3: Number of zones (grid cells) in the X1, X2, and X3 directions.
    N1M, N2M, N3M: Total size of the arrays, including ghost zones (extra cells used for boundary conditions).
    NG: Number of ghost zones per dimension.

Loops:

    ZLOOP: Loops over all active zones (grid cells) in the 3D domain.
    ZSLOOP: Loops over a custom subrange of grid indices.

Coordinate Systems:

    GR: Whether general relativity (GR) is used.
    BL: Whether Boyer-Lindquist (BL) coordinates are used.
    
*/
#define NG (3) /*Number of ghost cells per boundary. This means 3 layers of ghost zones around the physical domain*/
#define WRITEGHOST (0) // Whether to write ghost zones to output files (0 = no, 1 = yes)
/* how many cells near the poles to stabilize, choose 0 for no stabilization */
#define POLEFIX (2) // Number of zones near the poles (in θ-direction) where special stabilization (like smoothing) is applied: 0 means no stabilization; 2 means stabilize 2 cells near the pole

//whether to write out ener files. since so far no use of them, set to zero //in diag.c file
#define DOENER (0)

//which problem
#define MONOPOLE_PROBLEM_1D 1
#define MONOPOLE_PROBLEM_2D 2
#define BZ_MONOPOLE_2D 3
#define TORUS_PROBLEM 4
#define BONDI_PROBLEM_1D 5
#define BONDI_PROBLEM_2D 6
#define SHOCK_TEST 7
#define OT_TEST    8
#define SHOCK_ENTROPY_TEST 9
#define LIN_MODES         10
#define SNDWAVE_TEST 11
#define ENTWAVE_TEST 12
#define VTRAN_TEST 13
#define TURB 14
#define GRAD 15
#define STATCOND 16
#define ATM_TEST 17
#define BONDI_CON 18
#define STATCOND1D 19
#define HUBBLE 20
#define TORUS_GRB 21

#define WHICHPROBLEM TORUS_PROBLEM
#define DIS_DETECT 0
#define BETAHEAT 1
#define eCOND    0
#define eHEAT    1
#define DOKTOT   1
#define DOFLR    0
#define DONUCLEAR 0
#define DOCYLINDRIFYCOORDS 1
#define EVOLVEVPOT 1
#define DOPARTICLES 0
#define DOQDOTAVG 1

/** here are the few things that we change frequently **/

#if WHICHPROBLEM == MONOPOLE_PROBLEM_1D
#define N1       (2*386)      /* number of physical zones in X1-direction */
#define N2       (1)          /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == MONOPOLE_PROBLEM_2D
#define N1       (2*100) //386)      /* number of physical zones in X1-direction */
#define N2       (100)   //256)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == BZ_MONOPOLE_2D
#define N1       (128)      /* number of physical zones in X1-direction */
#define N2       (128)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == TORUS_PROBLEM
#define N1       (16)        /* number of physical zones in X1-direction */  //change back to 256x256
#define N2       (16)        /* number of physical zones in X2-direction */
#define N3       (20)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (3)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == TORUS_GRB
#define N1       (64)        /* number of physical zones in X1-direction */  //change back to 256x256
#define N2       (64)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (2)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
//#undef DONUCLEAR
//#define DONUCLEAR (0)         /* whether or not to keep track of nuclear and weak interactions in a GRB disk */
#elif WHICHPROBLEM == BONDI_PROBLEM_1D
#define N1       (256)        /* number of physical zones in X1-direction */
#define N2       (1)          /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == BONDI_PROBLEM_2D
#define N1       (100)        /* number of physical zones in X1-direction */
#define N2       (100)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (0)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == SHOCK_TEST
#define N1/*shock*/(20)         /* number of physical zones in X1-direction */
#define N2/*shock*/(1)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (1)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == SNDWAVE_TEST
#define N1/*sndwave*/(2000)         /* number of physical zones in X1-direction */
#define N2/*sndwave*/(1)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == ENTWAVE_TEST
#define N1/*ent*/(32)         /* number of physical zones in X1-direction */
#define N2       (1)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == TURB
#define N1/*turb*/(128)         /* number of physical zones in X1-direction */
#define N2/*turb*/(128)         /* number of physical zones in X1-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == OT_TEST
#define N1       (100)        /* number of physical zones in X1-direction */
#define N2       (100)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == SHOCK_ENTROPY_TEST
#define N1       (500)        /* number of physical zones in X1-direction */
#define N2       (1)        /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == LIN_MODES
#define N1/*mode*/(512)         /* number of physical zones in X1-direction */
#define N2/*mode*/(512)         /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == VTRAN_TEST
#define N1/*vtran*/(512)         /* number of physical zones in X1-direction */
#define N2/*vtran*/(1)         /* number of physical zones in X2-direction */
#define N3       (1)          /* number of physical zones in X3-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == GRAD
#define N1/*grad*/(500)         /* number of physical zones in X1-direction */
#define N2/*grad*/(1)         /* number of physical zones in X2-direction */
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == STATCOND
#define N1/*grad*/(250)         /* number of physical zones in X1-direction */
#define N2/*grad*/(250)         /* number of physical zones in X2-direction */
#define N3       (1)
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == ATM_TEST
#define N1/*atm*/(2000)         /* number of physical zones in X1-direction */
#define N2/*atm*/(1)         /* number of physical zones in X2-direction */
#define N3       (1)
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define DIRICHLET (1)
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#undef POLEFIX
#define POLEFIX (0)           /* disable POLEFIX for this test */
#elif WHICHPROBLEM == BONDI_CON
#define N1/*bondi*/(2048)         /* number of physical zones in X1-direction */
#define N2       (1)         /* number of physical zones in X2-direction */
#define N3       (1)
#define GR       (1)          /* whether or not to use GR */
#define BL       (1)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define DIRICHLET (1)
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == STATCOND1D
#define N1       (100)         /* number of physical zones in X1-direction */
#define N2       (1)         /* number of physical zones in X2-direction */
#define N3       (1)
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (1)          /* whether or not to use periodic boundary conditions */
#define OUTFLOW  (0)          /* whether or not to use outflow boundary conditions in all directions */
#elif WHICHPROBLEM == HUBBLE
#define N1/*hubble*/(100)         /* number of physical zones in X1-direction */
#define N2       (1)         /* number of physical zones in X2-direction */
#define N3       (1)
#define GR       (0)          /* whether or not to use GR */
#define BL       (0)          /* whether or not to use BL coords */
#define INFLOW   (1)          /* whether or not to allow inflow at boundaries */
#define PERIODIC (0)          /* whether or not to use periodic boundary conditions */
#define DIRICHLET (0)
#define OUTFLOW  (1)          /* whether or not to use outflow boundary conditions in all directions */
#endif


//only allocate memory for ghost cells for non-trivial dimensions
#define N1M ((N1>1)?(N1+2*NG):(1)) /*Total grid size including ghost cells. For example, N1M = N1 + 2*NG if N1 > 1, so like Total r cells = 16 + 2×3 = 22 */
#define N2M ((N2>1)?(N2+2*NG):(1)) /*Total grid size including ghost cells. For example, N1M = N1 + 2*NG if N1 > 1, so like Total θ cells = 16 + 2×3 = 22 */
#define N3M ((N3>1)?(N3+2*NG):(1)) /*Total grid size including ghost cells. For example, N1M = N1 + 2*NG if N1 > 1, so like Total φ cells = 20 + 2×3 = 26 */
/* Example (if N1 = 16 and NG = 3):

    Physical r-indices = 3 to 18

    Ghost cells:

        Inner = 0,1,2

        Outer = 19,20,21*/

#define N1G ((N1>1)?(NG):(0))
#define N2G ((N2>1)?(NG):(0))
#define N3G ((N3>1)?(NG):(0)) /*Ghost cell count per side (NG = 3 here). Used for indexing ghost zone*/

/*Flags indicating whether the domain spans >1 cell in a given direction
Used for logic like: if D1=0, no derivatives in x1-direction 
This just checks if the dimension is "active" (more than 1 cell).
    D1 = 1 → r is active
    D2 = 1 → θ is active
    D3 = 1 → φ is active 
Some parts of the code may run in lower dimension mode (like 1D tests).Example use:
if(D2) {
  // Do theta-related computations
}
*/
#define D1 (N1>1)
#define D2 (N2>1)
#define D3 (N3>1)

#if(DONUCLEAR)
#define NPR        (13)
#elif(eCOND || eHEAT)
#if(DOFLR)
#define NPR        (19)        /* number of primitive variables */
#else
#define NPR        (18)
#endif
#elif(DOKTOT)
#define NPR        (9)
#else
#define NPR        (8)
#endif
#define NDIM       (4)        /* number of total dimensions.  Never changes */
#define NPG        (5)        /* number of positions on grid for grid functions -> This is usually about locations like cell-center, face-center, corner, etc */
#define COMPDIM    (2)        /* number of non-trivial spatial dimensions used in computation */

#define NIMG       (5)        /* Number of types of diagnostics to save into fdump */






/* whether or not to use Font's  adiabatic/isothermal prim. var. inversion method: */
#define DO_FONT_FIX (1)

/* whether or not to rescale primitive variables before interpolating them for flux/BC's */
#define RESCALE     (0)

/* Reconstruction method */
#define LIN         (0)
#define PARA        (1)
#define RECONSTRUCT    (PARA)


/** FIXUP PARAMETERS, magnitudes of rho and u, respectively, in the floor : **/
#if DONUCLEAR
#define RHOMIN	(1.e-5)
#define UUMIN	(1.e-7)
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-22)
#else
#define RHOMIN	(1.e-6)
#define UUMIN	(1.e-8)
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#endif
#define POWRHO (2)

#define FLOORFACTOR (1.)
#define BSQORHOMAX (50.*FLOORFACTOR)
#define BSQOUMAX (2500.*FLOORFACTOR)
#define UORHOMAX (50.*FLOORFACTOR)

//add mass in the ZAMO frame (=1) instead of fluid frame (=0)
#define WHICH_FLOOR DRIFT_FRAME_FLOOR
#define DRIFT_FRAME_FLOOR (0)
#define ZAMO_FRAME_FLOOR (1)

/* A numerical convenience to represent a small non-zero quantity compared to unity:*/
#define SMALL	(1.e-20)

/* Max. value of gamma, the lorentz factor */
#define GAMMAMAX (50.)

/* maximum fractional increase in timestep per timestep */
#define SAFE	(1.3)


#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
#define SINGSMALL (1.E-20)


/* I/O format strings used herein : */
#define FMT_DBL_OUT "%32.22e"
#define FMT_INT_OUT "%10d"

#define MAXLEN (100)


/*************************************************************************
    MNEMONICS SECTION 
*************************************************************************/
/* boundary condition mnemonics */
//the below turns out not to be used
//#define OUTFLOW	(0)
//#define SYMM	(1)
//#define ASYMM	(2)
//#define FIXED	(3)

/* mnemonics for primitive vars; conserved vars */
#define RHO	(0)	
#define UU	(1)
#define U1	(2)
#define U2	(3)
#define U3	(4)
#define B1	(5)
#define B2	(6)
#define B3	(7)
#define KTOT    (8)

#define KEL4     (9)
#define KEL4A    (10)
#define KEL4B    (11)
#define KEL4C    (12)
#define KEL4D    (13)
#define KEL4E    (14)

#define KEL5    (15)
#define KELDIS  (16)
#define PHI     (17)
#define FLR     (18)

#if( DONUCLEAR )

#if(!DOKTOT || eCOND || eHEAT)
#error Can only use DONUCLEAR with DOKTOT but without eCOND and eHEAT
#endif

#define RHONP     (9)
#define RHOALPHA  (10)
#define RHOFLOOR  (11)
#define YE        (12)
#endif

//#define PHI2 (13)
//#define PHI3 (14)
//#define KEL42 (15)
//#define KEL43 (16)

/* mnemonics for dimensional indices */
// Index identifiers for coordinate directions
// TT=Time, RR=Radial, TH=Theta, PH=Phi
#define TT	(0)     
#define RR	(1)
#define TH	(2)
#define PH	(3)

/* mnemonics for centering of grid functions */
// Index identifiers for locations of variables
/*This tells where a variable is defined within a cell or between cells.
Name	Meaning	            Where variable is located
FACE1	On x1 (r) face	     r-interfaces (like B1)
FACE2	On x2 (θ) face	     θ-interfaces (like B2)
FACE3	On x3 (φ) face	     φ-interfaces (like B3)
CORN	At the corner of a cell	  (not common)
CENT	Cell center	         Usual for density, pressure, etc
EDGE1	Edge in x1	          Uncommon
EDGE2	Edge in x2	          Uncommon
EDGE3	Edge in x3	          Uncommon
Edges (like EDGE1, EDGE2, EDGE3) are actually defined too in your code: 
But — they are not part of the default NPG=5 locations — unless explicitly needed for something special like vector potential, reconstructions, or CT (constrained transport)*/
#define FACE1	(0)	// face between i and i+1 in X1 (r)
#define FACE2	(1) // face between j and j+1 in X2 (θ)
#define FACE3	(2) // face between k and k+1 in X3 (φ)
#define CORN	(3) // grid corner (i, j, k) -> Used for geometry quantities (interpolations etc)
#define CENT	(4) // cell center (i+½, j+½, k+½) -> Cell-centered quantities (primitive variables, conserved variables)
#define EDGE1   (5)
#define EDGE2   (6)
#define EDGE3   (7)

/* mnemonics for slope limiter */
#define MC	(0)
#define VANL	(1)
#define MINM	(2)

/* mnemonics for diagnostic calls */
#define INIT_OUT	(0)
#define DUMP_OUT	(1)
#define RDUMP_OUT	(2)
#define IMAGE_OUT	(3)
#define LOG_OUT		(4)
#define FINAL_OUT	(5)
#define DIVB_OUT        (6)
#define PDUMP_OUT       (7)

/* Directional Mnemonics */
// -------------> r
// |         3    
// |        1-0   
// |         2    
// v            
// theta      
#define X1UP    (0)
#define X1DN    (1)
#define X2UP    (2)
#define X2DN    (3)


/* failure modes */
#define FAIL_UTOPRIM        (1)
#define FAIL_VCHAR_DISCR    (2)
#define FAIL_COEFF_NEG	    (3)
#define FAIL_COEFF_SUP	    (4)
#define FAIL_GAMMA          (5)
#define FAIL_METRIC         (6)


/* For rescale() operations: */
#define FORWARD 1
#define REVERSE 2



/*************************************************************************
 MACROS
 *************************************************************************/
/* loop over all active zones */
#define ZLOOP for(i=0;i<N1;i++)for(j=0;j<N2;j++)for(k=0;k<N3;k++)

/* loop over all active zones */
#define IMAGELOOP for(k=0;k<N3;k++)for(j=0;j<N2;j++)for(i=0;i<N1;i++)

/* specialty loop */
extern int istart,istop,jstart,jstop,kstart,kstop ;
#define ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) \
for(i=istart;i<=istop;i++)\
for(j=jstart;j<=jstop;j++)\
for(k=kstart;k<=kstop;k++)

/* loop over Primitive variables */
#define PLOOP  for(m=0;m<NPR;m++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++) for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)
/* loop over all electron variables */
#define eLOOP for (m = KEL4; m<=KELDIS; m++)


extern double fval1,fval2;
#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))
#define MY_SIGN(fval) ( ((fval) <0.) ? -1. : 1. )

#define NMAX     MY_MAX(MY_MAX(N1,N2),N3) /* this sizes 1D slices */

#define delta(i,j) (((i) == (j)) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

//PARTICLES
#if(DOPARTICLES)
/* particles */
#define NPTOT (10000)
extern double xp[NPTOT][NDIM];
extern int myNp;
#else
#define NPTOT (0)
#endif


////////////////////////////////////////////////////////
//
//  MPI section
//
////////////////////////////////////////////////////////

typedef float dumptype;
typedef float sdumptype;
typedef float gdumptype;
typedef double gdump2type;
typedef double rdumptype;
typedef long long fdumptype;

#ifdef MPI
///////////////////////////////////////////////
// how to write dumps and gdumps in parallel
#define DO_PARALLEL_WRITE (1)

extern void *mpi_file_buffer;
extern dumptype *dump_buffer;
extern sdumptype *sdump_buffer;
extern gdumptype *gdump_buffer;
extern gdump2type *gdump2_buffer;
extern rdumptype *rdump_buffer;
extern fdumptype *fdump_buffer;

#define MPI_DUMP_TYPE MPI_FLOAT
#define MPI_SDUMP_TYPE MPI_FLOAT
#define MPI_GDUMP_TYPE MPI_FLOAT
#define MPI_GDUMP2_TYPE MPI_DOUBLE
#define MPI_RDUMP_TYPE MPI_DOUBLE
#define MPI_FDUMP_TYPE MPI_LONG_LONG_INT

#define DUMP_FILE (0)
#define SDUMP_FILE (1)
#define GDUMP_FILE (2)
#define GDUMP2_FILE (3)
#define RDUMP_FILE (4)
#define FDUMP_FILE (5)

extern MPI_Datatype gdump_file_type, gdump_cell_type;
extern MPI_Datatype gdump2_file_type, gdump2_cell_type;
extern MPI_Datatype dump_file_type, dump_cell_type;
extern MPI_Datatype sdump_file_type, sdump_cell_type;
extern MPI_Datatype rdump_file_type, rdump_cell_type;
extern MPI_Datatype fdump_file_type, fdump_cell_type;
//
/////////////////////////////////////////////////
extern int mpi_numtasks;
extern MPI_Comm mpi_cartcomm;
extern int mpi_reorder;
extern int mpi_nbrs[NDIM][2];
//primitives
extern MPI_Request mpi_reqs_send[NDIM][2];
extern MPI_Request mpi_reqs_recv[NDIM][2];
extern MPI_Status mpi_stat_send[NDIM][2];
extern MPI_Status mpi_stat_recv[NDIM][2];
extern double mpi_buf_send[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR+NPTOT];
extern double mpi_buf_recv[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR+NPTOT];
//pflag
extern MPI_Request mpi_reqs_send_pflag[NDIM][2];
extern MPI_Request mpi_reqs_recv_pflag[NDIM][2];
extern MPI_Status mpi_stat_send_pflag[NDIM][2];
extern MPI_Status mpi_stat_recv_pflag[NDIM][2];
extern int mpi_buf_send_pflag[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR];
extern int mpi_buf_recv_pflag[NDIM][2][(NMAX+2*NG)*(NMAX+2*NG)*NG*NPR];
#else
#define DO_PARALLEL_WRITE (0)
#endif
extern int mpi_periods[NDIM];
extern int mpi_rank;
extern int mpi_coords[NDIM];
extern int mpi_dims[NDIM];
extern int mpi_ntot[NDIM];
extern int mpi_ntile[NDIM];
extern int mpi_startn[NDIM];
extern int i_am_the_master;
extern double mpi_startx[NDIM], mpi_stopx[NDIM];
#define MPI_NDIM (NDIM-1)
#define MASTER (0)



/*************************************************************************
    GLOBAL ARRAY SECTION 
*************************************************************************/
//p[N1M][N2M][N3M][NPR]: Primitive variables (density, velocities, B-field, etc.) with ghost zone
/*    NPR → Number of primitive variables
Usually: | Index | Variable | |-------|-----------| | 0 | Density (ρ) | | 1 | Internal Energy (u) | | 2 | u^r velocity | | 3 | u^θ velocity | | 4 | u^φ velocity | | 5 | B^r | | 6 | B^θ | | 7 | B^φ |
(Exact order depends on the code setup)
In our code, NPR = 18, ρ (rest-mass density),u (internal energy),v^i (3-velocity),B^i (magnetic field components),extra variables if you include conduction, heating, nuclear reactions, etc */
//extern means that the actual memory allocation and definition of the a_p array occur in another part of your C program (in a different source file). This declaration simply tells the compiler that a variable with this name and type exists elsewhere.
extern double  a_p[N1M][N2M][N3M][NPR] ;	/* space for primitive vars */
extern double  a_dq[N1M][N2M][N3M][NPR] ;  /* slopes */
extern double  a_dsource[N1M][N2M][N3M] ;  /* sources */
extern double  a_duscon[N1M][N2M][N3M] ;  /* sources */
extern double  a_sour[N1M][N2M][N3M] ;  /* sources */
extern double  a_F1[N1M][N2M][N3M][NPR] ;	/* fluxes */
extern double  a_F2[N1M][N2M][N3M][NPR] ;	/* fluxes */
extern double  a_F3[N1M][N2M][N3M][NPR] ;	/* fluxes */
extern double  a_uelvar[N1M][N2M][N3M]; /* electron internal energy */
extern double  a_thetaold[N1M][N2M][N3M]; /* temperature at previous time step for guess */
extern double  a_qdot[N1M][N2M][N3M]; /* source term diagnostics */
#if(DOQDOTAVG)
extern double  a_qdotavg[N1M][N2M][N3M]; /* source term diagnostics */
extern double  a_feqdotavg[N1M][N2M][N3M]; /* source term diagnostics */
#endif
extern double  a_qisosq[N1M][N2M][N3M]; /* source term diagnostics */
extern double  a_qisodotb[N1M][N2M][N3M]; /* source term diagnostics */


extern double a_uelarray[N1M][N2M][N3M][KELDIS-KEL4 +1]; /* space for electron internal energy (used for failure points) */
extern double a_du_flr[N1M][N2M][N3M]; /* energy added by floor */





extern double complex a_delvx[N1M][N2M];
extern double complex a_delvy[N1M][N2M];
extern double  a_ph[N1M][N2M][N3M][NPR] ;	/* half-step primitives */
extern double  a_psave[N1M][N2M][N3M][NPR] ;
extern double  a_pbound[N1M][N2M][N3M][NPR] ;


extern int     a_pflag[N1M][N2M][N3M];	/* identifies failure points */

/* for debug */
//extern double fimage[NIMG][N1*N2*N3];
extern long long    failimage[N1][N2][N3][NIMG];

/* grid functions 
NPG = 5 → Usually (not sure):
    0 = FACE1
    1 = FACE2
    2 = FACE3
    3 = CORNER
    4 = CENTER 
These five locations are the only ones actually used for centering grid metric functions like a_gcon, a_conn, etc. 
So NPG = 5 refers specifically to grid position types relevant to geometry, not all possible edge locations 
So yes, the name NPG = 5 is code-specific: it reflects only the five centering types that this particular simulation framework uses for metric/grid-related quantities — not an exhaustive list of all centerings possible.*/
extern double a_conn[N1M][N2M][N3M][NDIM][NDIM][NDIM] ; //stores Christoffel symbols Γ^a_bc at every grid location NDIM = 4 (t, r, θ, φ) For each cell (i,j,k), it stores 4×4×4 Christoffel symbols needed for GR calculations
extern double a_gcon[N1M][N2M][N3M][NPG][NDIM][NDIM] ; //This stores inverse metric g^ab at each of NPG different positions in a grid cell.
extern double a_gcov[N1M][N2M][N3M][NPG][NDIM][NDIM] ;
extern double a_gdet[N1M][N2M][N3M][NPG] ;
extern double a_phys_coords[NDIM][N1M][N2M][N3M];
extern double emf[NDIM][N1+D1][N2+D2][N3+D3] ; //OPTMARK: could reduce NDIM to 1 in 2D and eliminate completely in 1D
#if(EVOLVEVPOT)
extern double vpot[NDIM][N1+D1][N2+D2][N3+D3] ; //OPTMARK: could reduce NDIM to 1 in 2D and eliminate completely in 1D
#endif

extern double (*   p)[N2M][N3M][NPR] ;
extern double (*   psave)[N2M][N3M][NPR] ;
extern double (*   pbound)[N2M][N3M][NPR];
extern double (*  dq)[N2M][N3M][NPR] ;
extern double (*  dsource)[N2M][N3M] ;
extern double (*  duscon)[N2M][N3M] ;
extern double (*  sour)[N2M][N3M] ;
extern double (*  F1)[N2M][N3M][NPR] ;
extern double (*  F2)[N2M][N3M][NPR] ;
extern double (*  F3)[N2M][N3M][NPR] ;
extern double (*  ph)[N2M][N3M][NPR] ;
extern int    (*  pflag)[N2M][N3M];
extern double (* conn)[N2M][N3M][NDIM][NDIM][NDIM] ;
extern double (* gcon)[N2M][N3M][NPG][NDIM][NDIM] ;
extern double (* gcov)[N2M][N3M][NPG][NDIM][NDIM] ;
extern double (* gdet)[N2M][N3M][NPG] ;
extern double (* phys_coords)[N1M][N2M][N3M];
extern double (* uelvar)[N2M][N3M];
extern double (* thetaold)[N2M][N3M];
extern double (* qdot)[N2M][N3M] ;
#if(DOQDOTAVG)
extern double (* qdotavg)[N2M][N3M] ; //pointer to a 2 dimesional array of double element; The parentheses around * qdotavg are crucial. They indicate that qdotavg itself is a pointer. If the parentheses were absent, double * qdotavg[N2M][N3M] would be interpreted as a 2D array of pointers to doubles, which is a different thing entirely.
extern double (* feqdotavg)[N2M][N3M] ;
#endif
extern double (* qisosq)[N2M][N3M] ;
extern double (* qisodotb)[N2M][N3M] ;
extern double complex (*delvx)[N2M];
extern double complex (*delvy)[N2M];

extern double (*   uelarray)[N2M][N3M][KELDIS-KEL4 +1] ;
extern double (*   du_flr)[N2M][N3M] ;




#if(DO_FONT_FIX)
extern double Katm[N1];
#endif

/*
Key Difference from the Previous Example:
In the previous example (extern double a_qdotavg[N1M][N2M][N3M];), a_qdotavg was the array itself. Here, (extern double (* qdotavg)[N2M][N3M]) qdotavg is a pointer that can point to such an array (specifically, a 2D array of doubles with dimensions N2M x N3M).
To use qdotavg, you would typically need to assign it the memory address of an actual 2D array of the correct type and dimensions that is defined elsewhere in your program.
*/
/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
extern double a ;
extern double gam ;
extern double game;
extern double game4;
extern double game5;
extern double mrat ;
extern double qosc ;
extern double kelmin;
extern double fel0;
extern double felfloor;
extern double rmax;
extern double rhomax;

/* numerical parameters */
extern double Rin,Rout,hslope,R0,fractheta,fracphi ;
extern double x1br, cpow2,npow2, rbr;
extern double cour ;
extern double dV,dx[NDIM],startx[NDIM],lenx[NDIM] ;
extern double dt ;
extern double t,tf ;
extern double x1curr,x2curr ;
extern int nstep ;
extern int threadid; //openmp id
extern int nthreads; //openmp number of threads

/* output parameters */
extern double DTd ;
extern double DTl ;
extern double DTi ;
extern double DTr ;
extern int    DTr01 ;
extern int    dump_cnt ;
extern int    image_cnt ;
extern int    rdump_cnt ;
extern int    rdump01_cnt ;
extern int    nstroke ;

/* global flags */
extern int failed ;
extern int lim ;
extern double defcon ;

/* diagnostics */
extern double mdot ;
extern double edot ;
extern double ldot ;

/* set global variables that indicate current local metric, etc. */
extern int icurr,jcurr,kcurr,pcurr ;
struct of_geom {
	double gcon[NDIM][NDIM] ;
	double gcov[NDIM][NDIM] ;
	double g ;
} ;

struct of_state {
	double ucon[NDIM] ;
	double ucov[NDIM] ;
	double bcon[NDIM] ;
	double bcov[NDIM] ;
} ;


/*************************************************************************
    FUNCTION DECLARATIONS 
*************************************************************************/
double bl_gdet_func(double r, double th, double ph) ;
double bsq_calc(double *pr, struct of_geom *geom) ;
int    FFT2D(complex double c[][N2M],int nx,int ny,int dir);
int    FFT(int dir,int m,double *x,double *y);
int    gamma_calc(double *pr, struct of_geom *geom, double *gamma) ;
void ut_calc_3vel(double *vcon, struct of_geom *geom, double *ut);
double gdet_func(double lgcov[][NDIM]) ;
double mink(int j, int k) ;
int Powerof2(int n,int *m,int *twopm);
double ranc(int seed) ;
double slope_lim(double y1, double y2, double y3) ;


int restart_init(void) ;

void area_map(int i, int j, int k, double prim[][N2M][N3M][NPR]) ;
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon) ;
void blgset(int i, int j, int k, struct of_geom *geom);
void bl_coord(double *X, double *r, double *th, double *ph) ;
void bl_coord_vec(double *X, double *V) ;
void bl_gcon_func(double r, double th, double ph, double gcov[][NDIM]) ;
void bl_gcov_func(double r, double th, double ph, double gcov[][NDIM]) ;
void bound_prim(double pr[][N2M][N3M][NPR]) ;

void dxdxp_func(double *X, double dxdxp[][NDIM]);
void conn_func(double *X, struct of_geom *geom, double lconn[][NDIM][NDIM]) ;
void coord(int i, int j, int k, int loc, double *X) ;
void diag(int call_code) ;
void diag_flux(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR]) ;
size_t dump(int dump_cnt, int issmall, int is_dry_run) ;
size_t gdump(int is_dry_run) ;
size_t gdump2(int is_dry_run) ;
void fdump(int dumpno);
void pdump(int dumpno, int issmall);


void fail(int fail_type) ;
void fixup(double (* pv)[N2M][N3M][NPR]) ;
void fixup1zone( int i, int j, int k, double prim[NPR] ) ;
void fixup_utoprim( double (*pv)[N2M][N3M][NPR] )  ;
void set_Katm(void);
int  get_G_ATM( double *g_tmp );
void fix_flux(double F1[][N2M][N3M][NPR], double F2[][N2M][N3M][NPR], double F3[][N2M][N3M][NPR]) ;
void flux_ct(double F1[][N2M][N3M][NPR],double F2[][N2M][N3M][NPR],double F3[][N2M][N3M][NPR]) ;
void gaussj(double **tmp, int n, double **b, int m) ;
void gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]) ;
void gcov_func(double *X, double lgcov[][NDIM]) ;
void get_geometry(int i, int j, int k, int loc, struct of_geom *geom) ;
void get_phys_coord_vec(int ii, int jj, int kk, double *V);
void get_phys_coord_r(int ii, int jj, int kk, double *r);
void get_phys_coord(int ii, int jj, int kk, double *r, double *theta, double *phi);

void get_state(double *pr, struct of_geom *geom, struct of_state *q) ;
void image_all( int image_count ) ;
void init(void) ;
void lower(double *a, struct of_geom *geom, double *b) ;
void ludcmp(double **a, int n, int *indx, double *d) ;
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd)  ;
void misc_source(double *ph, double *phxp1, double *phxm1, double *phyp1, double *phym1,int ii, int jj, int kk, struct of_geom *geom,
                 struct of_state *q, double *dU, double Dt);
void primtoflux(double *pa, struct of_state *q, int dir, struct of_geom *geom,
			double *fl) ;
void primtoU(double *p, struct of_state *q, struct of_geom *geom, double *U);
void weno(double x1, double x2, double x3, double x4,
          double x5, double *lout, double *rout) ;
void para(double x1, double x2, double x3, double x4,
          double x5, double *lout, double *rout) ;
void raise(double *v1, struct of_geom *geom, double *v2) ;
void rescale(double *pr, int which, int dir, int ii, int jj, int kk, int face,
			struct of_geom *geom) ;
void restart_write(int dump_cnt) ;
int restart_read(int dump_cnt);
void pdump_read(int dumpno, int issmall);
void set_arrays(void) ;
void set_grid(void) ;
void set_points(void) ;
void step_ch(double *ndt1, double *ndt2, double *ndt3) ;
void source(double *pa, struct of_geom *geom, int ii, int jj, int kk, double *Ua,double Dt) ;
void timestep(void) ;
void u_to_v(double *pr, int i, int j, int k) ;
void ucon_calc(double *pr, struct of_geom *geom, double *ucon) ;
void ucon_to_utcon(double *ucon,struct of_geom *geom, double *utcon);
void usrfun(double *pr,int n,double *beta,double **alpha) ;
void Utoprim(double *Ua, struct of_geom *geom, double *pa) ;
int Utoprim_2d(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR]);
int Utoprim_1dvsq2fix1(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR], double K );
int Utoprim_1dfix1(double U[NPR], double gcov[NDIM][NDIM], double gcon[NDIM][NDIM], 
               double gdet, double prim[NPR], double K );

void vchar(double *pr, struct of_state *q, struct of_geom *geom,
		int dir, double *cmax, double *cmin) ;

int invert_matrix( double A[][NDIM], double Ainv[][NDIM] );
int LU_decompose( double A[][NDIM], int permute[] );
void LU_substitution( double A[][NDIM], double B[], int permute[] );

void reconstruct_lr_lin(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]); 
void reconstruct_lr_par(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]);
void reconstruct_lr_weno(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]);

int mpi_init(int argc,char *argv[]);

int is_physical_bc( int dim, int isup );

void initialize_parallel_write(int stage);
void de_initialize_parallel_write();
void append_rank(char *name);

void parallel_readwrite(char *file_name, void *dump_buffer,
                        int type_of_file, int is_write, long long offset);
size_t write_to_dump( int is_dry_run, FILE *fp, dumptype *buf, dumptype val );
size_t write_to_gdump( int is_dry_run, FILE *fp, gdumptype *buf, gdumptype val );
size_t write_to_gdump2( int is_dry_run, FILE *fp, gdump2type *buf, gdump2type val );
size_t write_to_rdump( int is_dry_run, FILE *fp, rdumptype *buf, rdumptype val );

void getmax_densities(double (*prim)[N2M][N3M][NPR], double *rhomax, double *umax);
double get_maxprimvalrpow(double (*prim)[N2M][N3M][NPR], double rpow, int m );
int normalize_field_local_nodivb(double targbeta, double rhomax, double amax,
                                 double (*prim)[N2M][N3M][NPR],
                                 double (*A)[N2+D2][N3+D3], int dir);
double compute_rat(double (*prim)[N2M][N3M][NPR], double (*A)[N2+D2][N3+D3],
                   double rhomax, double amax, double targbeta, int loc, int i, int j, int k);


double compute_profile( double (*prim)[N2M][N3M][NPR], double amax, double aphipow, int loc, int i, int j, int k );
int compute_vpot_from_gdetB1( double (*prim)[N2M][N3M][NPR],  double (*A)[N2+D2][N3+D3] );
void get_rho_u_floor( double r, double th, double phi, double *rho_floor, double *u_floor );
void advance_particles(double ph[N1M][N2M][N3M][NPR], double Dt);
void init_particles();



extern double tdump,trdump,timage,tlog,tavgstart ;

////////////////////////////////
//SJETCOORDS
////////////////////////////////
extern double global_fracdisk;
extern double global_fracjet;
extern double global_disknu1;
extern double global_disknu2;
extern double global_jetnu1;
extern double global_jetnu2;
extern double global_rsjet;
extern double global_r0grid;
extern double global_r0jet;
extern double global_rjetend;
extern double global_r0disk;
extern double global_rdiskend;
extern double global_x10;
extern double global_x20;

//TORUS
extern double global_kappa;
extern double aphipow;

//DONUCLEAR
#if(DONUCLEAR)
double compute_temperature(double rho, double p, double Ye);
double compute_degeneracy(double rho, double T, double Ye);
#endif

/*
Item	        Meaning	                    Example/Value
NG	            Ghost cells	                    3
N1M, N2M, N3M	Total cells incl. ghost	    22, 22, 26
N1G, N2G, N3G	Ghost offset	                3 each
D1, D2, D3	    Is dimension active?	    (True for N1=16 >1)
NPG	            #positions in grid	        5 → FACE1, FACE2, FACE3, CORN, CENT
NPR	            #Primitive variables	     18 (in your case)

After analyzing the provided files, here's the comprehensive breakdown of variable locations in the staggered grid:
Cell Centered ( CENT):

Primitive variables (prim[NPR]): RHO (density)

 UU (internal energy)
UTCON1, UTCON2, UTCON3 (velocities)
Pressure (computed from EOS)
Conserved variables (U[NPR])
State variables (in struct of_state):ucon[NDIM], ucov[NDIM] (4-velocities), bcon[NDIM], bcov[NDIM] (magnetic 4-vectors)
Metric quantities (in struct of_geom):gcon[NDIM][NDIM], gcov[NDIM][NDIM] (metric tensors)
g (metric determinant)

Face Centered:
Magnetic field components:B1/BCON1 on FACE1 (r-faces)
B2/BCON2 on FACE2 (θ-faces)
B3/BCON3 on FACE3 (φ-faces)

Fluxes (F1, F2, F3 arrays)
gdet (metric determinant) at faces for flux calculations

Corner (CORN):
Vector potential components (A[N1+D1][N2+D2][N3+D3])
Grid coordinates at corners
Connection coefficients (conn[NDIM][NDIM][NDIM])

Edge Centered:
Electric fields/ or EMFs (EDGE1, EDGE2, EDGE3)
EMF (electromotive force) components for constrained transport
Area elements for flux calculations
In ideal MHD, the electric field is not an independent variable because of the ideal MHD condition:
E = -v × B (in special relativity), or more generally in general relativity:, E^μ = -ε^μναβ u_ν B_α u_β
The EMF (electromotive force) components emf[NDIM][N1+D1][N2+D2][N3+D3] used in the constrained transport scheme to evolve the magnetic fields while maintaining ∇·B = 0. This is evident from the code in  defs.h:
These EMF components are used in the numerical scheme for magnetic field evolution, but they are not physical electric fields in the sense of being independent variables. They are computed quantities used in the numerical method, specifically for the constrained transport algorithm.


Special notes:
The staggering is crucial for maintaining ∇·B = 0 constraint
Interpolation between different grid locations is handled by functions in  interp.c
The face-centered magnetic fields are evolved using constrained transport method
Metric terms are computed at different locations as needed for various calculations
This staggered arrangement is typical for GRMHD codes and follows the Constrained Transport approach for maintaining magnetic field divergence-free condition to machine precision.
*/