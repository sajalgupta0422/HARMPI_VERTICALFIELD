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
#include "defs.h"
#include "electron.h"

/*****************************************************************/
/*****************************************************************
   main():
   ------

     -- Initializes, time-steps, and concludes the simulation. 
     -- Handles timing of output routines;
     -- Main is main, what more can you say.  

-*****************************************************************/
int main(int argc,char *argv[])
{
    int i,j,k,m;
	int nfailed = 0 ;
    double ndt1, ndt2, ndt3;
    int is_restarted;
    
  
    mpi_init(argc,argv);
    //memory allocation for rdump
    initialize_parallel_write(0);

#ifdef _OPENMP
#pragma omp parallel default(none) shared(nthreads) private(threadid)
    {
        threadid = omp_get_thread_num();
        printf("tid = %d\n", threadid);
        if(threadid==0) {
          nthreads = omp_get_num_threads();
          printf("nthreads = %d\n",nthreads);
        }
    }
#else
    if(MASTER==mpi_rank) printf("_OPENMP not defined: running without OPENMP\n");
#endif
	nstep = 0 ;
    
	/* Perform Initializations, either directly or via checkpoint */
	if(MASTER==mpi_rank) system("mkdir -p dumps images");
    is_restarted = restart_init();
	if(!is_restarted) {
        init() ;
#if( eHEAT || eCOND || DOKTOT )
        if(WHICHPROBLEM != LIN_MODES && WHICHPROBLEM != ATM_TEST && WHICHPROBLEM != BONDI_CON){
          init_entropy(); //defined in electron.c 
        }
#endif
        tdump = t+DTd ;
        trdump = t+DTr ;
        timage = t+DTi ;
        tlog = t+DTl ;
    }
    tavgstart = t;

    //memory allocation for dump and gdump write buffers
    initialize_parallel_write(1);
  
    if (!is_restarted) {
        /* do initial diagnostics */
        diag(INIT_OUT) ;
    }
    else {
        diag(DIVB_OUT) ;
    }

    defcon = 1.;
	while(t < tf) {
        
        /* step variables forward in time */
        nstroke = 0 ;
        
        step_ch(&ndt1,&ndt2,&ndt3) ;
        
        //MPIMARK: need to add up nstroke's from all MPI processes AND also to *define* nstroke
        
        if(MASTER==mpi_rank) {
            fprintf(stderr,"%10.5g %10.5g (%10.5g,%10.5g,%10.5g) %8d %10.5g\n",t,dt,ndt1,ndt2,ndt3,nstep,
                    nstroke/(2.*mpi_ntot[1]*mpi_ntot[2]*mpi_ntot[3])) ;
        }
        
        //MPIMARK: need to sync up t's between MPI processes
        /* Handle output frequencies: */
        if(t >= tdump) {
            diag(DUMP_OUT) ;
            tdump += DTd ;
        }
        if(t >= timage) {
            diag(IMAGE_OUT) ;
            timage += DTi ;
        }
        if(t >= tlog) {
            diag(LOG_OUT) ;
            tlog += DTl ;
        }
        if(t >= trdump) {
            trdump += DTr ;
            diag(RDUMP_OUT) ;
        }
        
        /* restart dump */
        nstep++ ;
        if(nstep%DTr01 == 0) {
            //write out particles for restart
#if(DOPARTICLES)
            pdump(-1,0);
#endif
            restart_write(-1) ;
        }
        
        
        /* deal with failed timestep, though we usually exit upon failure */
        if(failed) {
            restart_init() ;
            failed = 0 ;
            nfailed = nstep ;
            defcon = 0.3 ;
        }
        if(nstep > nfailed + DTr*4.*(1 + 1./defcon)) defcon = 1. ;


	}
        //MPIMARK: need to add up all the steps from all MPI processes
	if(MASTER==mpi_rank)
          fprintf(stderr,"ns,ts: %d %d\n",nstep,nstep*mpi_ntot[1]*mpi_ntot[2]*mpi_ntot[3]) ;

	/* do final diagnostics */
	diag(FINAL_OUT) ;
//#ifdef _OPENMP
//     }
//#endif
    de_initialize_parallel_write();
#ifdef MPI
        MPI_Finalize();
#endif

	return(0) ;
}


/*****************************************************************/
/*****************************************************************
  set_arrays():
  ----------

       -- sets to zero all arrays, plus performs pointer trick 
          so that grid arrays can legitimately refer to ghost 
          zone quantities at positions  i = -2, -1, N1, N1+1 and 
          j = -2, -1, N2, N2+1

 *****************************************************************/
//SASMARK: zero out everything, including uelval
//SASMARK: remove excess variables, such as uelval perhaps?
void set_arrays()
{
  int i,j,k,m ;
  //Manipulating pointers to access "ghost zone" elements of grid arrays
  /*
  p: This is a pointer variable. Based on the cast on the RHS, it's a pointer to a 4-dimensional array of double with dimensions [N2M][N3M][NPR]. 
  Notice that the first dimension N1M is not present in the pointer type.
(double (*) [N2M][N3M][NPR]): This is a type cast. It tells the compiler to treat the memory address on the RHS as a pointer to a 4D array of double with the specified inner dimensions.
&(): This is the "address-of" operator. It takes the memory address of the expression inside the parentheses.
a_p[N1G][N2G][N3G][0]: This accesses a specific element within the a_p array, where a_p is declared as extern double a_p[N1M][N2M][N3M][NPR]. 
    N1G, N2G, and N3G are preprocessor macros representing indices along the first three dimensions of a_p. 
    These values are probably chosen such that they correspond to the start of the "active" computational domain, offset from the beginning of the allocated array to accommodate the ghost zones.
    [0] selects the first element along the last dimension (NPR)

    Let's visualize this for a simpler 1D case. Suppose you have an array data[N1M] and you want to have a ghost zone of size 2 on each side. 
    You might allocate an array of size N1M + 4. If N1G was 2, then &data[N1G] would point to the beginning of your "active" data, 
    and you could access the ghost cells using negative indices relative to this pointer (e.g., p[-1], p[-2]) 
    and indices beyond N1M - 1 (which would map to the upper ghost zones in the underlying data array).

  For example, if the code later accesses p[i][j][k][m], and if i ranges from -2 to N1 + 1, 
  this will actually correspond to accessing elements in the a_p array with indices ranging around N1G + i, 
  and similarly for the other dimensions. By choosing N1G, N2G, and N3G appropriately , the code can use a more 
  natural indexing scheme (-2 to N1+1, etc.) while still accessing the memory allocated for the ghost zones in the a_p array.
  */
  p     =  (double (*) [N2M][N3M][NPR])(&(  a_p[N1G][N2G][N3G][0])) ; 
  dq    =  (double (*) [N2M][N3M][NPR])(&( a_dq[N1G][N2G][N3G][0])) ;
  dsource =(double (*) [N2M][N3M])(&( a_dsource[N1G][N2G][N3G])) ;
  duscon =(double (*) [N2M][N3M])(&( a_duscon[N1G][N2G][N3G])) ;
  sour =(double (*) [N2M][N3M])(&( a_sour[N1G][N2G][N3G])) ;
  uelvar =(double (*) [N2M][N3M])(&( a_uelvar[N1G][N2G][N3G])) ;
  qisosq =(double (*) [N2M][N3M])(&( a_qisosq[N1G][N2G][N3G])) ;
  qisodotb =(double (*) [N2M][N3M])(&( a_qisodotb[N1G][N2G][N3G])) ;
  qdot =(double (*) [N2M][N3M])(&( a_qdot[N1G][N2G][N3G])) ;
#if(DOQDOTAVG)
    qdotavg =(double (*) [N2M][N3M])(&( a_qdotavg[N1G][N2G][N3G])) ;
    feqdotavg =(double (*) [N2M][N3M])(&( a_feqdotavg[N1G][N2G][N3G])) ;
#endif
  thetaold =(double (*) [N2M][N3M])(&( a_thetaold[N1G][N2G][N3G])) ;
  delvx =(double complex (*) [N2M])(&( a_delvx[N1G][N2G])) ;
  delvy =(double complex (*) [N2M])(&( a_delvy[N1G][N2G])) ;

  uelarray =(double (*) [N2M][N3M][KELDIS-KEL4+1])(&(  a_uelarray[N1G][N2G][N3G][0])) ;
  du_flr =(double (*) [N2M][N3M])(&(  a_du_flr[N1G][N2G][N3G])) ;

  
  
  F1    =  (double (*) [N2M][N3M][NPR])(&( a_F1[N1G][N2G][N3G][0])) ;
  F2    =  (double (*) [N2M][N3M][NPR])(&( a_F2[N1G][N2G][N3G][0])) ;
  F3    =  (double (*) [N2M][N3M][NPR])(&( a_F3[N1G][N2G][N3G][0])) ;
  ph    =  (double (*) [N2M][N3M][NPR])(&( a_ph[N1G][N2G][N3G][0])) ;
  psave =(double (*) [N2M][N3M][NPR])(&( a_psave[N1G][N2G][N3G][0])) ;
  pbound =(double (*) [N2M][N3M][NPR])(&( a_pbound[N1G][N2G][N3G][0])) ;

  
  pflag =  (int    (*) [N2M][N3M])(     &( a_pflag[N1G][N2G][N3G] )) ;

#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) default(none) collapse(3) shared(p,ph,psave,pbound,dq,F1,F2,F3,pflag,dsource,duscon,sour,qisosq,qisodotb,delvx,delvy,nthreads) private(i,j,k,m)
  /* everything must be initialized to zero */
  ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G-1) {
    PLOOP {
      p[i][j][k][m]   = 0. ;
      ph[i][j][k][m]  = 0. ;
      psave[i][j][k][m] = 0.;
      pbound[i][j][k][m] = 0.;

      dq[i][j][k][m]  = 0. ;
      
      F1[i][j][k][m]  = 0. ;
      F2[i][j][k][m]  = 0. ;
      F3[i][j][k][m]  = 0. ;
    }
    pflag[i][j][k] = 0 ;
    dsource[i][j][k] =0;
    duscon[i][j][k] = 0;
    sour[i][j][k] = 0;
    qisosq[i][j][k] = 0;
    qisodotb[i][j][k] = 0;
    delvx[i][j] = 0.;
    delvy[i][j] = 0.;
    
  }

#if(DOQDOTAVG)
#pragma omp parallel for schedule(static,N1M*N2M*N3M/nthreads) default(none) collapse(3) shared(qdotavg,feqdotavg,nthreads) private(i,j,k,m)
    /* everything must be initialized to zero */
    ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G-1) {
        PLOOP {
            qdotavg[i][j][k]  = 0. ;
            feqdotavg[i][j][k]  = 0. ;
        }
    }
#endif

    
#if(EVOLVEVPOT)
  //zero out vpot
#pragma omp parallel for schedule(static,(N1+D1)*(N2+D2)*(N3+D3)/nthreads) collapse(3) default(none) shared(vpot,nthreads) private(i,j,k)
  ZSLOOP(0,N1-1+D1,0,N2-1+D2,0,N3-1+D3) {
    vpot[0][i][j][k] = 0.;
    vpot[1][i][j][k] = 0.;
    vpot[2][i][j][k] = 0.;
    vpot[3][i][j][k] = 0.;
  }
#endif
  
  ZLOOP for(m = 0; m < NIMG; m++) {
    failimage[i][j][k][m] = 0L ;
  }
  
  
  /* grid functions */
  conn = (double (*) [N2M][N3M][NDIM][NDIM][NDIM])
  (& ( a_conn[N1G][N2G][N3G][0][0][0])) ;
  gcon = (double (*) [N2M][N3M][NPG][NDIM][NDIM])
  (& ( a_gcon[N1G][N2G][N3G][0][0][0])) ;
  gcov = (double (*) [N2M][N3M][NPG][NDIM][NDIM])
  (& ( a_gcov[N1G][N2G][N3G][0][0][0])) ;
  gdet = (double (*) [N2M][N3M][NPG])
  (& ( a_gdet[N1G][N2G][N3G][0])) ;
  phys_coords = (double (*) [N1M][N2M][N3M])
  (& ( a_phys_coords[0][N1G][N2G][N3G])) ;

}


/*****************************************************************/
/*****************************************************************
  set_grid():
  ----------

       -- calculates all grid functions that remain constant 
          over time, such as the metric (gcov), inverse metric 
          (gcon), connection coefficients (conn), and sqrt of 
          the metric's determinant (gdet).

 *****************************************************************/
void set_grid()
{
	int i,j,k,m,dim ;
	double X[NDIM], V[NDIM] ;
	struct of_geom geom ;

	/* set up boundaries, steps in coordinate grid */
	set_points() ;
	dV = dx[1]*dx[2]*dx[3] ;
        //printf("dx[1] = %g, dx[2] = %g, dx[3] = %g\n",dx[1],dx[2],dx[3]);
        if(MASTER == mpi_rank){
          printf("Starting grid pre-computation...");
          fflush(stdout);
        }

	DLOOPA X[j] = 0. ;

	ZSLOOP(-N1G,N1+N1G-1,-N2G,N2+N2G-1,-N3G,N3+N3G-1) {
		
		/* zone-centered */
		coord(i,j,k,CENT,X) ;
		gcov_func(X, gcov[i][j][k][CENT]) ;
		gdet[i][j][k][CENT] = gdet_func(gcov[i][j][k][CENT]) ;
		gcon_func(gcov[i][j][k][CENT], gcon[i][j][k][CENT]) ;

    bl_coord_vec(X, V);
    for(dim=0;dim<NDIM;dim++) phys_coords[dim][i][j][k] = V[dim];
          
		get_geometry(i,j,k,CENT,&geom) ;
		conn_func(X, &geom, conn[i][j][k]) ;

		/* corner-centered */
		coord(i,j,k,CORN,X) ;
		gcov_func(X,gcov[i][j][k][CORN]) ;
		gdet[i][j][k][CORN] = gdet_func(gcov[i][j][k][CORN]) ;

		gcon_func(gcov[i][j][k][CORN],gcon[i][j][k][CORN]) ;

		/* r-face-centered */
		coord(i,j,k,FACE1,X) ;
		gcov_func(X,gcov[i][j][k][FACE1]) ;
		gdet[i][j][k][FACE1] = gdet_func(gcov[i][j][k][FACE1]) ;
		gcon_func(gcov[i][j][k][FACE1],gcon[i][j][k][FACE1]) ;

		/* theta-face-centered */
		coord(i,j,k,FACE2,X) ;
		gcov_func(X,gcov[i][j][k][FACE2]) ;
		gdet[i][j][k][FACE2] = gdet_func(gcov[i][j][k][FACE2]) ;
		gcon_func(gcov[i][j][k][FACE2],gcon[i][j][k][FACE2]) ;

    /* phi-face-centered */
    coord(i,j,k,FACE3,X) ;
    gcov_func(X,gcov[i][j][k][FACE3]) ;
    gdet[i][j][k][FACE3] = gdet_func(gcov[i][j][k][FACE3]) ;
    gcon_func(gcov[i][j][k][FACE3],gcon[i][j][k][FACE3]) ;

	
        }

        if(MASTER == mpi_rank){
          printf(" done!\n");
          fflush(stdout);
        }
	/* done! */
}

