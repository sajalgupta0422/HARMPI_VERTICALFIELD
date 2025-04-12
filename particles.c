
#include "decs.h"

/*

 advance particle positions using fluid half-step primitives 

 revised cfg 11 apr 2016 to improve interpolation scheme 

*/
#if(DOPARTICLES)
void advance_particles(double ph[N1M][N2M][N3M][NPR], double Dt)
{
	int k, l, i, j;
	double ucon[NDIM], vel[NDIM];
	double f1,f2,f3;
	struct of_geom geom;

	for (l = 0; l < myNp; l++) {
#if(BL)
                //ensure that phi-periodicity is in place
                //SASMARK: only works for phi-periodic grids, so do this only for BL!=0
                if (xp[l][3] >= startx[3] + lenx[3]) {
                    xp[l][3] -= lenx[3];
                }
                else if (xp[l][3] < startx[3]) {
                    xp[l][3] += lenx[3];
                }
#endif
		/* don't update particles that are off-grid for this MPI process */
		if(xp[l][1] >= mpi_startx[1] && xp[l][2] >= mpi_startx[2] && xp[l][3] >= mpi_startx[3] &&
		   xp[l][1]  < mpi_stopx[1]  && xp[l][2]  < mpi_stopx[2]  && xp[l][3]  < mpi_stopx[3] ) {

			/* the four-velocities are zone-centered */
			f1 = (xp[l][1] - mpi_startx[1] + 0.5*dx[1]) / dx[1];
			f2 = (xp[l][2] - mpi_startx[2] + 0.5*dx[2]) / dx[2];
                        f3 = (xp[l][3] - mpi_startx[3] + 0.5*dx[3]) / dx[3];

		   	/* find nearest zone center */
			i = lround( f1 ) ;
			j = lround( f2 ) ;
                        k = lround( f3 ) ;

                        get_geometry(i, j, k, CENT, &geom);
                        ucon_calc(ph[i][j][k], &geom, ucon);
                        for (k = 1; k < NDIM; k++) vel[k] = ucon[k]/ucon[0];

			/* push particle forward */
			for (k = 1; k < NDIM; k++)
				xp[l][k] += Dt * vel[k];

		}

	}
	/* done! */
}


/*

 initialize Lagrangian tracer particles
 cfg 4 feb 09

 simplified 10 apr 2016 cfg

*/

void init_particles()
{
    int i, j, k, Np, newNp;
    double dmass, mass, Nexp, newNexp, sample_factor, X[NDIM];
    struct of_geom geom ;
    struct of_state q ;
    double U[NPR];
    int iglob, jglob, kglob;
    double rancval1, rancval2, rancval3;

    /* global variables */
    //pdump_cnt = 0;
    //DTp = 20.;

    /* assign particles according to restmass density */
    /* first find total rest-mass on grid */
    mass = 0.;
    ZLOOP {
        get_geometry(i, j, k, CENT, &geom) ;
        get_state(p[i][j][k], &geom, &q);
        primtoU(p[i][j][k], &q, &geom, U) ;

        dmass = U[RHO]*dx[1]*dx[2]*dx[3] ;
        mass += dmass;
    }
#if MPI
    //exchange the info between the MPI processes to get the total sum
    MPI_Allreduce(MPI_IN_PLACE,&mass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif


    /* normalization factor so that we get the # of
       particles we want */
    sample_factor = NPTOT / mass;

    /* now sample, keeping running total of particles
       created */
    Np = 0;
    Nexp = 0.;
    myNp = 0; //reset the particle counter
//        //ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
    for(iglob=0;iglob<mpi_ntot[1];iglob++) {
        for(jglob=0;jglob<mpi_ntot[2];jglob++) {
            for(kglob=0;kglob<mpi_ntot[3];kglob++) {
                i = iglob-mpi_startn[1];
                j = jglob-mpi_startn[2];
                k = kglob-mpi_startn[3];
                if(i<0 ||
                   j<0 ||
                   k<0 ||
                   i>=N1 ||
                   j>=N2 ||
                   k>=N3){
                    newNp = 0;
                    newNexp = 0.;

                    if (N1 - 1 == iglob % N1 || N2 - 1 == jglob % N2 || N3 - 1 == kglob % N3) {
#if MPI
                        //by design, this MPI process will be stuck until just before
                        //switching tiles (see below)
                        MPI_Allreduce(MPI_IN_PLACE,&newNp,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                        MPI_Allreduce(MPI_IN_PLACE,&newNexp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif
                    }
                    else {
                        continue;
                    }
                    //throw the dice the same number of times on all MPI processes to
                    //make the process independent of MPI partitioning
                    for( ; Np < newNp; Np++) {
                        rancval1 = ranc(0);
                        rancval2 = ranc(0);
                        rancval3 = ranc(0);
                    }
                    Np = newNp;
                    Nexp = newNexp;
                }
                else {
                    get_geometry(i, j, k, CENT, &geom) ;
                    get_state(p[i][j][k], &geom, &q);
                    primtoU(p[i][j][k], &q, &geom, U) ;

                    dmass = U[RHO]*dx[1]*dx[2]*dx[3] ;

                    Nexp += dmass * sample_factor;

                    /* assign particle to random position in cell */
                    while(Nexp >= 0.5) {
                        /* lower left corner of zone is here */
                        coord(i, j, k, CORN, X);

                        rancval1 = ranc(0);
                        rancval2 = ranc(0);
                        rancval3 = ranc(0);

                        xp[myNp][0] = (double) Np; //particle number is its tag
                        xp[myNp][1] = X[1] + rancval1 * dx[1];
                        xp[myNp][2] = X[2] + rancval2 * dx[2];
                        xp[myNp][3] = X[3] + rancval3 * dx[3];

                        Np++;
                        myNp++;
                        Nexp -= 1.;
                    }
#if MPI
                    //only do this just before switching tiles
                    //the rest of the processes will idly wait
                    if (N1 - 1 == iglob % N1 || N2 - 1 == jglob % N2 || N3 - 1 == kglob % N3) {
                        MPI_Allreduce(MPI_IN_PLACE,&Np,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
                        MPI_Allreduce(MPI_IN_PLACE,&Nexp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                    }
#endif
                }
            }
        }
    }

    /* report back! */
    fprintf(stderr, "[MPI rank %d] made %d particles of %d (%g). total made: %d\n", mpi_rank, myNp, NPTOT, Nexp, Np);

    /* done! */
}
#endif