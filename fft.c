//
//  fft.c
//  HARM2D
//
//

#include "decs.h"


/*-------------------------------------------------------------------------
 Perform a 2D FFT inplace given a complex 2D array
 The direction dir, 1 for forward, -1 for reverse
 The size of the array (nx,ny)
 Return false if there are memory problems or
 the dimensions are not powers of 2
 used in step_ch.c
 */
int FFT2D(complex double c[][N2M],int nx,int ny,int dir)
{
    int i,j;
    int m,twopm;
    double *real,*imag;
    
    /* Transform the rows */
    real = (double *)malloc(nx * sizeof(double));
    imag = (double *)malloc(nx * sizeof(double));
    if (real == NULL || imag == NULL){
        printf ("Returning Null\n");

        return(0);
    }
    
    if (!Powerof2(nx,&m,&twopm) || twopm != nx){
        printf ("nx fail: %d\n", nx);
        printf ("twopm: %d\n", twopm);


        return(0);
    }
    
    for (j=0;j<ny;j++) {
        for (i=0;i<nx;i++) {
            real[i] = creal(c[i][j]);
            imag[i] =  cimag(c[i][j]);
        }
        FFT(dir,m,real,imag);
        for (i=0;i<nx;i++) {
            c[i][j]= real[i]+I*imag[i];
        }
    }
    free(real);
    free(imag);
    
    /* Transform the columns */
    real = (double *)malloc(ny * sizeof(double));
    imag = (double *)malloc(ny * sizeof(double));
    if (real == NULL || imag == NULL){
        printf ("Returning Null: %d\n", i);

        return(0);
    }
    if (!Powerof2(ny,&m,&twopm) || twopm != ny){
        printf ("ny fail: %d\n", ny);

        return(0);
    }
    for (i=0;i<nx;i++) {
        for (j=0;j<ny;j++) {
            real[j] = creal(c[i][j]);
            imag[j] = cimag(c[i][j]);
        }
        FFT(dir,m,real,imag);
        for (j=0;j<ny;j++) {
            c[i][j] = real[j]+I*imag[j];
        }
    }
    free(real);
    free(imag);
    
    return(1);
}

/*-------------------------------------------------------------------------
 This computes an in-place complex-to-complex FFT
 x and y are the real and imaginary arrays of 2^m points.
 dir =  1 gives forward transform
 dir = -1 gives reverse transform
 
 Formula: forward
                N-1
                ---
        1      \          - j k 2 pi n / N
 X(n) = --      >   x(k) e                    = forward transform
        N      /                                n=0..N-1
                ---
                k=0
 
 Formula: reverse
            N-1
            ---
            \          j k 2 pi n / N
 X(n) =      >   x(k) e                    = forward transform
            /                                n=0..N-1
            ---
            k=0
 */
int FFT(int dir,int m,double *x,double *y)
{
    long nn,i,i1,j,k,i2,l,l1,l2;
    double c1,c2,tx,ty,t1,t2,u1,u2,z;
    
    /* Calculate the number of points */
    nn = 1;
    for (i=0;i<m;i++)
        nn *= 2;
    
    /* Do the bit reversal */
    i2 = nn >> 1;
    j = 0;
    for (i=0;i<nn-1;i++) {
        if (i < j) {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }
        k = i2;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    
    /* Compute the FFT */
    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for (l=0;l<m;l++) {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j=0;j<l1;j++) {
            for (i=j;i<nn;i+=l2) {
                i1 = i + l1;
                t1 = u1 * x[i1] - u2 * y[i1];
                t2 = u1 * y[i1] + u2 * x[i1];
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            z =  u1 * c1 - u2 * c2;
            u2 = u1 * c2 + u2 * c1;
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }
    
    /* Scaling for forward transform */
    if (dir == 1) {
        for (i=0;i<nn;i++) {
            x[i] /= (double)nn;
            y[i] /= (double)nn;
        }
    }
    
    return(1);
}

/*-------------------------------------------------------------------------
 Calculate the closest but lower power of two of a number
 twopm = 2**m <= n
 Return TRUE if 2**m == n
 */
int Powerof2(int n,int *m,int *twopm)
{
    
    if (n <= 1) {
        *m = 0;
        *twopm = 1;
        return(0);
    }
    
    *m = 1;
    *twopm = 2;
    do {
        (*m)++;
        (*twopm) *= 2;
    } while (2*(*twopm) <= n);
    

    if (*twopm != n)
        return(0);
    else
        return(1);
    
    

}


