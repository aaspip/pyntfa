
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "ntfa_triangle.h"
#include "ntfa_blas.h"
#include "ntfa_alloc.h"
#include "ntfa_conjgrad.h"

static int nw;
static float **w;

void tf_weight2_init(int nw1   /* number of components */, 
		     int n     /* model size */, 
		     float *ww /* weight [nw*n] */)
/*< initialize >*/
{
    int iw;

    nw = nw1;
    w = (float**) tf_alloc(nw,sizeof(float*));

    for (iw=0; iw < nw; iw++) {
	w[iw] = ww+iw*n;
    }
}

void tf_weight2_close(void)
/*< free allocated storage >*/
{
    free(w);
}

void tf_weight2_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy)
/*< linear operator >*/
{
    int i, iw;

    if (nw*ny != nx) printf("tf_weight2_lop size mismatch: %d*%d != %d\n",
			      nw,ny,nx);

    tf_adjnull (adj, add, nx, ny, xx, yy);
  
    if (adj) {
        for (iw=0; iw < nw; iw++) {
	    for (i=0; i < ny; i++) {
	        xx[i+iw*ny] += yy[i] * w[iw][i];
	    }
	}
    } else {
        for (iw=0; iw < nw; iw++) {
	    for (i=0; i < ny; i++) {
	        yy[i] += xx[i+iw*ny] * w[iw][i];
	    }
	}
    }
}