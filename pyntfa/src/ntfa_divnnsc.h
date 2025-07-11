
#ifndef _divnn_sc_h
#define _divnn_sc_h


void ntrianglen2_init(int nw      /* number of components */, 
		      int n       /* data size */,
		      int nf      /* first dimension */,
		      int *nbox   /* smoothing radius [nw] */);
/*< initialization >*/


void ntriangle2_close(void);
/*< free allocated storage >*/


void ntriangle2_lop (bool adj, bool add, int nx, int ny, float *x, float *y);
/*< combined linear operator (1-dimensional smoothing) >*/


void divnn_sc_init(int nw       /* number of components */, \
		     int ndim     /* number of dimensions */,
		     int n        /* data size */,
		     int *ndat    /* data dimensions [ndim] */, 
		     int *nbox    /* smoothing radius [nw] */,
		     float* den   /* denominator [nw*nd] */,
		     bool verb    /* verbosity flag */);
/*< initialize >*/


void divnn_sc_init2(int nw       /* number of components */, \
		     int ndim     /* number of dimensions */,
		     int n        /* data size */,
		     int *ndat    /* data dimensions [ndim] */, 
		     int *nbox    /* smoothing radius [nw] */,
		  	 float **rct /* triangle lengths [ndim][nd] */,
          	 int **sft /* triangle shifts [ndim][nd] */,
		     float* den   /* denominator [nw*nd] */,
		     bool verb    /* verbosity flag */);
/*< initialize >*/


void divnn_sc_close (void);
/*< free allocated storage >*/

void divnn_sc_close2 (void);
/*< free allocated storage >*/

void divnn_sc (float* num  /* numerator */, 
		   float* rat  /* ratio */, 
		   int niter   /* number of iterations */);
/*< smoothly divide num/rat >*/


void divnn_sc2 (float* num  /* numerator */, 
		   float* rat  /* ratio */, 
		   int niter   /* number of iterations */);
/*< smoothly divide num/rat >*/

#endif
