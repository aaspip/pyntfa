

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "ntfa_alloc.h"
#include "ntfa_conjgrad.h"
#include "ntfa_ntriangle.h"
#include "ntfa_ntrianglen.h"
#include "ntfa_decart.h"

static int *n, s[9], nd, dim, nrep;
static tf_ntriangle *tr;
static float *tmp, **tlen, **tsft;

void tf_ntrianglen_init (int ndim    /* number of dimensions */, 
		      int *nbox   /* triangle radius [ndim] */, 
		      int *ndat   /* data dimensions [ndim] */,
		      float **len /* triangle lengths [ndim][nd] */,
              float **sft /* triangle shifts [ndim][nd] */,
		      int repeat  /* repeated smoothing */)
/*< initialize >*/
{
    int i;

    n = ndat;
    dim = ndim;

    tr = (tf_ntriangle*) tf_alloc(dim,sizeof(tf_ntriangle));
    
    nd = 1;
    for (i=0; i < dim; i++) {
	tr[i] = (nbox[i] > 1)? tf_ntriangle_init (nbox[i],ndat[i]): NULL;
	s[i] = nd;
	nd *= ndat[i];
    }
    tlen = len; 
    tsft = sft;

    tmp = tf_floatalloc(nd);
    nrep = repeat;
}

void tf_ntrianglen_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i, j, i0, irep;

    if (nx != ny || nx != nd) 
	printf("tf_ntrianglen_lop: Wrong data dimensions: nx=%d, ny=%d, nd=%d\n",
		 nx,ny,nd);

    tf_adjnull (adj,add,nx,ny,x,y);
  
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    tmp[i] = x[i];
	}
    }

  
    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) {
	    for (j=0; j < nd/n[i]; j++) {
		i0 = tf_first_index (i,j,dim,n,s);

		for (irep=0; irep < nrep; irep++) {
		    tf_nsmooth (tr[i], i0, s[i], false, tlen[i], tsft[i], tmp);
		}
	    }
	}
    }
	
    if (adj) {
	for (i=0; i < nd; i++) {
	    x[i] += tmp[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    y[i] += tmp[i];
	}
    }    
}

void tf_ntrianglen_close(void)
/*< free allocated storage >*/
{
    int i;

    free (tmp);

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) tf_ntriangle_close (tr[i]);
    }

    free(tr);
}