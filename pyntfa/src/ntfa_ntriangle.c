/* Non-stationary triangle smoothing */
// #include "_bool.h"
// #include "alloc.h"

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "ntfa_alloc.h"

#include "ntfa_ntriangle.h"


tf_ntriangle tf_ntriangle_init (int nbox /* maximum triangle length */, 
			     int ndat /* data length */)
/*< initialize >*/
{
    tf_ntriangle tr;

    tr = (tf_ntriangle) tf_alloc(1,sizeof(*tr));

    tr->nx = ndat;
    tr->nb = nbox;
    tr->np = ndat + 2*nbox;
    
    tr->tmp = tf_floatalloc(tr->np);

    return tr;
}

static void fold (int o, int d, int nx, int nb, int np, 
		  const float *x, float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	tmp[i+nb] = x[o+i*d];
    
    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+(nx-1-i)*d];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    tmp[j+i] = x[o+i*d];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+i*d];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    tmp[j-1-i] = x[o+(nx-1-i)*d];
    }
}

static void fold2 (int o, int d, int nx, int nb, int np, 
		   float *x, const float* tmp)
{
    int i, j;

    /* copy middle */
    for (i=0; i < nx; i++) 
	x[o+i*d] = tmp[i+nb];

    /* reflections from the right side */
    for (j=nb+nx; j < np; j += nx) {
	for (i=0; i < nx && i < np-j; i++)
	    x[o+(nx-1-i)*d] += tmp[j+i];
	j += nx;
	for (i=0; i < nx && i < np-j; i++)
	    x[o+i*d] += tmp[j+i];
    }
    
    /* reflections from the left side */
    for (j=nb; j >= 0; j -= nx) {
	for (i=0; i < nx && i < j; i++)
	    x[o+i*d] += tmp[j-1-i];
	j -= nx;
	for (i=0; i < nx && i < j; i++)
	    x[o+(nx-1-i)*d] += tmp[j-1-i];
    }
}

static void doubint (int nx, float *xx, bool der)
{
    int i;
    float t;

    /* integrate backward */
    t = 0.;
    for (i=nx-1; i >= 0; i--) {
	t += xx[i];
	xx[i] = t;
    }

    if (der) return;

    /* integrate forward */
    t=0.;
    for (i=0; i < nx; i++) {
	t += xx[i];
	xx[i] = t;
    }
}

static void triple (int o, int d, int nx, int nb, 
		    const float* t, /* triangle radius */
		    const float* s, /* shift */
		    float* x, const float* tmp)
{
    int i, nt, nt1, ns;
    float tt, ss, wt, wt1, w, wp, wm, c, cp, cn;

    for (i=0; i < nx; i++) {
	tt = t[o+i*d]; 			/*t[i]->t[o+i*d]: corrected by Yangkang, Dec 13, 2021*/
	nt = floorf(tt);
	nt1 = nt+1;
	ss = s[o+i*d];			/*s[i]->s[o+i*d]: corrected by Yangkang, Dec 13, 2021*/
	ns = roundf(ss);
	ss -= ns;
	ns += nb;
	wt  = (nt1*nt1-tt*tt)/(nt*nt*(nt+nt1));
	wt1 = (tt*tt-nt*nt)/(nt1*nt1*(nt+nt1));
	c = 2*(wt+wt1)*tmp[i+ns] - 
	    (tmp[i+ns-nt1] + tmp[i+ns+nt1])*wt1 - 
	    (tmp[i+ns-nt]  + tmp[i+ns+nt ])*wt;
	if (ss == 0.0f) {
	    x[o+i*d] = c;
	} else {
	    cp = 2*(wt+wt1)*tmp[i+ns+1] - 
		(tmp[i+ns-nt1+1] + tmp[i+ns+nt1+1])*wt1 - 
		(tmp[i+ns-nt +1] + tmp[i+ns+nt +1])*wt;
	    cn = 2*(wt+wt1)*tmp[i+ns-1] - 
		(tmp[i+ns-nt1-1] + tmp[i+ns+nt1-1])*wt1 - 
		(tmp[i+ns-nt -1] + tmp[i+ns+nt -1])*wt;
	    w = (1.0f-ss)*(1.0f+ss);
	    wp = ss*(ss+1.0f)/2;
	    wm = ss*(ss-1.0f)/2;
	    x[o+i*d] = c*w+cp*wp+cn*wm;
	}	
    }
}


static void triple2 (int o, int d, int nx, int nb, 
		     const float* t, /* triangle radius */
		     const float* s, /* shift */
		     const float* x, float* tmp)
{
    int i, nt, nt1, ns;
    float tt, ss, wt, wt1, w, wp, wm;

    for (i=0; i < nx + 2*nb; i++) {
	tmp[i] = 0;
    }

    for (i=0; i < nx; i++) {
	tt = t[i];
	nt = floorf(tt);
	nt1 = nt+1;
	ss = s[i];
	ns = roundf(ss);
	ss -= ns;
	ns += nb;
	wt  = x[o+i*d]*(nt1*nt1-tt*tt)/(nt*nt*(nt+nt1));
	wt1 = x[o+i*d]*(tt*tt-nt*nt)/(nt1*nt1*(nt+nt1));
	if (ss == 0.0f) {
	    tmp[i+ns-nt1] -= wt1; 
	    tmp[i+ns-nt]  -= wt; 
	    tmp[i+ns]     += 2*(wt+wt1);
	    tmp[i+ns+nt]  -= wt;
	    tmp[i+ns+nt1] -= wt1;
	} else {
	    w = (1.0f-ss)*(1.0f+ss);
	    wp = ss*(ss+1.0f)/2;
	    wm = ss*(ss-1.0f)/2;
	    tmp[i+ns-nt1]   -= wt1*w; 
	    tmp[i+ns-nt]    -= wt*w; 
	    tmp[i+ns]       += 2*(wt+wt1)*w;
	    tmp[i+ns+nt]    -= wt*w;
	    tmp[i+ns+nt1]   -= wt1*w;
	    tmp[i+ns-nt1+1] -= wt1*wp; 
	    tmp[i+ns-nt+1]  -= wt*wp; 
	    tmp[i+ns+1]     += 2*(wt+wt1)*wp;
	    tmp[i+ns+nt+1]  -= wt*wp;
	    tmp[i+ns+nt1+1] -= wt1*wp;
	    tmp[i+ns-nt1-1] -= wt1*wm; 
	    tmp[i+ns-nt-1]  -= wt*wm; 
	    tmp[i+ns-1]     += 2*(wt+wt1)*wm;
	    tmp[i+ns+nt-1]  -= wt*wm;
	    tmp[i+ns+nt1-1] -= wt1*wm;
	}    
    }
}

static void double1 (int o, int d, int nx, int nb, 
        const float* t,
        const float* s,
        float* x, const float* tmp)
{
    int i, nt, nt1, ns;
    float tt, wt, wt1;

    for (i=0; i < nx; i++) {
  tt = t[o+i*d];      
  nt = floorf(tt); 
  nt1 = nt+1;     

  ns = nb + s[o+i*d];

  wt  = (nt1-tt)/(nt*nt); 
  wt1 = (tt-nt)/(nt1*nt1); 

  x[o+i*d] = -(tmp[i+ns-nt1] - tmp[i+ns+nt1])*wt1 
             +(tmp[i+ns-nt]  - tmp[i+ns+nt ])*wt;

    }
}

void tf_nsmooth (tf_ntriangle tr   /* smoothing object */, 
		 int o, int d   /* sampling */, 
		 bool der       /* derivative flag */, 
		 const float *t /* triangle lengths */, 
		 const float *s /* triangle shifts */,
		 float *x       /* data (smoothed in place) */)
/*< smooth >*/
{
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp); 
    doubint (tr->np,tr->tmp,der);
    triple (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
}

void tf_nsmooth2 (tf_ntriangle tr   /* smoothing object */, 
		  int o, int d   /* sampling */, 
		  bool der       /* derivative flag */, 
		  const float *t /* triangle lengths */,
		  const float *s /* triangle shifts */,
		  float *x       /* data (smoothed in place) */)
/*< alternative smooth >*/
{
    triple2 (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
    doubint (tr->np,tr->tmp,der);
    fold2 (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp);
}

void tf_ndsmooth (tf_ntriangle tr /* smoothing derivative object */, 
		  int o, int d /* sampling. o: starting index, d: stride in samples for 1/2/3rd dimension; all refer to a correct index in a 1D vector  */, 
		  bool der     /* derivative flag */, 
		  const float *t /* triangle lengths */, 
		  const float *s /* triangle shifts */,
		  float *x     /* data (smoothed in place) */)

/*< smooth derivative >*/
{
    fold (o,d,tr->nx,tr->nb,tr->np,x,tr->tmp); 
    doubint (tr->np,tr->tmp,der);
    double1 (o,d,tr->nx,tr->nb,t,s,x,tr->tmp);
}

void  tf_ntriangle_close(tf_ntriangle tr)
/*< free allocated storage >*/
{
    free (tr->tmp);
    free (tr);
}
