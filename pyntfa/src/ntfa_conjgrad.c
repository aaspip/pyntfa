
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "ntfa_blas.h"
#include "ntfa_alloc.h"
#include "ntfa_conjgrad.h"

static int np, nx, nr, nd;
static float *r, *sp, *sx, *sr, *gp, *gx, *gr;
static float eps, tol;
static bool verb, hasp0;

void tf_adjnull (bool adj /* adjoint flag */, 
		 bool add /* addition flag */, 
		 int nx   /* size of x */, 
		 int ny   /* size of y */, 
		 float* x, 
		 float* y) 
/*< Zeros out the output (unless add is true). 
  Useful first step for any linear operator. >*/
{
    int i;
    
    if(add) return;
    
    if(adj) {
	for (i = 0; i < nx; i++) {
	    x[i] = 0.;
	}
    } else {
	for (i = 0; i < ny; i++) {
	    y[i] = 0.;
	}
    }
}


void tf_conjgrad_init(int np1     /* preconditioned size */, 
		      int nx1     /* model size */, 
		      int nd1     /* data size */, 
		      int nr1     /* residual size */, 
		      float eps1  /* scaling */,
		      float tol1  /* tolerance */, 
		      bool verb1  /* verbosity flag */, 
		      bool hasp01 /* if has initial model */) 
/*< solver constructor >*/
{
    np = np1; 
    nx = nx1;
    nr = nr1;
    nd = nd1;
    eps = eps1*eps1;
    tol = tol1;
    verb = verb1;
    hasp0 = hasp01;

    r = tf_floatalloc(nr);  
    sp = tf_floatalloc(np);
    gp = tf_floatalloc(np);
    sx = tf_floatalloc(nx);
    gx = tf_floatalloc(nx);
    sr = tf_floatalloc(nr);
    gr = tf_floatalloc(nr);
}

void tf_conjgrad_close(void) 
/*< Free allocated space >*/
{
    free (r);
    free (sp);
    free (gp);
    free (sx);
    free (gx);
    free (sr);
    free (gr);
}

void tf_conjgrad(tf_operator prec  /* data preconditioning */, 
		 tf_operator oper  /* linear operator */, 
		 tf_operator shape /* shaping operator */, 
		 float* p          /* preconditioned model */, 
		 float* x          /* estimated model */, 
		 float* dat        /* data */, 
		 int niter         /* number of iterations */) 
/*< Conjugate gradient solver with shaping >*/
{
    double gn, gnp, alpha, beta, g0, dg, r0;
    float *d=NULL;
    int i, iter;
    
    if (NULL != prec) {
	d = tf_floatalloc(nd); 
	for (i=0; i < nd; i++) {
	    d[i] = - dat[i];
	}
	prec(false,false,nd,nr,d,r);
    } else {
	for (i=0; i < nr; i++) {
	    r[i] = - dat[i];
	}
    }
    
    if (hasp0) { /* initial p */
	shape(false,false,np,nx,p,x);
	if (NULL != prec) {
	    oper(false,false,nx,nd,x,d);
	    prec(false,true,nd,nr,d,r);
	} else {
	    oper(false,true,nx,nr,x,r);
	}
    } else {
	for (i=0; i < np; i++) {
	    p[i] = 0.;
	}
	for (i=0; i < nx; i++) {
	    x[i] = 0.;
	}
    } 
    
    dg = g0 = gnp = 0.;
    r0 = cblas_dsdot(nr,r,1,r,1);
    if (r0 == 0.) {
	if (verb) printf("zero residual: r0=%g\n",r0);
	return;
    }

    for (iter=0; iter < niter; iter++) {
	for (i=0; i < np; i++) {
	    gp[i] = eps*p[i];
	}
	for (i=0; i < nx; i++) {
	    gx[i] = -eps*x[i];
	}

	if (NULL != prec) {
	    prec(true,false,nd,nr,d,r);
	    oper(true,true,nx,nd,gx,d);
	} else {
	    oper(true,true,nx,nr,gx,r);
	}

	shape(true,true,np,nx,gp,gx);
	shape(false,false,np,nx,gp,gx);

	if (NULL != prec) {
	    oper(false,false,nx,nd,gx,d);
	    prec(false,false,nd,nr,d,gr);
	} else {
	    oper(false,false,nx,nr,gx,gr);
	}

	gn = cblas_dsdot(np,gp,1,gp,1);

	if (iter==0) {
	    g0 = gn;

	    for (i=0; i < np; i++) {
		sp[i] = gp[i];
	    }
	    for (i=0; i < nx; i++) {
		sx[i] = gx[i];
	    }
	    for (i=0; i < nr; i++) {
		sr[i] = gr[i];
	    }
	} else {
	    alpha = gn / gnp;
	    dg = gn / g0;

	    if (alpha < tol || dg < tol) {
		if (verb) 
		    printf(
			"convergence in %d iterations, alpha=%g, gd=%g\n",
			iter,alpha,dg);
		break;
	    }

	    cblas_saxpy(np,alpha,sp,1,gp,1);
	    cblas_sswap(np,sp,1,gp,1);

	    cblas_saxpy(nx,alpha,sx,1,gx,1);
	    cblas_sswap(nx,sx,1,gx,1);

	    cblas_saxpy(nr,alpha,sr,1,gr,1);
	    cblas_sswap(nr,sr,1,gr,1);
	}

	beta = cblas_dsdot(nr,sr,1,sr,1) + eps*(cblas_dsdot(np,sp,1,sp,1) - cblas_dsdot(nx,sx,1,sx,1));
	
	if (verb) printf("iteration %d res: %f grad: %f\n",
			     iter,cblas_snrm2(nr,r,1)/r0,dg);

	alpha = - gn / beta;

	cblas_saxpy(np,alpha,sp,1,p,1);
	cblas_saxpy(nx,alpha,sx,1,x,1);
	cblas_saxpy(nr,alpha,sr,1,r,1);

	gnp = gn;
    }

    if (NULL != prec) free (d);

}

