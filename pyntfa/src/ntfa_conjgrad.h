
#ifndef _tf_conjgrad_h
#define _tf_conjgrad_h

typedef void (*tf_operator)(bool,bool,int,int,float*,float*);
typedef void (*tf_solverstep)(bool,int,int,float*,
			   const float*,float*,const float*);
typedef void (*tf_weight)(int,const float*,float*);
/*^*/

void tf_adjnull (bool adj /* adjoint flag */, 
		 bool add /* addition flag */, 
		 int nx   /* size of x */, 
		 int ny   /* size of y */, 
		 float* x, 
		 float* y);
/*< Zeros out the output (unless add is true). 
  Useful first step for any linear operator. >*/

void tf_conjgrad_init(int np1     /* preconditioned size */, 
		      int nx1     /* model size */, 
		      int nd1     /* data size */, 
		      int nr1     /* residual size */, 
		      float eps1  /* scaling */,
		      float tol1  /* tolerance */, 
		      bool verb1  /* verbosity flag */, 
		      bool hasp01 /* if has initial model */);
/*< solver constructor >*/


void tf_conjgrad_close(void);
/*< Free allocated space >*/


void tf_conjgrad(tf_operator prec  /* data preconditioning */, 
		 tf_operator oper  /* linear operator */, 
		 tf_operator shape /* shaping operator */, 
		 float* p          /* preconditioned model */, 
		 float* x          /* estimated model */, 
		 float* dat        /* data */, 
		 int niter         /* number of iterations */);
/*< Conjugate gradient solver with shaping >*/


void tf_conjgrad_adj(bool adj /* adjoint flag */,
		     tf_operator oper  /* linear operator */, 
		     tf_operator shape /* shaping operator */, 
		     float* p          /* preconditioned model */, 
		     float* x          /* estimated model */, 
		     float* dat        /* data */, 
		     int niter         /* number of iterations */);
/*< Conjugate gradient solver with shaping and its adjoint. >*/

#endif
