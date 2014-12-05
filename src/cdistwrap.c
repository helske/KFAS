#include <R.h>
#include <Rmath.h>

void F77_NAME(dnormf)(double *x,double *mu, double *sigma, double *res);
void F77_NAME(dpoisf)(double *x,double *lambda, double *res);
void F77_NAME(dbinomf)(double *x, double *n, double *p, double *res);
void F77_NAME(dgammaf)(double *x, double *shape, double *scale, double *res);
void F77_NAME(dnbinomf)(double *x, double *size, double *mu, double *res);
int F77_NAME(finitex)(double *x);

void F77_SUB(dnormf)(double *x, double *mu, double *sigma, double *res){ 
 *res -= dnorm4(*x, *mu, *sigma, 1); }
 
void F77_SUB(dpoisf)(double *x, double *lambda, double *res){ 
 *res += dpois(*x, *lambda, 1); }
 
void F77_SUB(dbinomf)(double *x, double *n, double *p, double *res){ 
 *res += dbinom(*x, *n, *p, 1); }
 
void F77_SUB(dgammaf)(double *x, double *shape, double *scale, double *res){ 
 *res += dgamma(*x, *shape, *scale, 1); }
 
void F77_SUB(dnbinomf)(double *x, double *size, double *mu, double *res){ 
 *res += dnbinom_mu(*x, *size, *mu, 1); }

int F77_SUB(finitex)(double *x){
	return R_FINITE(*x); }
