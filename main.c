#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <float.h>
#include <gsl/gsl_integration.h>
#include <complex.h>
#include <fftw3.h>


#define EPS 3.0e-14 
#define H 1.054572e-27
#define PIM4 0.7511255444649425 
#define MAXIT 10 
#define _GNU_SOURCE
#define OMEGA 1
#define ALPHA 4.472135955
#define SIZE 20
#define NUM_POINTS 64
#define NPOINTS 600

int gauleg(double x1, double x2, double *x, double *w, int n);
double F(double t);
double IntegrateF(double a, double b, double *w, double *x, double n, int N);
double IntegrateP(double a, double b, double *w, double *x, double n);
double f (double x, void * params);
double A_evaluate(double t, void * params);
double A_integrate(double t1, double t);
double S_integrate(double t1, double t);
double S_integrate_func(double eps, void * params);
double S_integrate_func2(double eps, void * params);
double M(double eps, double t);
double FFT_M();

int main()
{
	double a_integrate = A_integrate(-100, 0);
	printf ("A_integration = % .18f\n", a_integrate);
	FILE *f;
	f = fopen("S.txt", "w");
	double t = 10;
	double tau = 0;
	double a_integrate_result;
	double temp = A_evaluate(5, (void *)0);
	/*
	 * FFTW BEGIN
	 * */
	fftw_complex *in, *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    
    fftw_complex signal[NUM_POINTS];
    fftw_complex result[NUM_POINTS];
    
    signal[0][0] = 1.0 * cos(10.0 ) + 0.5 * cos(25.0 );
    
    double tempT = 0;
	double tempEnd = 20;
	double step = tempEnd / NPOINTS;
	
	double complex underintegrateFunc;
	
	int i = 0;
	for(double tempT = 0; tempT < tempEnd; tempT+=step)
	{
		underintegrateFunc = (cexp(I * 1 * tempT)/H) 
			/ cpow(tempT, 3.0/2.0) 
			* (cexp(I * S_integrate(t-tempT, t) / H) - 1);
		in[i][0] = underintegrateFunc;
		i++;
	}
	
	
	
	//Old Cycle
	/*for(tau = 0; tau<20; tau+=0.1)
	{
															
		fprintf(f, "%f %f\n", tau, S_integrate(t-tau, t));
	}
	fclose(f);*/
	return 0;
}

double f (double x, void * params) 
{
  return exp(-(x*x)/(ALPHA*ALPHA))*cos(OMEGA*x);
}

double A_evaluate(double t, void * params)
{
	double t1 = -150;
	double result;
	double error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &f;
	F.params = &alpha;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 1e-13, 0, 10000, w, &result, &error); 
	gsl_integration_workspace_free (w);
	return result;
}

double A_integrate(double t1, double t)
{
	double result;
	double error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &A_evaluate;
	F.params = &alpha;
	gsl_integration_workspace * a = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 1e-13, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return result;
	
}

double S_integrate(double t1, double t)
{
	double a_integrate_result = A_integrate(t1, t)/( 2*(t1-t) );
	double result, error;
	gsl_function F;
	F.function = &S_integrate_func;
	F.params = &a_integrate_result;
	gsl_integration_workspace * a = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 1e-13, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return -result/2.0 + a_integrate_result;
}

double S_integrate_func(double eps, void * params)
{
	double temp = A_evaluate(eps, 0);
	return temp*temp;
}

double M_func(double t0, double eps, double t1)
{
	return exp(I * eps * (t1)/H) / pow(t1, 5/2) * ( exp(I*S_integrate(t0, t0-t1)/H) - 1 );
}

double FFT_M()
{
	fftw_complex *in, *out;
    fftw_plan p;
	int N = 200;
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_execute(p); /* repeat as needed */
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
   return 0;
}

