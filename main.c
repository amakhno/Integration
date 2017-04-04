#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#define EPS 3.0e-14 
#define H 1.054572e-27
#define PIM4 0.7511255444649425 
#define MAXIT 10 
#define _GNU_SOURCE
#define OMEGA 1.0
#define F0 1.0
#define ALPHA 4.472135955
#define SIZE 20
#define NUM_POINTS 64
#define NPOINTS 100	//
#define END_ROOT_FIND 30

int gauleg(double x1, double x2, double *x, double *w, int n);
double F(double t);
double IntegrateF(double a, double b, double *w, double *x, double n, int N);
double IntegrateP(double a, double b, double *w, double *x, double n);
double f (double x, void * params);
double A_evaluate(double t, void * params);
double A_integrate(double t1, double t);
double S_integrate(double t1, double t);
double S_integrate_func(double eps, void * params);
double Alpha2(double t1, double t);
double GetEvalSolve(double t);

struct quadratic_params
  {
    double t;
  };

double quadratic (double t0, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double t = p->t;

  return Alpha2(t0, t);
}

int FindRoots(double t, double* roots)
{
	gsl_function AlphaFunc;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	struct quadratic_params params = { t };
	AlphaFunc.function = &quadratic;
	AlphaFunc.params = &params;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);	
	//FILE *o;
	//o = fopen("out.txt", "w");
	
	double h = 0.1;
	double t1_hi = t, t1_lo = t - h;
	int countOfRoots = 0;
	while( t1_lo > -20 )
	{
		if ((Alpha2(t1_lo, t)) * (Alpha2(t1_hi, t)) < 0)
		{
			gsl_root_fsolver_set (s, &AlphaFunc, t1_lo, t1_hi);
			//printf("В промежутке : [%f, %f]\n", t1_lo, t1_lo+h);			
			int iter = 0;
			int status = 0;
			do
			{
				iter++;
				status = gsl_root_fsolver_iterate (s);
				roots[countOfRoots] = gsl_root_fsolver_root (s);
				t1_lo = gsl_root_fsolver_x_lower (s);
				t1_hi = gsl_root_fsolver_x_upper (s);
				status = gsl_root_test_interval (t1_lo, t1_hi,
                                       0, 0.001);
			}
			while (status == GSL_CONTINUE && iter < 1000);
			//printf("Решение: %f\n", roots[countOfRoots]);
			countOfRoots++;
			t1_hi = t1_lo;
			t1_lo -= h;
		}
		//fprintf(o, "%f %f\n", t1_lo, Alpha2(t1_lo, t) );
		t1_hi = t1_lo;
		t1_lo -= h;
	}
	
	gsl_root_fsolver_free (s);
	//fclose(o);
	return countOfRoots;
}

double GetAnalyticSolve(double t, int countOfRoots, double* roots)
{
	fftw_complex sum;
	double c = 3e9;
	double E = -2.0;
	for(int i = 0; i<countOfRoots; i++)
	{
		double d = Alpha2(roots[i], t) * ( f(roots[i], 0) + 1/c * ( 1/( (roots[i] - t)*(roots[i] - t) ) * A_integrate(t, roots[i]) + A_evaluate(roots[i], 0) / (t - roots[i])) );
		sum += cexp(I * S_integrate(roots[i], t) + E)/(d * pow( t - roots[i] , 3.0/2.0));
	}	
	return creal(sum)*creal(sum) + cimag(sum)*cimag(sum);;
}

int main()
{
	FILE *f;
	f = fopen("S.txt", "w");
	FILE *f1;
	f1 = fopen("F1.txt", "w");
	double t = 10;
	double a_integrate_result;
	
	double roots[50];
	
	for(double t = 10; t<15; t+=0.5)
	{
		
		int countOfRoots = FindRoots( t , roots );
		double m1 = GetAnalyticSolve( t , countOfRoots, roots );
		double m2 = GetEvalSolve( t );
		fprintf(f, "%f %f\n", t , m1);
		fprintf(f1, "%f %f\n", t, m2);
	}
	fclose(f);
	fclose(f1);
	return 0;
}

double f (double x, void * params) 
{
  return F0*exp(-(x*x)/(ALPHA*ALPHA))*cos(OMEGA*x);
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

double Alpha2(double t1, double t)
{
	return A_evaluate(t1, 0) - A_integrate(t1, t) / (t1 - t);
}

double S_integrate_func(double eps, void * params)
{
	double temp = A_evaluate(eps, 0);
	return temp*temp;
}

double GetEvalSolve(double t)
{
	fftw_complex *in, *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    
	double tempEnd = 50;
	double step = tempEnd / NPOINTS;
	double tempPower;	
	int i = 0;	
	
	for(double tempT = 0; tempT < tempEnd; tempT+=step)
	{
		in[i] = cexp(I * 1 * tempT) 
			/ cpow(tempT, 3.0/2.0) 
			* (cexp(I * S_integrate(t-tempT, t)) - 1);
			
		if(!tempT)
		{
			in[i] = 0;
		}
		i++;
	}
	
	p = fftw_plan_dft_1d(NPOINTS, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);    
	double result;
    
    for(int i = 0; i<NPOINTS; ++i)
	{
		tempPower = creal(out[i])*creal(out[i]) + cimag(out[i])*cimag(out[i]);													
	}
    
    result = creal(out[0])*creal(out[0]) + cimag(out[0])*cimag(out[0]);
	printf("Попытка очистить OUT\n");
	free(out);
	printf("Попытка очистить PLAN\n");
	fftw_destroy_plan(p);
	return result;
}
