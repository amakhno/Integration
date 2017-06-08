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
#include "omp.h"
 
#define TBEGIN (-20.0)
#define TEND (20.0)
#define TSTEP (0.1)
#define EPS 0.5
#define H 1.054572e-27
#define PIM4 0.7511255444649425 
#define PI 3.14159265359
#define _GNU_SOURCE
#define OMEGA 1
#define F0 1
#define ALPHA 4.472135955
#define NPOINTS (1024.0)	//


double f(double x, void * params);
double A_evaluate(double t, void * params);
double A_integrate(double t1, double t);
double S_integrate_new(double tau, double t);
double S_integrate_func_new(double eps, void * params);
double Alpha2(double t1, double t);
double GetEvalSolve(double t);
double GetEvalSolve2(double t);
double GetEvalSolve3(double t);

double countOfRootsPrev = 0;


//Структура для численнного решения уравнения
struct quadratic_params
  {
    double t;
  };

struct S_integrate_params
  {
    double Aint;
  };

//Структура для численнного решения уравнения
double quadratic (double t0, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double t = p->t;

  return Alpha2(t0, t)-sqrt(2*EPS);
}

double quadratic2 (double t0, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double t = p->t;

  return Alpha2(t0, t)+sqrt(2*EPS);
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
	double h = 0.5;
	double t1_hi = t, t1_lo = t - h;
	int countOfRoots = 0;
	while( t1_lo > -10 )
	{
		if ((quadratic(t1_lo, &params)) * (quadratic(t1_hi, &params)) < 0)
		{
			gsl_root_fsolver_set (s, &AlphaFunc, t1_lo, t1_hi);			
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
			countOfRoots++;
			t1_hi = t1_lo;
			t1_lo -= h;
		}
		t1_hi = t1_lo;
		t1_lo -= h;
	}
	
	
	AlphaFunc.function = &quadratic2;
	AlphaFunc.params = &params;
	h = 0.5;
	t1_hi = t, t1_lo = t - h;
	while( t1_lo > -10 )
	{
		if ((quadratic2(t1_lo, &params)) * (quadratic2(t1_hi, &params)) < 0)
		{
			gsl_root_fsolver_set (s, &AlphaFunc, t1_lo, t1_hi);			
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
			countOfRoots++;
			t1_hi = t1_lo;
			t1_lo -= h;
		}
		t1_hi = t1_lo;
		t1_lo -= h;
	}
	gsl_root_fsolver_free (s);
	return countOfRoots;
}



int main()
{
	FILE *f44;
	f44 = fopen("f.txt", "w");
	FILE *f1;
	f1 = fopen("a.txt", "w");
	FILE *f2;
	f2 = fopen("roots.txt", "w");
	double tBegin = TBEGIN, tEnd = TEND, tStep = TSTEP;
	int steps = (int)(tEnd-tBegin)/tStep + 1;
	double* m1 = (double*) malloc(sizeof(double) * steps);
	double* m2 = (double*) malloc(sizeof(double) * steps);
	
	//#pragma omp parallel for shared (m1, m2)
	for(int i = 0; i<steps; i++)
	{				
		//fprintf(f2, "%2.15f %2.15f\n", tBegin + i*tStep, GetEvalSolve2(tBegin + i*tStep)
		double t = tBegin + i*tStep;
		fprintf(f44, "%f %2.15f\n" ,t, f(t, 0));
		fprintf(f1, "%f %2.15f\n", t, A_evaluate(t, 0));
		/*printf("T = %f\n", t);
		double roots[20];
		int countOfRoots = FindRoots( tBegin + i*tStep , roots );
		
		for(int j = 0; j<countOfRoots; j++)
		{
			fprintf(f2, "%3.10f %3.10f\n", t, roots[j]);
		}
		
		double m1 = GetAnalyticSolve( t , countOfRoots, roots );		
		double m2 = GetEvalSolve3( t );		
		printf("M_analytic(eps, t) = %2.15f\n", m1);
		printf("M_numerical(eps, t) = %2.15f\n", m2);
		fprintf(f, "%2.15f %2.15f\n", t , m1);
		fprintf(f1, "%2.15f %2.15f\n", t, m2);*/
	}
	
	fclose(f44);
	fclose(f1);
	fclose(f2);	
	return 0;
}

double GetAnalyticSolve(double t, int countOfRoots, double* roots)
{
	complex sum = 0;

	int a = 0;
	
	for(int i = 0; i<countOfRoots; i++)
	{	
		double S_1 = S_integrate_new(roots[i], t) + EPS*(t-roots[i]);
		double d = Alpha2(roots[i], t) * (( f(roots[i], 0) ) -
			( 1/( (roots[i] - t)*(roots[i] - t) ) * 
			A_integrate(roots[i], t) + A_evaluate(roots[i], 0) / (t - roots[i])));
		sum += cexp(I * S_1)/(csqrt(d) * cpow( t - roots[i] , 3.0/2.0));	
	}	
	sum *= csqrt(2);
	return creal(sum)*creal(sum) + cimag(sum)*cimag(sum);
}

double f (double x, void * params) 
{
  return -F0*exp(-(x*x)/(ALPHA*ALPHA))*cos(OMEGA*x);
}

//A как интеграл от f
double A_evaluate(double t, void * params)
{
	double t1 = -50; 		//Типа -\infty
	double result;
	double error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &f;
	F.params = &alpha;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000000);
	gsl_integration_qags (&F, t1, t, 1e-7, 0, 1000000, w, &result, &error); 
	gsl_integration_workspace_free (w);
	return result;
}

//А как интеграл от A_eval
double A_integrate(double t1, double t)
{
	double result;
	double error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &A_evaluate;
	F.params = &alpha;
	gsl_integration_workspace * a = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 1e-7, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return result;
	
}


//Под интегралом функция S_integrate_func
double S_integrate_new(double tau, double t)
{
	if (tau >= t)
	{
		return 0;
	}
	double Aint = A_integrate(tau, t)/(t-tau);
	struct S_integrate_params params = { Aint };
	double result, error;
	gsl_function F;
	F.function = &S_integrate_func_new;
	F.params = &params;
	gsl_integration_workspace * b = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, tau, t, 1e-7, 0, 10000, b, &result, &error); 
	gsl_integration_workspace_free (b);
	return -result/2.0;
}

double S_integrate_func_new(double eps, void * params)
{	
	struct S_integrate_params *p 
    = (struct S_integrate_params *) params;

	double Aint = p->Aint;
	double temp = A_evaluate(eps, 0) - Aint;
	return temp*temp;
}


//Функця по которой решается уравнение
double Alpha2(double t1, double t)
{
	return A_evaluate(t1, 0) - A_integrate(t1, t) / (t - t1);
}

double S_integrate_func(double eps, void * params)
{
	double temp = A_evaluate(eps, 0);
	return temp*temp;
}

double GetEvalSolve3(double t)
{    
	double tempEnd = 350;
	double step = tempEnd / NPOINTS;
	double tempPower;	
	int i = 0;	
	complex sum = 0;
	
	for(int i = 1; i < (int)NPOINTS; i++)
	{
		sum += cexp(I * EPS * i*step) 
			/ cpow(i*step, 3.0/2.0) 
			* (cexp(I * S_integrate_new(t-i*step, t)) - 1);
	}

	sum *= step;
	sum *=1/csqrt(2*PI*I);
	double result = creal(sum)*creal(sum) + cimag(sum)*cimag(sum);
	return result;
}

double GetEvalSolve(double t)
{
	fftw_complex *in, *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    
	double tempEnd = 200;
	double step = tempEnd / NPOINTS;
	double tempPower;	
	int i = 0;	
	
	for(double tempT = 0; tempT < tempEnd; tempT+=step)
	{
		if(!tempT)
		{
			in[i] = 0;
			i++;
			continue;
		}
		in[i] = cexp(I * EPS * tempT) 
			/ cpow(tempT, 3.0/2.0) 
			* (cexp(I * S_integrate_new(t-tempT, t)) - 1);
		i++;	
		
	}
	
	p = fftw_plan_dft_1d(NPOINTS, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);    
	out[0]*=1.0/NPOINTS;
	out[0]*=tempEnd;
	out[0]*=1/csqrt(2*PI*I);
	double result = creal(out[0])*creal(out[0]) + cimag(out[0])*cimag(out[0]);
	free(out);
	fftw_destroy_plan(p);
	return result;
}

double GetEvalSolveS_0(double t)
{
	fftw_complex *in, *out;
    fftw_plan p;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * NPOINTS);
    
	double tempEnd = 90;
	double step = tempEnd / NPOINTS;
	double tempPower;	
		
	int i = 0;	
	printf("step = %f\n", step);
	for(double tempT = 0; tempT < tempEnd; tempT+=step)
	{
		if(!tempT)
		{
			in[i] = 0;
			i++;
			continue;
		}
		double S_0=-f(t,0)*f(t,0)/24*tempT*tempT*tempT;
		printf("double S_0 = %f\n", S_0);
		in[i] = cexp(I * EPS * tempT) 
			/ cpow(tempT, 3.0/2.0) 
			* (cexp(I * S_0) - 1);
		i++;			
	}
	
	p = fftw_plan_dft_1d(NPOINTS, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);    
	out[0]*=1.0/NPOINTS;
	out[0]*=tempEnd;
	double result = creal(out[0])*creal(out[0]) + cimag(out[0])*cimag(out[0]);
	free(out);
	fftw_destroy_plan(p);
	return result;
}
