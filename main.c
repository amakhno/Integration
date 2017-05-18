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

double f(double x, void * params);
double A_evaluate(double t, void * params);
double A_integrate(double t1, double t);
double S_integrate(double t1, double t);
double S_integrate_func(double eps, void * params);
double S_integrate_new(double tau, double t);
double S_integrate_func_new(double eps, void * params);
double Alpha2(double t1, double t);
double GetEvalSolve(double t);


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
	double h = 0.2;
	double t1_hi = t, t1_lo = t - h;
	int countOfRoots = 0;
	while( t1_lo > -30 )
	{
		if ((Alpha2(t1_lo, t)) * (Alpha2(t1_hi, t)) < 0)
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
	printf("\nНа промежутке [%f, %f]\n", t1_lo, t);
	for(int i = 0; i<countOfRoots; i++)
	{
		printf("Корень %d: %f\n", i, roots[i]);
	}
	gsl_root_fsolver_free (s);
	return countOfRoots;
}

double GetAnalyticSolve(double t, int countOfRoots, double* roots)
{
	printf("\nПоиска аналитеческого решения\n");
	printf("\nВсего корней при t = %f:\t%d\n", t, countOfRoots);
	fftw_complex sum = 0;
	double E = -2;
	
	for(int i = 0; i<countOfRoots; i++)
	{		
		double S_1 = S_integrate(roots[i], t) + E*(t-roots[i]);
		double d = - E * ( f(roots[i], 0) ); //-
			//( 1/( (roots[i] - t)*(roots[i] - t) ) * 
			//A_integrate(roots[i], t) - A_evaluate(roots[i], 0) / (t - roots[i])) );
		sum += cexp(I * S_1)/(csqrt(d) * cpow( t - roots[i] , 3.0/2.0));	
		//printf("d = %f\n", d);	
		//printf("S_integrate(roots[i], t) = %f\n", creal(S_integrate(roots[i], t)));	
		//printf("E*(t-roots[i]) = %f\n", E*(t-roots[i]));
		//printf("cexp(I * S_integrate(roots[i], t) + E*(t-roots[i])) = %f\n", creal(cexp(I * S_integrate(roots[i], t) + E*(t-roots[i]))));
	}	
	printf("Value = %f\n", creal(sum)*creal(sum) + cimag(sum)*cimag(sum));
	return creal(sum)*creal(sum) + cimag(sum)*cimag(sum);
}

int main()
{
	printf("A(5) = %f\n", A_evaluate(5, 0));
	printf("A_integrate(5, 6) = %f\n", A_integrate(5, 6));
	printf("S(5, 6) = %f\n", S_integrate_new(5, 6));
	FILE *f;
	f = fopen("analytical.txt", "w");
	FILE *f1;
	f1 = fopen("numerical.txt", "w");
	FILE *f2;
	f2 = fopen("roots.txt", "w");
	
	//----------------------PrintS
	FILE *f3;
	f3 = fopen("S.txt", "w");	
	for(double t1 = -6; t1<5.9; t1+=0.01)
	{
		printf("t1 = %f\n", t1);
		fprintf(f3, "%f %f\n", t1 , S_integrate_new(t1, 6.0));
	}
	fclose(f3);
	//----------------------PrintS
	
	double roots[50];
	
	for(double t = 0; t<10; t+=0.1)
	{		
		int countOfRoots = FindRoots( t , roots );
		for(int j = 0; j<countOfRoots; j++)
		{
			fprintf(f2, "%3.10f %3.10f\n", t, roots[j]);
		}
		
		double m1 = GetAnalyticSolve( t , countOfRoots, roots );
		double m2 = GetEvalSolve( t );
		fprintf(f, "%f %f\n", t , m1);
		fprintf(f1, "%f %f\n", t, m2);
	}
	fclose(f);
	fclose(f1);
	fclose(f2);	
	return 0;
}

double f (double x, void * params) 
{
  return -F0*exp(-(x*x)/(ALPHA*ALPHA))*cos(OMEGA*x);
}


//A как интеграл от f
double A_evaluate(double t, void * params)
{
	double t1 = -150; 		//Типа -\infty
	double result;
	double error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &f;
	F.params = &alpha;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000000);
	gsl_integration_qags (&F, t1, t, 1e-13, 0, 1000000, w, &result, &error); 
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
	gsl_integration_qags (&F, t1, t, 1e-13, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return result;
	
}

//Под интегралом функция S_integrate_func
double S_integrate(double t1, double t)
{
	double a_integrate_result = A_integrate(t1, t)/( 2*(t-t1) );
	double result, error;
	gsl_function F;
	F.function = &S_integrate_func;
	F.params = &a_integrate_result;
	gsl_integration_workspace * b = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 1e-13, 0, 10000, b, &result, &error); 
	gsl_integration_workspace_free (b);
	return -result/2.0 + a_integrate_result;
}

//Под интегралом функция S_integrate_func
double S_integrate_new(double tau, double t)
{
	double Aint = A_integrate(tau, t)/(t-tau);
	struct S_integrate_params params = { Aint };
	double result, error;
	gsl_function F;
	F.function = &S_integrate_func_new;
	F.params = &params;
	gsl_integration_workspace * b = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, tau, t, 1e-13, 0, 10000, b, &result, &error); 
	gsl_integration_workspace_free (b);
	return -result/2.0;
}

double S_integrate_func_new(double eps, void * params)
{	
	struct S_integrate_params *p 
    = (struct S_integrate_params *) params;

	double Aint = p->Aint;
	double temp = A_evaluate(eps, 0) - Aint;
	/*(A(eps) - intA(tau, t)/(t-tau))*/
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
		if(!tempT)
		{
			in[i] = 0;
			i++;
			continue;
		}
		
		in[i] = cexp(I * 1 * tempT) 
			/ cpow(tempT, 3.0/2.0) 
			* (cexp(I * S_integrate_new(t-tempT, t)) - 1);
		i++;	
		
	}
	
	p = fftw_plan_dft_1d(NPOINTS, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	fftw_execute(p);    
	double result = creal(out[0])*creal(out[0]) + cimag(out[0])*cimag(out[0]);
	free(out);
	fftw_destroy_plan(p);
	return result;
}
