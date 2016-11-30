#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <float.h>
#include <gsl/gsl_integration.h>
#include <complex.h>


#define EPS 3.0e-14 
#define PIM4 0.7511255444649425 
#define MAXIT 10 
#define _GNU_SOURCE
#define OMEGA 1
#define ALPHA 4.472135955

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

int main()
{
	int n = 10;
	double *X = (double*)malloc((n+1)*sizeof(double));
	double *W = (double*)malloc((n+1)*sizeof(double));
	complex cmplx;
	
	
	X[0] = 0.013046735741414;
	X[1] = 0.067468316655508;
	X[2] = 0.160295215850488;
	X[3] = 0.283302302935376;
	X[4] = 0.425562830509184;
	X[5] = 0.574437169490816;
	X[6] = 0.716697697064624;
	X[7] = 0.839704784149512;
	X[8] = 0.932531683344492;
	X[9] = 0.986953264258586;
	W[0] = 0.033335672154341;
	W[1] = 0.074725674575290;
	W[2] = 0.109543181257991;
	W[3] = 0.134633359654996;
	W[4] = 0.147762112357376;
	W[5] = 0.147762112357376;
	W[6] = 0.134633359654996;
	W[7] = 0.109543181257991;
	W[8] = 0.074725674575290;
	W[9] = 0.033335672154341;	
	double a = -50;
	double b = 0;
	double A = IntegrateF(a, b, W, X, n, 12);
	printf ("My A evaluate   = % .18f\n", A);
	//---------------------------------------------Gsl
	//gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  
	double result, error;
	double expected = 0.026704689271297935868873867110;
	double alpha = 1.0;

	gsl_function F;
	F.function = &f;
	F.params = &alpha;
	result = A_evaluate(0, 0);
	
	//gsl_integration_qags (&F, -50, 0, 1e-14, 0, 1000,
    //                    w, &result, &error); 
    //---------------------------------------------Gsl-end
	
	printf ("\nresult          = % .18f\n", result);
	printf ("exact result    = % .18f\n", expected);
	printf ("estimated error = % .18f\n", error);
	printf ("actual error    = % .18f\n", result - expected);
	double a_integrate = A_integrate(-100, 0);
	printf ("A_integration = % .18f\n", a_integrate);
	/*double t1 = -10;
	double t = 10;
	double eps = 2.1;*/
	FILE *f;
	f = fopen("S.txt", "w");
	double t = 10;
	double tau = 0;
	double a_integrate_result;
	double temp = A_evaluate(5, (void *)0);
	//return (temp - a_integrate_result) * (temp - a_integrate_result);
	for(tau = 0; tau<30; tau+=0.1)
	{
		//a_integrate_result =  A_integrate(t-tau, t);
		//fprintf(f, "%f %f\n", tau, (temp - a_integrate_result) * (temp - a_integrate_result));
		fprintf(f, "%f %f\n", tau, S_integrate(t-tau, t));
	}
	fclose(f);
	//printf ("intervals       = %zu\n", w->size);

	//gsl_integration_workspace_free (w);
	return 0;
}

int gauleg(double x1, double x2, double *x, double *w, int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) 
	{
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do 
		{
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) 
			{
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n+1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i]=w[i];
	}
	return 0;
}

double F(double t)
{
	return exp(-(t*t)/(ALPHA*ALPHA))*cos(OMEGA*t);
}


double IntegrateF(double a, double b, double *w, double *x, double n, int N)
{
	double sum=0;
	
	double Size = (b-a)/N;
	
	for(int i = 0; i<N; i++)
	{
		sum += IntegrateP(a + Size*i, a + Size * (1+i), w, x, n);
	}
	
	return sum;
}

/*double IntegrateF2(double a, double b, double *w, double *x, double n, int N)
{
	double sum=0;
	
	double Size = (b-a)/N;
	
	for(int i = 0; i<N; i++)
	{
		sum += IntegrateP(a + Size*i, a + Size * (1+i), w, x, n);
	}
	
	double At = sum;
	
	
	
	return sum;
}*/

double IntegrateP(double a, double b, double *w, double *x, double n)
{
	double sum=0;
	
	double scale = (b-a);
	for(int i = 0; i<=n; i++)
	{		
		sum+=w[i]*F(a + scale * x[i])*scale;
	}

	return sum;
}

double f (double x, void * params) 
{
  return exp(-(x*x)/(ALPHA*ALPHA))*cos(OMEGA*x);
}

double A_evaluate(double t, void * params)
{
	double t1 = -100;
	double result;
	double error;
	double alpha = 1.0;
	gsl_function F;
	F.function = &f;
	F.params = &alpha;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	gsl_integration_qags (&F, t1, t, 0.000001, 0, 1000, w, &result, &error); 
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
	gsl_integration_qags (&F, t1, t, 0.00000001, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return result;
	
}

double S_integrate(double t1, double t)
{
	/*double param = 1;
	double result, error;
	double a_integrate_result = A_integrate(t1, t)/t1-t;
	gsl_function F;
	F.function = &S_integrate_func;
	F.params = &a_integrate_result;
	gsl_integration_workspace * a = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 0.001, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return -result;*/
	double a_integrate_result = A_integrate(t1, t)/( 2*(t1-t) );
	double result, error;
	gsl_function F;
	F.function = &S_integrate_func2;
	F.params = &a_integrate_result;
	gsl_integration_workspace * a = gsl_integration_workspace_alloc (10000);
	gsl_integration_qags (&F, t1, t, 0.001, 0, 10000, a, &result, &error); 
	gsl_integration_workspace_free (a);
	return -result/2.0 + a_integrate_result;
}

/*double S_integrate_func(double eps, void * params)
{
	double temp = A_evaluate(eps, (void *)0);
	return (temp - *(double *)params) * (temp - *(double *)params);
}*/

double S_integrate_func2(double eps, void * params)
{
	double temp = A_evaluate(eps, 0);
	return temp*temp;
}
