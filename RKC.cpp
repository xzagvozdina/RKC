#include <iostream>
#include <fstream>
//#include <cmath>
#include <math.h>
#include <iomanip>
#include <float.h>
//#include <stdafx.h>
#include <cstdlib>
#include <string.h>   
#include <stdlib.h>
#include <sstream>
//#include <vector>

using namespace std;

typedef double(*function)(double x);
double T0(double x);
double T1(double x);
double T2(double x);
double T3(double x);
double T4(double x);
double T5(double x);
double dT0(double x);
double dT1(double x);
double dT2(double x);
double dT3(double x);
double dT4(double x);
double dT5(double x);
double ddT0(double x);
double ddT1(double x);
double ddT2(double x);
double ddT3(double x);
double ddT4(double x);
double ddT5(double x);
double funcm(int fl, double y1, double y2);
double exactm1(double x);
double exactm2(double x);
double maxerror(double *y, double *x, function f, int N);
void rkc3m(double *y1, double *y2, double *x, int N, double h);
void rkc5m(double *y1, double *y2, double *x, int N, double h);

int main()
{
	setlocale(LC_ALL,"Russian");

	int n;
	cout << "Количество точек:  ";
	cin >> n;
	//int l = 2 * n;
	cin.get();

	double a = 0.;   //
	double b = 2. * 6.2831;   //
	double *y1 = new double[n + 1];
	double *y2 = new double[n + 1];
	/*y1[0] = 0.;       //1.;
	y2[0] = 1.;       //
	double *yy1 = new double[l + 1];
	double *yy2 = new double[l + 1];
	yy1[0] = 0.;       //1.;
	yy2[0] = 1.;       //

	double h = (b - a) / n;
	double hh = h / 2;

	double *x = new double[n + 1];
	for (int i = 0; i <= n; ++i)
	{
		x[i] = a + i*h;
	}

	double *xx = new double[l + 1];
	for (int i = 0; i <= l; ++i)
	{
		xx[i] = a + i * hh;
	}*/

	ofstream errorm;
	errorm.open("errorm.dat");

	ofstream dotsm1;
	dotsm1.open("dots1m.dat");

	ofstream dotsm2;
	dotsm2.open("dots2m.dat");

	ofstream faz;
	faz.open("faz.dat");

	ofstream en;
	en.open("en.dat");

	ofstream errl;
	errl.open("errl.dat");

	for (int k = 0; k <= 3; ++k)
	{
		int l = pow(2, k) * n;
		double h = (b - a) / l;
		double *y1 = new double[l + 1];
		double *y2 = new double[l + 1];
		y1[0] = 1.;       //1.;
		y2[0] = 0.;
		
		en << h << '\n';

		double *x = new double[l + 1];
		for (int i = 0; i <= l; ++i)
		{
			x[i] = a + i*h;
		}

		dotsm1 << h << '\n';
		dotsm2 << h << '\n';
		faz << h << '\n';

		rkc5m(y1, y2, x, l, h);
		for (int i = 0; i <= l; ++i)
		{
			dotsm1 << x[i] << " " << y1[i] << '\n';
			dotsm2 << x[i] << " " << y2[i] << '\n';
			faz << y1[i] << " " << y2[i] << '\n';
		}

		double energy = pow(y2[0], 2) + pow(y1[0], 2) / 4.;
		double erren = pow(y2[1], 2) + pow(y1[1], 2) / 4. - energy;
		for (int i = 2; i <= l; ++i)
		{
			double erren1 = pow(y2[i], 2) + pow(y1[i], 2) / 4. - energy;
			if (erren1 > erren)
				erren = erren1;
		}

		double erl = 0.;
		for (int i = 0; i <= l; ++i)
			erl += h*fabs(exactm1(x[i])-y1[i]);

		en << erren << '\n';
		errl << erl << '\n';
		dotsm1 << '\n';
		dotsm2 << '\n';
		faz << '\n';
		en << '\n';

		errorm << h << " " << maxerror(y1, x, exactm1, n) << " " << maxerror(y2, x, exactm2, n) << '\n';
	}

	errorm.close();
	en.close();

	/*rkc3m(y1, y2, x, n, h);
	rkc3m(yy1, yy2, xx, l, hh);

	ofstream dotsm;
	dotsm.open("dotsm.dat");
	if (dotsm.is_open())
	{
		for (int i = 0; i <= n; ++i)
			dotsm << x[i] << " " << y1[i] << " " << y2[i] << '\n';
		cout << "dotsm.dat" << endl;
	}
	else
	{
		cout << "File Error" << endl;
	}
	dotsm.close();

	ofstream ddotsm;
	ddotsm.open("dots1m.dat");
	if (ddotsm.is_open())
	{
		for (int i = 0; i <= l; ++i)
			ddotsm << xx[i] << " " << yy1[i] << " " << yy2[i] << '\n';
		cout << "dots1m.dat" << endl;
	}
	else
	{
		cout << "File Error" << endl;
	}
	ddotsm.close();

	ofstream errorm;
	errorm.open("errorm.dat");
	if (errorm.is_open())
	{
		errorm << h << " " << maxerror(y1, x, exactm1, n) << " " << maxerror(y2, x, exactm2, n) << '\n';
		errorm << hh << " " << maxerror(yy1, xx, exactm1, l) << " " << maxerror(yy2, xx, exactm2, l) << '\n';
		cout << "errorm.dat" << endl;
	}
	else
	{
		cout << "File Error" << endl;
	}
	errorm.close();

	ofstream fazm;
	fazm.open("fazm.dat");
	if (fazm.is_open())
	{
		for (int i = 0; i <= n; ++i)
			fazm << y1[i] << " " << y2[i] << '\n';
		cout << "fazm.dat" << endl;
	}
	else
	{
		cout << "File Error" << endl;
	}
	fazm.close();

	ofstream fazm1;
	fazm1.open("fazm1.dat");
	if (fazm1.is_open())
	{
		for (int i = 0; i <= l; ++i)
			fazm1 << yy1[i] << " " << yy2[i] << '\n';
		cout << "fazm1.dat" << endl;
	}
	else
	{
		cout << "File Error" << endl;
	}
	fazm1.close();*/

	cin.get();

  	return 0;
}

double T0(double x)
{
	return 1.;
}

double T1(double x)
{
	return x;
}

double T2(double x)
{
	return 2 * x*x - 1;
}

double T3(double x)
{
	return 4 * x*x*x - 3 * x;
}

double T4(double x)
{
	return 1. - 8.*x*x + 8.*x*x*x*x;
}

double T5(double x)
{
	return 5.*x - 20.*x*x*x + 16.*x*x*x*x*x;
}

double dT0(double x)
{
	return 0.;
}

double dT1(double x)
{
	return 1.;
}

double dT2(double x)
{
	return 4.*x;
}

double dT3(double x)
{
	return 12.*x*x - 3.;
}

double dT4(double x)
{
	return -16.*x + 32.*x*x*x;
}

double dT5(double x)
{
	return 5. - 60.*x*x + 80.*x*x*x*x;
}

double ddT0(double x)
{
	return 0.;
}

double ddT1(double x)
{
	return 0.;
}

double ddT2(double x)
{
	return 4.;
}

double ddT3(double x)
{
	return 24.*x;
}

double ddT4(double x)
{
	return -16. + 96.*x*x;
}

double ddT5(double x)
{
	return -120.*x + 320.*x*x*x;
}

double funcm(int fl, double y1, double y2) //double func(double y)
{
	if (fl == 1)
		return y2;
	else
		return -1. / 4. * y1;
}

double exactm1(double x)
{
	//return sin(x);
	//return 2 * sin(x/2.);  //cos(x) + sin(x);
	return cos(x / 2.);
}

double exactm2(double x)
{
	//return cos(x);
	//return cos(x/2.);  //cos(x) - sin(x);
	return -1 / 2. * sin(x / 2.);
}

double maxerror(double *y, double *x, function f, int N)
{
	double *u = new double[N + 1];
	double *err = new double[N + 1];

	for (int i = 0; i <= N; ++i)
	{
		u[i] = f(x[i]);
		err[i] = fabs(y[i] - u[i]);
	}


	double max = err[0];

	for (int i = 0; i <= N; ++i)
	{
		if (max < err[i])
			max = err[i];
	}

	return max;
}

void rkc3m(double *y1, double *y2, double *x, int N, double h)
{
	int k = 4;
	double eps = 0.15;

	double w0 = 1 + eps / ((k - 1)*(k - 1));
	double w1 = dT3(w0) / ddT3(w0);

	double *b = new double[k];
	b[0] = b[1] = b[2] = ddT2(w0) / (dT2(w0)*dT2(w0));
	b[3] = ddT3(w0) / (dT3(w0)*dT3(w0));

	double *a = new double[k - 1];
	a[0] = 0;
	a[1] = 1 - b[1] * T1(w0);
	a[2] = 1 - b[2] * T2(w0);

	double *mu = new double[k];
	mu[0] = 0;
	mu[1] = 0.;
	for (int i = 2; i <= k - 1; ++i)
		mu[i] = 2. * w0 * b[i] / b[i - 1];

	double *nu = new double[k];
	nu[0] = 0;
	nu[1] = 0;
	for (int i = 2; i <= k - 1; ++i)
		nu[i] = -b[i] / b[i - 2];

	double *c = new double[k];
	c[0] = 0;
	c[1] = w1 * b[1] * dT1(w0);
	c[2] = w1 * b[2] * dT2(w0);
	c[3] = 1;

	double *kappa = new double[k];
	kappa[0] = 0;
	kappa[1] = c[1];
	for (int i = 2; i <= k - 1; ++i)
		kappa[i] = 2 * w1*b[i] / b[i - 1];

	double *g1 = new double[k];
	double *g2 = new double[k];

	for (int i = 0; i < N; ++i)
	{
		g1[0] = y1[i];
		g2[0] = y2[i];
		g1[1] = y1[i] + kappa[1] * h * funcm(1, g1[0], g2[0]); 
		g2[1] = y2[i] + kappa[1] * h * funcm(2, g1[0], g2[0]);
		g1[2] = y1[i] + mu[2] * (g1[1] - y1[i]) + kappa[2] * h * funcm(1, g1[1], g2[1]) - a[1] * kappa[2] * h * funcm(1, g1[0], g2[0]);
		g2[2] = y2[i] + mu[2] * (g2[1] - y2[i]) + kappa[2] * h * funcm(2, g1[1], g2[1]) - a[1] * kappa[2] * h * funcm(2, g1[0], g2[0]);
		g1[3] = y1[i] + mu[3] * (g1[2] - y1[i]) + nu[3] * (g1[1] - y1[i]) + kappa[3] * h * funcm(1, g1[2], g2[2]) - a[2] * kappa[3] * h * funcm(1, g1[0], g2[0]);
		g2[3] = y2[i] + mu[3] * (g2[2] - y2[i]) + nu[3] * (g2[1] - y2[i]) + kappa[3] * h * funcm(2, g1[2], g2[2]) - a[2] * kappa[3] * h * funcm(2, g1[0], g2[0]);
		y1[i + 1] = g1[3];
		y2[i + 1] = g2[3];
	}

	cin.get();
}

void rkc5m(double *y1, double *y2, double *x, int N, double h)
{
	int k = 6;
	double eps = 0.15;

	double w0 = 1 + eps / ((k - 1)*(k - 1));
	double w1 = dT5(w0) / ddT5(w0);

	double *b = new double[k];
	b[0] = b[1] = b[2] = ddT2(w0) / (dT2(w0)*dT2(w0));
	b[3] = ddT3(w0) / (dT3(w0)*dT3(w0));
	b[4] = ddT4(w0) / (dT4(w0)*dT4(w0));
	b[5] = ddT5(w0) / (dT5(w0)*dT5(w0));

	double *a = new double[k];
	a[0] = 0;
	a[1] = 1 - b[1] * T1(w0);
	a[2] = 1 - b[2] * T2(w0);
	a[3] = 1 - b[3] * T3(w0);
	a[4] = 1 - b[4] * T4(w0);
	a[5] = 1 - b[5] * T5(w0);

	double *mu = new double[k];
	mu[0] = 0;
	mu[1] = 0.;
	for (int i = 2; i <= k - 1; ++i)
		mu[i] = 2. * w0 * b[i] / b[i - 1];

	double *nu = new double[k];
	nu[0] = 0;
	nu[1] = 0;
	for (int i = 2; i <= k - 1; ++i)
		nu[i] = -b[i] / b[i - 2];

	double *c = new double[k];
	c[0] = 0;
	c[1] = w1 * b[1] * dT1(w0);
	c[2] = w1 * b[2] * dT2(w0);
	c[3] = w1 * b[3] * dT3(w0);
	c[4] = w1 * b[4] * dT4(w0);
	c[5] = 1.;

	double *kappa = new double[k];
	kappa[0] = 0;
	kappa[1] = c[1];
	for (int i = 2; i <= k - 1; ++i)
		kappa[i] = 2 * w1*b[i] / b[i - 1];

	double *g1 = new double[k];
	double *g2 = new double[k];

	for (int i = 0; i < N; ++i)
	{
		g1[0] = y1[i];
		g2[0] = y2[i];
		g1[1] = y1[i] + kappa[1] * h * funcm(1, g1[0], g2[0]);
		g2[1] = y2[i] + kappa[1] * h * funcm(2, g1[0], g2[0]);
		g1[2] = y1[i] + mu[2] * (g1[1] - y1[i]) + kappa[2] * h * funcm(1, g1[1], g2[1]) - a[1] * kappa[2] * h * funcm(1, g1[0], g2[0]);
		g2[2] = y2[i] + mu[2] * (g2[1] - y2[i]) + kappa[2] * h * funcm(2, g1[1], g2[1]) - a[1] * kappa[2] * h * funcm(2, g1[0], g2[0]);
		g1[3] = y1[i] + mu[3] * (g1[2] - y1[i]) + nu[3] * (g1[1] - y1[i]) + kappa[3] * h * funcm(1, g1[2], g2[2]) - a[2] * kappa[3] * h * funcm(1, g1[0], g2[0]);
		g2[3] = y2[i] + mu[3] * (g2[2] - y2[i]) + nu[3] * (g2[1] - y2[i]) + kappa[3] * h * funcm(2, g1[2], g2[2]) - a[2] * kappa[3] * h * funcm(2, g1[0], g2[0]);
		g1[4] = y1[i] + mu[4] * (g1[3] - y1[i]) + nu[4] * (g1[2] - y1[i]) + kappa[4] * h * funcm(1, g1[3], g2[3]) - a[3] * kappa[4] * h * funcm(1, g1[0], g2[0]);
		g2[4] = y2[i] + mu[4] * (g2[3] - y2[i]) + nu[4] * (g2[2] - y2[i]) + kappa[4] * h * funcm(2, g1[3], g2[3]) - a[3] * kappa[4] * h * funcm(2, g1[0], g2[0]);
		g1[5] = y1[i] + mu[5] * (g1[4] - y1[i]) + nu[5] * (g1[3] - y1[i]) + kappa[5] * h * funcm(1, g1[4], g2[4]) - a[4] * kappa[5] * h * funcm(1, g1[0], g2[0]);
		g2[5] = y2[i] + mu[5] * (g2[4] - y2[i]) + nu[5] * (g2[3] - y2[i]) + kappa[5] * h * funcm(2, g1[4], g2[4]) - a[4] * kappa[5] * h * funcm(2, g1[0], g2[0]);
		y1[i + 1] = g1[5];
		y2[i + 1] = g2[5];
	}

	cin.get();
}