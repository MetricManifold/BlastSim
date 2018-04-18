
#include "weno.h"

#define F(i) f[D + i]

double WENO::beta0(double *f)
{
	//return POW2(3 * F(0) - 4 * F(1) + F(2)) / 4 + POW2(F(0) - 2 * F(1) + F(2));
	return (13.0 / 12.0) * POW2(F(0) - 2 * F(1) + F(2)) + POW2(3 * F(0) - 4 * F(1) + F(2)) / 4.0;
}

double WENO::beta1(double *f)
{
	//return POW2(F(-1) - F(1)) / 4 + POW2(F(-1) - 2 * F(0) + F(-1));
	return (13.0 / 12.0) * POW2(F(-1) - 2 * F(0) + F(1)) + POW2(F(-1) - F(1)) / 4.0;
}

double WENO::beta2(double *f)
{
	//return POW2(F(-2) - 4 * F(-1) + 3 * F(0)) / 4 + POW2(F(-2) - 2 * F(-1) + F(0));
	return (13.0 / 12.0) * POW2(F(-2) - 2 * F(-1) + F(0)) + POW2(F(-2) - 4 * F(-1) + 3 * F(0)) / 4.0;
}

double WENO::flux(double *f, double d[3])
{
	//double l = abs(0.5 * F(-2) - 2 * F(-1) + 1.5 * F(0)) - abs(-1.5 * F(0) + 2 * F(1) - 0.5 * F(2));
	//double r = F(-2) - F(-1) - 3 * F(0) + 5 * F(1) - 2 * F(2);
	//double tau = l * l + r * r;

	//double d[] = { 0.3, 0.6, 0.1 };
	double a[] = {
		d[0] / POW2(EPS + WENO::beta0(f)),
		d[1] / POW2(EPS + WENO::beta1(f)),
		d[2] / POW2(EPS + WENO::beta2(f)) };

	//double tau = abs(beta0(f) - beta2(f));
	//double a[] = {
	//	d[0] * (1 + POW2(tau / (EPS + WENO::beta0(f)))),
	//	d[1] * (1 + POW2(tau / (EPS + WENO::beta1(f)))),
	//	d[2] * (1 + POW2(tau / (EPS + WENO::beta2(f)))) };

	double sum = a[0] + a[1] + a[2];
	double w[] = {
		a[0] / sum,
		a[1] / sum,
		a[2] / sum };

	double R[] = {
		(2 * F(0) + 5 * F(1) - F(2)) / 6.0,
		(-F(-1) + 5 * F(0) + 2 * F(1)) / 6.0,
		(2 * F(-2) - 7 * F(-1) + 11 * F(0)) / 6.0 };

	double flux = 0;
	for (size_t n = 0; n < 3; n++)
	{
		flux += w[n] * R[n];
	}
	
	return flux;
}

double WENO::grad(double f[6], double df[6], double u[6])
{
	double g[] = { f[6], f[5], f[4], f[3], f[2], f[1], f[0] };
	double v[] = { u[6], u[5], u[4], u[3], u[2], u[1], u[0] };
	double left, right;

	//for (double *iter = df; iter < df + 6; iter++) *iter = abs(*iter);
	//double a = *std::max_element(df, df + 6);

	double a = *df;
	double d[] = { 0.3, 0.6, 0.1 };

	//// if the flux is from left to right
	//if (*df > 0)
	//{
	//	right = WENO::flux(f, d);
	//	left = WENO::flux(g, d);
	//}
	//// flux if from right to left
	//else
	//{
	//	right = WENO::flux(f - 1, d);
	//	left = WENO::flux(g - 1, d);
	//}

	right = (WENO::flux(f, d) + WENO::flux(g - 1, d) - a * (WENO::flux(v - 1, d) - WENO::flux(u, d))) / 2.0;
	left = (WENO::flux(f - 1, d) + WENO::flux(g, d) - a * (WENO::flux(v, d) - WENO::flux(u - 1, d))) / 2.0;

	//right = (WENO::flux(f) + WENO::flux(g - 1) - *df * (WENO::flux(v - 1) - WENO::flux(u))) / 2.0;
	//left = (WENO::flux(f - 1) + WENO::flux(g) - *df * (WENO::flux(v) - WENO::flux(u - 1))) / 2.0;

	return (right - left) / H;
}

