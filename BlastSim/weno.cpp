
#include "weno.h"

#define F(i) f[D + i]


double WENO::beta0(double *f)
{
	//return POW2(F(-2) - 4 * F(-1) + 3 * F(0)) / 4 + POW2(F(-2) - 2 * F(-1) + F(0));
	return 13 / 12 * POW2(F(-2) - 2 * F(-1) + F(0)) + POW2(F(-2) - 4 * F(-1) + 3 * F(0)) / 4;
}

double WENO::beta1(double *f)
{
	//return POW2(F(-1) - 2 * F(1)) + POW2(F(-1) - 2 * F(0) + F(-1));
	return 13 / 12 * POW2(F(-1) - 2 * F(0) + F(1)) + POW2(F(-1) - F(1)) / 4;
}

double WENO::beta2(double *f)
{
	//return POW2(3 * F(0) - 4 * F(1) + F(2)) / 4 + POW2(F(0) - 2 * F(1) + F(2));
	return 13 / 12 * POW2(F(0) - 2 * F(1) + F(2)) + POW2(3 * F(0) - 4 * F(1) + 3 * F(2)) / 4;
}

double WENO::flux(double *f, bool b)
{
	double beta[] = {
		WENO::beta0(f),
		WENO::beta1(f),
		WENO::beta2(f) };

	double l = abs(0.5 * F(-2) - 2 * F(-1) + 1.5 * F(0)) - abs(-1.5 * F(0) + 2 * F(1) - 0.5 * F(2));
	double r = F(-2) - F(-1) - 3 * F(0) + 5 * F(1) - 2 * F(2);
	//double tau = l * l + r * r;
	double tau = abs(beta[0] - beta[2]);

	//double a[] = {
	//	0.1 * (1 + POW2(tau / (EPS + beta[0]))),
	//	0.6 * (1 + POW2(tau / (EPS + beta[1]))),
	//	0.3 * (1 + POW2(tau / (EPS + beta[2]))) };
	
	double d[] = { 0.1, 0.6, 0.3 };
	if (b) d[2] = 0.1, d[0] = 0.3;

	double a[] = {
		d[0] / (POW2((EPS + beta[0]))),
		d[1] / (POW2((EPS + beta[1]))),
		d[2] / (POW2((EPS + beta[2]))) };

	double sum = a[0] + a[1] + a[2];

	double w[] = {
		a[0] / sum,
		a[1] / sum,
		a[2] / sum };

	double R[] = {
		(2 * F(-2) - 7 * F(-1) + 11 * F(0)) / 6,
		(-F(-1) + 5 * F(0) + 2 * F(1)) / 6,
		(2 * F(0) + 5 * F(1) - F(2)) / 6 };

	double flux = 0;
	for (size_t n = 0; n < 3; n++)
	{
		flux += w[n] * R[n];
	}

	return flux;
}

double WENO::grad(double f[6], double df[6], double u[6])
{
	double g[] = { f[5], f[4], f[3], f[2], f[1], f[0] };
	double v[] = { u[5], u[4], u[3], u[2], u[1], u[0] };
	double a, left, right;

	for (double *iter = df; iter < df + 6; iter++) *iter = abs(*iter);
	a = *std::max_element(df, df + 6);

	//double roe = (f[D + 1] - f[D]) / (u[D + 1] - u[D]);
	
	// if the flux is from left to right
	if (df[D] > 0)
	{
		right = WENO::flux(f, false);
		left = WENO::flux(g-1, true);
	}
	// flux if from right to left
	else
	{
		right = WENO::flux(g, true);
		left = WENO::flux(f-1, false);
	}
	
	//right = (WENO::flux(f, false) + WENO::flux(g, false) - a * (WENO::flux(v, false) - WENO::flux(u, false))) / 2;
	//left = (WENO::flux(f - 1, true) + WENO::flux(g - 1, true) - a * (WENO::flux(v - 1, true) - WENO::flux(u - 1, true))) / 2;
	
	return -(right - left) / H;
}

