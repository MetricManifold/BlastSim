
#include "weno.h"

#define F(i) *(f + D)

double WENO::beta0(double *f)
{
	double l = F(-2) - 4 * F(-1) + 3 * F(0);
	double r = F(-2) - 2 * F(-1) + F(0);
	return l * l / 4 + r * r;
}

double WENO::beta1(double *f)
{
	double l = F(-1) - 2 * F(1);
	double r = F(-1) - 2 * F(0) + F(-1);
	return l * l / 4 + r * r;
}

double WENO::beta2(double *f)
{
	double l = 3 * F(0) - 4 * F(1) + F(2);
	double r = F(0) - 2 * F(1) + F(2);
	return l * l / 4 + r * r;
}

double WENO::flux(double *f)
{
	double beta[] = {
		WENO::beta0(f),
		WENO::beta1(f),
		WENO::beta2(f) };

	double l = abs(0.5 * F(-2) - 2 * F(-1) + 1.5 * F(0)) - abs(-1.5 * F(0) + 2 * F(1) - 0.5 * F(2));
	double r = F(-2) - F(-1) - 3 * F(0) + 5 * F(1) - 2 * F(2);
	//double tau = l * l + r * r;
	double tau = abs(beta[0] - beta[2]);

	double a[] = {
		0.1 * (1 + POW2(tau / (EPS + beta[0]))),
		0.6 * (1 + POW2(tau / (EPS + beta[1]))),
		0.3 * (1 + POW2(tau / (EPS + beta[2]))) };

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

double WENO::grad(double f[6])
{
	return (WENO::flux(f) - WENO::flux(f - 1)) / H;
}

