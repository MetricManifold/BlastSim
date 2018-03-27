
#include "weno.h"

#define F(i) f[i + 2]

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
	double R[] = {
		(2 * F(-2) - 7 * F(-1) + 11 * F(0)) / 6,
		(-F(-2) + 5 * F(-1) + 2 * F(0)) / 6,
		(2 * F(-2) + 5 * F(-1) - F(0)) / 6 };

	double beta[] = {
		WENO::beta0(f, i),
		WENO::beta1(f, i),
		WENO::beta2(f, i) };

	double l = abs(0.5 * F(-2) - 2 * F(-1) + 1.5 * F(0))
		- abs(-1.5 * F(0) + 2 * F(1) - 0.5 * F(2));
	double r = F(-2) - F(-1) - 3 * F(0) + 5 * F(1) - 2 * F(2);
	double tau = l * l + r * r;

	double a[] = {
		0.1 * (1 + POW2(tau / (EPS + beta[0]))),
		0.6 * (1 + POW2(tau / (EPS + beta[0]))),
		0.3 * (1 + POW2(tau / (EPS + beta[0]))) };

	double sum = a[0] + a[1] + a[2];

	double w[] = {
		a[0] / sum,
		a[1] / sum,
		a[2] / sum };

	double flux = 0;
	for (size_t l = 0; l < 2; l++)
	{
		flux += w[0] * R[0];
	}

	return flux;
}

double WENO::weno(double *f)
{
	return (WENO::flux(f, i) - WENO::flux(f, i - 1)) / H;
}

