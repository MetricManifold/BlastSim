
#include "weno.h"

#define F(i) f[i + D]

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

double WENO::flux(double *f, double w[3])
{
	double R[] = {
		(2 * F(-2) - 7 * F(-1) + 11 * F(0)) / 6,
		(-F(-1) + 5 * F(0) + 2 * F(1)) / 6,
		(2 * F(0) + 5 * F(1) - F(2)) / 6 };

	double flux = 0;
	for (size_t l = 0; l < 2; l++)
	{
		flux += w[l] * R[l];
	}

	return flux;
}

double WENO::weno(double f[6])
{
	double beta[] = {
		WENO::beta0(f + 1),
		WENO::beta1(f + 1),
		WENO::beta2(f + 1) };

	double l = abs(0.5 * F(-2) - 2 * F(-1) + 1.5 * F(0)) - abs(-1.5 * F(0) + 2 * F(1) - 0.5 * F(2));
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

	return (WENO::flux(f + 1, w) - WENO::flux(f, w)) / H;
}

