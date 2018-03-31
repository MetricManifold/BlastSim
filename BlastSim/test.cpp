
#include "test.h"

void testCurveMu()
{
	FILE *f;
	f = fopen("results_curveMu.txt", "w");

	double
		e0 = 0.6,
		en = 3.3,
		de = (en - e0) / 1000;

	for (auto &rho : {
		RHO0 * 1e1, RHO0 * 1e0,
		RHO0 * 1e-1, RHO0 * 1e-2,
		RHO0 * 1e-3, RHO0 * 1e-4,
		RHO0 * 1e-5 })
	{

		fprintf(f, "# %f\n", rho / RHO0);
		for (double loge = e0; loge < en; loge += de)
		{
			fprintf(f, "%f %f\n", loge, FIT::mu(pow(10, loge) * E0, rho) / MU0);
		}
		fprintf(f, "\n\n");
	}

	fclose(f);
	exit(0);
}


void testCurveP()
{
	FILE *f;
	f = fopen("results_curveP.txt", "w");

	double
		e0 = 0.6,
		en = 3.6,
		de = (en - e0) / 1000;

	for (auto &rho : {
		RHO0 * 1e-7, RHO0 * 1e-6,
		RHO0 * 1e-5, RHO0 * 1e-4,
		RHO0 * 1e-3, RHO0 * 1e-2,
		RHO0 * 1e-1, RHO0 * 1e0,
		RHO0 * 1e1, RHO0 * 1e2,
		RHO0 * 1e3 })
	{

		fprintf(f, "# %f\n", rho / RHO0);
		for (double loge = e0; loge < en; loge += de)
		{
			fprintf(f, "%f %f\n", log10(FIT::p(pow(10, loge) * E0, rho) / P0), loge);
		}
		fprintf(f, "\n\n");
	}

	fclose(f);
	exit(0);
}

void testCurveT()
{
	FILE *f;
	f = fopen("results_curveT.txt", "w");

	for (auto pack : std::initializer_list<double[3]>{
		{RHO0 * 1e-7, -6.8, -4.5}, {RHO0 * 1e-6, -5.75, -3.5},
		{RHO0 * 1e-5, -4.8, -3}, {RHO0 * 1e-4, -3.8, -1.75},
		{RHO0 * 1e-3, -2.8, -0.75}, {RHO0 * 1e-2, -1.9, 0.5},
		{RHO0 * 1e-1, -0.9, 1.5}, {RHO0 * 1e0, 0.1, 2.0},
		{RHO0 * 1e1, 1.1, 3.0 }, {RHO0 * 1e2, 2.1, 4.0},
		{RHO0 * 1e3, 3.2, 5.0 } })
	{
		double rho = pack[0];
		double p0 = pack[1];
		double pn = pack[2];
		double pe = (pn - p0) / 100;

		fprintf(f, "# %f\n", rho / RHO0);
		for (double logp = p0; logp < pn; logp += pe)
		{
			fprintf(f, "%f %f\n", logp, log10(FIT::T(pow(10, logp) * P0, rho) / T0));
		}
		fprintf(f, "\n\n");
	}

	fclose(f);
	exit(0);
}
