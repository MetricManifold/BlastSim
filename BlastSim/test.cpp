
#include "test.h"

namespace TEST
{
	void test()
	{
		testCurveMu();
		testCurveP();
		testCurveT();
		testCurveK();
		testCurveA();

		exit(0);
	}

	void testCurveMu()
	{
		FILE *f;
		f = fopen("results_curveMu.txt", "w");

		double
			e0 = 0.6,
			en = 3.3,
			de = (en - e0) / 1000;

		for (auto &rho : {
			RHO_CURVE * 1e1, RHO_CURVE * 1e0,
			RHO_CURVE * 1e-1, RHO_CURVE * 1e-2,
			RHO_CURVE * 1e-3, RHO_CURVE * 1e-4,
			RHO_CURVE * 1e-5 })
		{

			fprintf(f, "# %f\n", rho / RHO_CURVE);
			for (double loge = e0; loge < en; loge += de)
			{
				fprintf(f, "%f %f\n", loge, FIT::mu(pow(10, loge), rho/ RHO_CURVE));
			}
			fprintf(f, "\n\n");
		}

		fclose(f);
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
			RHO_CURVE * 1e-7, RHO_CURVE * 1e-6,
			RHO_CURVE * 1e-5, RHO_CURVE * 1e-4,
			RHO_CURVE * 1e-3, RHO_CURVE * 1e-2,
			RHO_CURVE * 1e-1, RHO_CURVE * 1e0,
			RHO_CURVE * 1e1, RHO_CURVE * 1e2,
			RHO_CURVE * 1e3 })
		{

			fprintf(f, "# %f\n", rho / RHO_CURVE);
			for (double loge = e0; loge < en; loge += de)
			{
				fprintf(f, "%f %f\n", log10(FIT::p(pow(10, loge), rho / RHO_CURVE)), loge);
			}
			fprintf(f, "\n\n");
		}

		fclose(f);
	}

	void testCurveT()
	{
		FILE *f;
		f = fopen("results_curveT.txt", "w");

		for (auto pack : std::initializer_list<double[3]>{
			{ RHO_CURVE * 1e-7, -6.8, -4.5 },{ RHO_CURVE * 1e-6, -5.75, -3.5 },
			{ RHO_CURVE * 1e-5, -4.8, -3 },{ RHO_CURVE * 1e-4, -3.8, -1.75 },
			{ RHO_CURVE * 1e-3, -2.8, -0.75 },{ RHO_CURVE * 1e-2, -1.9, 0.5 },
			{ RHO_CURVE * 1e-1, -0.9, 1.5 },{ RHO_CURVE * 1e0, 0.1, 2.0 },
			{ RHO_CURVE * 1e1, 1.1, 3.0 },{ RHO_CURVE * 1e2, 2.1, 4.0 },
			{ RHO_CURVE * 1e3, 3.2, 5.0 } })
		{
			double rho = pack[0];
			double p0 = pack[1];
			double pn = pack[2];
			double pe = (pn - p0) / 100;

			fprintf(f, "# %f\n", rho / RHO_CURVE);
			for (double logp = p0; logp < pn; logp += pe)
			{
				fprintf(f, "%f %f\n", logp, log10(FIT::T(pow(10, logp), rho / RHO_CURVE)));
			}
			fprintf(f, "\n\n");
		}

		fclose(f);
	}

	void testCurveA()
	{
		FILE *f;
		f = fopen("results_curveA.txt", "w");

		double
			e0 = 0.6,
			en = 3.0,
			de = (en - e0) / 1000;

		for (auto pack : std::initializer_list<double[3]>{
			{ RHO_CURVE * 1e-7, -6.8, -4.5 },{ RHO_CURVE * 1e-6, -5.75, -3.5 },
			{ RHO_CURVE * 1e-5, -4.8, -3 },{ RHO_CURVE * 1e-4, -3.8, -1.75 },
			{ RHO_CURVE * 1e-3, -2.8, -0.75 },{ RHO_CURVE * 1e-2, -1.9, 0.5 },
			{ RHO_CURVE * 1e-1, -0.9, 1.5 },{ RHO_CURVE * 1e0, 0.1, 2.0 },
			{ RHO_CURVE * 1e1, 1.1, 3.0 },{ RHO_CURVE * 1e2, 2.1, 4.0 },
			{ RHO_CURVE * 1e3, 3.2, 5.0 } })
		{
			double rho = pack[0];
			double p0 = pack[1];
			double pn = pack[2];
			double pe = (pn - p0) / 100;

			//for (double logp = p0; logp < pn; logp += pe)

			fprintf(f, "# %f\n", rho / RHO_CURVE);
			for (double loge = e0; loge < en; loge += de)
			{
				double p = FIT::p(pow(10, loge), rho / RHO_CURVE);
				fprintf(f, "%f %f\n", log10(p), log10(FIT::a(pow(10, loge), rho / RHO_CURVE) / A_CURVE));
			}
			fprintf(f, "\n\n");
		}

		fclose(f);
	}

	void testCurveK()
	{
		FILE *f;
		f = fopen("results_curveK.txt", "w");

		double
			e0 = 0.6,
			en = 3.3,
			de = (en - e0) / 1000;

		for (auto &rho : {
			RHO_CURVE * 1e-5, RHO_CURVE * 1e-4,
			RHO_CURVE * 1e-3, RHO_CURVE * 1e-2,
			RHO_CURVE * 1e-1, RHO_CURVE * 1e0,
			RHO_CURVE * 1e1})
		{

			fprintf(f, "# %f\n", rho / RHO_CURVE);
			for (double loge = e0; loge < en; loge += de)
			{
				fprintf(f, "%f %f\n", loge, FIT::k(pow(10, loge), rho / RHO_CURVE));
			}
			fprintf(f, "\n\n");
		}

		fclose(f);
	}
}