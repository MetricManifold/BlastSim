
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

	void testWENO(double(*fn)(double), double(*df)(double), bool(*cn)(double))
	{
		const size_t len = 200;
		double x0 = -1.0, xn = 1.0;
		double h = (xn - x0) / len, k = 0.5 * h;

		double u[len + D + E];
		double du[len];


		for (size_t i = 0; i < len; i++)
		{
			if (cn(i * h + x0))
			{
				//u[i + D] = -(i * h + x0);
				u[i + D] = -sin(PI * (i * h + x0));
				//u[i + D] = 1;
			}
			else
			{
				u[i + D] = 0;
			}

		}

		FILE *file;
		file = fopen("output.txt", "w");
		for (size_t i = 0; i < len; i++)
		{
			fprintf(file, "%f %f\n", i * h + x0, u[i + D]);
		}

		for (double time = 0; time < 0.55; time += k)
		{

			static auto boundary = [&]() {
				double a = 0;
				size_t index = 0;
				for (size_t i = D; i < len + D; i++)
				{
					double b = abs(df(u[i]));
					if (a < b)
					{
						a = b;
						index = i;
					}
				}
				//a = abs(df(u[index]));

				for (size_t i = 0; i < D; i++)
				{
					u[i] = u[i + len];
				}
				for (size_t i = len + D; i < D + len + E; i++)
				{
					u[i] = u[i - len];
				}

				return a;
			};

			double dk[4][len];
			static auto update = [&](size_t index, double dt) {
				double a = boundary();
				for (size_t i = 0; i < len; i++)
				{
					double fs[] = { fn(u[i]), fn(u[i + 1]), fn(u[i + 2]), fn(u[i + 3]), fn(u[i + 4]), fn(u[i + 5]), fn(u[i + 6]) };
					double ds[] = { df(u[i]), df(u[i + 1]), df(u[i + 2]), df(u[i + 3]), df(u[i + 4]), df(u[i + 5]), df(u[i + 6]) };
					double s[] = { u[i], u[i + 1], u[i + 2], u[i + 3], u[i + 4], u[i + 5], u[i + 6] };

					dk[index][i] = -WENO::grad(fs, &a, s) * H / h;
				}
				for (size_t i = 0; i < len; i++) u[i + D] = du[i] + dt * dk[index][i];
			};

			for (size_t i = 0; i < len; i++) du[i] = u[i + D];
			update(0, k / 2.0);
			update(1, k / 2.0);
			update(2, k);
			update(3, 0);
			for (size_t i = 0; i < len; i++) u[i + D] = du[i] + k / 6.0 * (dk[0][i] + 2.0 * dk[1][i] + 2.0 * dk[2][i] + dk[3][i]);

		}

		fprintf(file, "\n\n");
		for (size_t i = 0; i < len; i++)
		{
			fprintf(file, "%f %f\n", i * h + x0, u[i + D]);
		}
		fclose(file);

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
				fprintf(f, "%f %f\n", loge, FIT::mu(pow(10, loge), rho / RHO_CURVE));
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
			RHO_CURVE * 1e1 })
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