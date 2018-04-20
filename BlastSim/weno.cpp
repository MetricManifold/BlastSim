
#include "weno.h"

// where 2 defines the center of the stencil [0 1 2 3 4]
#define F(i) f[2 + i]

double WENO::beta0(double *f)
{
	return POW2(3 * F(0) - 4 * F(1) + F(2)) / 4 + POW2(F(0) - 2 * F(1) + F(2));
}

double WENO::beta1(double *f)
{
	return POW2(F(-1) - F(1)) / 4 + POW2(F(-1) - 2 * F(0) + F(-1));
}

double WENO::beta2(double *f)
{
	return POW2(F(-2) - 4 * F(-1) + 3 * F(0)) / 4 + POW2(F(-2) - 2 * F(-1) + F(0));
}

double WENO::grad(double f[7], double u[7], double a)
{
	double g[] = { f[6], f[5], f[4], f[3], f[2], f[1], f[0] };
	double v[] = { u[6], u[5], u[4], u[3], u[2], u[1], u[0] };
	double left, right;

	right = (WENO::flux(f + 1) + WENO::flux(g) - a * (WENO::flux(v) - WENO::flux(u + 1))) / 2.0;
	left = (WENO::flux(f) + WENO::flux(g + 1) - a * (WENO::flux(v + 1) - WENO::flux(u))) / 2.0;

	return (right - left) / H;
}

double WENO::flux(double *f)
{
	return flux(f, f);
}

double WENO::flux(double *f, double *b)
{
	double tau = POW2(abs(0.5 * F(-2) - 2 * F(-1) + 1.5 * F(0)) -
		abs(-1.5 * F(0) + 2 * F(1) - 0.5 * F(2))) +
		POW2(F(-2) - F(-1) - 3 * F(0) + 5 * F(1) - 2 * F(2));

	double a[] = {
		0.3 * (1 + POW2(tau / (EPS + WENO::beta0(b)))),
		0.6 * (1 + POW2(tau / (EPS + WENO::beta1(b)))),
		0.1 * (1 + POW2(tau / (EPS + WENO::beta2(b)))) };

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

//double *WENO::gradV(double Ft[4][7], double Ut[4][7], double **L[2], double **R[2], double a)
//{
//	double G[7][4], V[7][4], F[7][4], U[7][4];
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 7; j++)
//		{
//			F[j][i] = Ft[i][j];
//			U[j][i] = Ut[i][j];
//		}
//	}
//	for (size_t i = 0; i < 4; i++)
//	{
//		for (size_t j = 0; j < 7; j++)
//		{
//			G[j][i] = F[6 - j][i];
//			V[j][i] = U[6 - j][i];
//		}
//	}
//
//	double left[4][4], right[4][4];
//	for (size_t s = 0; s < 4; s++)
//	{
//		double *l = L[0][s];
//		double *ll = L[1][s];
//		double f[5], ff[5], g[5], gg[5], u[5], uu[5], v[5], vv[5];
//		for (size_t i = 0; i < 5; i++)
//		{
//			f[i] = HELP::dot(ll, F[i + 1]);
//			ff[i] = HELP::dot(l, F[i]);
//
//			g[i] = HELP::dot(l, G[i + 1]);
//			gg[i] = HELP::dot(ll, G[i]);
//
//			u[i] = HELP::dot(ll, U[i + 1]);
//			uu[i] = HELP::dot(l, U[i]);
//
//			v[i] = HELP::dot(l, V[i + 1]);
//			vv[i] = HELP::dot(ll, V[i]);
//		}
//
//		left[s][0] = WENO::flux(ff);
//		left[s][1] = WENO::flux(g);
//		left[s][2] = WENO::flux(v);
//		left[s][3] = WENO::flux(uu);
//
//		right[s][0] = WENO::flux(f);
//		right[s][1] = WENO::flux(gg);
//		right[s][2] = WENO::flux(vv);
//		right[s][3] = WENO::flux(u);
//	}
//
//	double LEFT[4]{ 0 }, RIGHT[4]{ 0 };
//	for (size_t i = 0; i < 4; i++)
//	{
//		double l[4]{ 0 }, r[4]{ 0 };
//		for (size_t j = 0; j < 4; j++)
//		{
//			for (size_t s = 0; s < 4; s++)
//			{
//				l[j] += left[s][j] * R[0][s][i];
//				r[j] += right[s][j] * R[1][s][i];
//			}
//		}
//
//		LEFT[i] = (l[0] + l[1] - a * (l[2] - l[3])) / 2.0;
//		RIGHT[i] = (r[0] + r[1] - a * (r[2] - r[3])) / 2.0;
//	}
//
//	double *result = new double[4];
//	for (size_t i = 0; i < 4; i++)
//	{
//		result[i] = (RIGHT[i] - LEFT[i]) / H;
//	}
//	
//	return result;
//}


double *WENO::gradV(double F[4][7], double U[4][7], double **L[2], double **R[2], double a)
{
	double G[4][7], V[4][7];// , F[7][4], U[7][4];
	for (size_t i = 0; i < 4; i++)
	{
		for (size_t j = 0; j < 7; j++)
		{
			G[i][j] = F[i][6 - j];
			V[i][j] = U[i][6 - j];
		}
	}

	double *result = new double[4];
	for (size_t c = 0; c < 4; c++)
	{
		double left[4][4][4], right[4][4][4];
		for (size_t s = 0; s < 4; s++)
		{
			double *l = L[0][s], *ll = L[1][s];
			double f[5], ff[5], g[5], gg[5], u[5], uu[5], v[5], vv[5];
			
			for (size_t n = 0; n < 4; n++)
			{
				HELP::dot(ll[n], F[c] + 1, f, 5);
				HELP::dot(l[n], F[c], ff, 5);

				HELP::dot(l[n], G[c] + 1, g, 5);
				HELP::dot(ll[n], G[c], gg, 5);

				HELP::dot(ll[n], U[c] + 1, u, 5);
				HELP::dot(l[n], U[c], uu, 5);

				HELP::dot(l[n], V[c] + 1, v, 5);
				HELP::dot(ll[n], V[c], vv, 5);

				left[0][s][n] = WENO::flux(ff);
				left[1][s][n] = WENO::flux(g);
				left[2][s][n] = WENO::flux(v);
				left[3][s][n] = WENO::flux(uu);

				right[0][s][n] = WENO::flux(f);
				right[1][s][n] = WENO::flux(gg);
				right[2][s][n] = WENO::flux(vv);
				right[3][s][n] = WENO::flux(u);
			}
		}


		double LEFT = 0, RIGHT = 0;
		for (size_t s = 0; s < 4; s++)
		{
			double l[4]{ 0 }, r[4]{ 0 };
			for (size_t i = 0; i < 4; i++)
			{
				l[i] += HELP::dot(left[i][s], R[0][s]);
				r[i] += HELP::dot(right[i][s], R[1][s]);
			}

			LEFT = (l[0] + l[1] - a * (l[2] - l[3])) / 2.0;
			RIGHT = (r[0] + r[1] - a * (r[2] - r[3])) / 2.0;
		}

		result[c] = (RIGHT - LEFT) / H;
	}

	return result;
}

