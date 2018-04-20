#include "roe.h"

using namespace EQN;

namespace ROE
{

	// enthalpy
	double h(size_t i, size_t j)
	{
		return rho[i][j] * (FIT::gamma(e[i][j], rho[i][j]) - 1) / (FIT::gamma(e[i][j], rho[i][j]) * p[i][j]);
	}

	double wbarX(size_t i, size_t j, size_t o, size_t r)
	{
		double w[] = {
			(rho[i][j] + rho[i + o][j]) / 2,
			(sqrt(rho[i][j]) * u[i][j] + sqrt(rho[i + o][j]) * u[i + o][j]) / 2,
			(sqrt(rho[i][j]) * v[i][j] + sqrt(rho[i + o][j]) * v[i + o][j]) / 2,
			(sqrt(rho[i][j]) * (e[i][j] + p[i][j] / rho[i][j])
			+ sqrt(rho[i + o][j]) * (e[i + o][j] + p[i + o][j] / rho[i + o][j])) / 2 };

		return w[r];
	}


	// specific heat
	double sh(size_t i, size_t j, size_t o)
	{
		//return (1 / (FIT::gamma(e[i][j], rho[i][j]) * p[i][j] * e[i][j] / rho[i][j] + 1) +
		//	1 / (FIT::gamma(e[i + o][j], rho[i + o][j]) * p[i + o][j] * e[i + o][j] / rho[i + o][j] + 1)) / 2;
		return  (FIT::gamma(e[i][j], rho[i][j]) + FIT::gamma(e[i + o][j], rho[i + o][j])) / 2;
	}

	double qq(size_t i, size_t j, size_t o)
	{
		return Ux * Ux + Vx * Vx;
	}

	double aa(size_t i, size_t j, size_t o)
	{
		return (sh(i, j, o) - 1) * (Hx - 0.5 * qq(i, j, o));
	}

	// each row corresponds to an eigenvector
	// left eigenvalue matrix corresponding to F
	double **lF(size_t i, size_t j, size_t o)
	{
		if (aa(i, j, o) < 0 || isnan(aa(i, j, o)))
		{
			double asdf = Hx;
			double asdfdd = aa(i, j, o);
			double sdf = 0;
		}
		double a = sqrt(aa(i, j, o));
		double d = (1 - sh(i, j, o)) / (-2 * aa(i, j, o));
		double r = Ux * Ux * a - 0.5 * qq(i, j, o) * (Ux + a) + Ux * Hx + Vx * a;
		double s = -Ux * Ux * a - 0.5 * qq(i, j, o) * (Ux - a) + Ux * Hx - Vx * a;
		double t = -0.5 * qq(i, j, o) + Hx;
		double u = 2 * Hx - 2 * Vx - 2 * Ux * Ux;

		double **l = new double*[4]{
			new double[4] { -d * s / a,		d * (t - Ux * a) / a,		-d,				d },
			new double[4] { d * u,			2 * Ux * d,					2 * d,			-2 * d },
			new double[4] { -1,				0,							1 / (EPS + Vx),	0 },
			new double[4] { d * r / a,		-d * (Ux * a + t) / a,		-d,				d } };

		return l;
	}

	// each row corresponds to an eigenvector
	// left eigenvalue matrix corresponding to F
	double **rF(size_t i, size_t j, size_t o)
	{
		double a = sqrt(aa(i, j, o));

		double **r = new double*[4]{
			new double[4] { 1,		Ux - a,		Vx,			Hx - Ux * a },
			new double[4] { 0,		0,			Vx,			Vx * Vx },
			new double[4] { 1,		Ux,			Vx,			0.5 * qq(i, j, o) },
			new double[4] { 1,		Ux + a,		Vx,			Hx + Ux * a } };
		return r;
	}

}