
#include "viscous.h"

double VISCOUS::mu[D + M + E][D + N + E];
using namespace EQN;

namespace VISCOUS
{
	double fnrhou(size_t i, size_t j)
	{
		static auto vy = [](size_t i, size_t j) { return mu[i][j] * v[i][j] / y; };

		return (type == INVISCID) ? 0 : (CS::gradx(xxp, i, j) + CS::grady(xy, i, j))
			+ A / y * (xy(i, j) - 2 / 3 * y * CS::gradx(vy, i, j));
	}

	double fnrhov(size_t i, size_t j)
	{
		static auto vy = [](size_t i, size_t j) { return mu[i][j] * v[i][j] / y; };

		return (type == INVISCID) ? 0 : (CS::gradx(xy, i, j) + CS::grady(yyp, i, j))
			+ A / y * (yyp(i, j) - tt(i, j) - 2 / 3 * vy(i, j) - y * 2 / 3 * CS::grady(vy, i, j));
	}

	double fnrhoe(size_t i, size_t j)
	{
		static auto Ev = [](size_t i, size_t j) { return u[i][j] * xxp(i, j) + v[i][j] * xy(i, j) - qx(i, j); };
		static auto Fv = [](size_t i, size_t j) { return u[i][j] * xy(i, j) + v[i][j] * yyp(i, j) - qy(i, j); };
		static auto vy = [](size_t i, size_t j) { return mu[i][j] * v[i][j] / y; };
		static auto vvy = [](size_t i, size_t j) { return mu[i][j] * v[i][j] / y; };
		static auto uvy = [](size_t i, size_t j) { return mu[i][j] * v[i][j] / y; };

		return (type == INVISCID) ? 0 : (CS::gradx(Ev, i, j) + CS::grady(Fv, i, j))
			+ A / y * (mu[i][j] * xy(i, j) + v[i][j] * yyp(i, j) - qy(i, j) -
				2 / 3 * vvy(i, j) - y * 2 / 3 * CS::grady(vvy, i, j) - 
				y * 2 / 3 * CS::gradx(uvy, i , j));
	}

	double xxp(size_t i, size_t j)
	{
		return (type == LAMINAR) ? 0 : mu[i][j] / 3 * (4 * CS::gradx(u, i, j) - 2 * CS::grady(v, i, j));
	}

	double yyp(size_t i, size_t j)
	{
		return (type == LAMINAR) ? 0 : mu[i][j] / 3 * (4 * CS::grady(v, i, j) - 2 * CS::gradx(u, i, j));
	}

	double xy(size_t i, size_t j)
	{
		return (type == LAMINAR) ? 0 : mu[i][j] * (CS::grady(u, i, j) + CS::gradx(v, i, j));
	}

	double tt(size_t i, size_t j)
	{
		return (type == LAMINAR) ? 0 : -2 / 3 * mu[i][j] * (CS::gradx(u, i, j) + 
			CS::grady(v, i, j)) + 4 / 3 * mu[i][j] * v[i][j] / y;
	}

	double qx(size_t i, size_t j)
	{
		static auto T = [](size_t i, size_t j) { return FIT::T(p[i][j], rho[i][j]); };
		return FIT::k(e[i][j], rho[i][j]) * CS::gradx(T, i, j);
	}

	double qy(size_t i, size_t j)
	{
		static auto T = [](size_t i, size_t j) { return FIT::T(p[i][j], rho[i][j]); };
		return FIT::k(e[i][j], rho[i][j]) * CS::grady(T, i, j);
	}
}
