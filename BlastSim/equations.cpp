
#include "equations.h"

double EQN::rho[D + M + E][D + N + E];
double EQN::drho[M][N];

double EQN::u[D + M + E][D + N + E];
double EQN::du[M][N];

double EQN::v[D + M + E][D + N + E];
double EQN::dv[M][N];

double EQN::e[D + M + E][D + N + E];
double EQN::de[M][N];

double EQN::p[D + M + E][D + N + E];

double EQN::fnrho(size_t i, size_t j)
{
	double rhou[] = SLICE_TR(rho);
	double us[] = SLICE_TR(u);
	STENCIL_MUL(rhou, us);

	double rhov[] = SLICE(rho);
	double vs[] = SLICE(u);
	STENCIL_MUL(rhov, vs);

	return WENO::weno(rhou) + WENO::weno(rhov) + A * (rhov[D]) / y;
}

double EQN::fnrhou(size_t i, size_t j)
{

	double rhoup[] = SLICE_TR(rho);
	{
		double ps[] = SLICE_TR(p);
		double usx[] = SLICE_TR(u);
		STENCIL_MUL(usx, usx);
		STENCIL_MUL(rhoup, usx);
		STENCIL_ADD(rhoup, ps);
	}

	double rhovu[] = SLICE(rho);
	{
		double vs[] = SLICE(v);
		double us[] = SLICE(u);
		STENCIL_MUL(rhovu, vs);
		STENCIL_MUL(rhovu, us);
	}

	return WENO::weno(rhoup) + WENO::weno(rhovu) + A * (rhovu[D]) / y;
}

double EQN::fnrhov(size_t i, size_t j)
{
	double rhouv[] = SLICE_TR(rho);
	{
		double us[] = SLICE_TR(u);
		double vs[] = SLICE_TR(v);
		STENCIL_MUL(rhouv, us);
		STENCIL_MUL(rhouv, vs);
	}

	double rhovp[] = SLICE(rho);
	{
		double ps[] = SLICE(p);
		double vs[] = SLICE(v);
		STENCIL_MUL(vs, vs);
		STENCIL_MUL(rhovp, vs);
		STENCIL_ADD(rhovp, ps);
	}

	return WENO::weno(rhovp) + WENO::weno(rhovp) + A * (rhovp[D] - p[i + D][j + D]) / y;
}

double EQN::fnrhoe(size_t i, size_t j)
{
	double rhopu[] = SLICE_TR(rho);
	{
		double es[] = SLICE_TR(e);
		double ps[] = SLICE_TR(p);
		double us[] = SLICE_TR(u);
		STENCIL_MUL(rhopu, es);
		STENCIL_ADD(rhopu, ps);
		STENCIL_MUL(rhopu, us);
	}

	double rhopv[] = SLICE(rho);
	{
		double es[] = SLICE(e);
		double ps[] = SLICE(p);
		double vs[] = SLICE(v);
		STENCIL_MUL(rhopu, es);
		STENCIL_ADD(rhopu, ps);
		STENCIL_MUL(rhopu, vs);
	}

	return WENO::weno(rhopu) + WENO::weno(rhopv) + A * (rhopv[D]) / y;
}