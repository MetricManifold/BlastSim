
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
	double rhou[] = SLICEX(rho);
	double us[] = SLICEX(u);
	STENCIL_MUL(rhou, us);

	double rhov[] = SLICEY(rho);
	double vs[] = SLICEY(u);
	STENCIL_MUL(rhov, vs);

	return -((WENO::weno(rhou) + WENO::weno(rhov)) + A / y * rhov[D]);
}

double EQN::fnrhou(size_t i, size_t j)
{
	double rhoup[] = SLICEX(rho);
	{
		double ps[] = SLICEX(p);
		double us[] = SLICEX(u);
		STENCIL_MUL(us, us);
		STENCIL_MUL(rhoup, us);
		STENCIL_ADD(rhoup, ps);
	}

	double rhovu[] = SLICEY(rho);
	{
		double vs[] = SLICEY(v);
		double us[] = SLICEY(u);
		STENCIL_MUL(rhovu, vs);
		STENCIL_MUL(rhovu, us);
	}

	return -((WENO::weno(rhoup) + WENO::weno(rhovu)) + 
		A / y * rhovu[D]) + VISCOUS::fnrhou(i, j);
}

double EQN::fnrhov(size_t i, size_t j)
{
	double rhouv[] = SLICEX(rho);
	{
		double us[] = SLICEX(u);
		double vs[] = SLICEX(v);
		STENCIL_MUL(rhouv, us);
		STENCIL_MUL(rhouv, vs);
	}

	double rhovp[] = SLICEY(rho);
	{
		double ps[] = SLICEY(p);
		double vs[] = SLICEY(v);
		STENCIL_MUL(vs, vs);
		STENCIL_MUL(rhovp, vs);
		STENCIL_ADD(rhovp, ps);
	}

	return -((WENO::weno(rhovp) + WENO::weno(rhovp)) + 
		A / y * (rhovp[D] - p[i + D][j + D])) + VISCOUS::fnrhov(i, j);
}

double EQN::fnrhoe(size_t i, size_t j)
{
	double rhopu[] = SLICEX(rho);
	{
		double es[] = SLICEX(e);
		double ps[] = SLICEX(p);
		double us[] = SLICEX(u);
		STENCIL_MUL(rhopu, es);
		STENCIL_ADD(rhopu, ps);
		STENCIL_MUL(rhopu, us);
	}

	double rhopv[] = SLICEY(rho);
	{
		double es[] = SLICEY(e);
		double ps[] = SLICEY(p);
		double vs[] = SLICEY(v);
		STENCIL_MUL(rhopu, es);
		STENCIL_ADD(rhopu, ps);
		STENCIL_MUL(rhopu, vs);
	}

	return -((WENO::weno(rhopu) + WENO::weno(rhopv)) + 
		A / y * rhopv[D]) + VISCOUS::fnrhoe(i, j);
}