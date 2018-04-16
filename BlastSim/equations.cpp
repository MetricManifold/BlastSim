
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
	double rhox[] = SLICEX(rho);
	double us[] = SLICEX(u);
	STENCIL_MUL(rhou, us);

	double rhov[] = SLICEY(rho);
	double rhoy[] = SLICEY(rho);
	double vs[] = SLICEY(v);
	STENCIL_MUL(rhov, vs);

	return -((WENO::grad(rhou, us, rhox) + 0 * WENO::grad(rhov, vs, rhoy)) + A / y * rhov[D]);
}

double EQN::fnrhou(size_t i, size_t j)
{
	double rhoup[] = SLICEX(rho);
	double rhox[] = SLICEX(rho);
	double us[] = SLICEX(u);
	{
		double ps[] = SLICEX(p);
		STENCIL_MUL(rhoup, us);
		STENCIL_MUL(rhoup, us);
		STENCIL_ADD(rhoup, ps);
		STENCIL_MUL(rhox, us);
	}

	double rhovu[] = SLICEY(rho);
	double rhoy[] = SLICEY(rho);
	double vs[] = SLICEY(v);
	{
		double us[] = SLICEY(u);
		STENCIL_MUL(rhovu, vs);
		STENCIL_MUL(rhovu, us);
		STENCIL_MUL(rhoy, us);
	}

	return -((WENO::grad(rhoup, us, rhox) + 0*WENO::grad(rhovu, vs, rhoy)) + 
		A / y * rhovu[D]) + VISCOUS::fnrhou(i, j) / Re;
}

double EQN::fnrhov(size_t i, size_t j)
{
	double rhouv[] = SLICEX(rho);
	double rhox[] = SLICEX(rho);
	double us[] = SLICEX(u);
	{
		double vs[] = SLICEX(v);
		STENCIL_MUL(rhouv, us);
		STENCIL_MUL(rhouv, vs);
		STENCIL_MUL(rhox, vs);
	}

	double rhovp[] = SLICEY(rho);
	double rhoy[] = SLICEY(rho);
	double vs[] = SLICEY(v);
	{
		double ps[] = SLICEY(p);
		STENCIL_MUL(vs, vs);
		STENCIL_MUL(rhovp, vs);
		STENCIL_ADD(rhovp, ps);
		STENCIL_MUL(rhoy, vs);
	}

	return -((WENO::grad(rhouv, us, rhox) + 0*WENO::grad(rhovp, vs, rhoy)) + 
		A / y * (rhovp[D] - p[i + D][j + D])) + VISCOUS::fnrhov(i, j) / Re;
}

double EQN::fnrhoe(size_t i, size_t j)
{
	double rhopu[] = SLICEX(rho);
	double rhox[] = SLICEX(rho);
	double us[] = SLICEX(u);
	{
		double es[] = SLICEX(e);
		double ps[] = SLICEX(p);
		STENCIL_MUL(rhopu, es);
		STENCIL_ADD(rhopu, ps);
		STENCIL_MUL(rhopu, us);
		STENCIL_MUL(rhox, es);
	}

	double rhopv[] = SLICEY(rho);
	double rhoy[] = SLICEY(rho);
	double vs[] = SLICEY(v);
	{
		double es[] = SLICEY(e);
		double ps[] = SLICEY(p);
		STENCIL_MUL(rhopv, es);
		STENCIL_ADD(rhopv, ps);
		STENCIL_MUL(rhopv, vs);
		STENCIL_MUL(rhoy, es);
	}

	return -((WENO::grad(rhopu, us, rhox) + 0*WENO::grad(rhopv, vs, rhoy)) + 
		A / y * rhopv[D]) + VISCOUS::fnrhoe(i, j) / Re;
}