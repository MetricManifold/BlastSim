
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

double maxX(size_t j)
{
	double a = 0;
	for (size_t i = D; i < D + M; i++)
	{
		double b = abs(EQN::u[i][j] + sqrt(ROE::aa(i, j, 0)));
		if (a < b)
		{
			a = b;
		}
	}

	for (size_t i = D; i < D + M; i++)
	{
		double b = abs(EQN::u[i][j] - sqrt(ROE::aa(i, j, 0)));
		if (a < b)
		{
			a = b;
		}
	}

	return a;
}

double maxY(size_t i)
{
	double a = 0;
	for (size_t j = D; j < D + N; j++)
	{
		double b = 2 * abs(EQN::v[i][j]);
		if (a < b)
		{
			a = b;
		}
	}

	return a;
}

double *EQN::system(size_t i, size_t j)
{
	double U[4][7] = { SLICEX(rho), SLICEX(u), SLICEX(v), SLICEX(e) };
	double F[4][7] = { SLICEX(rho), SLICEX(u), SLICEX(v), SLICEX(e) };
	double P[] = SLICEX(p);

	// F_1
	STENCIL_MUL(F[0], U[1]);

	// F_2
	STENCIL_MUL(F[1], U[0]);
	STENCIL_MUL(F[1], U[1]);
	STENCIL_ADD(F[1], P);

	// F_3
	STENCIL_MUL(F[2], U[0]);
	STENCIL_MUL(F[2], U[1]);
	STENCIL_MUL(F[2], U[2]);

	// F_4
	STENCIL_MUL(F[3], U[0]);
	STENCIL_ADD(F[3], P);
	STENCIL_MUL(F[3], U[1]);

	double **Lm = ROE::lF(i + D, j + D, -1);
	double **Lp = ROE::lF(i + D, j + D, 1);
	double **Rm = ROE::rF(i + D, j + D, -1);
	double **Rp = ROE::rF(i + D, j + D, 1);

	double **L[2] = { Lm, Lp };
	double **R[2] = { Rm, Rp };

	return WENO::gradV(F, U, L, R, maxX(j + D));

	for (size_t i = 0; i < 4; i++)
	{
		delete[] Lm[i];
		delete[] Lp[i];
		delete[] Rm[i];
		delete[] Rp[i];
	}

	delete[] Lm;
	delete[] Lp;
	delete[] Rm;
	delete[] Rp;
}

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

	return -((WENO::grad(rhou, rhox, maxX(j + D))));
	// + WENO::grad(rhov, rhoy, maxY(v, i + D))) + A / y * rhov[D]);
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

	return -((WENO::grad(rhoup, rhox, maxX(j + D))));
	// + WENO::grad(rhovu, rhoy, maxY(v, i + D))) +
		//A / y * rhovu[D]) + VISCOUS::fnrhou(i, j) / Re;
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

	return -((WENO::grad(rhouv, rhox, maxX(j + D))));
	// + WENO::grad(rhovp, rhoy, maxY(v, i + D))) +
		//A / y * (rhovp[D] - p[i + D][j + D])) + VISCOUS::fnrhov(i, j) / Re;
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

	return -((WENO::grad(rhopu, rhox, maxX(j + D))));
	// + WENO::grad(rhopv, rhoy, maxY(v, i + D))) +
		//A / y * rhopv[D]) + VISCOUS::fnrhoe(i, j) / Re;
}