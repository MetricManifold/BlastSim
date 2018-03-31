
#include <iostream>

#include "weno.h"
#include "equations.h"
#include "curve.h"
#include "test.h"

#define LOOP for (size_t i = 0; i < M; i++) for (size_t j = 0; j < N; j++)
#define LOOPALL for (size_t i = 0; i < D + M + E; i++) for (size_t j = 0; j < D + N + E; j++)
#define LOOPLEFT for (size_t i = 0; i < D; i++) for (size_t j = 0; j < D + N + E; j++)
#define LOOPRIGHT for (size_t i = D + M; i < D + M + E; i++) for (size_t j = 0; j < D + N + E; j++)
#define LOOPBOTTOM for (size_t i = 0; i < D + M + E; i++) for (size_t j = 0; j < D; j++)
#define LOOPTOP for (size_t i = 0; i < D + M + E; i++) for (size_t j = D + N; j < D + N + E; j++)

// compute k[E] for all values
#define COMPUTE(B) \
LOOP \
{ \
k[B][i][j] = {fnrho(i, j), fnrhou(i, j), fnrhov(i, j), fnrhoe(i, j) }; \
}

#define UPDATE(B, h) \
LOOP \
{ \
rho[i + D][j + D] = drho[i][j] + h * k[B][i][j].rho; \
u[i + D][j + D] = (du[i][j] * drho[i][j] + h * k[B][i][j].u) / rho[i + D][j + D]; \
v[i + D][j + D] = (dv[i][j] * drho[i][j] + h * k[B][i][j].v) / rho[i + D][j + D]; \
e[i + D][j + D] = (de[i][j] * drho[i][j] + h * k[B][i][j].e) / rho[i + D][j + D]; \
}

using namespace EQN;
using namespace TAU;

int main(int argc, char *argv[])
{
	const size_t end = atoi(argv[1]);

	struct { double rho, u, v, e; } k[4][M][N];

	testCurveT();

	// initial conditions
	LOOPALL
	{
		rho[i][j] = RHO0;
		u[i][j] = 0;
		v[i][j] = 0;
		e[i][j] = E0;
		p[i][j] = P0;
		mu[i][j] = MU0;
	}

	for (size_t time = 0; time < end; time++)
	{

		// save the current state of the system
		LOOP
		{
			drho[i][j] = rho[i + D][j + D];
			du[i][j] = u[i + D][j + D];
			dv[i][j] = v[i + D][j + D];
			de[i][j] = e[i + D][j + D];
		}

		COMPUTE(0);
		UPDATE(0, H / 2);

		COMPUTE(1);
		UPDATE(1, H / 2);

		COMPUTE(2);
		UPDATE(2, H);

		COMPUTE(3);
		
		// update the state of the system
		LOOP
		{
			rho[i + D][j + D] = drho[i][j] + H * (k[0][i][j].rho + 2 * k[1][i][j].rho + 2 * k[2][i][j].rho + k[3][i][j].rho) / 6;
			u[i + D][j + D] = du[i][j] + H * (k[0][i][j].u + 2 * k[1][i][j].u + 2 * k[2][i][j].u + k[3][i][j].u) / 6;
			v[i + D][j + D] = dv[i][j] + H * (k[0][i][j].v + 2 * k[1][i][j].v + 2 * k[2][i][j].v + k[3][i][j].v) / 6;
			e[i + D][j + D] = de[i][j] + H * (k[0][i][j].e + 2 * k[1][i][j].e + 2 * k[2][i][j].e + k[3][i][j].e) / 6;
		}

		// symmetry condition on x = 0
		LOOPLEFT
		{
			rho[i][j] = rho[D + D - 1 - i][j];
			u[i][j] = u[D + D - 1 - i][j];
			v[i][j] = v[D + D - 1 - i][j];
			e[i][j] = e[D + D - 1 - i][j];
		}

		// non reflecting
		LOOPRIGHT
		{

		}

		// non reflecting
		LOOPTOP
		{

		}

		// no-slip condition for viscous, or tangency for inviscid
		LOOPBOTTOM
		{

		}

		LOOPALL
		{
			mu[i][j] = FIT::mu(e[i][j], rho[i][j]);
		}
	}

	return 0;
}