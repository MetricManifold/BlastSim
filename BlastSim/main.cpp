
#include <iostream>

#include "weno.h"
#include "equations.h"

#define LOOP for (size_t i = 0; i < M; i++) for (size_t j = 0; j < N; j++)
#define LOOPALL for (size_t i = 0; i < M + D + D; i++) for (size_t j = 0; j < N + D + D; j++)

// compute k[E] for all values
#define COMPUTE(E) \
LOOP \
{ \
double r = fnrho(i, j); \
k[E][i][j] = {r, fnrhou(i, j) / r, fnrhov(i, j) / r, fnrhoe(i, j) / r }; \
}

#define UPDATE(E, h) \
LOOP \
{ \
rho[i + D][j + D] = drho[i][j] + h * k[E][i][j].rho; \
u[i + D][j + D] = du[i][j] + h * k[E][i][j].u; \
v[i + D][j + D] = dv[i][j] + h * k[E][i][j].v; \
e[i + D][j + D] = de[i][j] + h * k[E][i][j].e; \
}

using namespace EQN;

int main(int argc, char *argv[])
{
	const size_t end = atoi(argv[1]);

	struct { double rho, u, v, e; } k[4][M][N];

	LOOPALL
	{
		rho[i][j] = 1.225;
		u[i][j] = 0;
		v[i][j] = 0;
		e[i][j] = 51.3;
		p[i][j] = 1;
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
	}

	return 0;
}