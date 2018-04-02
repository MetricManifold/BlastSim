#pragma once

#include "equations.h"
#include "weno.h"

// compute k[E] for all values
#define COMPUTE(IND) \
	LOOP k[IND][i][j] = { fnrho(i, j), fnrhou(i, j), fnrhov(i, j), fnrhoe(i, j) };

#define UPDATE(IND, h) \
	LOOP rho[i + D][j + D] = drho[i][j] + h * k[IND][i][j].rho; \
	LOOP u[i + D][j + D] = (du[i][j] * drho[i][j] + h * k[IND][i][j].u) / rho[i + D][j + D]; \
	LOOP v[i + D][j + D] = (dv[i][j] * drho[i][j] + h * k[IND][i][j].v) / rho[i + D][j + D]; \
	LOOP e[i + D][j + D] = (de[i][j] * drho[i][j] + h * k[IND][i][j].e) / rho[i + D][j + D];

struct Data { double rho, u, v, e; };

namespace RK
{
	extern Data k[4][M][N];
	void rungekutta4();
}