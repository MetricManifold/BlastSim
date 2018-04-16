#pragma once

#include "equations.h"
#include "weno.h"
#include "cases.h"

// compute k[E] for all values
#define COMPUTE(IND) \
	LOOP k[IND][i][j] = { fnrho(i, j), fnrhou(i, j), fnrhov(i, j), fnrhoe(i, j) };

#define UPDATE(IND, dt) \
	LOOP rho[i + D][j + D] = drho[i][j] + dt * k[IND][i][j].rho; \
	LOOP u[i + D][j + D] = (du[i][j] * drho[i][j] + dt * k[IND][i][j].rhou) / rho[i + D][j + D]; \
	LOOP v[i + D][j + D] = (dv[i][j] * drho[i][j] + dt * k[IND][i][j].rhov) / rho[i + D][j + D]; \
	LOOP e[i + D][j + D] = (de[i][j] * drho[i][j] + dt * k[IND][i][j].rhoe) / rho[i + D][j + D];

#define FITS \
	if (type = INVISCID) LOOPALL VISCOUS::mu[i][j] = 0; \
	else LOOPIN VISCOUS::mu[i][j] = FIT::mu(e[i][j], rho[i][j]); \
	LOOPIN p[i][j] = FIT::p(e[i][j], rho[i][j]); \
	//LOOPIN u[i][j] /= FIT::a(e[i][j], rho[i][j]); \
	LOOPIN v[i][j] /= FIT::a(e[i][j], rho[i][j]);

struct Data { double rho, rhou, rhov, rhoe; };

namespace RK
{
	extern Data k[4][M][N];
	void rungekutta4();
}