#pragma once

#include "weno.h"
#include "viscous.h"

#define STENCIL_MUL(f, g) for (size_t n = 0; n < 6; f[n] *= g[n], n++);
#define STENCIL_ADD(f, g) for (size_t n = 0; n < 6; f[n] += g[n], n++);

#define SLICEX(f) { f[i][j + D], f[i + 1][j + D], f[i + 2][j + D], f[i + 3][j + D], f[i + 4][j + D], f[i + 5][j + D] }
#define SLICEY(f) { f[i + D][j], f[i + D][j + 1], f[i + D][j + 2], f[i + D][j + 3], f[i + D][j + 4], f[i + D][j + 5] }

namespace EQN
{
	extern double rho[D + M + E][D + N + E];
	extern double drho[M][N];

	extern double u[D + M + E][D + N + E];
	extern double du[M][N];

	extern double v[D + M + E][D + N + E];
	extern double dv[M][N];

	extern double e[D + M + E][D + N + E];
	extern double de[M][N];
	
	extern double p[D + M + E][D + N + E];

	double fnrho(size_t i, size_t j);
	double fnrhou(size_t i, size_t j);
	double fnrhov(size_t i, size_t j);
	double fnrhoe(size_t i, size_t j);
}

