#pragma once

#include "weno.h"
#include "viscous.h"

#define STENCIL_MUL(f, g) for (size_t n = 0; n < 5; f[n] *= g[n], n++);
#define STENCIL_ADD(f, g) for (size_t n = 0; n < 5; f[n] += g[n], n++);

#define SLICE_TR(f) { f[i][j], f[i + 1][j], f[i + 2][j], f[i + 3][j], f[i + 4][j] }
#define SLICE(f) { f[i][j], f[i][j + 1], f[i][j + 2], f[i][j + 3], f[i][j + 4] }

namespace EQN
{
	extern double rho[D + M + D][D + N + D];
	extern double drho[M][N];

	extern double u[D + M + D][D + N + D];
	extern double du[M][N];

	extern double v[D + M + D][D + N + D];
	extern double dv[M][N];

	extern double e[D + M + D][D + N + D];
	extern double de[M][N];
	
	extern double p[D + M + D][D + N + D];

	double fnrho(size_t i, size_t j);
	double fnrhou(size_t i, size_t j);
	double fnrhov(size_t i, size_t j);
	double fnrhoe(size_t i, size_t j);
}

