#pragma once

#include "constants.h"
#include <algorithm>

#define POW2(x) ((x) * (x))
#define EPS 1e-40

namespace WENO
{
	double beta0(double *f);
	double beta1(double *f);
	double beta2(double *f);

	double flux(double *f, double d[3]);
	double grad(double f[6], double df[6], double u[6]);
}
