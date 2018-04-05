#pragma once

#include "constants.h"

#define POW2(x) ((x) * (x))
#define EPS 1e-6

namespace WENO
{
	double beta0(double *f);
	double beta1(double *f);
	double beta2(double *f);

	double flux(double *f);
	double grad(double f[6]);
}
