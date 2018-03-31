#pragma once

#include "constants.h"

#define EPS 1e-6
#define POW2(x) ((x) * (x))

namespace WENO
{
	double beta0(double *f);
	double beta1(double *f);
	double beta2(double *f);

	double flux(double *f, double w[3]);
	double weno(double f[6]);
}
