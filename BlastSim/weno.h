#pragma once

#include "constants.h"
#include "roe.h"
#include "helper.h"

#define POW2(x) ((x) * (x))
#define EPS 1e-60

namespace WENO
{
	double beta0(double *f);
	double beta1(double *f);
	double beta2(double *f);

	double flux(double *f);
	double flux(double *f, double *w);
	double grad(double f[7], double u[7], double a);
	double *gradV(double F[4][7], double U[4][7], double **L[2], double **R[2], double a);
}
