#pragma once

#include "constants.h"
#include <cmath>

/*
* curve constants
*/
#define E_CURVE 78408.4		// m^2 s^-2
#define P_CURVE 101325		// N m^-2, Pa
#define RHO_CURVE 1.243		// kg m^-3
#define T_CURVE 273.15		// K
#define A_CURVE 331.23		// m s^-1

/*
 * conversions
 */
#define KCAL_TO_M2_PS2 4184
#define ATM_TO_N_PM2 101325

namespace FIT
{
	// viscosity
	double mu(double e, double rho);

	// pressure
	double p(double e, double rho);

	// speed of sound
	double a(double e, double rho);

	// temperature
	double T(double p, double rho);

	// coefficient of thermal conductivity
	double k(double e, double rho);
}
