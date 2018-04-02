#pragma once

#include "constants.h"
#include <cmath>

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
