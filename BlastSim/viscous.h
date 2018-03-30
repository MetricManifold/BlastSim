#pragma once

#include "equations.h"
#include "central.h"

namespace TAU
{
	extern double mu[D + M + E][D + N + E];

	double momx(size_t i, size_t j);
	double momy(size_t i, size_t j);
	double ene(size_t i, size_t j);
	
	double xxpdx(size_t i, size_t j);
	double xydx(size_t i, size_t j);
	double yypdy(size_t i, size_t j);
	double xydy(size_t i, size_t j);
	

}