#pragma once

#include "equations.h"
#include "central.h"

namespace TAU
{
	extern double mu[D + M + D][D + N + D];

	void momx(size_t i, size_t j);
	void momy(size_t i, size_t j);
	void ene(size_t i, size_t j);

	void xxpdx(size_t i, size_t j);
	void xydx(size_t i, size_t j);
	void yypdy(size_t i, size_t j);
	void xydy(size_t i, size_t j);
	

}