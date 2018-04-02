#pragma once

#include "equations.h"
#include "central.h"
#include "curve.h"

#define LOOPLOCAL for (size_t I = i - 1; I < i + 2; I++) for (size_t J = j - 1; J < j + 2; J++)


namespace VISCOUS
{
	extern double mu[D + M + E][D + N + E];

	double fnrhou(size_t i, size_t j);
	double fnrhov(size_t i, size_t j);
	double fnrhoe(size_t i, size_t j);
	
	double xxp(size_t i, size_t j);
	double yyp(size_t i, size_t j);
	double xy(size_t i, size_t j);
	double tt(size_t i, size_t j);
	double qx(size_t i, size_t j);
	double qy(size_t i, size_t j);
	

}