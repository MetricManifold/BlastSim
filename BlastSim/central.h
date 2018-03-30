#pragma once

#include "constants.h"

namespace CS
{
	double gradx(double f[D + M + E][D + N + E], size_t i, size_t j);
	double grady(double f[D + M + E][D + N + E], size_t i, size_t j);
	double laplx(double f[D + M + E][D + N + E], size_t i, size_t j);
	double laply(double f[D + M + E][D + N + E], size_t i, size_t j);
}