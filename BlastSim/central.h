#pragma once

#include "constants.h"

namespace CS
{
	double gradx(double f[D + M + E][D + N + E], size_t i, size_t j);
	double gradx(double(*f)(size_t, size_t), size_t i, size_t j);
	double grady(double f[D + M + E][D + N + E], size_t i, size_t j);
	double grady(double(*f)(size_t, size_t), size_t i, size_t j);
}