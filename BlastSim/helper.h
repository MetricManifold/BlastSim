#pragma once

#include "roe.h"

namespace HELP
{
	double dot(double *a, double *b, size_t n = 4);
	void dot(double a, double *b, double *out, size_t n = 4);
}