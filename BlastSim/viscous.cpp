
#include "viscous.h"

double TAU::mu[D + M + E][D + N + E];

namespace TAU
{
	double xxpdx(size_t i, size_t j)
	{
		return 4 / 3 * CS::gradx(mu, i, j);
	}
}
