
#include "viscous.h"

double TAU::mu[D + M + D][D + N + D];

namespace TAU
{
	void xxpdx(size_t i, size_t j)
	{
		4 / 3 * CS::gradx(mu, i, j);
	}
}
