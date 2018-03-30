#include "central.h"

namespace CS
{
	double gradx(double f[D + M + E][D + N + E], size_t i, size_t j)
	{
		return (f[i + 1][j] - f[i - 1][j]) / (H + H);
	}

	double grady(double f[D + M + E][D + N + E], size_t i, size_t j)
	{
		return (f[i][j + 1] - f[i][j - 1]) / (H + H);
	}

	double laplx(double f[D + M + E][D + N + E], size_t i, size_t j)
	{
		return (f[i + 1][j] - 2 * f[i][j] + f[i - 1][j]) / (H * H);
	}

	double laply(double f[D + M + E][D + N + E], size_t i, size_t j)
	{
		return (f[i][j + 1] - 2 * f[i][j] + f[i][j - 1]) / (H * H);
	}
}
