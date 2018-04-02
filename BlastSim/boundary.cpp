
#include "boundary.h"

using namespace EQN;

namespace BND
{

	void symmetryY0(double f[D + M + E][D + N + E])
	{
		LOOPBOTTOM f[i][j] = f[i][D + D - 1 - j];
	}

	void noslipXM(double f[D + M + E][D + N + E])
	{
		if (f == v || f == u)
		{
			LOOPRIGHT f[i][j] = 0;
			LOOPRIGHT f[i][j] = 0;
		}
		else if (f == p)
		{
			LOOPRIGHT f[i][j] = (v[i][j] - 2 * v[i - 1][j] + v[i - 2][j]) / (H * Re) + f[i - 1][j];
		}
		else if (f == rho)
		{
			LOOPRIGHT f[i][j] = f[D + M - 1][j];
		}
	}

	void noslipYN(double f[D + M + E][D + N + E])
	{
		if (f == v || f == u)
		{
			LOOPTOP f[i][j] = 0;
			LOOPTOP f[i][j] = 0;
		}
		else if (f == p)
		{
			LOOPTOP f[i][j] = (v[i][j] - 2 * v[i][j - 1] + v[i][j - 2]) / (H * Re) + f[i][j - 1];
		}
		else if (f == rho)
		{
			LOOPTOP f[i][j] = f[i][D + N - 1];
		}
	}

	void noslipX0(double f[D + M + E][D + N + E])
	{
		if (f == v || f == u)
		{
			LOOPLEFT f[i][j] = 0;
			LOOPLEFT f[i][j] = 0;
		}
		else if (f == p)
		{
			LOOPLEFT f[i][j] = (v[i + 2][j] - 2 * v[i + 1][j] + v[i][j]) / (H * Re) + f[i + 1][j];
		}
		else if (f == rho)
		{
			LOOPLEFT f[i][j] = f[D][j];
		}
	}

	void tangencyX0(double f[D + M + E][D + N + E])
	{
		if (f == v)
		{
			LOOPLEFT f[i][j] = f[D][j];
		}
		else if (f == u)
		{
			LOOPLEFT f[i][j] = 0;
		}
		else if (f == p)
		{
			LOOPLEFT f[i][j] = (v[i + 2][j] - 2 * v[i + 1][j] + v[i][j]) / (H * Re) + f[i + 1][j];
		}
		else if (f == rho)
		{
			LOOPRIGHT f[i][j] = f[D][j];
		}
	}

	void farfieldXM(double f[D + M + E][D + N + E])
	{
		LOOPRIGHT return;
	}

	void farfieldYN(double f[D + M + E][D + N + E])
	{
		LOOPTOP return;
	}

}