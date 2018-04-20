
#include "boundary.h"

using namespace EQN;

namespace BND
{

	void symmetryY0()
	{
		LOOPBOTTOM u[i][j] = u[i][D + D - 1 - j];
		LOOPBOTTOM v[i][j] = v[i][D + D - 1 - j];
		LOOPBOTTOM p[i][j] = p[i][D + D - 1 - j];
		LOOPBOTTOM e[i][j] = e[i][D + D - 1 - j];
		LOOPBOTTOM rho[i][j] = rho[i][D + D - 1 - j];
	}

	void noslipXM()
	{
		LOOPRIGHT u[i][j] = 0;
		LOOPRIGHT v[i][j] = 0;
		LOOPRIGHT p[i][j] = -(v[i][j] - 2 * v[i - 1][j] + v[i - 2][j]) / (H * Re) + p[i - 1][j];//p[D + M - 1][j];
		LOOPRIGHT rho[i][j] = rho[D + M - 1][j];
		LOOPRIGHT e[i][j] = e[D + M - 1][j];
	}

	void noslipYN()
	{
		LOOPTOP u[i][j] = 0;
		LOOPTOP v[i][j] = 0;
		LOOPTOP p[i][j] = -(u[i][j] - 2 * u[i][j - 1] + u[i][j - 2]) / (H * Re) + p[i][j - 1];//p[i][D + N - 1];
		LOOPTOP rho[i][j] = rho[i][D + N - 1];
		LOOPTOP e[i][j] = e[i][D + N - 1];
	}

	void noslipX0()
	{
		LOOPLEFT v[i][j] = 0;
		LOOPLEFT u[i][j] = 0;
		LOOPLEFT p[i][j] = -(v[i + 2][j] - 2 * v[i + 1][j] + v[i][j]) / (H * Re) + p[i + 1][j];//p[D][j];
		LOOPLEFT rho[i][j] = rho[D][j];
		LOOPLEFT e[i][j] = e[D][j];
	}

	void tangencyX0()
	{
		LOOPLEFT u[i][j] = u[D][j];
		LOOPLEFT v[i][j] = 0;
		LOOPLEFT p[i][j] = -(v[i + 2][j] - 2 * v[i + 1][j] + v[i][j]) / (H * Re) + p[i + 1][j];
		LOOPLEFT rho[i][j] = rho[D][j];
		LOOPLEFT e[i][j] = e[D][j];
	}

	void farfieldXM()
	{
		LOOPRIGHT return;
	}

	void farfieldYN()
	{
		LOOPTOP return;
	}

}