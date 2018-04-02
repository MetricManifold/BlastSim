
#include "rungekutta.h"

using namespace EQN;

namespace RK
{
	struct Data k[4][M][N];

	void rungekutta4()
	{
		// save the current state of the system
		LOOP
		{
			drho[i][j] = rho[i + D][j + D];
			du[i][j] = u[i + D][j + D];
			dv[i][j] = v[i + D][j + D];
			de[i][j] = e[i + D][j + D];
		}

		// perform runge kutta steps
		COMPUTE(0);
		UPDATE(0, H / 2);

		COMPUTE(1);
		UPDATE(1, H / 2);

		COMPUTE(2);
		UPDATE(2, H);

		COMPUTE(3);

		// update the state of the system
		LOOP
		{
			rho[i + D][j + D] = drho[i][j] + H * (k[0][i][j].rho + 2 * k[1][i][j].rho + 2 * k[2][i][j].rho + k[3][i][j].rho) / 6;
			u[i + D][j + D] = (du[i][j] * drho[i][j] + H * (k[0][i][j].u + 2 * k[1][i][j].u + 2 * k[2][i][j].u + k[3][i][j].u) / 6) / rho[i + D][j + D];
			v[i + D][j + D] = (dv[i][j] * drho[i][j] + H * (k[0][i][j].v + 2 * k[1][i][j].v + 2 * k[2][i][j].v + k[3][i][j].v) / 6) / rho[i + D][j + D];
			e[i + D][j + D] = (de[i][j] * drho[i][j] + H * (k[0][i][j].e + 2 * k[1][i][j].e + 2 * k[2][i][j].e + k[3][i][j].e) / 6) / rho[i + D][j + D];
		}
	}
}