
#include "rungekutta.h"

using namespace EQN;
using namespace TUBE;

#include <iostream>
namespace RK
{
	struct Data k[4][M][N];

	void rungekutta4()
	{
		// save the current state of the system
		LOOP drho[i][j] = rho[i + D][j + D];
		LOOP du[i][j] = u[i + D][j + D];
		LOOP dv[i][j] = v[i + D][j + D];
		LOOP de[i][j] = e[i + D][j + D];
		
		// perform runge kutta steps
 		COMPUTE(0);
		UPDATE(0, K / 2);
		FITS;
		boundaries();

		COMPUTE(1);
		UPDATE(1, K / 2);
		FITS;
		boundaries();
		
		COMPUTE(2);
		UPDATE(2, K);
		FITS;
		boundaries();
		
		COMPUTE(3);
		LOOP rho[i + D][j + D] = drho[i][j] + K * (k[0][i][j].rho + 2 * k[1][i][j].rho + 2 * k[2][i][j].rho + k[3][i][j].rho) / 6;
		LOOP u[i + D][j + D] = (du[i][j] + K * (k[0][i][j].u + 2 * k[1][i][j].u + 2 * k[2][i][j].u + k[3][i][j].u) / 6);
		LOOP v[i + D][j + D] = (dv[i][j] + K * (k[0][i][j].v + 2 * k[1][i][j].v + 2 * k[2][i][j].v + k[3][i][j].v) / 6);
		LOOP e[i + D][j + D] = (de[i][j] + K * (k[0][i][j].e + 2 * k[1][i][j].e + 2 * k[2][i][j].e + k[3][i][j].e) / 6);
		FITS;
		boundaries();

	}
}