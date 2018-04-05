
#include "cases.h"

using namespace EQN;

namespace TUBE
{
	double
		P4 = 100 / P0,
		RHO4 = 2.641 / RHO0,
		E4 = 5216 / E0,
		T4 = 9000 / T0,		// is this overdetermined??
		P1 = P0 / P0,
		RHO1 = RHO0 / RHO0,
		E1 = E0 / E0,
		T1 = T0 / T0;

	void initial()
	{
		LOOPIN
		{
			if (i < M / 2)
			{
				p[i][j] = P4;
				rho[i][j] = RHO4;
				e[i][j] = E4;
			}
			else
			{
				p[i][j] = P1;
				rho[i][j] = RHO1;
				e[i][j] = E1;
			}
		}

		boundaries();
	}

	void boundaries()
	{
		BND::noslipYN();
		BND::noslipXM();
		BND::noslipX0();
		BND::symmetryY0();
	}
}
