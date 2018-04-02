
#include "cases.h"

using namespace EQN;

namespace TUBE
{
	void initial()
	{
		double
			P4 = 100 * P0,
			RHO4 = 2.641,
			E4 = 21823744,
			T4 = 9000,		// is this overdetermined??
			P1 = P0,
			RHO1 = 1.174,
			E1 = 214764.72,
			T1 = 300;

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
	}

	void boundaries()
	{
		BND::noslipYN(rho);
		BND::noslipYN(p);
		BND::noslipYN(u);
		BND::noslipYN(v);

		BND::noslipXM(rho);
		BND::noslipXM(p);
		BND::noslipXM(u);
		BND::noslipXM(v);

		BND::noslipX0(rho);
		BND::noslipX0(p);
		BND::noslipX0(u);
		BND::noslipX0(v);

		BND::symmetryY0(rho);
		BND::symmetryY0(p);
		BND::symmetryY0(u);
		BND::symmetryY0(v);
		BND::symmetryY0(e);
	}
}