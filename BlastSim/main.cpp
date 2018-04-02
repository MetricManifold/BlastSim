
#include <iostream>

#include "rungekutta.h"
#include "boundary.h"
#include "test.h"
#include "initial.h"

using namespace VISCOUS;
using namespace EQN;

int main(int argc, char *argv[])
{
	const size_t end = atoi(argv[1]);

	if (strcmp(argv[2], "INVISCID") == 0) type = INVISCID;
	else if (strcmp(argv[2], "LAMINAR") == 0) type = LAMINAR;
	else if (strcmp(argv[2], "TURBULENT") == 0) type = TURBULENT;
	else exit(1);

	FILE *f;
	f = fopen("output.txt", "w");

	fprintf(f, "# dimensions are %zd,%zd\n", M, N);
	fprintf(f, "# Finite differences; k=%f, h=%f\n", K, H);
	fprintf(f, "# computing %zd time steps\n", end);
	fprintf(f, "# using %s flow\n", argv[2]);

	// initial conditions
	LOOPALL
	{
		rho[i][j] = RHO0;
		u[i][j] = 0;
		v[i][j] = 0;
		e[i][j] = E0;
		p[i][j] = P0;
		mu[i][j] = MU0;
	}

	INIT::tube();

	for (size_t time = 0; time < end; time++)
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

		RK::rungekutta4();

		LOOPIN p[i][j] = FIT::p(e[i][j], rho[i][j]);
		if (type = INVISCID) LOOPALL mu[i][j] = 0;
		else LOOPALL mu[i][j] = FIT::mu(e[i][j], rho[i][j]);

		/*
		if (type == INVISCID)
		{
			BND::tangencyY0(rho);
			BND::tangencyY0(u);
			BND::tangencyY0(v);
		}
		else
		{
			BND::noslipY0(rho);
			BND::noslipY0(u);
			BND::noslipY0(v);

			BND::noslipYN(rho);
			BND::noslipYN(u);
			BND::noslipYN(v);
		}*/


	}
	
	LOOPIN fprintf(f, "%f %f\n", x, rho[i][j] / 1.174);

	return 0;
}