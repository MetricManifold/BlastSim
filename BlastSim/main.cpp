
#include <iostream>

#include "rungekutta.h"
#include "boundary.h"
#include "test.h"
#include "cases.h"

using namespace VISCOUS;
using namespace EQN;

int main(int argc, char *argv[])
{
	const size_t end = atoi(argv[1]);

	if (strcmp(argv[2], "INVISCID") == 0) type = INVISCID;
	else if (strcmp(argv[2], "LAMINAR") == 0) type = LAMINAR;
	else if (strcmp(argv[2], "TURBULENT") == 0) type = TURBULENT;
	else exit(1);

	//TEST::test();

	FILE *f;
	f = fopen("output.txt", "w");

	fprintf(f, "# dimensions are %zd,%zd\n", M, N);
	fprintf(f, "# Finite differences; k=%f, h=%f\n", K, H);
	fprintf(f, "# computing %zd time steps\n", end);
	fprintf(f, "# using %s flow\n", argv[2]);

	// initial conditions
	LOOPALL
	{
		rho[i][j] = 1;
		u[i][j] = 0;
		v[i][j] = 0;
		e[i][j] = 1;
		p[i][j] = 1;
		mu[i][j] = 1;
	}

	TUBE::initial();

	for (size_t time = 0; time < end; time++)
	{
		RK::rungekutta4();

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
	
	LOOPIN fprintf(f, "%f %f\n", x, rho[i][j]);
	fprintf(f, "\n\n");
	LOOPIN fprintf(f, "%f %f\n", x, p[i][j]);

	return 0;
}