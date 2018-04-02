#pragma once

#include "constants.h"
#include "equations.h"


namespace BND
{
	void symmetryY0(double f[D + M + E][D + N + E]);
	void noslipX0(double f[D + M + E][D + N + E]);
	void noslipXM(double f[D + M + E][D + N + E]);
	void noslipYN(double f[D + M + E][D + N + E]);
	void tangencyX0(double f[D + M + E][D + N + E]);
	void farfieldXM(double f[D + M + E][D + N + E]);
	void farfieldYN(double f[D + M + E][D + N + E]);

}