#pragma once

#include "constants.h"
#include "equations.h"
#include "boundary.h"

namespace TUBE
{
	extern double P4, RHO4, E4, T4, P1, RHO1, E1, T1;

	void initial();
	void boundaries();
}