#pragma once

#include "constants.h"
#include "equations.h"

#define RHOW 0

namespace BND
{
	void symmetryY0();
	void noslipX0();
	void noslipXM();
	void noslipYN();
	void tangencyX0();
	void farfieldXM();
	void farfieldYN();

}