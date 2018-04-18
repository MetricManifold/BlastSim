#pragma once

#include <iostream>

#include "constants.h"
#include "equations.h"
#include "curve.h"

namespace TEST
{
	void test();
	void testCurveMu();
	void testCurveP();
	void testCurveT();
	void testCurveK();
	void testCurveA();
	void testWENO(double(*fn)(double), double(*df)(double), bool(*cn)(double));
}
