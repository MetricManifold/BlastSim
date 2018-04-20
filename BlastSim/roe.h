#pragma once

#include "equations.h"


#define Ux (wbarX(i, j, o, 1) / wbarX(i, j, o, 0))
#define Vx (wbarX(i, j, o, 2) / wbarX(i, j, o, 0))
#define Hx (wbarX(i, j, o, 3) / wbarX(i, j, o, 0))

namespace ROE
{

	double wbarX(size_t i, size_t j, size_t o, size_t r);
	double sh(size_t i, size_t j, size_t o);
	double qq(size_t i, size_t j, size_t o);
	double aa(size_t i, size_t j, size_t o);
	double **lF(size_t i, size_t j, size_t o);
	double **rF(size_t i, size_t j, size_t o);

}
