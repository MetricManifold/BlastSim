#include "helper.h"


double HELP::dot(double *a, double *b, size_t n)
{
	double dot = 0;
	for (size_t i = 0; i < n; i++)
	{
		dot += a[i] * b[i];
	}
	return dot;
}

void HELP::dot(double a, double *b, double *out, size_t n)
{
	for (size_t i = 0; i < n; i++)
	{
		out[i] = a * b[i];
	}
}