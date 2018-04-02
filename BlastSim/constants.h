#pragma once

#include <cmath>

/*
 * dimensions of the system
 */
#define X0 0
#define Y0 0
#define XM 3.52
#define YN 0.02
#define HOB ((XM - X0) / 2)
#define R0 0.3

/* 
 * finite difference method constants
 */
#define C 0.6			// courant number
#define K (H * C)		// time increment
#define H (0.02)	// spatial increment
#define D 3				// left boundary thickness
#define E 2				// right boundary thickness
#define M size_t((XM - X0) / H)
#define N size_t((YN - Y0) / H)

/*
 * gas constants
 */
#define Re 50			// Reynold's number
#define Rgas 287.06		// ideal gas constant

/*
 * initial values
 */
#define MU0 1.7894e-5	// kg m^-1 s^-1
#define RHO0 1.225		// kg m^-3
#define E0 78408.4		// m^2 s^-2
#define P0 101325		// N m^-2
#define T0 288.15		// K
#define K0 1.87915e-2	// J (K M S)^1 (thermal conductivity)
#define RHOWALL 10e5

 // value of y
#define y (j * H + Y0 + H)
#define x (i * H + X0 + H)

/*
 * the type of flow, determined by input
 */
static enum { INVISCID, LAMINAR, TURBULENT } type;


// 1 for 2d axisymmetric flow
// 0 for 2d planar flow
#define A 0

/*
 * loop algorithms
 */
#define LOOP for (size_t i = 0; i < M; i++) for (size_t j = 0; j < N; j++)
#define LOOPIN for (size_t i = D; i < M + D; i++) for (size_t j = D; j < N + D; j++)
#define LOOPALL for (size_t i = 0; i < D + M + E; i++) for (size_t j = 0; j < D + N + E; j++)
#define LOOPLEFT for (size_t i = D - 1; i + 1 > 0; i--) for (size_t j = D; j < D + N; j++)
#define LOOPRIGHT for (size_t i = D + M; i < D + M + E; i++) for (size_t j = D; j < D + N; j++)
#define LOOPBOTTOM for (size_t i = D; i < D + M; i++) for (size_t j = D - 1; j + 1 > 0; j--)
#define LOOPTOP for (size_t i = D; i < D + M; i++) for (size_t j = D + N; j < D + N + E; j++)

