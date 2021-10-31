/*
	\file:		driftFunctions.cpp
	\brief:		This file defines structures relevant to drift and speed
				information of the quasi-potential problem.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include <cmath>
#include "linalg.h"
#include "driftFunctions.h"

// scaling parameters for rescaling coordinates
// TODO: refactor this
double xfac = 1, yfac = 1;

// master switch of drift functions
PairDoub masterDrift(
	const PairDoub& z, 
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return PairDoub(
			-params->a * z.x - params->b * z.y, 
			params->a * params->b * z.x - z.y
		);
	case 3:
		return PairDoub(
			-params->a * z.x * z.x * z.x - params->b * z.y * z.y * z.y,
			params->a * params->b * z.x * z.x * z.x - z.y * z.y * z.y
		);
	case 5:
		return PairDoub(
			-2 * z.x - 0.75 * params->a * z.x * z.x - params->b * z.y, 
			2 * params->b * z.x + .75 * params->a * params->b * z.x * z.x - z.y
		);
	case 6:
		return PairDoub(
			-z.x * xfac * (params->a * params->a - 1 - params->a * z.x * xfac 
				+ 1.0 / 3.0 * z.x * xfac * z.x * xfac) - z.y * yfac, 
			params->b * z.x * xfac
		);
	case 8:
		return PairDoub(
			-2 * z.x - z.x * z.x * z.x - params->b * z.x * z.y * z.y 
				+ 3 * z.x * z.x + params->b * z.y * z.y, 
			-2 * z.y - z.x * z.x * z.y + 2 * z.x * z.y
		);

	}
}

// master switch of gradient of drift functions
MatrixDoub masterGradDrift(
	const PairDoub& z, 
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return MatrixDoub(
			PairDoub(
				-params->a, 
				-params->b
			), 
			PairDoub(
				params->a * params->b, 
				-1
			)
		);
	case 3:
		return MatrixDoub(
			PairDoub(
				-3 * params->a * z.x * z.x, 
				-3 * params->b * z.y * z.y
			),
			PairDoub(
				3 * params->a * params->b * z.x * z.x, 
				-3 * z.y * z.y
			)
		);
	case 5:
		return MatrixDoub(
			PairDoub(
				-2 - 1.5 * params->a * z.x, 
				-params->b
			), 
			PairDoub(
				2 * params->b + 1.5 * params->a * params->b * z.x, 
				-1
			)
		);
	case 6:
		return MatrixDoub(
			PairDoub(
				xfac * (1 - params->a * params->a + 2 * params->a * xfac * z.x 
					- xfac * z.x * xfac * z.x), 
				-yfac
			), 
			PairDoub(
				xfac * params->b, 
				0
			)
		);
	case 8:
		return MatrixDoub(
			PairDoub(
				-2 - 3 * z.x * z.x - params->b * z.y * z.y + 6 * z.x, 
				-2 * params->b * z.x * z.y + 2 * params->b * z.y
			), 
			PairDoub(
				-2 * z.x * z.y + 2 * z.y, 
				-2 - z.x * z.x + 2 * z.x
			)
		);
	}
}

// master switch of Hessian drift functions
TensorDoub masterHessDrift(
	const PairDoub& z, 
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;
	switch (params->switchKey) {
	case 0: 
		return TensorDoub();
	case 3:
		return TensorDoub(
			MatrixDoub(
				PairDoub(
					-6 * params->a * z.x, 
					0
				), 
				PairDoub(
					0, 
					-6 * params->b * z.y
				)
			),
			MatrixDoub(
				PairDoub(
					6 * params->a * params->b * z.x, 
					0
				), 
				PairDoub(
					0, 
					-6 * z.y
				)
			)
		);
	case 5:
		return TensorDoub(
			MatrixDoub(
				PairDoub(
					-1.5 * params->a, 
					0
				), 
				PairDoub(
					0, 
					0
				)
			), 
			MatrixDoub(
				PairDoub(
					1.5 * params->a * params->b, 
					0
				), 
				PairDoub(
					0, 
					0
				)
			)
		);
	case 6:
		return TensorDoub(
			MatrixDoub(
				PairDoub(
					xfac * xfac * (2 * params->a - 2 * xfac * z.x), 
					0
				), 
				PairDoub(
					0, 
					0
				)
			), 
			MatrixDoub(
				PairDoub(
					0, 
					0
				), 
				PairDoub(
					0, 
					0
				)
			)
		);
	case 8:
		return TensorDoub(
			MatrixDoub(
				PairDoub(
					-6 * z.x + 6, 
					-2 * params->b * z.y
				), 
				PairDoub(
					-2 * params->b * z.y, 
					-2 * params->b * z.x + 2 * params->b
				)
			), 
			MatrixDoub(
				PairDoub(
					-2 * z.y, 
					-2 * z.x + 2
				), 
				PairDoub(
					-2 * z.x + 2, 
					0
				)
			)
		);
	}
}

// master switch of divergence drift functions
double masterDivDrift(
	const PairDoub& z,
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;

	switch (params->switchKey) {
	case 0:
		return params->a * 2 + 2;
	case 5:
		return -3 - 3 * z.x;
	case 6:
		return 0;
	case 8:
		return 2 * (z.x * z.x 
			+ z.y * z.y 
			+ .25 * z.x * z.x * z.x * z.x 
			+ 0.5 * z.x * z.x * z.y * z.y 
			- z.x * z.x * z.x 
			- z.x * z.y * z.y);
	case 7:
		double B11 = xfac * (1.0 - params->a * params->a), 
			B12 = -yfac, 
			B21 = xfac * params->b, 
			B22 = 0.0;
		double aux1 = B21 - B12, 
			aux2 = B11 + B22, 
			aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return Q11 * z.x * z.x + (Q12 + Q21) * z.x * z.y + Q22 * z.y * z.y;
	}
}

// master switch of string functions
std::string masterStr(
	const DriftParams& spdPrms
) {
	return "";
}

// master switch of solution functions
double masterSol(
	const PairDoub& z, 
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;

	switch (params->switchKey) {
	default:
	case 0:
		return params->a * z.x * z.x + z.y * z.y;
	case 1:
		return sqrt(z.x * z.x + z.y * z.y);
	case 2:
		return z.x * z.x / 2 + 1 - cos(z.x + z.y);
	case 3:
		return 0.5 * params->a * z.x * z.x * z.x * z.x 
			+ 0.5 * z.y * z.y * z.y * z.y;
	case 4:
		return z.x * z.x + z.y * z.y + params->b * sin(z.x) * sin(z.x);
	case 5:
		return 2 * z.x * z.x + .5 * params->a * z.x * z.x * z.x + z.y * z.y;
	case 6:
		return 0;
	case 8:
		return 2 * (z.x * z.x + z.y * z.y + .25 * z.x * z.x * z.x * z.x 
			+ 0.5 * z.x * z.x * z.y * z.y - z.x * z.x * z.x - z.x * z.y * z.y);
	case 7:
		double B11 = xfac * (1.0 - params->a * params->a), 
			B12 = -yfac, 
			B21 = xfac * params->b, 
			B22 = 0.0;
		double aux1 = B21 - B12, 
			aux2 = B11 + B22, 
			aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return Q11 * z.x * z.x + (Q12 + Q21) * z.x * z.y + Q22 * z.y * z.y;

	}
}

// master switch of gradient solution functions
PairDoub masterGradSol(
	const PairDoub& z, 
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return PairDoub(2 * params->a * z.x, 2 * z.y);
	case 1:
		return PairDoub(z) / z.norm();
	case 2:
		return PairDoub(z.x + sin(z.x + z.y), sin(z.x + z.y));
	case 3:
		return PairDoub(2 * params->a * z.x * z.x * z.x, 2 * z.y * z.y * z.y);
	case 4:
		return PairDoub(2 * z.x + params->b * sin(2 * z.x), 2 * z.y);
	case 5:
		return PairDoub(4 * z.x + 1.5 * params->a * z.x * z.x, 2 * z.y);
	case 6:
		return PairDoub(0, 0);
	case 8:
		return PairDoub(2 * (2 * z.x + z.x * z.x * z.x + z.x * z.y * z.y 
			- 3 * z.x * z.x - z.y * z.y), 
			2 * (2 * z.y + z.x * z.x * z.y - 2 * z.x * z.y));
	case 7:
		double B11 = xfac * (1.0 - params->a * params->a), 
			B12 = -yfac, 
			B21 = xfac * params->b,
			B22 = 0.0;
		double aux1 = B21 - B12, 
			aux2 = B11 + B22, 
			aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return PairDoub(2 * Q11 * z.x + (Q12 + Q21) * z.y, 
			(Q12 + Q21) * z.y + 2 * Q22 * z.y);

	}
}

// master switch of laplacian solution
double masterLaplacianSol(
	const PairDoub& z, 
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;

	switch (params->switchKey) {
	default:
		return 0;
	case 0:
		return params->a * 2 + 2;
	case 5:
		return 6 + 3 * params->a * z.x;
	case 6:
		return 0;
	case 8:
		return 2 * (z.x * z.x + z.y * z.y + .25 * z.x * z.x * z.x * z.x 
			+ 0.5 * z.x * z.x * z.y * z.y - z.x * z.x * z.x - z.x * z.y * z.y);
	case 7:
		double B11 = xfac * (1.0 - params->a * params->a), 
			B12 = -yfac, 
			B21 = xfac * params->b,
			B22 = 0.0;
		double aux1 = B21 - B12, 
			aux2 = B11 + B22, 
			aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return Q11 * z.x * z.x + (Q12 + Q21) * z.x * z.y + Q22 * z.y * z.y;
	}
}


// master switch of lineariztion of drift functions
MatrixDoub masterDriftLinearization(
	const void* sp
) {
	const DriftParams* params = (DriftParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return MatrixDoub(PairDoub(-params->a, -params->b), 
			PairDoub(params->a * params->b, -1));
	case 1:	case 2:
		return MatrixDoub(PairDoub(INFINITY, INFINITY), 
			PairDoub(INFINITY, INFINITY));
	case 3:
		MatrixDoub(PairDoub(0, 0), PairDoub(0, 0));
	case 4:
		return MatrixDoub(PairDoub(-(1 + params->b), -1), 
			PairDoub(1 + params->b, -1));
	case 5:
		return MatrixDoub(PairDoub(-2, -params->b), 
			PairDoub(2 * params->b, -1));
	case 6:
		return MatrixDoub(PairDoub(xfac * (1 - params->a * params->a), -yfac), 
			PairDoub(xfac * params->b, 0));
	case 7:
		return MatrixDoub(PairDoub(xfac * (1 - params->a * params->a), -yfac), 
			PairDoub(xfac * params->b, 0));
	case 8:
		return MatrixDoub(PairDoub(-2, 0), PairDoub(-2, 0));
	}
}

