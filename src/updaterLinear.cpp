/*
	\file:		updaterLinear.cpp
	\brief:		This file defines update formulas for linear MAP
				approximations.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#include <cmath>
#include "linalg.h"
#include "updaterLinear.h"

// compute one-point update value of u and grad u
double UpdaterLinearBase::calcOnePtUpdateValues(
	const SolverParams1Pt& p,
	PairDoub& u1grad, 
	Flags& flags
) {
	// Compute the updated value of grad u
	// Namely, Du = b(z1) (z1-z2)/|z1-z2| - b(z1)

	PairDoub v( p.z1 - p.z2);
	PairDoub drift_z1 = p.si.drift(p.z1, &p.si.sp);
	u1grad = drift_z1.norm() * v / v.norm() - drift_z1;

	// Set flags
	flags.lambda = 0;
	flags.updateType = 1;
	return calcF1(p);
}

// compute triangle update value of u and grad u
// returns infinity if it cannot find an interior minimizer
double UpdaterLinearBase::calcTriangleUpdateValues(
	const SolverParams2Pt& p,
	PairDoub& u1grad,
	bool& interiorMinimizer, 
	Flags& flags
) {
	// Find stationary point of triangle update F2

	bool interiorSolution = false;
	double lambdaStat = calcStationaryPointF2(0, 1, p, interiorSolution);
	
	// Terminate whenever the solver doesn't return critical points in [0,1].
	if (!interiorSolution) {
		u1grad = PairDoub{ INFINITY,INFINITY };
		interiorMinimizer = false;
		return INFINITY;
	}
	double uTriangle = calcF2(lambdaStat, p);

	// Check if uTriangle is less than u values proposed by endpoints

	SolverParams1Pt z2EndpointParams{ p.z1,p.z2,p.u2,p.u2del,p.si },
		z3EndpointParams{ p.z1,p.z3,p.u3,p.u3del,p.si };
	double uEndpoints = std::min(
		calcF1(z2EndpointParams), 
		calcF1(z3EndpointParams)
	);
	if (uTriangle < uEndpoints) {
		interiorMinimizer = true;
		
		// Compute updated value of grad u
		PairDoub zLam{ (1 - lambdaStat) * p.z2 + lambdaStat * p.z3 };
		PairDoub v{ p.z1 - zLam };
		PairDoub b = p.si.drift(p.z1, &p.si.sp);
		u1grad = b.norm() * v / v.norm() - b;

		// Set flags
		flags.lambda = lambdaStat;
		flags.updateType = 2;
		return uTriangle;
	}
	else
		return INFINITY;
}

// Run Hybrid solver to find stationary points of triangle update formula F2
// See "Afternotes on Numerical Analysis" - G.W. Stewart, Lecture 5, p. 37
double UpdaterLinearBase::calcStationaryPointF2(
	double a0, 
	double a1, 
	const SolverParams2Pt& p,
	bool& interiorSolution
) const {
	double tol = 1e-8; // thLameshold for asserting convergence
	double smallVal = 1e-13; // value for determining sign of numbers
	int maxIter = 25;
	double b = a0, c = a1;
	double fb = calcF2FirstDeriv(b, p), fc = calcF2FirstDeriv(c, p);
	interiorSolution = false;
	if (fb * fc > -tol) {
		return INFINITY;
	}
	double a = 0, t = 0, d = 0;
	double fa = 0, fd = 0;
	double dm = 0, df = 0, ds = 0, dd = 0;
	int iter = 0;
	while (iter++ <= maxIter) { // kill3DNewtonSolver after maxIter iteration
		if (abs(fc) < abs(fb)) {
			t = c; c = b; b = t;
			t = fc; fc = fb; fb = t;
			a = c; fa = fc;
		}
		if (abs(b - c) <= tol) {
			interiorSolution = true; 
			return b;
		}
		dm = (c - b) / 2.0;
		df = (fa - fb);
		if (df == 0)
			ds = dm;
		else
			ds = -fb * (a - b) / df;
		if (ds * dm < -smallVal || abs(ds) > abs(dm))
			dd = dm;
		else
			dd = ds;
		if (abs(dd) < tol)
			dd = 0.5 * (static_cast<int64_t>(dm > 0) 
				- static_cast<int64_t>(dm < 0))* tol;
		d = b + dd;
		fd = calcF2FirstDeriv(d, p);
		if (fd == 0) {
			b = c = d; fb = fc = fd;
			break;
		}
		a = b; b = d;
		fa = fb; fb = fd;
		if (fb * fc > smallVal* smallVal) {
			c = a; fc = fa;
		}
	}
	interiorSolution = true;
	return b;
}

///////////////////////////////////////////////////////////////////////////////
// Endpoint quadrature formulas
///////////////////////////////////////////////////////////////////////////////


// compute one-point update F1
double UpdaterLinearEndpoint::calcF1(
	const SolverParams1Pt& p
) const {
	PairDoub v = p.z1 - p.z2;
	return p.si.calcSlowness(p.z1, v, &p.si.sp) * v.norm() + p.u2;
}

// compute two-point update F2 for given lambda value
// i.e. update z1 from (1-lambda)*z2 + lambda*z3
double UpdaterLinearEndpoint::calcF2(
	double lambda, 
	const SolverParams2Pt& p
) const {
	PairDoub zLam = (1 - lambda) * p.z2 + lambda * p.z3;
	PairDoub vLam = p.z1 - zLam;
	double u1 = p.si.calcSlowness(p.z1, vLam, &p.si.sp) * vLam.norm() 
		+ (1 - lambda) * p.u2 + lambda * p.u3;
	return u1;
}

// compute derivative of two-point update calcF2 with respect to lambda
double UpdaterLinearEndpoint::calcF2FirstDeriv(
	double lambda, 
	const SolverParams2Pt& p
) const {
	const void* sp = &p.si.sp;
	PairDoub zLam = (1 - lambda) * p.z2 + lambda * p.z3;
	PairDoub vLam = p.z1 - zLam;
	PairDoub dz = p.z3 - p.z2;
	double hLam = vLam.norm();
	double u1prime = -p.si.calcSlowness(p.z1, vLam, sp) * dot(vLam, dz) / hLam
		- hLam * dot(p.si.slowGradV(p.z1, vLam, sp), dz)
		+ p.u3 - p.u2;
	return u1prime;
}


///////////////////////////////////////////////////////////////////////////////
// Midpoint quadrature formulas
///////////////////////////////////////////////////////////////////////////////

// compute one-point update F1
double UpdaterLinearMidpoint::calcF1(
	const SolverParams1Pt& p
) const {
	PairDoub v = p.z1 - p.z2;
	return p.si.calcSlowness(0.5 * (p.z1 + p.z2), v, &p.si.sp) * v.norm() 
		+ p.u2;
}

// compute two-point update F2 for given lambda value
double UpdaterLinearMidpoint::calcF2(
	double lambda, 
	const SolverParams2Pt& p
) const {
	const void* spp = &p.si.sp;
	PairDoub zLam = (1 - lambda) * p.z2 + lambda * p.z3;
	PairDoub vLam = p.z1 - zLam;
	PairDoub drift_zLam = (1 - lambda) * p.si.drift(0.5 * (p.z1 + p.z2), spp) 
		+ lambda * p.si.drift(0.5 * (p.z1 + p.z3), spp);
	//PairDoub drift_zLam = p.si.drift(zLam, spp);
	double u1 = drift_zLam.norm() * (vLam).norm() - dot(drift_zLam, vLam) 
		+ (1 - lambda) * p.u2 + lambda * p.u3;
	return u1;
}

// compute derivative of two-point update calcF2 with respect to lambda
double UpdaterLinearMidpoint::calcF2FirstDeriv(
	double lambda, 
	const SolverParams2Pt& p
) const {
	const void* spp = &p.si.sp;
	PairDoub(*b)(const PairDoub&, const void*) = p.si.drift;
	PairDoub zLam = (1 - lambda) * p.z2 + lambda * p.z3;
	PairDoub vLam = p.z1 - zLam;
	PairDoub dz = p.z3 - p.z2;
	double hLam = vLam.norm();
	PairDoub zLamMid = (zLam + p.z1) / 2;

	PairDoub
		drift_zLam = (1 - lambda) * b(0.5 * (p.z1 + p.z2), spp)
		+ lambda * b(0.5 * (p.z1 + p.z3), spp),
		drift_z13Mid = b(0.5 * (p.z1 + p.z3), spp),
		drift_z12Mid = b(0.5 * (p.z1 + p.z2), spp);
	double u1_prime = 
		dot(drift_zLam, drift_z13Mid-drift_z12Mid) 
			* (p.z1 - zLam).norm() / drift_zLam.norm()
		+ drift_zLam.norm() * dot(p.z1 - zLam, p.z2 - p.z3) 
			/ (p.z1 - zLam).norm()
		- dot(drift_z13Mid - drift_z12Mid, p.z1 - zLam)
		- dot(drift_zLam, p.z2 - p.z3)
		+ p.u3 - p.u2;
	return u1_prime;
}