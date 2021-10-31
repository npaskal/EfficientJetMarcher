/*
	\file:		updaterHermite.cpp
	\brief:		This file defines update formulas for Hermite cubic MAP
				approximations.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#include <cassert>
#include <iostream>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "linalg.h"
#include "updaterBase.h"
#include "updaterHermite.h"

double gauss()
{
	double eps = 1.e-50, pi2 = 8 * atan(1);
	return sqrt(eps - 2. * log(eps + rand() / (double)RAND_MAX)) 
		* cos(rand() * pi2 / RAND_MAX);
}

// Auxiliary function declarations.
bool kill3DNewtonSolver(double a0, double a1, double r);
PairDoub my_rotate(const PairDoub& vec, double a);
inline double p0(double r) {
	return 1 - 3 * r * r + 2 * r * r * r;
}
inline double p1(double r) {
	return r * (1 - r) * (1 - r);
}
inline double p0p(double r) {
	return 6 * r * (r - 1);
}
inline double p1p(double r) {
	return 1 - 4 * r + 3 * r * r;
}
inline double p0pp(double r) {
	return 12 * r - 6;
}
inline double p1pp(double r) {
	return -4 + 6 * r;
}



/*
*******************************************************************************
*******************************************************************************
							UpdaterHermite Solver -- Quasipotential
*******************************************************************************
*******************************************************************************
For explanation of all the math, see the documentation.
*/

double Q(double a0, double a1, double r) {
	return a0 * r * (1 - r) * (1 - r) - a1 * r * r * (1 - r);
}

double Qp(double a0, double a1, double r) {
	return a0 * (3 * r * r - 4 * r + 1) + a1 * (3 * r * r - 2 * r);
}

// Given input vector v, returns [1,-a; a, 1] * v
PairDoub my_rotate(
	const PairDoub& vec, 
	double a
) {
	return rightMultiply(MatrixDoub(PairDoub(1, -a), PairDoub(a, 1)), vec);
}

// compute one-point update value of u and grad u
double UpdaterHermite::calcOnePtUpdateValues(
	const SolverParams1Pt& p, 
	PairDoub& u1grad,
	Flags& flags
) {
	bool converges = false;
	std::vector<double> outputs = calcStationaryPointF1(0, 0, p, converges);
	double a0, a1;

	// if solver does not converge, use linear MAP
	if (!converges) {
		a0 = 0; 
		a1 = 0;
	}
	else 
	{
		a0 = outputs[0], 
		a1 = outputs[1];
	}

	double updateVal = INFINITY;
	if (m_refinedQuadratureOnePt) 
		updateVal = calcF1Refined(a0, a1, p);
	else 
		updateVal = calcF1(a0, a1, p);

	// Set grad U direction
	PairDoub dir(p.z1 - p.z2);
	dir = my_rotate(dir, a1);
	dir = dir / dir.norm();
	u1grad = p.si.drift(p.z1, &p.si.sp).norm() * dir 
		- p.si.drift(p.z1, &p.si.sp);

	// Update flags
	flags.updateType = 1;
	flags.a0 = a0;
	flags.a1 = a1;
	flags.lambda = 0;

	return updateVal;
}

// compute triangle update value of u and grad u
double UpdaterHermite::calcTriangleUpdateValues(
	const SolverParams2Pt& p, 
	PairDoub& u1grad,
	bool& interiorMinimizer, 
	Flags& flags
) {
	double lambdaStart = 0.5; // TODO - implement guess.
	bool converges = false;
	std::vector<double> stationaryParams = 
		calcStationaryPointF2(0, 0, lambdaStart, p, converges);

	// if no convergence of Newton solver, propose no update value
	if (!converges) {
		interiorMinimizer = false;
		return INFINITY;
	}
	double lambda = stationaryParams[2], 
		a1 = stationaryParams[1], 
		a0 = stationaryParams[0];

	// kill update if the stationary point falls out of the relevant region.
	double angThreshold = 3; // Upper bound on acceptable |tan(a)| values.
	if (
		lambda > 1 
		|| lambda < 0 
		|| std::max(abs(a0), abs(a1)) > angThreshold
	) {
		interiorMinimizer = false;
		return INFINITY;
	}

	interiorMinimizer = true;
	// this point is only a stationary point, we do not check if its a min

	// compute grad u value
	PairDoub dir = (p.z1 - (1 - lambda) * p.z2 - lambda * p.z3);
	dir = my_rotate(dir, a1);
	dir = dir / dir.norm();
	u1grad = p.si.drift(p.z1, &p.si.sp).norm() * dir 
		- p.si.drift(p.z1, &p.si.sp);

	// update flags
	flags.updateType = 2;
	flags.a0 = a0;
	flags.a1 = a1;
	flags.lambda = lambda;

	double updateVal = INFINITY;
	if (m_refinedQuadratureTriangle) 
		updateVal = calcF2Refined(a0, a1, lambda, p);
	else 
		updateVal = calcF2(a0, a1, lambda, p);
	return updateVal;
}

// use Newton's method to find a stationary point of F2 
std::vector<double> UpdaterHermite::calcStationaryPointF1(
	double a0Start, // initial guess for a0
	double a1Start,
	const SolverParams1Pt& p, 
	bool& converges
) {
	m_numOnePtCalls++;

	double tolConv = 1e-4; // tolerance for signifying convergence
	double tolSingularity = 1e-5; // tolConv for evaluating hessian singularity
	int maxIter = 20; // max number of iterations

	std::vector<double> gradVec = calcF1FirstDeriv(a0Start, a1Start, p);
	std::vector<double> hessVec = calcF1SecondDeriv(a0Start, a1Start, p);

	Eigen::Vector2d b(-gradVec[0], -gradVec[1]), // rhs
		s(a0Start, a1Start), // current iter (a0,a1) vector
		sprev(s), // previous iter (a0,a1) vector
		sdiff(0, 0); // difference (a0,a1) vector
	Eigen::Matrix2d A; // Hessian matrix
	A << hessVec[0], hessVec[1], hessVec[2], hessVec[3];
	A.trace();
	int iter = 0;

	std::vector<double> sol; 
	while (++iter <= maxIter) {
		// TODO: come up with a better condition for checking singularity of A
		if (abs(A.determinant()) / pow(abs(A(0)) + abs(A(3)), 2) 
			< tolSingularity
		) {
			converges = false;
			m_numOnePtFailures++;
			return sol;
		}
		Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
		s += x;
		gradVec = calcF1FirstDeriv(s(0), s(1), p);
		hessVec = calcF1SecondDeriv(s(0), s(1), p);
		b(0) = -gradVec[0]; b(1) = -gradVec[1];
		for (int i = 0; i < 4; i++) 
			A(i) = hessVec[i];
		sdiff = s - sprev;
		sprev = s;
		if (sdiff.norm() < tolConv) {
			sol = { s(0),s(1) };
			converges = true;
			m_numOnePtSuccesses++;
			return sol;
		}
	}
	converges = false;
	m_numOnePtFailures++;
	return sol;
}

std::vector<double> UpdaterHermite::calcStationaryPointF2(
	double a0Start, 
	double a1Start, 
	double lambdaStart,
	const SolverParams2Pt& p, 
	bool& converges
) {
	m_numTriangleCalls++;
	double tolConv = 1e-5, // tolerance for signifying convergence
		tolSingularity = 1e-5; 	// tolConv for evaluating hessian singularity
	int numIter = 20;	// max number iterations

	std::vector<double> gradVec 
		= calcF2FirstDeriv(a0Start, a1Start, lambdaStart, p);
	std::vector<double> hessVec 
		= calcF2SecondDeriv(a0Start, a1Start, lambdaStart, p);
	assert(gradVec.size() == 3 && hessVec.size() == 9);
	Eigen::Vector3d b(-gradVec[0], -gradVec[1], -gradVec[2]),
		s(a0Start, a1Start, lambdaStart),
		sprev(s),
		sdiff(0, 0, 0);
	Eigen::Matrix3d A;
	for (int k = 0; k < 9; k++)
		A(k) = hessVec[k];
	int iter = 0;


	std::vector<double> sol;

	while (++iter <= numIter) {
		double a = A.determinant();
		if (abs(a) / pow(abs(A(0)) + abs(A(4)) + abs(A(8)), 3) 
				< tolSingularity
			|| kill3DNewtonSolver(s(0), s(1), s(2))
		) {
			converges = false;
			return sol;
		}
		Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
		s += x;
		gradVec = calcF2FirstDeriv(s(0), s(1), s(2), p);
		hessVec = calcF2SecondDeriv(s(0), s(1), s(2), p);
		b(0) = -gradVec[0]; b(1) = -gradVec[1]; b(2) = -gradVec[2];
		for (int k = 0; k < 9; k++)
			A(k) = hessVec[k];
		sdiff = s - sprev;
		sprev = s;

		if (sdiff.norm() < tolConv) {
			sol = { s(0),s(1), s(2) };
			converges = true;
			return sol;

		}
	}
	converges = false;

	return sol;
}

// condition at which to kill3DNewtonSolver 3D Newton solver
bool kill3DNewtonSolver(double a0, 
	double a1, 
	double lambda) {
	if (lambda > 1.01 || lambda < -.01) 
		return true;
	else 
		return false;
}

/*
*******************************************************************************
Beginning of Formulas
*******************************************************************************
*/


// compute one-point update F1, given entry and exit angles a0, a1
double UpdaterHermite::calcF1(
	double a0, 
	double a1, 
	const SolverParams1Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;

	PairDoub dir = z1 - z2;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + z2) + 0.125 * h * (a0 - a1) * w;
	return u2 + h / 6 * (
		p.si.calcSlowness(z2, my_rotate(dir, a0), &sp) * sqrt(1 + a0 * a0)
		+ 4 * p.si.calcSlowness(y, my_rotate(dir, -(a0 + a1) / 4), &sp)
		* sqrt(1 + (a0 + a1) * (a0 + a1) / 16)
		+ p.si.calcSlowness(z1, my_rotate(dir, a1), &sp) * sqrt(1 + a1 * a1)
		);
}

// compute one-point update using higher-order Simpsons quadrature
double UpdaterHermite::calcF1Refined(
	double a0, 
	double a1,
	const SolverParams1Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;

	PairDoub dir = z1 - z2;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + z2) + 0.125 * h * (a0 - a1) * w;

	double u_update = 0;
	int sumWeights = 0;
	for (unsigned k = 0; k < m_refinedNumOnePtNodes; k++) {
		double g = (double)k / (m_refinedNumOnePtNodes - (int)1);
		int weight = 0;
		if (k == 0 || k == m_refinedNumOnePtNodes - 1) weight = 1;
		else if (k % 2 == 0) weight = 2;
		else weight = 4;
		PairDoub z_curr = (1 - g) * z2 + g * z1 + h * Q(a0, a1, g) * w;
		double qp = Qp(a0, a1, g);
		u_update += weight * p.si.calcSlowness(z_curr, my_rotate(dir, qp), &sp) 
			* sqrt(1 + qp * qp);
		sumWeights += weight;
	}
	return u2 + h * u_update / sumWeights;
}

// compute first-order derivatives of F1 (2D vector)
std::vector<double> UpdaterHermite::calcF1FirstDeriv(
	double a0, 
	double a1,
	const SolverParams1Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) 
		= p.si.gradDrift;

	PairDoub dir = z1 - z2;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + z2) + 0.125 * h * (a0 - a1) * w;
	PairDoub by = drift(y, &sp);
	double bymag = by.norm();
	MatrixDoub dby = gradDrift(y, &sp);
	double sqrtTerm = sqrt(1 + (a0 + a1) * (a0 + a1) / 16);
	double dFda0 = h / 6 * (
		drift(z2, &sp).norm() * a0 / sqrt(1 + a0 * a0)
		- dot(drift(z2, &sp), w)
		+ h / (2 * bymag) * dot(by, rightMultiply(dby, w)) * sqrtTerm
		+ 0.25 * bymag * (a0 + a1) / sqrtTerm
		- 0.5 * dot(rightMultiply(dby, w), my_rotate(dir, -(a0 + a1) / 4))
		+ dot(by, w)
		);
	double dFda1 = h / 6 * (
		drift(z1, &sp).norm() * a1 / sqrt(1 + a1 * a1) - dot(drift(z1, &sp), w)
		- h / (2 * bymag) * dot(by, rightMultiply(dby, w)) * sqrtTerm
		+ 0.25 * bymag * (a0 + a1) / sqrtTerm
		+ 0.5 * dot(rightMultiply(dby, w), my_rotate(dir, -(a0 + a1) / 4))
		+ dot(by, w)
		);
	return { dFda0,dFda1 };

}

// compute second-order derivatives of F1 (as 4D vector)
std::vector<double> UpdaterHermite::calcF1SecondDeriv(
	double a0, 
	double a1,
	const SolverParams1Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) 
		= p.si.gradDrift;
	TensorDoub(*hessDrift) (const PairDoub & z, const void* sp) 
		= p.si.hessDrift;

	PairDoub dir = z1 - z2;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + z2) + 0.125 * h * (a0 - a1) * w;
	PairDoub by = drift(y, &sp);
	double bymag = by.norm();
	MatrixDoub dby = gradDrift(y, &sp);
	TensorDoub ddby = hessDrift(y, &sp);
	PairDoub dbyw = rightMultiply(dby, w);
	double sqrtTerm = sqrt(1 + (a0 + a1) * (a0 + a1) / 16);
	double dFda0a0 = h / 6 * (
		drift(z2, &sp).norm() / pow(1 + a0 * a0, 1.5)
		+ h * (a0 + a1) / (16 * sqrtTerm) 
			* dot(by, rightMultiply(dby, w)) / bymag
		+ h * h / 16 * sqrtTerm / bymag * (
			pow(dbyw.norm(), 2) + dot(by, ddby.sandwich(w, w)) 
			- pow(dot(by / bymag, dbyw), 2)
			)
		+ 1.0 / 4 * bymag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		+ h / 4 * dot(dbyw, w)
		);
	double dFda0a1 = h / 6 * (
		-h * h / 16 * sqrtTerm / bymag * (
			pow(dbyw.norm(), 2) + dot(by, ddby.sandwich(w, w)) 
			- pow(dot(by / bymag, dbyw), 2)
			)
		+ 1.0 / 4 * bymag / pow(sqrtTerm, 3)
		+ h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		);
	double dFda1a1 = h / 6 * (
		drift(z1, &sp).norm() / pow(1 + a1 * a1, 1.5)
		- h * (a0 + a1) / (16 * sqrtTerm) * dot(by, rightMultiply(dby, w)) 
			/ bymag
		+ h * h / 16 * sqrtTerm / bymag * (
			pow(dbyw.norm(), 2) + dot(by, ddby.sandwich(w, w)) 
			- pow(dot(by / bymag, dbyw), 2)
			)
		+ 1.0 / 4 * bymag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		- h / 4 * dot(dbyw, w)
		);
	return { dFda0a0, dFda0a1, dFda0a1, dFda1a1 };
}

// compute two-point update F2, given lambda, a0, a1
double UpdaterHermite::calcF2(
	double a0, 
	double a1, 
	double r,
	const SolverParams2Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2; double u3 = p.u3;
	PairDoub gradU2 = p.u2del, gradU3 = p.u3del;

	PairDoub zlam = (1 - r) * z2 + r * z3;
	PairDoub dir = z1 - zlam;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + zlam) + 0.125 * h * (a0 - a1) * w;

	double u0prime = dot(z3 - z2, gradU2), u1prime = dot(z3 - z2, gradU3);
	double H = u2 * p0(r) + u3 * p0(1 - r) + u0prime * p1(r)
		- u1prime * p1(1 - r);
	return H + h / 6 * (
		p.si.calcSlowness(zlam, my_rotate(dir, a0), &sp) * sqrt(1 + a0 * a0)
		+ 4 * p.si.calcSlowness(y, my_rotate(dir, -(a0 + a1) / 4), &sp)
		* sqrt(1 + (a0 + a1) * (a0 + a1) / 16)
		+ p.si.calcSlowness(z1, my_rotate(dir, a1), &sp) * sqrt(1 + a1 * a1)
		);
}

// compute triangle update using higher-order Simpsons quadrature
double UpdaterHermite::calcF2Refined(double a0, double a1, double  r,
	const SolverParams2Pt& p) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2; double u3 = p.u3;
	PairDoub gradU2 = p.u2del, gradU3 = p.u3del;

	PairDoub zlam = (1 - r) * z2 + r * z3;
	PairDoub dir = z1 - zlam;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + zlam) + 0.125 * h * (a0 - a1) * w;

	double u0prime = dot(z3 - z2, gradU2), u1prime = dot(z3 - z2, gradU3);
	double H =
		u2 * p0(r) + u3 * p0(1 - r) + u0prime * p1(r) - u1prime * p1(1 - r);

	double u_update = 0;
	int sumWeights = 0;
	for (unsigned k = 0; k < m_refinedNumTriangleNodes; k++) {
		double g = (double)k / (m_refinedNumTriangleNodes - (int)1);
		int weight = 0;
		if (k == 0 || k == m_refinedNumTriangleNodes - 1) weight = 1;
		else if (k % 2 == 0) weight = 2;
		else weight = 4;
		PairDoub z_curr = (1 - g) * zlam + g * z1 + h * Q(a0, a1, g) * w;
		double qp = Qp(a0, a1, g);
		double utemp = weight * 
			p.si.calcSlowness(
				z_curr, my_rotate(dir, qp), 
				&sp
			) 
			* sqrt(1 + qp * qp);
		u_update += weight * p.si.calcSlowness(z_curr, my_rotate(dir, qp), &sp) 
			* sqrt(1 + qp * qp);
		sumWeights += weight;
	}
	return H + h * u_update / sumWeights;

}

// compute first-order derivatives of F2 (3D vector)
std::vector<double> UpdaterHermite::calcF2FirstDeriv(
	double a0, 
	double a1, 
	double r,
	const SolverParams2Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2, u3 = p.u3;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) 
		= p.si.gradDrift;
	PairDoub gradU2 = p.u2del, gradU3 = p.u3del;

	PairDoub zlam = (1 - r) * z2 + r * z3;
	PairDoub dir = z1 - zlam;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + zlam) + 0.125 * h * (a0 - a1) * w;
	PairDoub by = drift(y, &sp);
	double bymag = by.norm();
	MatrixDoub dby = gradDrift(y, &sp);
	double sqrtTerm = sqrt(1 + (a0 + a1) * (a0 + a1) / 16);
	double u0prime = dot(z3 - z2, gradU2), u1prime = dot(z3 - z2, gradU3);

	double dFda0 = h / 6 * (
		drift(zlam, &sp).norm() * a0 / sqrt(1 + a0 * a0)
		- dot(drift(zlam, &sp), w)
		+ h / (2 * bymag) * dot(by, rightMultiply(dby, w)) * sqrtTerm
		+ 0.25 * bymag * (a0 + a1) / sqrtTerm
		- 0.5 * dot(rightMultiply(dby, w), my_rotate(dir, -(a0 + a1) / 4))
		+ dot(by, w)
		);
	double dFda1 = h / 6 * (
		drift(z1, &sp).norm() * a1 / sqrt(1 + a1 * a1) - dot(drift(z1, &sp), w)
		- h / (2 * bymag) * dot(by, rightMultiply(dby, w)) * sqrtTerm
		+ 0.25 * bymag * (a0 + a1) / sqrtTerm
		+ 0.5 * dot(rightMultiply(dby, w), my_rotate(dir, -(a0 + a1) / 4))
		+ dot(by, w)
		);
	PairDoub dx = z3 - z2, bzlam = drift(zlam, &sp);
	MatrixDoub dbzlam = gradDrift(zlam, &sp);
	double dFdr = u2 * p0p(r) - u3 * p0p(1 - r) 
		+ u0prime * p1p(r) + u1prime * p1p(1 - r)
		- dot(z1 - zlam, dx) / (6 * h) * sqrt(1 + a0 * a0) * bzlam.norm()
		+ h / 6 * sqrt(1 + a0 * a0) * dot(bzlam, rightMultiply(dbzlam, dx)) 
			/ bzlam.norm()
		- 1.0 / 6 * dot(rightMultiply(dbzlam, dx), my_rotate(dir, a0))
		+ 1.0 / 6 * dot(bzlam, my_rotate(dx, a0))
		- dot(z1 - zlam, dx) / (6 * h) * 4 * sqrtTerm * bymag
		+ h / 3 * sqrtTerm * 
			dot(
				by, 
				rightMultiply(dby, dx - (a0 - a1) / 4 * dx.perp())
			) / bymag
		- 1.0 / 3 * 
			dot(
				rightMultiply(
					dby, 
					dx - (a0 - a1) / 4 * dx.perp()
				), 
				my_rotate(dir, -(a0 + a1) / 4)
			)
		+ 2.0 / 3 * dot(by, my_rotate(dx, -(a0 + a1) / 4))
		- dot(z1 - zlam, dx) / (6 * h) * sqrt(1 + a1 * a1) 
		* drift(z1, &sp).norm()
		+ 1.0 / 6 * dot(drift(z1, &sp), my_rotate(dx, a1))
		;
	return { dFda0,dFda1,dFdr };
}

// compute second-order derivatives of F2 (as 9D vector)
std::vector<double> UpdaterHermite::calcF2SecondDeriv(
	double a0, 
	double a1, 
	double r,
	const SolverParams2Pt& p
) const {
	DriftParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2, u3 = p.u3;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) 
		= p.si.gradDrift;
	TensorDoub(*hessDrift) (const PairDoub & z, const void* sp) 
		= p.si.hessDrift;
	PairDoub gradU2 = p.u2del, gradU3 = p.u3del;

	PairDoub zr = (1 - r) * z2 + r * z3;
	PairDoub dir = z1 - zr;
	double h = dir.norm();

	PairDoub w = dir.perpUnit();
	PairDoub yr = 0.5 * (z1 + zr) + 0.125 * h * (a0 - a1) * w;
	PairDoub delx = z3 - z2,
		q = delx - (a0 - a1) / 4 * delx.perp();
	PairDoub byr = drift(yr, &sp),
		bzr = drift(zr, &sp);
	double byrmag = byr.norm(), bzrmag = bzr.norm();
	PairDoub dbyw = rightMultiply(gradDrift(yr, &sp), w),
		dbyq = rightMultiply(gradDrift(yr, &sp), q),
		dbzrdx = rightMultiply(gradDrift(zr, &sp), delx)
		;
	double hp = -dot(dir, delx) / h;
	double hpp = delx.norm() * delx.norm() / h - dot(dir, delx) 
		* dot(dir, delx) / (h * h * h);
	MatrixDoub dby = gradDrift(yr, &sp);
	TensorDoub ddby = hessDrift(yr, &sp),
		ddbzr = hessDrift(zr, &sp);
	double sqrtTerm = sqrt(1 + (a0 + a1) * (a0 + a1) / 16);
	double u0prime = dot(z3 - z2, gradU2), u1prime = dot(z3 - z2, gradU3);
	double dFda0a0 = h / 6 * (
		bzr.norm() / pow(1 + a0 * a0, 1.5)
		+ h * (a0 + a1) / (16 * sqrtTerm) * dot(byr, dbyw) / byrmag
		+ h * h / 16 * sqrtTerm / byrmag * (
			pow(dbyw.norm(), 2) + dot(byr, ddby.sandwich(w, w)) 
			- pow(dot(byr / byrmag, dbyw), 2)
			)
		+ 1.0 / 4 * byrmag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		+ h / 4 * dot(dbyw, w)
		);
	double dFda0a1 = h / 6 * (
		-h * h / 16 * sqrtTerm / byrmag * (
			pow(dbyw.norm(), 2) + dot(byr, ddby.sandwich(w, w)) 
			- pow(dot(byr / byrmag, dbyw), 2)
			)
		+ 1.0 / 4 * byrmag / pow(sqrtTerm, 3)
		+ h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		);
	double dFda0r =
		h / 6 * a0 / sqrt(1 + a0 * a0) * dot(bzr, dbzrdx) / bzrmag
		- 1.0 / 6 * dot(dbzrdx, dir.perp())
		+ 1.0 / 6 * (dot(bzr - byr, delx.perp()) 
			+ h / 2 * dot(dbyw, my_rotate(delx, -(a0 + a1) / 4)))
		+ hp / 6 * (a0 / sqrt(1 + a0 * a0) * bzrmag 
			+ 0.25 * (a0 + a1) / sqrtTerm * byrmag 
			+ h / 2 * sqrtTerm * dot(byr, dbyw) / byrmag)
		+ h / 48 * (a0 + a1) / sqrtTerm * dot(byr, dbyq) / byrmag
		+ h * h / 24 * sqrtTerm / byrmag * (dot(dbyw, dbyq) 
			+ dot(byr, ddby.sandwich(w, q)))
		- h / 12 * sqrtTerm * dot(byr, rightMultiply(dby, delx.perp()))/ byrmag
		- h * h / 24 * sqrtTerm * dot(byr, dbyq) * dot(byr, dbyw) 
		/ (byrmag * byrmag * byrmag)
		- h / 24 * dot(ddby.sandwich(w, q), my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 12 * dot(rightMultiply(dby, delx.perp()), 
			my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 12 * dot(dbyq, dir.perp())
		;
	double dFda1a1 = h / 6 * (
		drift(z1, &sp).norm() / pow(1 + a1 * a1, 1.5)
		- h * (a0 + a1) / (16 * sqrtTerm) 
		* dot(byr, rightMultiply(dby, w)) / byrmag
		+ h * h / 16 * sqrtTerm / byrmag * (
			pow(dbyw.norm(), 2) 
			+ dot(byr, ddby.sandwich(w, w)) - pow(dot(byr / byrmag, dbyw), 2)
			)
		+ 1.0 / 4 * byrmag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		- h / 4 * dot(dbyw, w)
		);
	double dFda1r =
		+1.0 / 6.0 * (-dot(byr, delx.perp()) 
			- h / 2 * dot(dbyw, my_rotate(delx, -(a0 + a1) / 4)))

		+ 1.0 / 6 * dot(drift(z1, &sp), delx.perp())
		+ hp / 6 * (a1 / sqrt(1 + a1 * a1) * drift(z1, &sp).norm() 
			+ 0.25 * (a0 + a1) / sqrtTerm * byrmag 
			- h / 2 * sqrtTerm * dot(byr, dbyw) / byrmag)
		+ h / 48 * (a0 + a1) / sqrtTerm * dot(byr, dbyq) / byrmag
		- h * h / 24 * sqrtTerm / byrmag * (dot(dbyw, dbyq) 
			+ dot(byr, ddby.sandwich(w, q)))
		+ h / 12 * sqrtTerm* dot(byr, rightMultiply(dby, delx.perp())) / byrmag
		+ h * h / 24 * sqrtTerm * dot(byr, dbyq) * dot(byr, dbyw) 
		/ (byrmag * byrmag * byrmag)
		+ h / 24 * dot(ddby.sandwich(w, q), my_rotate(dir, -(a0 + a1) / 4))
		- 1.0 / 12 * dot(rightMultiply(dby, delx.perp()), 
			my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 12 * dot(dbyq, dir.perp())
		;
	double dFdrr =
		u2 * p0pp(r) + u3 * p0pp(1 - r) 
		+ u0prime * p1pp(r) - u1prime * p1pp(1 - r)
		+ hp / 6 * sqrt(1 + a0 * a0) * dot(bzr, dbzrdx) / bzrmag
		+ h / 6 * sqrt(1 + a0 * a0) / bzrmag * (dbzrdx.norm() * dbzrdx.norm() 
			+ dot(bzr, ddbzr.sandwich(delx, delx)) -
			pow(dot(bzr, dbzrdx) / bzrmag, 2))
		- 1.0 / 6 * dot(ddbzr.sandwich(delx, delx), my_rotate(dir, a0))
		+ 1.0 / 6 * dot(dbzrdx, my_rotate(delx, a0))
		+ 1.0 / 6 * dot(dbzrdx, my_rotate(delx, a0))
		+ 1.0 / 3 * dot(dbyq, my_rotate(delx, -(a0 + a1) / 4))
		+ 1.0 / 6 * hpp * (sqrt(1 + a0 * a0) * bzrmag + 4 * sqrtTerm * byrmag 
			+ sqrt(1 + a1 * a1) * drift(z1, &sp).norm())
		+ 1.0 / 6 * hp * (sqrt(1 + a0 * a0) * dot(bzr, dbzrdx) / bzrmag 
			+ 2 * sqrtTerm * dot(byr, dbyq) / byrmag)
		+ hp / 3 * sqrtTerm * dot(byr, dbyq) / byrmag
		+ h / 6 * sqrtTerm / byrmag * (dbyq.norm() * dbyq.norm() 
			+ dot(byr, ddby.sandwich(q, q)) 
			- dot(byr, dbyq) * dot(byr, dbyq) / (byrmag * byrmag))
		- 1.0 / 6 * dot(ddby.sandwich(q, q), my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 3 * dot(dbyq, my_rotate(delx, -(a0 + a1) / 4))
		;
	return { dFda0a0, dFda0a1,dFda0r,dFda0a1,
		dFda1a1,dFda1r,dFda0r,dFda1r,dFdrr };
}

// double check math is correct. Check derivatives agree with finite diff
void UpdaterHermite::checkFormulas(
	const SpeedInfo& speeds, 
	const DriftParams& sp
) const {
	PairDoub z1(gauss(), gauss()), z2(gauss(), gauss()), z3(gauss(), gauss());
	double u2(abs(gauss())), u3(abs(gauss()));
	SolverParams1Pt pars1{ z1, z2, u2,  
		PairDoub(gauss(), gauss()) , speeds };
	SolverParams2Pt pars{ z1, z2, z3, u2, u3, 
		PairDoub(gauss(), gauss()), 
		PairDoub(gauss(), gauss()), speeds };
	double delta = 1e-5;
	std::vector<double> p1{ gauss(), gauss(),gauss() };
	std::vector<double> grad0_ = calcF1FirstDeriv(p1[0], p1[1], pars1),
		grad1_ = calcF1FirstDeriv(p1[0] + delta, p1[1], pars1),
		grad2_ = calcF1FirstDeriv(p1[0], p1[1] + delta, pars1);


	std::cout << "Function checks \t 2 Point Update: \n "
		<< "Finite difference value vs. derivative value\n"
		<< "dF/da0: " << calcF1(p1[0] + delta, p1[1], pars1) 
		- calcF1(p1[0], p1[1], pars1)
		<< " vs. " << grad0_[0] * delta << std::endl
		<< "dF/da1: " << calcF1(p1[0], p1[1] + delta, pars1) 
		- calcF1(p1[0], p1[1], pars1)
		<< " vs. " << grad0_[1] * delta << std::endl
		<< "d^2F/da0: " << grad1_[0] - grad0_[0]
		<< " vs. " << calcF1SecondDeriv(p1[0], p1[1], pars1)[0] * delta 
		<< std::endl
		<< "d^2F/da1: " << grad2_[1] - grad0_[1]
		<< " vs. " << calcF1SecondDeriv(p1[0], p1[1], pars1)[3] * delta 
		<< std::endl
		<< "d^2F/da0a1: " << grad1_[1] - grad0_[1]
		<< " vs. " << calcF1SecondDeriv(p1[0], p1[1], pars1)[1] * delta 
		<< std::endl 
		<< std::endl;

	std::vector<double>
		grad0 = calcF2FirstDeriv(p1[0], p1[1], p1[2], pars),
		grad1 = calcF2FirstDeriv(p1[0] + delta, p1[1], p1[2], pars),
		grad2 = calcF2FirstDeriv(p1[0], p1[1] + delta, p1[2], pars),
		grad3 = calcF2FirstDeriv(p1[0], p1[1], p1[2] + delta, pars);

	std::cout << "Function checks \t 2 Point Update: \n "
		<< "Finite difference value vs. derivative value\n"
		<< "dF/da0: " << calcF2(p1[0] + delta, p1[1], p1[2], pars) 
		- calcF2(p1[0], p1[1], p1[2], pars)
		<< " vs. " << grad0[0] * delta << std::endl
		<< "dF/da1: " << calcF2(p1[0], p1[1] + delta, p1[2], pars) 
		- calcF2(p1[0], p1[1], p1[2], pars)
		<< " vs. " << grad0[1] * delta << std::endl
		<< "dF/dr: " << calcF2(p1[0], p1[1], p1[2] + delta, pars) 
		- calcF2(p1[0], p1[1], p1[2], pars)
		<< " vs. " << grad0[2] * delta << std::endl
		<< "d^2F/da0: " << grad1[0] - grad0[0]
		<< " vs. " << calcF2SecondDeriv(p1[0], p1[1], p1[2], pars)[0] * delta 
		<< std::endl
		<< "d^2F/da1: " << grad2[1] - grad0[1]
		<< " vs. " << calcF2SecondDeriv(p1[0], p1[1], p1[2], pars)[4] * delta 
		<< std::endl
		<< "d^2F/dr: " << grad3[2] - grad0[2]
		<< " vs. " << calcF2SecondDeriv(p1[0], p1[1], p1[2], pars)[8] * delta 
		<< std::endl
		<< "d^2F/da0a1: " << grad2[0] - grad0[0]
		<< " vs. " << calcF2SecondDeriv(p1[0], p1[1], p1[2], pars)[1] * delta
		<< std::endl
		<< "d^2F/da0r: " << grad1[2] - grad0[2]
		<< " vs. " << calcF2SecondDeriv(p1[0], p1[1], p1[2], pars)[2] * delta 
		<< std::endl
		<< "d^2F/da1r: " << grad2[2] - grad0[2]
		<< " vs. " << calcF2SecondDeriv(p1[0], p1[1], p1[2], pars)[5] * delta 
		<< std::endl;

}