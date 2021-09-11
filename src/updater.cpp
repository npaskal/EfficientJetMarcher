#include "fastmarch.h"
#include "updater.h"
#include "structs.h"
#include <iostream>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cassert>

double gauss()
{
	double eps = 1.e-50, pi2 = 8 * atan(1);
	return sqrt(eps - 2. * log(eps + rand() / (double)RAND_MAX)) * cos(rand() * pi2 / RAND_MAX);
}
// Auxiliary function declarations.
bool kill(double a0, double a1, double r);
PairDoub my_rotate(const PairDoub& vec, double a);
inline double p0(double r) { return 1 - 3 * r * r + 2 * r * r * r; }
inline double p1(double r) { return r * (1 - r) * (1 - r); }
inline double p0p(double r) { return 6 * r * (r - 1); }
inline double p1p(double r) { return 1 - 4 * r + 3 * r * r; }
inline double p0pp(double r) { return 12 * r - 6; }
inline double p1pp(double r) { return -4 + 6 * r; }

double LinearChar::onePointUpdateValue(const SolverParams1Pt& p,
	PairDoub& u1grad, Flags& flags) const {
	// Note the gradient does not depend on the quadrature rule.
	PairDoub v{ p.z1 - p.z2 };
	PairDoub b = p.si.drift(p.z1, &p.si.sp);
	u1grad = b.norm() * v / v.norm() - b;
	flags.lambda = 0;
	flags.updateType = 1;
	return F1(p);
}

double LinearChar::twoPointUpdateValue(const SolverParams2Pt& p,
	PairDoub& u1grad, bool& interiorMinimizer, Flags& flags) const {
	bool interiorSolution = false;
	double rCand = hybridSolver(0, 1, p, interiorSolution);
	// Terminate whenever the solver doesn't return critical points in [0,1].
	if (!interiorSolution) {
		u1grad = PairDoub{ INFINITY,INFINITY };
		interiorMinimizer = false;
		return INFINITY;
	}
	double unew = F2(rCand, p);

	SolverParams1Pt z2params{ p.z1,p.z2,p.u2,p.u2del,p.si },
		z3params{ p.z1,p.z3,p.u3,p.u3del,p.si };

	// To identify if the critical point is a min, we check the endpoints.
	// TODO 2nd derivative test with finite difference. 

	double uEnds = std::min(F1(z2params), F1(z3params));
	if (unew < uEnds) {
		interiorMinimizer = true;
		PairDoub zr{ (1 - rCand) * p.z2 + rCand * p.z3 };
		PairDoub v{ p.z1 - zr };
		PairDoub b = p.si.drift(p.z1, &p.si.sp);
		u1grad = b.norm() * v / v.norm() - b;
		flags.lambda = rCand;
		flags.updateType = 2;
		return unew;
	}
	return INFINITY;
}

double LinearChar::hybridSolver(double a0, double a1, const SolverParams2Pt& p,
	bool& interiorSolution) const {
	double smallVal = 1e-13;
	double tol = 1e-8; int maxIter = 25;
	double b = a0, c = a1;
	double fb = F2p(b, p), fc = F2p(c, p);
	interiorSolution = false;
	if (fb * fc > -smallVal) {
		return INFINITY;
	}
	double a = 0, t = 0, d = 0;
	double fa = 0, fd = 0;
	double dm = 0, df = 0, ds = 0, dd = 0;
	int iter = 0;
	while (iter++ <= maxIter) { // kill after 20 iterations even if no convergence
		if (abs(fc) < abs(fb)) {
			t = c; c = b; b = t;
			t = fc; fc = fb; fb = t;
			a = c; fa = fc;
		}
		if (abs(b - c) <= tol) {
			interiorSolution = true; //indicating that the solver found a minimum.
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
			dd = 0.5 * (static_cast<int64_t>(dm > 0) - static_cast<int64_t>(dm < 0))* tol;
		d = b + dd;
		fd = F2p(d, p);
		if (fd == 0) {
			b = c = d; fb = fc = fd;
			break;
		}
		a = b; b = d;
		fa = fb; fb = fd;
		if (fb * fc > smallVal*smallVal) {
			c = a; fc = fa;
		}
	}
	interiorSolution = true;
	return b;
}

double EndpointLinear::F1(const SolverParams1Pt& p) const {
	PairDoub v = p.z1 - p.z2;
	return p.si.slow(p.z1, v, &p.si.sp) * v.norm() + p.u2;
}
double EndpointLinear::F2(double r, const SolverParams2Pt& p) const {
	double(*slow)(const PairDoub&, const PairDoub&, const void*) = p.si.slow;
	const void* spp = &p.si.sp;
	PairDoub zr = (1 - r) * p.z2 + r * p.z3;
	PairDoub vr = p.z1 - zr;
	double u1 = slow(p.z1, vr, spp) * vr.norm() + (1 - r) * p.u2 + r * p.u3;
	return u1;
}
double EndpointLinear::F2p(double r, const SolverParams2Pt& p) const {
	double(*slow)(const PairDoub&, const PairDoub&, const void*) = p.si.slow;
	const void* sp = &p.si.sp;
	PairDoub(*gradV)(const PairDoub&, const PairDoub&, const void*) = p.si.slowGradV;
	PairDoub zr = (1 - r) * p.z2 + r * p.z3;
	PairDoub vr = p.z1 - zr;
	PairDoub dz = p.z3 - p.z2;
	double hr = vr.norm();
	double u1prime = -slow(p.z1, vr, sp) * dot(vr, dz) / hr
		- hr * dot(gradV(p.z1, vr, sp), dz)
		+ p.u3 - p.u2;
	return u1prime;
}


double MidpointLinear::F1(const SolverParams1Pt& p) const {
	PairDoub v = p.z1 - p.z2;
	int num = 1;
	if (false) {
		return p.si.slow(0.5 * (p.z1 + p.z2), v, &p.si.sp) * v.norm() + p.u2;
	}
	else {
		double u1 =0;
		for (int i = 0; i < num; i++) {
			double rat = (2.0*i+1.0)/(2.0*num);
			u1 += p.si.slow(rat * p.z1 + (1 - rat) * p.z2, v, &p.si.sp) / num;
		}
		u1 *= v.norm();
		u1 += p.u2;
		return u1;
	}
}
double MidpointLinear::F2(double r, const SolverParams2Pt& p) const {
	PairDoub(*b)(const PairDoub&, const void*) = p.si.drift;
	const void* spp = &p.si.sp;
	PairDoub zr = (1 - r) * p.z2 + r * p.z3;
	PairDoub vr = p.z1 - zr;
	PairDoub byr = (1 - r) * b(0.5 * (p.z1 + p.z2), spp) + r * b(0.5 * (p.z1 + p.z3), spp);
	double u1 = byr.norm() * (vr).norm() - dot(byr, vr) + (1 - r) * p.u2 + r * p.u3;
	//return u1;

	int num = 1;
	if (false) {
		return u1;
	}
	else {
		double u1 = 0;
		for (int i = 0; i < num; i++) {
			double rat = (2.0 * i + 1.0) / (2.0 * num);
			u1 += p.si.slow(rat * p.z1 + (1 - rat) * zr, vr, &p.si.sp) / num;
		}
		u1 *= vr.norm();
		u1 += (1 - r) * p.u2 + r * p.u3;
		return u1;
	}

}




double MidpointLinear::F2p(double r, const SolverParams2Pt& p) const {
	double(*slow)(const PairDoub&, const PairDoub&, const void*) = p.si.slow;
	const void* spp = &p.si.sp;
	PairDoub(*gradV)(const PairDoub&, const PairDoub&, const void*) = p.si.slowGradV;
	PairDoub(*gradZ)(const PairDoub&, const PairDoub&, const void*) = p.si.slowGradZ;
	PairDoub(*b)(const PairDoub&, const void*) = p.si.drift;
	PairDoub zr = (1 - r) * p.z2 + r * p.z3;
	PairDoub vr = p.z1 - zr;
	PairDoub dz = p.z3 - p.z2;
	double hr = vr.norm();
	PairDoub yr = (zr + p.z1) / 2;

	PairDoub byr = (1 - r) * b(0.5 * (p.z1 + p.z2), spp) + r * b(0.5 * (p.z1 + p.z3), spp);
	double u1_prime = dot(byr, b(0.5 * (p.z1 + p.z3), spp) - b(0.5 * (p.z1 + p.z2), spp)) * (p.z1 - zr).norm() / byr.norm()
		+ byr.norm() * dot(p.z1 - zr, p.z2 - p.z3) / (p.z1 - zr).norm()
		- dot(b(0.5 * (p.z1 + p.z3), spp) - b(0.5 * (p.z1 + p.z2), spp), p.z1 - zr)
		- dot(byr, p.z2 - p.z3)
		+ p.u3 - p.u2;
	return u1_prime;
}

/*
*******************************************************************************
*******************************************************************************
							Hermite Solver -- Quasipotential
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

PairDoub my_rotate(const PairDoub& vec, double a) {
	return rightMultiply(MatrixDoub(PairDoub(1, -a), PairDoub(a, 1)), vec);
}

double Hermite::onePointUpdateValue(const SolverParams1Pt& p, PairDoub& u1grad, 
	Flags& flags) const {
	bool converges = false;
	std::vector<double> outputs = NewtonSolver2D( 0, 0, p, converges);
	double a0, a1;
	if (!converges) { a0 = 0; a1 = 0; }
	else { a0 = outputs[0], a1 = outputs[1]; }
	double updateVal = INFINITY;
	if (higherEval_1pt) updateVal = F1_higher(a0, a1, p);
	else updateVal = F1(a0, a1, p);

	PairDoub dir(p.z1 - p.z2);
	dir = my_rotate(dir, a1);
	dir = dir / dir.norm();
	u1grad = p.si.drift(p.z1, &p.si.sp).norm() * dir - p.si.drift(p.z1, &p.si.sp);
// 	flags = { 1,a0,a1,0 };
	flags.updateType = 1;
	flags.a0 = a0;
	flags.a1 = a1;
	flags.lambda = 0;
	return updateVal;
}

double Hermite::twoPointUpdateValue(const SolverParams2Pt& p, PairDoub& u1grad, 
	bool& interiorMinimizer, Flags& flags) const {
	double rguess = 0.5; // TODO - implement guess.
	bool converges = false;
	std::vector<double> outputs = NewtonSolver3D(0, 0, rguess, p, converges);
	if (!converges) {
		interiorMinimizer = false;
		return INFINITY;
	}	
	double r = outputs[2], a1 = outputs[1], a0 = outputs[0];

	// kill update if the critical point fell out of the relevant region.
	double angThresh = 3; // Upper bound on acceptable |tan(a)| values.
	if (r > 1 || r < 0 || std::max(abs(a0), abs(a1)) > angThresh) {
		interiorMinimizer = false;
		return INFINITY;
	}
	interiorMinimizer = true;
	// Note: this is only really a critical point.
	// Whether it's a minimizer is handled by calling function.

	PairDoub dir = (p.z1 - (1 - r) * p.z2 - r * p.z3);
	dir = my_rotate(dir, a1);
	dir = dir / dir.norm();
// 	flags = { 2,a0,a1,r };
	flags.updateType = 2;
	flags.a0 = a0;
	flags.a1 = a1;
	flags.lambda = r;
	u1grad = p.si.drift(p.z1, &p.si.sp).norm() * dir - p.si.drift(p.z1, &p.si.sp);
	double updateVal = INFINITY;
	if (higherEval_2pt) updateVal = F2_higher(a0, a1, r, p);
	else updateVal = F2(a0, a1, r, p);
	return updateVal;
}

// Note, using Eigen to solve the system was about twice as fast as the by-hand
// linear system solver.
extern int FAIL, SUCCESS;

std::vector<double> Hermite::NewtonSolver2D(double s0, double t0, 
	const SolverParams1Pt& p, bool &converges) const {
	std::vector<double> sol;
	double tol = 1e-4, smallVal = 1e-5; 	
	int numIter = 20;
	std::vector<double> gradVec = F1p(s0, t0, p);
	std::vector<double> hessVec = F1pp(s0, t0, p);
	assert(gradVec.size() == 2 && hessVec.size() == 4);
	Eigen::Vector2d b(-gradVec[0], -gradVec[1]),
		s(s0, t0),
		sprev(s), 
		sdiff(0, 0);
	Eigen::Matrix2d A;
	A << hessVec[0], hessVec[1], hessVec[2], hessVec[3];
	A.trace();
	int iter = 0;
	double norm = 0;
	while (++iter <= numIter) {
		// TODO: come up with a better condition for checking singularity of A
		if (abs(A.determinant()) / pow(abs(A(0)) + abs(A(3)), 2) < smallVal) {
			converges = false;
			FAIL++;
			return sol;
		}
		Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
		s += x;
		gradVec = F1p(s(0), s(1), p);
		hessVec = F1pp(s(0), s(1), p);
		b(0) = -gradVec[0]; b(1) = -gradVec[1];
		for (int i = 0; i < 4; i++) A(i) = hessVec[i];
		sdiff = s - sprev;
		sprev = s;
		norm = sdiff.norm();
		if (norm < tol) {
			sol = { s(0),s(1) };
			converges = true;
			SUCCESS++;
			return sol;
		}
	}
	converges = false;
	FAIL++;
	return sol;
}

std::vector<double> Hermite::NewtonSolver3D(double s0, double t0, double r0, 
	const SolverParams2Pt& p, bool& converges) const {
	std::vector<double> sol;
	double tol = 1e-5, smallVal = 1e-5; 	int numIter = 20;
	std::vector<double> gradVec = F2p(s0, t0, r0, p);
	std::vector<double> hessVec = F2pp(s0, t0, r0, p);
	assert(gradVec.size() == 3 && hessVec.size() == 9);
	Eigen::Vector3d b(-gradVec[0], -gradVec[1], -gradVec[2]), 
		s(s0, t0, r0), 
		sprev(s), 
		sdiff(0, 0, 0);
	Eigen::Matrix3d A;
	for (int k = 0; k < 9; k++)
		A(k) = hessVec[k];
	int iter = 0;
	double norm = 0;
	while (++iter <= numIter) {
		double a = A.determinant();
		if (abs(a) / pow(abs(A(0)) + abs(A(4)) + abs(A(8)), 3) < smallVal
			|| kill(s(0), s(1), s(2)) ){
			converges = false;
			
			return sol;
		}
		Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
		s += x;
		gradVec = F2p(s(0), s(1), s(2), p);
		hessVec = F2pp(s(0), s(1), s(2), p);
		b(0) = -gradVec[0]; b(1) = -gradVec[1]; b(2) = -gradVec[2];
		for (int k = 0; k < 9; k++)
			A(k) = hessVec[k];
		sdiff = s - sprev;
		sprev = s;
		norm = sdiff.norm();
		if (norm < tol) {
			sol = { s(0),s(1), s(2) };
			converges = true;
			return sol;

		}
	}
	converges = false;

	return sol;
}

bool kill(double a0, double a1, double r) {
	if (r > 1.01 || r < -.01) return true;
	else return false;
}

/*
*******************************************************************************
Beginning of Formulas
*******************************************************************************
*/

double Hermite::F1(double a0, double a1, const SolverParams1Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;
	double (*slow) (const PairDoub & z, const PairDoub & v, 
		const void* sp) = p.si.slow;

	PairDoub dir = z1 - z2;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + z2) + 0.125 * h * (a0 - a1) * w;
	return u2 + h / 6 * (
		slow(z2, my_rotate(dir, a0), &sp) * sqrt(1 + a0 * a0)
		+ 4 * slow(y, my_rotate(dir, -(a0 + a1) / 4), &sp) 
			* sqrt(1 + (a0 + a1) * (a0 + a1) / 16)
		+ slow(z1, my_rotate(dir, a1), &sp) * sqrt(1 + a1 * a1)
		);
}

// 1- point update using Hermite characteristic, but with a finer Simpson's rule
// approximation to the action, using 'numPoints' number of equally spaced
// points between z1 and z2.
double Hermite::F1_higher(double a0, double a1, 
	const SolverParams1Pt& p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;
	double (*slow) (const PairDoub & z, const PairDoub & v, 
		const void* sp) = p.si.slow;

	PairDoub dir = z1 - z2;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + z2) + 0.125 * h * (a0 - a1) * w;

	double u_update = 0;
	int sumWeights = 0;
	for (unsigned k = 0; k < evalPts_1pt; k++) {
		double g = (double)k / (evalPts_1pt - (int)1);
		int weight = 0;
		if (k == 0 || k == evalPts_1pt - 1) weight = 1;
		else if (k % 2 == 0) weight = 2;
		else weight = 4;
		PairDoub z_curr = (1 - g) * z2 + g * z1 + h * Q(a0, a1, g) * w;
		double qp = Qp(a0, a1, g);
		u_update += weight * slow(z_curr, my_rotate(dir, qp), &sp) * sqrt(1 + qp * qp);
		sumWeights += weight;
	}
	return u2 + h * u_update / sumWeights;
}

std::vector<double> Hermite::F1p(double a0, double a1, 
	const SolverParams1Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) = p.si.gradDrift;

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

std::vector<double> Hermite::F1pp(double a0, double a1, 
	const SolverParams1Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2;
	double u2 = p.u2;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) = p.si.gradDrift;
	TensorDoub(*hessDrift) (const PairDoub & z, const void* sp) = p.si.hessDrift;

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
		+ h * (a0 + a1) / (16 * sqrtTerm) * dot(by, rightMultiply(dby, w)) / bymag
		+ h * h / 16 * sqrtTerm / bymag * (
			pow(dbyw.norm(), 2) + dot(by, ddby.sandwich(w, w)) - pow(dot(by / bymag, dbyw), 2)
			)
		+ 1.0 / 4 * bymag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		+ h / 4 * dot(dbyw, w)
		);
	double dFda0a1 = h / 6 * (
		-h * h / 16 * sqrtTerm / bymag * (
			pow(dbyw.norm(), 2) + dot(by, ddby.sandwich(w, w)) - pow(dot(by / bymag, dbyw), 2)
			)
		+ 1.0 / 4 * bymag / pow(sqrtTerm, 3)
		+ h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		);
	double dFda1a1 = h / 6 * (
		drift(z1, &sp).norm() / pow(1 + a1 * a1, 1.5)
		- h * (a0 + a1) / (16 * sqrtTerm) * dot(by, rightMultiply(dby, w)) / bymag
		+ h * h / 16 * sqrtTerm / bymag * (
			pow(dbyw.norm(), 2) + dot(by, ddby.sandwich(w, w)) - pow(dot(by / bymag, dbyw), 2)
			)
		+ 1.0 / 4 * bymag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		- h / 4 * dot(dbyw, w)
		);
	return { dFda0a0, dFda0a1, dFda0a1, dFda1a1 };
}

double Hermite::F2(double a0, double a1, double r, 
	const SolverParams2Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2; double u3 = p.u3;
	PairDoub gradU2 = p.u2del, gradU3 = p.u3del;
	double (*slow) (const PairDoub & z, const PairDoub & v, 
		const void* sp) = p.si.slow;

	PairDoub zlam = (1 - r) * z2 + r * z3;
	PairDoub dir = z1 - zlam;
	double h = dir.norm();
	PairDoub w = dir.perpUnit();
	PairDoub y = 0.5 * (z1 + zlam) + 0.125 * h * (a0 - a1) * w;

	double u0prime = dot(z3 - z2, gradU2), u1prime = dot(z3 - z2, gradU3);
	double H = u2 * p0(r) + u3 * p0(1 - r) + u0prime * p1(r) 
		- u1prime * p1(1 - r);
	return H + h / 6 * (
		slow(zlam, my_rotate(dir, a0), &sp) * sqrt(1 + a0 * a0)
		+ 4 * slow(y, my_rotate(dir, -(a0 + a1) / 4), &sp) 
			* sqrt(1 + (a0 + a1) * (a0 + a1) / 16)
		+ slow(z1, my_rotate(dir, a1), &sp) * sqrt(1 + a1 * a1)
		);
}

double Hermite::F2_higher(double a0, double a1, double  r, 
	const SolverParams2Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2; double u3 = p.u3;
	PairDoub gradU2 = p.u2del, gradU3 = p.u3del;
	double (*slow) (const PairDoub & z, const PairDoub & v, const void* sp) = p.si.slow;

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
	for (unsigned k = 0; k < evalPts_2pt; k++) {
		double g = (double)k / (evalPts_2pt - (int)1);
		int weight = 0;
		if (k == 0 || k == evalPts_2pt - 1) weight = 1;
		else if (k % 2 == 0) weight = 2;
		else weight = 4;
		PairDoub z_curr = (1 - g) * zlam + g * z1 + h * Q(a0, a1, g) * w;
		double qp = Qp(a0, a1, g);
		double utemp = weight * slow(z_curr, my_rotate(dir, qp), &sp) * sqrt(1 + qp * qp);
		u_update += weight * slow(z_curr, my_rotate(dir, qp), &sp) * sqrt(1 + qp * qp);
		sumWeights += weight;
	}
	return H + h * u_update / sumWeights;

}

std::vector<double> Hermite::F2p(double a0, double a1, double r, 
	const SolverParams2Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2, u3 = p.u3;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) = p.si.gradDrift;
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
	double dFdr = u2 * p0p(r) - u3 * p0p(1 - r) + u0prime * p1p(r) + u1prime * p1p(1 - r)
		- dot(z1 - zlam, dx) / (6 * h) * sqrt(1 + a0 * a0) * bzlam.norm()
		+ h / 6 * sqrt(1 + a0 * a0) * dot(bzlam, rightMultiply(dbzlam, dx)) / bzlam.norm()
		- 1.0 / 6 * dot(rightMultiply(dbzlam, dx), my_rotate(dir, a0))
		+ 1.0 / 6 * dot(bzlam, my_rotate(dx, a0))
		- dot(z1 - zlam, dx) / (6 * h) * 4 * sqrtTerm * bymag
		+ h / 3 * sqrtTerm * dot(by, rightMultiply(dby, dx - (a0 - a1) / 4 * dx.perp())) / bymag
		- 1.0 / 3 * dot(rightMultiply(dby, dx - (a0 - a1) / 4 * dx.perp()), my_rotate(dir, -(a0 + a1) / 4))
		+ 2.0 / 3 * dot(by, my_rotate(dx, -(a0 + a1) / 4))
		- dot(z1 - zlam, dx) / (6 * h) * sqrt(1 + a1 * a1) * drift(z1, &sp).norm()
		+ 1.0 / 6 * dot(drift(z1, &sp), my_rotate(dx, a1))
		;
	return { dFda0,dFda1,dFdr };
}

std::vector<double> Hermite::F2pp(double a0, double a1, double r, 
	const SolverParams2Pt &p) const {
	SpeedParams sp = p.si.sp;
	PairDoub z1 = p.z1, z2 = p.z2, z3 = p.z3;
	double u2 = p.u2, u3 = p.u3;
	PairDoub(*drift) (const PairDoub & z, const void* sp) = p.si.drift;
	MatrixDoub(*gradDrift) (const PairDoub & z, const void* sp) = p.si.gradDrift;
	TensorDoub(*hessDrift) (const PairDoub & z, const void* sp) = p.si.hessDrift;
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
	double hpp = delx.norm() * delx.norm() / h - dot(dir, delx) * dot(dir, delx) / (h * h * h);
	MatrixDoub dby = gradDrift(yr, &sp);
	TensorDoub ddby = hessDrift(yr, &sp),
		ddbzr = hessDrift(zr, &sp);
	double sqrtTerm = sqrt(1 + (a0 + a1) * (a0 + a1) / 16);
	double u0prime = dot(z3 - z2, gradU2), u1prime = dot(z3 - z2, gradU3);
	double dFda0a0 = h / 6 * (
		bzr.norm() / pow(1 + a0 * a0, 1.5)
		+ h * (a0 + a1) / (16 * sqrtTerm) * dot(byr, dbyw) / byrmag
		+ h * h / 16 * sqrtTerm / byrmag * (
			pow(dbyw.norm(), 2) + dot(byr, ddby.sandwich(w, w)) - pow(dot(byr / byrmag, dbyw), 2)
			)
		+ 1.0 / 4 * byrmag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		+ h / 4 * dot(dbyw, w)
		);
	double dFda0a1 = h / 6 * (
		-h * h / 16 * sqrtTerm / byrmag * (
			pow(dbyw.norm(), 2) + dot(byr, ddby.sandwich(w, w)) - pow(dot(byr / byrmag, dbyw), 2)
			)
		+ 1.0 / 4 * byrmag / pow(sqrtTerm, 3)
		+ h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		);
	double dFda0r =
		h / 6 * a0 / sqrt(1 + a0 * a0) * dot(bzr, dbzrdx) / bzrmag
		- 1.0 / 6 * dot(dbzrdx, dir.perp())
		+ 1.0 / 6 * (dot(bzr - byr, delx.perp()) + h / 2 * dot(dbyw, my_rotate(delx, -(a0 + a1) / 4)))
		+ hp / 6 * (a0 / sqrt(1 + a0 * a0) * bzrmag + 0.25 * (a0 + a1) / sqrtTerm * byrmag + h / 2 * sqrtTerm * dot(byr, dbyw) / byrmag)
		+ h / 48 * (a0 + a1) / sqrtTerm * dot(byr, dbyq) / byrmag
		+ h * h / 24 * sqrtTerm / byrmag * (dot(dbyw, dbyq) + dot(byr, ddby.sandwich(w, q)))
		- h / 12 * sqrtTerm * dot(byr, rightMultiply(dby, delx.perp())) / byrmag
		- h * h / 24 * sqrtTerm * dot(byr, dbyq) * dot(byr, dbyw) / (byrmag * byrmag * byrmag)
		- h / 24 * dot(ddby.sandwich(w, q), my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 12 * dot(rightMultiply(dby, delx.perp()), my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 12 * dot(dbyq, dir.perp())
		;
	double dFda1a1 = h / 6 * (
		drift(z1, &sp).norm() / pow(1 + a1 * a1, 1.5)
		- h * (a0 + a1) / (16 * sqrtTerm) * dot(byr, rightMultiply(dby, w)) / byrmag
		+ h * h / 16 * sqrtTerm / byrmag * (
			pow(dbyw.norm(), 2) + dot(byr, ddby.sandwich(w, w)) - pow(dot(byr / byrmag, dbyw), 2)
			)
		+ 1.0 / 4 * byrmag / pow(sqrtTerm, 3)
		- h / 16 * dot(ddby.sandwich(w, w), my_rotate(dir, -(a0 + a1) / 4))
		- h / 4 * dot(dbyw, w)
		);
	double dFda1r =
		+1.0 / 6.0 * (-dot(byr, delx.perp()) - h / 2 * dot(dbyw, my_rotate(delx, -(a0 + a1) / 4)))

		+ 1.0 / 6 * dot(drift(z1, &sp), delx.perp())
		+ hp / 6 * (a1 / sqrt(1 + a1 * a1) * drift(z1, &sp).norm() + 0.25 * (a0 + a1) / sqrtTerm * byrmag - h / 2 * sqrtTerm * dot(byr, dbyw) / byrmag)
		+ h / 48 * (a0 + a1) / sqrtTerm * dot(byr, dbyq) / byrmag
		- h * h / 24 * sqrtTerm / byrmag * (dot(dbyw, dbyq) + dot(byr, ddby.sandwich(w, q)))
		+ h / 12 * sqrtTerm * dot(byr, rightMultiply(dby, delx.perp())) / byrmag
		+ h * h / 24 * sqrtTerm * dot(byr, dbyq) * dot(byr, dbyw) / (byrmag * byrmag * byrmag)
		+ h / 24 * dot(ddby.sandwich(w, q), my_rotate(dir, -(a0 + a1) / 4))
		- 1.0 / 12 * dot(rightMultiply(dby, delx.perp()), my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 12 * dot(dbyq, dir.perp())
		;
	double dFdrr =
		u2 * p0pp(r) + u3 * p0pp(1 - r) + u0prime * p1pp(r) - u1prime * p1pp(1 - r)
		+ hp / 6 * sqrt(1 + a0 * a0) * dot(bzr, dbzrdx) / bzrmag
		+ h / 6 * sqrt(1 + a0 * a0) / bzrmag * (dbzrdx.norm() * dbzrdx.norm() + dot(bzr, ddbzr.sandwich(delx, delx)) -
			pow(dot(bzr, dbzrdx) / bzrmag, 2))
		- 1.0 / 6 * dot(ddbzr.sandwich(delx, delx), my_rotate(dir, a0))
		+ 1.0 / 6 * dot(dbzrdx, my_rotate(delx, a0))
		+ 1.0 / 6 * dot(dbzrdx, my_rotate(delx, a0))
		+ 1.0 / 3 * dot(dbyq, my_rotate(delx, -(a0 + a1) / 4))
		+ 1.0 / 6 * hpp * (sqrt(1 + a0 * a0) * bzrmag + 4 * sqrtTerm * byrmag + sqrt(1 + a1 * a1) * drift(z1, &sp).norm())
		+ 1.0 / 6 * hp * (sqrt(1 + a0 * a0) * dot(bzr, dbzrdx) / bzrmag + 2 * sqrtTerm * dot(byr, dbyq) / byrmag)
		+ hp / 3 * sqrtTerm * dot(byr, dbyq) / byrmag
		+ h / 6 * sqrtTerm / byrmag * (dbyq.norm() * dbyq.norm() + dot(byr, ddby.sandwich(q, q)) - dot(byr, dbyq) * dot(byr, dbyq) / (byrmag * byrmag))
		- 1.0 / 6 * dot(ddby.sandwich(q, q), my_rotate(dir, -(a0 + a1) / 4))
		+ 1.0 / 3 * dot(dbyq, my_rotate(delx, -(a0 + a1) / 4))
		;
	return { dFda0a0, dFda0a1,dFda0r,dFda0a1,dFda1a1,dFda1r,dFda0r,dFda1r,dFdrr };
}

void Hermite::checkFormulas(const SpeedInfo& speeds, const SpeedParams& sp) const {
	PairDoub z1(gauss(), gauss()), z2(gauss(), gauss()), z3(gauss(), gauss());
	double u2(abs(gauss())), u3(abs(gauss()));
	SolverParams1Pt pars1{ z1, z2, u2,  PairDoub(gauss(), gauss()) , speeds };
	SolverParams2Pt pars{ z1, z2, z3, u2, u3, PairDoub(gauss(), gauss()), PairDoub(gauss(), gauss()), speeds };
	double delta = 1e-5;
	std::vector<double> p1{ gauss(), gauss(),gauss() };
	std::vector<double> grad0_ = F1p(p1[0], p1[1], pars1),
		grad1_ = F1p(p1[0] + delta, p1[1], pars1),
		grad2_ = F1p(p1[0], p1[1] + delta, pars1);


	std::cout << "Function checks \t 2 Point Update: \n Finite difference value vs. derivative value\n"
		<< "dF/da0: " << F1(p1[0] + delta, p1[1], pars1) - F1(p1[0], p1[1], pars1)
		<< " vs. " << grad0_[0] * delta << std::endl
		<< "dF/da1: " << F1(p1[0], p1[1] + delta, pars1) - F1(p1[0], p1[1], pars1)
		<< " vs. " << grad0_[1] * delta << std::endl
		<< "d^2F/da0: " << grad1_[0] - grad0_[0]
		<< " vs. " << F1pp(p1[0], p1[1], pars1)[0] * delta << std::endl
		<< "d^2F/da1: " << grad2_[1] - grad0_[1]
		<< " vs. " << F1pp(p1[0], p1[1], pars1)[3] * delta << std::endl
		<< "d^2F/da0a1: " << grad1_[1] - grad0_[1]
		<< " vs. " << F1pp(p1[0], p1[1], pars1)[1] * delta << std::endl << std::endl;

	std::vector<double>grad0 = F2p(p1[0], p1[1], p1[2], pars),
		grad1 = F2p(p1[0] + delta, p1[1], p1[2], pars),
		grad2 = F2p(p1[0], p1[1] + delta, p1[2], pars),
		grad3 = F2p(p1[0], p1[1], p1[2] + delta, pars);

	std::cout << "Function checks \t 2 Point Update: \n Finite difference value vs. derivative value\n"
		<< "dF/da0: " << F2(p1[0] + delta, p1[1], p1[2], pars) - F2(p1[0], p1[1], p1[2], pars)
		<< " vs. " << grad0[0] * delta << std::endl
		<< "dF/da1: " << F2(p1[0], p1[1] + delta, p1[2], pars) - F2(p1[0], p1[1], p1[2], pars)
		<< " vs. " << grad0[1] * delta << std::endl
		<< "dF/dr: " << F2(p1[0], p1[1], p1[2] + delta, pars) - F2(p1[0], p1[1], p1[2], pars)
		<< " vs. " << grad0[2] * delta << std::endl
		<< "d^2F/da0: " << grad1[0] - grad0[0]
		<< " vs. " << F2pp(p1[0], p1[1], p1[2], pars)[0] * delta << std::endl
		<< "d^2F/da1: " << grad2[1] - grad0[1]
		<< " vs. " << F2pp(p1[0], p1[1], p1[2], pars)[4] * delta << std::endl
		<< "d^2F/dr: " << grad3[2] - grad0[2]
		<< " vs. " << F2pp(p1[0], p1[1], p1[2], pars)[8] * delta << std::endl
		<< "d^2F/da0a1: " << grad2[0] - grad0[0]
		<< " vs. " << F2pp(p1[0], p1[1], p1[2], pars)[1] * delta << std::endl
		<< "d^2F/da0r: " << grad1[2] - grad0[2]
		<< " vs. " << F2pp(p1[0], p1[1], p1[2], pars)[2] * delta << std::endl
		<< "d^2F/da1r: " << grad2[2] - grad0[2]
		<< " vs. " << F2pp(p1[0], p1[1], p1[2], pars)[5] * delta << std::endl;

}