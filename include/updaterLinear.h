#pragma once
/*
	\file:		updaterLinear.h
	\brief:		This file declares class for updater using linear MAPs.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include "linalg.h"
#include "updaterBase.h"

// purely virtual base class for updaters using linear characteristics
// i.e. linear approximations to the MAP
class UpdaterLinearBase : public UpdaterBase {
public:
	// compute one-point update value of u and grad u
	double calcOnePtUpdateValues(
		const SolverParams1Pt& p,				// input parameters
		PairDoub& u1grad,						// output: gradient U
		Flags& flags							// output: other parameter vals
	)  override;

	// compute triangle update value of u and grad u
	double calcTriangleUpdateValues(
		const SolverParams2Pt& p,				// input parameters
		PairDoub& u1grad,						// output: gradient U
		bool& interiorMinimizer,				// output: flag for inter. min.
		Flags& flags							// outpur: other parameter vals
	)  override;

private:	
	// find a staionary point of F2 for lambda between x_left and x_right
	double calcStationaryPointF2(
		double x_left, 
		double x_right, 
		const SolverParams2Pt& p, 
		bool& interiorSolution
	) const;

	// compute one-point update F1
	virtual double calcF1(
		const SolverParams1Pt& p
	) const = 0;

	// compute two-point update F2 for given lambda value.
	// i.e. update z1 from (1-lambda)*z2 + lambda*z3
	virtual double calcF2(
		double lambda, 
		const SolverParams2Pt& p
	) const = 0;

	// compute derivative of two-point update F2 with respect to lambda
	virtual double calcF2FirstDeriv(
		double r, 
		const SolverParams2Pt& p
	) const = 0;
	
};



// Class for computing updates using linear characteristics and an
// ENDPOINT quadrature rule for computing the action integral
class UpdaterLinearEndpoint : public UpdaterLinearBase {
private:
	// compute one-point update F1
	double calcF1(
		const SolverParams1Pt& p
	) const override;

	// compute two-point update F2 for given lambda value
	// i.e. update z1 from (1-lambda)*z2 + lambda*z3
	double calcF2(
		double lambda, const SolverParams2Pt& p
	) const override;

	// compute derivative of two-point update calcF2 with respect to lambda
	double calcF2FirstDeriv(
		double lambda, 
		const SolverParams2Pt& p
	) const override;
};



// Class for computing updates using linear characteristics with a
// MIDPOINT quadrature rule for computing the action integral
class UpdaterLinearMidpoint : public UpdaterLinearBase {
private:
	// compute one-point update F1
	double calcF1(
		const SolverParams1Pt& p
	) const override;

	// compute two-point update F2 for given lambda value
	// i.e. update z1 from (1-lambda)*z2 + lambda*z3
	double calcF2(
		double lambda, 
		const SolverParams2Pt& p
	) const override;

	// compute derivative of two-point update calcF2 with respect to lambda
	double calcF2FirstDeriv(
		double lambda, 
		const SolverParams2Pt& p
	) const override;
};

