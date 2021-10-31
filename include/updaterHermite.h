#pragma once
/*
	\file:		updaterHermite.h
	\brief:		This file declares class for Hermite updater.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include <vector>
#include "updaterBase.h"
#include "driftFunctions.h"

// Updater class using Hermite cubic approximations to the MAP and
// Hermite cubic interpolations of U
class UpdaterHermite : public UpdaterBase {
public:
	// compute one-point update value of u and grad u
	double calcOnePtUpdateValues(const SolverParams1Pt& p, PairDoub& u1grad,
		Flags& flags)  override;

	// compute triangle update value of u and grad u
	double calcTriangleUpdateValues(const SolverParams2Pt& p, PairDoub& u1grad,
		bool& interiorMinimizer, Flags& flags)  override;

	// double check math is correct. Check derivatives agree with finite diff
	void checkFormulas(const SpeedInfo& speeds, const DriftParams& sp) const;

private:
	// compute one-point update F1, given entry and exit angles a0, a1
	double calcF1(
		double a0, 
		double a1, 
		const SolverParams1Pt& p
	) const;

	// compute first-order derivatives of F1 (2D vector)
	std::vector<double> calcF1FirstDeriv(
		double a0, 
		double a1,
		const SolverParams1Pt& p
	) const;

	// compute second-order derivatives of F1 (as 4D vector)
	std::vector<double> calcF1SecondDeriv(
		double a0, 
		double a1,
		const SolverParams1Pt& p
	) const;

	// compute two-point update F2, given lambda, a0, a1
	double calcF2(
		double a0, 
		double a1, 
		double lambda, 
		const SolverParams2Pt& p
	) const;

	// compute first-order derivatives of F2 (3D vector)
	std::vector<double> calcF2FirstDeriv(
		double a0, 
		double a1, 
		double r,
		const SolverParams2Pt& p
	) const;

	// compute second-order derivatives of F2 (as 9D vector)
	std::vector<double> calcF2SecondDeriv(double a0, 
		double a1, 
		double r,
		const SolverParams2Pt& p
	) const;

	// use Newton's method to find a stationary point of F1
	std::vector<double> calcStationaryPointF1(
		double s0, 
		double t0,
		const SolverParams1Pt& p,
		bool& converges
	);

	// use Newton's method to find a stationary point of F2 
	std::vector<double> calcStationaryPointF2(
		double a0, 
		double a1, 
		double lam,
		const SolverParams2Pt& p, 
		bool& converges
	);

	// compute one-point update using higher-order Simpsons quadrature
	double calcF1Refined(
		double a0, 
		double a1, 
		const SolverParams1Pt& p
	) const;

	// compute triangle update using higher-order Simpsons quadrature
	double calcF2Refined(
		double a0, 
		double a1, 
		double r,
		const SolverParams2Pt& p
	) const;

	   
	// Use higher order quadrature for final u calculations: one-pt updates
	bool m_refinedQuadratureOnePt = true;

	// Use higher order quadrature for final u calculations: triangle updates
	bool m_refinedQuadratureTriangle = true;

	// Number of nodes for final u calculations: one-pt updates
	unsigned m_refinedNumOnePtNodes = 13;

	// Nubmer of nodes for final u calculations: triangle updates
	unsigned m_refinedNumTriangleNodes = 13;
};

