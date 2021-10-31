#pragma once
/*
	\file:		updaterBase.h
	\brief:		This file declares base class for updater that computes
				one-point and triangle updates.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include "linalg.h"
#include "driftFunctions.h"

struct SolverParams1Pt;

struct SolverParams2Pt; 

struct Flags; 

// Purely virtual base class for computing one-point and triangle updates 
class UpdaterBase {
public:

	virtual ~UpdaterBase() {}

	// compute one-point update value of u and grad u
	virtual double calcOnePtUpdateValues(
		const SolverParams1Pt& params, 
		PairDoub& u1grad, Flags& flags)  
		=0;

	// compute triangle update value of u and grad u
	virtual double calcTriangleUpdateValues(
		const SolverParams2Pt& params, 
		PairDoub& u1grad, 
		bool& interiorMinimizer, 
		Flags& flags)  
		=0;

protected:

	int m_numOnePtCalls = 0; // number of one-point update calls

	int m_numTriangleCalls = 0; // number of triangle update calls

	int m_numOnePtFailures = 0; // number of one-point update failures

	int m_numOnePtSuccesses = 0; // number of one-point update successes
};

// Information required to perform a triangle update of point 
// z1 from points z2 and z3
struct SolverParams2Pt {
	PairDoub z1; // Coordinates of point z1 to be updated

	PairDoub z2; // Coordinates of updating point z2

	PairDoub z3; // Coordinates of updating point z3

	double u2; // Solution value at updating point z2

	double u3; // Solution value at updating point z3

	PairDoub u2del; // Gradient of calcSolution at updating point z2

	PairDoub u3del; // Gradient of calcSolution at updating point z3

	const SpeedInfo& si; // Object containing drift information
};

// Information required to perform a one-point update of point 
// z1 from point z2
struct SolverParams1Pt {
	PairDoub z1; // Coordinates of point z1 to be updated

	PairDoub z2; // Coordinates of updating point z2

	double u2; // Solution value at updating point z2

	PairDoub u2del; // Gradient of calcSolution at updating point z2

	const SpeedInfo& si; // Object containing drift information
};

// Structure containing output values of a one-point or triangle update
struct Flags {
	// Type of update
	// 0 = none, 1 = one-point, 2 = triangle
	int updateType{ 0 };

	double a0{ 0 }; // Parameter a0 of triangle updates, for UpdaterHermite only

	double a1{ 0 }; // Parameter a1, of triangle updates, for UpdaterHermite only

	double lambda{ 0 }; // Parameter lambda of triangle updates

};

