#pragma once
/*
	\file:		driftFunctions.h
	\brief:		This file declares structures relevant to drift and speed
				information of the quasi-potential problem.
	\author:	Nick Paskal
	\date:		10/27/2021
*/



#include <cmath>
#include <tuple>
#include <vector>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "linalg.h"


// struct containing drift parameters 
struct DriftParams {

	int switchKey; // index indicating switch in master list of drifts

	double a; // parameter 1

	double b; // parameter 2

	// constructor
	DriftParams(int key, double aIn, double bIn) :
		switchKey(key), 
		a(aIn), 
		b(bIn) 
	{}

	// copy constructor
	DriftParams(const DriftParams& spIn) :
		switchKey(spIn.switchKey), 
		a(spIn.a), 
		b(spIn.b) 
	{}
};


// Forward declarations for master drift functions and calcSolution functions
// For now, each is a switch.
// TODO: refactor this

PairDoub masterDrift(const PairDoub& z, const void* sp);
MatrixDoub masterGradDrift(const PairDoub& z, const void* sp);
TensorDoub masterHessDrift(const PairDoub& z, const void* sp);
double masterDivDrift(const PairDoub& z, const void* sp);
double masterSol(const PairDoub& z, const void* sp);
PairDoub masterGradSol(const PairDoub& z, const void* sp);
double masterLaplacianSol(const PairDoub& z, const void* sp);
std::string masterStr(const DriftParams& spdPrms);
MatrixDoub masterDriftLinearization(const void* sp);


// struct containing all relevant speed and drift function pointers
struct SpeedInfo {

	std::string str; // function text string

	const DriftParams& sp; // reference to drift parameters

	// function pointer to drift
	PairDoub(*drift)(
		const PairDoub& z, 
		const void* sp
		);

	// function pointer to grad drift 
	MatrixDoub(*gradDrift)(
		const PairDoub& z, 
		const void* sp
		);

	// function pointer to hessian dirft
	TensorDoub(*hessDrift)(
		const PairDoub& z, 
		const void* sp
		);

	// function pointer to exact calcSolution
	// returns 0 if no analytic calcSolution is available
	double (*calcSolution) (
		const PairDoub& z, 
		const void* sp
	);

	// function poitner to grad calcSolution
	// returns 0 vector if no analytic calcSolution is available
	PairDoub(*calcGradSolution)(
		const PairDoub& z, 
		const void* sp
	);

	// function pointer to laplacian calcSolution
	// returns 0 if no analytic calcSolution is available
	double (*laplaceSolution) (
		const PairDoub& z, 
		const void* sp
	);

	// function poitner to linearized drift around point attractor
	MatrixDoub(*driftLinearization)(
		const void* sp
		);

	// function poitner to divergence of drift
	double(*divDrift)(
		const PairDoub& z, 
		const void* sp
		);


	/*
	The slowness is given by
		s(z,v) = |b(z)| - <b(z),v/|v|>
	This has derivatives
		Ds/Dz = [b(z)/|b(z)| - v/|v| ] Db(z)
		Ds/Dv = -b(z)/|v| + v <b(z),v>/|v|^3

	*/


	// calculate slowness function
	double calcSlowness(
		const PairDoub& z, 
		const PairDoub& v, 
		const void* sp
	) const {
		return drift(z, sp).norm() - dot(drift(z, sp), v) / v.norm();
	}

	// calculate gradient of slowness with respect to current position z
	PairDoub slowGradZ(
		const PairDoub& z, 
		const PairDoub& v, 
		const void* sp
	) const {
		return leftMultiply(masterDrift(z, sp) / masterDrift(z, sp).norm() 
			- v / v.norm(), masterGradDrift(z, sp));
	}

	// calculate gradient of slowness in direction of motion v
	PairDoub slowGradV(
		const PairDoub& z, 
		const PairDoub& v, 
		const void* sp
	) const {
		return  v * dot(masterDrift(z, sp), v) / pow(v.norm(), 3) 
			- masterDrift(z, sp) / v.norm();
	}

	// constructor
	// currently, sets all function pointers equal to master drift
	// which is controllect by a switchKey
	SpeedInfo(
		const DriftParams& params
	) :
		sp(params),
		drift(masterDrift),
		calcSolution(masterSol),
		laplaceSolution(masterLaplacianSol),
		calcGradSolution(masterGradSol),
		gradDrift(masterGradDrift),
		str(masterStr(params)),
		hessDrift(masterHessDrift),
		driftLinearization(masterDriftLinearization),
		divDrift(masterDivDrift) {}
};

