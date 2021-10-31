#pragma once
/*
	\file:		fastMarchEJM.h
	\brief:		This file declares the fastMarchEJM class and information
				about the update stencil.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#include <cmath>
#include <vector>
#include "updaterBase.h"
#include "fastmarch.h"
#include "driftFunctions.h"


struct StencilData;

// class for Efficient Jet Marching method
class FastMarchEJM : public FastMarch {
public:
	const StencilData &stencilData; // stencil setup information

	// collection of stencils, indexed by drift angles
	std::vector<std::vector<PairInt> > stencilCollection; 

	// vector containing stencilCollection index for each mesh point
	std::vector<int> stencilIndex;

	// the only EJM constructor
	FastMarchEJM(
		const MeshInfo& mesh,
		UpdaterBase& upIn,
		const RunOptions& opIn,
		const SpeedInfo& spIn,
		const StencilData& stencilDataIn
	);

	// run a step of EJM unless the kill condition is met
	void runStep(bool& kill3DNewtonSolver) override;
	
	// construct the collection of stencils
	void createStencils();

	// update the neighbors of the newly minted Accepted point
	void updateNeighbors(int ind) override;

	// initialize the mesh points near the attractor
	void initializeSolution() override;

	// write the stencil corresponding to mesh point ind to txt file
	void writeStencilToTXT(int ind);

	// write solution and stencil to txt file
	void writeSolutionToTXT() override;

};

struct StencilData {
	int stencilType;	// algorithm for constructing stencil

	double alpha;		// parameter alpha in construction of Acute stencil

	int stencilCutoff;		// max radius of stencil (in multiples of h)

	int numThetaCells;	// number of theta bins

	// set to true to use neighbors of 8 pt neighborhood for z3 candidate
	// set to false to use stencil neighborhood as z3 candidate
	bool boxUpdate;

	// collection of drift angle thetas to use for stencil bins
	std::vector<double> theta_res;

	// Constructor
	StencilData(
		double alphaIn, 
		int maxStencIn, 
		int binsIn, 
		bool boxUpIn, 
		int stencilTypeIn) 
		:
		alpha(alphaIn), 
		stencilCutoff(maxStencIn), 
		numThetaCells(binsIn), 
		boxUpdate(boxUpIn),
		stencilType(stencilTypeIn)
	{
		createStencilResolution();
	}

	// Sets theta angles for BUBBLE stencil.
	// TODO - refactor, this should not be in Stencil Data
	void createStencilResolution() {
		// Note this puts the midpoint angle in between pi/2^k and pi/2^{k+1}
		int kmax = 8;
		for (int k = kmax; k >= 1; k--) {
			theta_res.push_back(M_PI / pow(2, k));
			if (k > 1) 
				theta_res.push_back(M_PI * 3.0 / pow(2, k + 1));
		}
		for (int k = 1; k <= kmax; k++) {
			if (k > 1) 
				theta_res.push_back(2 * M_PI - M_PI * 3.0 / pow(2, k + 1));
			theta_res.push_back(2 * M_PI - M_PI / pow(2, k));

		}

	}
};


// enume for algorithm for creating stencils
enum StencilType {
	ACUTE = 0,  // as in Mirebeau paper
	BUBBLE = 1,
	OBLONG = 2	// as in our paper
};