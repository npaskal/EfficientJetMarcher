#pragma once
/*
	\file:		fastMarchOLIM.h
	\brief:		This file declares the fastMarchOLIM class.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include <cmath>
#include <vector>
#include "fastmarch.h"
#include "updaterBase.h"
#include "driftFunctions.h"



// class for OLIM methods for computing quasi-potential
class FastMarchOLIM : public FastMarch {
public:
	int updateRadius; // search radius for neighbors

	// number of Accepted neighbors for each mesh point
	std::vector<int> acceptedNeighborCount; 
	
	// constructor
	FastMarchOLIM(
		const MeshInfo& mesh,  
		UpdaterBase& upIn, 
		const RunOptions& opIn, 
		const SpeedInfo& spIn, 
		int update_radius_in
	) :
		FastMarch(mesh, upIn, opIn, spIn), 
		updateRadius(update_radius_in),
		acceptedNeighborCount(std::vector<int>(mesh.gridLength(),0)) 
	{
		for (int ind = 0; ind < mesh.gridLength(); ind++) {
			if (group[ind] == ACCEPTED) {
				group[ind] = FRONT;
			}
		}
		const std::vector<int> neighborIndices{ 
			-mesh.nx - 1,
			-mesh.nx,
			-mesh.nx + 1,
			-1,
			1,
			mesh.nx - 1,
			mesh.nx,
			mesh.nx + 1 };
		for (int ind = 0; ind < mesh.gridLength(); ind++) {
			for (const int k : neighborIndices) {
				if (!isInvalidIndex(ind + k) 
					&& (group[ind + k] == FRONT 
						|| group[ind + k] == ACCEPTED) ) {
					++acceptedNeighborCount[ind];
				}
			}
			if (acceptedNeighborCount[ind] == 8)
				group[ind] = ACCEPTED;
			
		}
	}

	// initialize the mesh points near the attractor
	void initializeSolution() override;

	// run a step of Fast March OLIM unless the kill condition is met
	void runStep(bool& kill3DNewtonSolver) override;

	// update the neighbors of the newly minted Accepted point
	void updateNeighbors(int ind) override;

	// check if the given mesh point should be initialized
	bool isInitialMeshPoint(int index) override;
};



