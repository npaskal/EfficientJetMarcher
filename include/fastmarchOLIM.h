#ifndef _FASTMARCHOLIM_H
#define _FASTMARCHOLIM_H
#include <cmath>
#include <tuple>
#include <vector>
#include <string>
#include "structs.h"
#include "fastmarch.h"
#include "binarytree.h"
#include "updater.h"

// Better to save stencils as (xi,yi) than as indexes




class FastMarchOLIM : public FastMarch {
public:
	int update_radius;
	std::vector<int> acceptedNeighborCount;



	FastMarchOLIM(const GridInfo& grid, const Updater& upIn, 
		const RunOptions& opIn, const SpeedInfo& spIn, int update_radius_in) :
		FastMarch(grid, upIn, opIn, spIn), update_radius(update_radius_in),
		acceptedNeighborCount(std::vector<int>(grid.gridLength(),0)) {
		for (int ind = 0; ind < grid.gridLength(); ind++) {
			if (group[ind] == ACCEPTED) {
				group[ind] = FRONT;
			}
		}
		const std::vector<int> neighborIndices{ -grid.nx - 1,-grid.nx,-grid.nx + 1,-1,1,grid.nx - 1,grid.nx,grid.nx + 1 };
		for (int ind = 0; ind < grid.gridLength(); ind++) {
			for (const int k : neighborIndices) {
				if (!invalidIndex(ind + k) && (group[ind + k] == FRONT || group[ind + k] == ACCEPTED) ) {
					++acceptedNeighborCount[ind];
				}
			}
			if (acceptedNeighborCount[ind] == 8)
				group[ind] = ACCEPTED;
			
		}
	}
	void initializeSolution() override;
	void runStep(bool& kill) override;
	void updateNeighbors(int ind) override;
	bool initialPoint(int index) override;
};




#endif
