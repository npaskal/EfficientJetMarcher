#ifndef _FASTMARCHASR_H
#define _FASTMARCHASR_H
#include <cmath>
#include <tuple>
#include <vector>
#include <string>
#include "structs.h"
#include "fastmarch.h"
#include "binarytree.h"
#include "updater.h"






class FastMarchASR : public FastMarch{
public:
	// stencil information
	const StencilData &stencilData;
	std::vector<std::vector<PairInt> > stencil;
	std::vector<int> stencilIndex;

	// necessary information for 2-pt update rules.



	FastMarchASR(const GridInfo& grid, 
		const Updater& upIn, 
		const RunOptions& opIn, 
		const SpeedInfo &spIn, 
		const StencilData& stencilDataIn) :
		FastMarch(grid, upIn,opIn,spIn), 
		stencilData(stencilDataIn), 
		stencilIndex(std::vector<int>(grid.gridLength(), INT_MAX))

	{
		if (stencilData.bubble) createStencils_bubble_new();
		else createStencils_Mirebeau();		
	} 
	void runStep(bool& kill) override;
	void createStencils_Mirebeau();
	void createStencils_bubble();
	void createStencils_bubble_new();

	void updateNeighbors(int ind) override;
	//void twoPointUpdate(int ind1, int ind2, int ind3) override;
	//void onePointUpdate(int ind1, int ind2) override;
	void initializeSolution() override;
	void writeStencil(int ind);
	void writeToTXT() override;

};




#endif
