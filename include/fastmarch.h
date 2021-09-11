#ifndef FASTMARCH_H
#define FASTMARCH_H
#include <cmath>
#include <vector>
#include <string>
#include "structs.h"
#include "binarytree.h"
#include "updater.h"




class FastMarch {
public:
	const GridInfo &grid;					
	const Updater &updater;																
	const RunOptions &options;				
	const SpeedInfo& speeds;



	std::vector<UpdateInfo> debugInfo;		
	SummaryStats stats;

	std::vector<double> u;
	std::vector<PairDoub> uGrad;
	std::vector<double> laplaceU;
	std::vector<char> group;	

	// These vectors are for filtering purpose.
	std::vector<double> lastUpdateDist;
	std::vector<double> minOnePtUpdateVal;
	std::vector<int> updateType;

	std::vector<int> indexInConsidered;
	GridPointTree considered;
	std::vector<int> accepted;	


	// Temporary items for overcoming structural issues in algorithm design.
	
	std::vector<PairDoub> udel_check;
	std::vector<double> uOne;

	int failSafe = 0;
	int trackInd = INT_MAX;


	// Temporary for debugging.
	std::vector<PairDoub> z1, z2;
	std::vector<double> lam, a0, a1;

public:
	FastMarch(const GridInfo& grid, const Updater& upIn,  
		const RunOptions &opIn, const SpeedInfo & spIn) :
		grid(grid),
		updater(upIn),
		options(opIn),
		speeds(spIn),
		considered(&u, &indexInConsidered), 
		u(std::vector<double>(grid.gridLength(), INFINITY)),
		uGrad(std::vector<PairDoub>(grid.gridLength(), PairDoub( INFINITY, INFINITY)) ),
		laplaceU(std::vector<double>(grid.gridLength(), INFINITY)),
		group(std::vector<char>(grid.gridLength(), FAR)), 
		indexInConsidered(std::vector<int>(grid.gridLength(), INT_MAX)), 		
		uOne(std::vector<double>(grid.gridLength(), INFINITY)),
		udel_check(std::vector<PairDoub>(grid.gridLength(), PairDoub(INFINITY, INFINITY))),
		lastUpdateDist(std::vector<double>(grid.gridLength(), INFINITY)),
		minOnePtUpdateVal(std::vector<double>(grid.gridLength(), INFINITY)),
		updateType(std::vector<int>(grid.gridLength(), INT_MAX)),
	//	z1(std::vector<PairDoub>(grid.gridLength(), PairDoub(INFINITY, INFINITY))),
	//	z2(std::vector<PairDoub>(grid.gridLength(), PairDoub(INFINITY, INFINITY))),
	//	lam(std::vector<double>(grid.gridLength(), INFINITY)),
	//	a0(std::vector<double>(grid.gridLength(), INFINITY)),
	//	a1(std::vector<double>(grid.gridLength(), INFINITY)),
		stats()
	{
		if (opIn.debugMode) debugInfo = std::vector<UpdateInfo>(grid.gridLength(), UpdateInfo());
	}		


	double getXValue(int index) const { return grid.getXValue(index); }
	double getYValue(int index) const { return grid.getYValue(index); }
	int getIndex(const PairInt& indexPair) const { return indexPair.xi + indexPair.yi * grid.ny; } 
	PairDoub getCoords(int index) const { return PairDoub{ getXValue(index),getYValue(index) }; }
	PairInt getIndexPair(int index) const { return PairInt{ getXIndex(index),getYIndex(index) }; }
	int getXIndex(int index) const { return grid.getXIndex(index); }
	int getYIndex(int index) const { return grid.getYIndex(index); }
	double getSol(int index) const { return u[index]; }
	bool inBoundary(int index) const { return grid.inBoundary(index); }
	bool invalidIndex(int index) const { return grid.invalidIndex(index); }
	bool invalidIndexPair(PairInt indexPair) const {
		return indexPair.xi < 0 || indexPair.xi >= grid.nx || indexPair.yi < 0 || indexPair.yi >= grid.ny;
	}
	// test function
	void runMarch();
	virtual void runStep(bool &kill);
	virtual void updateNeighbors(int ind);
	virtual void initializeSolution();
	virtual bool initialPoint(int index);
	double initialValue(int index);
	bool consideredEmpty() { return (considered.isEmpty()); }
	double slow(int index) { return speeds.slow(getCoords(index),getCoords(index),&speeds.sp); }
	double solution(int index) { return speeds.solution(getCoords(index), &speeds.sp); }
	PairDoub gradSolution(int index) { return speeds.gradSolution(getCoords(index), &speeds.sp); }
	double laplaceSolution(int index) { return speeds.laplaceSolution(getCoords(index), &speeds.sp); }
	double divDrift(int index) { return speeds.divDrift(getCoords(index), &speeds.sp); }
	PairDoub drift(int index) { return speeds.drift(getCoords(index), &speeds.sp); }
	PairDoub drift(PairDoub z) { return speeds.drift(z, &speeds.sp); }
	virtual void writeToTXT();
	double computeSupError();
	double computeRMSError();
	double computeGradSupError();
	double computeGradRMSError();
	void computeMashaStats(std::vector<double> *stats);
	virtual void twoPointUpdate(int ind1, int ind2, int ind3);
	virtual void onePointUpdate(int ind1, int ind2);
	PairDoub initialGradU(int index);
	virtual void printIntro();
	void printStats();
	void writeDebug();
	void computeStats(double procTime);
	void fixFailures(int ind, int num);
	void fixFailures_soupedUp(int ind, int num);

	void writeMAP(const PairDoub& z);

	double shootMaps(PairDoub x0);
	PairDoub interpolate(PairDoub xcur);
	PairDoub interpolate_cub(PairDoub xcur);
	void printDebugConsole();
	void computeLaplacian();

	double interpLaplacian(const PairDoub &pair);



	
};


#endif
