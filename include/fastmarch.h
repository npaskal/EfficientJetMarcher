#pragma once
/*
	\file:		fastMarch.h
	\brief:		This file declares base class for Fast Marching methods
				applied to quasi-potential problem
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#include <cmath>
#include <vector>
#include <string>
#include "binarytree.h"
#include "updaterBase.h"
#include "driftFunctions.h"



struct MeshInfo;

struct RunOptions;

struct UpdateInfo;

// Summary statistics for FastMarch object
struct SummaryStats {
	double err_sup; // sup error of solution U

	double err_rms; // RMS error of solution U

	double err_grad_sup; // sup error of grad U

	double err_grad_rms; // RMS error of grad U

	double run_time;	// total runtime of FastMarch

	int num_acc;		// number of Accepted points

	int num_acc_1pt;	// number of Accepted points from 1-pt udpates

	int num_acc_2pt;	// number of Accepted poitns from triangle updates

	// default constructor
	SummaryStats() : 
		err_sup(INFINITY), 
		err_rms(INFINITY), 
		err_grad_sup(INFINITY), 
		err_grad_rms(INFINITY),
		run_time(INFINITY), 
		num_acc(0), 
		num_acc_1pt(0), 
		num_acc_2pt(0) {}

};



// Class for running fast march method for quasi-potential problem
// TODO: generalize to anisotropic eikonal solvers
class FastMarch {
public:
	///////////////////////////////////////////////////////////////////////////
	// Core Problem Information:
	///////////////////////////////////////////////////////////////////////////


	// Information about computation mesh
	const MeshInfo &mesh;			

	// Collection of one-point update and triangle update functions
	UpdaterBase &updater;											

	// Struct of miscellaneous settings
	const RunOptions &options;			

	// Collection of drift and derivative formulas.
	// Contains exact quasi-potential formulas, if available
	const SpeedInfo& speeds;


	///////////////////////////////////////////////////////////////////////////
	// Primary Computed Information:
	///////////////////////////////////////////////////////////////////////////


	// Vector of computed quasi-potential (U) values
	std::vector<double> u;

	// Vector of computed quasi-potential gradient (DU) values
	std::vector<PairDoub> uGrad;

	// Vector of computed Laplacian quasi-potential values
	std::vector<double> laplaceU;

	// Vector containing current group (Unk., Cons., Acc.) of each mesh point
	std::vector<char> group;

	// Vector of detailed update information for each mesh point
	std::vector<UpdateInfo> debugInfo;		

	// Final summary statistics
	SummaryStats stats;


	///////////////////////////////////////////////////////////////////////////
	// Supplementary Information:
	///////////////////////////////////////////////////////////////////////////
	
	// Binary tree containing Considered points
	BinaryTree considered;

	// List of Accepted points
	std::vector<int> accepted;

	// List of each mesh points current index in the Considered binary tree
	std::vector<int> indexInConsidered;

	// Vector of distance of most recent update
	std::vector<double> lastUpdateDist;

	// Minimum length of all processed one-point updates for each mesh point
	std::vector<double> minOnePtUpdateVal;

	// Type of most recent update (1 or 2) for each Cons. or Acc. mesh point
	std::vector<int> updateType;

	
	///////////////////////////////////////////////////////////////////////////
	// Debugging Information:
	///////////////////////////////////////////////////////////////////////////


	// Index of mesh point to provide thorough debugging tracking
	int trackInd = INT_MAX;


	int failSafeCalls = 0;




	// Temporary items for overcoming structural issues in algorithm design.
	
	std::vector<PairDoub> udel_check;

	std::vector<double> uOne;

	std::vector<PairDoub> z1, z2;

	std::vector<double> lam, a0, a1;

public:

	// only constructor
	FastMarch(
		const MeshInfo& mesh, 
		UpdaterBase& upIn,
		const RunOptions& opIn, 
		const SpeedInfo& spIn
	);


	///////////////////////////////////////////////////////////////////////////
	// Mesh related getters
	///////////////////////////////////////////////////////////////////////////


	// get x coordinate from mesh point index
	double getXValue(int index) const;

	// get y coordinate from mesh point index
	double getYValue(int index) const;

	// get mesh point index from mesh point PairIndex
	int getIndex(const PairInt& indexPair) const;

	// get (x,y) coordinates from mesh point index
	PairDoub getCoords(int index) const;

	// get mesh point PairIndex from mesh point index
	PairInt getIndexPair(int index) const;

	// get mesh point xIndex
	int getXIndex(int index) const;

	// get mesh point yIndex
	int getYIndex(int index) const;

	// get computed calcSolution at mesh point index
	double getSol(int index) const;

	// check if given mesh point is in boundary
	bool isInBoundary(int index) const;

	// check if given mesh point index is valid
	bool isInvalidIndex(int index) const;

	// check if given PairIndex is valid
	bool isInvalidIndexPair(PairInt indexPair) const;


	///////////////////////////////////////////////////////////////////////////
	// Primary Fast March functions
	///////////////////////////////////////////////////////////////////////////


	// run the Fast March -- main function
	void runMarch();

	// run a step of EJM unless the kill condition is met
	virtual void runStep(bool &killCondition);

	// update the neighbors of the newly minted Accepted point
	virtual void updateNeighbors(int ind);

	// update ind1 via one-point update from ind2
	virtual void onePointUpdate(int ind1, int ind2);

	// update ind1 via triangle update from ind2 and ind3
	virtual void triangleUpdate(int ind1, int ind2, int ind3);

	// initialize the mesh points near the attractor
	virtual void initializeSolution();

	// performs exhaustive triangle update search if last update is one-pt
	void fixFailures(int ind, int num);

	// performs exhaustive triangle update search if last update is one-pt.
	// This is the slower and more thorough version.
	// Use only for debugging purposes
	void fixFailures_soupedUp(int ind, int num);

	// check if the given mesh point should be initialized
	virtual bool isInitialMeshPoint(int index);

	// compute the value to initialize the given mesh point with
	double computeInitialValue(int index);

	// compute value to initialize gradient of solution
	PairDoub computeInitialGradValue(int index);

	// check if the Considered list is empty
	bool isConsideredEmpty();


	///////////////////////////////////////////////////////////////////////////
	// Speed-related function shortcuts
	///////////////////////////////////////////////////////////////////////////


	// calculates the slowness at given mesh point index
	double calcSlowness(int index);

	// calculates exact solution at given mesh point index
	double calcSolution(int index);

	// calculates exact gradient of solution at given mesh point index
	PairDoub calcGradSolution(int index);

	// calculates exact Laplacian of solution at given mesh point index
	double laplaceSolution(int index);
	
	// calculate drift of divergence at given index
	double divDrift(int index);

	// calculate drift at given index
	PairDoub drift(int index);

	// calculate drift at given point
	PairDoub drift(PairDoub z);


	///////////////////////////////////////////////////////////////////////////
	// Output related functions
	///////////////////////////////////////////////////////////////////////////


	// writes solution to txt file
	virtual void writeSolutionToTXT();

	// write detailed point by point last update information to txt file
	// for Accepted points
	void writeDebugInfotToTXT();

	// compute Sup error of solution
	double computeSupError();

	// compute RMS error of solution
	double computeRMSError();

	// compute Sup error of gradient solution
	double computeGradSupError();

	// compute RMS error of gradient solution
	double computeGradRMSError();

	// compute statistics in format comparable with results in OLIM paper
	void computeStatsOLIMComparison(std::vector<double> *stats);

	// print method specific info to cout
	virtual void printIntro();

	// print summary statistics
	void printStats();

	// compute default stats
	void computeStats(double procTime);

	// compute finite difference Laplacian for each mesh point
	void computeLaplacian();


	///////////////////////////////////////////////////////////////////////////
	// Auxiliary functions
	///////////////////////////////////////////////////////////////////////////


	//ARCHIVED -- DO NOT USE
	// shoot maps from attractor to x0
	// TODO - cleanup and export to Matlab
	double shootMaps(PairDoub x0);

	// ARCHIVED -- DO NOT USE
	PairDoub interpolate(PairDoub xcur);
	
	// ARCHIVED -- DO NOT USE
	// compute the Laplacian U value at input point by interpolating solution
	double interpLaplacian(const PairDoub &pair);	
};

// Object containing information about mesh and boundary
// Assumes coordinates are shifted so that attractor is at (0,0)
// For now, ONLY uniform rectangular mesh supported
struct MeshInfo {
public:
	int nx; // number of mesh points in x direction

	double xl; // x coordinate of rectangle left edge

	double xr; // x coordinate of rectangle right edge

	int ny; // number of mesh points in y direction

	double yl; // y coordinate of rectangle bottom edge

	double yr; // y coordinate of rectangle top edge

	double hx; // uniform x spacing

	double hy; // uniform y spacing

	// constant vector containing integer index shifts to get to
	// 8 point nearest neighbors
	const std::vector<PairInt> eightPtNeighborShifts;

	// constant vector containing integer index shifts to get to 
	// 4 point nearest neighbors
	const std::vector<PairInt> fourPtNeighborShifts;


	// list constructor
	MeshInfo(
		int nx_i, 
		double xl_i, 
		double xr_i, 
		int ny_i, 
		double yl_i, 
		double yr_i) :
		nx(nx_i), 
		xl(xl_i), 
		xr(xr_i), 
		ny(ny_i), 
		yl(yl_i), 
		yr(yr_i),
		hx(double(xr_i - xl_i) / (static_cast<double>(nx_i) - 1)),
		hy(double(yr_i - yl_i) / (static_cast<double>(ny_i) - 1)),
		eightPtNeighborShifts(
			{ PairInt(1,0),PairInt(1,1),	
			PairInt(0,1), PairInt(-1,1),	
			PairInt(-1,0),PairInt(-1,-1),
			PairInt(0,-1),PairInt(1,-1) 
			}
		),
		fourPtNeighborShifts(
			{ PairInt(1,0),PairInt(0,1),
			PairInt(-1,0),PairInt(0,-1) 
			}
		) 
	{}

	// copy constructor
	MeshInfo(
		const MeshInfo& in
	) : 
		nx(in.nx), 
		xl(in.xl), 
		xr(in.xr), 
		ny(in.ny), 
		yl(in.yl), 
		yr(in.yr),
		hx(in.hx), 
		hy(in.hy),
		eightPtNeighborShifts(in.eightPtNeighborShifts),
		fourPtNeighborShifts(in.fourPtNeighborShifts) 
	{}

	// get x coordinate from index
	double getXValue(int index) const {
		return xl + hx * (index % ny);
	} 

	// get y coordinate from index
	double getYValue(int index) const {
		return yl + hy * (index / ny);
	}

	// get x index from index
	int getXIndex(int index) const {
		return index % ny;
	}

	// get y index from index
	int getYIndex(int index) const {
		return index / ny;
	}

	// compute number of mesh points
	int gridLength() const {
		return nx * ny;
	}

	// check if mesh point is in boundary
	bool isInBoundary(int index) const {
		return (
			(getXIndex(index) == 0) 
			|| (getXIndex(index) == nx - 1) 
			|| (getYIndex(index) == 0) 
			|| (getYIndex(index) == ny - 1)
			);
	}

	// check if mesh point index is valis
	bool isInvalidIndex(int index) const {
		return ((index < 0) 
			|| (index >= gridLength())
			);
	}
};

// struct containing miscellaneous options for Fast March
struct  RunOptions {
	bool verbose;	// flag for printing cout output

	// flag for storing detailed debug informaion
	// set to true, will take much more memory
	bool debugMode;		

	// set to true if triangle updates are to be accepted on the basis of their
	// hypotenuse length, rather than their proposed value
	bool twoPtUpdateMinDist;

	// set to true to run fail-safe on Accepted points without triangle updates
	bool failSafe;

	// set to true to require triangle updates to propose values of U smaller
	// than all previous one-point updates, in order to adopt them
	bool fakeFilter;

	// switch key for specifying how to initialize
	int initType;

	// set to true if you want one-point update to be conducted for z2 and z3
	// set to false if you only want it to be conducted for z2
	bool z3one = false;
};


// struct containing information about a particular update
// to be used as detailed Debug information
struct UpdateInfo {
	int ind1, // index of point z1 in question

		ind2, // index of point z2 of final-update of z1

		ind3; // index of point z3 of final-update of z1, set to z2 if one-pt

	PairDoub z1, // coordinates of point in question

		z2, // coordinates of point z2 of final-update of z1

		z3; // coords of point z3 of final-update of z1, set to z2 if one-pt

	double u1, // computed solution of z1

		u2,		// computed solution of z2

		u3;		// comptued solution of z3

	PairDoub u1del, // computed grad solution at z1
		
		u2del,	// computed grad solution at z2
		
		u3del;		// computed grad solution at z3

	double lam,		// lambda value of final update of z1

		a0,			// a0 value of final update of z1

		a1;			// a1 value of final update of z1

	int lastUpdateType; // type of last update of z1, i.e. 1 or 2

	int rank = 0;		// rank in the Accepted list

	// update parameter values 
	void updateFlags(const Flags& flags_in) {
		lastUpdateType = flags_in.updateType;
		lam = flags_in.lambda;
		a0 = flags_in.a0;
		a1 = flags_in.a1;
	}

	//update z2 and z3 indices
	void updateInd(int ind2_in, int ind3_in) {
		ind2 = ind2_in;
		ind3 = ind3_in;
	}

};


// Algorithm enum
enum algorithm {
	FM = 0,
	EJM = 1,
	OLIM = 2
};

// Group list enum
enum group {
	FAR = 0,
	CONSIDERED = 1,
	FRONT = 2,
	ACCEPTED = 3
};

// initialization shape enum
enum initShape {
	ELLIPSE = 0,
	BOX = 1
};



