/*
	\file:		fastMarch.cpp
	\brief:		This file defines base class for Fast Marching methods
				applied to quasi-potential problem.
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
#include <sstream>
#include "binarytree.h"
#include "fastmarch.h"


// constructor
FastMarch::FastMarch(
	const MeshInfo& mesh, 
	UpdaterBase& upIn,		// non-constant because incrementer. TODO: change
	const RunOptions& opIn, 
	const SpeedInfo& spIn
) :
	mesh(mesh),
	updater(upIn),
	options(opIn),
	speeds(spIn),
	considered(&u, &indexInConsidered),
	u(std::vector<double>(mesh.gridLength(), INFINITY)),
	uGrad(std::vector<PairDoub>(
		mesh.gridLength(), 
		PairDoub(INFINITY, INFINITY))
	),
	laplaceU(std::vector<double>(
		mesh.gridLength(), 
		INFINITY)
	),
	group(std::vector<char>(
		mesh.gridLength(), 
		FAR)
	),
	indexInConsidered(std::vector<int>(
		mesh.gridLength(), 
		INT_MAX)
	),
	uOne(std::vector<double>(
		mesh.gridLength(), 
		INFINITY)
	),
	udel_check(std::vector<PairDoub>(
		mesh.gridLength(), 
		PairDoub(INFINITY, INFINITY))
	),
	lastUpdateDist(std::vector<double>(
		mesh.gridLength(), 
		INFINITY)
	),
	minOnePtUpdateVal(std::vector<double>(
		mesh.gridLength(), 
		INFINITY)
	),
	updateType(std::vector<int>(
		mesh.gridLength(), 
		INT_MAX)
	),
	stats()
{
	if (opIn.debugMode) debugInfo = std::vector<UpdateInfo>(
		mesh.gridLength(), 
		UpdateInfo()
		);
}

double FastMarch::getXValue(int index) const {
	return mesh.getXValue(index);
}

double FastMarch::getYValue(int index) const {
	return mesh.getYValue(index);
}

int FastMarch::getIndex(const PairInt& indexPair) const {
	return indexPair.xi + indexPair.yi * mesh.ny;
}

PairDoub FastMarch::getCoords(int index) const {
	return PairDoub{ getXValue(index),getYValue(index) };
}

PairInt FastMarch::getIndexPair(int index) const {
	return PairInt{ getXIndex(index),getYIndex(index) };
}

int FastMarch::getXIndex(int index) const {
	return mesh.getXIndex(index);
}

int FastMarch::getYIndex(int index) const {
	return mesh.getYIndex(index);
}

double FastMarch::getSol(int index) const {
	return u[index];
}

bool FastMarch::isInBoundary(int index) const {
	return mesh.isInBoundary(index);
}

bool FastMarch::isInvalidIndex(int index) const {
	return mesh.isInvalidIndex(index);
}

bool FastMarch::isInvalidIndexPair(PairInt indexPair) const {
	return indexPair.xi < 0 
		|| indexPair.xi >= mesh.nx 
		|| indexPair.yi < 0 
		|| indexPair.yi >= mesh.ny;
}

bool FastMarch::isConsideredEmpty() {
	return (considered.isEmpty());
}

double FastMarch::calcSlowness(int index) {
	return speeds.calcSlowness(getCoords(index), getCoords(index), &speeds.sp);
}

double FastMarch::calcSolution(int index) {
	return speeds.calcSolution(getCoords(index), &speeds.sp);
}

PairDoub FastMarch::calcGradSolution(int index) {
	return speeds.calcGradSolution(getCoords(index), &speeds.sp);
}

double FastMarch::laplaceSolution(int index) {
	return speeds.laplaceSolution(getCoords(index), &speeds.sp);
}

double FastMarch::divDrift(int index) {
	return speeds.divDrift(getCoords(index), &speeds.sp);
}

PairDoub FastMarch::drift(int index) {
	return speeds.drift(getCoords(index), &speeds.sp);
}

PairDoub FastMarch::drift(PairDoub z) {
	return speeds.drift(z, &speeds.sp);
}

// initialize the mesh points near the attractor
void FastMarch::initializeSolution() {
	std::vector<std::pair<double, int>> pairs;
	for (int ind = 0; ind < mesh.gridLength(); ind++) {
		if (isInitialMeshPoint(ind)) {
			u[ind] = computeInitialValue(ind);
			uGrad[ind] = computeInitialGradValue(ind);
			group[ind] = ACCEPTED;
			pairs.push_back(std::make_pair(u[ind], ind));
		}
	}
	/* To enforce causality, we sort the initialized points by u, which may be
	   time-consuming if the initialization region is large.
	   Then we perform the update procedure on the sorted list of initialized
	   Accepted points. 
	*/
	std::sort(pairs.begin(), pairs.end());
	for (auto& p : pairs) {
		accepted.push_back(p.second);
		updateNeighbors(p.second);
	}
}

// check if the given mesh point should be initialized
bool FastMarch::isInitialMeshPoint(int index) {
	double ellipse_umax = .1;
	double boxWidth = 1;
	double xbound = INFINITY;

	// default to 0 or 1, to set initialization region as 8 pt rectangle
	switch (options.initType) {
	default:
	case 0: case 1:
		return (abs(getXValue(index)) < (boxWidth + .5) * mesh.hx 
			&& abs(getYValue(index)) < (boxWidth + .5) * mesh.hy);
	case 2:		// this seems wrong, correct this
		return (computeInitialValue(index) < ellipse_umax 
			&& (
				abs(getXValue(index)) < xbound 
				&& abs(getYValue(index)) < xbound
				)
			|| abs(getXValue(index)) < (boxWidth + .5) * mesh.hx 
			&& abs(getYValue(index)) < (boxWidth + .5) * mesh.hy);
	}
}

// compute the value to initialize the given mesh point with
double FastMarch::computeInitialValue(int index) {
	bool linearApprox(true);
	if (options.initType == 0) {
		// return linearized solution
		MatrixDoub linMat = speeds.driftLinearization(&speeds.sp);
		double B11 = linMat.row1.x, 
			B12 = linMat.row1.y,
			B21 = linMat.row2.x, 
			B22 = linMat.row2.y;
		double aux1 = B21 - B12, 
			aux2 = B11 + B22, 
			aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		double x = getXValue(index), y = getYValue(index);
		return Q11 * x * x + (Q12 + Q21) * x * y + Q22 * y * y;
	}
	else if (options.initType == 1)
		return speeds.calcSolution(getCoords(index), &speeds.sp);
	else
		return 0;
}

// compute value to initialize gradient of solution
PairDoub FastMarch::computeInitialGradValue(int index) {
	bool linearApprox(true);
	if (options.initType == 0) {
		// return linearized solution
		MatrixDoub linMat = speeds.driftLinearization(&speeds.sp);
		double B11 = linMat.row1.x, 
			B12 = linMat.row1.y, 
			B21 = linMat.row2.x, 
			B22 = linMat.row2.y;
		double aux1 = B21 - B12, 
			aux2 = B11 + B22, 
			aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		double x = getXValue(index), y = getYValue(index);
		return PairDoub(
			2 * Q11 * x + (Q12 + Q21) * y, 
			(Q12 + Q21) * y + 2 * Q22 * y
		);
	}
	else if (options.initType == 1)
		return speeds.calcGradSolution(getCoords(index), &speeds.sp);
	else
		return PairDoub(0, 0);
}

// print method specific info to cout
void FastMarch::printIntro() {
	std::cout << "Running fast march with ";
	if (true) std::cout << "8"; else std::cout << "4";
	std::cout << "-point nearest neighbor updates\n\n"
		<< "Speed: \t\ts(x,y) = " << speeds.str << "\n\n"
		<< "Grid: \t \t[" << mesh.xl << ", " << mesh.xr << "] x [" 
		<< mesh.yl << ", " << mesh.yr << "]\n"
		<< "Spacing: \t(nx,ny) = (" << mesh.nx << ", " << mesh.ny << ")\n"
		<< "Quadrature: \t";
	/*switch (updater.quadType) {
	case 'e':
		std::cout 
		<< "Linear characteristics with endpoint quadrature\numNeighbors";
		break;
	case 'm': 
		std::cout 
		<< "Linear characteristics with midpoint rule quadrature\numNeighbors";
		break;
	case 's':
		std::cout 
		<< "Linear characteristics with Simpsons rule quadrature\numNeighbors";
		break;
	case 'h':
		std::cout << "Cubic hermite polynomial characteristics "
		<<< "with Simpson's rule quadrature\numNeighbors";
		break;
	}*/
	std::cout << "\n\n";
}

// compute default stats
void FastMarch::computeStats(double procTime) {
	stats.run_time = procTime;
	stats.err_sup = computeSupError();
	stats.err_rms = computeRMSError();
	stats.err_grad_sup = computeGradSupError();
	stats.err_grad_rms = computeGradRMSError();
	for (int k = 0; k < accepted.size(); k++) {
		int ind = accepted[k];
		if (isInvalidIndex(ind)) continue;
		if (options.debugMode) {
			if (debugInfo[ind].lastUpdateType == 1) stats.num_acc_1pt++;
			else if (debugInfo[ind].lastUpdateType == 2) stats.num_acc_2pt++;
		}
		stats.num_acc++;
	}
}

// print stats to cout
void FastMarch::printStats() {
	std::cout << "Sup error for u:\t\t " << stats.err_sup << std::endl;
	std::cout << "RMS error for u:\t\t " << stats.err_rms << std::endl;
	std::cout << "Sup error for Du:\t\t " << stats.err_grad_sup << std::endl;
	std::cout << "RMS error for Du:\t\t " << stats.err_grad_rms << std::endl;
	std::cout << "Computation time:\t\t " << stats.run_time << std::endl;
	std::cout << "Total 1-point updates:\t " << stats.num_acc_1pt << std::endl;
	std::cout << "Total 2-point updates:\t " << stats.num_acc_2pt << std::endl;
}


// run the Fast March -- main function
void FastMarch::runMarch() {
	initializeSolution();
	int k = 0;
	bool kill = false;
	if (options.verbose) {
		printIntro();
		std::cout << "Completion progress: ";
	}
	clock_t time_s = clock();
	while (!isConsideredEmpty()) {
		runStep(kill);
		if (kill == true)
			break;
		if (options.verbose) {
			if ((k++) % (mesh.gridLength() / 20) == 0)
				std::cout << 5 * k / (mesh.gridLength() / 20) << "%...";
		}
	}
	clock_t time_f = clock();
	double procTime = (double(time_f) - time_s) / CLOCKS_PER_SEC;
	computeStats(procTime);
	computeLaplacian();
	if (options.verbose) {
		std::cout << " Completed.\n\n";
		printStats();
	}	
}

// run a step of Fast March unless the kill condition is met
void FastMarch::runStep(bool &kill3DNewtonSolver) {
	int ind = considered.takeMin();	
	if (isInBoundary(ind)) {
		kill3DNewtonSolver = true;
		return;	
	}	
	group[ind] = ACCEPTED;
	accepted.push_back(ind);
	updateNeighbors(ind);
}

// update the neighbors of the newly minted Accepted point
void FastMarch::updateNeighbors(int index) {
	const std::vector<PairInt>& neighborShifts =  mesh.eightPtNeighborShifts;
	int numNeighbors = neighborShifts.size();
	for (int neighbor = 0; neighbor < numNeighbors; neighbor++) {
		PairInt indexPairUpdate(getIndexPair(index)+ neighborShifts[neighbor]);
		if (isInvalidIndexPair(indexPairUpdate)) 
			continue;
		int index_update = getIndex(indexPairUpdate);
		if (group[index_update] == ACCEPTED)  
			continue;
		/* Remark:
		We choose to perform the 1-point update first before checking for 
		triangle updates.	
		*/	
		onePointUpdate(index_update, index); 
		for (auto m : { -1,1 }) {
			//TODO: explain this
			int ind_updater =  getIndex( getIndexPair(index_update) 
				+ neighborShifts[(numNeighbors+numNeighbors/2+neighbor+m)
				%numNeighbors]);
			/* Remark:
			The 1st numNeighbors eliminates issue with modular arithmetic.
			The numNeighbors/2 flips the shifts to be measured from 
			index_update, whereas originally they were measured in the opposite 
			direction from index. The m (either 1 or -1) ensures that only
			adjacent vertices can be considered for triangles.
			*/
			if (isInvalidIndex(ind_updater)) continue;
			if (group[ind_updater] == ACCEPTED)
				triangleUpdate(index_update, index, ind_updater);
		}
		if (group[index_update] == FAR) {
			group[index_update] = CONSIDERED;
			considered.tree.push_back(index_update);
			indexInConsidered[index_update] = considered.tree.size() - 1;
		}
		considered.reassembleFromTarget(indexInConsidered[index_update]);
	}
}

// update z1_i via one-point update from z2_i
void FastMarch::onePointUpdate(
	int z1_i, 
	int z2_i
) {
	// TODO: get rid of brace initialization
	const SolverParams1Pt params{ 
		getCoords(z1_i),
		getCoords(z2_i),
		u[z2_i],
		uGrad[z2_i],
		speeds 
	};

	PairDoub u1grad(INFINITY,INFINITY);
	Flags flags_new;
	double unew = updater.calcOnePtUpdateValues(params,u1grad, flags_new);
	if (z1_i == trackInd) {
		std::cout << std::endl << "New update: 1 pt\tu=" << unew << "\t";
		(getIndexPair(z2_i) - getIndexPair(z1_i)).print();
		if (unew < u[z1_i]) std::cout << "\tSuccess";
		else std::cout << "\tFail";
		std::cout << std::endl;
	}
	if (unew < u[z1_i]) {
		u[z1_i] = unew;
		uGrad[z1_i] = u1grad;
		updateType[z1_i] = 1;
		if (options.debugMode) {
			debugInfo[z1_i].updateFlags(flags_new);
			debugInfo[z1_i].updateInd(z2_i, z2_i);
		}
	}
	if (unew < minOnePtUpdateVal[z1_i]) 
		minOnePtUpdateVal[z1_i] = unew;
}

// update z1_i via triangle update from z2_i and z3_i
void FastMarch::triangleUpdate(int z1_i, int z2_i, int z3_i) {
	// Checks if this two point update would result in an update distance 
	// larger than the update distance from the previous
	// successful 2-point update.
	double newDist = std::min((getCoords(z1_i) - getCoords(z2_i)).norm(),
		(getCoords(z1_i) - getCoords(z3_i)).norm());
	bool tooFarAway = options.twoPtUpdateMinDist
		&& updateType[z1_i] >= 2
		&& (newDist >= lastUpdateDist[z1_i]);
	if (tooFarAway) 
		return;

	// Run the minimization.
	const SolverParams2Pt params{ getCoords(z1_i),getCoords(z2_i),
		getCoords(z3_i),u[z2_i],u[z3_i], uGrad[z2_i],uGrad[z3_i],speeds };
	bool interiorMinimizer = false;
	PairDoub u1grad(INFINITY, INFINITY);
	Flags flags_new;
	double unew = updater.calcTriangleUpdateValues(params, u1grad,
		interiorMinimizer, flags_new);

	// Checks if the new value is less than a previous 1-point update value.
	// This is designed to filter out other local minima, not corresponding to 
	// the MAP.
	bool wrongLocalMin = options.fakeFilter && unew > minOnePtUpdateVal[z1_i];

	// For debugging.
	if (z1_i == trackInd) {
		std::cout << std::endl << "New update: 2 pt\tu=" << unew << "\t";
		(getIndexPair(z2_i) - getIndexPair(z1_i)).print();
		std::cout << "\tand\t";
		(getIndexPair(z3_i) - getIndexPair(z1_i)).print();
		if (unew < u[z1_i]) std::cout << "\tSuccess";
		else std::cout << "\tFail";
		std::cout << std::endl;
	}
	if ( (interiorMinimizer && !wrongLocalMin && options.twoPtUpdateMinDist) 
	|| (interiorMinimizer && !options.twoPtUpdateMinDist && unew < u[z1_i]  ) 
		) {
		u[z1_i] = unew;
		uGrad[z1_i] = u1grad;
		updateType[z1_i] = 2;
		if (options.debugMode) {
			debugInfo[z1_i].updateFlags(flags_new);
			debugInfo[z1_i].updateInd(z2_i, z3_i);
		}

		// Compute linear distance between z1 and zr
		lastUpdateDist[z1_i] = (getCoords(z1_i) - ( (1- flags_new.lambda) 
			*getCoords(z2_i) + flags_new.lambda*getCoords(z3_i))).norm();
	}
}

// compute finite difference Laplacian for each mesh point
void FastMarch::computeLaplacian() {
	for (unsigned x = 0; x < mesh.nx ; x++)
		for (unsigned y = 0; y < mesh.ny ; y++) {
			int index = y * mesh.nx + x;
			if (group[index] != ACCEPTED 
				&& group[index] != FRONT) 
				continue;
			int index_right = y * mesh.nx + x + 1;
			int index_left = y * mesh.nx + x - 1;
			int index_up = (y + 1) * mesh.nx + x;
			int index_down = (y - 1) * mesh.nx + x;
			double dxx = 0, dyy = 0;

			// compute dxx part of Laplacian
			if (!isInvalidIndex(index_right) 
				&& !isInvalidIndex(index_left) 
				&& group[index_right] >= FRONT 
				&& group[index_left] >= FRONT
				)	// interior point
				dxx = (uGrad[index_right].x - uGrad[index_left].x) 
					/ (2 * mesh.hx);
			else if (!isInvalidIndex(index_right) 
				&&  group[index_right >= FRONT]) // edge case
				dxx = (uGrad[index_right].x - uGrad[index].x) / (mesh.hx);
			else if (!isInvalidIndex(index_left) 
				&& group[index_left] >= FRONT) // edge case
				dxx = (uGrad[index].x - uGrad[index_left].x) / mesh.hx;
			else
				dxx = INFINITY;

			// compute dyy part of Laplacian
			if (!isInvalidIndex(index_up) 
				&& !isInvalidIndex(index_down) 
				&& group[index_up] >= FRONT 
				&& group[index_down] >= FRONT)
				dyy = (uGrad[index_up].y - uGrad[index_down].y) 
				/ (2 * mesh.hy);
			else if (!isInvalidIndex(index_up) && group[index_up >= FRONT])
				dyy = (uGrad[index_up].y - uGrad[index].y) / mesh.hy;
			else if (!isInvalidIndex(index_down) && group[index_down >= FRONT])
				dyy = (uGrad[index].y - uGrad[index_down].y) / mesh.hy;
			else
				dyy = INFINITY;
			laplaceU[index] = dxx + dyy;
		}
}

// writes solution to txt file
void FastMarch::writeSolutionToTXT() {
	// text file containing all solution values (Considered & Accepted)
	std::ofstream outFileAllVals("outputs/solutionAll.txt");
	bool partial = false;
	outFileAllVals 
		<< "x,y,U,Err(U),Err(DU),Last Update Type, DU_x,DU_y" 
		<< std::endl;
	for (unsigned x=0; x < mesh.nx; x++) 
		for (unsigned y = 0; y < mesh.ny; y++) {
			int index = y * mesh.nx + x;
			outFileAllVals
				<< getXValue(index) << ","
				<< getYValue(index) << ","
				<< u[index] << ","
				<< uGrad[index].x << ","
				<< uGrad[index].y << ","
				<< abs(u[index] - calcSolution(index)) << ","
				<< (uGrad[index] - calcGradSolution(index)).norm() << ","
				<< debugInfo[index].lastUpdateType << std::endl;
		}
	// tack on some problem information to end of txt file
	// useful for Matlab scripts
	outFileAllVals 
		<< mesh.nx << "," 
		<< mesh.ny << "," 
		<< mesh.nx*mesh.ny << "," 
		<< speeds.sp.switchKey << "," 
		<< speeds.sp.a << "," 
		<< speeds.sp.b 
		<< std::endl;
	outFileAllVals.close();

	// text file containing ONLY Accepted solution values
	std::ofstream outFileAccepted("outputs/solutionAccepted.txt");
	outFileAccepted 
		<< "x,y,U,Err(U),Err(DU),Last Update Type, DU_x,DU_y" 
		<< std::endl;
	for (unsigned x = 0; x < mesh.nx; x++)
		for (unsigned y = 0; y < mesh.ny; y++) {
			int index = y * mesh.nx + x;
			if (group[index] == ACCEPTED || group[index] == FRONT) {
				outFileAccepted
					<< getXValue(index) << ","
					<< getYValue(index) << ","
					<< u[index] << ","
					<< uGrad[index].x << ","
					<< uGrad[index].y << ","
					<< abs(u[index] - calcSolution(index)) << ","
					<< (uGrad[index] - calcGradSolution(index)).norm() << ","
					<< debugInfo[index].lastUpdateType
					//<< "," 	<< laplaceU[index] << "," 
					//<< getCoords( debugInfo[index].ind2).x << "," 
					//<< getCoords(debugInfo[index].ind2).y << "," 
					//<< getCoords(debugInfo[index].ind3).x << "," 
					//<< getCoords(debugInfo[index].ind3).y << ","
					//<< debugInfo[index].lam << "," 
					//<< debugInfo[index].a0 << "," 
					//<< debugInfo[index].a1 << "," 
					//<< debugInfo[index].rank 
					<< std::endl;
			}
			else
				outFileAccepted
				<< getXValue(index) << ","
				<< getYValue(index) << ","
				<< "nan,"
				<< "nan,"
				<< "nan,"
				<< "nan,"
				<< "nan,"
				<< "nan"
				<< std::endl;
		}
	outFileAccepted 
		<< mesh.nx << "," 
		<< mesh.ny << "," 
		<< mesh.nx * mesh.ny << "," 
		<< speeds.sp.switchKey << "," 
		<< speeds.sp.a << "," 
		<< speeds.sp.b 
		<< std::endl;
	outFileAccepted.close();
}

// compute statistics in format comparable with results in OLIM paper
void FastMarch::computeStatsOLIMComparison(
	std::vector<double>* stats
) {
	*stats = std::vector<double>(8, 0);
	double max_err = 0, 
		max_grad_err = 0, 
		max_u = 0;
	double rms_sum = 0, 
		rms_grad_sum = 0, 
		rms_u_sum = 0;
	double max_laplace_err = 0, 
		rms_laplace_sum = 0;
	unsigned num_points = 0, 
		num_laplace = 0;
	for (int x = 0; x < mesh.nx; x++) {
		for (int y = 0; y < mesh.ny; y++) {
			int index = getIndex(PairInt(x,y));
			if (group[index] == ACCEPTED) {
				double abs_err = abs(u[index] - calcSolution(index));
				double abs_grad_err = 
					(uGrad[index] - calcGradSolution(index)).norm();
				double abs_laplace_err = 
					abs(laplaceU[index] - laplaceSolution(index));
				max_err = std::max(abs_err, max_err);
				max_grad_err = std::max(abs_grad_err, max_grad_err);

				max_u = std::max(u[index], max_u);
				rms_sum += abs_err * abs_err;
				rms_grad_sum += abs_grad_err * abs_grad_err;
				rms_u_sum += u[index] * u[index];
				if (abs_laplace_err < INFINITY) {
					max_laplace_err = 
						std::max(
							abs_laplace_err, 
							max_laplace_err
						);
					rms_laplace_sum += abs_laplace_err * abs_laplace_err;
					num_laplace++;
				}


				num_points++;				
			}
		}
	}
	(*stats)[0] = max_err;
	(*stats)[1] = sqrt(rms_sum / num_points);
	(*stats)[2] = max_err / max_u;
	(*stats)[3] = sqrt(rms_sum) / sqrt(rms_u_sum);
	(*stats)[4] = max_grad_err;
	(*stats)[5] = sqrt(rms_grad_sum / num_points);	
	(*stats)[6] = max_laplace_err;
	(*stats)[7] = sqrt(rms_laplace_sum / num_points);
}

// compute Sup error of solution
double FastMarch::computeSupError() {
	double max_error = 0;
	for (int x = 0; x < mesh.nx; x++) {
		for (int y = 0; y < mesh.ny; y++) {
			int index = y * mesh.ny + x;
			if ((group[index] == ACCEPTED)
				&& abs(u[index] - calcSolution(index)) > max_error
				) {
				max_error = abs(u[index] - calcSolution(index));
			}
		}
	}
	return max_error;
}

// compute RMS error of solution
double FastMarch::computeRMSError() {
	double rms_sum = 0;
	unsigned num_points = 0; 
	for (int ind = 0; ind < mesh.gridLength() ; ind++) {
		double u_exact = calcSolution(ind);
		if (group[ind] == ACCEPTED) {
			rms_sum += (u[ind] - u_exact) * (u[ind] - u_exact);
			num_points++;
		}
	}
	return sqrt(rms_sum / num_points);
}

// compute Sup error of gradient solution
double FastMarch::computeGradSupError() {
	double max_error = 0;
	for (int x = 0; x < mesh.nx; x++) {
		for (int y = 0; y < mesh.ny; y++) {
			int index = y * mesh.ny + x;
			double curr_err = (uGrad[index] - calcGradSolution(index)).norm();
			if ( (
				group[index] == ACCEPTED)	
				&& curr_err > max_error
				) 
				max_error = curr_err;
		}
	}
	return max_error;
}



double FastMarch::computeGradRMSError() {
	double sumSquares = 0;
	unsigned num_points = 0;
	for (int ind = 0; ind < mesh.gridLength(); ind++) {
		if (group[ind] == ACCEPTED) {
			sumSquares += (uGrad[ind] - calcGradSolution(ind)).normsq();
			num_points++;
		}
	}
	return sqrt(sumSquares / num_points);
}

// write detailed point by point last update information to txt file
// for Accepted points
void FastMarch::writeDebugInfotToTXT() {
	std::ofstream outDebug("outputs/debug.txt");
	outDebug << std::setprecision(12); 

	outDebug
		<< "rank,i1,e1,e2,e3,ge1,x1,y1,u1,gu1x,gu1y,i2,x2,y2,u2,gu2x,gu2y,i3,"
		<< "x3,y3,u3,gu3x,gu3y,type,lambda,a0,a1,dist" 
		<< std::endl;
	int rank = 0;
	int sizeLimit = mesh.gridLength();
	sizeLimit = 50000;
	for (unsigned k = 0; k < accepted.size(); k++) {


		int ind1 = accepted[k];
		if (isInitialMeshPoint(ind1)) continue;
		if (rank > sizeLimit) break;
		rank++;
		UpdateInfo inf = debugInfo[ind1];
		int ind2 = inf.ind2, ind3 = inf.ind3;
		outDebug
			<< rank << ","
			<< ind1 << ","
			<< u[ind1] - calcSolution(ind1) << ","
			<< u[ind2] - calcSolution(ind2) << ","
			<< u[ind3] - calcSolution(ind3) << ","
			<< (uGrad[ind1] -calcGradSolution(ind1) ).norm() << ","
			<< getCoords(ind1).x << ","
			<< getCoords(ind1).y << ","
			<< u[ind1] << ","
			<< uGrad[ind1].x << ","
			<< uGrad[ind1].y << ","
			<< inf.ind2 << ","
			<< getCoords(ind2).x << ","
			<< getCoords(ind2).y << ","
			<< u[ind2] << ","
			<< uGrad[ind2].x << ","
			<< uGrad[ind2].y << ","
			<< inf.ind3 << ","
			<< getCoords(ind3).x << ","
			<< getCoords(ind3).y << ","
			<< u[ind3] << ","
			<< uGrad[ind3].x << ","
			<< uGrad[ind3].y << ","
			<< inf.lastUpdateType<< ","
			<< inf.lam << ","
			<< inf.a0 << ","
			<< inf.a1 << ","
			<< (getCoords(ind1) - getCoords(ind2)).norm() / mesh.hx
			<< std::endl;
	}
	// tack on some problem information to end of txt file
// useful for Matlab scripts
	outDebug << mesh.xl << "," 
		<< mesh.xr << "," 
		<< mesh.nx << "," 
		<< mesh.yl << "," 
		<< mesh.yr << "," 
		<< mesh.ny << "," 
		<< speeds.sp.a << "," 
		<< speeds.sp.b << "," 
		<< speeds.sp.switchKey <<
		std::endl;
	outDebug.close();


}



// performs exhaustive triangle update search if last update is one-pt
void FastMarch::fixFailures(int ind, int num) {
	PairInt ipair = getIndexPair(ind);
	int x0 = -num, y0 = 0;
	int xprev = -num + 1, yprev = 1;
	int xinc = 1, yinc = -1;

	for (int k = 0; k < 4 * num; k++) {
		PairInt p2{ x0,y0 };
		PairInt p3{ xprev,yprev };
		int i2 = getIndex(ipair + p2);
		int i3 = getIndex(ipair + p3);


		if (!mesh.isInvalidIndex(i2) 
			&& !mesh.isInvalidIndex(i3) 
			&& group[i2] == ACCEPTED 
			&& group[i3] == ACCEPTED)
			triangleUpdate(ind, i2, i3);
		xprev = x0, yprev = y0;
		x0 += xinc; y0 += yinc;
		if (abs(x0) == num)
			xinc *= -1;
		else if (abs(y0) == num)
			yinc *= -1;
	}
}

// performs exhaustive triangle update search if last update is one-pt.
// This is the slower and more thorough version.
// Use only for debugging purposes
void FastMarch::fixFailures_soupedUp(int ind, int num) {
	PairInt ipair = getIndexPair(ind);
	int x0 = -num, y0 = 0;
	int xprev = -num + 1, yprev = 1;
	int xinc = 1, yinc = -1;

	for (int k = 0; k < 4 * num; k++) {
		PairInt p2{ x0,y0 };
		PairInt p3{ xprev,yprev };
		int i2 = getIndex(ipair + p2);


		for (const PairInt& shift : mesh.eightPtNeighborShifts) {
			int i3 = getIndex(getIndexPair(i2) + shift);
			if (!mesh.isInvalidIndex(i2) 
				&& !mesh.isInvalidIndex(i3) 
				&& group[i2] == ACCEPTED 
				&& group[i3] == ACCEPTED
				)
				triangleUpdate(ind, i2, i3);
		}
		xprev = x0, yprev = y0;
		x0 += xinc; y0 += yinc;
		if (abs(x0) == num)
			xinc *= -1;
		else if (abs(y0) == num)
			yinc *= -1;
	}
}

//ARCHIVED -- DO NOT USE
// shoot maps from attractor to x0
// TODO - cleanup and export to Matlab
double FastMarch::shootMaps(PairDoub x0) {
	PairDoub origin(0, 0);
	double step = mesh.hx/10; // Set to h/10.

	PairDoub xcur(x0);
	PairDoub gucur(0, 0);
	int runs = 10000;
	std::ofstream sol_out("Outputs/mapShooter.csv");
	int k = 0;
	std::vector<PairDoub> path;
	path.push_back(xcur);
	bool rk = true;
	bool print = true;
	if (rk)
		while ((xcur - origin).norm() > 2*step) // Set to radius of termination
		{
			gucur = interpolate(xcur);
			xcur = xcur - step * (gucur + drift(xcur)) 
				/ (gucur + drift(xcur)).norm();
			path.push_back(xcur);
			sol_out << xcur.x << "," << xcur.y << std::endl;
			k++;
			if (k >= runs)
				break;
		}
	else 
		while ((xcur - origin).norm() > 2 * step) {
			PairDoub gu1 = interpolate(xcur);
			PairDoub k1 = (gu1 + drift(xcur)) / (gu1 + drift(xcur)).norm();
			PairDoub gu2 = interpolate(xcur - k1 * step / 2);
			PairDoub k2 = (gu2 + drift(xcur)) / (gu2 + drift(xcur)).norm();
			PairDoub gu3 = interpolate(xcur - k2 * step / 2);
			PairDoub k3 = (gu3 + drift(xcur)) / (gu3 + drift(xcur)).norm();
			PairDoub gu4 = interpolate(xcur - k3 * step);
			PairDoub k4 = (gu4 + drift(xcur)) / (gu4 + drift(xcur)).norm();
			xcur = xcur - step / 6.0 * (k1 + 2*k2 + 2*k3 + k4);
			path.push_back(xcur);
			sol_out << xcur.x << "," << xcur.y << std::endl;

			k++;
			if (k >= runs)
				break;
		}
	for (int j = 0; j < 100; j++) {
	//	path[j].print(); std::cout << std::endl;
	}
	sol_out.close();
	std::ifstream infile;
	infile.open("SourceFiles/master.csv");
	std::vector<PairDoub> points;
	std::string line;

	while (std::getline(infile, line)) {
		std::size_t pos = line.find(",");     
		PairDoub point(
			std::stod(line.substr(0, pos)), 
			std::stod(line.substr(pos + 1, line.size()))
		); 
		points.push_back(point);
	}

	std::vector<double> errs;
	double minerr = 0;
	for (int i = 0; i < path.size(); i++) {
		double curmin = INFINITY;
		int jmin = 0;
		for (int j = 0; j < points.size(); j++) {
			double curval = (path[i] - points[j]).normsq();
			if ( curval< curmin) {
				curmin = curval;
				jmin = j;
			}
		}
		if (jmin == 0 || jmin == points.size() - 1) {
			errs.push_back(0);
			continue;
		}
		
		PairDoub p2_cur = points[jmin], 
			p2_next = points[jmin + 1], 
			p2_prev = points[jmin - 1];
		PairDoub p1 = path[i];

		if (p1.norm()  < 0.1) {
			errs.push_back(0);
			continue;
		}


		double dotprod = dot(p1 - p2_prev, p2_cur - p2_prev);
		double d12 = (p1 - p2_prev).normsq() 
			- dotprod * dotprod / (p2_cur - p2_prev).normsq();
		d12 = sqrt(d12);
		double dotprod2 = dot(p1 - p2_cur, p2_next - p2_cur);
		double d23 = (p1 - p2_cur).normsq() 
			- dotprod2 * dotprod2 / (p2_next - p2_cur).normsq();
		d23 = sqrt(d23);
		errs.push_back(std::min(d12, d23));
		minerr = std::max(minerr, std::min(d12, d23));
	}

	// Compute pre-factor
	double F = 0,  Fprev= 0;
	double sum = 0;
	double tot = 0;
	int numFin = 0;
	for (int i = path.size()-1; i >1 ; i--) {
		PairDoub mid = (path[i - 1] + path[i]) / 2;
		F = speeds.divDrift(mid, &speeds.sp);
		F += 0.5 * interpLaplacian(mid); 
		if (F > 100000) F = Fprev;
		else {
			Fprev = F;
			tot += abs(F);
			numFin++;
		}
		sum += F * (path[i-1] - path[i ]).norm();
	}

	//return minerr;
	return exp(sum)-1;
	//return tot / numFin;

}

// ARCHIVED -- DO NOT USE
PairDoub FastMarch::interpolate(PairDoub xcur) {
	PairInt bl(
		floor((xcur.x - mesh.xl) / mesh.hx), 
		floor((xcur.y - mesh.yl) / mesh.hy)
	);
	PairInt br(bl.xi + 1, bl.yi), 
		tl(bl.xi, bl.yi + 1), 
		tr(bl.xi + 1, bl.yi + 1);
	int bli = getIndex(bl), 
		bri = getIndex(br), 
		tli = getIndex(tl), 
		tri = getIndex(tr);
	PairDoub gu_bl(uGrad[bli]), 
		gu_br(uGrad[bri]), 
		gu_tl(uGrad[tli]), 
		gu_tr(uGrad[tri]);
	PairDoub gu(0, 0);
	double gux = 0, guy = 0;
	MatrixDoub xmat(
		PairDoub(gu_bl.x, gu_tl.x), 
		PairDoub(gu_br.x, gu_tr.x)
	);
	MatrixDoub ymat(
		PairDoub(gu_bl.y, gu_tl.y), 
		PairDoub(gu_br.y, gu_tr.y)
	);
	double x1 = getCoords(bli).x, 
		x2 = getCoords(bri).x;
	double y1 = getCoords(bli).y, 
		y2 = getCoords(tli).y;
	PairDoub vec1(x2 - xcur.x, xcur.x - x1);
	PairDoub vec2(y2 - xcur.y, xcur.y - y1);
	double fac = 1.0 / ((x2 - x1) * (y2 - y1));
	gux = fac * dot(vec1, rightMultiply(xmat, vec2));
	guy = fac * dot(vec1, rightMultiply(ymat, vec2));
	gu.x = gux;
	gu.y = guy;
	//std::cout << "xi:\t" << bli;
	return gu;

}


// ARCHIVED -- DO NOT USE
// compute the Laplacian U value at input point by interpolating solution
double FastMarch::interpLaplacian(const PairDoub &xcur) {
	PairInt bl(
		floor((xcur.x - mesh.xl) / mesh.hx), 
		floor((xcur.y - mesh.yl) / mesh.hy)
	);
	PairInt br(bl.xi + 1, bl.yi), 
		tl(bl.xi, bl.yi + 1), 
		tr(bl.xi + 1, bl.yi + 1);
	int bli = getIndex(bl), 
		bri = getIndex(br), 
		tli = getIndex(tl), 
		tri = getIndex(tr);
	double gu_bl(laplaceU[bli]), 
		gu_br(laplaceU[bri]), 
		gu_tl(laplaceU[tli]), 
		gu_tr(laplaceU[tri]);
	double gu(0);
	MatrixDoub xmat(
		PairDoub(gu_bl, gu_tl), 
		PairDoub(gu_br, gu_tr)
	);
	double x1 = getCoords(bli).x, x2 = getCoords(bri).x;
	double y1 = getCoords(bli).y, y2 = getCoords(tli).y;
	PairDoub vec1(x2 - xcur.x, xcur.x - x1);
	PairDoub vec2(y2 - xcur.y, xcur.y - y1);
	double fac = 1.0 / ((x2 - x1) * (y2 - y1));
	gu = fac * dot(vec1, rightMultiply(xmat, vec2));
	return gu;
}