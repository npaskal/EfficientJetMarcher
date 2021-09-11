
#include <cmath>
#include <tuple>
#include <vector>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "structs.h"
#include "binarytree.h"
#include "fastmarch.h"
#include <sstream>




//	Run the initialization procedure.
//	All initialized points are added to ACCEPTED and their neighbors updated.
void FastMarch::initializeSolution() {
	std::vector<std::pair<double, int>> pairs;
	for (int ind = 0; ind < grid.gridLength(); ind++) {
		if (initialPoint(ind)) {
			u[ind] = initialValue(ind);
			uGrad[ind] = initialGradU(ind);
			group[ind] = ACCEPTED;
			pairs.push_back(std::make_pair(u[ind], ind));
		}
	}
	// To enforce causality, we sort the initialized points by u, which may be
	// time-consuming if the initialization region is large.
	// Then we perform the update procedure on the sorted list of initialized
	// Accepted points.
	std::sort(pairs.begin(), pairs.end());
	for (auto& p : pairs) {
		accepted.push_back(p.second);
		updateNeighbors(p.second);
	}
}

/*
	Return true if given index corresponds to a mesh point that should be
	initialized. 
	This depends on value of options.initType
*/
bool FastMarch::initialPoint(int index) {
	//double ellipse_umax = .00005;
	double ellipse_umax = .1;
	double boxWidth = 1;
	double xbound = INFINITY;
	switch (options.initType) {
	default:
	case 0: case 1:
		return (abs(getXValue(index)) < (boxWidth + .5) * grid.hx && abs(getYValue(index)) < (boxWidth + .5) * grid.hy);
	case 2:
		//return  (speeds.solution(getCoords(index), &speeds.sp) < ellipse_umax);
		return (initialValue(index) < ellipse_umax && (abs(getXValue(index)) < xbound && abs(getYValue(index)) < xbound)
			|| abs(getXValue(index)) < (boxWidth + .5) * grid.hx && abs(getYValue(index)) < (boxWidth + .5) * grid.hy); // / (grid.nx / 17.0));
	}
}

/*
	Return initialization value of u at index.
*/
double FastMarch::initialValue(int index) {
	bool linearApprox(true);
	if (options.initType == 0) {
		MatrixDoub linMat = speeds.driftLinearization(&speeds.sp);
		double B11 = linMat.row1.x, B12 = linMat.row1.y, B21 = linMat.row2.x, B22 = linMat.row2.y;
		double aux1 = B21 - B12, aux2 = B11 + B22, aux = aux1 * aux1 + aux2 * aux2;
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
		return speeds.solution(getCoords(index), &speeds.sp);
	else
		return 0;
}

/*
	Return initialization value of grad u at index.
*/
PairDoub FastMarch::initialGradU(int index) {
	bool linearApprox(true);
	if (options.initType == 0) {
		MatrixDoub linMat = speeds.driftLinearization(&speeds.sp);
		double B11 = linMat.row1.x, B12 = linMat.row1.y, B21 = linMat.row2.x, B22 = linMat.row2.y;
		double aux1 = B21 - B12, aux2 = B11 + B22, aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		double x = getXValue(index), y = getYValue(index);
		return PairDoub(2 * Q11 * x + (Q12 + Q21) * y, (Q12 + Q21) * y + 2 * Q22 * y);
	}
	else if (options.initType == 1)
		return speeds.gradSolution(getCoords(index), &speeds.sp);
	else
		return PairDoub(0, 0);
}

/*
	Set cout introduction.
*/
void FastMarch::printIntro() {
	std::cout << "Running fast march with ";
	if (true) std::cout << "8"; else std::cout << "4";
	std::cout << "-point nearest neighbor updates\n\n"
		<< "Speed: \t\ts(x,y) = " << speeds.str << "\n\n"
		<< "Grid: \t \t[" << grid.xl << ", " << grid.xr << "] x [" << grid.yl << ", " << grid.yr << "]\n"
		<< "Spacing: \t(nx,ny) = (" << grid.nx << ", " << grid.ny << ")\n"
		<< "Quadrature: \t";
	/*switch (updater.quadType) {
	case 'e':
		std::cout << "Linear characteristics with endpoint quadrature\n";
		break;
	case 'm': 
		std::cout << "Linear characteristics with midpoint rule quadrature\n";
		break;
	case 's':
		std::cout << "Linear characteristics with Simpson's rule quadrature\n";
		break;
	case 'h':
		std::cout << "Cubic hermite polynomial characteristics with Simpson's rule quadrature\n";
		break;
	}*/
	std::cout << "\n\n";
}

/*
	Set run stats for a Fast March run. 
	This should be done after runMarch is called.
*/
void FastMarch::computeStats(double procTime) {
	stats.run_time = procTime;
	stats.err_sup = computeSupError();
	stats.err_rms = computeRMSError();
	stats.err_grad_sup = computeGradSupError();
	stats.err_grad_rms = computeGradRMSError();
	for (int k = 0; k < accepted.size(); k++) {
		int ind = accepted[k];
		if (invalidIndex(ind)) continue;
		if (options.debugMode) {
			if (debugInfo[ind].lastUpdateType == 1) stats.num_acc_1pt++;
			else if (debugInfo[ind].lastUpdateType == 2) stats.num_acc_2pt++;
		}
		stats.num_acc++;
	}
}

/*
	Print stats to console.
*/
void FastMarch::printStats() {
	std::cout << "Sup error for u:\t\t " << stats.err_sup << std::endl;
	std::cout << "RMS error for u:\t\t " << stats.err_rms << std::endl;
	std::cout << "Sup error for Du:\t\t " << stats.err_grad_sup << std::endl;
	std::cout << "RMS error for Du:\t\t " << stats.err_grad_rms << std::endl;
	std::cout << "Computation time:\t\t " << stats.run_time << std::endl;
	std::cout << "Total 1-point updates:\t " << stats.num_acc_1pt << std::endl;
	std::cout << "Total 2-point updates:\t " << stats.num_acc_2pt << std::endl;
}



void FastMarch::runMarch() {
	initializeSolution();
	int k = 0;
	bool kill = false;
	if (options.verbose) {
		printIntro();
		std::cout << "Completion progress: ";
	}
	clock_t time_s = clock();
	while (!consideredEmpty()) {
		runStep(kill);
		if (kill == true)
			break;
		if (options.verbose) {
			if ((k++) % (grid.gridLength() / 20) == 0)
				std::cout << 5 * k / (grid.gridLength() / 20) << "%...";
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

/*
	Run a single step.
*/
void FastMarch::runStep(bool &kill) {
	int ind = considered.takeMin();	
	if (inBoundary(ind)) {
		kill = true;
		return;	
	}	
	group[ind] = ACCEPTED;
	accepted.push_back(ind);
	updateNeighbors(ind);
}

/*
	Update neighbors of mesh point index. OVERLOADED function.
*/
void FastMarch::updateNeighbors(int index) {
	const std::vector<PairInt>& neighborShifts =  grid.eightPtNeighborShifts;
	int n = neighborShifts.size();
	for (int k = 0; k < n; k++) {
		PairInt indexPairUpdate{ getIndexPair(index)+ neighborShifts[k] };
		if (invalidIndexPair(indexPairUpdate)) continue;
		int index_update = getIndex(indexPairUpdate);
		if (group[index_update] == ACCEPTED)  continue;
		/* Remark:
		We choose to perform the 1-point update first before checking for 2-point updates.
		Note that 2-point updates only update the values of u if there is an interior
			minimizer, and not if the minimizer occurs at the endpoint. Putting the 1-point
			update in front of the 2-point updates handles that.		
		*/	
		onePointUpdate(index_update, index); 
		for (auto m : { -1,1 }) {
			int ind_updater =  getIndex( getIndexPair(index_update) + neighborShifts[(n+n/2+k+m)%n]);
			/* Remark:
			The 1st n eliminates issue with modular arithmetic and negative numbers.
			The n/2 flips the shifts to be measured from index_update, whereas
				originally they were measured in the opposite direction from index.			
			The m (either 1 or -1) ensures that only adjacent vertices can be considered for triangles.
			*/
			if (index_update == 2093053) {
				(getIndexPair(index) - getIndexPair(index_update)).print();
				std::cout << "\t";
				(getIndexPair(ind_updater) - getIndexPair(index_update)).print();
				std::cout << std::endl;
			}
			if (invalidIndex(ind_updater)) continue;
			if (group[ind_updater] == ACCEPTED)
				twoPointUpdate(index_update, index, ind_updater);
		}
		if (group[index_update] == FAR) {
			group[index_update] = CONSIDERED;
			considered.tree.push_back(index_update);
			indexInConsidered[index_update] = considered.tree.size() - 1;
		}
		considered.reassembleFromTarget(indexInConsidered[index_update]);
	}
}
/*
	Update mesh point z1_i from z2_i.
*/
void FastMarch::onePointUpdate(int z1_i, int z2_i) {
	const SolverParams1Pt params{ getCoords(z1_i),getCoords(z2_i),u[z2_i],uGrad[z2_i],speeds };
	PairDoub u1grad(INFINITY,INFINITY);
	Flags flags_new;
	double unew = updater.onePointUpdateValue(params,u1grad, flags_new);
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

/*
	Update mesh point z1_i via triangle update from z2_i and z3_i.
*/
void FastMarch::twoPointUpdate(int z1_i, int z2_i, int z3_i) {
	// Checks if this two point update would result in an update distance 
	// larger than the update distance from the previous successful 2-point update.
	double dist_new = std::min((getCoords(z1_i) - getCoords(z2_i)).norm(),
		(getCoords(z1_i) - getCoords(z3_i)).norm());
	bool tooFarAway = options.twoPtUpdateMinDist
		&& updateType[z1_i] >= 2
		&& (dist_new >= lastUpdateDist[z1_i]);
	if (tooFarAway) 
		return;

	// Run the minimization.
	const SolverParams2Pt params{ getCoords(z1_i),getCoords(z2_i),
		getCoords(z3_i),u[z2_i],u[z3_i], uGrad[z2_i],uGrad[z3_i],speeds };
	bool interiorMinimizer = false;
	PairDoub u1grad(INFINITY, INFINITY);
	Flags flags_new;
	double unew = updater.twoPointUpdateValue(params, u1grad,
		interiorMinimizer, flags_new);

	// Checks if the new value is less than a previous 1-point update value.
	// This is designed to filter out other local minima, not corresponding to the
	// action characteristic.
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

/*
	Write solution to csv file.
*/

void FastMarch::computeLaplacian() {
	for (unsigned x = 0; x < grid.nx ; x++)
		for (unsigned y = 0; y < grid.ny ; y++) {
			int index = y * grid.nx + x;
			if (group[index] != ACCEPTED && group[index] != FRONT) continue;
			int index_right = y * grid.nx + x + 1;
			int index_left = y * grid.nx + x - 1;
			int index_up = (y + 1) * grid.nx + x;
			int index_down = (y - 1) * grid.nx + x;
			double dxx = 0, dyy = 0;
			if (!invalidIndex(index_right) && !invalidIndex(index_left) && group[index_right] >= FRONT && group[index_left] >= FRONT)
				dxx = (uGrad[index_right].x - uGrad[index_left].x) / (2 * grid.hx);
			else if (!invalidIndex(index_right) &&  group[index_right >= FRONT])
				dxx = (uGrad[index_right].x - uGrad[index].x) / (grid.hx);
			else if (!invalidIndex(index_left) && group[index_left] >= FRONT)
				dxx = (uGrad[index].x - uGrad[index_left].x) / grid.hx;
			else
				dxx = INFINITY;

			if (!invalidIndex(index_up) && !invalidIndex(index_down) && group[index_up] >= FRONT && group[index_down] >= FRONT)
				dyy = (uGrad[index_up].y - uGrad[index_down].y) / (2 * grid.hy);
			else if (!invalidIndex(index_up) && group[index_up >= FRONT])
				dyy = (uGrad[index_up].y - uGrad[index].y) / grid.hy;
			else if (!invalidIndex(index_down) && group[index_down >= FRONT])
				dyy = (uGrad[index].y - uGrad[index_down].y) / grid.hy;
			else
				dyy = INFINITY;
			laplaceU[index] = dxx + dyy;
		}
}

void FastMarch::writeToTXT() {
	std::ofstream sol_out("Outputs/solution.csv");
	bool partial = false;
	//sol_out << "x,y,U,Err(U),Err(DU),Last Update Type, DU_x,DU_y" << std::endl;
	if (partial == true) {
		// Compute this on the center 1/4 of the grid.
		unsigned xs = grid.nx / 4, xf = 3 * grid.nx / 4, ys = grid.ny / 4, yf = 3 * grid.ny / 4;
		int totXN = xf - xs, totYN = yf - ys;
		int XNMAX = 500, YNMAX = 500;
		int num = 0, numX = 0, numY = 0, Ycounter = 0;
		// In order to control the size of the output file, we limit the number of points saved to 500 by 500.
		//sol_out << "x,y,u,e,eg,type";
		for (unsigned x = xs; x < xf; x++) {
			Ycounter = 0;
			for (unsigned y = ys; y < yf; y++) {
				int index = y * grid.nx + x;
				if (((x * XNMAX / totXN - (x - 1) * XNMAX / totXN) >= 1) && ((y * YNMAX / totYN - (y - 1) * YNMAX / totYN) >= 1)) {// note integer division is required. 
					sol_out << getXValue(index) << "," << getYValue(index) << "," << u[index] << "," << abs(u[index] - solution(index)) << "," << (uGrad[index] - gradSolution(index)).norm() << ","
						<< debugInfo[index].lastUpdateType << std::endl;
					Ycounter++;
					numY = std::max(numY, Ycounter);
					num++;
				}
			}
		}
		numX = num / numY;
		sol_out << numX << "," << numY << "," << num << std::endl;
		sol_out.close();
	} 
	else {
		for (unsigned x=0; x < grid.nx; x++) 
			for (unsigned y = 0; y < grid.ny; y++) {
				int index = y * grid.nx + x;
				sol_out << getXValue(index) << "," << getYValue(index) << "," << u[index] << "," << abs(u[index] - solution(index)) << "," << (uGrad[index] - gradSolution(index)).norm() << ","
					<< debugInfo[index].lastUpdateType << "," << uGrad[index].x << "," << uGrad[index].y << std::endl;
			}
		sol_out << grid.nx << "," << grid.ny << "," << grid.nx*grid.ny << "," << speeds.sp.switchKey << "," << speeds.sp.a << "," << speeds.sp.b << std::endl;
		sol_out.close();

	}
	std::ofstream sol_out2("Outputs/solutionAccepted.csv");
	//sol_out2 << "x,y,U,Err(U),Err(DU),Last Update Type, DU_x,DU_y" << std::endl;
	for (unsigned x = 0; x < grid.nx; x++)
		for (unsigned y = 0; y < grid.ny; y++) {
			int index = y * grid.nx + x;
			if (group[index] == ACCEPTED) {
				sol_out2 << getXValue(index) << "," << getYValue(index) << "," << u[index] << "," << abs(u[index] - solution(index)) << "," << (uGrad[index] - gradSolution(index)).norm() << ","
					<< debugInfo[index].lastUpdateType << "," << uGrad[index].x << "," << uGrad[index].y << "," << laplaceU[index] << "," 
					<< getCoords( debugInfo[index].ind2).x << "," << getCoords(debugInfo[index].ind2).y << "," 
					<< getCoords(debugInfo[index].ind3).x << "," << getCoords(debugInfo[index].ind3).y << ","
					<< debugInfo[index].lam << "," << debugInfo[index].a0 << "," << debugInfo[index].a1 << "," << debugInfo[index].rank 
					<<  std::endl;
			}
			else
				sol_out2 << getXValue(index) << "," << getYValue(index) << "," << "nan" << "," << abs(u[index] - solution(index)) << "," << (uGrad[index] - gradSolution(index)).norm() << ","
				<< -1<< "," << uGrad[index].x << "," << uGrad[index].y << "," << laplaceU[index] << ","
				<< debugInfo[index].z1.x << "," << debugInfo[index].z1.y << "," << debugInfo[index].z2.x << "," << debugInfo[index].z2.y << ","
				<< 0<< "," << 0 << "," << 0 << std::endl;
		}
	sol_out2 << grid.nx << "," << grid.ny << "," << grid.nx * grid.ny << "," << speeds.sp.switchKey << "," << speeds.sp.a << "," << speeds.sp.b << std::endl;
	sol_out2.close();
}

/*
	Compute stats used in Masha's OLIM files.
*/
void FastMarch::computeMashaStats(std::vector<double>* stats) {
	*stats = std::vector<double>(8, 0);
	double max_err = 0, max_grad_err = 0, max_u = 0;
	double rms_sum = 0, rms_grad_sum = 0, rms_u_sum = 0;
	double max_laplace_err = 0, rms_laplace_sum = 0;
	unsigned num_points = 0, num_laplace = 0;
	for (int x = 0; x < grid.nx; x++) {
		for (int y = 0; y < grid.ny; y++) {
			int index = getIndex(PairInt(x,y));
			if (group[index] == ACCEPTED) {
				double abs_err = abs(u[index] - solution(index));
				double abs_grad_err = (uGrad[index] - gradSolution(index)).norm();
				double abs_laplace_err = abs(laplaceU[index] - laplaceSolution(index));
				max_err = std::max(abs_err, max_err);
				max_grad_err = std::max(abs_grad_err, max_grad_err);

				max_u = std::max(u[index], max_u);
				rms_sum += abs_err * abs_err;
				rms_grad_sum += abs_grad_err * abs_grad_err;
				rms_u_sum += u[index] * u[index];
				if (abs_laplace_err < INFINITY) {
					max_laplace_err = std::max(abs_laplace_err, max_laplace_err);
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

double FastMarch::computeSupError() {
	double max_error = 0;
	for (int x = 0; x < grid.nx; x++) {
		for (int y = 0; y < grid.ny; y++) {
			int index = y * grid.ny + x;
			if ((group[index] == ACCEPTED)
				&& abs(u[index] - solution(index)) > max_error) {
				max_error = abs(u[index] - solution(index));
			}
		}
	}
	return max_error;
}

double FastMarch::computeRMSError() {
	double rms_sum = 0;
	unsigned num_points = 0; 
	for (int ind = 0; ind < grid.gridLength() ; ind++) {
		double u_exact = solution(ind);
		if (group[ind] == ACCEPTED) {
			rms_sum += (u[ind] - u_exact) * (u[ind] - u_exact);
			num_points++;
		}
	}
	return sqrt(rms_sum / num_points);
}
double FastMarch::computeGradSupError() {
	double max_error = 0;
	for (int x = 0; x < grid.nx; x++) {
		for (int y = 0; y < grid.ny; y++) {
			int index = y * grid.ny + x;
			double curr_err = (uGrad[index] - gradSolution(index)).norm();
			if ( (group[index] == ACCEPTED)	&& curr_err > max_error) max_error = curr_err;
		}
	}
	return max_error;
}
double FastMarch::computeGradRMSError() {
	double sumSquares = 0;
	unsigned num_points = 0;
	for (int ind = 0; ind < grid.gridLength(); ind++) {
		if (group[ind] == ACCEPTED) {
			sumSquares += (uGrad[ind] - gradSolution(ind)).normsq();
			num_points++;
		}
	}
	return sqrt(sumSquares / num_points);
}

/*
	Writes point by point information about ACCEPTED mesh points to file, in 
	order of which they were added to accepted.
	Default size limit is 50,000 otherwise file gets too big. 
	Can change this within. 
*/
void FastMarch::writeDebug() {
	std::ofstream sol_out("Outputs/debug.csv");
	sol_out << std::setprecision(12); 

	sol_out
		<< "rank,i1,e1,e2,e3,ge1,x1,y1,u1,gu1x,gu1y,i2,x2,y2,u2,gu2x,gu2y,i3,x3,y3,u3,gu3x,gu3y,type,lambda,a0,a1,dist" << std::endl;
	int rank = 0;
	int sizeLimit = grid.gridLength();
	sizeLimit = 50000;
	for (unsigned k = 0; k < accepted.size(); k++) {


		int ind1 = accepted[k];
		if (initialPoint(ind1)) continue;
		if (rank > sizeLimit) break;
		rank++;
		UpdateInfo inf = debugInfo[ind1];
		int ind2 = inf.ind2, ind3 = inf.ind3;
		sol_out
			<< rank << ","
			<< ind1 << ","
			<< u[ind1] - solution(ind1) << ","
			<< u[ind2] - solution(ind2) << ","
			<< u[ind3] - solution(ind3) << ","
			<< (uGrad[ind1] -gradSolution(ind1) ).norm() << ","
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
			<< (getCoords(ind1) - getCoords(ind2)).norm() / grid.hx
			<< std::endl;
	}
	sol_out << grid.xl << "," << grid.xr << "," << grid.nx << "," << grid.yl << "," << grid.yr << "," << grid.ny << "," << speeds.sp.a << "," << speeds.sp.b << "," << speeds.sp.switchKey <<
		std::endl;

	sol_out.close();


}


void FastMarch::printDebugConsole() {
	std::cout << std::setprecision(4);

	int rank = 0;
	int sizeLimit = grid.gridLength();
	sizeLimit = 100;
	for (unsigned k = 0; k < accepted.size(); k++) {


		int ind1 = accepted[k];
		if (rank > sizeLimit) break;
		rank++;
		UpdateInfo inf = debugInfo[ind1];
		int ind2 = inf.ind2, ind3 = inf.ind3;
		std::cout << "Rank:\t" << rank << "\tIndex:\t" << ind1 << "\tU:\t" << u[ind1] << "\tType:\t" << inf.lastUpdateType << "\ti2:\t" << ind2 << "\ti3:\t" << ind3 <<std::endl;

	}


}

/*
	Run num'th iteration of the failsafe (i.e. ell^1 shell of radius num)
*/
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


		if (!grid.invalidIndex(i2) && !grid.invalidIndex(i3) && group[i2] == ACCEPTED && group[i3] == ACCEPTED)
			twoPointUpdate(ind, i2, i3);
		xprev = x0, yprev = y0;
		x0 += xinc; y0 += yinc;
		if (abs(x0) == num)
			xinc *= -1;
		else if (abs(y0) == num)
			yinc *= -1;
	}
}

/*
	More thorough fail-safe...don't use this.
*/
void FastMarch::fixFailures_soupedUp(int ind, int num) {
	PairInt ipair = getIndexPair(ind);
	int x0 = -num, y0 = 0;
	int xprev = -num + 1, yprev = 1;
	int xinc = 1, yinc = -1;

	for (int k = 0; k < 4 * num; k++) {
		PairInt p2{ x0,y0 };
		PairInt p3{ xprev,yprev };
		int i2 = getIndex(ipair + p2);


		for (const PairInt& shift : grid.eightPtNeighborShifts) {
			int i3 = getIndex(getIndexPair(i2) + shift);
			if (!grid.invalidIndex(i2) && !grid.invalidIndex(i3) && group[i2] == ACCEPTED && group[i3] == ACCEPTED)
				twoPointUpdate(ind, i2, i3);
		}
		xprev = x0, yprev = y0;
		x0 += xinc; y0 += yinc;
		if (abs(x0) == num)
			xinc *= -1;
		else if (abs(y0) == num)
			yinc *= -1;
	}
}

void FastMarch::writeMAP(const PairDoub& z) {
	double tol = 1e-4;
	PairDoub zcurr(z), zprev(z);
	PairDoub gu(INFINITY, INFINITY);
	while (zcurr.norm() > tol) {
		zprev = zcurr;
		double d = std::min(zcurr.norm() / 2.0, this->grid.hx);
		
		zcurr = zcurr;

	}

}



double FastMarch::shootMaps(PairDoub x0) {
	PairDoub origin(0, 0);
	double step = grid.hx/10; // Set to h/10.

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
			xcur = xcur - step * (gucur + drift(xcur)) / (gucur + drift(xcur)).norm();
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
		std::size_t pos = line.find(",");      // position of the end of the name of each one in the respective string
		PairDoub point(std::stod(line.substr(0, pos)), std::stod(line.substr(pos + 1, line.size()))); // convert string age to a double
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
		
		PairDoub p2_cur = points[jmin], p2_next = points[jmin + 1], p2_prev = points[jmin - 1];
		PairDoub p1 = path[i];

		if (p1.norm()  < 0.1) {
			errs.push_back(0);
			continue;
		}


		double dotprod = dot(p1 - p2_prev, p2_cur - p2_prev);
		double d12 = (p1 - p2_prev).normsq() - dotprod * dotprod / (p2_cur - p2_prev).normsq();
		d12 = sqrt(d12);
		double dotprod2 = dot(p1 - p2_cur, p2_next - p2_cur);
		double d23 = (p1 - p2_cur).normsq() - dotprod2 * dotprod2 / (p2_next - p2_cur).normsq();
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

PairDoub FastMarch::interpolate(PairDoub xcur) {
	PairInt bl(floor((xcur.x - grid.xl) / grid.hx), floor((xcur.y - grid.yl) / grid.hy));
	PairInt br(bl.xi + 1, bl.yi), tl(bl.xi, bl.yi + 1), tr(bl.xi + 1, bl.yi + 1);
	int bli = getIndex(bl), bri = getIndex(br), tli = getIndex(tl), tri = getIndex(tr);
	PairDoub gu_bl(uGrad[bli]), gu_br(uGrad[bri]), gu_tl(uGrad[tli]), gu_tr(uGrad[tri]);
	PairDoub gu(0, 0);
	double gux = 0, guy = 0;
	MatrixDoub xmat(PairDoub(gu_bl.x, gu_tl.x), PairDoub(gu_br.x, gu_tr.x));
	MatrixDoub ymat(PairDoub(gu_bl.y, gu_tl.y), PairDoub(gu_br.y, gu_tr.y));
	double x1 = getCoords(bli).x, x2 = getCoords(bri).x;
	double y1 = getCoords(bli).y, y2 = getCoords(tli).y;
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

double FastMarch::interpLaplacian(const PairDoub &xcur) {
	PairInt bl(floor((xcur.x - grid.xl) / grid.hx), floor((xcur.y - grid.yl) / grid.hy));
	PairInt br(bl.xi + 1, bl.yi), tl(bl.xi, bl.yi + 1), tr(bl.xi + 1, bl.yi + 1);
	int bli = getIndex(bl), bri = getIndex(br), tli = getIndex(tl), tri = getIndex(tr);
	double gu_bl(laplaceU[bli]), gu_br(laplaceU[bri]), gu_tl(laplaceU[tli]), gu_tr(laplaceU[tri]);
	double gu(0);
	MatrixDoub xmat(PairDoub(gu_bl, gu_tl), PairDoub(gu_br, gu_tr));
	double x1 = getCoords(bli).x, x2 = getCoords(bri).x;
	double y1 = getCoords(bli).y, y2 = getCoords(tli).y;
	PairDoub vec1(x2 - xcur.x, xcur.x - x1);
	PairDoub vec2(y2 - xcur.y, xcur.y - y1);
	double fac = 1.0 / ((x2 - x1) * (y2 - y1));
	gu = fac * dot(vec1, rightMultiply(xmat, vec2));
	return gu;

}