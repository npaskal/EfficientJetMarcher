#define _USE_MATH_DEFINES 
#include <algorithm>
#include <iostream>
#include "fastmarchOLIM.h"

int round2 = 0;
void FastMarchOLIM::runStep(bool& kill) {
	int ind = considered.takeMin();
	PairInt indPair = getIndexPair(ind);
	int num = 1;
	while (options.failSafe && updateType[ind] == 1 && num < 20) {
		//fixFailures(ind, num);
		num += 1;
	}
	group[ind] = FRONT;
	accepted.push_back(ind);
	// Kill process when a boundary point is accepted.
	if (inBoundary(ind)
		|| indPair.xi <=1 || indPair.xi >= grid.nx - 2 || indPair.yi <=1 || indPair.yi >= grid.ny - 2    ) {
		kill = true;
		return;
	}
	// Change any FRONT points to ACCEPTED if 'ind' is their 8th accepted neighbor.
	for ( PairInt iNeigh : grid.eightPtNeighborShifts) {
		int indNeigh = getIndex(getIndexPair(ind) + iNeigh);
		if ( !invalidIndex(indNeigh)  && (++acceptedNeighborCount[indNeigh] == 8) && (group[indNeigh] == FRONT))
			group[indNeigh] = ACCEPTED;
	}
	//if (debugInfo[ind].lastUpdateType == 1) num1Pt++;
	//else if (debugInfo[ind].lastUpdateType >= 2) num2Pt++;
	updateNeighbors(ind);
}


void FastMarchOLIM::updateNeighbors(int index) {
	/*
	Algorithm:
	1. Separate 8 neighbors of 'index' into FRONT points and FAR points. Disregard any neighbors that are in CONSIDERED or ACCEPTED.
	2. Loop over every grid point 'index_update' within the K radius of 'index'.
		a. Disregard 'index_update' unless it is in CONSIDERED.
		b. Perform a 2-point update of 'index_update' from 'index' and each ind_updater of 'index' that is also in FRONT.			
		c. Reassemble the CONSIDERED list to account for changes.
	3. Loop over every ind_updater 'index_update' of 'index' that is in FAR.
		a. For each 'index_update', loop over every grid point 'ind_cand' within the K radius of 'index_update'.
			i. Disregard 'ind_cand' unless it is in FRONT.
			ii. Perform a 1-point update of 'index_update' from 'ind_cand'.
			iii. Let 'index_min' be the the corresponding 'ind_cand' that led to the lowest update value of 'index_update'.
		b. Loop over every ind_updater 'ind_updater' of 'index_min'.
			i. Disregard 'ind_updater' unless it is in FRONT.																		
			ii. Perform a 2-point update of 'index_update' from 'index_min' and 'ind_updater'.
		c. Add 'index_update' to the CONSIDERED list, and reorder the CONSIDERED list.

	*/

	// Step 1.
	int xi = getXIndex(index), yi = getYIndex(index);
	std::vector<int> frontNeighbors, farNeighbors;
	for (  const PairInt &shift : grid.eightPtNeighborShifts) {
		int ind_new = getIndex(getIndexPair(index) + shift);
		if (group[ind_new] == FRONT)
			frontNeighbors.push_back(ind_new);
		if (group[ind_new] == FAR)
			farNeighbors.push_back(ind_new);
	}

	// Step 2.
	int K2 = update_radius * update_radius;
	for (int xi_new = std::max(xi - update_radius, 0); xi_new <= std::min(xi +update_radius, grid.nx - 1); xi_new++) {
		int yi_diff_max = int(sqrt(abs(K2 - (xi_new - xi) * (xi_new - xi))));
		for (int yi_new = std::max(yi - yi_diff_max, 0); yi_new <= std::min(yi + yi_diff_max, grid.ny - 1); yi_new++) {
			int index_update = yi_new * grid.ny + xi_new;
			double utemp = u[index];
			if (group[index_update] == CONSIDERED) {
				// first do a 2-point update of "index_update" from "ind" and each of the accepted front neighbors of "ind", stored in "frontNeighbors"
				if (index == 70 && index_update == 38) {
					std::cout << "stop" << std::endl;
				}
				onePointUpdate(index_update, index);
				for (const int& ind_updater : frontNeighbors) 
					if (!invalidIndex(ind_updater))	twoPointUpdate(index_update, index, ind_updater);
				considered.reassembleFromTarget(indexInConsidered[index_update]);
			}
		}
	}

	// Step 3.
	PairDoub emp(0, 0);
	for (const int& index_update : farNeighbors) {
		// first identify the accepted Front point in a K ball that minimizes the 1-point update, this is onePointMinimizer
		int onePointMinimizer = INT_MAX;
		u[index_update] = INFINITY;

		double umin = INFINITY, umin2 = INFINITY;
		int a1xi = getXIndex(index_update), a1yi = getYIndex(index_update);
		for (int xi_cand = std::max(a1xi - update_radius, 0); xi_cand <= std::min(a1xi + update_radius, grid.nx - 1); xi_cand++) {
			int yi_diff_max = int(sqrt(abs(K2 - (xi_cand - a1xi) * (xi_cand - a1xi))));
			for (int yi_cand = std::max(a1yi - yi_diff_max, 0); yi_cand <= std::min(a1yi + yi_diff_max, grid.ny - 1); yi_cand++) {
				int ind_cand = yi_cand * grid.ny + xi_cand;
				double x_cand = getXValue(ind_cand), y_cand = getYValue(ind_cand);
				double utemp = INFINITY;
				Flags flags;
				if (group[ind_cand] == FRONT ){//|| group[ind_cand] == ACCEPTED) {
					utemp = updater.onePointUpdateValue(SolverParams1Pt{ getCoords(index_update),getCoords(ind_cand), u[ind_cand], uGrad[ind_cand], speeds },emp, flags);
					if (utemp < umin) {
						umin = utemp;
						onePointMinimizer = ind_cand;
						if ((getIndexPair(ind_cand) - getIndexPair(index_update)).norm() < 2) {
							umin2 = utemp;
						}
					}
				}
			}
		}
		uOne[index_update] = umin2;
		// next do a 2-point update over all accepted front points that are neighbors of onePointMinimizer
		onePointUpdate(index_update, onePointMinimizer);
		for (const PairInt& shift : grid.eightPtNeighborShifts) {
			int ind_updater = getIndex(getIndexPair(onePointMinimizer) + shift);  
			if (ind_updater < 0 || ind_updater >= grid.nx * grid.ny)
				continue;
			// if (group[ind_updater] == FRONT || group[ind_updater] == ACCEPTED) {		// Masha just considers front neighbors, not sure why
			if (group[ind_updater] == FRONT ) {
				twoPointUpdate(index_update, onePointMinimizer, ind_updater);
			}
		}

		group[index_update] = CONSIDERED;
		considered.tree.push_back(index_update);
		indexInConsidered[index_update] = considered.tree.size() - 1;
		considered.reassembleFromTarget(indexInConsidered[index_update]);
	}
}


void FastMarchOLIM::initializeSolution() {
	std::vector<std::pair<double, int>> pairs;
	for (int ind = 0; ind < grid.gridLength(); ind++) {
		if (initialPoint(ind)) {
			u[ind] = initialValue(ind);
			uGrad[ind] = initialGradU(ind);
			group[ind] = FRONT;
			pairs.push_back(std::make_pair(u[ind], ind));
		}
	}
	std::sort(pairs.begin(), pairs.end());
	for (auto& p : pairs) {
		accepted.push_back(p.second);
		updateNeighbors(p.second);
		int ind = p.second;
		for (PairInt iNeigh : grid.eightPtNeighborShifts) {
			int indNeigh = getIndex(getIndexPair(ind) + iNeigh);
			if (!invalidIndex(indNeigh) && (++acceptedNeighborCount[indNeigh] == 8) && (group[indNeigh] == FRONT))
				group[indNeigh] = ACCEPTED;
		}
	}

}


bool FastMarchOLIM::initialPoint(int index) {
	//double ellipse_umax = .00005;
	double ellipse_umax = .01;
	double boxWidth = 1;
	double xbound = INFINITY;
	switch (options.initType) {
	case 1:
		//return  (updater.speeds.solution(getCoords(index), &updater.speeds.sp) < ellipse_umax);
		return (initialValue(index) < ellipse_umax && (abs(getXValue(index)) < xbound && abs(getYValue(index)) < xbound)
			|| abs(getXValue(index)) < (boxWidth + .5) * grid.hx && abs(getYValue(index)) < (boxWidth + .5) * grid.hy); // / (grid.nx / 17.0));
	default:
	case 0:
		return (abs(getXValue(index)) < (boxWidth + .5) * grid.hx && abs(getYValue(index)) < (boxWidth + .5) * grid.hy);
			//&& getXValue(index) >=-.5*grid.hx && getYValue(index) >= -.5*grid.hy );
	}
}


