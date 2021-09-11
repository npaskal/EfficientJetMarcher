#define _USE_MATH_DEFINES 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "fastmarchASR.h"

void FastMarchASR::initializeSolution() {
	std::vector<std::pair<double, int>> pairs;
	switch (options.initType) {
	case 0: case 2: { // Default initialization
		for (int ind = 0; ind < grid.gridLength(); ind++) {
			if (initialPoint(ind)) {
				u[ind] = initialValue(ind);
				uGrad[ind] = initialGradU(ind);
				group[ind] = ACCEPTED;
				pairs.push_back(std::make_pair(u[ind], ind));
			}
		}
		break;
	}
	case 1: { // Runs solver on finer grid around equilibria.
		if (options.verbose)
			std::cout << "Initializing:\n";
		int numSteps = 64; // must be multiple of 2 as currently configured.
		double tol_deg = 15;
		double tol = M_PI / 180.0 * tol_deg;
		double usqrt = 1.0 / tol * grid.hx * sqrt(1 + speeds.sp.b * speeds.sp.b) / sqrt(.5);
		usqrt = std::min(usqrt, .2);
		int xi_left = (1.0 - usqrt) / grid.hx - 1;
		int xi_right = (grid.nx - 1 - xi_left);
		PairDoub coord_low = getCoords(getIndex(PairInt(xi_left, xi_left)));
		PairDoub coord_high = getCoords(getIndex(PairInt(xi_right, xi_right)));
		double time_frac_allowed = 0.25;
		int per_point = (double)grid.nx / (double)(xi_right - xi_left) * sqrt(time_frac_allowed);
		int n_new = per_point * (xi_right - xi_left) + 1;

		//const GridInfo gridInit((grid.nx + 1) / 2, coord_low.x, 
		//	coord_high.x, (grid.ny + 1) / 2, coord_low.y, coord_high.y);
		n_new = 129;
		coord_low.x = -.002; coord_high.x = .002; coord_low.y = -.002; coord_high.y = .002;
		const GridInfo gridInit(n_new, coord_low.x,
			coord_high.x, n_new, coord_low.y, coord_high.y);
		RunOptions optsInit(options);
		optsInit.verbose = false;
		optsInit.initType = 0;
		FastMarchASR marchInit(gridInit, updater, optsInit, speeds, stencilData);
		marchInit.runMarch();

		for (int xi = xi_left; xi <= xi_right; xi++)
			for (int yi = xi_left; yi <= xi_right; yi++) {
				int ind = getIndex(PairInt(xi, yi));
				int xi_new = (xi - xi_left) * std::round(grid.hx / gridInit.hx);
				int yi_new = (yi - xi_left) * std::round(grid.hy / gridInit.hy);
				int ind_new = marchInit.getIndex(PairInt(xi_new, yi_new));
				if (marchInit.group[ind_new] == ACCEPTED) {
					u[ind] = marchInit.u[ind_new];
					uGrad[ind] = marchInit.uGrad[ind_new];
					group[ind] = ACCEPTED;
					pairs.push_back(std::make_pair(u[ind], ind));
				}
			}
		for (int xi = xi_left; xi <= xi_right; xi++)
			for (int yi = xi_left; yi <= xi_right; yi++) {
				int ind = getIndex(PairInt(xi, yi));
				int xi_new = (xi - xi_left) * std::round(grid.hx / gridInit.hx);
				int yi_new = (yi - xi_left) * std::round(grid.hy / gridInit.hy);
				int ind_new = marchInit.getIndex(PairInt(xi_new, yi_new));
				if (marchInit.group[ind_new] == ACCEPTED) {
					u[ind] = marchInit.u[ind_new];
					uGrad[ind] = marchInit.uGrad[ind_new];
					group[ind] = ACCEPTED;
					pairs.push_back(std::make_pair(u[ind], ind));
				}
			}
		break;
	}
	}

	// To enforce causality, we sort the initialized points by u, which may be
	// time-consuming if the initialization region is large.
	std::sort(pairs.begin(), pairs.end());
	for (auto& p : pairs) {
		accepted.push_back(p.second);
		updateNeighbors(p.second);
	}
}

void FastMarchASR::createStencils_Mirebeau() {
	// Loop through the grid of angular directions of the drift field b.
	stencil.reserve(stencilData.num_theta_cells);
	//std::vector<PairInt> M0{ PairInt(1,0),PairInt(0,-1),PairInt(-1,0),PairInt(0,1) }, L0{ PairInt(1,0) };
	std::vector<PairInt> M0{ PairInt(1,0),PairInt(1,-1),PairInt(0,-1),PairInt(-1,-1),PairInt(-1,0),PairInt(-1,1),PairInt(0,1),PairInt(1,1) }, L0{ PairInt(1,0) };
	for (int i = 0; i < stencilData.num_theta_cells; i++) { 
		double theta = 2 * M_PI * double(i) / double(stencilData.num_theta_cells);
		// Theta is the NEGATIVE angle of drift b. The stencil will be refined more finely in this direction.
		// NOTE the stencil is the set of points that affect you. This should be centered on the points
		// coming inwards, since they expend no energy following the drift inwards, hence you look outwards.
		PairDoub bhat{ cos(theta),sin(theta) }; 
		std::vector<PairInt> L{ L0 }, M{ M0 };
		while (!M.empty()) {
			PairInt u = M.back(), v = L.back();
			PairDoub uhat{ (double)u.xi/u.norm(),(double)u.yi/u.norm() }; 
			PairDoub vhat{ (double)v.xi/v.norm(),(double)v.yi/v.norm() }; 
			double uhvh = dot(uhat, vhat);
			double temp2 = stencilData.alpha * std::max(dot(uhat, bhat), dot(vhat, bhat));
			if ((dot(uhat,vhat) < stencilData.alpha * std::max(dot(uhat,bhat),dot(vhat,bhat))) 
				&& (abs(u.xi + v.xi) + abs(u.yi + v.yi) <= stencilData.stencil_cutoff)) {
				M.push_back(u+v);
			}
			else {
				L.push_back(M.back());
				M.pop_back();
			}
		}
		L.pop_back();
		stencil.push_back(L);
	}
	
	// This next loop precomputes the direction of b for all grid points
	// and saves the z2_i of the appropriate stencil.
	// h = the width of a theta bin = 2*PI/n.
	// Bins are [ (2n-1) h/2,h/2),[h/2,3h/2],...,[(2n-3) h/2,(2n-1)h/2], where
	// convention will be the interval is inclusive on the left.
	double h = 2 * M_PI / stencilData.num_theta_cells;
	for (int index = 0; index < grid.gridLength(); index++) {
		PairDoub b = drift(index);
		double theta = atan2(-b.y,- b.x);
		if (theta < 0)
			theta += 2 * M_PI;
		stencilIndex[index] = int((theta + h / 2) / h) % stencilData.num_theta_cells;
	}
}


void FastMarchASR::createStencils_bubble() {
	// first initialize all stencils to the 4-point stencil
	stencil = std::vector<std::vector<PairInt>>(stencilData.num_theta_cells,
		std::vector<PairInt>{PairInt(1, 0), PairInt(0, 1), PairInt(-1, 0), PairInt(0, -1)});
	double delta = .1;
	int stencSize = stencilData.theta_res.size();
	std::vector<int> indices_reordered(stencSize,0);
	int midInd = stencSize / 2; // supposed to be truncated
	for (int i = 0; i < stencilData.num_theta_cells; i++) {
		stencil.push_back(std::vector<PairInt>());
		double theta_bneg = 2 * M_PI * double(i) / double(stencilData.num_theta_cells);
		double prev_dist = 0;
		for (int k = 0; k < stencSize; k++) {
			int ind = INT_MAX;
			if (k <= midInd) ind = midInd - k;
			else ind = k;
			if (k == midInd + 1) prev_dist = 0; // reset for the 2nd half of points
			double theta = std::fmod(stencilData.theta_res[ind] + theta_bneg,2*M_PI);
			int quad = int(theta / (M_PI / 2));
			PairInt init(0, 0), fin(0, 0);
			switch (quad) {
			case 0:
				init = PairInt(0, 1); fin = PairInt(1, 0);
				break;
			case 1:
				init = PairInt(-1, 0); fin = PairInt(0, 1);
				break;
			case 2:
				init = PairInt(0, -1); fin = PairInt(-1, 0);
				break;
			case 3:
				init = PairInt(1, 0); fin = PairInt(0, -1);
				break;
			}
			int rad = 1;
			PairInt increment = fin - init;
			bool found = false;

			while (found == false) {
				PairInt init_new(rad * init), fin_new(rad * fin);
				PairInt curr(init_new);
				for (int j = 0; j < rad + 1; j++) {
					double dist = abs(dot(curr, PairDoub(cos(theta + M_PI / 2), sin(theta + M_PI / 2))));
					if (dist < delta && (curr.norm() > std::min(prev_dist,(double)stencilData.stencil_cutoff))) {
						stencil[i].push_back(curr);
						prev_dist = curr.norm();
						found = true;
						break;
					}
					else if (dist < delta) {
						stencil[i].push_back(curr);
					}
					curr = curr+ increment;
				}
				rad++;
			}

		}
		for (auto& pair : grid.eightPtNeighborShifts)
			stencil[i].push_back(pair);
		//std::cout << "\n Stencil Size is " << stencil[i].size();
	}

	double h = 2 * M_PI / stencilData.num_theta_cells;
	for (int index = 0; index < grid.gridLength(); index++) {
		PairDoub b = drift(index);
		double theta = atan2(-b.y, -b.x);
		if (theta < 0)
			theta += 2 * M_PI;
		stencilIndex[index] = int((theta + h / 2) / h) % stencilData.num_theta_cells;
	}
}

void FastMarchASR::createStencils_bubble_new() {
	double L = 12.0;
	double ell = 3.0;
	int N = 5;
	std::vector<PairDoub> coords;
	
	for (int i = 0; i < 2 * N - 1; i++) {
		double legLength = 0;
		double theta = 0;
		if (i < N - 1) {
			theta = M_PI / 2.0 * pow(0.5, i + 1);
			legLength = (i + 1) * (L / static_cast<double>( N));
		}
		else if (i == N - 1) {
			theta = 0.0;
			legLength = (i + 1) * (L / static_cast<double>(N));
		} 
		else {
			theta = -M_PI / 2.0 * pow(0.5, 2 * N - (i + 1));
			legLength = (2 * N - (i + 1)) * (L / static_cast<double>(N));
		}
		int numSamples = ceil( legLength / ell);
		for (int j = 1; j <= numSamples; j++) {
			double curLen = static_cast<double>(j)* legLength / static_cast<double>(numSamples);
			coords.push_back(PairDoub(curLen * cos(theta), curLen * sin(theta)));
		}
	}
	for (int i = 0; i < stencilData.num_theta_cells; i++) {
		double theta_bneg = 2 * M_PI * double(i) / double(stencilData.num_theta_cells);
		std::vector<PairDoub> coords_rot;
		std::vector<PairInt> coords_rounded;
		for (int j = 0; j < coords.size(); j++) {
			coords_rot.push_back(coords[j].rotate(theta_bneg));
			coords_rounded.push_back(PairInt(round(coords_rot[j].x), round(coords_rot[j].y)));

		}
		for (auto& pair : grid.eightPtNeighborShifts)
			coords_rounded.push_back(pair);
		std::sort(coords_rounded.begin(), coords_rounded.end());
		coords_rounded.erase(std::unique(coords_rounded.begin(), coords_rounded.end()), coords_rounded.end());
		stencil.push_back(coords_rounded);
	}

	double h = 2 * M_PI / stencilData.num_theta_cells;
	for (int index = 0; index < grid.gridLength(); index++) {
		PairDoub b_neg = -1.0*drift(index);
		double theta = b_neg.angle();
		stencilIndex[index] = int((theta + h / 2) / h) % stencilData.num_theta_cells;
	}
	std::cout << std::endl;
}


void FastMarchASR::updateNeighbors(int z2_i) {
	int theta_i = stencilIndex[z2_i]; // theta bin z2_i
	PairInt z2_ip = getIndexPair(z2_i);

	int ind_track = INT_MAX, ind_track2 = INT_MAX;
	for (int i = 0; i < static_cast<int64_t>(stencil.at(theta_i).size());i++) {
		PairInt z1_ip{ getIndexPair(z2_i) - stencil[theta_i][i]};
		int z1_i = getIndex(z1_ip); 
		if (invalidIndexPair(z1_ip) || invalidIndex(z1_i) 
			|| group[z1_i] == ACCEPTED)	continue;
		const std::vector<PairInt> *shifts{ nullptr };
		std::vector<PairInt> shifts_noBox;
		if (stencilData.boxUpdate)
			shifts = &grid.eightPtNeighborShifts;
		else {
			shifts_noBox.push_back(z1_ip - z2_ip 
				+ stencil[theta_i][(i + (int)1) % stencil[theta_i].size()]);
			shifts_noBox.push_back(z1_ip - z2_ip + stencil[theta_i]
				[(i + stencil[theta_i].size() - 1) % stencil[theta_i].size()]);		
			shifts = &shifts_noBox;
		}
		onePointUpdate(z1_i, z2_i); // Precautionary 1-pt update
		for (const PairInt& shift : *shifts) {
			int z3_i = getIndex(z2_ip + shift);
			int ind_track = INT_MAX;
			if (z1_i == ind_track && z2_i == ind_track2) {
				(getIndexPair(z2_i) - getIndexPair(z1_i)).print();
				std::cout << "\t";
				(getIndexPair(z3_i) - getIndexPair(z1_i)).print();
				std::cout << std::endl;
			}


			if (invalidIndex(z3_i))	
				continue;
			else if (group[z3_i] == ACCEPTED ) {
				if (options.z3one) onePointUpdate(z1_i, z3_i); // Precautionary 1-pt update
				twoPointUpdate(z1_i, z2_i, z3_i);
			}
		}
		if (group[z1_i] == FAR) {
			group[z1_i] = CONSIDERED;
			considered.tree.push_back(z1_i);
			indexInConsidered[z1_i] = considered.tree.size() - 1;
		}
		considered.reassembleFromTarget(indexInConsidered[z1_i]);
	}
}

void FastMarchASR::runStep(bool& kill) {
	int ind = considered.takeMin();
	if (inBoundary(ind)) {
		kill = true;
		// return;
	}
	// In the event that 'ind' was last updated with a 1-point update, we now
	// cycle outwards through shells in L^1, until we get to a successful 2-point update.
	int num = 1;
	if (updateType[ind] == 1) failSafe++;
	while ( options.failSafe && updateType[ind] == 1 && num < 20) {
		fixFailures(ind,num);

		num += 1;
	}
	num = 1;
	while (options.failSafe && updateType[ind] == 1 && num < 20) {
		//fixFailures_soupedUp(ind, num);

		num += 1;
	}
	group[ind] = ACCEPTED;

	accepted.push_back(ind);
	debugInfo[ind].rank = (int)accepted.size();
	if (kill) {
		return;
	}
	updateNeighbors(ind);
}


/*
Cycles through neighboring pairs (z2,z3) of the perimeter of the ball B(z1, num*h)
in L^1. It then attempts to perform a 2-point update of z1 from z2 and z3.
*/


void FastMarchASR::writeStencil(int index) {
	std::ofstream sol_out("Outputs/stencil.csv");
	sol_out << "index,xi,yi,x,y\n";
	sol_out << index << "," << getXIndex(index) << "," << getYIndex(index) << ","
		<< getXValue(index) << "," << getYValue(index) << std::endl;
	int si = stencilIndex[index];
	PairInt coords = getIndexPair(index);
	for (const PairInt& shift : stencil[si]) {
		PairInt neighbor = coords + shift;
		if (invalidIndexPair(neighbor)) continue;
		int neighborIndex = getIndex(neighbor);
		sol_out
			<< neighborIndex << ","
			<< neighbor.xi << ","
			<< neighbor.yi << ","
			<< getCoords(neighborIndex).x << ","
			<< getCoords(neighborIndex).y
			<< std::endl;
	}
	sol_out.close();
}

void FastMarchASR::writeToTXT() {
	int index = 7;
	FastMarch::writeToTXT();
	writeStencil(index);
}

