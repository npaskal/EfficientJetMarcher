/*
	\file:		fastMarchEJM.cpp
	\brief:		This file defines the fastMarchEJM class.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#define _USE_MATH_DEFINES 
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "fastmarchEJM.h"
#include "driftFunctions.h"

// the only EJM constructor
FastMarchEJM::FastMarchEJM(
	const MeshInfo& mesh,
	UpdaterBase& upIn,
	const RunOptions& opIn,
	const SpeedInfo& spIn,
	const StencilData& stencilDataIn
) :
	FastMarch(mesh, upIn, opIn, spIn),
	stencilData(stencilDataIn),
	stencilIndex(std::vector<int>(mesh.gridLength(), INT_MAX))
{
	createStencils();
}

// initialize the mesh points near the attractor
void FastMarchEJM::initializeSolution() {
	std::vector<std::pair<double, int>> pairs;
	switch (options.initType) {
	case 0: case 2: { // Default initialization
		for (int ind = 0; ind < mesh.gridLength(); ind++) {
			if (isInitialMeshPoint(ind)) {
				u[ind] = computeInitialValue(ind);
				uGrad[ind] = computeInitialGradValue(ind);
				group[ind] = ACCEPTED;
				pairs.push_back(std::make_pair(u[ind], ind));
			}
		}
		break;
	}
	// TODO: Clean-up case 1, i.e. running FM on submesh near attractor
	case 1: { // Runs solver on finer mesh around equilibria.
		if (options.verbose)
			std::cout << "Initializing:\n";
		int numSteps = 64; // must be multiple of 2 as currently configured.
		double tol_deg = 15;
		double tol = M_PI / 180.0 * tol_deg;
		double usqrt = 1.0 / tol * mesh.hx 
			* sqrt(1 + speeds.sp.b * speeds.sp.b) / sqrt(.5);
		usqrt = std::min(usqrt, .2);
		int xi_left = (1.0 - usqrt) / mesh.hx - 1;
		int xi_right = (mesh.nx - 1 - xi_left);
		PairDoub coord_low = getCoords(getIndex(PairInt(xi_left, xi_left)));
		PairDoub coord_high = getCoords(getIndex(PairInt(xi_right, xi_right)));
		double time_frac_allowed = 0.25;
		int per_point = (double)mesh.nx 
			/ (double)(xi_right - xi_left) * sqrt(time_frac_allowed);
		int n_new = per_point * (xi_right - xi_left) + 1;

				coord_low.x = -.002; 
		coord_high.x = .002; 
		coord_low.y = -.002; 
		coord_high.y = .002;

		// construct subgrid around attracotr
		// above mesh is for figuring out size of subgrid to use
		const MeshInfo gridInit(
			n_new, 
			coord_low.x,
			coord_high.x, 
			n_new, 
			coord_low.y, 
			coord_high.y
		);
		RunOptions optsInit(options);
		optsInit.verbose = false;
		optsInit.initType = 0;

		// construct EJM object for sub-grid
		FastMarchEJM marchInit(
			gridInit, 
			updater, 
			optsInit, 
			speeds, 
			stencilData);
		marchInit.runMarch();

		for (int xi = xi_left; xi <= xi_right; xi++) {
			for (int yi = xi_left; yi <= xi_right; yi++) {
				int ind = getIndex(PairInt(xi, yi));
				int xi_new = (xi - xi_left) 
					* std::round(mesh.hx / gridInit.hx);
				int yi_new = (yi - xi_left) 
					* std::round(mesh.hy / gridInit.hy);
				int ind_new = marchInit.getIndex(PairInt(xi_new, yi_new));
				if (marchInit.group[ind_new] == ACCEPTED) {
					u[ind] = marchInit.u[ind_new];
					uGrad[ind] = marchInit.uGrad[ind_new];
					group[ind] = ACCEPTED;
					pairs.push_back(std::make_pair(u[ind], ind));
				}
			}
		}
		break;
	}
	}

	// To enforce causality, we sort the initialized points by u, which may be
	// time-consuming if the initialization region is large
	std::sort(pairs.begin(), pairs.end());
	for (auto& p : pairs) {
		accepted.push_back(p.second);
		updateNeighbors(p.second);
	}
}

// construct the collection of stencils
void FastMarchEJM::createStencils() {
	switch (stencilData.stencilType) {

		// This is stencil using (modified) algorithm in Mirebeau 
		// paper "Efficient Fast Marching with Finsler metrics" 
	case ACUTE: {

		// Loop through the mesh of angular directions of the drift field b.
		stencilCollection.reserve(stencilData.numThetaCells);
		// use 8-pt neighbor as base stencil (can also try 4)
		std::vector<PairInt> M0{
			PairInt(1,0),
			PairInt(1,-1),
			PairInt(0,-1),
			PairInt(-1,-1),
			PairInt(-1,0),
			PairInt(-1,1),
			PairInt(0,1),
			PairInt(1,1)
		},
			L0{ PairInt(1,0) };
		for (int i = 0; i < stencilData.numThetaCells; i++) {
			double theta = 2 * M_PI * double(i) 
				/ double(stencilData.numThetaCells);
			// Theta is the NEGATIVE angle of drift b
			// stencil will be refined more finely in this direction
			PairDoub bhat{
				cos(theta),
				sin(theta)
			};
			std::vector<PairInt>
				L{ L0 },
				M{ M0 };
			while (!M.empty()) {
				PairInt u = M.back(),
					v = L.back();
				PairDoub uhat{
					(double)u.xi / u.norm(),
					(double)u.yi / u.norm()
				};
				PairDoub vhat{
					(double)v.xi / v.norm(),
					(double)v.yi / v.norm()
				};
				double uhvh = dot(uhat, vhat);
				double temp2 = stencilData.alpha *
					std::max(
						dot(uhat, bhat),
						dot(vhat, bhat)
					);
				if ((dot(uhat, vhat) < stencilData.alpha *
					std::max(
						dot(uhat, bhat),
						dot(vhat, bhat)
					)
					)
					&& (abs(u.xi + v.xi) + abs(u.yi + v.yi)
						<= stencilData.stencilCutoff)) {
					M.push_back(u + v);
				}
				else {
					L.push_back(M.back());
					M.pop_back();
				}
			}
			L.pop_back();
			stencilCollection.push_back(L);
		}
		break;
	}

	// ARCHIVED -- do not use
	// original version of modifed causal stencils
	// instead, use OBLONG
	case BUBBLE: {
		// first initialize all stencils to the 4-point stencilCollection
		stencilCollection = std::vector<std::vector<PairInt>>(
			stencilData.numThetaCells,
			std::vector<PairInt>{
				PairInt(1, 0), 
				PairInt(0, 1), 
				PairInt(-1, 0), 
				PairInt(0, -1)
			}
		);
		double delta = .1;
		int stencSize = stencilData.theta_res.size();
		std::vector<int> indices_reordered(stencSize, 0);
		int midInd = stencSize / 2; // supposed to be truncated
		for (int i = 0; i < stencilData.numThetaCells; i++) {
			stencilCollection.push_back(std::vector<PairInt>());
			double theta_bneg = 2 * M_PI * double(i) 
				/ double(stencilData.numThetaCells);
			double prev_dist = 0;
			for (int k = 0; k < stencSize; k++) {
				int ind = INT_MAX;
				if (k <= midInd) 
					ind = midInd - k;
				else 
					ind = k;
				if (k == midInd + 1) 
					prev_dist = 0; 
				// reset for the 2nd half of points
				double theta = std::fmod(
					stencilData.theta_res[ind] + theta_bneg, 
					2 * M_PI
				);
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
						double dist = abs(
							dot(
								curr, 
								PairDoub(
									cos(theta + M_PI / 2), 
									sin(theta + M_PI / 2)
								)
							)
						);
						if (dist < delta 
							&& (curr.norm() > std::min(
								prev_dist, 
								(double)stencilData.stencilCutoff
							))) {
							stencilCollection[i].push_back(curr);
							prev_dist = curr.norm();
							found = true;
							break;
						}
						else if (dist < delta) {
							stencilCollection[i].push_back(curr);
						}
						curr = curr + increment;
					}
					rad++;
				}

			}
			for (auto& pair : mesh.eightPtNeighborShifts)
				stencilCollection[i].push_back(pair);
			//std::cout << "\n Stencil Size is "<< stencilCollection[i].size();
		}
		break;
	}

	// oblong causal stencils, say paper "An efficient jet marcher for 
	// computing the quasi-potential
	case OBLONG: {
		double L = 12.0;
		double ell = 3.0;
		int N = 5;
		std::vector<PairDoub> coords;

		for (int i = 0; i < 2 * N - 1; i++) {
			double legLength = 0;
			double theta = 0;
			if (i < N - 1) {
				theta = M_PI / 2.0 * pow(0.5, i + 1);
				legLength = (i + 1) * (L / static_cast<double>(N));
			}
			else if (i == N - 1) {
				theta = 0.0;
				legLength = (i + 1) * (L / static_cast<double>(N));
			}
			else {
				theta = -M_PI / 2.0 * pow(0.5, 2 * N - (i + 1));
				legLength = (2 * N - (i + 1)) * (L / static_cast<double>(N));
			}
			int numSamples = ceil(legLength / ell);
			for (int j = 1; j <= numSamples; j++) {
				double curLen = static_cast<double>(j)* legLength 
					/ static_cast<double>(numSamples);
				coords.push_back(
					PairDoub(
						curLen * cos(theta), 
						curLen * sin(theta)
					)
				);
			}
		}
		for (int i = 0; i < stencilData.numThetaCells; i++) {
			double theta_bneg = 2 * M_PI * double(i) 
				/ double(stencilData.numThetaCells);
			std::vector<PairDoub> coords_rot;
			std::vector<PairInt> coords_rounded;
			for (int j = 0; j < coords.size(); j++) {
				coords_rot.push_back(coords[j].rotate(theta_bneg));
				coords_rounded.push_back(
					PairInt(
						round(coords_rot[j].x), 
						round(coords_rot[j].y)
					)
				);

			}
			for (auto& pair : mesh.eightPtNeighborShifts)
				coords_rounded.push_back(pair);
			std::sort(coords_rounded.begin(), coords_rounded.end());
			coords_rounded.erase(
				std::unique(
					coords_rounded.begin(), 
					coords_rounded.end()
				), 
				coords_rounded.end()
			);
			stencilCollection.push_back(coords_rounded);
		}
		break;
	}
	}
	// OUT OF SWITCH

	// This next loop precomputes the direction of b for all mesh points
	// and saves the z2_i of the appropriate stencilCollection.
	// h = the width of a theta bin = 2*PI/n.
	// Bins are [ (2n-1) h/2,h/2),[h/2,3h/2],...,[(2n-3) h/2,(2n-1)h/2], where
	// convention will be the interval is inclusive on the left.
	double h = 2 * M_PI / stencilData.numThetaCells;
	for (int index = 0; index < mesh.gridLength(); index++) {
		PairDoub b_neg = -1.0 * drift(index);
		double theta = b_neg.angle();
		stencilIndex[index] 
			= int((theta + h / 2) / h) % stencilData.numThetaCells;
	}
	std::cout << std::endl;
}

// update the neighbors of the newly minted Accepted point
void FastMarchEJM::updateNeighbors(
	int z2_i
) {
	int theta_i = stencilIndex[z2_i]; // stencil theta bin of z2_i
	PairInt z2_ip = getIndexPair(z2_i);

	for (int i = 0; 
		i < static_cast<int64_t>(stencilCollection.at(theta_i).size());
		i++) 
	{
		PairInt z1_ip{ getIndexPair(z2_i) - stencilCollection[theta_i][i]};
		int z1_i = getIndex(z1_ip); 
		if (isInvalidIndexPair(z1_ip) 
			|| isInvalidIndex(z1_i) 
			|| group[z1_i] == ACCEPTED)	
			continue;
		const std::vector<PairInt> *shifts{ nullptr };
		std::vector<PairInt> shifts_noBox;
		if (stencilData.boxUpdate) // use z3 as 8 pt neighbor shifts of z2 
			shifts = &mesh.eightPtNeighborShifts;
		else {	// use z3 as stencil neighbor of z2
			shifts_noBox.push_back(
				z1_ip 
				- z2_ip 
				+ stencilCollection[theta_i]
					[(i + (int)1) % stencilCollection[theta_i].size()]);
			shifts_noBox.push_back(
				z1_ip 
				- z2_ip 
				+ stencilCollection[theta_i]
					[(i + stencilCollection[theta_i].size() - 1) 
					% stencilCollection[theta_i].size()
					]
			);		
			shifts = &shifts_noBox;
		}
		onePointUpdate(z1_i, z2_i); // Precautionary 1-pt update
		for (const PairInt& shift : *shifts) {
			int z3_i = getIndex(z2_ip + shift);
			if (z1_i ==trackInd) {
				(getIndexPair(z2_i) - getIndexPair(z1_i)).print();
				std::cout << "\t";
				(getIndexPair(z3_i) - getIndexPair(z1_i)).print();
				std::cout << std::endl;
			}


			if (isInvalidIndex(z3_i))	
				continue;
			else if (group[z3_i] == ACCEPTED ) {
				if (options.z3one) // if settings wants z3 1-pt updates ran
					onePointUpdate(z1_i, z3_i); 
				triangleUpdate(z1_i, z2_i, z3_i);
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

// run a step of EJM unless the kill condition is met
void FastMarchEJM::runStep(bool& killCondition) {
	int ind = considered.takeMin();
	if (isInBoundary(ind)) {
		killCondition = true;
	}
	// In the event that 'ind' was last updated with a 1-point update, we now
	// run the fail-safe to try to get a successful triangle update.
	int iterFS = 1;
	if (updateType[ind] == 1) 
		failSafeCalls++;
	while ( options.failSafe 
		&& updateType[ind] == 1 
		&& iterFS < 20) 
	{
		fixFailures(ind,iterFS);
		iterFS += 1;
	}
	group[ind] = ACCEPTED;

	accepted.push_back(ind);
	if (options.debugMode)
		debugInfo[ind].rank = (int)accepted.size();
	if (killCondition) {
		return;
	}
	updateNeighbors(ind);
}

// write the stencil corresponding to mesh point ind to txt file
void FastMarchEJM::writeStencilToTXT(int index) {
	std::ofstream outFile("outputs/stencil.txt");
	outFile << "index,xi,yi,x,y\n";
	outFile << index << "," 
		<< getXIndex(index) << "," 
		<< getYIndex(index) << ","
		<< getXValue(index) << "," 
		<< getYValue(index) << std::endl;
	int si = stencilIndex[index];
	PairInt coords = getIndexPair(index);
	for (const PairInt& shift : stencilCollection[si]) {
		PairInt neighbor = coords + shift;
		if (isInvalidIndexPair(neighbor)) continue;
		int neighborIndex = getIndex(neighbor);
		outFile
			<< neighborIndex << ","
			<< neighbor.xi << ","
			<< neighbor.yi << ","
			<< getCoords(neighborIndex).x << ","
			<< getCoords(neighborIndex).y
			<< std::endl;
	}
	outFile.close();
}

// write solution and stencil to txt file
void FastMarchEJM::writeSolutionToTXT() {
	int index = 7;
	FastMarch::writeSolutionToTXT();
	writeStencilToTXT(index);
}

