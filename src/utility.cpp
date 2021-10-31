/*
	\file:		utility.cpp
	\brief:		This file defines run and diagnostic functions for main.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include <iostream>
#include "utility.h"
#include "fastmarch.h"
#include "fastmarchEJM.h"
#include "fastmarchOLIM.h"
#include "updaterBase.h"
#include "updaterLinear.h"
#include "updaterHermite.h"

// run Fast March object for single mesh and set of settings 
void runSingle(
	const MeshInfo& mesh,  
	UpdaterBase& upIn,
	const RunOptions& opIn, 
	const StencilData& stIn,
	const SpeedInfo& spIn, 
	int olimK, 
	char algKey)
{
	std::unique_ptr<FastMarch> fastMarcher(nullptr);
	switch (algKey) {
	case FM:
		fastMarcher.reset(new FastMarch(mesh, upIn, opIn, spIn));
		break;
	case EJM:
		fastMarcher.reset(new FastMarchEJM(mesh, upIn, opIn, spIn, stIn));
		break;
	case OLIM:
		fastMarcher.reset(new FastMarchOLIM(mesh, upIn, opIn, spIn, olimK));
		break;
	}
	fastMarcher->runMarch();
	if (opIn.debugMode) {
		fastMarcher->writeSolutionToTXT();
		fastMarcher->writeDebugInfotToTXT();
		std::cout << "Fail-Safe Calls:\t" 
			<< fastMarcher->failSafeCalls 
			<< std::endl;
		//fastMarcher->shootMaps(PairDoub(-1.333, 0));
	}
}

// run accuracy diagnostic to determine error convergence rate for single 
// set of settings
void convergenceRate(
	UpdaterBase& upIn, 
	RunOptions opIn,
	const StencilData& stIn, 
	const SpeedInfo& spIn, 
	int olimK, 
	int algKey,
	bool bulkWrite, 
	char quadType)
{
	opIn.debugMode = false;
	opIn.verbose = false;
	std::cout << "Error statistics for varying step size n:\n";

	switch (algKey) {
	case OLIM:
		std::cout << "OLIM algorithm ";
		break;
	case EJM:
		std::cout << "ASR algorithm ";
		break;
	}
	if (quadType == 'h') {
		std::cout << "with Hermite cubic quadrature\n";
	}
	else if (quadType == 'm') {
		std::cout << "with Midpoint linear quadrature\n";
	}
	std::cout << "Drift key: " << spIn.sp.switchKey << " with parameters a = "
		<< spIn.sp.a << ", b = " << spIn.sp.b << "\n\n";


	//std::vector<int> ns{ 17, 33, 65, 129, 257, 513, 1025, 2049};
	std::vector<int> ns{ 17, 33, 65, 129, 257, 513, 1025, 2049, 4097 };
	std::vector<int>olim_k{ 5,7,8,8,9,11,13,16,20 };
	//std::vector<int> ns{ 17, 33, 65, 129};
	std::string out_file_name = "Outputs/";
	if (bulkWrite) out_file_name += "ConvergencePlots/";


	if (algKey == EJM && quadType == 'h') 
		out_file_name += "EJM_";
	else if (algKey == OLIM && quadType == 'm') 
		out_file_name += "MidpointOLIM_";
	else if (algKey == OLIM && quadType == 'h') 
		out_file_name += "HermiteOLIM_";
	else if (algKey == EJM && quadType == 'm') 
		out_file_name += "MidpointASR_";
	out_file_name += std::to_string(spIn.sp.switchKey) 
		+ "_" + std::to_string((int)spIn.sp.b);
	out_file_name += ".txt";

	std::ofstream sol_out(out_file_name);
	// DATA COLUMNS: (1) n, (2) ErrMax, (3) ErrRMS, (4) ErrMax/UMax, 
	// (5) ErrRMS/URMS (6) ErrGradMax, (7) ErrGradRMS, 
	// (8) CPU time, (9) ?, (10) ?
	std::cout.precision(2); std::cout << std::scientific;
	std::cout << "n\tErrMax\t\tErrRMS\t\tErrMaxGrad\t"
		<< "ErrRMSGrad\tErrMaxLap\tErrRMSLap\tProcTime\n";
	int iteration = 0;
	for (int n : ns) {
		MeshInfo mesh(n, -4.0 / 3, 4.0 / 3, n, -2, 2);
		clock_t time_s = clock();
		std::vector<double> stats;
		//double shootError = 0;
		switch (algKey) {
		case FM: {
			FastMarch fastMarcher(mesh, upIn, opIn, spIn);
			fastMarcher.runMarch();
			fastMarcher.computeStatsOLIMComparison(&stats);
			break; }
		case EJM: {
			FastMarchEJM fastMarcher(mesh, upIn, opIn, spIn, stIn);
			fastMarcher.runMarch();
			fastMarcher.computeStatsOLIMComparison(&stats);
			fastMarcher.computeLaplacian();
			//shootError = fastMarcher.shootMaps(PairDoub(-1.333, 0));
			int xi = 1, yi = (n - 1) / 2;
			int ind = xi + mesh.nx * yi;
			double lp = fastMarcher.laplaceU
				[fastMarcher.getIndex(PairInt(xi, yi))];
			std::cout << abs(lp - fastMarcher.speeds.laplaceSolution(
				fastMarcher.getCoords(ind), 
				&fastMarcher.speeds.sp)
			) 
				<< std::endl;
			break; }
		case OLIM: {
			//olimK = 10 + 4 * (static_cast<int>(std::log2(mesh.nx)) - 7);
			//olimK = olim_k[iteration];
			FastMarchOLIM fastMarcher(mesh, upIn, opIn, spIn, olimK);
			fastMarcher.runMarch();
			fastMarcher.computeStatsOLIMComparison(&stats);
			//shootError = fastMarcher.shootMaps(
						// PairDoub(0.9769818, 0.02923383));
			//shootError = fastMarcher.shootMaps(PairDoub(-1.333, 0));

			break; }
		}
		clock_t time_f = clock();
		double procTime = (double(time_f) - time_s) / CLOCKS_PER_SEC;
		std::cout << n << "\t" 
			<< stats[0] << "\t" 
			<< stats[1] << "\t" 
			<< stats[4] << "\t" 
			<< stats[5] << "\t" 
			<< stats[6] << "\t" 
			<< stats[7] << "\t" 
			<< procTime  
			<< std::endl;

		sol_out << n << "\t" 
			<< stats[0] << "\t" 
			<< stats[1] << "\t" 
			<< stats[2] << "\t" 
			<< stats[3] << "\t" 
			<< stats[4] << "\t" 
			<< stats[5] << "\t" 
			<< procTime << "\t"
			<< std::endl;

		iteration++;
	}
	sol_out.close();

}

// run accuracy diagnostic to determine error convergence rate for several 
// sets of settings
void fullConvRate(
	UpdaterBase& upIn, 
	const RunOptions& opIn, 
	const StencilData& stIn, 
	int olimK
) {
	std::vector<DriftParams> sp_vec{ 
		DriftParams(5,2,.1), 
		DriftParams(5,2,1), 
		DriftParams(5,2,10) 
	};
	bool quasi(true), 
		verbose(true), 
		twoPtUpdateMinDist(true),
		failSafe(true), 
		denseInit(true), 
		debugMode(true);
	bool initType(false); // int initType(0)
	RunOptions op_ASR_Herm{ true,false,false,twoPtUpdateMinDist,
		failSafe,initType };
	RunOptions op_OLIM_Herm{ true,false,false,twoPtUpdateMinDist,
		failSafe,initType };
	RunOptions op_ASR_End{ true,false,false,!twoPtUpdateMinDist,
		!failSafe,initType };
	op_ASR_End.z3one = true;
	RunOptions op_OLIM_Mid{ true,false,false,!twoPtUpdateMinDist,
		!failSafe,initType };
	UpdaterLinearMidpoint mid;
	UpdaterLinearEndpoint end;
	UpdaterHermite herm;
	for (auto& sp : sp_vec) {
		const SpeedInfo speeds(sp);
		convergenceRate(
			end, 
			op_ASR_End, 
			stIn, 
			speeds, 
			olimK, 
			EJM, 
			true, 
			'e');
		convergenceRate(
			mid, 
			op_OLIM_Mid, 
			stIn, 
			speeds, 
			olimK, 
			OLIM, 
			true, 
			'm');

	}
}

// run OLIM for a variety of different k values
void olimKTest(
	const MeshInfo& mesh,  
	UpdaterBase& upIn, 
	const SpeedInfo& spIn, 
	RunOptions& opIn
) {
	opIn.debugMode = false;
	std::string out_file_name = "Outputs/OLIM_K/olimTestK_";
	out_file_name += std::to_string(spIn.sp.switchKey);
	out_file_name += "_";
	std::cout.precision(2); std::cout << std::scientific;
	std::cout << "Running K Test for n = " << mesh.nx 
		<< std::endl 
		<< std::endl;
	out_file_name += std::to_string(mesh.nx);
	std::cout << "K\tErrMax\t\tErrRMS\t\tErrMax/UMax\t"
		<< "ErrRMS/URMS\tErrMaxGrad\tErrRMSGrad\tProcTime\n";
	for (int k = 1; k < 42; k++) {
		out_file_name += "_" + std::to_string(k) + ".csv";
		std::ofstream sol_out(out_file_name);
		clock_t time_s = clock();
		FastMarchOLIM fastMarcher(mesh, upIn, opIn, spIn, k);
		fastMarcher.runMarch();
		std::vector<double> stats;
		fastMarcher.computeStatsOLIMComparison(&stats);
		clock_t time_f = clock();
		double procTime = (double(time_f) - time_s) / CLOCKS_PER_SEC;
		std::cout 
			<< k << "\t" 
			<< stats[0] << "\t" 
			<< stats[1] << "\t" 
			<< stats[2]	<< "\t" 
			<< stats[3] << "\t" 
			<< stats[4] << "\t" 
			<< stats[5] << "\t" 
			<< procTime 
			<< std::endl;
		sol_out 
			<< k << "\t" 
			<< stats[0] << "\t" 
			<< stats[1] << "\t" 
			<< stats[2] << "\t" 
			<< stats[3] << "\t" 
			<< stats[4] << "\t" 
			<< stats[5] << "\t" 
			<< procTime 
			<< std::endl;
		sol_out.close();
	}




}




//void checkSpeedDerivsCorrectness(const DriftParams* sp, bool quasi) {
//	PairDoub z(gauss(), gauss()), v(gauss(), gauss());
//	double u2(gauss()), u3(gauss());
//
//	double delta = 1e-7;
//	PairDoub e1(1, 0), e2(0, 1);
//
//	if (!quasi) {
//		std::cout << "Eikonal equation slowness derivative checks:\n\n";
//		std::cout << "Finite difference value vs. derivative value\n"
//			<< "ds/dx: " << masterSlow(z + delta * e1, v, sp) - masterSlow(z, v, sp)
//			<< " vs. " << masterGradZ(z, v, sp).x * delta << std::endl
//			<< "ds/dy: " << masterSlow(z + delta * e2, v, sp) - masterSlow(z, v, sp)
//			<< " vs. " << masterGradZ(z, v, sp).y * delta << std::endl;
//
//	}
//	else {
//		std::cout << "Quasipotential drift derivative checks:\n\n";
//		std::cout << "Finite difference value vs. derivative value\n"
//			<< "db1/dx: " << masterDrift(z + delta * e1, sp).x - masterDrift(z, sp).x
//			<< " vs. " << masterGradDrift(z, sp).row1.x * delta << std::endl
//			<< "db1/dy: " << masterDrift(z + delta * e2, sp).x - masterDrift(z, sp).x
//			<< " vs. " << masterGradDrift(z, sp).row1.y * delta << std::endl
//			<< "db2/dx: " << masterDrift(z + delta * e1, sp).y - masterDrift(z, sp).y
//			<< " vs. " << masterGradDrift(z, sp).row2.x * delta << std::endl
//			<< "db2/dy: " << masterDrift(z + delta * e2, sp).y - masterDrift(z, sp).y
//			<< " vs. " << masterGradDrift(z, sp).row2.y * delta << std::endl
//			<< "\nSecond Order Derivatives\n"
//			<< "db1/dxx: " << masterGradDrift(z + delta * e1, sp).row1.x - masterGradDrift(z, sp).row1.x
//			<< " vs. " << masterHessDrift(z, sp).mat1.row1.x * delta << std::endl
//			<< "db1/dxy: " << masterGradDrift(z + delta * e2, sp).row1.x - masterGradDrift(z, sp).row1.x
//			<< " vs. " << masterHessDrift(z, sp).mat1.row1.y * delta << std::endl
//			<< "db1/dyx: " << masterGradDrift(z + delta * e1, sp).row1.y - masterGradDrift(z, sp).row1.y
//			<< " vs. " << masterHessDrift(z, sp).mat1.row2.x * delta << std::endl
//			<< "db1/dyy: " << masterGradDrift(z + delta * e2, sp).row1.y - masterGradDrift(z, sp).row1.y
//			<< " vs. " << masterHessDrift(z, sp).mat1.row2.y * delta << std::endl
//			<< "db2/dxx: " << masterGradDrift(z + delta * e1, sp).row2.x - masterGradDrift(z, sp).row2.x
//			<< " vs. " << masterHessDrift(z, sp).mat2.row1.x * delta << std::endl
//			<< "db2/dxy: " << masterGradDrift(z + delta * e2, sp).row2.x - masterGradDrift(z, sp).row2.x
//			<< " vs. " << masterHessDrift(z, sp).mat2.row1.y * delta << std::endl
//			<< "db2/dyx: " << masterGradDrift(z + delta * e1, sp).row2.y - masterGradDrift(z, sp).row2.y
//			<< " vs. " << masterHessDrift(z, sp).mat2.row2.x * delta << std::endl
//			<< "db2/dyy: " << masterGradDrift(z + delta * e2, sp).row2.y - masterGradDrift(z, sp).row2.y
//			<< " vs. " << masterHessDrift(z, sp).mat2.row2.y * delta << std::endl;
//
//
//
//	}
//
//}


// Pre-saved drift and slowness functions.
/*
	For s(z,v) = |b(z)| - <b(z),v> / |v|, we calculate the following.
		ds/dv = v <b,v> /|v|^3 - b/|v|
		d^2s/dv^2 = TODO
		ds/dz = b Db/|b| - v Db
		d^2s/dz^2 = TODO
*/



