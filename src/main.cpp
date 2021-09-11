#include "fastmarch.h"
#include "fastmarchASR.h"
#include "fastmarchOLIM.h"
#include "structs.h"
#include <cmath>
#include <tuple>
#include <vector>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>



void convergenceRate(const Updater& upIn,  RunOptions opIn, const StencilData& stIn, const SpeedInfo &spIn, int olimK, int algKey,bool bulkWrite, char quadType);
void runSingle(const GridInfo& grid, const Updater& upIn, const RunOptions& opIn, const StencilData& stIn, const SpeedInfo& spIn, int olimK, char algKey);
void checkSpeedDerivsCorrectness(const SpeedParams* sp, bool quasi);
void fullConvRate(const Updater& upIn, const RunOptions& opIn, const StencilData& stIn, int olimK);
void olimKTest(const GridInfo& grid, const Updater& upIn, const SpeedInfo& spIn,  RunOptions& opIn);

// Temporary global variables that I should delete later.
int FAIL = 0, SUCCESS = 0;
double xfac = 1, yfac = 1;

int main() {
	/* Options description:

	Run Type Key:
			's' = single run
			'c' = run and write analysis data for varying N values 
					and fixed other settings
			'z' = run and write analysis data for varying N values 
					and varying other settings
			'b' = analyze efficacy of different types of stencils for EJM
			'k' = analyze effectiveness of different K values for OLIM

	Algorithm Key: (variable algKey)
			FM = Traditional fast march
			ASR = Efficient jet marcher
			OLIM = Ordered line integral method

	Speed Key:		(with parameters a & b)
			0: (Q)		drift = [- a*x -b*y, a*b*x - y]
			1:			slowness = 1
			2:			slowness = 1 + x^2/2 - cos(x+y)
			3: (Q)		drift = [-ax^3 - by^3, abx^3 - y^3]
			4: (Q)		drift = [-x-y-b*sin(x)*cos(x), x - y + b*sin(x)*cos(x)]
			5: (Q)		NONLINEAR DIAGNOSTIC DRIFT
							drift = [-2x - 0.75*a*x^2 - by, 
							2*b*x + 0.75*a*b*x^2 - y]
			6: (Q)		FITZHUGH-NAGUMO
							drift = [-x *(a^2-1 - a*x + 1/3*x^2) - y, bx]
			7: (Q)		LINEARIZED FITZHUGH-NAGUMO
							drift = [-x*(a^2-1) - y, bx]
			8: (Q)		MAIER-STEIN
							drift = [-2x - x^3 - b*x*y^2 + 3x^2 + by^2, 
							-2y - x^2*y + 2x*y]
		*Note: Q denotes quasi-potential only.

	Quadrature Key:
			'e' = endpoint quadrature
			'm' = midpoint rule
			's' = Simpson's rule
			'h' = Hermite interpolation with Simpsons' rule quadrature

	Updating Settings: (should set all to true for EJM and all to false for OLIM
			'twoPtUpdateMinDist':
				true -- 2 pt updates are accepted if new hypotenuse is smaller
				false -- 2 pt updates are accepted if new value is smaller
			'failSafe': turn on to use failSafe
			'fakeFilter': 
				true -- new updates are accepted only if they propose smaller
						values than	all previous 1 pt updates of the point
				false -- usual fast march rules

	Initialization Settings: Possible values of variable initType
			0 -- Initialize u to linearized solution on a box around attractor.
				Default is 8 pt neighborhood, to change modify functions
				FastMarch::initialPoint and FastMarch::initialValue
			1 -- Initialize u to exact solution on a box around attractor.
			2 -- Initialize u to linearized solution 


	*/

	
	
	char runType = 's';
	int algKey = ASR;

	//********** GRID SPECIFICATIONS FOR SINGLE RUN **********
	int n = 257;
	const GridInfo grid(n, -1, 1, n, -1, 1);
	//const GridInfo grid(n, -4.0/3, 4.0/3, n, -2, 2);
	// uniform rectangular mesh is specified by 
	//		grid(nx, xleft, xright, ny, ybottom, ytop)



	//********** ADDITIONAL OPTIONS **********
	bool quasi = true; // Turn off for eikonal examples
	bool verbose = true; // Turn off to suppress console output
	bool debugMode = true; // Turn off to not store detailed debug info (saves memory)

	// Updating Settings
	bool twoPtUpdateMinDist = false;
	bool failSafe = false;
	bool fakeFilter = false;
	if (algKey == ASR) {
		twoPtUpdateMinDist = true;
		failSafe = true;
		fakeFilter = true;
	}

	// Initialization Settings
	// 0 = linear initalize values on 8 neighbor grid of (0,0) 
	// 1 = run fast march on refined neighborhood around (0,0)
	int initType = 0;

	RunOptions opt{quasi, verbose, debugMode, twoPtUpdateMinDist,
		failSafe, fakeFilter, initType };


	//********** SPEED & ALGORITHM SETTINGS FOR SINGLE RUN AND DIAGNOSTIC **********
	char quadKey = 'h';
	if (algKey == ASR) quadKey = 'h';
	std::unique_ptr<Updater> updater(nullptr);
	switch (quadKey) {
	case 'e':
		updater.reset(new EndpointLinear());
		break;
	case 'm':
		updater.reset(new MidpointLinear());
		break;
	case 'h':
		updater.reset(new Hermite());
		break;
	}

	int speedKey =8;
	double a = 1.05, b = .1;
	a = 2.0, b = 3;
	const SpeedParams spdPrms(speedKey, a, b);
	const SpeedInfo speeds(spdPrms, quasi);	


	//********** ASR ADDITIONAL SETTINGS **********
	double alpha = 0.999;
	int maxStencilSize = 20;
	int numStencilThetaBins = 1000;
	bool boxUpdateASR = true;
	bool bubble = false;
	const StencilData stencilData{ alpha,maxStencilSize,numStencilThetaBins,boxUpdateASR, false };
	const StencilData stencilDataBubble{ alpha,maxStencilSize,numStencilThetaBins,boxUpdateASR, true };
	const StencilData stencilDataBubbleNew{ alpha,maxStencilSize,numStencilThetaBins,boxUpdateASR, true };

	   	 
	//********** OLIM ADDITIONAL SETTINGS **********
	int olimK =26;	

	// Uncomment next line to check if drift derivative formulas are correct
	//checkSpeedDerivsCorrectness(&spdPrms,true); 


	switch (runType) {
	case 's': {
		opt.verbose = true;
		runSingle(grid, *updater, opt, stencilData, speeds, olimK, algKey);
		break; }
	case 'c': {
		convergenceRate(*updater, opt, stencilData, speeds, olimK, algKey,true, quadKey);
		break; }
	case 'z': {
		fullConvRate(*updater, opt, stencilData, olimK);
		break;
	}
	case 'b': { // bubble
		std::ofstream sol_out("Outputs/testBubble.csv");
		sol_out << "n, 1-pointers Mirebeau, 1-pointers balloon" << std::endl;
		std::cout << "n\t 1-pointers Mirebeau\t 1-pointers balloon" << std::endl;
		int kmax = 10, kmin = 6, j = 1;
		for (int k = kmin; k <= kmax; k++) {
			opt.verbose = false;
			int n = pow(2, k) + 1;
			const GridInfo grid(n, -1, 1, n, -1, 1);
			FastMarchASR bub(grid, *updater, opt, speeds, stencilDataBubble);
			FastMarchASR nobub(grid, *updater, opt, speeds, stencilData);
			sol_out << n << ",";
			std::cout << n << "\t,";
			nobub.runMarch();
			sol_out << nobub.stats.num_acc_1pt << ",";
			std::cout << nobub.stats.num_acc_1pt << "\t";
			bub.runMarch();
			sol_out << bub.stats.num_acc_1pt << "\n";
			std::cout << bub.stats.num_acc_1pt << "\n";
			j++;
		}

		sol_out.close();
		break; 
	}
	case 'k': {
		opt.verbose = false;
		olimKTest(grid, *updater, speeds, opt);
		break;
	}
	}
} 

/*
*******************************************************************************
Define Run Types
*******************************************************************************
*/
void runSingle(const GridInfo& grid, const Updater &upIn, 
	const RunOptions &opIn, const StencilData &stIn, 
	const SpeedInfo &spIn, int olimK, char algKey) 
{
	std::unique_ptr<FastMarch> fastMarcher(nullptr);
	switch (algKey) {
	case FM:
		fastMarcher.reset(new FastMarch(grid, upIn, opIn, spIn));
		break;
	case ASR:
		fastMarcher.reset(new FastMarchASR(grid, upIn, opIn, spIn, stIn));
		break;
	case OLIM:
		fastMarcher.reset(new FastMarchOLIM(grid, upIn, opIn, spIn, olimK));
		break;
	}
	fastMarcher->runMarch();
	if (opIn.debugMode) {
		fastMarcher->writeToTXT();
		fastMarcher->writeDebug();
		std::cout << "Fail-Safe Calls:\t" << fastMarcher->failSafe << std::endl;
		//fastMarcher->shootMaps(PairDoub(0.9984693, 0.00195196));
		//fastMarcher->shootMaps(PairDoub(0.9769818, 0.02923383));
		//fastMarcher->shootMaps(PairDoub(-.01758, 0.232422));
		fastMarcher->shootMaps(PairDoub(-1.333, 0));
		//fastMarcher->printDebugConsole();
	}
}

void convergenceRate(const Updater& upIn, RunOptions opIn, 
	const StencilData& stIn, const SpeedInfo &spIn, int olimK, int algKey, 
	bool bulkWrite, char quadType) 
{
	opIn.debugMode = false;
	opIn.verbose = false;
	std::cout << "Error statistics for varying step size n:\n";

	switch (algKey) {
	case OLIM:
		std::cout << "OLIM algorithm ";
		break;
	case ASR:
		std::cout << "ASR algorithm ";
		break;
	}
	if (quadType == 'h') { std::cout << "with Hermite cubic quadrature\n"; }
	else if (quadType == 'm') {	std::cout << "with Midpoint linear quadrature\n"; }
	std::cout << "Drift key: " << spIn.sp.switchKey << " with parameters a = " 
		<< spIn.sp.a << ", b = " << spIn.sp.b << "\n\n";


	//std::vector<int> ns{ 17, 33, 65, 129, 257, 513, 1025, 2049};
	std::vector<int> ns{ 17, 33, 65, 129, 257, 513, 1025, 2049, 4097 };
	std::vector<int>olim_k{ 5,7,8,8,9,11,13,16,20 };
	//std::vector<int> ns{ 17, 33, 65, 129};
	std::string out_file_name = "Outputs/";
	if (bulkWrite) out_file_name += "ConvergencePlots/";


	if (algKey == ASR && quadType == 'h') out_file_name += "EJM_";
	else if (algKey == OLIM && quadType == 'm') out_file_name += "MidpointOLIM_";
	else if (algKey == OLIM && quadType == 'h') out_file_name += "HermiteOLIM_";
	else if (algKey == ASR && quadType == 'm') out_file_name += "MidpointASR_";
	out_file_name += std::to_string(spIn.sp.switchKey) + "_" + std::to_string((int)spIn.sp.b);
	out_file_name += ".txt";

	std::ofstream sol_out(out_file_name);
	// DATA COLUMNS: (1) n, (2) ErrMax, (3) ErrRMS, (4) ErrMax/UMax, (5) ErrRMS/URMS
	// (6) ErrGradMax, (7) ErrGradRMS, (8) CPU time, (9) ?, (10) ?
	std::cout.precision(2); std::cout << std::scientific;
	std::cout << "n\tErrMax\t\tErrRMS\t\tErrMaxGrad\tErrRMSGrad\tErrMaxLap\tErrRMSLap\tProcTime\n";
	int iteration = 0;
	for (int n : ns) {
		GridInfo grid(n, -4.0/3, 4.0/3, n, -2, 2);
		clock_t time_s = clock();
		std::vector<double> stats;
		double shootError = 0;
		switch (algKey) {
		case FM: {
			FastMarch fastMarcher(grid, upIn, opIn, spIn);
			fastMarcher.runMarch();
			fastMarcher.computeMashaStats(&stats);
			break; }
		case ASR: {
			FastMarchASR fastMarcher(grid, upIn, opIn, spIn, stIn);
			fastMarcher.runMarch();
			fastMarcher.computeMashaStats(&stats);
			fastMarcher.computeLaplacian();
			//shootError = fastMarcher.shootMaps(PairDoub(0.9769818, 0.02923383));
			shootError = fastMarcher.shootMaps(PairDoub(-1.333, 0));
			int xi = 1, yi = (n - 1) / 2 ;
			int ind = xi + grid.nx * yi;
			double lp = fastMarcher.laplaceU[fastMarcher.getIndex(PairInt( xi, yi))];
			std::cout << abs(lp - fastMarcher.speeds.laplaceSolution(fastMarcher.getCoords(ind), &fastMarcher.speeds.sp)) << std::endl;
			break; }
		case OLIM: {
			//olimK = 10 + 4 * (static_cast<int>(std::log2(grid.nx)) - 7);
			//olimK = olim_k[iteration];
			FastMarchOLIM fastMarcher(grid, upIn, opIn, spIn, olimK);
			fastMarcher.runMarch();
			fastMarcher.computeMashaStats(&stats);
			//shootError = fastMarcher.shootMaps(PairDoub(0.9769818, 0.02923383));
			shootError = fastMarcher.shootMaps(PairDoub(-1.333, 0));

			break; }
		}
		clock_t time_f = clock();
		double procTime = (double(time_f) - time_s) / CLOCKS_PER_SEC;
		std::cout << n << "\t" << stats[0] << "\t" << stats[1] <<  "\t" << stats[4] << "\t" << stats[5] << "\t" << stats[6] << "\t" << stats[7] << "\t" <<  procTime << "\t" << shootError << std::endl;

		sol_out << n << "\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2]
			<< "\t" << stats[3] << "\t" << stats[4] << "\t" << stats[5] << "\t" << procTime << "\t" 
			<< shootError << std::endl;

		iteration++;
	}
	sol_out.close();

}

void olimKTest(const GridInfo& grid, const Updater& upIn, const SpeedInfo &spIn,  RunOptions& opIn) {
	opIn.debugMode = false;
	std::string out_file_name = "Outputs/OLIM_K/olimTestK_";
	out_file_name += std::to_string(spIn.sp.switchKey);
	out_file_name += "_";
	std::cout.precision(2); std::cout << std::scientific;
	std::cout << "Running K Test for n = " << grid.nx << std::endl << std::endl;
	out_file_name += std::to_string(grid.nx);
	std::cout << "K\tErrMax\t\tErrRMS\t\tErrMax/UMax\tErrRMS/URMS\tErrMaxGrad\tErrRMSGrad\tProcTime\n";
	for (int k = 1; k < 42; k++) {
		out_file_name += "_" + std::to_string(k) + ".csv";
		std::ofstream sol_out(out_file_name);
		clock_t time_s = clock();
		FastMarchOLIM fastMarcher(grid, upIn, opIn, spIn, k);
		fastMarcher.runMarch();
		std::vector<double> stats;
		fastMarcher.computeMashaStats(&stats);
		clock_t time_f = clock();
		double procTime = (double(time_f) - time_s) / CLOCKS_PER_SEC;
		std::cout << k << "\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2]
			<< "\t" << stats[3] << "\t" << stats[4] << "\t" << stats[5] << "\t" << procTime << std::endl;
		sol_out << k << "\t" << stats[0] << "\t" << stats[1] << "\t" << stats[2]
			<< "\t" << stats[3] << "\t" << stats[4] << "\t" << stats[5] << "\t" << procTime << std::endl;
		sol_out.close();
	}


	

}


void fullConvRate(const Updater& upIn, const RunOptions& opIn, const StencilData& stIn, int olimK) {
	//std::vector<SpeedParams> sp_vec{ SpeedParams(0,2,1), SpeedParams(0,2,10), SpeedParams(3,2,1), SpeedParams(5,2,10), SpeedParams(5,2,1) };
	std::vector<SpeedParams> sp_vec{ SpeedParams(5,2,.1), SpeedParams(5,2,1), SpeedParams(5,2,10)};
	bool quasi(true), verbose(true), twoPtUpdateMinDist(true), 
		failSafe(true), denseInit(true), debugMode(true);
	bool initType(false); // int initType(0)
	RunOptions op_ASR_Herm{ true,false,false,twoPtUpdateMinDist,
		failSafe,initType};
	RunOptions op_OLIM_Herm{ true,false,false,twoPtUpdateMinDist,
		failSafe,initType };
	RunOptions op_ASR_End{ true,false,false,!twoPtUpdateMinDist,
		!failSafe,initType };
	op_ASR_End.z3one = true;
	RunOptions op_OLIM_Mid{ true,false,false,!twoPtUpdateMinDist,
		!failSafe,initType };
	const MidpointLinear mid;
	const EndpointLinear end;
	const Hermite herm;
	for (auto& sp : sp_vec) {
		const SpeedInfo speeds(sp, true);
		//convergenceRate(herm, op_ASR_Herm, stIn, speeds, olimK, ASR,true, 'h');
		convergenceRate(end, op_ASR_End, stIn, speeds, olimK, ASR, true, 'e');
		convergenceRate(mid, op_OLIM_Mid, stIn, speeds, olimK, OLIM, true, 'm');
		//convergenceRate(herm, op_OLIM_Herm, stIn, speeds, olimK, OLIM, true, 'h');

	}
}

void checkSpeedDerivsCorrectness(const SpeedParams* sp, bool quasi) {
	PairDoub z(gauss(), gauss()), v(gauss(), gauss());
	double u2(gauss()), u3(gauss());

	double delta = 1e-7;
	PairDoub e1(1, 0), e2(0, 1);

	if (!quasi) {
		std::cout << "Eikonal equation slowness derivative checks:\n\n";
		std::cout << "Finite difference value vs. derivative value\n"
			<< "ds/dx: " << masterSlow(z + delta * e1, v, sp) - masterSlow(z, v, sp)
			<< " vs. " << masterGradZ(z, v, sp).x * delta << std::endl
			<< "ds/dy: " << masterSlow(z + delta * e2, v, sp) - masterSlow(z, v, sp)
			<< " vs. " << masterGradZ(z, v, sp).y * delta << std::endl
			<< "d^2s/dx^2: " << masterGradZ(z + delta * e1, v, sp).x - masterGradZ(z, v, sp).x
			<< " vs. " << masterHessZ(z, v, sp).row1.x * delta << std::endl
			<< "d^2s/dy^2: " << masterGradZ(z + delta * e2, v, sp).y - masterGradZ(z, v, sp).y
			<< " vs. " << masterHessZ(z, v, sp).row2.y * delta << std::endl
			<< "d^2s/dxdy: " << masterGradZ(z + delta * e2, v, sp).x - masterGradZ(z, v, sp).x
			<< " vs. " << masterHessZ(z, v, sp).row2.x * delta << std::endl;
	}
	else {
		std::cout << "Quasipotential drift derivative checks:\n\n";
		std::cout << "Finite difference value vs. derivative value\n"
			<< "db1/dx: " << masterDrift(z + delta * e1, sp).x - masterDrift(z, sp).x
			<< " vs. " << masterGradDrift(z, sp).row1.x * delta << std::endl
			<< "db1/dy: " << masterDrift(z + delta * e2, sp).x - masterDrift(z, sp).x
			<< " vs. " << masterGradDrift(z, sp).row1.y * delta << std::endl
			<< "db2/dx: " << masterDrift(z + delta * e1, sp).y - masterDrift(z, sp).y
			<< " vs. " << masterGradDrift(z, sp).row2.x * delta << std::endl
			<< "db2/dy: " << masterDrift(z + delta * e2, sp).y - masterDrift(z, sp).y
			<< " vs. " << masterGradDrift(z, sp).row2.y * delta << std::endl
			<< "\nSecond Order Derivatives\n"
			<< "db1/dxx: " << masterGradDrift(z + delta * e1, sp).row1.x - masterGradDrift(z, sp).row1.x
			<< " vs. " << masterHessDrift(z, sp).mat1.row1.x * delta << std::endl
			<< "db1/dxy: " << masterGradDrift(z + delta * e2, sp).row1.x - masterGradDrift(z, sp).row1.x
			<< " vs. " << masterHessDrift(z, sp).mat1.row1.y * delta << std::endl
			<< "db1/dyx: " << masterGradDrift(z + delta * e1, sp).row1.y - masterGradDrift(z, sp).row1.y
			<< " vs. " << masterHessDrift(z, sp).mat1.row2.x * delta << std::endl
			<< "db1/dyy: " << masterGradDrift(z + delta * e2, sp).row1.y - masterGradDrift(z, sp).row1.y
			<< " vs. " << masterHessDrift(z, sp).mat1.row2.y * delta << std::endl
			<< "db2/dxx: " << masterGradDrift(z + delta * e1, sp).row2.x - masterGradDrift(z, sp).row2.x
			<< " vs. " << masterHessDrift(z, sp).mat2.row1.x * delta << std::endl
			<< "db2/dxy: " << masterGradDrift(z + delta * e2, sp).row2.x - masterGradDrift(z, sp).row2.x
			<< " vs. " << masterHessDrift(z, sp).mat2.row1.y * delta << std::endl
			<< "db2/dyx: " << masterGradDrift(z + delta * e1, sp).row2.y - masterGradDrift(z, sp).row2.y
			<< " vs. " << masterHessDrift(z, sp).mat2.row2.x * delta << std::endl
			<< "db2/dyy: " << masterGradDrift(z + delta * e2, sp).row2.y - masterGradDrift(z, sp).row2.y
			<< " vs. " << masterHessDrift(z, sp).mat2.row2.y * delta << std::endl;



	}

}


// Pre-saved drift and slowness functions.
/*
	For s(z,v) = |b(z)| - <b(z),v> / |v|, we calculate the following.
		ds/dv = v <b,v> /|v|^3 - b/|v|
		d^2s/dv^2 = TODO
		ds/dz = b Db/|b| - v Db
		d^2s/dz^2 = TODO
*/



PairDoub masterDrift(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return PairDoub(-params->a * z.x - params->b * z.y, params->a * params->b * z.x - z.y);
	case 1:	case 2:
		return PairDoub(INFINITY, INFINITY);
	case 3:
		return PairDoub(-params->a * z.x*z.x*z.x - params->b * z.y*z.y*z.y, 
			params->a * params->b * z.x*z.x*z.x - z.y*z.y*z.y);
	case 4:
		return PairDoub(-z.x - z.y -params->b*sin(z.x)*cos(z.x), z.x-z.y+params->b*sin(z.x)*cos(z.x));
	case 5:
		return PairDoub(-2 * z.x - 0.75 * params->a * z.x * z.x - params->b * z.y, 2 * params->b * z.x + .75 * params->a * params->b * z.x * z.x - z.y);
	case 6:
		return PairDoub(-z.x*xfac * (params->a * params->a - 1 - params->a * z.x*xfac + 1.0 / 3.0 * z.x*xfac * z.x*xfac) - z.y*yfac, params->b * z.x*xfac);
	case 7:
		return PairDoub(-z.x*xfac * (params->a * params->a - 1 ) - z.y*yfac, params->b * z.x*xfac);
	case 8:
		return PairDoub(-2 * z.x - z.x * z.x * z.x - params->b * z.x * z.y * z.y + 3 * z.x * z.x + params->b * z.y * z.y, -2 * z.y - z.x * z.x * z.y + 2 * z.x * z.y);

	}
}

MatrixDoub masterGradDrift(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0: 
		return MatrixDoub(PairDoub(-params->a, -params->b), PairDoub(params->a * params->b, -1));
	case 1: case 2:
		return MatrixDoub(PairDoub(INFINITY, INFINITY), PairDoub(INFINITY, INFINITY));
	case 3:
		return MatrixDoub(PairDoub(-3 * params->a * z.x * z.x, -3 * params->b * z.y * z.y),
			PairDoub(3 * params->a * params->b * z.x * z.x, -3*z.y * z.y));
	case 4:
		return MatrixDoub(PairDoub(-1 - params->b * cos(2 * z.x), -1), PairDoub(1 + params->b * cos(2 * z.x), -1));
	case 5:
		return MatrixDoub(PairDoub(-2 - 1.5 * params->a * z.x, -params->b), PairDoub(2 * params->b + 1.5 * params->a * params->b * z.x, -1));
	case 6:
		return MatrixDoub(PairDoub(xfac*(1 - params->a * params->a + 2 * params->a * xfac*z.x - xfac*z.x * xfac* z.x), -yfac), PairDoub(xfac*params->b, 0));
	case 7:
		return MatrixDoub(PairDoub(xfac*(1 - params->a * params->a) , -yfac), PairDoub(xfac*params->b, 0));
	case 8 :
		return MatrixDoub(PairDoub(-2 - 3 * z.x * z.x - params->b * z.y * z.y + 6 * z.x, -2 * params->b * z.x * z.y + 2 * params->b * z.y), PairDoub(-2 * z.x * z.y + 2 * z.y, -2 - z.x * z.x + 2 * z.x));
	}
}

TensorDoub masterHessDrift(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0: case 1: case 2: 
		return TensorDoub();
	case 3:
		return TensorDoub(MatrixDoub(PairDoub(-6*params->a*z.x,0), PairDoub(0,-6*params->b*z.y)), 
			MatrixDoub(PairDoub(6*params->a*params->b*z.x,0), PairDoub(0,-6*z.y)));
	case 4:
		return TensorDoub(MatrixDoub(PairDoub(2 * params->b * sin(2 * z.x), 0), PairDoub(0, 0)),
			MatrixDoub(PairDoub(-2 * params->b * sin(2 * z.x), 0), PairDoub(0, 0)));
	case 5:
		return TensorDoub(MatrixDoub(PairDoub(-1.5 * params->a, 0), PairDoub(0, 0)), MatrixDoub(PairDoub(1.5 * params->a * params->b, 0), PairDoub(0, 0)));
	case 6:
		return TensorDoub(MatrixDoub(PairDoub(xfac*xfac*(2 * params->a - 2 * xfac*z.x), 0), PairDoub(0, 0)), MatrixDoub(PairDoub(0, 0), PairDoub(0, 0)));
	case 7:
		return TensorDoub();
	case 8:
		return TensorDoub(MatrixDoub(PairDoub(-6 * z.x + 6, -2 * params->b * z.y), PairDoub(-2 * params->b * z.y, -2 * params->b * z.x + 2 * params->b)), MatrixDoub(PairDoub(-2 * z.y, -2 * z.x + 2), PairDoub(-2 * z.x + 2, 0)));
	}
}

double masterDivDrift(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;

	switch (params->switchKey) {
	default:
		return 0;
	case 0:
		return params->a * 2 + 2;
	case 5:
		return -3-3*z.x;
	case 6:
		return 0;
	case 8:
		return 2 * (z.x * z.x + z.y * z.y + .25 * z.x * z.x * z.x * z.x + 0.5 * z.x * z.x * z.y * z.y - z.x * z.x * z.x - z.x * z.y * z.y);
	case 7:
		double B11 = xfac * (1.0 - params->a * params->a), B12 = -yfac, B21 = xfac * params->b, B22 = 0.0;
		double aux1 = B21 - B12, aux2 = B11 + B22, aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return Q11 * z.x * z.x + (Q12 + Q21) * z.x * z.y + Q22 * z.y * z.y;
	}
}

double masterSol(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;

	switch (params->switchKey) {
	default:
	case 0: 
		return params->a * z.x * z.x + z.y * z.y;
	case 1:
		return sqrt(z.x * z.x + z.y * z.y);
	case 2:
		return z.x*z.x/2 + 1 -cos(z.x+z.y);
	case 3:
		return 0.5 * params->a * z.x * z.x * z.x * z.x + 0.5 * z.y * z.y * z.y * z.y;
	case 4:
		return z.x * z.x  +z.y * z.y + params->b * sin(z.x) * sin(z.x);
	case 5:
		return 2 * z.x * z.x + .5 * params->a * z.x * z.x * z.x + z.y * z.y;
	case 6: 
		return 0;
	case 8:
		return 2 * (z.x * z.x + z.y * z.y + .25 * z.x * z.x * z.x * z.x + 0.5 * z.x * z.x * z.y * z.y - z.x * z.x * z.x - z.x * z.y * z.y);
	case 7:
		double B11 =xfac*( 1.0 - params->a*params->a), B12 = -yfac, B21 = xfac* params->b, B22 = 0.0; 
		double aux1 = B21 - B12, aux2 = B11 + B22, aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return Q11 * z.x * z.x + (Q12 + Q21) * z.x * z.y + Q22 * z.y * z.y;

	}
}
PairDoub masterGradSol(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return PairDoub(2*params->a * z.x, 2* z.y);
	case 1:
		return PairDoub(z)/z.norm();
	case 2:
		return PairDoub(z.x + sin(z.x+z.y),sin(z.x+z.y));
	case 3:
		return PairDoub(2 * params->a * z.x * z.x * z.x, 2 * z.y * z.y * z.y);
	case 4:
		return PairDoub(2 * z.x + params->b * sin(2 * z.x), 2 * z.y);
	case 5:
		return PairDoub(4 * z.x + 1.5 * params->a * z.x * z.x, 2 * z.y);
	case 6: 
		return PairDoub(0, 0); 
	case 8:
		return PairDoub(2 * (2 * z.x + z.x * z.x * z.x + z.x * z.y * z.y - 3 * z.x * z.x - z.y * z.y), 2 * (2 * z.y + z.x * z.x * z.y - 2 * z.x * z.y));
	case 7:
		double B11 = xfac*(1.0 - params->a * params->a), B12 = -yfac, B21 = xfac*params->b, B22 = 0.0;
		double aux1 = B21 - B12, aux2 = B11 + B22, aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return PairDoub(2 * Q11 * z.x + (Q12 + Q21) * z.y, (Q12 + Q21) * z.y + 2 * Q22 * z.y);

	}
}

double masterLaplacianSol(const PairDoub& z, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;

	switch (params->switchKey) {
	default:
		return 0;
	case 0:
		return params->a * 2 + 2;
	case 5:
		return 6+ 3 * params->a * z.x ;
	case 6:
		return 0;
	case 8:
		return 2 * (z.x * z.x + z.y * z.y + .25 * z.x * z.x * z.x * z.x + 0.5 * z.x * z.x * z.y * z.y - z.x * z.x * z.x - z.x * z.y * z.y);
	case 7:
		double B11 = xfac * (1.0 - params->a * params->a), B12 = -yfac, B21 = xfac * params->b, B22 = 0.0;
		double aux1 = B21 - B12, aux2 = B11 + B22, aux = aux1 * aux1 + aux2 * aux2;
		aux1 *= aux2 / aux;
		aux2 *= aux2 / aux;
		double Q11 = -(B11 * aux2 + B21 * aux1);
		double Q12 = -(B12 * aux2 + B22 * aux1);
		double Q21 = Q12;
		double Q22 = -(B22 * aux2 - B12 * aux1);
		return Q11 * z.x * z.x + (Q12 + Q21) * z.x * z.y + Q22 * z.y * z.y;
	}
}

/*
The first case of the slowness here is for all quasipotential cases.
The slowness is given by
	s(z,v) = |b(z)| - <b(z),v/|v|>
This has derivatives
	Ds/Dz = [b(z)/|b(z)| - v/|v| ] Db(z)
	Ds/Dv = -b(z)/|v| + v <b(z),v>/|v|^3

*/
double masterSlow(const PairDoub& z, const PairDoub& v, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0: case 3: case 4: case 5: case 6: case 7: case 8:
		return masterDrift(z, sp).norm() - dot(masterDrift(z, sp), v) / v.norm();
	case 1:
		return 1;
	case 2:
		return sqrt((z.x + sin(z.x+z.y))* (z.x + sin(z.x + z.y)) + sin(z.x+z.y)* sin(z.x+z.y));
	}
}
PairDoub masterGradZ(const PairDoub& z, const PairDoub& v, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0: case 3: case 4: case 5: case 6: case 7: case 8:
		return leftMultiply(masterDrift(z, sp) / masterDrift(z, sp).norm() - v / v.norm(), masterGradDrift(z, sp));
	case 1:
		return PairDoub(0, 0);
	case 2: {
		double th = z.x + z.y;
		return PairDoub( (z.x+sin(th))*(1+cos(th)) + sin(th)*cos(th),(z.x+sin(th))*cos(th) + sin(th)*cos(th)) / masterSlow(z,v,sp);
	}
	}
}
MatrixDoub masterHessZ(const PairDoub& z, const PairDoub& v, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	double temp = 0;
	switch (params->switchKey) {
	default:
	case 0: case 3: case 4: case 5: case 6: case 7: case 8:
		return MatrixDoub(PairDoub(INFINITY, INFINITY), PairDoub(INFINITY, INFINITY));
	case 1:
		return MatrixDoub(PairDoub(0, 0), PairDoub(0, 0));
	case 2: {
		double th = z.x + z.y;
		PairDoub grad = masterGradZ(z, v, sp);
		double A11 =((1 + cos(th))*(1+cos(th)) - (z.x + sin(th)) * sin(th) + cos(th) * cos(th) - sin(th) * sin(th) - grad.x*grad.x) / masterSlow(z, v, sp);
		double A12 = (cos(th) * (1 + cos(th) ) - (z.x + sin(th)) * sin(th) + cos(th) * cos(th) - sin(th) * sin(th) - grad.x*grad.y) / masterSlow(z, v, sp);
		double A22 = (cos(th) * cos(th) - sin(th) * (z.x + sin(th)) + cos(th) * cos(th) - sin(th) * sin(th) - grad.y*grad.y) / masterSlow(z, v, sp);
		return MatrixDoub(PairDoub(A11,A12),PairDoub(A12,A22));
	}
	}
}


PairDoub masterGradV(const PairDoub& z, const PairDoub& v, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0: case 3: case 4: case 5: case 6: case 7: case 8:
		return  v * dot(masterDrift(z, sp), v) / pow(v.norm(), 3)- masterDrift(z, sp) / v.norm();
	case 1: case 2: 
		return PairDoub(0, 0);
	}
}

// do I use this?
MatrixDoub masterHessV(const PairDoub& z, const PairDoub& v, const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0: 
		return MatrixDoub(PairDoub(INFINITY, INFINITY), PairDoub(INFINITY, INFINITY));
	case 1:
		return MatrixDoub(PairDoub(0, 0), PairDoub(0, 0));
	}
}

std::string masterStr(const SpeedParams& spdPrms) {
	std::string str;
	switch (spdPrms.switchKey) {
	default:
	case 0:
		return str;
	case 1:
		str = "1";
		return str;
	case 2:
		str = "sqrt(sin^2(x + y) + (x + sin(x + y) )^2)";
		return str;
	}

}

MatrixDoub masterDriftLinearization(const void* sp) {
	const SpeedParams* params = (SpeedParams*)sp;
	switch (params->switchKey) {
	default:
	case 0:
		return MatrixDoub(PairDoub(-params->a, -params->b), PairDoub(params->a * params->b, -1)); 
	case 1:	case 2:
		return MatrixDoub(PairDoub(INFINITY, INFINITY), PairDoub(INFINITY, INFINITY));
	case 3:
		MatrixDoub(PairDoub(0, 0), PairDoub(0, 0));
	case 4:
		return MatrixDoub(PairDoub(-(1 + params->b), -1), PairDoub(1 + params->b, -1)); 
	case 5:
		return MatrixDoub(PairDoub(-2, -params->b), PairDoub(2* params->b, -1));
	case 6:
		return MatrixDoub(PairDoub(xfac*(1 - params->a * params->a), -yfac), PairDoub(xfac*params->b, 0));
	case 7:
		return MatrixDoub(PairDoub(xfac*(1 - params->a * params->a), -yfac), PairDoub( xfac*params->b, 0));
	case 8:
		return MatrixDoub(PairDoub(-2, 0), PairDoub(-2, 0));
	}
}

