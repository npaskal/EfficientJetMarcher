/*
	\file:		main.cpp
	\brief:		This file defines startup function main.cpp for EJM project.
	\author:	Nick Paskal
	\date:		10/27/2021
*/
#include <iostream>
#include "fastmarch.h"
#include "fastmarchEJM.h"
#include "fastmarchOLIM.h"
#include "updaterBase.h"
#include "updaterLinear.h"
#include "updaterHermite.h"
#include "utility.h"
#include "driftFunctions.h"






int main() {
	/* Primary options description:

	Run Type Key:
			's' = single run
			'c' = run and write analysis data for varying N values 
					and fixed other settings
			'z' = run and write analysis data for varying N values 
					and varying other settings
			'k' = analyze effectiveness of different K values for OLIM

	Algorithm Key: (variable algKey)
			FM = Traditional fast march
			EJM = Efficient jet marcher
			OLIM = Ordered line integral method

	Speed Key:		(with parameters a & b)
			0: 			drift = [ - a*x -b*y, 
								    a*b*x - y  ]
			3: 			drift = [  -a*x^3 - b*y^3, 
									a*b*x^3 - y^3 ]
			5: 			NONLINEAR DIAGNOSTIC DRIFT
							drift = [ -2x - 0.75*a*x^2 - by, 
									  2*b*x + 0.75*a*b*x^2 - y ]
			6:  		FITZHUGH-NAGUMO
							drift = [ - x *(a^2-1 - a*x + 1/3*x^2) - y, 
										b *x ]
			8:  		MAIER-STEIN
							drift = [-2x - x^3 - b*x*y^2 + 3x^2 + by^2, 
									 -2y - x^2*y + 2x*y ]

	Quadrature Key:
			'e' = LINEAR MAP -- endpoint quadrature
			'm' = LINEAR MAP -- midpoint rule
			'h' = CUBIC MAP -- Hermite interpolation with Simpsons' rule quad
	*/


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// PRIMARY SETTINGS: -- user should modify these
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	#include <fstream>
	
	// run type (see description above)
	char runType = 's';

	// algorithym (see descrition above)
	int algKey = EJM;

	// set to true to use default settings for algorithm (recommended)
	bool useDefaultSettings = true; 

	int n =257; // number of mesh points

	// uniform rectangular mesh constructor:
	//		(nx, xleft, xright, ny, yleft, yright)
	const MeshInfo mesh(n, -1.0, 1.0, n, -1.0, 1.0);

	// quadrature/updater key (see description above)
	char quadKey = 'h'; 

	
	// drift information -- see masterDrift function in driftFunctions.cpp

	// switch key for drift functions in masterDrift -- see description above
	int speedKey = 8;

	double 
		a = 2.0,			// parameter 1 of drift function
		b = 1.0;			// parameter 2 of drift function

	// set to false to suppress console output
	bool verbose = true; 

	// set to true to store detailed debug info about each Accepted point
	// Note: memory requirements may be large when set to true
	bool debugMode = true; 



	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	// SECONDARY SETTINGS: -- may not need to modify 
	/////////////////////  -- default values already set to optimal  
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////


	bool twoPtUpdateMinDist = true;
	bool failSafe = true;
	bool fakeFilter = true; 

	if (quadKey != 'h' && useDefaultSettings) {
		twoPtUpdateMinDist = false;
		failSafe = false;
		fakeFilter = false;
	}
		
	int initType = 0;

	RunOptions opt{ verbose, debugMode, twoPtUpdateMinDist,
		failSafe, fakeFilter, initType };

	std::unique_ptr<UpdaterBase> updater(nullptr);
	switch (quadKey) {
	case 'e':
		updater.reset(new UpdaterLinearEndpoint());
		break;
	case 'm':
		updater.reset(new UpdaterLinearMidpoint());
		break;
	case 'h':
		updater.reset(new UpdaterHermite());
		break;
	}


	const DriftParams spdPrms(speedKey, a, b);
	const SpeedInfo speeds(spdPrms);	


	//********** EJM ADDITIONAL SETTINGS **********
	double alpha = 0.999;
	int maxStencilSize = 20;
	int numStencilThetaBins = 1000;
	bool boxUpdateASR = true;
	int stencilType = ACUTE;
	const StencilData stencilData(alpha, maxStencilSize, 
		numStencilThetaBins, boxUpdateASR, stencilType);

	   	 
	//********** OLIM ADDITIONAL SETTINGS **********
	int olimK =26;	

	switch (runType) {
	case 's': {
		runSingle(mesh, *updater, opt, stencilData, speeds, olimK, algKey);
		break; }
	case 'c': {
		convergenceRate(*updater, opt, stencilData, speeds, olimK, algKey,
			true, quadKey);
		break; }
	case 'z': {
		fullConvRate(*updater, opt, stencilData, olimK);
		break;
	}
	case 'k': {
		opt.verbose = false;
		olimKTest(mesh, *updater, speeds, opt);
		break;
	}
	}
} 

