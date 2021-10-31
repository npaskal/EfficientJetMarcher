#pragma once
/*
	\file:		utility.h
	\brief:		This file declares run and diagnostic functions for main.
	\author:	Nick Paskal
	\date:		10/27/2021
*/

#include "fastmarch.h"
#include "fastMarchEJM.h"

// run Fast March object for single mesh and set of settings 
void runSingle(
	const MeshInfo& mesh,
	UpdaterBase& upIn,
	const RunOptions& opIn,
	const StencilData& stIn,
	const SpeedInfo& spIn,
	int olimK,
	char algKey
);

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
	char quadType
);

// run accuracy diagnostic to determine error convergence rate for several 
// sets of settings
void fullConvRate(
	UpdaterBase& upIn, 
	const RunOptions& opIn, 
	const StencilData& stIn, 
	int olimK
);

// run OLIM for a variety of different k values
void olimKTest(
	const MeshInfo& mesh, 
	UpdaterBase& upIn, 
	const SpeedInfo& spIn, 
	RunOptions& opIn
);




//void checkSpeedDerivsCorrectness(const DriftParams* sp, bool quasi);