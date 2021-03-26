#pragma once

/* This is the main file of project that is the base for .so object 
	which you can open with ROOT after compilation stage.

	Add functions to header and implement them in .cpp file.
	Run make in CherenkovRootProject directory to compile.
	After compilation, you will be able to call these functions in ROOT running 
	root thisFileName.so

	Don't forget to change Linkdef.h file after adding new functionality.
*/

void runSimulation2D();

void runSimulation3D();

void testSmearing();

void testSmearingAndScattering();

void simulationFromGeantTrajectory(double stepLength, double qMultiplier, int imputNumber);

void simulationFromOwnTrajectory(double stepLength, double qMultiplier, int inputNumber);

void simulationSingleScattering();