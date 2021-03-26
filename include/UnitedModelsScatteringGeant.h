#pragma once

/*
    This is the main file of project that is the base for .so object 
    which you can open with ROOT after compilation stage.

    Add functions to header and implement them in .cpp file.
    Run make in CherenkovRootProject directory to compile.
    After compilation, you will be able to call next functions in ROOT running 
    root thisFileName.so

    Don't forget to change Linkdef.h file after adding new functionality.
*/

/**
  Simple 2D simulation with one step
    */
void runSimulation2D();

/**
  Deprecated and removed from .cpp
  */
void runSimulation3D();

/**
  Deprecated and removed from .cpp
  Used q-parameter for smearing with mulptiple scattering theory
  */
void testSmearing();

/**
  Deprecated and removed from .cpp
  */
void testSmearingAndScattering();

/**
  Deprecated and removed from .cpp
  */
void simulationFromGeantTrajectory(double stepLength, double qMultiplier, int imputNumber);

/**
  Deprecated and removed from .cpp
  */
void simulationFromOwnTrajectory(double stepLength, double qMultiplier, int inputNumber);

/**
  Function to run simulation using Geant4 trajectory.
  Name of the input and output file as the number of bins and tracks are hard-coded.
  */
void simulationSingleScattering();