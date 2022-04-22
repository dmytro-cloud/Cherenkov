#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

// includes all header files
// #include "CherenkovRadiationModels.h"
// #include <TComlex.h>

// All classes
#pragma link C++ class CherenkovRadiationModels+;
#pragma link C++ class ParticleScattering+;

// all functions
#pragma link C++ function runSimulation2D;
#pragma link C++ function runSimulation3D;
#pragma link C++ function testSmearing;
#pragma link C++ function testSmearingAndScattering;
#pragma link C++ function simulationFromGeantTrajectory;
#pragma link C++ function simulationFromOwnTrajectory;
#pragma link C++ function simulationSingleScattering;
#pragma link C++ function calculate_light_output;
#pragma link C++ function runSimulation2D_default;

#endif