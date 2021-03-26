#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>

#include <TRandom3.h>
#include <TMath.h>
#include <TVector3.h>

#include "CherenkovRadiationModels.h"

struct geantStep
{
  TVector3 coordinates = TVector3(0, 0, 0);
  double p = 0;
};

typedef std::vector<geantStep> electronTrajectory;

/**
  Now only for Theta angle distribution
  Maybe better to unify Geant4 and this class output
  This class should calculte trajectory
  And wrtie it in file in the same style
  As Geant4 does
**/

class ParticleScattering
{
public:
  ParticleScattering(CherenkovRadiationModels *ch_);
  ~ParticleScattering();

  void SetMeanAngle(double meanAngle_) {
    meanAngle = meanAngle_;
  }

  double GetMeanAngle() const {
    return meanAngle;
  }

  void SetSeed(unsigned int seed_);

  unsigned int GetSeed() const;

  void SetQMultiplier(double qMultiplier_) {
    qMultiplier = qMultiplier_;
  }

  double GetQMultiplier() const {
    return qMultiplier;
  }

  void SetStepLength(double stepLength_) {
    stepLength = stepLength_;
  }

  double GetStepLength() const {
    return stepLength;
  }

  double GetDeltaAngle();

  void ParseTrajectory(std::string filename, int numberOfElectrons=100);
  
  double transfromPtoE(double p) const {
    return TMath::Sqrt(electronMass * electronMass + p * p);
  }

  double calculateQ();

  double GetAngleSquaredMean();

  TVector3 GenerateStep(TVector3 begPoint, TVector3 endPoint);

  TVector3 GenerateStepRand(TVector3 begPoint);

  //FINISH!
  void WriteTrajectoryToFile(std::string filename="default",
       std::string path="Trajectories");

  std::vector<electronTrajectory> allElectrons;

private:

  TRandom3 *angleGenerator;
  CherenkovRadiationModels* ch;

  double meanAngle = 0;
  double stepLength = 1e-4; // In mkm in gauss
  double qMultiplier = 1;

  const double electronMass = 0.511; //MeV

};