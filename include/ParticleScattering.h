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

/**
  Structure for each geant4 step
  */
struct geantStep
{
  TVector3 coordinates = TVector3(0, 0, 0);
  double p = 0;
};

typedef std::vector<geantStep> electronTrajectory;

/** 
  Class that should have been smear the trajectory
  But now it parses and works with Geant4 trajectories
  */
class ParticleScattering
{
public:
  ParticleScattering(CherenkovRadiationModels *ch_);
  ~ParticleScattering();

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  void SetMeanAngle(double meanAngle_) {
    meanAngle = meanAngle_;
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double GetMeanAngle() const {
    return meanAngle;
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  void SetSeed(unsigned int seed_);

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  unsigned int GetSeed() const;

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  void SetQMultiplier(double qMultiplier_) {
    qMultiplier = qMultiplier_;
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double GetQMultiplier() const {
    return qMultiplier;
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  void SetStepLength(double stepLength_) {
    stepLength = stepLength_;
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double GetStepLength() const {
    return stepLength;
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double GetDeltaAngle();

  /**
    Parses geant4 traectory from the file
    */
  void ParseTrajectory(std::string filename, int numberOfElectrons=100);
  
  /** 
    Transforms momentum to full energy
    */
  double transfromPtoE(double p) const {
    return TMath::Sqrt(electronMass * electronMass + p * p);
  }

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double calculateQ();

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double GetAngleSquaredMean();

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  TVector3 GenerateStep(TVector3 begPoint, TVector3 endPoint);

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  TVector3 GenerateStepRand(TVector3 begPoint);

  //FINISH!
  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  void WriteTrajectoryToFile(std::string filename="default",
       std::string path="Trajectories");

  /**
    Vector of electron trajectories
    */
  std::vector<electronTrajectory> allElectrons;

private:
  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  TRandom3 *angleGenerator;
  CherenkovRadiationModels* ch;

  /** 
    Deprecated. Used for multiple scattering smearing.
    */
  double meanAngle = 0;
  double stepLength = 1e-4; // In mkm in gauss
  double qMultiplier = 1;

  const double electronMass = 0.511; //MeV

};