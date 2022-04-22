#pragma once

#include <iostream>
#include <string>
#include <fstream>
// #include <sys/stat.h>

#include <TMath.h>
#include <TComplex.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TVector3.h>

class CherenkovRadiationModels {
public:
  /**
    Constructor
    */
  CherenkovRadiationModels() {};
  ~CherenkovRadiationModels() {};

  /**
    Get momentum in MeV/c
    */
  double GetMomentum() const {
    return momentum;
  }

  /**
    Set momentum in MeV/c
    */
  void SetMomentum(double momentum_) {
    momentum = momentum_;
  }

  /**
   Set mass in MeV
   */
  void SetMass(double mass_) {
    mass = mass_;
  }

  /**
   Get mass in MeV
   */
  double GetMass() const {
    return mass;
  }

  // Set WaveLength in cm
  void SetWaveLenght(double wavelength_) {
    lambda = wavelength_;
  }

  /**
   Get WaveLength in cm
   */
  double GetWaveLength() const {
    return lambda;
  }

  /**
   Set refractive index
   */
  void SetRefractiveIndex(double n_) {
    n = n_;
  }

  /**
   Get refractive index
   */
  double GetRefractiveIndex() const {
    return n;
  }

  /** 
    Set sphere radius 
    */
  void SetSphereR(double r_) {
    sphereR = r_;
  }

  /** 
    Get Sphere radius
    */
  double GetSphereR() const {
    return sphereR;
  }

  /**
   Set number of bins in hist (number of cells that circle is divided)
   */
  void SetNBins(double nBins_) {
    nBins = nBins_;
  }

  /**
   Get number of bins in hist
   */
  double GetNBins() const {
    return nBins;
  }

  /**
    Get frequency omega
    */
  double GetOmega() const;

  /** 
    Get wave vector k
    */
  double GetWaveVector() const;

  /**
    Get relativistic ratio beta = v/c
    */
  double GetBeta() const;

  /** 
    Get particle velocity that calculates from energy
    */
  double GetVelocity() const;

  /**
    Get Cherenkov cosinus as c / (v * n)
    */
  double GetCherenkovCos() const;

  /**
    Get amplitude for 2D scattering of partice
    */
  double GetAmplitude2D() const {
    return amplitude2D;
  }

  /**
    Set amplitude for 2D scattering of partice
    */
  void SetAmplitude2D(double amplitude2D_) {
    amplitude2D = amplitude2D_;
  }

  /** 
    Get frequency of interactions of electron with matter for L-distance 
    */
  double GetOmegaY() const {
    return omegaY;
  }

  /** 
    Set frequency of interactions of electron with matter for L-distance 
    */
  void SetOmegaY(double omegaY_){
    omegaY = omegaY_;
  }

  /** 
    Get light speed in cm/s
    */
  static double GetLightSpeed() {
    return c; // cm/s
  }

  /**
    Get electron charge in Fr
    */
  static double GetElectronCharge() {
    return electronCharge; // Fr
  }

  /** 
    Calculates distance to te point (L) and theta in the particle's system.
    Default: 1D with phi=0
    */
  std::pair<double, double> AngleTransform(double angle, TVector3 initialPosition,
     TVector3 nextPosition, double phi=0.) const;

   /** 
    As in Dedrick
    */
  std::pair<double, double> AngleTransform_2(double angle, TVector3 initialPosition,
     TVector3 nextPosition, double phi=0.) const;

  /**
    Actually, almost the same that Geant4 does but taking into account the finite length of step
    */
  double CoherentMyModel(double cos, TVector3 initialPosition, TVector3 nextPosition,
   double phi=0., bool returnSquared=true, bool SI=false);

  /**
    CURRENTLY IN USE
    Modified Dedrick's analitical formula to calcualte differential density of power
    Irradiated by charge particle that moves in the volume. Considers the waves itself.
    Returns complex value for the wave created during the step.
    */
  TComplex CoherentDedricksModel(double cos, TVector3 initialPosition, TVector3 nextPosition, double currentTime,
   double phi=0., bool returnSquared = false);

  /** 
    Formula with sinus aproximation of particle movement for scattering.
    First approximation
    */
  double CoherentSinusoidalModel(double cos, TVector3 initialPosition, TVector3 nextPosition,
   double phi=0., bool returnSquared = true);

  /**
   Save the currently drawn histogram, in png and .C formats. The plot
   Will be placed in the appropriate directory.
   */
  void PrintDrawnHistogram(TCanvas* gCanvas, const char* basename);

  /**
    It parses Geant4 output
    TODO: IMPROVE PARSE AND GEANT4 OUTPUT
    FIX: already implemented in ParticleScattering.h
  */
  std::vector<TVector3> ParseGeantFile(std::ifstream& stream);

private:
  static constexpr double c = 29979245800; // cm/s
  static constexpr double electronCharge = 1.602e-19*c/10; // Fr (statC) 

    /*  Unit system: Gauss */
  // double energy = 5; // MeV
  double momentum = 5; // MeV / c
  double mass = 0.511; // Mass of particle in MeV

  // Fix omega for now
  double lambda = 450 * 1e-7; // Cherenkov radiation wavelength (in cm)

  // Fix n for now (depends on omega)
  double n = 1.33;

  // Deprecated. Some variables for sin trajectory idea
  /* Amplitude for 2D scattering of partice */
  double amplitude2D = GetWaveLength() / 10;
  /* Frequency of interactions of electron with matter for L-distance */
  double omegaY = GetOmega() / 20;

  // Observed position on a sphere at 8.5 m?
  double sphereR = 850; // сm
  int nBins = 10000;
};

