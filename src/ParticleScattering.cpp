#include "ParticleScattering.h"

#include <string>
#include <sys/stat.h>

#include <TMath.h>

/* Started to implement class that can
  follow geant trajectory  */

ParticleScattering::ParticleScattering(CherenkovRadiationModels* ch_) :
 ch(ch_) {
  angleGenerator = new TRandom3;
}

ParticleScattering::~ParticleScattering() {
  delete angleGenerator;
}

void ParticleScattering::SetSeed(unsigned int seed_) {
  angleGenerator->SetSeed(seed_);
}

unsigned int ParticleScattering::GetSeed() const {
  return angleGenerator->GetSeed();
}

double ParticleScattering::calculateQ() {

  double Z = 10;
  double c = 3e8; // Speed of light
  double re = 2.8e-15; //m
  double m_e = 0.511; // MeV/c^2
  double nElectrons = 10; // Electrons per 1 molecule of water
  double n_A = 6.02e23;  // Avogadro's number
  double m_w = 0.02;  //kg/mol heavy water
  double eps_0 = 8.85e-12;
  double charge = 1.6e-19;
  // In SI system kg, m^3 etc.

  double ne = 1000/m_w * n_A * nElectrons; // N/m^3
  double wp = TMath::Sqrt( charge * charge * ne / 9.1e-31 / eps_0  );

  double energyDependentPart = TMath::Power(ch->GetEnergy() / m_e , 2) * 
    ch->GetBeta() * TMath::Log(137 * TMath::Power(Z, -1./3)) / 
    TMath::Power( TMath::Power(1 + ch->GetEnergy()/ m_e , 2) - 1, 2);
    double result = Z * wp * wp * re * energyDependentPart / c;
    // std::cout << "Q: " << result << std::endl;
  return qMultiplier * result;
}

double ParticleScattering::GetAngleSquaredMean() {
   return calculateQ() / 2 * GetStepLength() / ch->GetVelocity();
}

double ParticleScattering::GetDeltaAngle() {
  return angleGenerator->Gaus(GetMeanAngle(), TMath::Sqrt(GetAngleSquaredMean()));
}

void ParticleScattering::ParseTrajectory(std::string filename, int numberOfElectrons) {
  std::ifstream stream;
  stream.open(filename);

  int counter = 0;

  electronTrajectory electron;
  double itrack = 0, x = 0, y = 0, z = 0, p = 0;
  std::string firstRaw;
  getline(stream, firstRaw);
  double itrackOld = 0;

  while (stream >> itrack && stream >> x && stream >> y
   && stream >> z && stream >> p){
    if (itrackOld != itrack) {
      allElectrons.push_back(electron);
      electron.clear();
      itrackOld = itrack;

      ++counter;
      if (counter >= numberOfElectrons) break;
    }

    geantStep stepInfo;
    TVector3 tmp(0, 0, 0);
    // From mm to cm
    tmp.SetX(x / 10);
    tmp.SetY(y / 10);  
    tmp.SetZ(z / 10);
    stepInfo.coordinates = tmp;
    stepInfo.p = p; //transfromPtoE(p);
    electron.push_back(stepInfo);
  }

  allElectrons.push_back(electron);
  stream.close();
}

TVector3 ParticleScattering::GenerateStep(TVector3 begPoint, TVector3 endPoint) {
  TVector3 deltaVector = endPoint - begPoint;
  if (deltaVector.Mag() < GetStepLength()) return endPoint;
  double theta = deltaVector.Theta();
  double phi = deltaVector.Phi();
  double deltaTheta = angleGenerator->Gaus(GetMeanAngle(), TMath::Sqrt(GetAngleSquaredMean()));
  double rotationPhi = angleGenerator->Rndm() * 2 * TMath::Pi();

  // double deltaPhi = angleGenerator->Gaus(GetMeanAngle(), TMath::Sqrt(GetAngleSquaredMean()));
  // Should be rotated!
  // deltaVector.SetPhi(phi + deltaPhi);
  TVector3 rotationAxes = deltaVector;
  deltaVector.SetTheta(theta + deltaTheta);
  deltaVector.SetMag(GetStepLength());
  deltaVector.Rotate(rotationPhi, deltaVector);
  return begPoint + deltaVector;
}

TVector3 ParticleScattering::GenerateStepRand(TVector3 begPoint) {
  TVector3 deltaVector = begPoint;

  double theta = deltaVector.Theta();
  double deltaTheta = angleGenerator->Gaus(GetMeanAngle(), TMath::Sqrt(GetAngleSquaredMean()));
  double rotationPhi = angleGenerator->Rndm() * 2 * TMath::Pi();

  // double deltaPhi = angleGenerator->Gaus(GetMeanAngle(), TMath::Sqrt(GetAngleSquaredMean()));
  // Should be rotated!
  // deltaVector.SetPhi(phi + deltaPhi);
  deltaVector.SetTheta(theta + deltaTheta);
  deltaVector.SetMag(GetStepLength());
  deltaVector.Rotate(rotationPhi, begPoint);
  return begPoint + deltaVector;
}

// FINISH!
// void ParticleScattering::WriteTrajectoryToFile(std::string filename,
//      std::string path) {

//   std::ofstream file;
//   std::string fileNameCoordinates = (filename == "default") ? path + "/CherenkovScatteredTrack_angleSigma_" 
//   + std::to_string(GetAngleSquaredMean) + ".txt" : path + filename;
//   std::cout << "Coordinates saved to " << fileNameCoordinates << std::endl;
//   mkdir(path.data(), 0777);
//   file.open(fileNameCoordinates.data());
//   file << "X\tY\tZ\n";
// }