#include "CherenkovRadiationModels.h"
#include "UnitedModelsScatteringGeant.h"
#include "ParticleScattering.h"

#include <TH2F.h>
#include <TRandom3.h>
#include <TComplex.h>
#include <TFile.h>

#include <ctime>

double calculateEnergyLoses(CherenkovRadiationModels *ch) {

double c = 3e8; // Speed of light
double re = 2.8e-15; //m
double m_e = 9.1e-31; // MeV/c^2
double nElectrons = 10; // Electrons per 1 molecule of water
double n_A = 6.02e23;  // Avogadro's number
double m_w = 0.018;  //kg/mol
double eps_0 = 8.85e-12;
double I = 75e-6; // MeV exictation value for water
double charge = 1.6e-19;

// In SI system kg, m^3 etc.
double ne = 1000/m_w * n_A * nElectrons; // N/m^3

  double logFactor = TMath::Log( 2 * ch->GetMass() * ch->GetBeta() * ch->GetBeta() / 
    (I * (1 - ch->GetBeta() * ch->GetBeta()) ) ) - ch->GetBeta() * ch->GetBeta();
  double coef1 = 4 * TMath::Pi() / (m_e * c * c) * ne / TMath::Power(ch->GetBeta(), 2);
  double coef2 = TMath::Power(charge * charge / (4 * TMath::Pi() * eps_0), 2);

  double result = coef1 * coef2 * logFactor * 6.242e+18 / 100 / 1e6; // positive dE/dx in MeV / cm
  return result;
}


// void runSimulation2D() {
//   CherenkovRadiationModels cherenkov;
//   //Show beta
//   std::cout << "Beta: " << cherenkov.GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at cos(theta) = " << cherenkov.GetCherenkovCos() << std::endl;
//   double binsn;
//   std::cin >> binsn;
//   cherenkov.SetNBins(binsn);
//   double currentTime = 0;

//   // File to write coordinates of particle
//   std::ifstream file;
//   std::string fileNameCoordinates = "geant4coordinates/" + std::to_string(int(cherenkov.GetEnergy())) + "MeVElectron.txt";
//   bool test = true;
//   if ( test ) {
//     std::string testFileName = "geant4coordinates/" + std::to_string(0) + "MeVElectronTest.txt";
//     fileNameCoordinates = testFileName;
//   }

//   std::cout << "File name: " << fileNameCoordinates << std::endl;
//   file.open(fileNameCoordinates.data());
//   std::vector<TVector3> allCoordinates = cherenkov.ParseGeantFile(file);
//   int nsteps = allCoordinates.size(); // In fact number of steps is n - 1
//   std::cout << "Steps: " << nsteps - 1 << std::endl;
//   // Set initial position
//   TVector3 particlePosition = allCoordinates.at(0);
  
//   TH1F* hPowerRadiated2 = new TH1F("hPowerRadiated2", "Power Radiated RE; cos(#theta); Power (normalized)" , cherenkov.GetNBins(), 0, 3.14);
//   TH1F* hPowerRadiated = new TH1F("hPowerRadiated", "Power Radiated IM; cos(#theta); Power (normalized)" , cherenkov.GetNBins(), 0, 3.14);

//   // double partitionOfLambda, partitionOfOmega;
//   // std::cout << "Input a and omegaY" << std::endl;
//   // std::cin >> partitionOfLambda >> partitionOfOmega;
//   // cherenkov.SetAmplitude2D(cherenkov.GetWaveLength() * partitionOfLambda);
//   // cherenkov.SetOmegaY(cherenkov.GetOmega() * partitionOfOmega);

//  // Now it is moving with a stepLength and without energy loses
//  // light emission is considered to be from the begining point of each line
//   for (int j = 1; j < nsteps; j++) {
//     TVector3 nextParticlePosition = allCoordinates.at(j);
    
//     // Loop through angles of observation
//     for (int i = 0; i < cherenkov.GetNBins(); i++) {
//       double theta = hPowerRadiated->GetXaxis()->GetBinCenter(i+1);

//       TComplex resultPower = cherenkov.CoherentKarinsModel(theta, particlePosition, nextParticlePosition, currentTime, 0., false);
//       hPowerRadiated->Fill(theta, resultPower.Im());
//       // if (cos > 0.755 && cos < 0.7552 && cos > 2) {
//       //   std::cout << "Cos++: " << cos + 2./1000000 << " PowerIm: " << cherenkov.CoherentKarinsModel(cos + 2./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
//       //   std::cout << "Cos+: " << cos + 1./1000000 << " PowerIm: " << cherenkov.CoherentKarinsModel(cos + 1./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
//       //   std::cout << "Cos: " << cos << " PowerIm: " << resultPower.Im() << "\n";
//       //   std::cout << "Cos-: " << cos - 1./1000000 << " PowerIm: " << cherenkov.CoherentKarinsModel(cos - 1./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
//       //   std::cout << "Cos--: " << cos - 2./1000000 << " PowerIm: " << cherenkov.CoherentKarinsModel(cos - 2./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
//       //   std::cout << "###" << std::endl;
//       // }
//       // double resultPower2 = cherenkov.CoherentMyModel(cos, particlePosition, nextParticlePosition);
//       hPowerRadiated2->Fill(theta, resultPower.Re());  

//     } // End of loop through angles

//     particlePosition = nextParticlePosition;
//   }

//   TCanvas* gCanvas = nullptr;
//   // Make plots square
//   Double_t w = 1000;
//   Double_t h = 800;
//   gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   // gCanvas->SetLogy();
//   gStyle->SetOptTitle(kFALSE);
//   gStyle->SetOptStat(0);
//   gStyle->SetPalette(kDeepSea);

//   // Double_t scale = 1. / hPowerRadiated->Integral();
//   // hPowerRadiated->Scale(scale);

//   // Scaling. Maybe there is a problem with error bars after scaling.
//   // scale = 1. / hPowerRadiated2->Integral();
//   // hPowerRadiated2->Scale(scale);
//   std::cout << "Integral: " << hPowerRadiated2->Integral();

//   TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
//    "Power Radiated CoherentMyModel; #theta; dP/d#omega" , cherenkov.GetNBins(), 0, 3.14);// TMath::Pi());

//     // Loop through angles of observation and square the result
//   for (int i = 0; i < cherenkov.GetNBins(); i++) {
//     double theta = hPowerRadiatedSquared->GetXaxis()->GetBinCenter(i+1);  

//     double im = hPowerRadiated->GetBinContent(i+1);
//     double re = hPowerRadiated2->GetBinContent(i+1);
//     double resultIntegral = re * re + im * im;
//     hPowerRadiatedSquared->Fill(theta, resultIntegral);      
//   } // End of loop through angles

//   // hPowerRadiated->SetMarkerStyle(kFullTriangleDown);
//   // // hPowerRadiated->Draw("PLC PMC");
//   // hPowerRadiated->SetLineColor(kRed);
//   // hPowerRadiated->SetMarkerStyle(kFullCircle);
//   // hPowerRadiated->Draw("HIST"); 
//   // hPowerRadiated2->Draw("SAMEHIST"); 
//   hPowerRadiatedSquared->Draw("HIST");

//   gPad->BuildLegend();

//   std::string fileName = "CherenkovScatteredTrackGeantSinusoidal_" + std::to_string(cherenkov.GetEnergy())
//    + "MeV_" + std::to_string(cherenkov.GetNBins()) + "bins_" + std::to_string(nsteps - 1) + "Steps";
//   cherenkov.PrintDrawnHistogram(gCanvas, fileName.data());
// }

// void runSimulation3D() {

//   CherenkovRadiationModels cherenkov;
//   double numberOfBins;
//   std::cout << "Input number of bins: " << std::endl;
//   std::cin >> numberOfBins;
//   cherenkov.SetNBins(numberOfBins);

//   //Show beta
//   std::cout << "Beta: " << cherenkov.GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at cos(theta) = " << cherenkov.GetCherenkovCos() << std::endl;

//   // File to write coordinates of particle
//   std::ifstream file;
//   std::string fileNameCoordinates = "geant4coordinates/" + std::to_string(0/*int(cherenkov.GetEnergy())*/) + "MeVElectron.txt";
//   std::cout << "File name: " << fileNameCoordinates << std::endl;
//   file.open(fileNameCoordinates.data());
//   std::vector<TVector3> allCoordinates = cherenkov.ParseGeantFile(file);
//   int nsteps = allCoordinates.size(); // In fact number of steps is n - 1
//   std::cout << "Steps: " << nsteps - 1 << std::endl;
//   // Set initial position
//   TVector3 particlePosition = allCoordinates.at(0);
  
//   TH2F* hPowerRadiated2 = new TH2F("hPowerRadiated2", "Power Radiated CoherentMyModel; #theta; #phi; Power (normalized)",
//    cherenkov.GetNBins(), 0, TMath::Pi(), cherenkov.GetNBins(), 0, TMath::Pi() * 2);
//   TH2F* hPowerRadiated = new TH2F("hPowerRadiated", "Power Radiated CoherentSinusoidalModel; #theta; #phi; Power (normalized)",
//    cherenkov.GetNBins(), 0, TMath::Pi(), cherenkov.GetNBins(), 0, TMath::Pi() * 2);

//  // Now it is moving with a stepLength and without energy loses
//  // light emission is considered to be from the begining point of each line
//   for (int j = 1; j < nsteps; j++) {
//     TVector3 nextParticlePosition = allCoordinates[j];

//     // Loop through angles of observation
//     for (int i = 0; i < cherenkov.GetNBins(); i++) {
//       for (int j = 0; j < cherenkov.GetNBins(); j++) {

//         double theta = hPowerRadiated->GetXaxis()->GetBinCenter(i+1);
//         double phi = hPowerRadiated->GetYaxis()->GetBinCenter(j+1);

//         double resultPower = cherenkov.CoherentSinusoidalModel(TMath::Cos(theta), particlePosition, nextParticlePosition, phi);
//         hPowerRadiated->Fill(theta, phi, resultPower);     

//         double resultPower2 = cherenkov.CoherentMyModel(TMath::Cos(theta), particlePosition, nextParticlePosition, phi);
//         hPowerRadiated2->Fill(theta, phi, resultPower2); 

//       } 
//     } // End of loop through angles

//     particlePosition = nextParticlePosition;
//   }

//   TCanvas* gCanvas = nullptr;
//   // Make plots square
//   Double_t w = 1000;
//   Double_t h = 800;
//   gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   gCanvas->SetLogz();
//   gStyle->SetOptTitle(kFALSE);
//   gStyle->SetOptStat(0);
//   gStyle->SetPalette(kDeepSea);

//   Double_t scale = 5. / hPowerRadiated->Integral();
//   hPowerRadiated->Scale(scale);

//   // Scaling. Maybe there is a problem with error bars after scaling.
//   // scale = 1. / hPowerRadiated2->Integral();
//   // hPowerRadiated2->Scale(scale);

//   // hPowerRadiated->SetMarkerStyle(kFullTriangleDown);
//   // hPowerRadiated->Draw("COLZ");  
//   hPowerRadiated2->SetLineColor(kRed);
//   hPowerRadiated2->SetMarkerStyle(kFullCircle);
//   hPowerRadiated2->Draw("COLZ"); 

//   gPad->BuildLegend();

//   std::string fileName = "CherenkovScattered3D" + std::to_string(cherenkov.GetEnergy())
//    + "MeV_" + std::to_string(cherenkov.GetNBins()) + "bins_" + std::to_string(nsteps - 1) + "Steps";
//   cherenkov.PrintDrawnHistogram(gCanvas, fileName.data());

// }

// void testSmearing() {

//   CherenkovRadiationModels cherenkov;

//   //Show beta
//   std::cout << "Beta: " << cherenkov.GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at cos(theta) = " << cherenkov.GetCherenkovCos() << std::endl;
//   double binsn, stepLength;
//   std::cout << "Enter number of bins \n";
//   std::cin >> binsn;
//   cherenkov.SetNBins(binsn);

//   std::cout << "Enter Path Length (cm) \n";
//   double pathLength; std::cin >> pathLength;
//   std::cout << "Enter step length (mkm) \n";
//   std::cin  >> stepLength;
//   stepLength *= 10e-5;
//   // For 1 cm path 
//   double nsteps = pathLength / stepLength;
//   double currentTime = 0; // Time in seconds

//   std::cout << "Number of steps: " << nsteps << std::endl;

//   // Theta angle distributions
//   TRandom3 *angleGenerator = new TRandom3();
//   angleGenerator->SetSeed(1);
//   auto seed = angleGenerator->GetSeed();
//   std::cout << "Seed: " << seed << std::endl;
//   double angleSquaredMean = 0.001; // In radians
//   double meanAngle = 0;


//   // File to write coordinates of particle
//   std::ofstream file;
//   std::string fileNameCoordinates = "coordinates/Cherenkov_withoutScatter_stepLength_" + std::to_string(stepLength) + ".txt";
//   std::cout << fileNameCoordinates << std::endl;
//   // mkdir("coordinates", 0777);
//   file.open(fileNameCoordinates.data());
//   file << "X\tY\tZ\n";

//   TVector3 particlePosition(0.,0., 0.000001);
//   TVector3 nextParticlePosition(0.,0.,0.000001);
  
//   TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
//    "Integral Result sum CoherentKarinsModel Real; #theta; Integral result real part" , cherenkov.GetNBins(), 0, 3.14);// TMath::Pi());
//   TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
//    "Integral Result sum CoherentKarinsModel Im; #theta; Integral result im part" , cherenkov.GetNBins(), 0, 3.14);// TMath::Pi());
//  // It moves with a stepLength and without energy loses
//  // light emission is considered to be from the begining point of each line
//   for (int j = 1; j < nsteps; j++) {
//     double deltaAngle = 0; // angleGenerator->Gaus(meanAngle, TMath::Sqrt(angleSquaredMean));
//     file << particlePosition.X() << "\t" << particlePosition.Y() << "\t" << particlePosition.Z() << "\n";

//     // nextParticlePosition.SetMag(finalMag);

//     nextParticlePosition.SetZ(particlePosition.Z() + stepLength);

//     // Loop through angles of observation
//     for (int i = 0; i < cherenkov.GetNBins(); i++) {
//       double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);  

//       TComplex resultPower = cherenkov.CoherentKarinsModel(theta,
//        particlePosition, nextParticlePosition, currentTime, 0., false);
//       hPowerRadiatedRe->Fill(theta, resultPower.Re());        
//       hPowerRadiatedIm->Fill(theta, resultPower.Im());      
//     } // End of loop through angles

//     currentTime += stepLength / cherenkov.GetVelocity();
//     particlePosition = nextParticlePosition;
//   }

//   TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
//    "Power Radiated CoherentMyModel; #theta; Power (normalized)" , cherenkov.GetNBins(), 0, 3.14);// TMath::Pi());

//     // Loop through angles of observation and square the result
//   for (int i = 0; i < cherenkov.GetNBins(); i++) {
//     double theta = hPowerRadiatedSquared->GetXaxis()->GetBinCenter(i+1);  

//     double im = hPowerRadiatedIm->GetBinContent(i+1);
//     double re = hPowerRadiatedRe->GetBinContent(i+1);
//     double resultIntegral = re * re + im * im;
//     hPowerRadiatedSquared->Fill(theta, resultIntegral);      
//   } // End of loop through angles

//   TCanvas* gCanvas = nullptr;
//   // Make plots square
//   Double_t w = 1000;
//   Double_t h = 800;
//   gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   // gCanvas->SetLogy();
//   gStyle->SetOptTitle(kFALSE);
//   gStyle->SetOptStat(0);
//   gStyle->SetPalette(kDeepSea);

//   // Double_t scale = 1. / hPowerRadiatedSquared->Integral();
//   // hPowerRadiatedSquared->Scale(scale);


//   // hPowerRadiatedSquared->SetMarkerStyle(kFullTriangleDown);
//   // hPowerRadiatedSquared->Draw("HIST");

//   // hPowerRadiatedRe->Draw("HIST");
//   // hPowerRadiatedIm->SetLineColor(kRed);
//   // hPowerRadiatedIm->Draw("SAME");
//   hPowerRadiatedSquared->SetLineColor(kGreen);
//   hPowerRadiatedSquared->Draw("HIST");

//   gPad->BuildLegend();

//   std::string fileName = "TestStepLength_" + std::to_string(int(cherenkov.GetEnergy())) 
//   + "MeV_" + std::to_string(cherenkov.GetNBins()) + "bins_" + std::to_string(stepLength) + "Length";
//   cherenkov.PrintDrawnHistogram(gCanvas, fileName.data());
// }

// void testSmearingAndScattering() {

//   CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
//   auto scattering = ParticleScattering(cherenkov);
//   scattering.ParseTrajectory("electronTrajectories/test.txt");
//   cherenkov->SetEnergy(5.);
//   double qMultiplier; std::cin >> qMultiplier;
//   scattering.SetQMultiplier(qMultiplier);

//   //Show beta
//   std::cout << "Beta: " << cherenkov->GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at theta = " << TMath::ACos(cherenkov->GetCherenkovCos()) << std::endl;
//   double binsn = 10000 , stepLength;

//   cherenkov->SetNBins(binsn);

//   cherenkov->SetWaveLenght(450 * 1e-7);

//   std::cout << "Enter step length (mkm) \n";
//   std::cin  >> stepLength;
//   stepLength *= 10e-5;    // translation to mkm in gauss system
//   // For 1 cm path 
//   double trackLength = 5;
//   double nsteps = trackLength / stepLength;
//   double currentTime = 0; // Time in seconds


//   // Theta angle distributions
//   TRandom3 *angleGenerator = new TRandom3();
//   angleGenerator->SetSeed(1);
//   auto seed = angleGenerator->GetSeed();
//   std::cout << "Seed: " << seed << std::endl;
//   double angleSquaredMean = 0.001; // In radians
//   double meanAngle = 0;

//   // Theta angle distributions
//   std::cout << "Seed: " << seed << std::endl;

//  // In radians
//   std::cout << "Well, estimated sigma^2 = " << scattering.GetAngleSquaredMean()  << std::endl;

//   // File to write coordinates of particle
//   std::ofstream file;
//   std::string fileNameCoordinates = "coordinates/Cherenkov_Scatter_stepLength_" + std::to_string(stepLength) + ".txt";
//   std::cout << fileNameCoordinates << std::endl;
//   // mkdir("coordinates", 0777);
//   file.open(fileNameCoordinates.data());
//   file << "X\tY\tZ\n";

//   TVector3 particlePosition(0.,0., 0.00001);
//   TVector3 nextParticlePosition(0.,0.,0.00001);

//   TVector3 geantParticlePosition(0.,0., 0.00001);
//   TVector3 geantNextParticlePosition(0.,0.,0.00001);

  
//   TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
//    "Integral Result sum CoherentKarinsModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
//   TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
//    "Integral Result sum CoherentKarinsModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());
//   TH1D* hPowerRadiatedSumAsGeantPhotons = new TH1D("hPowerRadiatedSumAsGeantPhotons",
//    "Sum as in geant; #theta; Sum result" , cherenkov->GetNBins(), 0, TMath::Pi());

//  // Now it is moving with a stepLength and without energy loses
//  // light emission is considered to be from the begining point of each line
//   int geantStepsFrac = 55;
//   for (int j = 1; j < nsteps; j++) {
//     double energyLoses = calculateEnergyLoses(cherenkov) * stepLength;

//     double q = scattering.calculateQ();
//     angleSquaredMean =  scattering.calculateQ() / 2 * stepLength / cherenkov->GetVelocity();
//     double deltaAngle = angleGenerator->Gaus(meanAngle, TMath::Sqrt(angleSquaredMean));
//     double deltaAngle2 = angleGenerator->Gaus(meanAngle, TMath::Sqrt(angleSquaredMean));
    
//     TVector3 deltaPosition(stepLength * deltaAngle, stepLength * deltaAngle2, stepLength);

//     nextParticlePosition += deltaPosition;
//     geantNextParticlePosition += deltaPosition;

//     // nextParticlePosition.SetPhi(deltaAngle / 2);
//     file << nextParticlePosition.X() << "\t" << nextParticlePosition.Y() << "\t" << nextParticlePosition.Z() << "\n";

//     // Loop through angles of observation
//     for (int i = 0; i < cherenkov->GetNBins(); i++) {
//       double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);  

//       TComplex resultPower = cherenkov->CoherentKarinsModel(TMath::Cos(theta),
//        particlePosition, nextParticlePosition, currentTime, 0., false);
//       hPowerRadiatedRe->Fill(theta, resultPower.Re());      
//       hPowerRadiatedIm->Fill(theta, resultPower.Im());

//       // Case for long geant steps
//       if (j % geantStepsFrac == 0){
//         hPowerRadiatedSumAsGeantPhotons->Fill(theta, cherenkov->Geant4like(TMath::Cos(theta),
//        geantParticlePosition, geantNextParticlePosition, 0));
//       }

//     } // End of loop through angles
      
//       // Case for long geant steps
//       if (j % geantStepsFrac == 0) {
//         geantParticlePosition = geantNextParticlePosition;
//         std::cout << "Well, estimated sigma^2 = " << angleSquaredMean  << std::endl;
//         std::cout << "Loses = " << energyLoses  << std::endl;
//         std::cout << geantParticlePosition.Z() << std::endl;
//       }

//     currentTime += (nextParticlePosition - particlePosition).Mag() / cherenkov->GetVelocity();
//     particlePosition = nextParticlePosition;
//     cherenkov->SetEnergy(cherenkov->GetEnergy() - energyLoses);
//     if (cherenkov->GetBeta() * cherenkov->GetRefractiveIndex() < 1) break;

//   }

//   std::cout << "Last energy " << cherenkov->GetEnergy() << std::endl;

//   double suaredHistBins = 200;

//   TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
//    "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

//   TH1D* hPowerRadiatedSumAsGeantPhotonsFin = new TH1D("hPowerRadiatedSumAsGeantPhotons",
//    "Sum as in geant; #theta; Sum result" , suaredHistBins, 0, TMath::Pi());

//     // Loop through angles of observation and square the result
//   for (int i = 0; i < cherenkov->GetNBins(); i++) {
//     double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);  

//     double im = hPowerRadiatedIm->GetBinContent(i+1);
//     double re = hPowerRadiatedRe->GetBinContent(i+1);
//     double resultIntegral = re * re + im * im;
//     hPowerRadiatedSquared->Fill(theta, resultIntegral);
//     hPowerRadiatedSumAsGeantPhotonsFin->Fill(theta, hPowerRadiatedSumAsGeantPhotons->GetBinContent(i+1));      
//   } // End of loop through angles

//   TCanvas* gCanvas = nullptr;
//   // Make plots square
//   Double_t w = 1000;
//   Double_t h = 800;
//   gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   // gCanvas->SetLogy();
//   gStyle->SetOptTitle(kFALSE);
//   gStyle->SetOptStat(0);
//   gStyle->SetPalette(kDeepSea);

//   // Double_t scale = 1. / hPowerRadiatedSquared->Integral();
//   // hPowerRadiatedSquared->Scale(scale);

//   double scale = 1. / hPowerRadiatedSquared->Integral();
//   hPowerRadiatedSquared->Scale(scale);

//   scale = 1. / hPowerRadiatedSumAsGeantPhotonsFin->Integral();
//   hPowerRadiatedSumAsGeantPhotonsFin->Scale(scale);

//   // hPowerRadiatedSquared->SetMarkerStyle(kFullTriangleDown);
//   // hPowerRadiatedRe->Draw("HIST");
//   // hPowerRadiatedIm->SetLineColor(kRed);
//   // hPowerRadiatedIm->Draw("SAMEHIST");


//   hPowerRadiatedSumAsGeantPhotonsFin->SetLineColor(kRed);
//   hPowerRadiatedSumAsGeantPhotonsFin->Draw("HIST");
//   hPowerRadiatedSquared->SetLineColor(kGreen);
//   hPowerRadiatedSquared->Draw("HISTSAME");


//   gPad->BuildLegend();

//   std::string fileName = "TestStepLength_" + std::to_string(int(cherenkov->GetEnergy())) 
//   + "MeV_" + std::to_string(nsteps) + "steps_" + std::to_string(stepLength) + "Length";
//   cherenkov->PrintDrawnHistogram(gCanvas, fileName.data());
// }

// // 
// void testSimulationFromGeantTrajectory(double stepLength, double qMultiplier) {

//   CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
//   auto scattering = ParticleScattering(cherenkov);
//   scattering.ParseTrajectory("electronTrajectories/test.txt", 100);
//   scattering.SetQMultiplier(qMultiplier);

//   //Show beta
//   std::cout << "Beta: " << cherenkov->GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at theta = " << TMath::ACos(cherenkov->GetCherenkovCos()) << std::endl;
//   double binsn = 5000;

//   cherenkov->SetNBins(binsn);

//   cherenkov->SetWaveLenght(410 * 1e-7);

//   stepLength *= 10e-5;    // translation to mkm in gauss system
//   scattering.SetStepLength(stepLength);

//   // Theta angle distributions
//   auto seed = scattering.GetSeed();
//   std::cout << "Seed: " << seed << std::endl;

//  // In radians
//   std::cout << "Well, estimated sigma^2 = " << scattering.GetAngleSquaredMean()  << std::endl;

//   TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
//    "Integral Result sum CoherentKarinsModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
//   TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
//    "Integral Result sum CoherentKarinsModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());
//   TH1D* hPowerRadiatedSumAsGeantPhotons = new TH1D("hPowerRadiatedSumAsGeantPhotons",
//    "Sum as in geant; #theta; Sum result" , cherenkov->GetNBins(), 0, TMath::Pi());

//   double suaredHistBins = 200;

//   TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
//    "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

//   int track_counter = 0;

//   for (auto track:scattering.allElectrons){
//     track_counter++;
//     // 
//     std::cout << "Track points: " << track.size() << " number: " << track_counter << std::endl;
//     if (track.size() < 2) continue;

//     time_t start, end;
//     time(&start);

//     TRandom3 waveLengthGenerator;
//     waveLengthGenerator.SetSeed(track_counter);

//     for (size_t iGeantStep = 0; iGeantStep < track.size() - 1; ++iGeantStep) {
//       // double waveLengthGenerated = waveLengthGenerator.Gaus(450 * 1e-7, 50 * 1e-7);
//       // cherenkov->SetWaveLenght(waveLengthGenerated);
//       // std::cout << waveLengthGenerated << std::endl;

//       double currentTime = 0;
//       geantStep stepBeg = track.at(iGeantStep);
//       geantStep stepEnd = track.at(iGeantStep + 1);
//       // Threshold for Cherenkov radiation
//       if (stepBeg.energy < 0.2) break;
//       cherenkov->SetEnergy(stepBeg.energy);

//       for (int i = 0; i < cherenkov->GetNBins(); i++) {
//         double theta = hPowerRadiatedSumAsGeantPhotons->GetXaxis()->GetBinCenter(i+1);  

//         hPowerRadiatedSumAsGeantPhotons->Fill(theta, cherenkov->Geant4like(TMath::Cos(theta),
//          stepBeg.coordinates, stepEnd.coordinates, 0));
//         }

//       double geantStepLength = (stepEnd.coordinates - stepBeg.coordinates).Mag();
//       // int numSmallSteps = (geantStepLength % scattering.GetStepLength() == 0) ? geantStepLength / scattering.GetStepLength() : geantStepLength / scattering.GetStepLength() + 1;
//       // May be mistake due to step number
//       int numSmallSteps = geantStepLength / scattering.GetStepLength() + 1;
//       double averagedEnergyLosesPerStep = (stepBeg.energy - stepEnd.energy) / numSmallSteps;
//       geantStep tmpBeg = stepBeg;

//       for (int iSmallStep = 0; iSmallStep < numSmallSteps; ++iSmallStep) {
//         TVector3 nextPoint = scattering.GenerateStep(tmpBeg.coordinates, stepEnd.coordinates);
//         // Loop through angles of observation
//         for (int i = 0; i < cherenkov->GetNBins(); i++) {
//           double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);  

//           TComplex resultPower = cherenkov->CoherentKarinsModel(TMath::Cos(theta),
//            tmpBeg.coordinates, nextPoint, currentTime, 0., false);
//           hPowerRadiatedRe->Fill(theta, resultPower.Re());      
//           hPowerRadiatedIm->Fill(theta, resultPower.Im());
//         }

//         currentTime += (nextPoint - tmpBeg.coordinates).Mag() / cherenkov->GetVelocity();
//         tmpBeg.coordinates = nextPoint;
//         // double energyLoses = calculateEnergyLoses(cherenkov) * scattering.GetStepLength();
//         cherenkov->SetEnergy(cherenkov->GetEnergy() - averagedEnergyLosesPerStep);

//       }

//         // Loop through angles of observation and square the result
//       for (int i = 0; i < cherenkov->GetNBins(); i++) {
//         double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);  

//         double im = hPowerRadiatedIm->GetBinContent(i+1);
//         double re = hPowerRadiatedRe->GetBinContent(i+1);
//         double resultIntegral = re * re + im * im;
//         hPowerRadiatedSquared->Fill(theta, resultIntegral);
//       } // End of loop through angles

//     hPowerRadiatedIm->Reset();
//     hPowerRadiatedRe->Reset();

//     }

//     time(&end);
//     double seconds = difftime(end, start);
//     std::cout << "Time for one track: " << seconds << std::endl; // искомое время

//   }


//   std::cout << "Last energy " << cherenkov->GetEnergy() << std::endl;

//   TH1D* hPowerRadiatedSumAsGeantPhotonsFin = new TH1D("hPowerRadiatedSumAsGeantPhotons",
//    "Sum as in geant; #theta; Sum result" , suaredHistBins, 0, TMath::Pi());

//   // Loop through angles of observation and square the result
//   for (int i = 0; i < cherenkov->GetNBins(); i++) {
//     double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);  

//     double im = hPowerRadiatedIm->GetBinContent(i+1);
//     double re = hPowerRadiatedRe->GetBinContent(i+1);
//     double resultIntegral = re * re + im * im;
//     hPowerRadiatedSquared->Fill(theta, resultIntegral);
//     hPowerRadiatedSumAsGeantPhotonsFin->Fill(theta, hPowerRadiatedSumAsGeantPhotons->GetBinContent(i+1));      
//   } // End of loop through angles

//   TCanvas* gCanvas = nullptr;
//   // Make plots square
//   Double_t w = 1000;
//   Double_t h = 800;
//   gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   // gCanvas->SetLogy();
//   gStyle->SetOptTitle(kFALSE);
//   gStyle->SetOptStat(0);
//   gStyle->SetPalette(kDeepSea);

//   // Double_t scale = 1. / hPowerRadiatedSquared->Integral();
//   // hPowerRadiatedSquared->Scale(scale);

//   // double scale = 1. / hPowerRadiatedSquared->Integral();
//   // hPowerRadiatedSquared->Scale(scale);

//   // scale = 1. / hPowerRadiatedSumAsGeantPhotonsFin->Integral();
//   // hPowerRadiatedSumAsGeantPhotonsFin->Scale(scale);

//   // hPowerRadiatedSquared->SetMarkerStyle(kFullTriangleDown);
//   // hPowerRadiatedRe->Draw("HIST");
//   // hPowerRadiatedIm->SetLineColor(kRed);
//   // hPowerRadiatedIm->Draw("SAMEHIST");


//   // hPowerRadiatedSumAsGeantPhotonsFin->SetLineColor(kRed);
//   // hPowerRadiatedSumAsGeantPhotonsFin->Draw("HIST");
//   // hPowerRadiatedSquared->SetLineColor(kGreen);
//   // hPowerRadiatedSquared->Draw("HISTSAME");

//   gPad->BuildLegend();

//   std::string outpotFilename = "test.root";
//   TFile* outputFile = new TFile(outpotFilename.c_str(), "recreate");
//   outputFile->cd();

//   hPowerRadiatedSquared->Write();

//   outputFile->Close();

//   // std::string fileName = "TestStepLength_" + std::to_string(int(cherenkov->GetEnergy())) 
//   // + "MeV_" + std::to_string(stepLength) + "Length";
//   // cherenkov->PrintDrawnHistogram(gCanvas, fileName.data());
// }

// // 
// void simulationFromGeantTrajectory(double stepLength, double qMultiplier, int inputNumber) {

//   CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
//   auto scattering = ParticleScattering(cherenkov);
//   scattering.ParseTrajectory("electronTrajectories/test.txt", 1000);
//   scattering.SetQMultiplier(qMultiplier);

//   //Show beta
//   std::cout << "Beta: " << cherenkov->GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at theta = " << TMath::ACos(cherenkov->GetCherenkovCos()) << std::endl;
//   double binsn = 10000;

//   cherenkov->SetNBins(binsn);

//   cherenkov->SetWaveLenght(410 * 1e-7);
//   cherenkov->SetRefractiveIndex(1.33668);

//   stepLength *= 10e-5;    // translation to mkm in gauss system
//   scattering.SetStepLength(stepLength);

//   // Theta angle distributions
//   auto seed = scattering.GetSeed();
//   std::cout << "Seed: " << seed << std::endl;

//  // In radians
//   std::cout << "Estimated sigma^2 = " << scattering.GetAngleSquaredMean()  << std::endl;

//   TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
//    "Integral Result sum CoherentKarinsModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
//   TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
//    "Integral Result sum CoherentKarinsModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());

//   double suaredHistBins = 10000;

//   TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
//    "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

//   std::ofstream scattered_traj;
//   scattered_traj.open("scattered_traj" + std::to_string(qMultiplier) + ".txt", std::ios::out);
//   scattered_traj << "TrackID\tX\tY\tZ\n";

//   int track_counter = 0;

//   for (auto track:scattering.allElectrons){
//     track_counter++;
//     // 
//     std::cout << "Track points: " << track.size() << " number: " << track_counter << std::endl;
//     if (track.size() < 2) continue;

//     time_t start, end;
//     time(&start);

//     TRandom3 waveLengthGenerator;
//     waveLengthGenerator.SetSeed(track_counter);

//     for (size_t iGeantStep = 0; iGeantStep < track.size() - 1; ++iGeantStep) {
//       // double waveLengthGenerated = waveLengthGenerator.Gaus(450 * 1e-7, 50 * 1e-7);
//       // cherenkov->SetWaveLenght(waveLengthGenerated);
//       // std::cout << waveLengthGenerated << std::endl;

//       double currentTime = 0;
//       geantStep stepBeg = track.at(iGeantStep);
//       geantStep stepEnd = track.at(iGeantStep + 1);
//       // Threshold for Cherenkov radiation
//       if (stepBeg.energy < 0.2) break;
//       cherenkov->SetEnergy(stepBeg.energy);

//       double geantStepLength = (stepEnd.coordinates - stepBeg.coordinates).Mag();
//       // int numSmallSteps = (geantStepLength % scattering.GetStepLength() == 0) ? geantStepLength / scattering.GetStepLength() : geantStepLength / scattering.GetStepLength() + 1;
//       // May be mistake due to step number
//       int numSmallSteps = geantStepLength / scattering.GetStepLength() + 1;
//       double averagedEnergyLosesPerStep = (stepBeg.energy - stepEnd.energy) / numSmallSteps;
//       geantStep tmpBeg = stepBeg;

//       for (int iSmallStep = 0; iSmallStep < numSmallSteps; ++iSmallStep) {
//         TVector3 nextPoint = scattering.GenerateStep(tmpBeg.coordinates, stepEnd.coordinates);

//         scattered_traj << iGeantStep << "\t" << nextPoint.X() << "\t" << nextPoint.Y()<< "\t" << nextPoint.Z() << "\n";

//         // Loop through angles of observation
//         for (int i = 0; i < cherenkov->GetNBins(); i++) {
//           double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);  

//           TComplex resultPower = cherenkov->CoherentKarinsModel(theta,
//            tmpBeg.coordinates, nextPoint, currentTime, 0., false);
//           hPowerRadiatedRe->Fill(theta, resultPower.Re());      
//           hPowerRadiatedIm->Fill(theta, resultPower.Im());
//         }

//         currentTime += (nextPoint - tmpBeg.coordinates).Mag() / cherenkov->GetVelocity();
//         tmpBeg.coordinates = nextPoint;
//         // double energyLoses = calculateEnergyLoses(cherenkov) * scattering.GetStepLength();
//         cherenkov->SetEnergy(cherenkov->GetEnergy() - averagedEnergyLosesPerStep);

//       }

//         // Loop through angles of observation and square the result
//       for (int i = 0; i < cherenkov->GetNBins(); i++) {
//         double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);  

//         double im = hPowerRadiatedIm->GetBinContent(i+1);
//         double re = hPowerRadiatedRe->GetBinContent(i+1);
//         double resultIntegral = re * re + im * im;

//         // Correct-factor S~R^2
//         hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta) * TMath::Sin(theta));
//       } // End of loop through angles

//     hPowerRadiatedIm->Reset();
//     hPowerRadiatedRe->Reset();

//     }

//     time(&end);
//     double seconds = difftime(end, start);
//     std::cout << "Time for one track: " << seconds << std::endl; // искомое время

//   }


//   std::cout << "Last energy " << cherenkov->GetEnergy() << std::endl;
//   std::cout << "Q " << scattering.calculateQ() * stepLength / cherenkov->GetVelocity()<< std::endl;
//   std::cout << "q " << scattering.calculateQ() << std::endl;
//   std::cout << "Omega " << cherenkov->GetOmega() << std::endl;
//   std::cout << "Estimated from theory peak width: " << TMath::Power(scattering.calculateQ() /cherenkov->GetOmega(), 1./3) *
//   TMath::Power(TMath::Sqrt(1 - cherenkov->GetCherenkovCos() * cherenkov->GetCherenkovCos()) * cherenkov->GetRefractiveIndex() * cherenkov->GetBeta(), -1./3) << std::endl;

//   // TCanvas* gCanvas = nullptr;
//   // Make plots square
//   // Double_t w = 1000;
//   // Double_t h = 800;
//   // gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   // gCanvas->SetLogy();
//   // gStyle->SetOptTitle(kFALSE);
//   // gStyle->SetOptStat(0);
//   // gStyle->SetPalette(kDeepSea);

//   // Double_t scale = 1. / hPowerRadiatedSquared->Integral();
//   // hPowerRadiatedSquared->Scale(scale);

//   // hPowerRadiatedSquared->SetMarkerStyle(kFullTriangleDown);
//   // hPowerRadiatedRe->Draw("HIST");
//   // hPowerRadiatedIm->SetLineColor(kRed);
//   // hPowerRadiatedIm->Draw("SAMEHIST");

//   // hPowerRadiatedSquared->SetLineColor(kGreen);
//   // hPowerRadiatedSquared->Draw("HISTSAME");

//   // gPad->BuildLegend();

//   // std::string outpotFilename = "histograms_output/CherenkovSimulation_wCorrection_q" + std::to_string(qMultiplier) + "_part_" + std::to_string(inputNumber) + ".root";
//   std::string outputFilename = "Test_for_q_rot_" + std::to_string(qMultiplier) + ".root";
//   TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");
//   outputFile->cd();

//   hPowerRadiatedSquared->Write();
//   outputFile->Close();
//   scattered_traj.close();


//   // std::string fileName = "TestStepLength_" + std::to_string(int(cherenkov->GetEnergy())) 
//   // + "MeV_" + std::to_string(stepLength) + "Length";
//   // cherenkov->PrintDrawnHistogram(gCanvas, fileName.data());
// }


// // 
// void  simulationFromOwnTrajectory(double stepLength, double qMultiplier, int nTracks) {

//   CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
//   auto scattering = ParticleScattering(cherenkov);
//   scattering.ParseTrajectory("electronTrajectories/test.txt", 1000);
//   scattering.SetQMultiplier(qMultiplier);

//   //Show beta
//   std::cout << "Beta: " << cherenkov->GetBeta() << std::endl;

//   // The Cherenkov Radiation Description.
//   std::cout << "The expected Cherenkov angle should be at theta = " << TMath::ACos(cherenkov->GetCherenkovCos()) << std::endl;
//   double binsn = 1000;

//   cherenkov->SetNBins(binsn);

//   cherenkov->SetWaveLenght(410 * 1e-7);
//   cherenkov->SetRefractiveIndex(1.33668);

//   stepLength *= 10e-5;    // translation to mkm in gauss system
//   scattering.SetStepLength(stepLength);

//   // Theta angle distributions
//   auto seed = scattering.GetSeed();
//   std::cout << "Seed: " << seed << std::endl;

//  // In radians
//   std::cout << "Estimated sigma^2 = " << scattering.GetAngleSquaredMean()  << std::endl;

//   TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
//    "Integral Result sum CoherentKarinsModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
//   TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
//    "Integral Result sum CoherentKarinsModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());

//   double suaredHistBins = 1000;

//   TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
//    "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

//   std::ofstream scattered_traj;
//   scattered_traj.open("scattered_traj_my_" + std::to_string(qMultiplier) + ".txt", std::ios::out);
//   scattered_traj << "TrackID\tX\tY\tZ\n";

//   int track_counter = 0;


//   time_t start, end;
//   time(&start);

//   for (int track_counter = 0; track_counter < nTracks; ++track_counter){

//     TRandom3 waveLengthGenerator;
//     waveLengthGenerator.SetSeed(track_counter);

//     for (size_t iGeantStep = 0; iGeantStep < 1; ++iGeantStep) {

//       double currentTime = 0;
//       geantStep stepBeg = geantStep();
//       stepBeg.coordinates = TVector3(0, 0, 1);

//       cherenkov->SetEnergy(5.02604);

//       double geantStepLength = 1;
//       // int numSmallSteps = (geantStepLength % scattering.GetStepLength() == 0) ? geantStepLength / scattering.GetStepLength() : geantStepLength / scattering.GetStepLength() + 1;
//       // May be mistake due to step number
//       int numSmallSteps = geantStepLength / scattering.GetStepLength() + 1;
//       geantStep tmpBeg = stepBeg;

//       for (int iSmallStep = 0; iSmallStep < numSmallSteps; ++iSmallStep) {
//         // TVector3 nextPoint = scattering.GenerateStep(tmpBeg.coordinates, stepEnd.coordinates);
//         TVector3 nextPoint = scattering.GenerateStepRand(tmpBeg.coordinates);

//         scattered_traj << iGeantStep << "\t" << nextPoint.X() << "\t" << nextPoint.Y()<< "\t" << nextPoint.Z() << "\n";

//         // Loop through angles of observation
//         for (int i = 0; i < cherenkov->GetNBins(); i++) {
//           double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);  

//           TComplex resultPower = cherenkov->CoherentKarinsModel(theta,
//            tmpBeg.coordinates, nextPoint, currentTime, 0., false);
//           hPowerRadiatedRe->Fill(theta, resultPower.Re());      
//           hPowerRadiatedIm->Fill(theta, resultPower.Im());
//         }

//         currentTime += (nextPoint - tmpBeg.coordinates).Mag() / cherenkov->GetVelocity();
//         tmpBeg.coordinates = nextPoint;
//         // double energyLoses = calculateEnergyLoses(cherenkov) * scattering.GetStepLength();
//         // cherenkov->SetEnergy(cherenkov->GetEnergy() - averagedEnergyLosesPerStep);

//       }

//         // Loop through angles of observation and square the result
//       for (int i = 0; i < cherenkov->GetNBins(); i++) {
//         double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);  

//         double im = hPowerRadiatedIm->GetBinContent(i+1);
//         double re = hPowerRadiatedRe->GetBinContent(i+1);
//         double resultIntegral = re * re + im * im;

//         // Correct-factor S~R^2
//         hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta) * TMath::Sin(theta));
//       } // End of loop through angles

//     hPowerRadiatedIm->Reset();
//     hPowerRadiatedRe->Reset();

//     }


//   }

//   time(&end);
//   double seconds = difftime(end, start);
//   std::cout << "Running time: " << seconds << std::endl; // искомое время

//   std::cout << "Last energy " << cherenkov->GetEnergy() << std::endl;
//   std::cout << "Q " << scattering.calculateQ() * stepLength / cherenkov->GetVelocity()<< std::endl;
//   std::cout << "q " << scattering.calculateQ() << std::endl;
//   std::cout << "Omega " << cherenkov->GetOmega() << std::endl;
//   std::cout << "Estimated from theory peak width: " << TMath::Power(scattering.calculateQ() /cherenkov->GetOmega(), 1./3) *
//   TMath::Power(TMath::Sqrt(1 - cherenkov->GetCherenkovCos() * cherenkov->GetCherenkovCos()) * cherenkov->GetRefractiveIndex() * cherenkov->GetBeta(), -1./3) << std::endl;

//   // TCanvas* gCanvas = nullptr;
//   // Make plots square
//   // Double_t w = 1000;
//   // Double_t h = 800;
//   // gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
//   // gCanvas->SetLogy();
//   // gStyle->SetOptTitle(kFALSE);
//   // gStyle->SetOptStat(0);
//   // gStyle->SetPalette(kDeepSea);

//   // Double_t scale = 1. / hPowerRadiatedSquared->Integral();
//   // hPowerRadiatedSquared->Scale(scale);

//   // hPowerRadiatedSquared->SetMarkerStyle(kFullTriangleDown);
//   // hPowerRadiatedRe->Draw("HIST");
//   // hPowerRadiatedIm->SetLineColor(kRed);
//   // hPowerRadiatedIm->Draw("SAMEHIST");

//   // hPowerRadiatedSquared->SetLineColor(kGreen);
//   // hPowerRadiatedSquared->Draw("HISTSAME");

//   // gPad->BuildLegend();

//   // std::string outpotFilename = "histograms_output/CherenkovSimulation_wCorrection_q" + std::to_string(qMultiplier) + "_part_" + std::to_string(inputNumber) + ".root";
//   std::string outputFilename = "Test_for_q_rot_ownTraj_" + std::to_string(qMultiplier) + ".root";
//   TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");
//   outputFile->cd();

//   hPowerRadiatedSquared->Write();
//   outputFile->Close();
//   scattered_traj.close();


//   // std::string fileName = "TestStepLength_" + std::to_string(int(cherenkov->GetEnergy())) 
//   // + "MeV_" + std::to_string(stepLength) + "Length";
//   // cherenkov->PrintDrawnHistogram(gCanvas, fileName.data());
// }


// 
void simulationSingleScattering() {

  CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
  auto scattering = ParticleScattering(cherenkov);
  scattering.ParseTrajectory("electronTrajectories/Trajectory_coordinates_1MeV_1000.txt", 200);

  //Show beta
  std::cout << "Beta: " << cherenkov->GetBeta() << std::endl;

  // The Cherenkov Radiation Description.
  std::cout << "The expected Cherenkov angle should be at theta = " << TMath::ACos(cherenkov->GetCherenkovCos()) << std::endl;
  double binsn = 100;

  cherenkov->SetNBins(binsn);

  cherenkov->SetWaveLenght(410 * 1e-7);
  cherenkov->SetRefractiveIndex(1.33668);

  TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
   "Integral Result sum CoherentKarinsModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
  TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
   "Integral Result sum CoherentKarinsModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());

  double suaredHistBins = 1000;

  TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
   "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

  int track_counter = 0;

  for (auto track:scattering.allElectrons){
    track_counter++;
    // 
    std::cout << "Track points: " << track.size() << " number: " << track_counter << std::endl;
    if (track.size() < 2) continue;

    time_t start, end;
    time(&start);

    for (size_t iGeantStep = 0; iGeantStep < track.size() - 1; ++iGeantStep) {
      // double waveLengthGenerated = waveLengthGenerator.Gaus(450 * 1e-7, 50 * 1e-7);
      // cherenkov->SetWaveLenght(waveLengthGenerated);
      // std::cout << waveLengthGenerated << std::endl;

      double currentTime = 0;
      geantStep stepBeg = track.at(iGeantStep);
      geantStep stepEnd = track.at(iGeantStep + 1);
      // Threshold for Cherenkov radiation
      cherenkov->SetMomentum(stepBeg.p);
      if (cherenkov->GetBeta() < 1./cherenkov->GetRefractiveIndex()) break;

      // Loop through angles of observation
      for (int i = 0; i < cherenkov->GetNBins(); i++) {

        double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);

        TComplex resultPower = cherenkov->CoherentKarinsModel(theta,
         stepBeg.coordinates, stepEnd.coordinates, currentTime, 0., false);

        hPowerRadiatedRe->Fill(theta, resultPower.Re());      
        hPowerRadiatedIm->Fill(theta, resultPower.Im());

      }


      currentTime += (stepEnd.coordinates - stepBeg.coordinates).Mag() / cherenkov->GetVelocity();
      // double energyLoses = calculateEnergyLoses(cherenkov) * scattering.GetStepLength();


        // Loop through angles of observation and square the result
      for (int i = 0; i < cherenkov->GetNBins(); i++) {
        double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);  

        double im = hPowerRadiatedIm->GetBinContent(i+1);
        double re = hPowerRadiatedRe->GetBinContent(i+1);
        double resultIntegral = re * re + im * im;

        // Correct-factor S~R^2
        // Only sin
        // hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta) * TMath::Sin(theta));
        hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta));
      } // End of loop through angles

    hPowerRadiatedIm->Reset();
    hPowerRadiatedRe->Reset();

    }

    time(&end);
    double seconds = difftime(end, start);
    std::cout << "Time for one track: " << seconds << std::endl; // искомое время

  }


  std::cout << "Last beta " << cherenkov->GetBeta() << std::endl;

  // std::string outpotFilename = "histograms_output/CherenkovSimulation_wCorrection_q" + std::to_string(qMultiplier) + "_part_" + std::to_string(inputNumber) + ".root";
  std::string outputFilename = "Single_scattering_1MeV_1000tracks_full_100bins_sin_trsh.root";
  TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");  
  outputFile->cd();

  hPowerRadiatedSquared->Write();
  outputFile->Close();

}