#include "CherenkovRadiationModels.h"
#include "UnitedModelsScatteringGeant.h"
#include "ParticleScattering.h"

#include <TH2F.h>
#include <TRandom3.h>
#include <TComplex.h>
#include <TFile.h>

#include <ctime>

/**
  Depricated since Geant4 is used for enery loses calculations
*/
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

/**
  Simple 2D simulation with given number of steps
  */
void runSimulation2D() {
  CherenkovRadiationModels cherenkov;
  //Show beta
  std::cout << "Beta: " << cherenkov.GetBeta() << std::endl;
  cherenkov.SetMomentum(0.629761);
  cherenkov.SetWaveLenght(400 * 1e-7);
  cherenkov.SetRefractiveIndex(1.33668);

  // The Cherenkov Radiation Description.
  std::cout << "The expected Cherenkov angle should be at cos(theta) = " << cherenkov.GetCherenkovCos() << std::endl;
  double binsn;
  std::cin >> binsn;
  cherenkov.SetNBins(binsn);
  double currentTime = 0;

  // File to write coordinates of particle
  std::ifstream file;
  std::string fileNameCoordinates = "geant4coordinates/" + std::to_string(int(cherenkov.GetMomentum())) + "MeVElectron.txt";
  bool test = true;
  if ( test ) {
    std::string testFileName = "electronTrajectories/testing2.txt";
    fileNameCoordinates = testFileName;
  }

  std::cout << "File name: " << fileNameCoordinates << std::endl;
  file.open(fileNameCoordinates.data());
  std::vector<TVector3> allCoordinates = cherenkov.ParseGeantFile(file);
  int nsteps = allCoordinates.size(); // In fact number of steps is n - 1
  std::cout << "Steps: " << nsteps - 1 << std::endl;
  // Set initial position
  TVector3 particlePosition = allCoordinates.at(0);
  
  TH1F* hPowerRadiated2 = new TH1F("hPowerRadiated2", "Power Radiated RE; cos(#theta); Power (normalized)" , cherenkov.GetNBins(), 0, 3.14);
  TH1F* hPowerRadiated = new TH1F("hPowerRadiated", "Power Radiated IM; cos(#theta); Power (normalized)" , cherenkov.GetNBins(), 0, 3.14);

  // double partitionOfLambda, partitionOfOmega;
  // std::cout << "Input a and omegaY" << std::endl;
  // std::cin >> partitionOfLambda >> partitionOfOmega;
  // cherenkov.SetAmplitude2D(cherenkov.GetWaveLength() * partitionOfLambda);
  // cherenkov.SetOmegaY(cherenkov.GetOmega() * partitionOfOmega);

 // Now it is moving with a stepLength and without energy loses
 // light emission is considered to be from the begining point of each line
  for (int j = 1; j < nsteps; j++) {
    TVector3 nextParticlePosition = allCoordinates.at(j);
    
    // Loop through angles of observation
    for (int i = 0; i < cherenkov.GetNBins(); i++) {
      double theta = hPowerRadiated->GetXaxis()->GetBinCenter(i+1);

      TComplex resultPower = cherenkov.CoherentDedricksModel(theta, particlePosition, nextParticlePosition, currentTime, 0., false);
      hPowerRadiated->Fill(theta, resultPower.Im());
      // if (cos > 0.755 && cos < 0.7552 && cos > 2) {
      //   std::cout << "Cos++: " << cos + 2./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos + 2./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "Cos+: " << cos + 1./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos + 1./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "Cos: " << cos << " PowerIm: " << resultPower.Im() << "\n";
      //   std::cout << "Cos-: " << cos - 1./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos - 1./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "Cos--: " << cos - 2./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos - 2./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "###" << std::endl;
      // }
      // double resultPower2 = cherenkov.CoherentMyModel(cos, particlePosition, nextParticlePosition);
      hPowerRadiated2->Fill(theta, resultPower.Re());

    } // End of loop through angles

    particlePosition = nextParticlePosition;
  }

  TCanvas* gCanvas = nullptr;
  // Make plots square
  Double_t w = 1000;
  Double_t h = 800;
  gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
  // gCanvas->SetLogy();
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kDeepSea);

  std::cout << "Integral: " << hPowerRadiated2->Integral();

  TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
   "Power Radiated CoherentMyModel; #theta; dP/d#omega" , cherenkov.GetNBins(), 0, TMath::Pi());// 

    // Loop through angles of observation and square the result
  for (int i = 0; i < cherenkov.GetNBins(); i++) {
    double theta = hPowerRadiatedSquared->GetXaxis()->GetBinCenter(i+1);  

    double im = hPowerRadiated->GetBinContent(i+1);
    double re = hPowerRadiated2->GetBinContent(i+1);
    double resultIntegral = re * re + im * im;
    hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta));      
  } // End of loop through angles

  // hPowerRadiated->SetMarkerStyle(kFullTriangleDown);
  // // hPowerRadiated->Draw("PLC PMC");
  // hPowerRadiated->SetLineColor(kRed);
  // hPowerRadiated->SetMarkerStyle(kFullCircle);
  // hPowerRadiated->Draw("HIST"); 
  // hPowerRadiated2->Draw("SAMEHIST"); 
  hPowerRadiatedSquared->Draw("HIST");

  gPad->BuildLegend();

  // std::string fileName = "Cherenkov_test_2D";
  // cherenkov.PrintDrawnHistogram(gCanvas, fileName.data());
  std::string outputFilename = "Cherenkov_test_2D.root";
  TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");
  outputFile->cd();

  hPowerRadiatedSquared->Write();
  outputFile->Close();
}

/**
  Simple 2D simulationSingleScattering with given number of steps
  */
void runSimulation2D_default() {
  CherenkovRadiationModels cherenkov;
  //Show beta
  cherenkov.SetMomentum(0.629761);
  cherenkov.SetWaveLenght(400 * 1e-7);
  cherenkov.SetRefractiveIndex(1.33668);
  std::cout << "Beta: " << cherenkov.GetBeta() << std::endl;

  // The Cherenkov Radiation Description.
  std::cout << "The expected Cherenkov angle should be at cos(theta) = " << cherenkov.GetCherenkovCos() << std::endl;
  double binsn;
  std::cin >> binsn;
  cherenkov.SetNBins(binsn);
  double currentTime = 0;

  // File to write coordinates of particle
  std::ifstream file;
  std::string fileNameCoordinates = "electronTrajectories/" + std::to_string(int(cherenkov.GetMomentum())) + "MeVElectron.txt";
  bool test = true;
  if ( test ) {
    std::string testFileName = "electronTrajectories/testing2.txt";
    fileNameCoordinates = testFileName;
  }

  std::cout << "File name: " << fileNameCoordinates << std::endl;
  file.open(fileNameCoordinates.data());
  std::vector<TVector3> allCoordinates = cherenkov.ParseGeantFile(file);
  int nsteps = allCoordinates.size(); // In fact number of steps is n - 1
  std::cout << "Steps: " << nsteps - 1 << std::endl;
  // Set initial position
  TVector3 particlePosition = allCoordinates.at(0);

  TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
 "Power Radiated CoherentMyModel; #theta; dP/d#omega" , cherenkov.GetNBins(), 0, TMath::Pi());
  // double partitionOfLambda, partitionOfOmega;
  // std::cout << "Input a and omegaY" << std::endl;
  // std::cin >> partitionOfLambda >> partitionOfOmega;
  // cherenkov.SetAmplitude2D(cherenkov.GetWaveLength() * partitionOfLambda);
  // cherenkov.SetOmegaY(cherenkov.GetOmega() * partitionOfOmega);

 // Now it is moving with a stepLength and without energy loses
 // light emission is considered to be from the begining point of each line
  for (int j = 1; j < nsteps; j++) {
    TVector3 nextParticlePosition = allCoordinates.at(j);
    
    // Loop through angles of observation
    for (int i = 0; i < cherenkov.GetNBins(); i++) {
      double theta = hPowerRadiatedSquared->GetXaxis()->GetBinCenter(i+1);

      TComplex resultPower = cherenkov.CoherentMyModel(theta, particlePosition, nextParticlePosition, 0., true, false);
      hPowerRadiatedSquared->Fill(theta, resultPower * TMath::Sin(theta));
      // if (cos > 0.755 && cos < 0.7552 && cos > 2) {
      //   std::cout << "Cos++: " << cos + 2./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos + 2./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "Cos+: " << cos + 1./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos + 1./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "Cos: " << cos << " PowerIm: " << resultPower.Im() << "\n";
      //   std::cout << "Cos-: " << cos - 1./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos - 1./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "Cos--: " << cos - 2./1000000 << " PowerIm: " << cherenkov.CoherentDedricksModel(cos - 2./1000000, particlePosition, nextParticlePosition, currentTime, 0., false).Im() << "\n";   
      //   std::cout << "###" << std::endl;
      // }
      // double resultPower2 = cherenkov.CoherentMyModel(cos, particlePosition, nextParticlePosition);

    } // End of loop through angles

    particlePosition = nextParticlePosition;
  }

  TCanvas* gCanvas = nullptr;
  // Make plots square
  Double_t w = 1000;
  Double_t h = 800;
  gCanvas = new TCanvas("gCanvas", "gCanvas", w, h);
  // gCanvas->SetLogy();
  gStyle->SetOptTitle(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kDeepSea);


  // hPowerRadiated->SetMarkerStyle(kFullTriangleDown);
  // // hPowerRadiated->Draw("PLC PMC");
  // hPowerRadiated->SetLineColor(kRed);
  // hPowerRadiated->SetMarkerStyle(kFullCircle);
  // hPowerRadiated->Draw("HIST"); 
  // hPowerRadiated2->Draw("SAMEHIST"); 
  hPowerRadiatedSquared->Draw("HIST");

  gPad->BuildLegend();
  std::string outputFilename = "Cherenkov_test_2D_square.root";
  TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");
  outputFile->cd();

  hPowerRadiatedSquared->Write();
  outputFile->Close();
}



/**
  Function to run simulation using Geant4 trajectory.
  Name of the input and output file as the number of bins and tracks are hard-coded.
  */
void simulationSingleScattering() {

  // Creating the CherenkovRadiationModels instance and tie it with scattering class
  CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
  auto scattering = ParticleScattering(cherenkov);
  scattering.ParseTrajectory("electronTrajectories/testing.txt", 3000);

  // scattering.ParseTrajectory("electronTrajectories/Trajectory_coordinates_SS_400nm_03MeV.txt", 3000);

  //Show beta
  std::cout << "Beta: " << cherenkov->GetBeta() << std::endl;

  // The Cherenkov Radiation Description.
  std::cout << "The expected Cherenkov angle should be at theta = " << TMath::ACos(cherenkov->GetCherenkovCos()) << std::endl;
  
  // Set up number of bins
  double binsn = 10000;
  cherenkov->SetNBins(binsn);

  cherenkov->SetWaveLenght(400 * 1e-7);
  cherenkov->SetRefractiveIndex(1.33668);
  std::cout << cherenkov->GetOmega() << std::endl;

  // Creating histograms for real and complex part
  TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
   "Integral Result sum CoherentDedricksModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
  TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
   "Integral Result sum CoherentDedricksModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());

  // Defining number if bins in squared (final power) histogram
  double suaredHistBins = 10000;
  TH1D* hPowerRadiatedSquared = new TH1D("hPowerRadiatedSquared",
   "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

  int track_counter = 0;

  // Loop over all parsed tracks
  for (auto track:scattering.allElectrons){
    track_counter++;

    // Output to see progress and measure the time of running
    std::cout << "Track points: " << track.size() << " number: " << track_counter << std::endl;
    if (track.size() < 2) continue;

    time_t start, end;
    time(&start);

    // Time variable for the simulation
    double currentTime = 0;

    // Loop over track's steps
    for (size_t iGeantStep = 0; iGeantStep < track.size() - 1; ++iGeantStep) {
      // double waveLengthGenerated = waveLengthGenerator.Gaus(450 * 1e-7, 50 * 1e-7);
      // cherenkov->SetWaveLenght(waveLengthGenerated);
      // std::cout << waveLengthGenerated << std::endl;

      geantStep stepBeg = track.at(iGeantStep);
      geantStep stepEnd = track.at(iGeantStep + 1);

      // Threshold for Cherenkov radiation
      cherenkov->SetMomentum(stepBeg.p);
      if (cherenkov->GetBeta() < 0.4) break;

      // if (cherenkov->GetBeta() < 1./cherenkov->GetRefractiveIndex()) break;

      // Loop through angles of observation
      for (int i = 0; i < cherenkov->GetNBins(); i++) {

        double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);

        TComplex resultPower = cherenkov->CoherentDedricksModel(theta,
         stepBeg.coordinates, stepEnd.coordinates, currentTime, 0., false);

        hPowerRadiatedRe->Fill(theta, resultPower.Re());
        hPowerRadiatedIm->Fill(theta, resultPower.Im());

      }

      // Adding a time for the step
      currentTime += (stepEnd.coordinates - stepBeg.coordinates).Mag() / cherenkov->GetVelocity();
      // double energyLoses = calculateEnergyLoses(cherenkov) * scattering.GetStepLength();

      } 

    std::cout << "Last beta " << cherenkov->GetBeta() << std::endl;

    // Loop through angles of observation and square the result for the track
    for (int i = 0; i < cherenkov->GetNBins(); i++) {
      double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);

      double im = hPowerRadiatedIm->GetBinContent(i+1);
      double re = hPowerRadiatedRe->GetBinContent(i+1);
      double resultIntegral = re * re + im * im;

      // Correct-factor S~R^2
      // Only sin
      // hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta) * TMath::Sin(theta));
      hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta));
    }

    // Reseting the histograms with Re and Im part
    hPowerRadiatedIm->Reset();
    hPowerRadiatedRe->Reset();

    time(&end);
    double seconds = difftime(end, start);
    std::cout << "Time for one track: " << seconds << std::endl;

  }

  std::cout << "Last beta " << cherenkov->GetBeta() << std::endl;

  hPowerRadiatedSquared->Draw("HIST");


  // std::string outpotFilename = "histograms_output/CherenkovSimulation_wCorrection_q" + std::to_string(qMultiplier) + "_part_" + std::to_string(inputNumber) + ".root";
  // std::string outputFilename = "Single_scattering_03MeV_full_100bins_400nm.root";
  std::string outputFilename = "Cherenkov_test_2D_model.root";

  TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");
  outputFile->cd();

  hPowerRadiatedSquared->Write();
  outputFile->Close();

}

void calculate_light_output(){
  // Tables used in SNO+ and in my geant4 simulation

  std::vector<double> wavelengthVector = {200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0,
   340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0,
    540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0, 720.0,
     740.0, 760.0, 780.0, 800.0}; // in nm

  std::vector<double> refractiveIndexVector = {1.41615, 1.39727, 1.38395, 1.37414, 1.36667, 1.36082, 1.35615,
   1.35233, 1.34916, 1.3465, 1.34423, 1.34227, 1.34058, 1.33909, 1.33778, 1.33661, 1.33557, 1.33463,
   1.33378, 1.33301, 1.33231, 1.33167, 1.33108, 1.33054, 1.33004,
   1.32957, 1.32914, 1.32874, 1.32836, 1.32801, 1.32768};

  // std::string outputFilename = "SS_03MeV_1500tracks_allwavelengths_fixk.root";
  std::string outputFilename = "all_wavelengths_05MeV_1000bins.root";

  TFile* outputFile = new TFile(outputFilename.c_str(), "recreate");
  outputFile->cd();

     // Creating the CherenkovRadiationModels instance and tie it with scattering class
  CherenkovRadiationModels *cherenkov = new CherenkovRadiationModels;
  auto scattering = ParticleScattering(cherenkov);
  // scattering.ParseTrajectory("electronTrajectories/Trajectory_coordinates_03MeV_1500_P.txt", 3000);
  scattering.ParseTrajectory("electronTrajectories/Trajectory_coordinates_05MeV_1500_P.txt", 3000);

  
  // Set up number of bins
  double binsn = 1000;
  cherenkov->SetNBins(binsn);

  for (size_t iRI = 0; iRI < wavelengthVector.size(); ++iRI){

    cherenkov->SetWaveLenght(wavelengthVector.at(iRI) * 1e-7);
    cherenkov->SetRefractiveIndex(refractiveIndexVector.at(iRI));

    // Creating histograms for real and complex part
    TH1D* hPowerRadiatedRe = new TH1D("hPowerRadiatedReal",
     "Integral Result sum CoherentDedricksModel Real; #theta; Integral result real part" , cherenkov->GetNBins(), 0, TMath::Pi());
    TH1D* hPowerRadiatedIm = new TH1D("hPowerRadiatedComplex",
     "Integral Result sum CoherentDedricksModel Im; #theta; Integral result im part" , cherenkov->GetNBins(), 0, TMath::Pi());

    // Defining number if bins in squared (final power) histogram
    double suaredHistBins = 10000;
    std::string outputHistName = "hPowerRadiatedSquared_" + std::to_string(wavelengthVector.at(iRI));
    TH1D* hPowerRadiatedSquared = new TH1D(outputHistName.data(),
     "Power Radiated CoherentMyModel; #theta; Power (normalized)" , suaredHistBins, 0, TMath::Pi());

    int track_counter = 0;


    time_t start, end;
    time(&start);

    // Loop over all parsed tracks
    for (auto track:scattering.allElectrons){
      track_counter++;

      // Output to see progress and measure the time of running
      // std::cout << "Track points: " << track.size() << " number: " << track_counter << std::endl;
      if (track.size() < 2) continue;

      // Time variable for the simulation
      double currentTime = 0;

      // Loop over track's steps
      for (size_t iGeantStep = 0; iGeantStep < track.size() - 1; ++iGeantStep) {
        // double waveLengthGenerated = waveLengthGenerator.Gaus(450 * 1e-7, 50 * 1e-7);
        // cherenkov->SetWaveLenght(waveLengthGenerated);
        // std::cout << waveLengthGenerated << std::endl;

        geantStep stepBeg = track.at(iGeantStep);
        geantStep stepEnd = track.at(iGeantStep + 1);

        // Threshold for Cherenkov radiation
        cherenkov->SetMomentum(stepBeg.p);
        if (cherenkov->GetBeta() < 0.4) break;

        // if (cherenkov->GetBeta() < 1./cherenkov->GetRefractiveIndex()) break;

        // Loop through angles of observation
        for (int i = 0; i < cherenkov->GetNBins(); i++) {

          double theta = hPowerRadiatedRe->GetXaxis()->GetBinCenter(i+1);

          TComplex resultPower = cherenkov->CoherentDedricksModel(theta,
           stepBeg.coordinates, stepEnd.coordinates, currentTime, 0., false);

          hPowerRadiatedRe->Fill(theta, resultPower.Re());
          hPowerRadiatedIm->Fill(theta, resultPower.Im());

        }

        // Adding a time for the step
        currentTime += (stepEnd.coordinates - stepBeg.coordinates).Mag() / cherenkov->GetVelocity();
        // double energyLoses = calculateEnergyLoses(cherenkov) * scattering.GetStepLength();

        } 

      // Loop through angles of observation and square the result for the track
      for (int i = 0; i < cherenkov->GetNBins(); i++) {
        double theta = hPowerRadiatedIm->GetXaxis()->GetBinCenter(i+1);

        double im = hPowerRadiatedIm->GetBinContent(i+1);
        double re = hPowerRadiatedRe->GetBinContent(i+1);
        double resultIntegral = re * re + im * im;

        // Correct-factor S~R^2
        // Only sin
        // hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta) * TMath::Sin(theta));
        hPowerRadiatedSquared->Fill(theta, resultIntegral * TMath::Sin(theta) * suaredHistBins / binsn);
      }

      // Reseting the histograms with Re and Im part
      hPowerRadiatedIm->Reset();
      hPowerRadiatedRe->Reset();

    }
    time(&end);
    double seconds = difftime(end, start);
    std::cout << "Time for " + std::to_string(wavelengthVector.at(iRI)) + " wavelength: " << seconds << std::endl;
    hPowerRadiatedSquared->Write();

    delete hPowerRadiatedRe;
    delete hPowerRadiatedIm;
    delete hPowerRadiatedSquared;

    std::cout << "Last beta " << cherenkov->GetBeta() << std::endl;
  } // End of loop through refractive indecies 


  // std::string outpotFilename = "histograms_output/CherenkovSimulation_wCorrection_q" + std::to_string(qMultiplier) + "_part_" + std::to_string(inputNumber) + ".root";

  outputFile->Close();

}