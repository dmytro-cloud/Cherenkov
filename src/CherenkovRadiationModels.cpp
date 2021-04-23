#include "CherenkovRadiationModels.h"
#include <cmath>

//Overloading operator << for TVector3 for convinient writing into the stream
std::ofstream& operator<<(std::ofstream& os, const TVector3& v) {
  os << v.X() << "\t" << v.Y() << "\t" << v.Z() << "\n";
  return os;
}  

std::ostream& operator<<(std::ostream& os, const TVector3& v) {
  os << v.X() << "\t" << v.Y() << "\t" << v.Z() << "\n";
  return os;
}  


double CherenkovRadiationModels::GetOmega() const {
  return 2*TMath::Pi() * GetLightSpeed() / GetWaveLength();
}

double CherenkovRadiationModels::GetWaveVector() const {
  return GetOmega() / GetLightSpeed(); // сm^-1
}

// double CherenkovRadiationModels::GetBeta() const {
//   return TMath::Sqrt(GetEnergy() * GetEnergy() + 2 * GetEnergy() * GetMass()) / (GetEnergy() + GetMass()); // unitless
// }

double CherenkovRadiationModels::GetBeta() const {
  return TMath::Sqrt(1 - 1 / (1 + TMath::Power(GetMomentum() / GetMass(), 2)) ); // unitless
}

double CherenkovRadiationModels::GetVelocity() const {
  return GetBeta() * GetLightSpeed(); // сm/s
}

double CherenkovRadiationModels::GetCherenkovCos() const {
  return GetLightSpeed() / (GetVelocity() * GetRefractiveIndex());
}

/**
  Some vector algebra to find angle between track direction and observation point
  */
std::pair<double, double> CherenkovRadiationModels::AngleTransform(double angle, TVector3 initialPosition,
   TVector3 nextPosition, double phi) const{

  // if (initialPosition == TVector3(0., 0., 0.)){
  // initialPosition = TVector3(0., 0., 0.001);
  // }
  // double angle = std::acos(cos);
  auto particleMovementVector = nextPosition - initialPosition;

  // Define vector of the dot position on the sphere with the center in (0; 0; 0)
  TVector3 spherePosition(0., 0., 0.001);
  spherePosition.SetMag(sphereR);
  spherePosition.SetTheta(angle);
  spherePosition.SetPhi(phi);
  // Radiation comes from the beginning of step
  auto deltaVector = spherePosition - initialPosition;
  double obsL = deltaVector.Mag();

  double angleInParticleSystem = deltaVector.Angle(particleMovementVector);
  // if (angle == 3.14 / 10000 * 9999.5) //std::cout << angleInParticleSystem << " ";
  return std::make_pair(obsL, angleInParticleSystem);
}

/**
  Actually, almost the same that Geant4 does but taking into account the finite length of step
  */
double CherenkovRadiationModels::CoherentMyModel(double theta, TVector3 initialPosition, TVector3 nextPosition,
       double phi /* = 0 default*/, bool returnSquared /* = true default*/, bool SI) {
  
  auto [obsL, angleInParticleSystem] = AngleTransform(theta, initialPosition, nextPosition, phi);
  double stepLength = (nextPosition - initialPosition).Mag();
  double SICoef = 10e-7; // erg  = 10e-7 W

  // Build up an expression
  double bracketsExpr = 1/GetVelocity() - n * TMath::Cos(angleInParticleSystem)/c;
  double sinAndDenominator = TMath::Sin(angleInParticleSystem) * TMath::Sin(GetOmega() * stepLength * bracketsExpr / 2) / bracketsExpr;
  double coefficient = electronCharge * n  / (TMath::Pi() * TMath::Power(c, 3) * TMath::Power(obsL, 2) );     
  double resultPower = returnSquared ? coefficient * sinAndDenominator * sinAndDenominator : TMath::Sqrt(coefficient) * sinAndDenominator; // Divide by 2pi???
  if (SI == true && returnSquared == true) resultPower *= SICoef;

  return resultPower;
}

/**
  CURRENTLY IN USE
  Modified Dedrick's analitical formula to calcualte differential density of power
  Irradiated by charge particle that moves in the volume. Considers the waves itself.
  Returns complex value for the wave created during the step.
  */
TComplex CherenkovRadiationModels::CoherentDedricksModel(double theta, TVector3 initialPosition, TVector3 nextPosition, double currentTime,
        double phi /* = 0 default*/, bool returnSquared /* = false default*/) {
  double time = currentTime; // time of start of step
  auto [obsL, angleInParticleSystem] = AngleTransform(theta, initialPosition, nextPosition, phi);
  double observedAngle = theta;
  double stepLength = (nextPosition - initialPosition).Mag();
  
  // Build up expression for the integral
  double normalization = electronCharge * TMath::Sin(angleInParticleSystem) / 2 / TMath::Pi();
  double exp1 = GetOmega() * time - GetRefractiveIndex() * GetWaveVector() * ( initialPosition.X() * TMath::Sin(observedAngle) * TMath::Cos(phi) + 
    initialPosition.Y() * TMath::Sin(observedAngle) * TMath::Sin(phi) + initialPosition.Z() * TMath::Cos(theta) );//position.Dot(obsDir.Unit());
  double exp2 = GetOmega() * stepLength * (1/GetVelocity() - n * TMath::Cos(angleInParticleSystem)/ GetLightSpeed());
  double denominator = GetOmega() * (1/GetVelocity() - n * TMath::Cos(angleInParticleSystem)/ GetLightSpeed());

  TComplex compExp1 = TComplex(TMath::Cos(exp1), TMath::Sin(exp1));
  TComplex compExp2 = TComplex(TMath::Sin(exp2), -TMath::Cos(exp2));

  TComplex I = normalization * compExp1 * (compExp2 + TComplex(0, 1)) / denominator;
  TComplex I2 = I * TComplex::Conjugate(I);
  double realI2 = I2.Re();
  TComplex power = TComplex(-999, -999);

 // Remove L dependency (?) (should be before integrals sum)
  if (returnSquared){
    power = GetWaveVector()*GetWaveVector() / (2*TMath::Pi()*c*obsL*obsL) * I2; // Needs n?
  } else {
    power = GetRefractiveIndex() * GetWaveVector() / ( TMath::Sqrt(2 * TMath::Pi() * c) * obsL ) * I;
  }

  return power;
}


double CherenkovRadiationModels::CoherentSinusoidalModel(double theta, TVector3 initialPosition, TVector3 nextPosition,
      double phi /* = 0 default*/, bool returnSquared /* = true default*/) {

  auto [obsL, angleInParticleSystem] = AngleTransform(theta, initialPosition, nextPosition, phi);
  double stepLength = (nextPosition - initialPosition).Mag();

  // Build up an expression
  double inBracketsFirst = stepLength * GetOmega() *
   (GetLightSpeed() - GetRefractiveIndex() * GetVelocity() * TMath::Cos(angleInParticleSystem))
    / (2 * GetLightSpeed() * GetVelocity() );
  double numeratorFirst = 2 * GetVelocity() * TMath::Sin(angleInParticleSystem) * TMath::Sin(inBracketsFirst);
  double denominatorFirst = GetLightSpeed() * GetOmega() - GetRefractiveIndex() * GetOmega()
   * GetVelocity() * TMath::Cos(angleInParticleSystem);
  double firstTerm = numeratorFirst / denominatorFirst;

  // Common factor for 2-nd and 3-rd terms
  double commonFactor = GetAmplitude2D() * GetOmegaY() * TMath::Cos(angleInParticleSystem);

  double inBracketsSecond = stepLength *
   (GetLightSpeed() * (GetOmega() - GetOmegaY()) - GetRefractiveIndex() * GetVelocity() * GetOmega() * TMath::Cos(angleInParticleSystem)) 
   / (2 * GetLightSpeed() * GetVelocity() );
   double secondTerm = TMath::Sin(inBracketsSecond)
   / (GetLightSpeed() * (GetOmega() - GetOmegaY()) - GetRefractiveIndex() * GetOmega() * GetVelocity() * TMath::Cos(angleInParticleSystem));

   double inBracketsThird = stepLength * 
   (GetLightSpeed() *  (GetOmega() + GetOmegaY()) - GetRefractiveIndex() * GetVelocity() * GetOmega() * TMath::Cos(angleInParticleSystem)) 
   / (2 * GetLightSpeed() * GetVelocity() );
   double thirdTerm = TMath::Sin(inBracketsThird)
   / (GetLightSpeed() *  (GetOmega() + GetOmegaY()) - GetRefractiveIndex() * GetVelocity() * GetOmega() * TMath::Cos(angleInParticleSystem));

   double preSquaredExpressionCoefficient = GetRefractiveIndex() * GetWaveVector() * GetWaveVector()
    * electronCharge * electronCharge / (2 * TMath::Pi() * obsL * obsL * GetLightSpeed() );
   double allExpression = GetLightSpeed() * GetVelocity() * (firstTerm + commonFactor * (secondTerm + thirdTerm));
   double squaredExpression = allExpression * allExpression; // In fact should be conjugated but already real

   double resultPower = returnSquared ? preSquaredExpressionCoefficient * squaredExpression : TMath::Sqrt(preSquaredExpressionCoefficient) * allExpression;

  return resultPower;
}

void CherenkovRadiationModels::PrintDrawnHistogram(TCanvas* gCanvas, const char* basename){
  // mkdir("plots", 0777);
  std::string plotName = "plots/scattered";
  
  plotName += basename;
  
  std::string pngName = plotName;
  pngName += ".png";
  std::string pdfName = plotName;
  pdfName += ".pdf";
  std::string cName = plotName;
  cName += ".C";
  
  // gCanvas->Print(pdfName.c_str(),"pdf");
  gCanvas->Print(pngName.c_str(),"png");
  // gCanvas->Print(cName.c_str());
}

/* Deprecated but still needed by 2D simulation */
std::vector<TVector3> CherenkovRadiationModels::ParseGeantFile(std::ifstream& stream) {
  std::vector<TVector3> allCoordinates;
  double x, y, z;
  std::string firstRaw;
  getline(stream, firstRaw);

  while (stream >> x && stream >> y && stream >> z){
    TVector3 tmp(0, 0, 0);
    tmp.SetX(x);
    /// !!! Now only for 2D
    tmp.SetY(0);  
    /// !!!
    tmp.SetZ(z);
    allCoordinates.push_back(tmp);
  }
  return allCoordinates;
}
