#include <MuScleFit/Calibration/interface/MuScleFitCorrector.h>


/// Method to get correct pt
double MuScleFitCorrector::getCorrectPt( const TLorentzVector & lorentzVector , const int & chg) {
  
  // Loop on all the functions and apply them iteratively on the pt corrected by the previous function.
  double pt = lorentzVector.Pt();
  pt = ( scaleFunction_->scale( pt, lorentzVector.Eta(), lorentzVector.Phi(), chg, scaleParArray_) );
  return pt;
}
//----------------------------//


// Return the squared difference of relative pT resolutions data-MC
double MuScleFitCorrector::getSigmaPtDiffSquared(const double & pt, const double & eta){
  
  double sigmaPtData = resolutionFunction_->sigmaPt(pt, eta, resolDataParArray_);
  double sigmaPtMC = resolutionFunction_->sigmaPt(pt, eta, resolMCParArray_);
    return (sigmaPtData*sigmaPtData - sigmaPtMC*sigmaPtMC);
    
}
//----------------------------//


// Method to get correct pt
double MuScleFitCorrector::getSmearedPt( const TLorentzVector & lorentzVector , const int & chg, bool fake=false) {
  
  double pt = lorentzVector.Pt();
  double eta = lorentzVector.Eta();
  double squaredDiff = getSigmaPtDiffSquared(pt,eta);
  if (squaredDiff < 0) return pt;
  double Cfact = 0.;
  if (fileName_.Contains("smearPrompt")) Cfact = 0.85; // smear Summer12 against 2012 Prompt data
  else if (fileName_.Contains("smearReReco")) Cfact = 0.75; // smear Summer12 against 2012 ReReco data
  else if (fileName_.Contains("2011_MC")) Cfact = 0.8; // smear Fall11 against 2011 ReReco data
  else {
    std::cout<<"Are you sure you want to smear data??"<<std::endl;
    exit(1);
  }
  double sPar = Cfact*sqrt(squaredDiff);
  double curv = ((double)chg/pt);
  double normSmearFact = gRandom_->Gaus(0,sPar);
  if (fake) normSmearFact = sPar;
  double smearedCurv = curv + fabs(curv)*normSmearFact;
  
  return 1./((double)chg*smearedCurv);
  
}
//----------------------------//

// Method to apply correction to a TLorentzVector
void MuScleFitCorrector::applyPtCorrection( TLorentzVector& lorentzVector, const int & chg ){

  double eta = lorentzVector.Eta();
  double phi = lorentzVector.Phi();
  double m  = lorentzVector.M(); 
  double corrPt = getCorrectPt(lorentzVector, chg);
  lorentzVector.SetPtEtaPhiM(corrPt,eta,phi,m);
  
}
//----------------------------//


// Method to apply smearing to a TLorentzVector
void MuScleFitCorrector::applyPtSmearing(TLorentzVector& lorentzVector, const int & chg, bool fakeSmear=false){
  double eta = lorentzVector.Eta();
  double phi = lorentzVector.Phi();
  double m  = lorentzVector.M(); 
  double smearedPt = getSmearedPt(lorentzVector, chg, fakeSmear);
  lorentzVector.SetPtEtaPhiM(smearedPt,eta,phi,m);  
}
//----------------------------//




