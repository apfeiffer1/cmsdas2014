
//
// $Id: MuScleFitMuonCorrector.cc,v 1.10 2013/05/22 12:46:15 scasasso Exp $
//

/**
  \class    modules::MuScleFitMuonCorrectorT MuScleFitMuonCorrectorT.h 
  \brief    Applies MuScleFit corrections to muons            
  \author   Giovanni Petrucciani (modified by Stefano Casasso)
  \version  $Id: MuScleFitMuonCorrector.cc,v 1.10 2013/05/22 12:46:15 scasasso Exp $
*/


#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "MuScleFit/Calibration/interface/MuScleFitCorrector.h"


namespace modules {

  template<typename T>
  class MuScleFitMuonCorrectorT : public edm::EDProducer {
    public:
      explicit MuScleFitMuonCorrectorT(const edm::ParameterSet & iConfig);
      virtual ~MuScleFitMuonCorrectorT() { }

      virtual void produce(edm::Event & iEvent, const edm::EventSetup & iSetup);

    private:
      /// Labels for input collections
      edm::InputTag src_;
      std::string identifier_;
      bool debug_;
      bool applySmearing_;
      bool fakeSmearing_;

      // MuScleFit corrector
      MuScleFitCorrector* corrector_;
      MuScleFitCorrector* corrector2012D_;
  };

} // namespace


template<typename T>
modules::MuScleFitMuonCorrectorT<T>::MuScleFitMuonCorrectorT(const edm::ParameterSet & iConfig) :
  src_(iConfig.getParameter<edm::InputTag>("src")),
  identifier_(iConfig.getParameter<std::string>("identifier")),
  debug_(iConfig.getParameter<bool>("debug")),
  applySmearing_(iConfig.getParameter<bool>("applySmearing")),
  fakeSmearing_(iConfig.getParameter<bool>("fakeSmearing"))
{

  TString fileName = "";
  TString fileName2012D = "";

  if (identifier_=="Summer12_DR53X_smearPrompt"){ // MC 2012 smeared against Prompt data
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2012_MC_53X_smearPrompt.txt");
  }
  if (identifier_=="Summer12_DR53X_smearReReco"){ // MC 2012 smeared against ReReco data
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2012_MC_53X_smearReReco.txt");
  }
  else if (identifier_=="Data2012_53X"){ // DATA 2012 (Prompt)
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2012ABC_DATA_53X.txt");
    fileName2012D = "MuScleFit/Calibration/data/MuScleFit_2012D_DATA_53X.txt";
  }
  else if (identifier_=="Data2012_53X_ReReco"){ // DATA 2012 (ReReco)
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2012ABC_DATA_ReReco_53X.txt");
    fileName2012D = "MuScleFit/Calibration/data/MuScleFit_2012D_DATA_ReReco_53X.txt";
  }
  else if (identifier_=="Fall11_START42"){ // MC 2011 (42X)
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2011_MC_42X.txt");
  }
  else if (identifier_=="Fall11_START44"){ // MC 2011 (44X)
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2011_MC_44X.txt");
  }
  else if (identifier_=="Data2011_42X"){ // DATA 2011 (42X)
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2011_DATA_42X.txt");
  }
  else if (identifier_=="Data2011_44X"){ // DATA 2011 (44X)
    fileName.Append("MuScleFit/Calibration/data/MuScleFit_2011_DATA_44X.txt");
  }
  else {
    std::cout<<"%MuScleFitCorrector% wrong identifier, choose among:"<<std::endl;
    std::cout<<"  data: 'Data2012_53X', 'Data2012_53X_ReReco', 'Data2011_42X', 'Data2011_44X'"<<std::endl;
    std::cout<<"   MC : 'Summer12_DR53X', 'Fall11_START44', 'Fall11_START42'"<<std::endl;
    exit(1);
  }

  edm::FileInPath fileWithFullPath(fileName.Data());
  corrector_ = new MuScleFitCorrector(fileWithFullPath.fullPath());

  if (identifier_.find("Data2012")!=std::string::npos){
    edm::FileInPath file2012DWithFullPath(fileName2012D.Data());
    corrector2012D_ = new MuScleFitCorrector(file2012DWithFullPath.fullPath());
  }
  
  produces<std::vector<T> >(); 
}

template<typename T>
void 
modules::MuScleFitMuonCorrectorT<T>::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {

  using namespace edm;
  using namespace std;
  
  
  Handle<View<T> > src;
  iEvent.getByLabel(src_, src);
  
  unsigned int nsrc = src->size();
  auto_ptr<vector<T> > out(new vector<T>());
  out->reserve(nsrc);
  
  unsigned int event = (unsigned int)iEvent.id().event(); 
  unsigned int run = (unsigned int)iEvent.id().run(); 
  
  for (unsigned int i = 0; i < nsrc; ++i) {
    T mu = (*src)[i];
    double chg = mu.charge();
    TLorentzVector* p4 = new TLorentzVector(mu.px(),mu.py(),mu.pz(),mu.energy());
    
    
    if(debug_ && event%100==0) {
      cout<<"-- MuScleFitCorrector, debug mode --"<<endl;
      cout<<"   Muon pT (RAW) =        "<<p4->Pt()<<endl;
    }

    if (run < 203773){
      corrector_->applyPtCorrection(*p4,chg);
    } else {
      corrector2012D_->applyPtCorrection(*p4,chg);
    }
    
    
    if(debug_ && event%100==0) {
      cout<<"   Muon pT (CORR) =       "<<p4->Pt()<<endl;
    }
    
    
     if (applySmearing_){
       corrector_->applyPtSmearing(*p4,chg,fakeSmearing_);
       
       if(debug_ && event%100==0) {
	 cout<<"   Muon pT (CORR+SMEAR) = "<<p4->Pt()<<endl;
       }
     }

     if (debug_ && event%100==0) cout<<endl;
     
     
     math::XYZTLorentzVector newP4(p4->Px(),p4->Py(),p4->Pz(),p4->Energy());
     mu.setP4(newP4);
     
     
     out->push_back(mu);
     
     
  }
  
  iEvent.put(out);
}


namespace modules {
  //typedef modules::MuScleFitMuonCorrectorT<reco::Muon>  MuScleFitMuonCorrector;
  typedef modules::MuScleFitMuonCorrectorT<pat::Muon>   MuScleFitPATMuonCorrector;
}

#include "FWCore/Framework/interface/MakerMacros.h"
using namespace modules;
//DEFINE_FWK_MODULE(MuScleFitMuonCorrector);
DEFINE_FWK_MODULE(MuScleFitPATMuonCorrector);
