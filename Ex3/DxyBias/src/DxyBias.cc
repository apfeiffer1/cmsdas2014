// -*- C++ -*-
//
// Package:    DxyBias
// Class:      Vertices
//
/**\class DxyBias DxyBias.cc Ex3/DxyBias/src/DxyBias.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/


// system include files
#include <algorithm>
#include <memory>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TProfile.h"

//
// class declaration
//

class DxyBias : public edm::EDAnalyzer {
 public:
  explicit DxyBias(const edm::ParameterSet&);
  ~DxyBias();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


 private:
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
                                    edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&,
                                  edm::EventSetup const&);

  TProfile * dxy_vs_phi_000;
  TProfile * dxy_vs_phi_beamspot;

  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DxyBias::DxyBias(const edm::ParameterSet& iConfig) {
  edm::Service<TFileService> fs;   //now do what ever initialization is needed

  dxy_vs_phi_000 = fs->make<TProfile>("dxy_vs_phi_000",
                                      "dxy_vs_phi_000",
                                      100, -3.14, 3.14);
  dxy_vs_phi_beamspot =fs->make<TProfile>("dxy_vs_phi_beamspot",
                                          "dxy_vs_phi_beamspot",
                                          100, -3.14, 3.14) ;
}


DxyBias::~DxyBias() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // We need to force a write here explicitely *before* deleting the
  // histograms, otherwise we will end up with and empty file. The
  // TFileService, by default, will save histograms in its destructor,
  // which is guaranteed to be called *after* the current one.
  edm::Service<TFileService> fs;
  fs->file().Write();

  delete dxy_vs_phi_000;
  delete dxy_vs_phi_beamspot;

}



// ------------ method called for each event  ------------
void
DxyBias::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // General Tracks + BeamSpot
  edm::Handle<std::vector<reco::Track> > tracks;
  iEvent.getByLabel("generalTracks", tracks);
  edm::Handle<reco::BeamSpot> beamspot;
  iEvent.getByLabel("offlineBeamSpot", beamspot );

  for(  std::vector<reco::Track>::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
    //compute dxy in relation to origin
    dxy_vs_phi_000->Fill(itTrack->phi(), itTrack->dxy(ROOT::Math::XYZPoint(0, 0, 0)));
    //compute dxy in relation to beamspot
    dxy_vs_phi_beamspot->Fill(itTrack->phi(), itTrack->dxy( *beamspot.product()));
  }

  dxy_vs_phi_000->SetAxisRange(-0.2, 0.2, "Y");
  dxy_vs_phi_beamspot->SetAxisRange(-0.2, 0.2, "Y");

}


// ------------ method called once each job just before starting event loop  ------------
void
DxyBias::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void
DxyBias::endJob() {}

// ------------ method called when starting to processes a run  ------------
void
DxyBias::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void
DxyBias::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void
DxyBias::beginLuminosityBlock(edm::LuminosityBlock const&,
                               edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void
DxyBias::endLuminosityBlock(edm::LuminosityBlock const&,
                             edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DxyBias::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DxyBias);
