// -*- C++ -*-
//
// Package:    TrackProperties
// Class:      TrackProperties
//
/**\class TrackProperties TrackProperties.cc Ex1/TrackProperties/src/TrackProperties.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Marco Rovere,40 1-B02,+41227671637,
//         Created:  Mon Okt 28 15:38:59 CET 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TH1.h"

//
// class declaration
//

class TrackProperties : public edm::EDAnalyzer {
 public:
  explicit TrackProperties(const edm::ParameterSet&);
  ~TrackProperties();

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

  void diMuonAnalysis(const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  float  muon_mass_;
  TH1F * h_track_signedpt;
  TH1F * h_track_phi;
  TH1F * h_track_eta;
  TH1F * h_track_dxy;
  TH1F * h_track_dz;
  TH1F * h_track_quality;
  TH1F * h_track_highpurity_signedpt;
  TH1F * h_track_highpurity_phi;
  TH1F * h_track_highpurity_eta;
  TH1F * h_track_highpurity_dxy;
  TH1F * h_track_highpurity_dz;
  TH1F * h_total_px;
  TH1F * h_total_py;
  TH1F * h_total_pz;
  TH1F * h_dimuon;
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
TrackProperties::TrackProperties(const edm::ParameterSet& iConfig)
    : muon_mass_(0.1056) {
  edm::Service<TFileService> fs;   //  now do whatever initialization is needed
  h_track_signedpt = fs->make<TH1F>("Track_SignedPt",
                                    "Track_SignedPt",
                                    200, -100., 100.);
  h_track_phi = fs->make<TH1F>("Track_Phi",
                               "Track_Phi",
                               100, -3.2, 3.2);
  h_track_eta = fs->make<TH1F>("Track_Eta",
                               "Track_Eta",
                               100, -3., 3.);
  h_track_dxy = fs->make<TH1F>("Track_dxy",
                               "Track_dxy",
                               100, -3., 3.);
  h_track_dz = fs->make<TH1F>("Track_dz",
                              "Track_dz",
                              100, -40., 40.);
  h_track_quality = fs->make<TH1F>("Track_quality",
                                   "Track_quality",
                                   9, -1.5, 7.5);
  // Attach proper labels to bins
  TAxis * t = h_track_quality->GetXaxis();
  for (int q = 0; q < reco::TrackBase::qualitySize; ++q) {
    t->SetBinLabel(t->FindBin(q),
                   reco::TrackBase::qualityNames[q].c_str());
  }
  t->SetBinLabel(t->FindBin(-1), "Any");
  h_track_highpurity_signedpt = fs->make<TH1F>("Track_Highpurity_SignedPt",
                                               "Track_Highpurity_SignedPt",
                                               200, -100., 100.);
  h_track_highpurity_phi = fs->make<TH1F>("Track_Highpurity_Phi",
                                          "Track_Highpurity_Phi",
                                          100, -3.2, 3.2);
  h_track_highpurity_eta = fs->make<TH1F>("Track_Highpurity_Eta",
                                          "Track_Highpurity_Eta",
                                          100, -3., 3.);
  h_track_highpurity_dxy = fs->make<TH1F>("Track_Highpurity_dxy",
                                          "Track_Highpurity_dxy",
                                          100, -3., 3.);
  h_track_highpurity_dz = fs->make<TH1F>("Track_Highpurity_dz",
                                         "Track_Highpurity_dz",
                                         100, -40., 40.);
  h_total_px = fs->make<TH1F>("Total_Tracks_px",
                              "Total_Tracks_px",
                              200, -100., 100.);
  h_total_py = fs->make<TH1F>("Total_Tracks_py",
                              "Total_Tracks_py",
                              200, -100., 100.);
  h_total_pz = fs->make<TH1F>("Total_Tracks_pz",
                              "Total_Tracks_pz",
                              200, -100., 100.);
  h_dimuon = fs->make<TH1F>("diMuon",
                            "diMuon",
                            100, 0., 5.);
}


TrackProperties::~TrackProperties() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // We need to force a write here explicitely *before* deleting the
  // histograms, otherwise we will end up with and empty file. The
  // TFileService, by default, will save histograms in its destructor,
  // which is guaranteed to be called *after* the current one.
  edm::Service<TFileService> fs;
  fs->file().Write();

  delete h_track_signedpt;
  delete h_track_phi;
  delete h_track_eta;
  delete h_track_dxy;
  delete h_track_dz;
  delete h_track_quality;
  delete h_track_highpurity_signedpt;
  delete h_track_highpurity_phi;
  delete h_track_highpurity_eta;
  delete h_track_highpurity_dxy;
  delete h_track_highpurity_dz;
  delete h_total_px;
  delete h_total_py;
  delete h_total_pz;
  delete h_dimuon;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackProperties::analyze(const edm::Event& iEvent,
                         const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("generalTracks", tracks);

  reco::TrackCollection::const_iterator iti = tracks->begin();
  reco::TrackCollection::const_iterator ite = tracks->end();

  float total_px, total_py, total_pz;
  total_px = total_py = total_pz = 0.;
  for (int j = 0; iti != ite; ++iti, ++j) {
    h_track_signedpt->Fill(iti->charge()*iti->pt());
    h_track_phi->Fill(iti->phi());
    h_track_eta->Fill(iti->eta());
    h_track_dxy->Fill(iti->dxy());
    h_track_dz->Fill(iti->dz());
    // test all possible quality flags
    for (int q = 0; q < reco::TrackBase::qualitySize; ++q) {
      if (iti->quality(iti->qualityByName(reco::TrackBase::qualityNames[q])))
        h_track_quality->Fill(q);
    }
    h_track_quality->Fill(-1);

    // Fill Histograms only for highPurity tracks
    if (iti->quality(iti->qualityByName("highPurity"))) {
      h_track_highpurity_signedpt->Fill(iti->charge()*iti->pt());
      h_track_highpurity_phi->Fill(iti->phi());
      h_track_highpurity_eta->Fill(iti->eta());
      h_track_highpurity_dxy->Fill(iti->dxy());
      h_track_highpurity_dz->Fill(iti->dz());
    }

    total_px += iti->px();
    total_py += iti->py();
    total_pz += iti->pz();

    if (j < 1)
      std::cout << "    Track " << j
                << " " << iti->charge()*iti->pt()
                << " " << iti->phi()
                << " " << iti->eta()
                << " " << iti->dxy()
                << " " << iti->dz() << std::endl;
  }
  h_total_px->Fill(total_px);
  h_total_py->Fill(total_py);
  h_total_pz->Fill(total_pz);

  // Perform diMuon analyse
  diMuonAnalysis(iEvent, iSetup);
}


void TrackProperties::diMuonAnalysis(const edm::Event & iEvent,
                                     const edm::EventSetup&) {
  using namespace edm;

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel("globalMuons", tracks);

  // Analyse only events with exactly 2 global muons.
  if (tracks->size() != 2)
    return;

  const reco::Track & mu1 = tracks->at(0);
  const reco::Track & mu2 = tracks->at(1);
  float e1 = sqrt(muon_mass_*muon_mass_ + mu1.p()*mu1.p());
  float e2 = sqrt(muon_mass_*muon_mass_ + mu2.p()*mu2.p());
  math::XYZTLorentzVectorF dimuon(mu1.px()+mu2.px(),
                                  mu1.py()+mu2.py(),
                                  mu1.pz()+mu2.pz(),
                                  e1+e2);
  h_dimuon->Fill(dimuon.mass());
}


// --- method called once each job just before starting event loop  ---
void
TrackProperties::beginJob() {}

// --- method called once each job just after ending the event loop  ---
void
TrackProperties::endJob() {}

// --- method called when starting to processes a run  ---
void
TrackProperties::beginRun(edm::Run const&, edm::EventSetup const&) {}

// --- method called when ending the processing of a run  ---
void
TrackProperties::endRun(edm::Run const&, edm::EventSetup const&) {}

// --- method called when starting to processes a luminosity block  ---
void
TrackProperties::beginLuminosityBlock(edm::LuminosityBlock const&,
                                      edm::EventSetup const&) {}

// --- method called when ending the processing of a luminosity block  ---
void
TrackProperties::endLuminosityBlock(edm::LuminosityBlock const&,
                                    edm::EventSetup const&) {}

// --- method fills 'descriptions' with the allowed parameters for the module  ---
void
TrackProperties::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(TrackProperties);
