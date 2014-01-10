// -*- C++ -*-
//
// Package:    Vertices
// Class:      Vertices
//
/**\class Vertices Vertices.cc Ex2/Vertices/src/Vertices.cc

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
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//

class Vertices : public edm::EDAnalyzer {
 public:
  explicit Vertices(const edm::ParameterSet&);
  ~Vertices();

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

  void kshortPeak(const char * label, TH1F * h, const edm::Event &iEvent, bool useCosAngle = false,
                          TH1F * h_cosAngle = NULL, // can be NULL
                          TH1F * h_cosAngleZoom = NULL,// can be NULL
                          TH1F * h_cosBestAngle = NULL, // can be NULL
                          TH1F * h_cosBestAngleZoom = NULL, // can be NULL
                          TH1F * h_mass_cosBestAngle = NULL // can be NULL
                          );

  TH2F * h_secondary_vertices_map;
  TH1F * h_global_loose_invmass;
  TH1F * h_global_tight_invmass;
  TH1F * h_global_hp_invmass;
  TH1F * h_global_cosAngle;
  TH1F * h_global_cosAngle_zoom;

  TH1F * h_global_cosBestAngle;
  TH1F * h_global_cosBestAngle_zoom;
  TH1F * h_global_invmass_cosBestAngle;

  TH2F * h_primary_vertices_map;
  TH1F * h_primary_vertices_distance;
  TH1F * h_pileup;
  TH1F * h_track_per_pv;
  TH2F * h_tracks_vs_pv;

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
Vertices::Vertices(const edm::ParameterSet& iConfig) {
  edm::Service<TFileService> fs;   //now do what ever initialization is needed
  h_secondary_vertices_map  = fs->make<TH2F>("Secondary_vertices_map",
                                             "Secondary_vertices_map",
                                             1200, -300., 300.,
                                             280, 0., 120. );
  h_global_loose_invmass = fs->make<TH1F>("Global_Loose_InvMass",
                                          "Global_Loose_InvMass",
                                          100, 0.4, 0.6);
  h_global_tight_invmass = fs->make<TH1F>("Global_Tight_InvMass",
                                          "Global_Tight_InvMass",
                                          100, 0.4, 0.6);
  h_global_hp_invmass = fs->make<TH1F>("Global_HP_InvMass",
                                       "Global_HP_InvMass",
                                       100, 0.4, 0.6);

  h_global_invmass_cosBestAngle = fs->make<TH1F>("Global_HP_InvMass_Best",
                                       "Global_HP_InvMass_Best",
                                       100, 0.4, 0.6);

  h_global_cosAngle = fs->make<TH1F>("Global_cosAngle",
                                       "Global_cosAngle",
                                       100, -1.0, 1.0);
  h_global_cosAngle_zoom = fs->make<TH1F>("Global_cosAngle_zoom",
                                       "Global_cosAngle_zoom",
                                       100, 0.99, 1.0);

  h_global_cosBestAngle = fs->make<TH1F>("Global_cosBestAngle",
                                       "Global_cosBestAngle",
                                       100, -1.0, 1.0);

  h_global_cosBestAngle_zoom = fs->make<TH1F>("Global_cosBestAngle_zoom",
                                       "Global_cosBestAngle_zoom",
                                       100, 0.99, 1.0);

  h_primary_vertices_map  = fs->make<TH2F>("Primary_vertices_map",
                                           "Primary_vertices_map",
                                           120, -30., 30.,
                                           400, 0., 2. );
  h_primary_vertices_distance = fs->make<TH1F>("Primary_vertices_distance",
                                               "Primary_vertices_distance",
                                               2500, -25., 25);
  h_pileup = fs->make<TH1F>("PileUp",
                            "PileUp",
                            50, -0.5, 49.5);
  h_track_per_pv = fs->make<TH1F>("Tracks_per_PV",
                                  "Tracks_per_PV",
                                  120, -0.5, 119.5);
  h_tracks_vs_pv = fs->make<TH2F>("NumTracks_vs_NumPV",
                                  "NumTracks_vs_NumPV",
                                  50, -0.5, 49.5,
                                  100, 0., 2000.);
}


Vertices::~Vertices() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // We need to force a write here explicitely *before* deleting the
  // histograms, otherwise we will end up with and empty file. The
  // TFileService, by default, will save histograms in its destructor,
  // which is guaranteed to be called *after* the current one.
  edm::Service<TFileService> fs;
  fs->file().Write();

  delete h_secondary_vertices_map;
  delete h_global_loose_invmass;
  delete h_global_tight_invmass;
  delete h_global_hp_invmass;
  delete h_global_cosAngle;
  delete h_global_cosAngle_zoom;
  delete h_global_cosBestAngle;
  delete h_global_cosBestAngle_zoom;
  delete h_global_invmass_cosBestAngle;

  delete h_primary_vertices_map;
  delete h_primary_vertices_distance;
  delete h_pileup;
  delete h_track_per_pv;
  delete h_tracks_vs_pv;
}


//
// member functions
//

void Vertices::kshortPeak(const char * label,
                          TH1F * h,
                          const edm::Event &iEvent,
                          bool useCosAngle,
                          TH1F * h_cosAngle, // can be NULL
                          TH1F * h_cosAngleZoom, // can be NULL
                          TH1F * h_cosBestAngle, // can be NULL
                          TH1F * h_cosBestAngleZoom, // can be NULL
                          TH1F * h_mass_cosBestAngle // can be NULL
                          ) {
  using namespace edm;

  edm::Handle<reco::VertexCompositeCandidateCollection> secondary_vertices;
  iEvent.getByLabel(label, "Kshort", secondary_vertices);

  edm::Handle<std::vector<reco::Vertex> > primary_vertices;
  iEvent.getByLabel("offlinePrimaryVertices", primary_vertices);

  reco::VertexCompositeCandidateCollection::const_iterator ivi =
      secondary_vertices->begin();
  reco::VertexCompositeCandidateCollection::const_iterator ive =
      secondary_vertices->end();

  const double cutBestAngle = 0.99;

  // iterate over secondary vertices
  for (; ivi != ive; ++ivi){

    double bestAngle = -1.1f;
    double bestMass = 0.0f;

    // plot mass
    h->Fill(ivi->mass());

    std::vector<reco::Vertex>::const_iterator ipvi =
      primary_vertices->begin();
    std::vector<reco::Vertex>::const_iterator ipve =
      primary_vertices->end();

    // iterate over primary vertices
    for (; ipvi != ipve; ipvi ++ ) {
      // use vertex() on the VertexCompositeCandidateCollection
      // and position on the reco::Vertex
      math::XYZVectorD displacement = ivi->vertex() - ipvi->position();

      if ( useCosAngle ) {
        // compute cosAngle between displacement vector and secondary vertex momentum
        // -> normalized dot product
        const double cosAngle = displacement.Dot( ivi->momentum() ) /
                                ( sqrt( displacement.Mag2() ) * sqrt( ivi->momentum().Mag2() ) );

        h_cosAngle->Fill ( cosAngle );
        h_cosAngleZoom -> Fill ( cosAngle );
        if ( cosAngle > bestAngle ) {
           bestMass = ivi->mass();
           bestAngle = cosAngle;
        }
      }
    } // end: iterate over primary vertices

    // make sure there actually was a nested iteration
    if (useCosAngle && ( primary_vertices->size() > 0 )) {
      h_cosBestAngle->Fill ( bestAngle );
      h_cosBestAngleZoom -> Fill ( bestAngle );
      if ( bestAngle > cutBestAngle ) {
        h_mass_cosBestAngle->Fill ( bestMass );
      }
    }

  } // end: iterate over secondary vertices
}

// ------------ method called for each event  ------------
void
Vertices::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<reco::VertexCompositeCandidateCollection> secondary_vertices;
  iEvent.getByLabel("SecondaryVerticesFromLooseTracks", "Kshort", secondary_vertices);

  reco::VertexCompositeCandidateCollection::const_iterator ivi =
      secondary_vertices->begin();
  reco::VertexCompositeCandidateCollection::const_iterator ive =
      secondary_vertices->end();

  for (; ivi != ive; ++ivi)
    h_secondary_vertices_map->Fill(ivi->vz(), ivi->p4().R());

  kshortPeak("SecondaryVerticesFromLooseTracks", h_global_loose_invmass, iEvent);
  kshortPeak("SecondaryVerticesFromTightTracks", h_global_tight_invmass, iEvent);
  // plot high-purity k-short and apply the cosAngle method
  kshortPeak("SecondaryVerticesFromHighPurityTracks", h_global_hp_invmass,
                                                        iEvent,
                                                        true,
                                                        h_global_cosAngle,
                                                        h_global_cosAngle_zoom,
                                                        h_global_cosBestAngle,
                                                        h_global_cosBestAngle_zoom,
                                                        h_global_invmass_cosBestAngle );

  // Primary Vertices+General Tracks
  edm::Handle<std::vector<reco::Vertex> > primary_vertices;
  iEvent.getByLabel("offlinePrimaryVertices", primary_vertices);
  edm::Handle<std::vector<reco::Track> > tracks;
  iEvent.getByLabel("generalTracks", tracks);

  std::vector<reco::Vertex>::const_iterator ipvi =
      primary_vertices->begin();
  std::vector<reco::Vertex>::const_iterator ipve =
      primary_vertices->end();

  h_tracks_vs_pv->Fill(primary_vertices->size(),
                       tracks->size());

  // print information about the event's primary vertices
  std::cout << "Pile-up: " << primary_vertices->size() << std::endl;

  // iterate over primary vertices
  for (; ipvi != ipve; ++ipvi) {

    // print information about this vertex
    std::cout << "   Vertex has nTracks: " << ipvi->nTracks()
              << " and tracksSize:" << ipvi->tracksSize() << std::endl;

    if ( !ipvi->isFake() ) {
      h_pileup->Fill(primary_vertices->size());
      h_track_per_pv->Fill(ipvi->tracksSize());
      h_primary_vertices_map->Fill(ipvi->position().z(),
                                   ipvi->position().Rho());
    }
    // iterate over the i+1 remaining vertices
    std::vector<reco::Vertex>::const_iterator ipvn = ipvi+1;
    for (; ipvn != ipve; ++ipvn)
      // compute the distance between the vertices
      h_primary_vertices_distance->Fill(ipvi->z()-ipvn->z());
  }


}


// ------------ method called once each job just before starting event loop  ------------
void
Vertices::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void
Vertices::endJob() {}

// ------------ method called when starting to processes a run  ------------
void
Vertices::beginRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when ending the processing of a run  ------------
void
Vertices::endRun(edm::Run const&, edm::EventSetup const&) {}

// ------------ method called when starting to processes a luminosity block  ------------
void
Vertices::beginLuminosityBlock(edm::LuminosityBlock const&,
                               edm::EventSetup const&) {}

// ------------ method called when ending the processing of a luminosity block  ------------
void
Vertices::endLuminosityBlock(edm::LuminosityBlock const&,
                             edm::EventSetup const&) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Vertices::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Vertices);
