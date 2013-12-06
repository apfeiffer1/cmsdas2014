#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "UserCode/sm_cms_das/interface/SMEventSummary.h"
#include "UserCode/sm_cms_das/interface/SMUtils.h"

#include "TH1F.h"

using namespace std;
using namespace edm;
using namespace reco;

//
class SMDataAnalyzer : public edm::EDAnalyzer 
{

public:
  SMDataAnalyzer(const edm::ParameterSet &iConfig);
  void beginLuminosityBlock(const edm::LuminosityBlock &iLumi, const edm::EventSetup & iSetup );
  void endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup);
  virtual void beginRun(const edm::Run & iRun, edm::EventSetup const & iSetup); 
  virtual void analyze( const edm::Event &iEvent, const edm::EventSetup &iSetup) ;
 
private:

  //monitoring
  TH1F *obsPU_h, *truePU_h, *cutflow_h, *trigger_h;
  
  //selection configuration
  edm::ParameterSet analysisCfg_;
  
  //tool to trace prescale changes
  HLTConfigProvider hltConfig_;

  //event summary
  SMEventSummary summary_;
};

using namespace std;


//
SMDataAnalyzer::SMDataAnalyzer(const edm::ParameterSet &iConfig) : obsPU_h(0), truePU_h(0)
{
  //configure selection
  analysisCfg_ = iConfig.getParameter<edm::ParameterSet>("cfg");

  //init monitoring tools
  edm::Service<TFileService> fs;
  summary_.create(  fs->make<TTree>("data","Event Summary") );

  obsPU_h   = fs->make<TH1F>( "pileup",     ";Pileup; Events",      100,-0.5,99.5);
  truePU_h  = fs->make<TH1F>( "pileuptrue", ";True pileup; Events", 100,-0.5,99.5);

  cutflow_h = fs->make<TH1F>( "cutflow",    ";Cutflow; Events",     2,0,2);
  string selSteps[]={"RECO","#geq 1 vertex"};
  for(size_t istep=0; istep<sizeof(selSteps)/sizeof(string); istep++)  cutflow_h->GetXaxis()->SetBinLabel(istep+1,selSteps[istep].c_str());

  std::vector<string> trigs=analysisCfg_.getParameter<std::vector<string> >("triggerPaths");
  trigger_h = fs->make<TH1F>( "trigger",    ";Trigger path;",       trigs.size(),0,trigs.size());
  for(size_t itrig=0; itrig<trigs.size(); itrig++) trigger_h->GetXaxis()->SetBinLabel(itrig+1,trigs[itrig].c_str());
}

//
void SMDataAnalyzer::beginLuminosityBlock(const edm::LuminosityBlock&lumi, const edm::EventSetup & setup ) 
{
  edm::Handle<LumiSummary> l;
  lumi.getByLabel("lumiProducer", l); 
  if (!l.isValid())  return;
  summary_.instLumi = l->avgInsDelLumi();
}

//
void SMDataAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup)
{
  try{
    edm::Handle<edm::MergeableCounter> ctrHandle;
    iLumi.getByLabel("startCounter", ctrHandle);
    if(ctrHandle.isValid()){
      cutflow_h->Fill(0.,ctrHandle->value);
    }
  }catch(std::exception){
  }
}

//
void SMDataAnalyzer::beginRun(const edm::Run & iRun, const edm::EventSetup & iSetup) 
{
  bool changed(true);
  hltConfig_.init(iRun, iSetup,"HLT",changed);
}

//
void SMDataAnalyzer::analyze(const edm::Event &event, const edm::EventSetup &iSetup) 
{
  bool isData=event.isRealData();

  summary_.reset();

  //
  // EVENT HEADER
  //
  summary_.run    = event.id().run();
  summary_.lumi   = event.luminosityBlock();
  summary_.event  = event.id().event();

  //
  // GENERATOR LEVEL
  //
  int nGenVectorBosons(0);
  if(!isData){
    
    //pileup
    edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
    event.getByType(puInfoH);
    summary_.ngenITpu    = 0;
    summary_.ngenOOTpu   = 0;
    summary_.ngenOOTpum1 = 0;
    summary_.ngenTruepu  = 0;
    if(puInfoH.isValid())
      {
	for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++)
	  {
	    if(it->getBunchCrossing()==0)      { summary_.ngenITpu += it->getPU_NumInteractions();   summary_.ngenTruepu  = it->getTrueNumInteractions(); }
	    else if(it->getBunchCrossing()<0)  { summary_.ngenOOTpum1 += it->getPU_NumInteractions(); }
	    else                               { summary_.ngenOOTpu   += it->getPU_NumInteractions(); }
	  }
      }
    obsPU_h->Fill(summary_.ngenITpu);
    truePU_h->Fill(summary_.ngenTruepu);

    //pdf info
    edm::Handle<GenEventInfoProduct> genEventInfoProd;
    event.getByType( genEventInfoProd );
    if(genEventInfoProd.isValid())
      {
	summary_.genWeight = genEventInfoProd->weight();
	summary_.qscale = genEventInfoProd->qScale();
	if(genEventInfoProd->pdf())
	  {
	    summary_.qscale = genEventInfoProd->pdf()->scalePDF;
	    summary_.x1  = genEventInfoProd->pdf()->x.first;
	    summary_.x2  = genEventInfoProd->pdf()->x.second;
	    summary_.id1 = genEventInfoProd->pdf()->id.first;
	    summary_.id2 = genEventInfoProd->pdf()->id.second;
	  }
	if(genEventInfoProd->binningValues().size()>0) summary_.pthat = genEventInfoProd->binningValues()[0];
      }

    //matrix element info
    Handle<LHEEventProduct> lheH;
    event.getByType(lheH);
    if(lheH.isValid()) summary_.nup=lheH->hepeup().NUP;

    //generator level event
    summary_.mcn=0;
    Handle<View<Candidate> > genParticlesH;
    event.getByLabel(analysisCfg_.getParameter<edm::InputTag>("genSource"), genParticlesH);

    //analyze hard process (PYTHIA-like status)
    for(size_t i = 0; i < genParticlesH->size(); ++ i)
      {
	const reco::GenParticle & p = dynamic_cast<const reco::GenParticle &>( (*genParticlesH)[i] );
	bool isHardProc(p.status()==3);
	if(!isHardProc) continue;
	if(abs(p.pdgId())==23 || abs(p.pdgId())==24) nGenVectorBosons++;
	summary_.mc_id[summary_.mcn]=p.pdgId();
	summary_.mc_status[summary_.mcn]=p.status();
	summary_.mc_px[summary_.mcn]=p.px();
	summary_.mc_py[summary_.mcn]=p.py();
	summary_.mc_pz[summary_.mcn]=p.pz();
	summary_.mc_en[summary_.mcn]=p.energy();
	summary_.mcn++;
      }
  }
  
  //
  // TRIGGER
  //
  bool triggerHasFired(false);
  edm::InputTag trigSource              = analysisCfg_.getParameter<edm::InputTag>("triggerSource");
  std::vector<std::string> triggerPaths = analysisCfg_.getParameter<std::vector<std::string> >("triggerPaths");
  summary_.tn = triggerPaths.size();  
  edm::Handle<edm::TriggerResults> triggerBitsH;
  event.getByLabel( trigSource, triggerBitsH);
  const edm::TriggerNames &triggerNames = event.triggerNames( *triggerBitsH );
  summary_.t_bits=0;
  for(int i=0; i<summary_.tn; i++) { summary_.t_prescale[i]=0; }
  for (size_t itrig = 0; itrig != triggerBitsH->size(); ++itrig)
    {
      if( !triggerBitsH->wasrun(itrig) ) continue;
      if( triggerBitsH->error(itrig) ) continue;
      if( !triggerBitsH->accept(itrig) ) continue;
      std::string trigName = triggerNames.triggerName(itrig);
      for(size_t it=0; it<triggerPaths.size(); it++)
	{
	  if(trigName.find(triggerPaths[it]) == std::string::npos) continue;
	  summary_.t_bits |= (1 << it);
	  summary_.t_prescale[it]=hltConfig_.prescaleValue(event, iSetup, trigName);
	  triggerHasFired=true;
	  break;
	}
    }
  if(!isData) triggerHasFired=true;

  //
  // VERTEX, BEAM SPOT, ENERGY DENSITY
  //
  edm::Handle<reco::BeamSpot> beamSpotH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("beamSpotSource"), beamSpotH);
  edm::Handle<reco::VertexCollection> vtxH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("vtxSource"), vtxH);
  summary_.nvtx=vtxH->size();
  edm::Handle< double > rho;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("rhoSource"),rho);
  summary_.rho = *rho;
  if(summary_.nvtx==0) return;
  reco::VertexRef primVtx(vtxH,0);
  if(primVtx.isNull()) return;
  cutflow_h->Fill(1);

  //
  // CHARGED LEPTONS
  //
  summary_.ln=0; 

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
  edm::Handle<View<Candidate> > muH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("muonSource"),     muH);
  edm::Handle<View<reco::Track> > tracksH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("trkSource"), tracksH);
  int nMuons(0);
  for(size_t imu=0; imu< muH->size(); ++imu)
    {
      reco::CandidatePtr muonPtr    = muH->ptrAt(imu);
      const pat::Muon *muon         = dynamic_cast<const pat::Muon *>( muonPtr.get() );

      //apply minimal pre-selection
      if(muon->pt()<15 || fabs(muon->eta())>2.5) continue;
      
      //store information
      summary_.ln_id[summary_.ln]                         = -13*muon->charge();
      summary_.ln_px[summary_.ln]                         = muon->px();
      summary_.ln_py[summary_.ln]                         = muon->py();
      summary_.ln_pz[summary_.ln]                         = muon->pz();
      summary_.ln_en[summary_.ln]                         = muon->energy();
      summary_.ln_ecalIso[summary_.ln]                    = muon->isolationR03().emEt;
      summary_.ln_hcalIso[summary_.ln]                    = muon->isolationR03().hadEt;
      summary_.ln_trkIso[summary_.ln]                     = muon->isolationR03().sumPt;
      summary_.ln_gIso[summary_.ln]                       = muon->pfIsolationR03().sumPhotonEt;
      summary_.ln_chIso[summary_.ln]                      = muon->pfIsolationR03().sumChargedHadronPt;
      summary_.ln_puchIso[summary_.ln]                    = muon->pfIsolationR03().sumPUPt;
      summary_.ln_nhIso[summary_.ln]                      = muon->pfIsolationR03().sumNeutralHadronEt;
      if(!muon->innerTrack().isNull())
	{
	  std::pair<bool,Measurement1D> ip3dRes               = mytools::getImpactParameter<reco::TrackRef>(muon->innerTrack(), primVtx, iSetup, true);
	  summary_.ln_ip3d[summary_.ln]                       = ip3dRes.second.value();
	  summary_.ln_ip3dsig[summary_.ln]                    = ip3dRes.second.significance();
	}
      else
	{
	  summary_.ln_ip3d[summary_.ln]                       = -999999999;
	  summary_.ln_ip3dsig[summary_.ln]                    = -999999999;
	}

      //add trigger match
      int TrigSum(0);
      for(size_t it=0; it<triggerPaths.size(); it++)
	{
	  string tempTrigName = triggerPaths[it] + "*";
	  if ( muon->triggerObjectMatchesByPath(tempTrigName).size() > 0 ) TrigSum |= (1<<it); 
	}
      summary_.ln_Tbits[summary_.ln] = TrigSum ;

      //add gen match index
      const reco::Candidate *genLep = muon->genLepton();      
      if(genLep){
	for(Int_t imcn=0; imcn<summary_.mcn; imcn++)
	  {
	    if(summary_.mc_id[imcn]!=genLep->pdgId()) continue;
	    LorentzVector p4(summary_.mc_px[imcn],summary_.mc_py[imcn],summary_.mc_pz[imcn],summary_.mc_en[imcn]);
	    float dR=deltaR(p4.eta(),p4.phi(),genLep->eta(),genLep->phi());
	    if(dR>0.1) continue;
	    summary_.ln_genid[summary_.ln]=imcn;
	    break;
	  }
      }

      //check if it is matched to a track
      bool matchesCtfTrack(false);
      for(size_t itk=0; itk< tracksH->size(); ++itk)
	{
	  Float_t dR=deltaR(muon->eta(),muon->phi(),tracksH->ptrAt(itk)->eta(),tracksH->ptrAt(itk)->phi());
	  if(dR>0.05) continue;
	  matchesCtfTrack=true;
	  break;
	}

      //add id bits
      bool isPF( muon->isPFMuon() );
      bool isGlobal( muon->isGlobalMuon() );
      bool isTracker( muon->isTrackerMuon() );
      Float_t innerTrackChi2 = isTracker && !muon->innerTrack().isNull() ? 
	muon->innerTrack()->normalizedChi2() :                       
	0.;
      Float_t trkchi2 = isGlobal  && !muon->globalTrack().isNull() ? 
	muon->globalTrack()->normalizedChi2() : 
	innerTrackChi2;
      Float_t trkValidPixelHits = isGlobal && !muon->globalTrack().isNull() ? 
	muon->globalTrack()->hitPattern().numberOfValidPixelHits() : 
	!muon->innerTrack().isNull() ?
	muon->innerTrack()->hitPattern().numberOfValidPixelHits():
	-99999;
      Float_t d0  = !muon->muonBestTrack().isNull() ?
	fabs(muon->muonBestTrack()->dxy(primVtx->position())) :
	-999999.;
      Float_t dZ = !muon->muonBestTrack().isNull() ?
	fabs(muon->muonBestTrack()->dz(primVtx->position())) :
	-999999.;
      Float_t validMuonHits = isGlobal &&  !muon->globalTrack().isNull()?
	muon->globalTrack()->hitPattern().numberOfValidMuonHits() : 
	0.;
      Float_t nMatchedStations           = muon->numberOfMatchedStations();
      Float_t trkLayersWithMeasurement   = !muon->track().isNull() ?
	muon->track()->hitPattern().trackerLayersWithMeasurement() :
	0.;
      Float_t pixelLayersWithMeasurement = isTracker && !muon->innerTrack().isNull() ? 
	muon->innerTrack()->hitPattern().pixelLayersWithMeasurement() : 
	0.;
      bool isLoose( isPF && (isGlobal || isTracker) );
      bool isTight( isPF                                 
		    && isGlobal              
		    && fabs(d0)<0.2 	 
		    && fabs(dZ)<0.5         
		    && trkValidPixelHits>0         
		    && trkchi2<10. 
		    && validMuonHits>0.     
		    && nMatchedStations>1   
		    && trkLayersWithMeasurement>5 );
      bool isSoft(isTracker 
		  && muon->muonID("TMOneStationTight") 
		  && fabs(d0)<3.  
		  && fabs(dZ)<30.
		  && trkLayersWithMeasurement>5 
		  && pixelLayersWithMeasurement>1  
		  && innerTrackChi2 < 1.8 );
      bool isHighNew(false);
      try{
	isHighNew=muon::isHighPtMuon(dynamic_cast<const reco::Muon &>(*muon), dynamic_cast<const reco::Vertex &> (*primVtx),reco::improvedTuneP) ;
      }catch(...){
      }
      summary_.ln_idbits[summary_.ln]                     = 
	( (int(muon->muonID("GlobalMuonPromptTight")) & 0x1)   << 0)
	| ( (int(muon->muonID("TMLastStationLoose")) & 0x1)    << 1)
	| ( (int(muon->muonID("TMLastStationTight")) & 0x1)    << 2)
	| ( (int(muon->muonID("TMLastStationAngTight")) & 0x1) << 3)
	| ( (int(muon->muonID("TMOneStationTight")) & 0x1)     << 4)
	| ( isTracker                                          << 5)
	| ( isGlobal                                           << 6)
	| ( isPF                                               << 7)
	| ( isLoose                                            << 8)
	| ( isSoft                                             << 9)
	| ( isTight                                            << 10)
	| ( isHighNew                                          << 11)
	| ( matchesCtfTrack                                    << 12);
      
      summary_.ln++;
      nMuons++;
    }

  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolation
  int nElecs(0);
  edm::Handle<View<Candidate> > eH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("electronSource"), eH);
  edm::Handle<reco::ConversionCollection> convH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("conversionSource"), convH);
  for(size_t iele=0; iele< eH->size(); ++iele)
    {
      reco::CandidatePtr elePtr       = eH->ptrAt(iele);
      const pat::Electron *ele        = dynamic_cast<const pat::Electron *>( elePtr.get() );
      const reco::GsfElectron *gsfEle = dynamic_cast<const reco::GsfElectron *>(ele);
      if(ele->pt()<15 || !(ele->isEB() || ele->isEE()) )                        continue;
      if(ele->gsfTrack().isNull() || ele->superCluster().isNull() || gsfEle==0) continue;
      
      //store information
      summary_.ln_id[summary_.ln]                         = -11*ele->charge();
      summary_.ln_px[summary_.ln]                         = ele->px();
      summary_.ln_py[summary_.ln]                         = ele->py();
      summary_.ln_pz[summary_.ln]                         = ele->pz();
      summary_.ln_en[summary_.ln]                         = ele->energy();
      if(!ele->gsfTrack().isNull())
	{
	  std::pair<bool,Measurement1D> ip3dRes               = mytools::getImpactParameter<reco::GsfTrackRef>(ele->gsfTrack(), primVtx, iSetup, true);
	  summary_.ln_ip3d[summary_.ln]                       = ip3dRes.second.value();
	  summary_.ln_ip3dsig[summary_.ln]                    = ip3dRes.second.significance();
	}
      else
	{
	  summary_.ln_ip3d[summary_.ln]                       = -99999999999;
	  summary_.ln_ip3dsig[summary_.ln]                    = -99999999999;
	}
      summary_.ln_ecalIso[summary_.ln]                    = ele->dr03EcalRecHitSumEt();
      summary_.ln_hcalIso[summary_.ln]                    = ele->dr03HcalTowerSumEt();
      summary_.ln_trkIso[summary_.ln]                     = ele->dr03TkSumPt();
      summary_.ln_gIso[summary_.ln]                       = ele->photonIso();
      summary_.ln_chIso[summary_.ln]                      = ele->chargedHadronIso();
      summary_.ln_puchIso[summary_.ln]                    = ele->puChargedHadronIso();
      summary_.ln_nhIso[summary_.ln]                      = ele->neutralHadronIso();

      //add trigger match
      int TrigSum(0);
      for(size_t it=0; it<triggerPaths.size(); it++)
	{
	  string tempTrigName = triggerPaths[it] + "*";
	  if ( ele->triggerObjectMatchesByPath(tempTrigName).size() > 0 ) TrigSum |= (1<<it); 
	}
      summary_.ln_Tbits[summary_.ln] = TrigSum ;

      //add mc match
      const reco::Candidate *genLep = ele->genLepton();      
      if(genLep){
	for(Int_t imcn=0; imcn<summary_.mcn; imcn++)
	  {
	    if(summary_.mc_id[imcn]!=genLep->pdgId()) continue;
	    LorentzVector p4(summary_.mc_px[imcn],summary_.mc_py[imcn],summary_.mc_pz[imcn],summary_.mc_en[imcn]);
	    float dR=deltaR(p4.eta(),p4.phi(),genLep->eta(),genLep->phi());
	    if(dR>0.1) continue;
	    summary_.ln_genid[summary_.ln]=imcn;
	    break;
	  }
      }


      //save id summary
      Float_t sceta=ele->superCluster()->eta();
      Float_t eop_in=ele->eSuperClusterOverP(); 
      Float_t f_brem=ele->fbrem();
      Float_t deta_in=ele->deltaEtaSuperClusterTrackAtVtx();
      Float_t dphi_in=ele->deltaPhiSuperClusterTrackAtVtx();
      Float_t s_ihih=ele->sigmaIetaIeta();
      Float_t hoe=ele->hadronicOverEm();
      Float_t ooemoop=(1.0/ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy());
      Bool_t isConv  = ConversionTools::hasMatchedConversion(*gsfEle,convH,beamSpotH->position());
      Float_t trkLostInnerHits=ele->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
      Float_t dZ=fabs(ele->gsfTrack()->dz(primVtx->position()));
      Float_t d0=fabs(ele->gsfTrack()->dxy(primVtx->position()));
      summary_.ln_idbits[summary_.ln] = 
	(ele->ecalDrivenSeed()                                       << 0)
	| ( ele->trackerDrivenSeed()                                 << 1)
	| ( EgammaCutBasedEleId::PassEoverPCuts(sceta,eop_in,f_brem) << 2);
      for(size_t iid=0; iid<4; iid++)
	{
	  int id(EgammaCutBasedEleId::VETO);
	  if(iid==1) id=EgammaCutBasedEleId::LOOSE;
	  if(iid==2) id=EgammaCutBasedEleId::MEDIUM;
	  if(iid==3) id=EgammaCutBasedEleId::TIGHT;
	  bool hasId=EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::WorkingPoint(id), ele->isEB(), ele->pt(), ele->eta(), 
						  deta_in,dphi_in,s_ihih,hoe,ooemoop,
						  d0,dZ,
						  0., 0., 0., isConv, trkLostInnerHits, *rho);
	  summary_.ln_idbits[summary_.ln] |=  (hasId << (3+iid));
	}
      
      //increment counters
      summary_.ln++;
      nElecs++;
    }
  
  //
  // JETS
  //
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetID
  summary_.jn=0;
  Handle<pat::JetCollection> jetH;
  event.getByLabel( analysisCfg_.getParameter<edm::InputTag>("jetSource"), jetH);
  for(size_t ijet=0; ijet<jetH->size(); ++ijet)
    {
      edm::RefToBase<pat::Jet> jetRef(edm::Ref<pat::JetCollection>(jetH,ijet));
      const pat::Jet *jet              = &((*jetH)[ijet]);
      const reco::Candidate *genParton = jet->genParton();
      const reco::GenJet *genJet       = jet->genJet();

      //pre-selection (note: raw jet energy must be used)
      float rawJetEn( jet->correctedJet("Uncorrected").energy() );
      float nhf( (jet->neutralHadronEnergy() + jet->HFHadronEnergy())/rawJetEn );
      float nef( jet->neutralEmEnergy()/rawJetEn );
      float cef( jet->chargedEmEnergy()/rawJetEn );
      float chf( jet->chargedHadronEnergy()/rawJetEn );
      float nch    = jet->chargedMultiplicity();
      float nconst = jet->numberOfDaughters();
      bool passLooseId(nhf<0.99  && nef<0.99 && nconst>1);
      bool passMediumId(nhf<0.95 && nef<0.95 && nconst>1);
      bool passTightId(nhf<0.90  && nef<0.90 && nconst>1);
      if(fabs(jet->eta())<2.4) {
	passLooseId  &= (chf>0 && nch>0 && cef<0.99);
	passMediumId &= (chf>0 && nch>0 && cef<0.99);
	passTightId  &= (chf>0 && nch>0 && cef<0.99);
      }
      if(jet->pt()<10 || fabs(jet->eta())>4.7 || !passLooseId) continue;
      
      //save information
      summary_.jn_px[summary_.jn]      = jet->px();
      summary_.jn_py[summary_.jn]      = jet->py();
      summary_.jn_pz[summary_.jn]      = jet->pz();
      summary_.jn_en[summary_.jn]      = jet->energy();
      summary_.jn_torawsf[summary_.jn] = jet->correctedJet("Uncorrected").pt()/jet->pt();
      summary_.jn_genflav[summary_.jn] = jet->partonFlavour();
      summary_.jn_genid[summary_.jn]   = genParton ? genParton->pdgId() : 0;
      summary_.jn_genpx[summary_.jn]   = genJet    ? genJet->px()       : 0;
      summary_.jn_genpy[summary_.jn]   = genJet    ? genJet->py()       : 0;
      summary_.jn_genpz[summary_.jn]   = genJet    ? genJet->pz()       : 0;
      summary_.jn_genen[summary_.jn]   = genJet    ? genJet->energy()   : 0;
      summary_.jn_csv[summary_.jn]     = jet->bDiscriminator("combinedSecondaryVertexBJetTags");
      summary_.jn_area[summary_.jn]    = jet->jetArea();

      //a summary of the id bits
      summary_.jn_idbits[summary_.jn] =
	(passLooseId << 0)
	| (passMediumId << 1)
	| (passTightId << 2);
      summary_.jn++;
    }

  //
  // MISSING TRANSVERSE ENERGY
  //
  std::vector<edm::InputTag> metSources=analysisCfg_.getParameter<std::vector<edm::InputTag> >("metSource");
  summary_.metn=0;
  for(size_t imet=0; imet<metSources.size(); imet++)
    {
      Handle<View<reco::PFMET> > metH;
      event.getByLabel(metSources[imet], metH);
      summary_.met_pt[summary_.metn]    = metH.isValid() ? metH->ptrAt(0)->pt() : 0; 
      summary_.met_phi[summary_.metn]   = metH.isValid() ? metH->ptrAt(0)->phi() : 0; 
      summary_.met_sigx2[summary_.metn] = metH.isValid() ? metH->ptrAt(0)->getSignificanceMatrix()(0,0) : 0;
      summary_.met_sigxy[summary_.metn] = metH.isValid() ? metH->ptrAt(0)->getSignificanceMatrix()(0,1) : 0;
      summary_.met_sigy2[summary_.metn] = metH.isValid() ? metH->ptrAt(0)->getSignificanceMatrix()(1,1) : 0;
      Float_t significance(0.);
      if(metH.isValid() && summary_.met_sigx2[summary_.metn]<1.e10 && summary_.met_sigy2[summary_.metn]<1.e10) significance = metH->ptrAt(0)->significance();
      summary_.met_sig[summary_.metn] = significance;
      summary_.met_sumet[summary_.metn] = metH.isValid() ? metH->ptrAt(0)->sumEt() : 0;
      summary_.met_chsumet[summary_.metn] = metH.isValid() ? metH->ptrAt(0)->chargedHadronEt()+metH->ptrAt(0)->electronEt()+metH->ptrAt(0)->muonEt() : 0;
      summary_.metn++;
    }    

  //
  // SUPER CLUSTERS (for tag and probe)
  //
  summary_.scn=0;
  edm::Handle<reco::SuperClusterCollection> superClustersH;
  event.getByLabel(analysisCfg_.getParameter<edm::InputTag>("scSource"),superClustersH );
  for(reco::SuperClusterCollection::const_iterator scIt = superClustersH->begin(); scIt != superClustersH->end(); scIt++)
    {
      if(scIt->energy()<15 || fabs(scIt->eta())>2.5)  continue;
      summary_.scn_e[summary_.scn]  =scIt->energy();
      summary_.scn_eta[summary_.scn]=scIt->eta();
      summary_.scn_phi[summary_.scn]=scIt->phi();
      summary_.scn++;
    }
  
  //all done here: fill if at least one supercluster or lepton is found
  if(!isData && nGenVectorBosons==1) summary_.fill();
  else if(triggerHasFired && (summary_.scn>0 || summary_.ln>0)) summary_.fill();
}




DEFINE_FWK_MODULE(SMDataAnalyzer);

