#ifndef smeventsummary_h
#define smeventsummary_h

#include "TTree.h"
#include "TMath.h"
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"

#define MAXDATAOBJECTS 1000

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class SMEventSummary
{
 public:

  Int_t run, lumi, event;
  
  //trigger information
  Int_t tn, t_bits;
  Int_t t_prescale[MAXDATAOBJECTS];

  //pileup information
  Int_t nvtx;
  Float_t instLumi;
  Float_t rho;

  //gen information
  Int_t ngenITpu, ngenOOTpu, ngenOOTpum1, ngenTruepu;
  Float_t pthat, genWeight, qscale, x1,x2;
  Int_t id1, id2, nup;
  Int_t mcn, mc_id[MAXDATAOBJECTS], mc_status[MAXDATAOBJECTS];
  Float_t mc_px[MAXDATAOBJECTS],mc_py[MAXDATAOBJECTS],mc_pz[MAXDATAOBJECTS],mc_en[MAXDATAOBJECTS], mc_lxy[MAXDATAOBJECTS]; 

  //leptons
  Int_t ln;
  Int_t ln_id[MAXDATAOBJECTS],          ln_idbits[MAXDATAOBJECTS],    ln_genid[MAXDATAOBJECTS],   ln_Tbits[MAXDATAOBJECTS];
  Float_t ln_px[MAXDATAOBJECTS],        ln_py[MAXDATAOBJECTS],        ln_pz[MAXDATAOBJECTS],      ln_en[MAXDATAOBJECTS];
  Float_t ln_ecalIso[MAXDATAOBJECTS],   ln_hcalIso[MAXDATAOBJECTS],   ln_trkIso[MAXDATAOBJECTS];
  Float_t ln_gIso[MAXDATAOBJECTS],      ln_chIso[MAXDATAOBJECTS],     ln_puchIso[MAXDATAOBJECTS], ln_nhIso[MAXDATAOBJECTS]; 
  Float_t ln_ip3d[MAXDATAOBJECTS],      ln_ip3dsig[MAXDATAOBJECTS];
  
  //jets
  Int_t jn, jn_idbits[MAXDATAOBJECTS];
  Float_t jn_px[MAXDATAOBJECTS],    jn_py[MAXDATAOBJECTS],      jn_pz[MAXDATAOBJECTS],          jn_en[MAXDATAOBJECTS], jn_torawsf[MAXDATAOBJECTS], jn_csv[MAXDATAOBJECTS], jn_area[MAXDATAOBJECTS];
  Int_t   jn_genflav[MAXDATAOBJECTS], jn_genid[MAXDATAOBJECTS];
  Float_t jn_genpx[MAXDATAOBJECTS], jn_genpy[MAXDATAOBJECTS], jn_genpz[MAXDATAOBJECTS], jn_genen[MAXDATAOBJECTS];
  
  //met 
  Int_t metn;
  Float_t met_pt[MAXDATAOBJECTS], met_phi[MAXDATAOBJECTS],met_sig[MAXDATAOBJECTS],met_sigx2[MAXDATAOBJECTS],met_sigxy[MAXDATAOBJECTS],met_sigy2[MAXDATAOBJECTS],met_sumet[MAXDATAOBJECTS],met_chsumet[MAXDATAOBJECTS];


  //superclusters
  Int_t  scn;
  Float_t scn_e[MAXDATAOBJECTS], scn_eta[MAXDATAOBJECTS], scn_phi[MAXDATAOBJECTS];


  SMEventSummary() { } 
  ~SMEventSummary() { }

  //
  inline void reset()
  {
    run=0;    
    lumi=0;   
    event=0;  
    tn=0;
    mcn=0;
    ln=0; 
    jn=0;
    metn=0;
    scn=0;
  }

  //
  inline void fill() { for(size_t i=0; i<m_trees.size(); i++) m_trees[i]->Fill(); }
  
  //
  inline void attach(TTree *t)
  {
    if(t==0) return;

    //event info
    t->SetBranchAddress("run",      &this->run);
    t->SetBranchAddress("lumi",     &this->lumi);
    t->SetBranchAddress("event",    &this->event);
    
    //trigger bit
    t->SetBranchAddress("tn",        &this->tn);
    t->SetBranchAddress("tbits",     &this->t_bits);
    t->SetBranchAddress("tprescale", this->t_prescale);
    
    //pileup related observables
    t->SetBranchAddress("nvtx",     &this->nvtx);
    t->SetBranchAddress("instLumi", &this->instLumi);
    t->SetBranchAddress("rho",      &this->rho );

    //generator level info
    t->SetBranchAddress("ngenITpu",    &this->ngenITpu);
    t->SetBranchAddress("ngenOOTpu",   &this->ngenOOTpu);
    t->SetBranchAddress("ngenOOTpum1", &this->ngenOOTpum1);
    t->SetBranchAddress("ngenTruepu",  &this->ngenTruepu);
    t->SetBranchAddress("pthat",       &this->pthat);
    t->SetBranchAddress("genWeight",   &this->genWeight);
    t->SetBranchAddress("qscale",      &this->qscale);
    t->SetBranchAddress("x1",          &this->x1);
    t->SetBranchAddress("x2",          &this->x2);
    t->SetBranchAddress("id1",         &this->id1);
    t->SetBranchAddress("id2",         &this->id2);
    t->SetBranchAddress("nup",         &this->nup);
    t->SetBranchAddress("mcn",         &this->mcn);
    t->SetBranchAddress("mc_id",       this->mc_id);
    t->SetBranchAddress("mc_status",   this->mc_status);
    t->SetBranchAddress("mc_px",       this->mc_px);
    t->SetBranchAddress("mc_py",       this->mc_py);
    t->SetBranchAddress("mc_pz",       this->mc_pz);
    t->SetBranchAddress("mc_en",       this->mc_en);
    t->SetBranchAddress("mc_lxy",      this->mc_lxy);  
    
    //selected leptons
    t->SetBranchAddress("ln",            &this->ln);
    t->SetBranchAddress("ln_id",         this->ln_id);
    t->SetBranchAddress("ln_idbits",     this->ln_idbits);
    t->SetBranchAddress("ln_px",         this->ln_px);
    t->SetBranchAddress("ln_py",         this->ln_py);
    t->SetBranchAddress("ln_pz",         this->ln_pz);
    t->SetBranchAddress("ln_en",         this->ln_en);
    t->SetBranchAddress("ln_ecalIso",    this->ln_ecalIso);
    t->SetBranchAddress("ln_hcalIso",    this->ln_hcalIso);
    t->SetBranchAddress("ln_trkIso",     this->ln_trkIso);
    t->SetBranchAddress("ln_gIso",       this->ln_gIso);
    t->SetBranchAddress("ln_chIso",      this->ln_chIso);
    t->SetBranchAddress("ln_puchIso",    this->ln_puchIso);
    t->SetBranchAddress("ln_nhIso",      this->ln_nhIso);
    t->SetBranchAddress("ln_ip3d",       this->ln_ip3d);
    t->SetBranchAddress("ln_ip3dsig",    this->ln_ip3dsig);
    t->SetBranchAddress("ln_genid",      this->ln_genid);
    t->SetBranchAddress("ln_Tbits",      this->ln_Tbits);
   
    //selected jets
    t->SetBranchAddress("jn",             &this->jn);
    t->SetBranchAddress("jn_idbits",      this->jn_idbits);
    t->SetBranchAddress("jn_px",          this->jn_px);
    t->SetBranchAddress("jn_py",          this->jn_py);
    t->SetBranchAddress("jn_pz",          this->jn_pz);
    t->SetBranchAddress("jn_en",          this->jn_en);
    t->SetBranchAddress("jn_torawsf",     this->jn_torawsf);
    t->SetBranchAddress("jn_area",        this->jn_area);
    t->SetBranchAddress("jn_csv",         this->jn_csv);
    t->SetBranchAddress("jn_genflav",     this->jn_genflav);
    t->SetBranchAddress("jn_genid",       this->jn_genid);
   
    //MET
    t->SetBranchAddress("metn",       &this->metn);
    t->SetBranchAddress("met_phi",     this->met_phi);
    t->SetBranchAddress("met_pt",      this->met_pt);
    t->SetBranchAddress("met_sig",     this->met_sig);
    t->SetBranchAddress("met_sigx2",   this->met_sigx2);
    t->SetBranchAddress("met_sigxy",   this->met_sigxy);
    t->SetBranchAddress("met_sigy2",   this->met_sigy2);
    t->SetBranchAddress("met_sumet",   this->met_sumet);
    t->SetBranchAddress("met_chsumet",   this->met_chsumet);
  }

  //
  inline void create(TTree *t)
  {
    //event info
    t->Branch("run",      &this->run,      "run/I");
    t->Branch("lumi",     &this->lumi,     "lumi/I");
    t->Branch("event",    &this->event,    "event/I");

    //trigger bit
    t->Branch("tn",         &this->tn,        "tn/I");
    t->Branch("tbits",      &this->t_bits,    "tbits/I");
    t->Branch("tprescale",  this->t_prescale, "tprescale[tn]/I");
    
    //pileup related observables
    t->Branch("nvtx",     &this->nvtx,       "nvtx/I");
    t->Branch("instLumi", &this->instLumi,   "instLumi/F");
    t->Branch("rho",      &this->rho,        "rho/F");
    
    //generator level info
    t->Branch("ngenITpu",    &this->ngenITpu,    "ngenITpu/I");
    t->Branch("ngenOOTpu",   &this->ngenOOTpu,   "ngenOOTpu/I");
    t->Branch("ngenOOTpum1", &this->ngenOOTpum1, "ngenOOTpum1/I");
    t->Branch("ngenTruepu",  &this->ngenTruepu,  "ngenTruepu/I");
    t->Branch("pthat",       &this->pthat,       "pthat/F");
    t->Branch("genWeight",   &this->genWeight,   "genWeight/F");
    t->Branch("qscale",      &this->qscale,      "qscale/F");
    t->Branch("x1",          &this->x1,          "x1/F");
    t->Branch("x2",          &this->x2,          "x2/F");
    t->Branch("id1",         &this->id1,         "id1/I");
    t->Branch("id2",         &this->id2,         "id2/I");
    t->Branch("nup",         &this->nup,         "nup/I");
    
    //mc truth
    t->Branch("mcn",         &this->mcn,         "mcn/I");
    t->Branch("mc_id",       this->mc_id,        "mc_id[mcn]/I");
    t->Branch("mc_status",   this->mc_status,    "mc_status[mcn]/I");
    t->Branch("mc_px",       this->mc_px,        "mc_px[mcn]/F");
    t->Branch("mc_py",       this->mc_py,        "mc_py[mcn]/F");
    t->Branch("mc_pz",       this->mc_pz,        "mc_pz[mcn]/F");
    t->Branch("mc_en",       this->mc_en,        "mc_en[mcn]/F");

    //selected leptons
    t->Branch("ln",            &this->ln,          "ln/I");
    t->Branch("ln_id",         this->ln_id,        "ln_id[ln]/I");
    t->Branch("ln_idbits",     this->ln_idbits,    "ln_idbits[ln]/I");
    t->Branch("ln_genid",      this->ln_genid,     "ln_genid[ln]/I");
    t->Branch("ln_Tbits",      this->ln_Tbits,     "ln_Tbits[ln]/I");
    t->Branch("ln_px",         this->ln_px,        "ln_px[ln]/F");
    t->Branch("ln_py",         this->ln_py,        "ln_py[ln]/F");
    t->Branch("ln_pz",         this->ln_pz,        "ln_pz[ln]/F");
    t->Branch("ln_en",         this->ln_en,        "ln_en[ln]/F");
    t->Branch("ln_ecalIso",    this->ln_ecalIso,   "ln_ecalIso[ln]/F");
    t->Branch("ln_hcalIso",    this->ln_hcalIso,   "ln_hcalIso[ln]/F");
    t->Branch("ln_trkIso",     this->ln_trkIso,    "ln_trkIso[ln]/F");
    t->Branch("ln_gIso",       this->ln_gIso,      "ln_gIso[ln]/F");
    t->Branch("ln_chIso",      this->ln_chIso,     "ln_chIso[ln]/F");
    t->Branch("ln_puchIso",    this->ln_puchIso,   "ln_puchIso[ln]/F");
    t->Branch("ln_nhIso",      this->ln_nhIso,     "ln_nhIso[ln]/F");
    t->Branch("ln_ip3d",       this->ln_ip3d,      "ln_ip3d[ln]/F");
    t->Branch("ln_ip3dsig",    this->ln_ip3dsig,   "ln_ip3dsig[ln]/F");

    //selected jets
    t->Branch("jn",             &this->jn,         "jn/I");
    t->Branch("jn_idbits",      this->jn_idbits,   "jn_idbits[jn]/I");
    t->Branch("jn_px",          this->jn_px,       "jn_px[jn]/F");
    t->Branch("jn_py",          this->jn_py,       "jn_py[jn]/F");
    t->Branch("jn_pz",          this->jn_pz,       "jn_pz[jn]/F");
    t->Branch("jn_en",          this->jn_en,       "jn_en[jn]/F");
    t->Branch("jn_torawsf",     this->jn_torawsf,  "jn_torawsf[jn]/F");
    t->Branch("jn_area",        this->jn_area,     "jn_area[jn]/F");
    t->Branch("jn_csv",         this->jn_csv,      "jn_csv[jn]/F");
    t->Branch("jn_genflav",     this->jn_genflav,  "jn_genflav[jn]/I");
    t->Branch("jn_genid",       this->jn_genid,    "jn_genid[jn]/I");
    t->Branch("jn_genpx",       this->jn_genpx,    "jn_genpx[jn]/F");
    t->Branch("jn_genpy",       this->jn_genpy,    "jn_genpy[jn]/F");
    t->Branch("jn_genpz",       this->jn_genpz,    "jn_genpz[jn]/F");
    t->Branch("jn_genen",       this->jn_genen,    "jn_genen[jn]/F");

    //MET
    t->Branch("metn",        &this->metn,       "metn/I");
    t->Branch("met_phi",     this->met_phi,     "met_phi[metn]/F");
    t->Branch("met_pt",      this->met_pt,      "met_pt[metn]/F");
    t->Branch("met_sig",     this->met_sig,     "met_pt[metn]/F");
    t->Branch("met_sigx2",   this->met_sigx2,   "met_sigx2[metn]/F");
    t->Branch("met_sigxy",   this->met_sigxy,   "met_sigxy[metn]/F");
    t->Branch("met_sigy2",   this->met_sigy2,   "met_sigy2[metn]/F");
    t->Branch("met_sumet",   this->met_sumet,   "met_sumet[metn]/F");
    t->Branch("met_chsumet",   this->met_chsumet,   "met_chsumet[metn]/F");

    //superclusters
    t->Branch("scn",            &this->scn,          "scn/I");
    t->Branch("scn_e",          this->scn_e,         "scn_e[scn]/F");
    t->Branch("scn_eta",        this->scn_eta,       "scn_eta[scn]/F");
    t->Branch("scn_phi",        this->scn_phi,       "scn_phi[scn]/F");
    
    m_trees.push_back(t);
  }
  

 private:
  std::vector<TTree *> m_trees;

};

#endif
