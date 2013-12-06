import FWCore.ParameterSet.Config as cms
import os,sys

isMC=True
gtag="START53_V23::All"
fileList=cms.untracked.vstring(#'/store/mc/Summer12/DYToMuMu_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/8696787E-F298-E111-BA4A-001A928116BC.root'
    #'/store/mc/Summer12/WplusToMuNu_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/E84FC5D5-5198-E111-B935-003048679150.root'
    # '/store/mc/Summer12/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/FA765CC2-22A1-E111-8CB5-003048D4771E.root',
    '/store/mc/Summer12/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/FC1008B3-FEA0-E111-9E50-001A64789D70.root',
    '/store/mc/Summer12/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S8_START52_V9-v2/0000/FEC2E1A9-FCA0-E111-A929-002481E14D84.root'
    )
execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/sm_cms_das/test/runSManalyzer_cfg.py'))

