import FWCore.ParameterSet.Config as cms
import os,sys

isMC=False
gtag="FT_53_V21_AN4::All"
fileList=cms.untracked.vstring('/store/data//Run2012A/SingleMu/AOD//22Jan2013-v1/20000/FEC16A3C-8A70-E211-9AC8-1CC1DE040FE8.root')

execfile( os.path.expandvars('${CMSSW_BASE}/src/UserCode/sm_cms_das/test/runSManalyzer_cfg.py'))

