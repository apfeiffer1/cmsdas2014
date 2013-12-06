import FWCore.ParameterSet.Config as cms

process = cms.Process("SManalysis")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = gtag

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        ) 

#input
process.source = cms.Source("PoolSource", fileNames = fileList )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#output (we won't use it)
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('drop *', *patEventContent )
                               )


if(isMC) : process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.TFileService = cms.Service("TFileService", fileName = cms.string('SManalysis.root'))

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True

#apply a good vertex selector and filter out scraping
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
                                                   filterParams = pvSelector.clone( minNdof = cms.double(4.0),
                                                                                    maxZ = cms.double(24.0),
                                                                                    maxd0 = cms.double(2.0)
                                                                                    ),
                                                   src=cms.InputTag('offlinePrimaryVertices')
                                                   )
process.goodVertexFilter = cms.EDFilter("GoodVertexFilter",
                                        vertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
                                        minimumNDOF = cms.uint32(4),
                                        maxAbsZ = cms.double(24),
                                        maxd0 = cms.double(2)
                                        )

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                                  )

# optional MET filters : should add more? should run in tagging mode?
# cf.https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.metFilteringTaggers = cms.Sequence(process.HBHENoiseFilter)

#PF2PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *

postfix = "PFlow"
jetAlgo="AK5"
jecLevels=['L1FastJet', 'L2Relative', 'L3Absolute']
if(not isMC): jecLevels.append('L2L3Residual')
usePF2PAT(process,
          runPF2PAT=True,
          jetAlgo=jetAlgo,
          runOnMC=isMC,
          postfix=postfix,
          typeIMetCorrections=True,
          jetCorrections=('AK5PFchs', cms.vstring(jecLevels)),
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          )

#setup trigger matching
from UserCode.sm_cms_das.triggerMatching_cff import *
addTriggerMatchingTo(process)

#custom electrons
useGsfElectrons(process,postfix=postfix,dR="03")

#custom muons
process.patMuonsPFlow.pfMuonSource = cms.InputTag("pfSelectedMuonsPFlow")
process.muonMatchPFlow.src = cms.InputTag("pfSelectedMuonsPFlow")

#custom jets for CHS
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpIsoPFlow.checkClosestZVertex = cms.bool(False)
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False     # to use muon-cleaned electron collection set to True (check iso)
getattr(process,"pfNoElectron"+postfix).enable = False # to use electron-cleaned tau collection set to True (check iso)
getattr(process,"pfNoTau"+postfix).enable = False      # to use tau-cleaned jet collection set to True (check what is a tau)
getattr(process,"pfNoJet"+postfix).enable = True       # this i guess it's for photons...      

#compute rho from central pf candidates only
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsCentral = kt4PFJets.clone( rParam = cms.double(0.6),
                                            doAreaFastjet = cms.bool(True),
                                            doRhoFastjet = cms.bool(True),
                                            Rho_EtaMax = cms.double(2.5),
                                            Ghost_EtaMax = cms.double(2.5) )

# cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag( cms.InputTag('pfMETcorrType0'),
                                                                 cms.InputTag('pfJetMETcorr', 'type1')
                                                                 )

# superclusters
process.superClusterMerger =  cms.EDProducer("EgammaSuperClusterMerger",
                                             src = cms.VInputTag(cms.InputTag('correctedHybridSuperClusters'),
                                                                 cms.InputTag('correctedMulti5x5SuperClustersWithPreshower'))
                                             )

process.smDataAnalyzer = cms.EDAnalyzer( "SMDataAnalyzer",
                                         cfg=cms.PSet( triggerSource = cms.InputTag("TriggerResults::HLT"),
                                                       triggerPaths = cms.vstring('HLT_Ele22_CaloIdL_CaloIsoVL','HLT_Ele17_CaloIdL_CaloIsoVL',
                                                                                  'HLT_Mu15_eta2p1','HLT_Mu15',
                                                                                  'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL','HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL'
                                                                                  ),
                                                       genSource       = cms.InputTag("genParticles"),
                                                       vtxSource       = cms.InputTag("goodOfflinePrimaryVertices"),
                                                       trkSource       = cms.InputTag("generalTracks"),
                                                       beamSpotSource  = cms.InputTag("offlineBeamSpot"),
                                                       rhoSource       = cms.InputTag("kt6PFJets:rho"),
                                                       muonSource      = cms.InputTag("selectedPatMuonsWithTriggerMatch"),
                                                       electronSource  = cms.InputTag("selectedPatElectronsWithTriggerMatch"),
                                                       conversionSource= cms.InputTag("allConversions"),
                                                       jetSource       = cms.InputTag("selectedPatJetsPFlow"),
                                                       metSource       = cms.VInputTag("pfMETPFlow","pfMet","pfType1CorrectedMet","pfType1p2CorrectedMet"),
                                                       scSource        = cms.InputTag("superClusterMerger")
                                                       )
                                         )
#counters for specific filters
process.startCounter = cms.EDProducer("EventCountProducer")
process.p = cms.Path( process.startCounter
                      *process.noscraping
                      *process.goodOfflinePrimaryVertices
                      *process.goodVertexFilter
                      *process.metFilteringTaggers
                      *getattr(process,"patPF2PATSequence"+postfix)
                      *process.kt6PFJetsCentral
                      *process.type0PFMEtCorrection*process.producePFMETCorrections
                      *process.selectedPatElectronsWithTriggerMatch
                      *process.selectedPatMuonsWithTriggerMatch
                      *process.superClusterMerger
                      *process.smDataAnalyzer
                      )
