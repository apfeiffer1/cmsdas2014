import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi import *
from FWCore.GuiBrowsers.ConfigToolBase import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerL1RefsEventContent

"""
add all trigger matches to the leptons
"""
def addTriggerMatchingTo(process) :

    #triggers of interest
    pathTrigger_MuSingle_1  = 'path("HLT_Mu15_v*")'
    pathTrigger_MuSingle_2  = 'path("HLT_Mu15_eta2p1_v*")'
    pathTrigger_EleSingle_1 = 'path("HLT_Ele17_CaloIdL_CaloIsoVL_v*")'
    pathTrigger_EleSingle_2 = 'path("HLT_Ele22_CaloIdL_CaloIsoVL_v*")'
    pathTrigger_MuEle_1     = 'path("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")'
    pathTrigger_MuEle_2     = 'path("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*")'

    #muon trigger matching
    process.muonTriggerMatchHLTMuSingle1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                            src     = cms.InputTag( 'selectedPatMuonsPFlow' ) ,
                                                            matched = cms.InputTag( 'patTrigger' ),    # selections of trigger objects ,
                                                            matchedCuts = cms.string( pathTrigger_MuSingle_1 ),    # selection of matches ,
                                                            maxDPtRel   = cms.double( 0.5 ), 
                                                            maxDeltaR   = cms.double( 0.3 ) ,
                                                            resolveAmbiguities    = cms.bool( True ) ,
                                                            resolveByMatchQuality = cms.bool( True ) 
                                                            )
    process.muonTriggerMatchHLTMuSingle2  = process.muonTriggerMatchHLTMuSingle1.clone( matchedCuts = cms.string( pathTrigger_MuSingle_2 ) )
    process.muonTriggerMatchHLTMuEle1     = process.muonTriggerMatchHLTMuSingle1.clone( matchedCuts = cms.string( pathTrigger_MuEle_1 ) )
    process.muonTriggerMatchHLTMuEle2     = process.muonTriggerMatchHLTMuSingle1.clone( matchedCuts = cms.string( pathTrigger_MuEle_2 ) )

    process.selectedPatMuonsWithTriggerMatch = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
                                                               src     = cms.InputTag( "selectedPatMuonsPFlow" ),
                                                               matches = cms.VInputTag('muonTriggerMatchHLTMuSingle1', 'muonTriggerMatchHLTMuSingle2',
                                                                                       'muonTriggerMatchHLTMuEle1',    'muonTriggerMatchHLTMuEle2'
                                                                                       )
                                                               )
    switchOnTriggerMatchEmbedding(process,
                                  triggerMatchers = [ 'muonTriggerMatchHLTMuSingle1', 'muonTriggerMatchHLTMuSingle2',
                                                      'muonTriggerMatchHLTMuEle1',    'muonTriggerMatchHLTMuEle2'],
                                  sequence='patDefaultSequencePFlow')

    #electron trigger matching
    process.eleTriggerMatchHLTEleSingle1 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "selectedPatElectronsPFlow" ),
                                                           matched = cms.InputTag( "patTrigger"),
                                                           matchedCuts = cms.string(pathTrigger_EleSingle_1),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )
    process.eleTriggerMatchHLTEleSingle2 = process.eleTriggerMatchHLTEleSingle1.clone( matchedCuts = cms.string(pathTrigger_EleSingle_2) )
    process.eleTriggerMatchHLTMuEle1     = process.eleTriggerMatchHLTEleSingle1.clone( matchedCuts = cms.string(pathTrigger_MuEle_1) )
    process.eleTriggerMatchHLTMuEle2     = process.eleTriggerMatchHLTEleSingle1.clone( matchedCuts = cms.string(pathTrigger_MuEle_2) )
    
    switchOnTriggerMatching( process,
                             ['eleTriggerMatchHLTEleSingle1','eleTriggerMatchHLTEleSingle2',
                              'eleTriggerMatchHLTMuEle1','eleTriggerMatchHLTMuEle2'],
                             sequence ='patDefaultSequencePFlow',
                             hltProcess = '*' )
    
    process.selectedPatElectronsWithTriggerMatch = cms.EDProducer("PATTriggerMatchElectronEmbedder",
                                                                  src     = cms.InputTag("selectedPatElectronsPFlow"),
                                                                  matches = cms.VInputTag( cms.InputTag('eleTriggerMatchHLTEleSingle1'),
                                                                                           cms.InputTag('eleTriggerMatchHLTEleSingle2'),
                                                                                           cms.InputTag('eleTriggerMatchHLTMuEle1'),
                                                                                           cms.InputTag('eleTriggerMatchHLTMuEle2')
                                                                                           )
                                                                  )
    
    removeCleaningFromTriggerMatching( process, sequence='patDefaultSequencePFlow' )
    process.patTrigger.processName    = cms.string( "*" )
    process.patTrigger.onlyStandAlone = cms.bool( False )
