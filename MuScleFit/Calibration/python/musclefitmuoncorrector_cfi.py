import FWCore.ParameterSet.Config as cms

MuScleFit = cms.EDProducer("MuScleFitPATMuonCorrector", 
                         src = cms.InputTag("patMuonsWithTrigger"), 
                         debug = cms.bool(True), 
                         identifier = cms.string("Summer12_DR53X"), 
                         applySmearing = cms.bool(True), 
                         fakeSmearing = cms.bool(False)
                         )
