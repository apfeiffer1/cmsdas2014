import FWCore.ParameterSet.Config as cms
process = cms.Process("TEST")

IsMC = True


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'root://lxcms00//data3/2012/HZZ_cmgTuple/synchHCP2/H125_53X_V5100.root' # V5_10_0 version
                                )
                            )

process.MuScleFit = cms.EDProducer("MuScleFitPATMuonCorrector",
                                   src = cms.InputTag("patMuonsWithTrigger"),
                                   debug = cms.bool(True),
                                   identifier = cms.string("Summer12_DR53X"),
                                   applySmearing = cms.bool(True),
                                   fakeSmearing = cms.bool(False),
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.maxEvents.input = -1

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('test.root')
)



process.p = cms.EndPath( process.MuScleFit )
process.pOut = cms.EndPath( process.out )


