import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:../../tracks_and_vertices.root'
    )
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('dxyBias.root')
                                   )

process.dxyBias = cms.EDAnalyzer('DxyBias'
)


process.p = cms.Path(process.dxyBias)
