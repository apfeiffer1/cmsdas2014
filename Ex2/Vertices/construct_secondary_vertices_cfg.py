# Use the tracks_and_vertices.root file as input.
import FWCore.ParameterSet.Config as cms

process = cms.Process("KSHORTS")

# Use the tracks_and_vertices.root file as input.
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:../../tracks_and_vertices.root"))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Suppress messages that are less important than ERRORs.
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

# Load part of the CMSSW reconstruction sequence to make vertexing possible.
# We'll need the CMS geometry and magnetic field to follow the true, non-helical shapes of tracks through the detector.
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "FT_53_V6_AN1::All"
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# Copy most of the vertex producer's parameters, but accept tracks with progressively more strict quality.
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")

# loose
process.SecondaryVerticesFromLooseTracks = process.generalV0Candidates.clone(
    trackRecoAlgorithm = cms.InputTag("generalTracks"),
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(False),
    trackQualities = cms.vstring("loose"),
    innerHitPosCut = cms.double(-1.),
    )

# tight
process.SecondaryVerticesFromTightTracks = process.generalV0Candidates.clone(
    trackRecoAlgorithm = cms.InputTag("generalTracks"),
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(False),
    trackQualities = cms.vstring("tight"),
    innerHitPosCut = cms.double(-1.),
    )

# highPurity
process.SecondaryVerticesFromHighPurityTracks = process.generalV0Candidates.clone(
    trackRecoAlgorithm = cms.InputTag("generalTracks"),
    selectKshorts = cms.bool(True),
    selectLambdas = cms.bool(False),
    trackQualities = cms.vstring("highPurity"),
    innerHitPosCut = cms.double(-1.),
    )

# Run all three versions of the algorithm.
process.path = cms.Path(process.SecondaryVerticesFromLooseTracks *
                        process.SecondaryVerticesFromTightTracks *
                        process.SecondaryVerticesFromHighPurityTracks)

# Writer to a new file called output.root.  Save only the new K-shorts
# and the primary vertices (for later exercises).

process.output = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring("drop *",
                                           "keep *_generalTracks_*_RECO",
                                           "keep *_globalMuons_*_RECO",
                                           "keep *_*_*_KSHORTS",
                                           "keep *_offlineBeamSpot_*_*",
                                           "keep *_offlinePrimaryVertices_*_*",
                                           "keep *_offlinePrimaryVerticesWithBS_*_*",
    ),
    fileName = cms.untracked.string("secondaryVertices.root"))
process.endpath = cms.EndPath(process.output)
