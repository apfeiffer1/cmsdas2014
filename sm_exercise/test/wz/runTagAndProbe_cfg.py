import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

InputFileName = "results/SingleMu.root"
OutputFilePrefix = "efficiency-"

EfficiencyBinningSpecification = cms.PSet(  UnbinnedVariables = cms.vstring('mass'),
                                            BinnedVariables = cms.PSet( pt=cms.vdouble(20,25,40,55,500),
                                                                        eta=cms.vdouble(-2.1,-1.2,-0.8,0,0.8,1.2,2.1)
                                                                        ),
                                            BinToPDFmap = cms.vstring('pdfSplusB')
                                            )


####
# Muon -> id+iso efficiency
####
process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                                 InputFileNames = cms.vstring(InputFileName),
                                                 InputDirectoryName = cms.string("tandp"),
                                                 InputTreeName = cms.string("tandp"),
                                                 OutputFileName = cms.string(OutputFilePrefix+"MuonSelection.root"),
                                                 NumCPU = cms.uint32(1),
                                                 SaveWorkspace = cms.bool(True),
                                                 floatShapeParameters = cms.bool(True),
                                                 Variables = cms.PSet( mass = cms.vstring("Mass", "60.0", "120.0", "[GeV]"),
                                                                       pt   = cms.vstring("Transverse momentum", "0", "1000", "[GeV]"),
                                                                       eta  = cms.vstring("Pseudo-rapidity", "-2.5", "2.5", "")
                                                                       ),
                                                 Categories = cms.PSet( passLoose    = cms.vstring("passLoose",    "dummy[pass=1,fail=0]"),
                                                                        passLooseIso = cms.vstring("passLooseIso", "dummy[pass=1,fail=0]"),
                                                                        passTight    = cms.vstring("passTight",    "dummy[pass=1,fail=0]"),
                                                                        passTightIso = cms.vstring("passTightIso", "dummy[pass=1,fail=0]")
                                                                        ),
                                                 PDFs = cms.PSet( pdfSplusB = cms.vstring( 'Gaussian::signal(mass, mean[91.2,89.0,93.0], sigma[2.3,0.5,20])',
                                                                                           'RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])',
                                                                                           'RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])',
                                                                                           'efficiency[0.8,0,1]',
                                                                                           'signalFractionInPassing[0.9]'     
                                                                                           ),
                                                                  ),
                                                 Efficiencies = cms.PSet(  Loose = cms.PSet( EfficiencyBinningSpecification,
                                                                                             EfficiencyCategoryAndState = cms.vstring("passLoose","pass","passLooseIso","pass")
                                                                                             ),
                                                                           Tight = cms.PSet( EfficiencyBinningSpecification,
                                                                                             EfficiencyCategoryAndState = cms.vstring("passTight","pass","passTightIso","pass")
                                                                                             )
                                                                           )                              
                                                 )

# run the analyizer 
process.p = cms.Path( process.TagProbeFitTreeAnalyzer )
