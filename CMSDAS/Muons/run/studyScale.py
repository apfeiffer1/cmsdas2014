from CMSDAS.Muons.Analyzer import Analyzer
import ROOT
ROOT.gROOT.ProcessLine('.x tdrstyle.C')

analyzer = Analyzer()

#declare plots here like this:
zDiffU = ROOT.TH1F("zDiffU","",50,-15,15)
zDiffC = ROOT.TH1F("zDiffC","",50,-15,15)




#Put the analysis code in the method below .
#The inputs are TLorentzVectors
def fillPlots(muPos,muNeg,muPosCorrected,muNegCorrected):
    zDiffU.Fill((muPos+muNeg).M()-91.1876)
    zDiffC.Fill((muPosCorrected+muNegCorrected).M()-91.1876)

analyzer.processFunc = fillPlots
analyzer.run()



