from CMSDAS.Muons.Analyzer import Analyzer
import ROOT
ROOT.gROOT.ProcessLine('.x tdrstyle.C')

analyzer = Analyzer()



def softID(muon,vertex):
    return muon.isSoftMuon(vertex)

def detIso(muon,vertex):
    return (muon.isolationR03().sumPt+muon.isolationR03().emEt+muon.isolationR03().hadEt)/muon.pt()


analyzer.addSelection('softMuonSelection',softID,20,ROOT.kRed)
analyzer.addMultipleSelection('detIso',detIso,20,ROOT.kGreen,20,0,0.5)


analyzer.run()
canvas,g=analyzer.summarize()
canvas.Update()


