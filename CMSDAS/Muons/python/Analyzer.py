from DataFormats.FWLite import Events, Handle,Lumis
import itertools
import ROOT
import os
import json
import math
ROOT.gSystem.Load("MuScleFitCalibration")

def getFullPath(path):
    return os.environ['CMSSW_BASE']+"/src/"+path

        
class Analyzer (object):
    def __init__(self):
        self.vertexHandle  = Handle ('std::vector<reco::Vertex>')
        self.muonHandle    = Handle('std::vector<pat::Muon>')
        self.selections = {}
        self.signal=0
        self.background=0
        self.corrector = ROOT.MuScleFitCorrector(getFullPath("MuScleFit/Calibration/data/MuScleFit_2012D_DATA_53X.txt"))
#        self.corrector = ROOT.rochcor2012()
#        self.corrector = ROOT.MuScleFitCorrector(getFullPath("MuScleFit/Calibration/data/MuScleFit_2012ABC_DATA_53X.txt"))
        self.processFunc = None
        
    def readCollections(self,event):
        event.getByLabel('offlinePrimaryVertices',self.vertexHandle)
        event.getByLabel('selectedPatMuons',self.muonHandle)
        self.muons = self.muonHandle.product()
        self.vertex = self.vertexHandle.product()[0]  


    def addSelection(self,name,selection,markerStyle,markerColor):
        self.selections[name] = {'function':selection,'style':markerStyle,'color':markerColor,'signal':0.0,'background':0.0,'value':-1}

    def addMultipleSelection(self,name,selection,markerStyle,markerColor,points,min,max):
        offset = float(max-min)/points
        for i in range(0,points):
            self.selections[name+'_'+str(i*offset)] = {'function':selection,'style':markerStyle,'color':markerColor,'signal':0.0,'background':0.0,'value':i*offset}

    def analyze(self):
        for mu1,mu2 in itertools.combinations(self.muons,2):
            #signal selection
            if mu1.isTightMuon(self.vertex) and \
               mu1.pt()>20 and mu2.pt()>20 and \
               mu1.chargedHadronIso()<0.15 and \
               mu1.charge()+mu2.charge() ==0 and \
               (mu1.p4()+mu2.p4()).M()>80 and (mu1.p4()+mu2.p4()).M()<120.:
                self.signal=self.signal+1
                
                if self.processFunc is not None and mu2.isTightMuon(self.vertex) and mu2.chargedHadronIso()/mu2.pt()<0.1:
                    #correct
                    vector1 = ROOT.TLorentzVector(mu1.px(),mu1.py(),mu1.pz(),mu1.energy())
                    vector2 = ROOT.TLorentzVector(mu2.px(),mu2.py(),mu2.pz(),mu2.energy())
                    vectorC1 = ROOT.TLorentzVector(mu1.px(),mu1.py(),mu1.pz(),mu1.energy())
                    vectorC2 = ROOT.TLorentzVector(mu2.px(),mu2.py(),mu2.pz(),mu2.energy())
                    self.corrector.applyPtCorrection(vectorC1,mu1.charge())
                    self.corrector.applyPtCorrection(vectorC2,mu2.charge())
#                    self.corrector.momcor_data(vectorC1, mu1.charge(), 1);
#                    self.corrector.momcor_data(vectorC2, mu1.charge(), 1);
                    if mu1.charge()>0:
                        self.processFunc(vector1,vector2,vectorC1,vectorC2)
                    else:
                        self.processFunc(vector2,vector1,vectorC2,vectorC1)

                for name,selection in self.selections.iteritems():
                    if selection['value']<0:
                        if selection['function'](mu2,self.vertex):
                            selection['signal'] = selection['signal']+1
                    else:        
                        if selection['function'](mu2,self.vertex)<selection['value']:
                            selection['signal'] = selection['signal']+1

            if mu1.pt()>20 and mu2.pt()>20  \
               and mu1.chargedHadronIso()>5.0 and (mu1.p4()+mu2.p4()).M()<80. \
               and mu1.charge()+mu2.charge() !=0 :
                self.background=self.background+1
                for name,selection in self.selections.iteritems():
                    if selection['value']<0:
                        if selection['function'](mu2,self.vertex):
                            selection['background'] = selection['background']+1
                    else:        
                        if selection['function'](mu2,self.vertex)<selection['value']:
                            selection['background'] = selection['background']+1



    def summarize(self):
        graphs=[]
        c = ROOT.TCanvas("c")
        h = c.DrawFrame(0.0,0.0,1.0,1.0)
        h.GetXaxis().SetTitle("background fractionpassing selection")
        h.GetYaxis().SetTitle("signal fraction passing selection")
        graphs.append(h)
        for name,selection in self.selections.iteritems():
            selection['signal'] = selection['signal']/self.signal
            selection['background'] = selection['background']/self.background
            print 'Selection:'+str(name),'Signal (%):'+str( selection['signal']),'background (%):'+str( selection['background'])
            g=ROOT.TGraph()
            g.SetName(name+'g')
            g.SetPoint(0,selection['background'],selection['signal'])
            g.SetMarkerColor(selection['color'])
            g.SetMarkerStyle(selection['style'])
            g.Draw("Psame")
            graphs.append(g)
        c.Draw()
        return c,graphs

    def calculateSignificance(self,s,b):
        for name,selection in self.selections.iteritems():
            if selection['background']>0:
                value = selection['signal']*s/math.sqrt(selection['background']*b)
                print 'Selection:'+str(name),'s/sqrt(b):'+str(value)
            else:
                print 'Selection:'+str(name),'background is zero , signal is :'+str(selection['signal']*s)
                

    def run(self):
        files =[
         'root://eoscms//eos/cms//store/cmst3/group/das2014/Muons/input2.root'
            ]
        events=Events(files)
        for event in events:
            self.readCollections(event)
            self.analyze()

