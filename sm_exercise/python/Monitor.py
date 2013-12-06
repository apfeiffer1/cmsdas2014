from ROOT import TFile, TLorentzVector, TH1F, TH2F, TGraph, TObjArray, TVectorD

"""
An histogram and object container
"""
class Monitor:
    def __init__(self,outUrl):
        self.outUrl=outUrl
        self.values={}
        self.allHistos={}
        self.extra=TObjArray()
    def addValue(self,value,name):
        self.values[name]=value
    def addToMonitor(self,h,tag):
        h.SetDirectory(0)
        name=h.GetName()
        if tag != 'base' : h.SetName(tag+'_'+name) 
        else :             h.Sumw2()
        if name in self.allHistos :
            self.allHistos[name].update({tag:h})
        else:
            self.allHistos[name]={'base':h}
    def addObject(self,obj):
        self.extra.Add(obj)
    def addHisto(self,name,title,nx,xmin,xmax) :
        h=TH1F(name,title,nx,xmin,xmax)
        self.addToMonitor(h,'base')
    def add2DHisto(self,name,title,nx,xmin,xmax,ny,ymin,ymax):
        h=TH2F(name,title,nx,xmin,xmax,ny,ymin,ymax)
        self.addToMonitor(h,'base')
    def fill(self,name,tags,valx,valy,valz=None):
        if name in self.allHistos:
            for t in tags:
                if not (t in self.allHistos[name]) :
                    self.addToMonitor(self.allHistos[name]['base'].Clone(),t)
                h=self.allHistos[name][t]
                if valz is None : h.Fill(valx,valy)
                else :            h.Fill(valx,valy,valz)
    def close(self):
        fOut=TFile.Open(self.outUrl,'RECREATE')
        outDir=fOut.mkdir('histos')
        outDir.cd()
        for key in self.allHistos:
            for tag in self.allHistos[key] :
                self.allHistos[key][tag].SetDirectory(outDir)
                self.allHistos[key][tag].Write()
        fOut.cd()
        for i in xrange(0,self.extra.GetEntriesFast()):
            obj=self.extra.At(i)
            outDir=fOut.mkdir(obj.GetName())
            outDir.cd()
            obj.SetDirectory(outDir)
            obj.Write()
            fOut.cd()
        for v in self.values :
            vVec=TVectorD(1)
            vVec[0]=self.values[v]
            vVec.Write(v)
        fOut.Close()
