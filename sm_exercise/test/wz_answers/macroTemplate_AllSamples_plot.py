import math,ROOT
from UserCode.sm_cms_das.Tools import *

#declare the base histogram and the stack before the loop
hpt=ROOT.TH1F('leg1pt',';Transverse momentum [GeV];Events',50,0,250)
hpt.Sumw2()
mcstack=ROOT.THStack('mcstack','mcstack')

#get all samples of interest for the analysis
processes=getProcesses(jsonUrl="test/wz/wz_samples.json")

#we shall use idx to assign colors to different processes
for idx,p in enumerate(processes):
    files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=1
    xsec=0
    totalSelected=0
    totalGenerated=0

    #replicate the base histogram for a new process
    hproc=hpt.Clone(p+'_leg1pt')
    hproc.SetTitle(p)
    hproc.SetFillStyle(1001)
    hproc.SetFillColor(idx+1)
    hproc.SetDirectory(0)

    #loop over files
    for f in files:
        inF=ROOT.TFile.Open(f)
        tree=inF.Get("data/data")

        #increment the total selected and generated
        totalSelected=totalSelected+tree.GetEntries("abs(cat)==13")
        totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
        xsec=inF.Get("crossSection")[0]
         
        #project the tree to the base histogram and add it to the total expectations for this process
        hpt.SetDirectory(inF)
        tree.Draw('leg1_pt>>leg1pt','abs(cat)==13','goff')
        hproc.Add(hpt)
        hpt.SetDirectory(0)
        
        inF.Close()

    #re-scale the yields according to the cross section
    nExpected=xsec*lumi*totalSelected/totalGenerated
    nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)

    #re-scale the histograms according to cross section
    hproc.Scale(xsec*lumi/totalGenerated)
    mcstack.Add(hproc,"hist")

#display the result in a canvas
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas('c','c',500,500)
mcstack.Draw()
mcstack.GetXaxis().SetTitle(hpt.GetXaxis().GetTitle())
mcstack.GetYaxis().SetTitle(hpt.GetYaxis().GetTitle())
c.BuildLegend()
c.SaveAs('stackleg1pt.png')
