import math,ROOT
from UserCode.sm_cms_das.Tools import *

#declare the base histogram and the stack before the loop
hpt=ROOT.TH1F('leg1pt',';Transverse momentum [GeV];Events',50,0,250)
hpt.Sumw2()

#open a ROOT file to contain the histograms
fOut=ROOT.TFile('plotter.root','RECREATE')

#get all samples of interest for the analysis: both data and MC
processes=getProcesses(jsonUrl="test/wz/wz_samples.json",getMC=True,getData=True)

#loop over the different processes
for idx,p in enumerate(processes):
    files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=19  #update the luminosity to the one in data
    xsec=0
    totalSelected=0
    totalGenerated=0

    #replicate the base histogram for a new process
    hproc=hpt.Clone(p+'_leg1pt')
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

    if xsec>=0 :
        #re-scale the yields according to the cross section
        nExpected=xsec*lumi*totalSelected/totalGenerated
        nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
        #re-scale the histograms according to cross section
        hproc.Scale(xsec*lumi/totalGenerated)
    else:
        nExpected=totalSelected
        nExpectedErr=0
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)

    #create a directory with the name of the process and store the histogram
    fOut.mkdir(p).cd()
    hproc.Write()

#close the output file
fOut.Close()
