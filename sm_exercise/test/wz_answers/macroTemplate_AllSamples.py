import math,ROOT
from UserCode.sm_cms_das.Tools import *

#get all samples of interest for the analysis
processes=getProcesses(jsonUrl="test/wz/wz_samples.json")

#iterate over all processes
for p in processes:
    files=getFilesForProcess(jsonUrl="test/wz/wz_samples.json",tag=p,inDir="/store/cmst3/user/psilva/CMSDAS_v2/summary/")
    lumi=1
    xsec=0
    totalSelected=0
    totalGenerated=0

    #loop over files
    for f in files:
        inF=ROOT.TFile.Open(f)
        tree=inF.Get("data/data")

        #increment the total selected and generated
        totalSelected=totalSelected+tree.GetEntries("abs(cat)==13")
        totalGenerated=totalGenerated+inF.Get("iniEvents")[0]
        xsec=inF.Get("crossSection")[0]

        inF.Close()

    #re-scale the yields according to the cross section
    nExpected=xsec*lumi*totalSelected/totalGenerated
    nExpectedErr=xsec*lumi*math.sqrt(totalSelected)/totalGenerated
    print "N(%s)=%3.3f +/- %3.3f"%(p,nExpected,nExpectedErr)
