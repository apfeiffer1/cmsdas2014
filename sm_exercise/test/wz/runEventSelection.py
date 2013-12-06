#!/usr/bin/env python
import optparse
import os
import sys
import getopt
import math
import array
import random
import commands
from ROOT import gSystem,gInterpreter
from ROOT import TTree, TFile, TLorentzVector, TH1F, TH2F, TGraph, TObjArray, TNtuple

#Some wrappers to help the analysis
from UserCode.sm_cms_das.Monitor import *
from UserCode.sm_cms_das.Candidates import *

#JET/MET tools
gSystem.Load("libFWCoreFWLite.so");
from ROOT import AutoLibraryLoader
AutoLibraryLoader.enable()
gSystem.Load("libCondFormatsJetMETObjects.so")
from ROOT import JetCorrectorParameters, JetCorrectionUncertainty, FactorizedJetCorrector

"""
Decodes the bits in the trigger word written in the tree
"""
def decodeTriggerWord(trigBits) :
    eFire     = (((trigBits >> 0) & 0x1)>0) or (((trigBits >> 1) & 0x1)>0)
    mFire     = (((trigBits >> 2) & 0x1)>0) or (((trigBits >> 3) & 0x1)>0)
    emFire    = (((trigBits >> 4) & 0x3)>0)
    return eFire,mFire,emFire

"""
Performs the selection of a good lepton
"""
def selectLepton(id, idbits, gIso, chIso, nhIso, puchIso, pt) :
    isLoose=False
    isTight=False
    isLooseIso=False
    isTightIso=False

    relIso=9999.
    if abs(id)==11 :
        isLoose = ((idbits >> 3) & 0x1)
        isTight = ((idbits >> 4) & 0x1)
        relIso=(chIso+nhIso+gIso)/pt
        if relIso<0.15: isLooseIso=True
        if relIso<0.1:  isTightIso=True

    if abs(id)==13:
        isLoose = ((idbits >> 8) & 0x1)
        isTight = ((idbits >> 9) & 0x1)
        relIso=(chIso+nhIso+gIso)/pt
        if relIso<0.20: isLooseIso=True
        if relIso<0.12: isTightIso=True

    return relIso, isLoose, isLooseIso, isTight, isTightIso

"""
Returns a vector boson candidate depending on the kinematics, id and isolation of its decay leg candidates
"""
def buildVcand(eFire,mFire,emFire,leptonCands,met) :

    vCand=None
    if len(leptonCands)==0 : return vCand

    tightLeptons=[]
    tightNonIsoLeptons=[]
    vetoLeptons=[]
    for l in leptonCands :
        if   l.passTight and l.passTightIso     : tightLeptons.append(l)
        elif l.passLoose and l.passLooseIso     : vetoLeptons.append(l)
        elif l.passTight and not l.passLooseIso : tightNonIsoLeptons.append(l)

    # test first hypotheses: muon channels -> require muon trigger
    # a) 2 tight muons = Z->mm
    # b) 1 tight muon and no loose lepton = W->mv
    if vCand is None and mFire :
        if len(tightLeptons)>=2 and tightLeptons[0].id*tightLeptons[1].id==-13*13 :
            vCand = VectorBosonCand(-13*13,'mumu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])
        elif len(tightLeptons)==1 and abs(tightLeptons[0].id)==13 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(tightLeptons[0].id,'mu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(met)
        elif len(tightNonIsoLeptons)==1 and abs(tightNonIsoLeptons[0].id)==13 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(100*tightNonIsoLeptons[0].id,'munoniso')
            vCand.addLeg(tightNonIsoLeptons[0])
            vCand.addLeg(met)

    # test second hypotheses: electron channels -> require electron trigger
    # a) 2 tight electrons = Z->ee
    # b) 1 tight electron and no loose lepton = W->ev
    if vCand is None and eFire:
        if len(tightLeptons)>=2 and tightLeptons[0].id*tightLeptons[1].id==-11*11 :
            vCand = VectorBosonCand(-11*11,'ee')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])
        elif len(tightLeptons)==1 and abs(tightLeptons[0].id)==11 and len(vetoLeptons)==0:
            vCand = VectorBosonCand(tightLeptons[0].id,'e')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(met)
        elif len(tightNonIsoLeptons)==1 and abs(tightNonIsoLeptons[0].id)==11 and len(vetoLeptons)==0 :
            vCand = VectorBosonCand(100*tightNonIsoLeptons[0].id,'enoniso')
            vCand.addLeg(tightNonIsoLeptons[0])
            vCand.addLeg(met)

    # test third hypothesis with the emu trigger
    # a) 1 tight electron, 1 tight muon = Z->tt
    if vCand is None and emFire:
        if len(tightLeptons)>=2 and tightLeptons[0].id*tightLeptons[1].id==-11*13 :
            vCand = VectorBosonCand(-11*13,'emu')
            vCand.addLeg(tightLeptons[0])
            vCand.addLeg(tightLeptons[1])

    return vCand

"""
Runs the event selection over the trees
"""
def selectEvents(fileName,saveProbes=False,saveSummary=False,outputDir='./',xsec=-1,correctionsMap={}):

    gSystem.ExpandPathName(fileName)
    file=TFile.Open(fileName)

    #exclusivity of triggers per PD
    eTriggersOnly  = ('SingleEle' in fileName)
    muTriggersOnly = ('SingleMu'  in fileName)

    #normalizations and corrections
    origEvents=1.0
    puWeightsGr=None
    if xsec>0 :
        origEvents=file.Get('smDataAnalyzer/cutflow').GetBinContent(1)
        if origEvents==0 :
            print '[Warning] 0 initial events ?'

        #derive pileup weights
        origPileup=file.Get('smDataAnalyzer/pileup')
        try:
            dataPileupFile=TFile.Open(correctionsMap['pu'])
            dataPileup=dataPileupFile.Get('pileup')
            normF=origPileup.Integral()/dataPileup.Integral()
            if normF>0 :
                puWeightsGr=TGraph()
                for xbin in xrange(1,origPileup.GetXaxis().GetNbins()+1) :
                    iweight=1.0
                    if origPileup.GetBinContent(xbin)>0 :
                        iweight=normF*dataPileup.GetBinContent(xbin)/origPileup.GetBinContent(xbin)
                        puWeightsGr.SetPoint( puWeightsGr.GetN(), origPileup.GetXaxis().GetBinCenter(xbin), iweight )
            dataPileupFile.Close()
        except :
            print 'No data pileup file provided or other error occurred. If you wish add -w pu,pu_file.root'

    jecCorrector=None
    jecUncertainty=None
    try:
        prefix='Data'
        if xsec>0 : prefix='MC'
        jecDir=correctionsMap['jec']

        gSystem.ExpandPathName(jecDir)
        jetCorLevels='L1FastJet'
        jetCorFiles=jecDir+'/'+prefix+'_L1FastJet_AK5PFchs.txt'
        jetCorLevels=jetCorLevels+':L2Relative'
        jetCorFiles=jetCorFiles+':'+jecDir+'/'+prefix+'_L2Relative_AK5PFchs.txt'
        jetCorLevels=jetCorLevels+':L3Absolute'
        jetCorFiles=jetCorFiles+':'+jecDir+'/'+prefix+'_L3Absolute_AK5PFchs.txt'
        #if prefix=='Data':
        #    jetCorLevels=jetCorLevels+':L2L3Residual'
        #    jetCorFiles=jetCorFiles+':'+jecDir+'/'+prefix+'_L2L3Residual_AK5PFchs.txt'
        jecCorrector=FactorizedJetCorrector(jetCorLevels,jetCorFiles)
        print 'Jet energy corrector initialized with levels ',jetCorLevels,' for ',prefix

        if prefix=='MC':
            jecUncertainty=JetCorrectionUncertainty(jecDir+"/"+prefix+"_Uncertainty_AK5PFchs.txt")
            print 'Jet uncertainty is ',jecUncertainty

    except Exception as e:
        print '[Error]',e



    tree=file.Get("smDataAnalyzer/data")
    nev = tree.GetEntries()

    outUrl=outputDir+'/'+os.path.basename(fileName)
    monitor=Monitor(outUrl)

    #same the initial normalization and cross section
    monitor.addValue(origEvents,'iniEvents')
    monitor.addValue(xsec,'crossSection')

    #some basic histograms
    monitor.addHisto('nvtx',    ';Vertices;Events',                       50,0,50)
    monitor.addHisto('nvtxraw', ';Vertices;Events',                       50,0,50)
    monitor.addHisto('vmass',   ';Mass [GeV];Events',                     50,0,250)
    monitor.addHisto('vmt',     ';Transverse mass [GeV];Events',          50,0,250)
    monitor.addHisto('vpt',     ';Boson transverse momentum [GeV];Events',50,0,250)
    monitor.addHisto('leg1pt',  ';Transverse momentum [GeV];Events',      50,0,250)
    monitor.addHisto('leg2pt',  ';Transverse momentum [GeV];Events',      50,0,250)
    monitor.addHisto('leg1iso', ';Relative isolation;Events',             50,0,0.5)
    monitor.addHisto('leg2iso', ';Relative isolation;Events',             50,0,0.5)

    #save a summary ntuple for analysis
    summaryTuple=None
    if saveSummary :
        varList='cat:weight:nvtx:njets'
        varList=varList+':v_mass:v_mt:v_pt:genv_mass:genv_pt'
        varList=varList+':leg1_pt:leg1_eta:leg1_phi:genleg1_pt'
        varList=varList+':leg2_pt:leg2_eta:leg2_phi:genleg2_pt'
        varList=varList+':sumEt:ht'
        summaryTuple=TNtuple('data','summary',varList)
        summaryTuple.SetDirectory(0)
        monitor.addObject(summaryTuple)

    #save a dedicated ntuple for Tag and Probe
    probesTuple=None
    probesId   = array.array( 'f', [ 0 ] )
    probesPt   = array.array( 'f', [ 0 ] )
    probesEta  = array.array( 'f', [ 0 ] )
    probesPhi  = array.array( 'f', [ 0 ] )
    probesNvtx = array.array( 'f', [ 0 ] )
    probesMass = array.array( 'f', [ 0 ] )
    probesIsMatched = array.array( 'i', [0] )
    probesPassLoose = array.array( 'i', [ 0 ] )
    probesPassTight = array.array( 'i', [ 0 ] )
    probesFireTrigger = array.array( 'i', [ 0 ] )
    if saveProbes :
        probesTuple=TTree('tandp','summary for tandp')
        probesTuple.Branch( 'id', probesId, 'id/F' )
        probesTuple.Branch( 'pt', probesPt, 'pt/F' )
        probesTuple.Branch( 'eta', probesEta, 'eta/F' )
        probesTuple.Branch( 'phi', probesPhi, 'phi/F' )
        probesTuple.Branch( 'nvtx', probesNvtx, 'nvtx/F' )
        probesTuple.Branch( 'mass', probesMass, 'mass/F' )
        probesTuple.Branch( 'isMatched', probesIsMatched, 'isMatched/I' )
        probesTuple.Branch( 'passLoose', probesPassLoose, 'passLoose/I' )
        probesTuple.Branch( 'passTight', probesPassTight, 'passTight/I' )
        probesTuple.Branch( 'fireTrigger', probesFireTrigger, 'fireTrigger/I' )
        probesTuple.SetDirectory(0)
        monitor.addObject(probesTuple)

    #
    # LOOP OVER THE EVENTS
    #
    for iev in xrange(0,nev):
        tree.GetEntry(iev)
        if iev%10000 == 0 :
            sys.stdout.write("\r[ %d/100 ] completed" %(100.*iev/nev))
            sys.stdout.flush()

        #check mc truth (select V bosons from the hard process
        genBosonP4=TLorentzVector(0,0,0,0)
        genNeutP4=TLorentzVector(0,0,0,0)
        for g in xrange(0,tree.mcn):
            if tree.mc_status[g]!=3 : continue
            genP4=TLorentzVector(tree.mc_px[g],tree.mc_py[g],tree.mc_pz[g],tree.mc_en[g])
            if abs(tree.mc_id[g])==12 or abs(tree.mc_id[g])==14 or abs(tree.mc_id[g])==14 : genNeutP4=genNeutP4+genP4
            if abs(tree.mc_id[g])!=23 and abs(tree.mc_id[g])!=24 : continue
            genBosonP4=genP4
        
        #get triggers that fired
        eFire,mFire,emFire=decodeTriggerWord(tree.tbits)
        if eTriggersOnly :  mFire=False
        if muTriggersOnly : eFire=False

        #select the leptons
        leptonCands=[]
        validTags=[]
        lepSums=[TLorentzVector(0,0,0,0)]*3
        lepFlux=TLorentzVector(0,0,0,0)
        for l in xrange(0,tree.ln) :
            lep=LeptonCand(tree.ln_id[l],tree.ln_px[l],tree.ln_py[l],tree.ln_pz[l],tree.ln_en[l])
            if lep.p4.Pt()<20 : continue
            if abs(tree.ln_id[l])==11 :
                if math.fabs(lep.p4.Eta())>2.5 : continue
                if math.fabs(lep.p4.Eta())>1.4442 and math.fabs(lep.p4.Eta())<1.566 : continue
            if abs(tree.ln_id[l])==13 :
                if math.fabs(lep.p4.Eta())>2.1 : continue
            relIso, isLoose, isLooseIso, isTight, isTightIso = selectLepton(tree.ln_id[l],tree.ln_idbits[l],tree.ln_gIso[l],tree.ln_chIso[l],tree.ln_nhIso[l],tree.ln_puchIso[l],lep.p4.Pt())
            lep.selectionInfo(relIso,isLoose, isLooseIso, isTight, isTightIso)
            lep.triggerInfo(tree.ln_Tbits[l])

            #check the generator level information
            genMatchIdx=tree.ln_genid[l]
            if genMatchIdx < tree.mcn :
                lep.genMatch(tree.mc_id[genMatchIdx],tree.mc_px[genMatchIdx],tree.mc_py[genMatchIdx],tree.mc_pz[genMatchIdx],tree.mc_en[genMatchIdx])
            else :
                lep.genMatch(0,0,0,0,0)
            leptonCands.append(lep)

            if not saveProbes: continue
            if not isTight or not isTightIso or lep.Tbits==0 : continue
            if abs(lep.id)==11 and not eFire: continue
            if abs(lep.id)==13 and not mFire: continue

            validTags.append( len(leptonCands)-1 )
            lepSums[1]=lepSums[1]+lep.getP4('lesup')-lep.p4
            lepSums[2]=lepSums[2]+lep.getP4('lesdown')-lep.p4
            lepFlux=lepFlux+lep.p4

        #check if probes tree should be saved
        if saveProbes and len(validTags)>0:

            # choose a random tag
            tagIdx=random.choice(validTags)
            tag=leptonCands[tagIdx]

            #find probe
            probe=None
            for l in xrange(0,len(leptonCands)) :
                if l==tagIdx: continue
                if abs(tag.id)!=abs(leptonCands[l].id) : continue
                probe=leptonCands[l]
                break

            #for electrons save superclusters if probe is not found
            matchToEle=1
            #if abs(tag.id)==11 and probe is None :
            #    matchToEle=0
            #    for sc in xrange(0,tree.scn) :
            #        sc_en=tree.scn_e[sc]
            #        sc_eta=tree.scn_eta[sc]
            #        sc_phi=tree.scn_phi[sc]
            #        sc_pt=sc_en/math.cosh(sc_eta)
            #        sc_p4=TLorentzVector(0,0,0,0)
            #        sc_p4.SetPtEtaPhiE(sc_pt,sc_eta,sc_phi,sc_en)
            #        lscp4=tag.p4+sc_p4
            #        if math.fabs(lscp4.M()-91)>30 : continue
            #        scCand=LeptonCand(tag.id,sc_p4.Px(),sc_p4.Py(),sc_p4.Pz(),sc_p4.E())
            #        scCand.selectionInfo(0,0,0,0,0)
            #        scCand.triggerInfo(0)
            #        probe=scCand
            #        break
            if abs(tag.id)==13 : matchToEle=0

            #save info
            if probe is not None:
                tpp4=tag.p4+probe.p4
                if math.fabs(tpp4.M()-91)<30 :
                    probesId[0]=probe.id
                    probesPt[0]=probe.p4.Pt()
                    probesEta[0]=probe.p4.Eta()
                    probesPhi[0]=probe.p4.Phi()
                    probesNvtx[0]=tree.nvtx
                    probesMass[0]=tpp4.M()
                    probesIsMatched[0]=(probe.genId!=0)
                    probesPassLoose[0]=(probe.passLoose and probe.passLooseIso)
                    probesPassTight[0]=(probe.passTight and probe.passTightIso)
                    probesFireTrigger[0]=(probe.Tbits>0)
                    probesTuple.Fill()

        #jets
        selJets=[]
        jetSums=[TLorentzVector(0,0,0,0)]*5
        jetFlux=TLorentzVector(0,0,0,0)
        ht=0
        for j in xrange(0,tree.jn) :
            jet=JetCand(tree.jn_px[j],tree.jn_py[j],tree.jn_pz[j],tree.jn_en[j],tree.jn_area[j],tree.jn_torawsf[j])

            #cross clean with loose isolated leptons
            overlapFound=False
            for l in leptonCands:
                if not l.passLoose or not l.passLooseIso : continue
                dR=jet.p4.DeltaR(l.p4)
                if dR>0.4 : continue
                overlapFound=True
                break
            if overlapFound: continue

            #very loose kinematics cuts
            if math.fabs(jet.p4.Eta())>4.7 or jet.p4.Pt()<10 : continue

            #save it
            jet.genMatch(tree.jn_genpx[j],tree.jn_py[j],tree.jn_pz[j],tree.jn_en[j],tree.jn_genid[j],tree.jn_genflav[j])
            jet.updateJEC(jecCorrector,jecUncertainty,tree.rho,tree.nvtx)
            selJets.append(jet)

            #account for all the corrections you have applied
            jetSums[0]=jetSums[0] + jet.getCorrectedJet()          - jet.getCorrectedJet('raw')
            jetSums[1]=jetSums[1] + jet.getCorrectedJet('jesup')   - jet.getCorrectedJet()
            jetSums[2]=jetSums[2] + jet.getCorrectedJet('jesdown') - jet.getCorrectedJet()
            jetSums[3]=jetSums[3] + jet.getCorrectedJet('jerup')   - jet.getCorrectedJet()
            jetSums[4]=jetSums[4] + jet.getCorrectedJet('jerdown') - jet.getCorrectedJet()
            jetFlux=jetFlux+jet.p4
            ht=ht+jet.p4.Pt()

        # met
        metCand=METCand(tree.met_pt[0]*math.cos(tree.met_phi[0]),tree.met_pt[0]*math.sin(tree.met_phi[0]),0,tree.met_pt[0])
        metCand.genMatch(genNeutP4.Px(),genNeutP4.Py(),genNeutP4.Pz(),genNeutP4.E())
        metCand.addSumEts(tree.met_sumet[0], tree.met_chsumet[0])
        metCand.addJetCorrections(jetSums)
        metCand.addLeptonCorrections(lepSums)
        unclFlux=-(metCand.p4+lepFlux+jetFlux)
        unclSums=[TLorentzVector(0,0,0,0),unclFlux*0.10,unclFlux*(-0.10)]
        metCand.addUnclusteredCorrections(unclSums)

        #build the candidate
        vCand=buildVcand(eFire,mFire,emFire,leptonCands,metCand)
        if vCand is None : continue

        #prepare to save
        weight=1.0
        if puWeightsGr is not None:
            weight=puWeightsGr.Eval(tree.ngenITpu)

        #show isolations
        for ileg in [0,1]:
            hname='leg'+str(ileg+1)+'iso'
            lid=''
            if abs(vCand.m_legs[ileg].id)==11 :   lid='e'
            elif abs(vCand.m_legs[ileg].id)==13 : lid='mu'
            else : continue
            monitor.fill(hname,[lid],vCand.m_legs[ileg].relIso,weight)

        tags=[vCand.tag]
        monitor.fill('nvtxraw',tags, tree.nvtx,               1.0)
        monitor.fill('nvtx',   tags, tree.nvtx,               weight)
        monitor.fill('vmass',  tags, vCand.p4.M(),            weight)
        monitor.fill('vpt',    tags, vCand.p4.Pt(),           weight)
        monitor.fill('leg1pt', tags, vCand.m_legs[0].p4.Pt(), weight)
        monitor.fill('leg2pt', tags, vCand.m_legs[1].p4.Pt(), weight)

        for var in ['','lesup','lesdown','jesup','jesdown','jerup','jerdown','umetup','umetdown']:
            mtVar=vCand.computeMt(var)
            monitor.fill('vmt', [vCand.tag+var], mtVar, weight)


        if saveSummary :
            values=[
                vCand.id, weight, tree.nvtx, len(selJets),
                vCand.p4.M(), vCand.mt, vCand.p4.Pt(), genBosonP4.M(), genBosonP4.Pt(),
                vCand.m_legs[0].p4.Pt(),vCand.m_legs[0].p4.Eta(),vCand.m_legs[0].p4.Phi(), vCand.m_legs[0].genP4.Pt(),
                vCand.m_legs[1].p4.Pt(),vCand.m_legs[1].p4.Eta(),vCand.m_legs[1].p4.Phi(), vCand.m_legs[1].genP4.Pt(),
                metCand.sumet,ht
                ]
            summaryTuple.Fill(array.array("f",values))

    file.Close()
    monitor.close()


def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-i', '--inputFile'  ,    dest='input'  ,     help='Name of the local input file',                                    default=None)
    parser.add_option('-p', '--saveProbes' ,    dest='saveProbes',  help='Save probes tree for tag and probe',                              default=False, action="store_true")
    parser.add_option('-t', '--saveSummary' ,   dest='saveSummary', help='Save summary with selected events',                               default=False, action="store_true")
    parser.add_option('-o', '--outputDir'  ,    dest='outputDir'  , help='Name of the local output directory (default=./)',                 default='./')
    parser.add_option('-w', '--weights',        dest='weights',     help='CSV weights to apply (e.g. pu,pu_file.root,eff,eff_file.root)',   default='')
    parser.add_option('-x', '--xsec',           dest='xsec'       , help='If given is used to normalize the final histograms (default=-1)', default=-1.0, type='float')
    (opt, args) = parser.parse_args()

    if opt.input is None :
        parser.print_help()
        sys.exit(1)

    #decode weights to apply
    correctionsMap={}
    if len(opt.weights)>0 :
        tkns=opt.weights.split(',')
        for itkn in xrange(0,len(tkns)+1):
            try :
                key=tkns[itkn]
                wgtFile=tkns[itkn+1]
                if wgtFile.find('/store/')== 0:
                    wgtFile=commands.getstatusoutput('cmsPfn ' + wgtFile)[1]
                correctionsMap[key] = wgtFile
                itkn=itkn+1
            except:
                continue

    print '[runEventSelection] analyzing %s'%opt.input
    selectEvents(fileName=opt.input, saveProbes=opt.saveProbes, saveSummary=opt.saveSummary, outputDir=opt.outputDir, xsec=opt.xsec, correctionsMap=correctionsMap)


if __name__ == "__main__":
    main()
