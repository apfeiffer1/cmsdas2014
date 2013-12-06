import math
from ROOT import TLorentzVector
               
"""
Wrapper for a lepton candidate object
"""
class LeptonCand:
    def __init__(self,id,px,py,pz,en):
        self.id=id
        self.p4=TLorentzVector(px,py,pz,en)
        self.Tbits=0
        self.relIso=99999.
        self.passLoose=False
        self.passLooseIso=False
        self.passTight=False
        self.passTightIso=False
        self.p4Vars={'':self.p4}
        lesUnc=0.0
        if abs(id)==13 :
            lesUnc=0.01
        elif abs(id)==11:
            if math.fabs(self.p4.Eta())<1.5 :
                lesUnc=0.02
            else:
                lesUnc=0.05
        self.p4Vars['lesup']=self.p4*(1.0+lesUnc)
        self.p4Vars['lesdown']=self.p4*(1.0-lesUnc)        
    def triggerInfo(self,Tbits):
        self.Tbits=Tbits
    def selectionInfo(self,relIso,passLoose,passLooseIso,passTight,passTightIso):
        self.relIso=relIso
        self.passLoose=passLoose
        self.passLooseIso=passLooseIso
        self.passTight=passTight
        self.passTightIso=passTightIso
    def genMatch(self,genId, genPx, genPy, genPz, genEn) :
        self.genId=genId
        self.genP4=TLorentzVector(genPx,genPy,genPz,genEn)
    def getP4(self,var=''):
        try:
            return self.p4Vars[var]
        except:
            return self.p4
        
"""
Wrapper for a vector boson candidate
"""
class VectorBosonCand:
    def __init__(self, id,tag):
        self.id=id
        self.tag=tag
        self.m_legs=[]
    def addLeg(self,LeptonCand):
        if len(self.m_legs) >2 : return
        self.m_legs.append(LeptonCand)
        if len(self.m_legs)<2: return
        p1=self.m_legs[0].p4
        p2=self.m_legs[1].p4
        self.p4=p1+p2
        dphi=p1.DeltaPhi(p2)
        self.mt=math.sqrt(2*p1.Pt()*p2.Pt()*(1-math.cos(dphi)))
    def getTag(self) :
        return self.tag
    def computeMt(self,var):
        p1=self.m_legs[0].getP4(var)
        p2=self.m_legs[1].getP4(var)
        varP4=p1+p2
        dphi=p1.DeltaPhi(p2)
        return math.sqrt(2*p1.Pt()*p2.Pt()*(1-math.cos(dphi)))
        

"""
Wrapper for a jet object
"""
class JetCand:
    def __init__(self,px,py,pz,en,area,toRawSF):
        self.p4=TLorentzVector(px,py,pz,en)
        self.area=area
        self.genP4=TLorentzVector(0,0,0,0)
        self.partonId=0
        self.flavId=0
        self.corrSF=[toRawSF]
    def genMatch(self,px,py,pz,en,partonId,flavId):
        self.genP4=TLorentzVector(px,py,pz,en)
        self.partonId=partonId
        self.flavId=flavId
    def getCorrectedJet(self,type=''):
        if type==''          : return self.p4
        elif type=='raw'     : return self.p4*self.corrSF[0]
        elif type=='jesup'   : return self.p4*self.corrSF[1]
        elif type=='jesdown' : return self.p4*self.corrSF[2]
        elif type=='jerup'   : return self.p4*self.corrSF[3]
        elif type=='jerdown' : return self.p4*self.corrSF[4]
    def updateJEC(self,jecCorrector, jecUncertainty, rho, nvtx):

        eta=math.fabs(self.p4.Eta())

        #update the JEC
        newJECSF=self.corrSF[0]
        rawJet=self.getCorrectedJet('raw')
        if jecCorrector is not None:
            jecCorrector.setJetPt(rawJet.Pt())
            jecCorrector.setJetEta(eta)
            jecCorrector.setJetA(self.area)
            jecCorrector.setRho(rho);
            jecCorrector.setNPV(nvtx)
            newJECSF=jecCorrector.getCorrection()
                
        #check if resolution needs smearing (MC only)
        jerSFs=[1.0,1.0,1.0]
        genPt=self.genP4.Pt()
        if genPt>0 :
            ptSF=1.0
            ptSF_err=0.06;
            if eta<0.5 :
                ptSF=1.052;
                ptSF_err=math.sqrt(pow(0.012,2)+pow(0.5*(0.062+0.061),2))
            elif eta>=0.5 and eta<1.1:
                ptSF=1.057
                ptSF_err=math.sqrt(pow(0.012,2)+pow(0.5*(0.056+0.055),2))
            elif eta>=1.1 and eta<1.7:
                ptSF=1.096
                ptSF_err=math.sqrt(pow(0.017,2)+pow(0.5*(0.063+0.062),2))
            elif eta>=1.7 and eta<2.3:
                ptSF=1.134
                ptSF_err=math.sqrt(pow(0.035,2)+pow(0.5*(0.087+0.085),2))
            elif eta>=2.3 and eta<5.0:
                ptSF=1.288
                ptSF_err=math.sqrt(pow(0.127,2)+pow(0.5*(0.155+0.153),2))

            recPt=newJECSF*self.p4.Pt()
            jerSFs[0]=max(0.,(genPt+ptSF*(recPt-genPt)))/recPt
            jerSFs[1]=max(0.,(genPt+(ptSF+ptSF_err)*(recPt-genPt)))/recPt
            jerSFs[2]=max(0.,(genPt+(ptSF-ptSF_err)*(recPt-genPt)))/recPt

        if newJECSF==0 or jerSFs[0]==0:
            print 'Will not correct JER @ eta=%f SF_JEC=%f SF_JER=%f pT(gen)=%f pT(rec)=%f'%(eta,newJECSF,jerSFs[0],genPt,recPt)
            jerSFs[0]=1
            jerSFs[1]=1
            jerSFs[2]=1

        #JES uncertainties
        jecUncShift=0
        if jecUncertainty is not None:
            jecUncertainty.setJetPt(self.p4.Pt());
            jecUncertainty.setJetEta(eta);
            jecUncShift=math.fabs(jecUncertainty.getUncertainty(True))
        
        # the new corrected jet
        newJECSF=newJECSF*jerSFs[0]
        rawJet *= newJECSF;
        self.p4.SetPxPyPzE(rawJet.Px(),rawJet.Py(),rawJet.Pz(),rawJet.E())
        self.corrSF[0]=1./newJECSF
        self.corrSF.append(1+jecUncShift)        #+JES
        self.corrSF.append(1-jecUncShift)        #-JES
        self.corrSF.append(jerSFs[1]/jerSFs[0])  #+JER
        self.corrSF.append(jerSFs[2]/jerSFs[0])  #-JER

"""
Wrapper for a MET object
"""
class METCand:
    def __init__(self,px,py,pz,en):
        self.p4=TLorentzVector(px,py,pz,en)
        self.p4Vars={}
        self.id=0
        self.genP4=TLorentzVector(0,0,0,0)
    def genMatch(self,genPx,genPy,genPz,genEn):
        self.genId=0
        self.genP4=TLorentzVector(genPx,genPy,genPz,genEn)
    def addSumEts(self,sumet,chsumet):
        self.sumet=sumet
        self.chsumet=chsumet
    def addJetCorrections(self,jetSums):
        #make sure corrections are only in the transverse plane
        for i in xrange(0,len(jetSums)) :
            jetSums[i].SetPz(0)
            jetSums[i].SetE(jetSums[i].Pt())
        #jet/jer corrected
        self.p4=self.p4-jetSums[0]
        #uncertainties on jet/jer
        self.p4Vars['']        = self.p4
        self.p4Vars['jesup']   = self.p4-jetSums[1]
        self.p4Vars['jesdown'] = self.p4-jetSums[2]
        self.p4Vars['jerup']   = self.p4-jetSums[3]
        self.p4Vars['jerdown'] = self.p4-jetSums[4]
    def addLeptonCorrections(self,lepSums):
        #make sure corrections are only in the transverse plane
        for i in xrange(0,len(lepSums)) :
            lepSums[i].SetPz(0)
            lepSums[i].SetE(lepSums[i].Pt())
        #energy corrections
        self.p4=self.p4-lepSums[0]
        #uncertainties on energy scale
        self.p4Vars['']        = self.p4
        self.p4Vars['lesup']   = self.p4-lepSums[1]
        self.p4Vars['lesdown'] = self.p4-lepSums[2]
    def addUnclusteredCorrections(self,uncSums):
        #make sure corrections are only in the transverse plane
        for i in xrange(0,len(uncSums)) :
            uncSums[i].SetPz(0)
            uncSums[i].SetE(uncSums[i].Pt())
        #energy corrections
        self.p4=self.p4-uncSums[0]
        #uncertainties on energy scale
        self.p4Vars['']        = self.p4
        self.p4Vars['umetup']   = self.p4-uncSums[1]
        self.p4Vars['umetdown'] = self.p4-uncSums[2]
    def getP4(self,var=''):
        try:
            return self.p4Vars[var]
        except:
            return self.p4
