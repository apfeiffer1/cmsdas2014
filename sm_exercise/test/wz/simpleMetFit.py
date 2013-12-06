import math,ROOT
from ROOT import *
from UserCode.sm_cms_das.Tools import *
from UserCode.sm_cms_das.PlotUtils import *

"""
Displays the results of the fit
"""
def showFitResults(w,cat) :

    gROOT.SetBatch(True)
    customROOTstyle()

    c=TCanvas('c'+cat,'c'+cat,500,500)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetLeftMargin(0.15)
 
    frame=w.var('x').frame()
    data=w.data('roohist_data_'+cat)
    data.plotOn(frame, RooFit.Name('data'))
    pdf=w.pdf('model_'+cat)
    #pdfSubSet=RooArgSet(w.pdf('pdf_other_%s'%(cat) ))
    pdfSubSet=RooArgSet(w.pdf('pdf_signal_%s'%(cat) ))
    pdf.plotOn(frame,
               RooFit.Components( pdfSubSet ),
               RooFit.MoveToBack(), RooFit.FillColor(592), RooFit.DrawOption('lf'), RooFit.Name('other') )
    pdfSubSet=RooArgSet(w.pdf('pdf_other_%s'%(cat) ), w.pdf('pdf_signal'))
    #pdfSubSet=RooArgSet(w.pdf('pdf_other_%s'%(cat) ), w.pdf('pdf_qcd'))
    pdf.plotOn(frame,
               RooFit.Components( pdfSubSet ),
               RooFit.MoveToBack(), RooFit.FillColor(17), RooFit.DrawOption('lf'), RooFit.Name('mj') )
#    pdfSubSet=RooArgSet(w.pdf('pdf_signal_%s'%(cat) ), w.pdf('pdf_other_%s'%(cat) ), w.pdf('pdf_qcd'))
#    pdf.plotOn(frame,
#               RooFit.Components( pdfSubSet ),
#               RooFit.MoveToBack(), RooFit.FillColor(614), RooFit.DrawOption('lf'), RooFit.Name('signal') )
    frame.Draw()
    frame.GetYaxis().SetTitle('Events')
    frame.GetXaxis().SetTitle('Missing transverse energy [GeV]')
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(1.5)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleOffset(0.8)

    #the CMS header
    pt = TPaveText(0.12,0.96,0.9,1.0,"brNDC")
    pt.SetBorderSize(0)
    pt.SetFillColor(0)
    pt.SetFillStyle(0)
    pt.SetTextAlign(12)
    pt.AddText("CMS preliminary, #sqrt{s}=8 TeV")
    pt.Draw()

    #region header
    regpt = TPaveText(0.8,0.96,0.9,1.0,"brNDC")
    regpt.SetBorderSize(0)
    regpt.SetFillColor(0)
    regpt.SetFillStyle(0)
    regpt.SetTextAlign(12)
    regpt.SetTextFont(42)
    if cat == '' : regpt.AddText('[SR]')
    else         : regpt.AddText('[CR]')    
    regpt.Draw()

    #only for main category
    if cat =='':
        #make a legend
        leg=TLegend(0.5,0.76,0.85,0.9)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextAlign(12)
        leg.SetLineColor(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry("data","data","lp")
        leg.AddEntry("signal","W#rightarrow l#nu","f")
        leg.AddEntry("mj","Multi-jets (data)","f")
        leg.AddEntry("other","EWK","f")
        leg.SetNColumns(2)
        leg.Draw()
        
        fitpt = TPaveText(0.6,0.76,0.9,0.45,"brNDC")
        fitpt.SetBorderSize(0)
        fitpt.SetFillColor(0)
        fitpt.SetFillStyle(0)
        fitpt.SetTextAlign(12)
        fitpt.SetTextFont(42)
        fitpt.AddText('Fit results')
        for var in [['mu','#mu'],['alpha','#alpha'],['beta','#beta']]:
            fitpt.AddText('%s=%3.3f#pm%3.3f'%(var[1],w.var(var[0]).getVal(),w.var(var[0]).getError()))
        for var in [['N_qcd_','N_{SR}(QCD)'],['N_qcd_noniso','N_{CR}(QCD)']] :
            fitpt.AddText('%s=%d#pm%d'%(var[1],w.var(var[0]).getVal(),w.var(var[0]).getError()))
        fitpt.Draw()
    else :
        c.SetLogy()
        
    c.Modified()
    c.Update()
    c.SaveAs(cat+'_fit.png')
    

"""
Auxiliary function: fill the histograms which contain data and MC in the different regions from file
"""
def importPdfsAndNormalizationsFrom(url,w):
      
    #containers for the histos
    baseHisto='leg2pt'
    baseCategories={'':None,'noniso':None}
    histos={'data':baseCategories.copy(),'signal':baseCategories.copy(),'other':baseCategories.copy()}

    #fill from file
    fIn=ROOT.TFile(url)
    for key in [key.GetName() for key in fIn.GetListOfKeys()] :
        proc='other'
        if key[0]=='W' and key.find('To')>0     : proc='signal'
        if key in ['SingleMu','SingleElectron'] : proc='data'
        for name,hist in histos[proc].iteritems():
            h=fIn.Get(key+'/'+key+'_'+baseHisto+name)
            if h is None: continue
            if hist is None:
                histos[proc][name]=h.Clone(proc+'_'+name)
                histos[proc][name].SetDirectory(0)
            else:
                hist.Add(h)
    fIn.Close()

    #import histograms to the workspace => convert to RooDataHists
    for proc,dict in histos.iteritems() :
        for name,hist in dict.iteritems() :
            rooHist=RooDataHist('roohist_'+proc+'_'+name, proc, RooArgList(w.var('x')), RooFit.Import(hist))
            if proc is 'data' :
                getattr(w,'import')( rooHist )
            else:
                w.factory('N_%s_%s[%f]'%(proc,name,hist.Integral()))
                getattr(w,'import')( RooHistPdf('pdf_'+proc+'_'+name, proc, RooArgSet(w.var('x')), rooHist ) )

"""
Instantiate a workspace and perform a combined fit
"""
def main() :

    #start a workspace to store all information
    w=RooWorkspace("w")

    #variable in which the data is counted
    w.factory('x[50,0,100]')
       
    #import histos
    importPdfsAndNormalizationsFrom(url='plotter.root',w=w)

    #define the QCD pdf expression 
    w.factory('alpha[10,0,100]')
    w.factory('beta[0,0,1]')
    w.factory("EXPR::pdf_qcd('@0*TMath::Exp(-0.5*TMath::Power(@0/(@1+@2*@0),2))',{x,alpha,beta})") 
 
    #define the signal strength
    w.factory('mu[1.0,0,10]')

    #fit individidually the signal and control regions
    nllSet=RooArgSet()
    for cat in ['','noniso'] :
    
        w.factory('N_qcd_%s[0,0,9999999999999.]'%(cat))
        w.factory("FormulaVar::N_signal_postfit_%s('@0*@1',{mu,N_signal_%s})"%(cat,cat))
        w.factory("SUM::model_%s( N_signal_postfit_%s*pdf_signal_%s, N_qcd_%s*pdf_qcd, N_other_%s*pdf_other_%s)"%(cat,cat,cat,cat,cat,cat) )
        
        pdf=w.pdf("model_%s"%(cat))
        data=w.data("roohist_data_%s"%(cat))
        pdf.fitTo(data,RooFit.Extended())

        nllSet.add( pdf.createNLL(data) )

    #now perform a combined fit by summing up the log likelihoods
    nll=RooAddition('nll','nll',nllSet)
    minuit=RooMinuit(nll)
    minuit.migrad()
    minuit.hesse() 
    getattr(w,'import')(nll)

    #save for posterity
    w.writeToFile('MetFitWorkspace.root')

    #show the plots
    for cat in ['','noniso'] :
        showFitResults(w=w,cat=cat) 


if __name__ == "__main__":
    main()
            

