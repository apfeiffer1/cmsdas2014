void plainStyle()
{
  gStyle->Reset();
  gStyle->SetFillColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetOptStat(1111);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetNumberContours(100);
  gStyle->SetFrameFillColor(10);
  gROOT->ForceStyle();
}

void processOneHisto(TH1* h,
                     TF1 * f,
                     TH1* compare_signal,
                     TH1* compare_background,
                     TCanvas * c,
                     const char * name,
                     int num)
{
  c->cd();
  TPad * one = 0;
  if (num == 1)
    one = new TPad(name, name, 0, 0.5, 0.5, 1);
  if (num == 2)
    one = new TPad(name, name, 0.5, 0.5, 1, 1);
  if (num == 3)
    one = new TPad(name, name, 0, 0, 0.5, 0.5);
  one->Draw();
  one->cd();

  h->Fit(f);
  compare_signal->SetBinContent(num,
                                h->GetFunction("peak")
                                ->GetParameter(0) / h->GetBinWidth(1));
  compare_background->SetBinContent(num,
                                    h->GetFunction("peak"
                                    )->GetParameter(3) * (0.56 - 0.44) / h->GetBinWidth(1));
}


void fit()
{
  plainStyle();

  TCanvas * c1 = new TCanvas("C1", "c1", 1024, 768);
  TF1 * peak = new TF1("peak", "[0]*exp(-(x - [1])**2/(2.*[2]**2))/(sqrt(6.28)*[2]) + [3]", 0.44, 0.56);

  peak->SetParameter(0, 3.);
  peak->SetParameter(1, 0.5);
  peak->SetParameter(2, 0.004);
  peak->SetParameter(3, 30.);

  TH1F * compare_signal = new TH1F("compare_signal", "", 3, 0.5, 3.5);
  TH1F * compare_background = new TH1F("compare_background", "", 3, 0.5, 3.5);

  TH1F * tmp = (TH1F*)_file0->Get("vertexproperties/Global_Loose_InvMass");
  processOneHisto(tmp, peak, compare_signal, compare_background, c1, "Loose", 1);

  tmp = (TH1F*)_file0->Get("vertexproperties/Global_Tight_InvMass");
  processOneHisto(tmp, peak, compare_signal, compare_background, c1, "Tight", 2);

  tmp = (TH1F*)_file0->Get("vertexproperties/Global_HP_InvMass");
  processOneHisto(tmp, peak, compare_signal, compare_background, c1, "HP", 3);

  compare_signal->SetStats(false);

  TLegend * tlegend = new TLegend(0.55, 0.85, 0.95, 0.95);
  tlegend->SetFillColor(kWhite);
  tlegend->AddEntry(compare_signal, "Signal yield", "l");
  tlegend->AddEntry(compare_background, "Background yield", "l");

  c1->cd();
  TPad * results = new TPad("Results", "Results", 0.5, 0., 1, 0.5);
  results->Draw();
  results->cd();
  compare_signal->SetAxisRange(0, 11000, "Y");
  compare_signal->Draw();
  compare_background->SetLineStyle(2);
  compare_background->SetLineColor(kRed);
  compare_background->Draw("same");
  tlegend->Draw();
  c1->cd();
}
