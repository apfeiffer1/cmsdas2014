/*
Run by calling either in the shell:
  > root macroTemplate.C

or directly in root:
  root[0] .x macroTemplate.C
*/

{
	// Open the file
	TFile *file = TFile::Open("root://eoscms//eos/cms/store/cmst3/user/psilva/CMSDAS_v2/summary/SingleMu_0.root", "READ");
	// Get the tree
	TTree *t = (TTree *)file->Get("data/data");

	// Make a new TCanvas
	TCanvas *c = new TCanvas("canvas", "C", 600, 600);
	c->cd();

	// Draw the variable
	t->Draw("leg1_pt", "abs(cat)==13");

	// Save the plot
	c->SaveAs("plot.png");

}