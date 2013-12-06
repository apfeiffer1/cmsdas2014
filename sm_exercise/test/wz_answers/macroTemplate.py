#!/usr/bin/env python
import ROOT
from ROOT import TFile, TCanvas

## Run by calling in the shell:
## python macroTemplate.py

# Open the file
file = TFile.Open("root://eoscms//eos/cms/store/cmst3/user/psilva/CMSDAS_v2/summary/SingleMu_0.root", "READ")
# Get the tree
tree = file.Get("data/data")

# Make a new TCanvas
c = TCanvas("canvas", "C", 600, 600)
c.cd()

# Draw the variable
tree.Draw("leg1_pt", "abs(cat)==13")

# Save the plot
c.SaveAs("plot.png")

