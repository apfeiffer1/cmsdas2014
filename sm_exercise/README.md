SM exercise for CMS DAS
=======================

Note
----

These notes are mostly meant for the developers but apply for the analysis.
If doing the exercise please consider following the instructions at

https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolZandWinclusiveExercise


Before starting
---------------

register to git at https://github.com/ 

git config --global user.github &lt;your github username&gt;

Installation
------------

export SCRAM_ARCH=slc5_amd64_gcc462

cmsrel CMSSW_5_3_12_patch2

cd CMSSW_5_3_12_patch2/src/

cmsenv

wget -q -O - --no-check-certificate https://raw.github.com/pfs/sm_cms_das/master/TAGS.txt | sh

scram b -j 9

Creating ntuples
----------------

cmsRun test/runSManalyzer_data_cfg.py

cmsRun test/runSManalyzer_mc_cfg.py

or use the crab/multicrab files under test/grid

Analysis scripts
----------------

To analyze a single file

python test/wz/runEventSelection.py -i file.root

To analyze a set of files for full analysis

python test/runFullAnalysis.py -e test/wz/runEventSelection.py -i /store/cmst3/user/psilva/CMSDAS_v2 -j test/wz/wz_samples.json -o ./results -p "-p -t -w pu,/store/cmst3/user/psilva/CMSDAS_v2/LowLumiRuns_pileup.root,jec,${CMSSW_BASE}/src/UserCode/sm_cms_das/data" -s 8nh 

To launch the computation of the PDF weights

python test/runFullAnalysis.py -e ${CMSSW_BASE}/test/${SCRAM_ARCH}/computePDFvariations -i /store/cmst3/user/psilva/CMSDAS_v2 -j test/wz/wz_samples.json -o ./results/pdf -p "-t smDataAnalyzer/data -p CT10.LHgrid" -s 1nd 

Make some control plots

python test/runPlotter.py -i results -j test/wz/wz_samples.json -l 18

Run tag and probe

cmsRun test/wz/runTagAndProbe_cfg.py 


Authors
-------
