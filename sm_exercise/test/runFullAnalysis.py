#!/usr/bin/env python
import os,sys
import json
import optparse
import commands

from UserCode.sm_cms_das.Tools import *

"""
Loop over the inputs and launch jobs
"""
def runFullAnalysis(myExec, inDir, jsonUrl, params, outdir, onlyTag, queue='') :

    curWorkDir=os.getcwd()
    if not os.path.isabs(myExec) :
        myExec=curWorkDir+'/'+myExec
    if myExec.find('.py')>0 :
        myExec='python %s'%myExec
    if not os.path.isabs(outdir) :
        outdir=curWorkDir+'/'+outdir
    
    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()

    for proc in procList :
    
        for desc in proc[1] :

            isData=getByLabel(desc,'isdata',False)
            data = desc['data']
            for d in data :
                dtag = getByLabel(d,'dtag','')
                if onlyTag!='' and dtag.find(onlyTag)<0 : continue
                split=getByLabel(d,'split',1)
                xsec = -1
                if not isData : xsec=getByLabel(d,'xsec',-1)

                for segment in range(0,split) :
                    eventsFile=dtag
                    if split>1:
                        eventsFile=dtag + '_' + str(segment)
                    eventsFileUrl=inDir+'/'+eventsFile+'.root'
                    if(eventsFileUrl.find('/store/')==0)  :
                        eventsFileUrl = commands.getstatusoutput('cmsPfn ' + eventsFileUrl)[1]

                    cmd='%s -i %s -o %s -x %f %s'%(myExec,eventsFileUrl,outdir,xsec,params)
                    if queue=='' :
                        os.system(cmd)
                    else:
                        script='%s/%s.sh'%(outdir,eventsFile)
                        f = open(script,'w')
                        f.write('#!/bin/bash\n')
                        f.write('cd %s/src\n'%os.environ['CMSSW_BASE'])
                        f.write('export SCRAM_ARCH=%s\n'%os.environ['SCRAM_ARCH'])
                        f.write('eval `scram r -sh`\n')
                        f.write(cmd+'\n')
                        f.close()
                        os.system('chmod u+x %s'%script)
                        os.system('bsub -q %s %s'%(queue,script))
                        
def main():

    usage = 'usage: %prog [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-e', '--exe'        ,    dest='exe'                , help='Excecutable'                            , default=None)
    parser.add_option('-i', '--in'         ,    dest='inDir'              , help='Input directory'                        , default=None)
    parser.add_option('-s', '--sub'        ,    dest='queue'              , help='Batch queue (optional)'                 , default='')
    parser.add_option('-j', '--json'       ,    dest='json'               , help='A json file with the samples to analyze', default=None)
    parser.add_option('-o', '--out'        ,    dest='outdir'             , help='Output directory'                       , default='./results')
    parser.add_option('-t', '--tag'        ,    dest='onlyTag'            , help='Only matching'                          , default='')
    parser.add_option('-p', '--pars'       ,    dest='params'             , help='Extra parameters for the job'           , default='')
    (opt, args) = parser.parse_args()

    if opt.inDir is None or opt.exe is None or opt.json is None:
        parser.print_help()
        sys.exit(1)

    os.system('mkdir -p %s'%opt.outdir)
    runFullAnalysis(myExec=opt.exe, inDir=opt.inDir, jsonUrl=opt.json, params=opt.params, outdir=opt.outdir, onlyTag=opt.onlyTag, queue=opt.queue)

if __name__ == "__main__":
    main()
