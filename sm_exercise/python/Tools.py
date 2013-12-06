import json
import commands

"""
Gets the value of a given item
(if not available a default value is returned)
"""
def getByLabel(desc,key,defaultVal=None) :
    try :
        return desc[key]
    except KeyError:
        return defaultVal

"""
"""
def getProcessesTree(inDir,jsonUrl,tag=None):
    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()
    for proc in procList :
        for desc in proc[1] :
            data = desc['data']
            for d in data :
                dtag = getByLabel(d,'dtag','')
                if tag is not None:
                    if dtag!=tag :
                        continue
                dtagFiles=[]
                split=getByLabel(d,'split',1)
                for segment in range(0,split) :
                    eventsFile=dtag
                    if split>1:
                        eventsFile=dtag + '_' + str(segment)
                    eventsFileUrl=inDir+'/'+eventsFile+'.root'
                    if(eventsFileUrl.find('/store/')==0)  :
                        eventsFileUrl = commands.getstatusoutput('cmsPfn ' + eventsFileUrl)[1]
                    dtagFiles.append(eventsFileUrl)
                d['files']=dtagFiles
    return procList
    

"""
Parse the json file looking for a given process and build the list of files to process
"""
def getFilesForProcess(jsonUrl, tag, inDir):
    procList=getProcessesTree(inDir,jsonUrl,tag)
    for proc in procList :
        for desc in proc[1] :
            data = desc['data']
            for d in data :
                dtag = getByLabel(d,'dtag','')
                if dtag!=tag : continue
                return getByLabel(d,'files')
    return None

"""
Parse the json file and retrieve all processes
"""
def getProcesses(jsonUrl,getMC=True,getData=False,getDataDriven=False):
    toReturn=[]
    jsonFile = open(jsonUrl,'r')
    procList=json.load(jsonFile,encoding='utf-8').items()
    for proc in procList :
        for desc in proc[1] :
            isDataDriven=getByLabel(desc,'isdatadriven',False)
            if isDataDriven and not getDataDriven : continue
            isData=getByLabel(desc,'isdata',False)
            if isData and not getData : continue
            if not isData and not getMC : continue
            data = desc['data']
            for d in data :
                toReturn.append( getByLabel(d,'dtag') )
    return toReturn
