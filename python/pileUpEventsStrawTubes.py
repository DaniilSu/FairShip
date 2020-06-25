from __future__ import print_function
from __future__ import division
from argparse import ArgumentParser
from decorators import *

import ROOT,os,sys
import global_variables
import rootUtils as ut
import shipunit as u
import shipRoot_conf

shipRoot_conf.configure()

parser = ArgumentParser()

parser.add_argument("-f", "--inputFile", dest="inputFile", help="Input file", required=True)
parser.add_argument("-n", "--nEvents",   dest="nEvents",   help="Number of events to reconstruct", required=False,  default=999999,type=int)
parser.add_argument("-g", "--geoFile",   dest="geoFile",   help="ROOT geofile", required=True)
parser.add_argument("--withT0",          dest="withT0", help="simulate arbitrary T0 and correct for it", required=False, action="store_true")
parser.add_argument("-T", "--deltaT",   dest="deltaT",   help="The time window size", required=False,  default=3.0,type=float)

options = parser.parse_args()
print('configured to process ', options.nEvents, ' events from ', options.inputFile)
if not options.inputFile.find('_piled_up.root') < 0:
  outFile   = options.inputFile
  options.inputFile = outFile.replace('_piled_up.root','.root')
else:
  outFile = options.inputFile.replace('.root','_piled_up.root')
# outfile should be in local directory
  tmp = outFile.split('/')
  outFile = tmp[len(tmp)-1]
  if options.inputFile[:7]=="root://" : os.system('xrdcp '+options.inputFile+' '+outFile)
  else :       os.system('cp '+options.inputFile+' '+outFile)

if not options.geoFile:
  tmp = options.inputFile.replace('ship.','geofile_full.')
  options.geoFile = tmp.replace('_piled_up','')

fgeo = ROOT.TFile.Open(options.geoFile)
geoMat =  ROOT.genfit.TGeoMaterialInterface()  # if only called in ShipDigiReco -> crash, reason unknown

from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
#load Shipgeo dictionary
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')
ecalGeoFile = ShipGeo.ecal.File

import shipDet_conf
run = ROOT.FairRunSim()
run.SetName("TGeant4")  # Transport engine
run.SetOutputFile(ROOT.TMemFile('output', 'recreate'))  # Output file
run.SetUserConfig("g4Config_basic.C") # geant4 transport not used, only needed for creating VMC field
rtdb = run.GetRuntimeDb()
# -----Create geometry----------------------------------------------
modules = shipDet_conf.configure(run,ShipGeo)
# run.Init()
fgeo.FAIRGeom
import geomGeant4

if hasattr(ShipGeo.Bfield,"fieldMap"):
  fieldMaker = geomGeant4.addVMCFields(ShipGeo, '', True,withVirtualMC = False)

# make global variables
global_variables.fieldMaker = fieldMaker
global_variables.withT0 = options.withT0
global_variables.ecalGeoFile = ecalGeoFile
global_variables.ShipGeo = ShipGeo
global_variables.modules = modules
global_variables.iEvent = 0
global_variables.deltaT = options.deltaT

fn = ROOT.TFile.Open(options.inputFile,"read")
sTree = fn.cbmsim
options.nEvents   = min(sTree.GetEntries(),options.nEvents)
piledUpf = ROOT.TFile(outFile,"recreate")
newTree = sTree.CloneTree(0)
#  check that all containers are present, otherwise create dummy version
dummyContainers={}
branch_class = {"vetoPoint":"vetoPoint","ShipRpcPoint":"ShipRpcPoint","TargetPoint":"TargetPoint",\
                  "strawtubesPoint":"strawtubesPoint","EcalPointLite":"ecalPoint",\
                  "TimeDetPoint":"TimeDetPoint","muonPoint":"muonPoint","UpstreamTaggerPoint":"UpstreamTaggerPoint"}
for x in branch_class:
  if not newTree.GetBranch(x):
    dummyContainers[x+"_array"] = ROOT.TClonesArray(branch_class[x])
    dummyContainers[x] = newTree.Branch(x,dummyContainers[x+"_array"],32000,-1)
    setattr(newTree,x,dummyContainers[x+"_array"])
    dummyContainers[x].Fill()
#
#newTree.SetBranchStatus("strawtubesPoint*",0)
#newTree.SetBranchStatus("MCTrack*",0)
timeStamp = 0.0025 # us
rateSST = 2.2 # MHz
fPileUp = ROOT.TF1("fPileUp","[0]*[1]*TMath::Exp(-[0]*x)",0,global_variables.deltaT)
fPileUp.SetParameter(0,rateSST)
fPileUp.SetParameter(1,timeStamp)
global_time = 0 # us
frame_time = global_variables.deltaT * 1000 # ns
unitedStrawtubesArray = ROOT.TClonesArray("strawtubesPoint")
unitedMCTrackArray = ROOT.TClonesArray("ShipMCTrack")
unitedStrawtubesBranch = newTree.Branch("strawtubesPointUnited",unitedStrawtubesArray,32000,-1)
unitedMCTrackBranch = newTree.Branch("MCTrackUnited",unitedMCTrackArray,32000,-1)
bufferOfHits = ROOT.TClonesArray("strawtubesPoint",0)
bufferOfMCTracks = ROOT.TClonesArray("ShipMCTrack",0)
trackIDshift = 0
index = 0
index_buffer = 0
index_MCTrack_buffer = 0
remove_buffer = []
used_trids_buffer = []
newTreeMainEventCounter = 0
insideTheFrameEventCounter = 0
unitedMCTrackIndex = 0
# main loop
for global_variables.iEvent in range(0, options.nEvents):
  if global_variables.iEvent % 1000 == 0:
    print('event ', global_variables.iEvent)
  if global_time >= frame_time:
    if bufferOfHits.GetSize() > 0:
      i = 0
      for hit in bufferOfHits:
        if hit.GetTime() < frame_time:
          if unitedStrawtubesArray.GetSize() == index:
            unitedStrawtubesArray.Expand(index+1000)
          unitedStrawtubesArray[index] = hit
          trid = hit.GetTrackID()
          if trid >= 0:
            hit.SetTrackID(unitedMCTrackIndex+1)
            if unitedMCTrackArray.GetSize() <= unitedMCTrackIndex+1:
              unitedMCTrackArray.Expand(unitedMCTrackIndex+1+1000)
            unitedMCTrackArray.AddAt(bufferOfMCTracks[i],unitedMCTrackIndex+1)
            unitedMCTrackIndex += 1
          remove_buffer.append(i)
          index += 1
        i += 1
      if len(remove_buffer) > 0:
        #print("size before removing",bufferOfHits.GetSize())
        #print(remove_buffer)
        #bufferOfHits.Print()
        for j in range(len(remove_buffer)):
          bufferOfHits.RemoveAt(remove_buffer[j])
          bufferOfMCTracks.RemoveAt(remove_buffer[j])
        bufferOfHits.Compress()
        bufferOfMCTracks.Compress()
        bufferOfHits.Expand(bufferOfHits.GetSize()-len(remove_buffer))
        bufferOfMCTracks.Expand(bufferOfMCTracks.GetSize()-len(remove_buffer))
        #bufferOfHits.Print()
        del remove_buffer[:]
        #print("size after removing",bufferOfHits.GetSize())
        index_buffer = bufferOfHits.GetSize()
        #print("index_buffer, buffer loop",index_buffer)
    frame_time += global_variables.deltaT * 1000
    unitedStrawtubesArray.Compress()
    #unitedMCTrackArray.Compress()
    newTree.Fill()
    newTreeMainEventCounter += 1
    insideTheFrameEventCounter = 0
    #unitedStrawtubesBranch.Fill()
    #unitedMCTrackBranch.Fill()
    index = 0
    #unitedMCTrackBranch.Clear()
  rc = sTree.GetEvent(global_variables.iEvent)
  insideTheFrameEventCounter += 1
  #sTree.strawtubesPoint.Dump()
  #sTree.MCTrack.Dump()
  for aMCPoint in sTree.strawtubesPoint:
    aMCPoint.SetTime(aMCPoint.GetTime() + global_time)
    trid = aMCPoint.GetTrackID()
    #print("PDGcode: ",aMCPoint.PdgCode())
    #print("TrackID: ",trid)
    if aMCPoint.GetTime() >= frame_time:
      bufferOfHits.Expand(bufferOfHits.GetSize()+1)
      bufferOfMCTracks.Expand(bufferOfMCTracks.GetSize()+1)
      #print(bufferOfHits.GetSize())
      #print("index_buffer, strawtubes point loop",index_buffer)
      bufferOfHits[index_buffer] = aMCPoint
      if trid >= 0:
        bufferOfMCTracks[index_buffer] = sTree.MCTrack[trid]
      else:
        bufferOfMCTracks[index_buffer] = ShipMCTrack()
      index_buffer += 1
    else:
      if unitedStrawtubesArray.GetSize() == index:
        unitedStrawtubesArray.Expand(index+1000)
      if trid >= 0: 
        aMCPoint.SetTrackID(insideTheFrameEventCounter * 1000 + trid)
        unitedStrawtubesArray[index] = aMCPoint
        unitedMCTrackIndex = aMCPoint.GetTrackID()
        if not unitedMCTrackIndex in used_trids_buffer: 
          unitedMCTrackIndex = aMCPoint.GetTrackID()
          if unitedMCTrackArray.GetSize() <= unitedMCTrackIndex:
            unitedMCTrackArray.Expand(unitedMCTrackIndex+1000)
          unitedMCTrackArray[unitedMCTrackIndex] = sTree.MCTrack[trid]
          used_trids_buffer.append(unitedMCTrackIndex)
      else:
        unitedStrawtubesArray[index] = aMCPoint
      index += 1
  del used_trids_buffer[:]
  global_time += fPileUp.GetRandom() * 1000
 # memory monitoring
 # mem_monitor()
# end loop over events
newTree.AutoSave()
print("sTree entries: ",sTree.GetEntries())
print("newTree entries: ",newTree.GetEntries())
fn.Close()
piledUpf.Close()
