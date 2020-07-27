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
parser.add_argument("-T", "--deltaT",   dest="deltaT",   help="The time window size in us", required=False,  default=3.0,type=float)
parser.add_argument("-R", "--rate",   dest="rate",   help="The SST hit rate in MHz", required=False,  default=2.2,type=float)
parser.add_argument("-s", "--timeS",   dest="timeS",   help="The minimal time stamp for the TDC in us", required=False,  default=0.0025,type=float)

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
#modules = shipDet_conf.configure(run,ShipGeo)
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
#global_variables.modules = modules
global_variables.iEvent = 0
global_variables.deltaT = options.deltaT
global_variables.rate = options.rate
global_variables.timeS = options.timeS

fn = ROOT.TFile.Open(options.inputFile,"read")
sTree = fn.cbmsim
options.nEvents   = min(sTree.GetEntries(),options.nEvents)
piledUpf = ROOT.TFile(outFile,"recreate")
sTree.SetBranchStatus("strawtubesPoint*",0)
sTree.SetBranchStatus("MCTrack*",0)
newTree = sTree.CloneTree(0)
sTree.SetBranchStatus("strawtubesPoint*",1)
sTree.SetBranchStatus("MCTrack*",1)
#
timeStamp = global_variables.timeS # us
rateSST = global_variables.rate # MHz
fPileUp = ROOT.TF1("fPileUp","[0]*[1]*TMath::Exp(-[0]*x)",0,global_variables.deltaT)
fPileUp.SetParameter(0,rateSST)
fPileUp.SetParameter(1,timeStamp)
global_time = 0 # us
frame_time = global_variables.deltaT * 1000 # ns
unitedStrawtubesArray = ROOT.TClonesArray("strawtubesPoint")
unitedMCTrackArray = ROOT.TClonesArray("ShipMCTrack")
unitedStrawtubesBranch = newTree.Branch("strawtubesPoint",unitedStrawtubesArray,32000,-1)
unitedMCTrackBranch = newTree.Branch("MCTrack",unitedMCTrackArray,32000,-1)
bufferOfHits = ROOT.TClonesArray("strawtubesPoint",1)
bufferOfMCTracks = ROOT.TClonesArray("ShipMCTrack",1)
index = 0
index_buffer = 0
index_MC_buffer = 0
index_MCTrack = 0
remove_buffer = []
remove_MC_buffer = []
shift_TrackID_buffer = []
trids_buffer = []
used_trids_buffer = {}
used_trids_buffer_pile = {}
used_trids_buffer_remove = {}
insideTheFrameEventCounter = 0
# main loop
for global_variables.iEvent in range(0, options.nEvents):
  if global_variables.iEvent % 1000 == 0:
    print('event ', global_variables.iEvent)
  if global_time >= frame_time:
    if index_buffer > 0:
      i = 0
      #bufferOfHits.Dump()
      #bufferOfMCTracks.Dump()
      #print("size ",bufferOfHits.GetSize())
      for hit in bufferOfHits:
        if hit.GetTime() < frame_time:
          #print("I am inside !!!")
          if unitedStrawtubesArray.GetSize() == index:
            unitedStrawtubesArray.Expand(index+1000)
          unitedStrawtubesArray[index] = ROOT.strawtubesPoint(hit)
          trid = bufferOfHits[i].GetTrackID()
          #print("trid ", trid)
          if trid >= 0:
            if not trid in used_trids_buffer_remove:
              if unitedMCTrackArray.GetSize() <= index_MCTrack+1:
                unitedMCTrackArray.Expand(index_MCTrack+1000)
              unitedStrawtubesArray[index].SetTrackID(index_MCTrack)
              unitedMCTrackArray[index_MCTrack] = ROOT.ShipMCTrack(bufferOfMCTracks[trid])
              remove_MC_buffer.append(trid)
              used_trids_buffer_remove[trid] = index_MCTrack
              motherID = unitedMCTrackArray[index_MCTrack].GetMotherId()
              index_MCTrack += 1
              if motherID >= 0:
                if not motherID in used_trids_buffer_remove:
                  unitedMCTrackArray[index_MCTrack] = ROOT.ShipMCTrack(bufferOfMCTracks[motherID])
                  unitedMCTrackArray[index_MCTrack-1].SetMotherId(index_MCTrack)
                  used_trids_buffer_remove[motherID] = index_MCTrack
                  index_MCTrack += 1
                else:
                  unitedMCTrackArray[index_MCTrack-1].SetMotherId(used_trids_buffer_remove[motherID])
            else:
              unitedStrawtubesArray[index].SetTrackID(used_trids_buffer_remove[trid])
          remove_buffer.append(i)
          index += 1
        i += 1
      used_trids_buffer_remove.clear()
      #print("remove_buffer",remove_buffer)
      #print("remove_MC_buffer",remove_MC_buffer)
      if len(remove_MC_buffer) > 0:
        for j in range(bufferOfHits.GetSize()):
          trid = bufferOfHits[j].GetTrackID()
          if not j in remove_buffer and trid in remove_MC_buffer:
            shift_TrackID_buffer.append(trid)
          if not trid in trids_buffer:
            trids_buffer.append(trid)
            motherID = bufferOfMCTracks[trid].GetMotherId()
            if motherID >= 0 and not motherID in trids_buffer:
              trids_buffer.append(motherID)

        #print("shift_TrackID_buffer ",shift_TrackID_buffer)

        for i in range(bufferOfMCTracks.GetSize()):
          if not i in trids_buffer:
            remove_MC_buffer.append(i)
        remove_MC_buffer.sort()
        #print("remove_MC_buffer",remove_MC_buffer)
        for j in range(len(remove_MC_buffer)):
          if not remove_MC_buffer[j] in shift_TrackID_buffer:
            for k in range(remove_MC_buffer[j]+1,bufferOfMCTracks.GetSize()):
              motherID = bufferOfMCTracks[k].GetMotherId()
              if motherID > remove_MC_buffer[j]:
                bufferOfMCTracks[k].SetMotherId(motherID - 1)
                #print("new motherID after removing MCTrack", bufferOfMCTracks[k].GetMotherId())
            bufferOfMCTracks.RemoveAt(remove_MC_buffer[j])
            #print("MCTrack removed: ", remove_MC_buffer[j])
            for i in range(bufferOfHits.GetSize()):
              trid = bufferOfHits[i].GetTrackID()
              if trid > remove_MC_buffer[j]-j:
                bufferOfHits[i].SetTrackID(trid - 1)
                #print("new trackID after removing MCTrack", bufferOfHits[i].GetTrackID())
        bufferOfMCTracks.Compress()
        if len(remove_MC_buffer) == bufferOfMCTracks.GetSize():
          bufferOfMCTracks.Expand(bufferOfMCTracks.GetSize()-len(remove_MC_buffer)+1)
          index_MC_buffer = 0
        else:
          bufferOfMCTracks.Expand(bufferOfMCTracks.GetSize()-len(remove_MC_buffer))
          index_MC_buffer = bufferOfMCTracks.GetSize()
        del remove_MC_buffer[:]
        del shift_TrackID_buffer[:]

      if len(remove_buffer) > 0:
        for j in range(len(remove_buffer)):
          bufferOfHits.RemoveAt(remove_buffer[j])
        bufferOfHits.Compress()
        if len(remove_buffer) == bufferOfHits.GetSize():
          bufferOfHits.Expand(bufferOfHits.GetSize()-len(remove_buffer)+1)
          index_buffer = 0
        else:
          bufferOfHits.Expand(bufferOfHits.GetSize()-len(remove_buffer))
          index_buffer = bufferOfHits.GetSize()
        del remove_buffer[:]
        #print("size after removing",bufferOfHits.GetSize())
        #print("index_buffer, buffer loop",index_buffer)
    frame_time += global_variables.deltaT * 1000
    unitedStrawtubesArray.Compress()
    unitedMCTrackArray.Compress()
    newTree.Fill()
    insideTheFrameEventCounter = 0
    index = 0
    index_MCTrack = 0
  rc = sTree.GetEvent(global_variables.iEvent)
  #print("insideTheFrameEventCounter ",insideTheFrameEventCounter)
  insideTheFrameEventCounter += 1
  #sTree.strawtubesPoint.Dump()
  #sTree.MCTrack.Dump()
  for aMCPoint in sTree.strawtubesPoint:
    newPoint = ROOT.strawtubesPoint(aMCPoint)
    newPoint.SetTime(aMCPoint.GetTime() + global_time)
    trid = newPoint.GetTrackID()
    #print("PDGcode: ",aMCPoint.PdgCode())
    #print("TrackID: ",trid)
    if newPoint.GetTime() >= frame_time:
      bufferOfHits.Expand(index_buffer+1)
      #print(bufferOfHits.GetSize())
      #print("index_buffer, strawtubes point loop",index_buffer)
      bufferOfHits[index_buffer] = ROOT.strawtubesPoint(newPoint)
      if trid >= 0:
        if not trid in used_trids_buffer_pile:
          bufferOfMCTracks.Expand(index_MC_buffer+1)
          #print("index_MC_buffer",index_MC_buffer)
          bufferOfHits[index_buffer].SetTrackID(index_MC_buffer)
          #print("Hit new TrackID: ",bufferOfHits[index_buffer].GetTrackID())
          bufferOfMCTracks[index_MC_buffer] = ROOT.ShipMCTrack(sTree.MCTrack[trid])
          used_trids_buffer_pile[trid] = index_MC_buffer
          motherID = bufferOfMCTracks[index_MC_buffer].GetMotherId()
          index_MC_buffer += 1
          if motherID >= 0:
            if not motherID in used_trids_buffer_pile:
              bufferOfMCTracks.Expand(index_MC_buffer+1)
              bufferOfMCTracks[index_MC_buffer] = ROOT.ShipMCTrack(sTree.MCTrack[motherID])
              bufferOfMCTracks[index_MC_buffer-1].SetMotherId(index_MC_buffer)
              used_trids_buffer_pile[motherID] = index_MC_buffer
              index_MC_buffer += 1
            else:
              bufferOfMCTracks[index_MC_buffer-1].SetMotherId(used_trids_buffer_pile[motherID])
        else:
          bufferOfHits[index_buffer].SetTrackID(used_trids_buffer_pile[trid])
          #print("Hit new TrackID: ",bufferOfHits[index_buffer].GetTrackID())
      index_buffer += 1
    else:
      if unitedStrawtubesArray.GetSize() == index:
        unitedStrawtubesArray.Expand(index+1000)
      unitedStrawtubesArray[index] = ROOT.strawtubesPoint(newPoint)
      if trid >= 0:
        if not trid in used_trids_buffer:
          if unitedMCTrackArray.GetSize() <= index_MCTrack+1:
            unitedMCTrackArray.Expand(index_MCTrack+1000)
          unitedStrawtubesArray[index].SetTrackID(index_MCTrack) 
          unitedMCTrackArray[index_MCTrack] = ROOT.ShipMCTrack(sTree.MCTrack[trid])
          used_trids_buffer[trid] = index_MCTrack
          motherID = unitedMCTrackArray[index_MCTrack].GetMotherId()
          index_MCTrack += 1
          if motherID >= 0:
            if not motherID in used_trids_buffer:
              unitedMCTrackArray[index_MCTrack] = ROOT.ShipMCTrack(sTree.MCTrack[motherID])
              unitedMCTrackArray[index_MCTrack-1].SetMotherId(index_MCTrack)
              used_trids_buffer[motherID] = index_MCTrack
              index_MCTrack += 1
            else:
              unitedMCTrackArray[index_MCTrack-1].SetMotherId(used_trids_buffer[motherID])
        else:
          unitedStrawtubesArray[index].SetTrackID(used_trids_buffer[trid])
      index += 1
  #print("used_trids_buffer ",used_trids_buffer)
  #print("used_trids_buffer_pile ",used_trids_buffer_pile)
  used_trids_buffer.clear()
  used_trids_buffer_pile.clear()
  global_time += fPileUp.GetRandom() * 1000 # ns
 # memory monitoring
 # mem_monitor()
# end loop over events
newTree.AutoSave()
print("sTree entries: ",sTree.GetEntries())
print("newTree entries: ",newTree.GetEntries())
fn.Close()
piledUpf.Close()
