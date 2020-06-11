from __future__ import print_function
from __future__ import division
from argparse import ArgumentParser

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
global_variables.pidProton = pidProton
global_variables.withT0 = options.withT0
global_variables.ecalGeoFile = ecalGeoFile
global_variables.ShipGeo = ShipGeo
global_variables.modules = modules
global_variables.iEvent = 0
global_variables.deltaT = options.deltaT

fn = ROOT.TFile.Open(inputFile,"read")
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
newTree.SetBranchStatus("strawtubesPoint*",0)
newTree.SetBranchStatus("MCTrack*",0)
timeStamp = 0.0025 # us
rateSST = 2.2 # MHz
r = ROOT.TUnuran()
fPileUp = ROOT.TF1("fPileUp","[0]*[1]*TMath::Exp(-[0]*x)",0,global_variables.deltaT)
fPileUp.SetParameter(0,rateSST)
fPileUp.SetParameter(1,timeStamp)
r.Init(TUnuranDistrCont(fPileUp))
global_time = 0 # us
frame_time = global_variables.deltaT
unitedStrawtubesBranch = ROOT.TClonesArray("strawtubesPoint")
unitedMCTrackBranch = ROOT.TClonesArray("MCTrack")
bufferOfHits = []
bufferOfMCTracks = []
trackIDshift = 0
# main loop
for global_variables.iEvent in range(0, options.nEvents):
  if global_variables.iEvent % 1000 == 0:
    print('event ', global_variables.iEvent)
  if global_time < frame_time:
    rc = sTree.GetEvent(global_variables.iEvent)
    if len(bufferOfHits) > 0:
      for hit,track in bufferOfHits,bufferOfMCTracks:
        if hit.GetTime() < frame_time:
          unitedStrawtubesBranch.AddLast(hit)
          unitedMCTrackBranch.AddLast(track)
          bufferOhHits.remove(hit)
          bufferOfMCTracks.remove(track)
          trackIDshift += 1
    for aMCPoint in sTree.strawtubesPoint:
      aMCPoint.SetTime(aMCPoint.GetTime() + global_time)
      trid = aMCPoint.GetTrackID()
      if aMCPoint.GetTime() >= frame_time:
        aMCPoint.SetTrackID(len(bufferOfHits))
        bufferOfHits.append(aMCPoint)
        bufferOfMCTracks.append(sTree.MCTrack[trid])
      else:
        aMCPoint.SetTrackID(trackIDshift + aMCPoint.GetTrackID())
        unitedStrawtubesBranch.AddLast(aMCPoint)
        unitedMCTrackBranch.AddLast(sTree.MCTrack[trid])
    trackIDshift = 0
  else:
    frame_time += global_variables.deltaT
    newTree.Fill()
    newTree.strawtubesPoint = unitedStrawtubesBranch
    newTree.MCTrack = unitedMCTrackBranch
    unitedStrawtubesBranch.Clear()
    unitedMCTrackBranch.Clear()
  global_time += r.Sample()
 # memory monitoring
 # mem_monitor()
# end loop over events
newTree.AutoSave()
fn.Close()
piledUpf.Close()
