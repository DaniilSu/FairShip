#!/usr/bin/env python
inputFile = 'ship.conical.Pythia8-TGeant4.root'
geoFile   = None
debug = False
EcalDebugDraw = False
withNoStrawSmearing = None # True   for debugging purposes
nEvents    = 999999
firstEvent = 0
withHists = True
vertexing = True
dy  = None
saveDisk  = False # remove input file
pidProton = False # if true, take truth, if False fake with pion mass
realPR = ''
realPROptions=["FH", "AR", "TemplateMatching"]
withT0 = False

import resource
def mem_monitor():
 # Getting virtual memory size 
    pid = os.getpid()
    with open(os.path.join("/proc", str(pid), "status")) as f:
        lines = f.readlines()
    _vmsize = [l for l in lines if l.startswith("VmSize")][0]
    vmsize = int(_vmsize.split()[1])
    #Getting physical memory size  
    pmsize = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print "memory: virtuell = %5.2F MB  physical = %5.2F MB"%(vmsize/1.0E3,pmsize/1.0E3)

import ROOT,os,sys,getopt
import __builtin__ as builtin
import rootUtils as ut
import shipunit as u
import shipRoot_conf
from ROOT import TGraph, TFile

shipRoot_conf.configure()

try:
        opts, args = getopt.getopt(sys.argv[1:], "o:D:FHPu:n:f:g:c:hqv:sl:A:Y:i:",\
           ["ecalDebugDraw","inputFile=","geoFile=","nEvents=","noStrawSmearing","noVertexing","saveDisk","realPR=","withT0"])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter --inputFile=  --geoFile= --nEvents=  --firstEvent=,'
        print ' noStrawSmearing: no smearing of distance to wire, default on'
        print ' outputfile will have same name with _rec added'  
        print ' --realPR= defines track pattern recognition. Possible options: ',realPROptions, "if no option given, fake PR is used."
        print ' Options description:'
        print '      FH                        : Hough transform.'
        print '      AR                        : Artificial retina.'
        print '      TemplateMatching          : Tracks are searched for based on the template: track seed + hits within a window around the seed.'
        sys.exit()
for o, a in opts:
        if o in ("noVertexing",):
            vertexing = False
        if o in ("noStrawSmearing",):
            withNoStrawSmearing = True
        if o in ("--withT0",):
            withT0 = True
        if o in ("-f", "--inputFile",):
            inputFile = a
        if o in ("-g", "--geoFile",):
            geoFile = a
        if o in ("-n", "--nEvents=",):
            nEvents = int(a)
        if o in ("-Y",):
            dy = float(a)
        if o in ("--ecalDebugDraw",):
            EcalDebugDraw = True
        if o in ("--saveDisk",):
            saveDisk = True
        if o in ("--realPR",):
            realPR = a
            if not realPR in realPROptions:
              print "wrong option given for realPR,",a," should be one of ",realPROptions
              exit(1)
if EcalDebugDraw: ROOT.gSystem.Load("libASImage")

# need to figure out which geometry was used
if not dy:
  # try to extract from input file name
  tmp = inputFile.split('.')
  try:
    dy = float( tmp[1]+'.'+tmp[2] )
  except:
    dy = None
realPRoption = realPR
if realPR=='':realPRoption='No' 
print 'configured to process ',nEvents,' events from ' ,inputFile, \
      ' starting with event ',firstEvent, ' with option Yheight = ',dy,' with vertexing',vertexing,' and real pattern reco ',realPRoption
if not inputFile.find('_rec.root') < 0: 
  outFile   = inputFile
  inputFile = outFile.replace('_rec.root','.root') 
else:
  outFile = inputFile.replace('.root','_rec.root') 
# outfile should be in local directory
  tmp = outFile.split('/')
  outFile = tmp[len(tmp)-1]
  if inputFile[:7]=="root://" : os.system('xrdcp '+inputFile+' '+outFile)
  elif saveDisk: os.system('mv '+inputFile+' '+outFile)
  else :       os.system('cp '+inputFile+' '+outFile)

if not geoFile:
 tmp = inputFile.replace('ship.','geofile_full.')
 geoFile = tmp.replace('_rec','')

fgeo = ROOT.TFile.Open(geoFile)
geoMat =  ROOT.genfit.TGeoMaterialInterface()  # if only called in ShipDigiReco -> crash, reason unknown

from ShipGeoConfig import ConfigRegistry
from rootpyPickler import Unpickler
#load Shipgeo dictionary
upkl    = Unpickler(fgeo)
ShipGeo = upkl.load('ShipGeo')
ecalGeoFile = ShipGeo.ecal.File

h={}
log={}
if withHists:
 ut.bookHist(h,'distu','distance to wire',100,0.,2.)
 ut.bookHist(h,'distv','distance to wire',100,0.,2.)
 ut.bookHist(h,'disty','distance to wire',100,0.,2.)
 ut.bookHist(h,'nmeas','nr measuerements',100,0.,50.)
 ut.bookHist(h,'chi2','Chi2/DOF',100,0.,20.)
 ut.bookHist(h,'vshape','Drift Time vs distance to wire; Distance, cm; Time, ns',1000,0.,0.,10000,0.,0.)
 ut.bookHist(h,'vshape_original','Drift Time vs distance to wire; Distance, cm; Time, ns',1000,0.,0.,10000,0.,0.)
 ut.bookHist(h,'recoDist','Reco dist vs distance to wire; Reco Distance, cm; Distance to the wire, cm',250,0.,2.5,100,0.,1.)
 ut.bookHist(h,'TDC','TDC',1000,0.,0.)
 ut.bookHist(h,'driftTime','driftTime',500,0.,1500)
 ut.bookHist(h,'residuals','residuals',2000,-1,1)

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
builtin.debug    = debug
builtin.fieldMaker = fieldMaker
builtin.pidProton = pidProton
builtin.withT0 = withT0
builtin.realPR = realPR
builtin.vertexing = vertexing
builtin.ecalGeoFile = ecalGeoFile
builtin.ShipGeo = ShipGeo
builtin.modules = modules
builtin.EcalDebugDraw  = EcalDebugDraw
builtin.withNoStrawSmearing = withNoStrawSmearing
builtin.h    = h
builtin.log  = log
iEvent = 0
builtin.iEvent  = iEvent

# import reco tasks
import shipDigiReco

SHiP = shipDigiReco.ShipDigiReco(outFile,fgeo)
nEvents   = min(SHiP.sTree.GetEntries(),nEvents)

# for iEvent in range(firstEvent, nEvents):
#     SHiP.setupDriftTimeHist()
# graph = TGraph()
# ROOT.strawtubesDigi.Instance().d2w_dtRelation(h['driftTime'], graph)
# graphOUT = TFile("graph.root", "RECREATE")
# graphOUT.cd()
# graph.Write()
# graphOUT.Close()
# main loop
for iEvent in range(firstEvent, nEvents):
 if iEvent%100 == 0 or debug: print 'event ',iEvent
 rc    = SHiP.sTree.GetEvent(iEvent) 
 SHiP.digitize()
 SHiP.reconstruct()
 # memory monitoring
 # mem_monitor() 
# end loop over events
# h['vshape_original'] = ROOT.strawtubesDigi.Instance().initialVShape.Clone()
# h['residuals'] = ROOT.strawtubesDigi.Instance().residualsInStraw.Clone()
# print ROOT.strawtubesDigi.Instance().counter
# ut.bookCanvas(h,key='dist',title='dist',nx=1200,ny=600,cx=3,cy=1)
# cv=h['dist'].cd(1)
# h['disty'].Draw()
# cv=h['dist'].cd(2)
# h['distu'].Draw()
# cv=h['dist'].cd(3)
# h['distv'].Draw()
# h['dist'].Print('dist2Wire.gif')
SHiP.finish()
