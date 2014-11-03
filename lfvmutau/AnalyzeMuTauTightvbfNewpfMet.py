'''

Run LFV H->MuTau analysis in the mu+tau channel.

Authors: Maria Cepeda, Aaron Levine, Evan K. Friis, UW

'''

import MuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT
import math
import array
#import argparse
#parser = argparse.ArgumentParser()
#args = parser.parse_args()
#print args
optimized = bool('true' in os.environ['optimized'])
preselection = bool('preselection' in os.environ['selection'])
intermediate = bool('intermediate' in os.environ['selection'])
newTauID = bool('newTaus' in os.environ['selection'])
if preselection == True or intermediate == True:
	optimized=False
Isiso = bool('true' in os.environ['iso'])
wjets = bool('true' in os.environ['wjets'])
otherJets = bool('true' in os.environ['otherJets'])#select non Ztt events
Twomu = False
Twojets = True
if preselection == False:
	Twojets = True
vbfMassCut500 = True
#if vbfMassCut500 == True:
#	vbfMassCutstr = "500 GeV"
#else:
#	vbfMassCutstr = "400 GeV"
is7TeV = bool('7TeV' in os.environ['jobid'])
isPU1signal = bool ('false' in os.environ['PU'])
isEmbed = bool ('true' in os.environ['embed'])#use EmbedPt Weight
isZtt = bool ('true' in os.environ['Ztt']) #add Ztt variable
isZjetsM50 = bool('true' in os.environ['ZjetsM50']) #zjets MC
treeSkim = bool('true' in os.environ['treeSkim']) ##useing skimmed trees
print "Is ZjetsM50:", isZjetsM50
print "Is OtherZJets:", otherJets
print "Is wjets", wjets
print "Is 7TeV:", is7TeV
print "Is PU1 signal:", isPU1signal
print "Is Embed: ", isEmbed
print "Is Ztt: ", isZtt
print "Is ZjetsM50: ", isZjetsM50
print "Preselection: " + str(preselection)
print "Two Muon Selection: " + str(Twomu)
print "Is Isolation applied: " + str(Isiso)
#print "vbfMassCut is: " + vbfMassCutstr
print "Two Jets required for vbf Preselection: " + str(Twojets)
print "Is Optimized: " + str(optimized)
###### Because I need to add a bunch more branches to the ntuple...
from math import sqrt, pi


def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI

def MH(row):
	taupx = row.tPt*math.cos(row.tPhi)
        taupy = row.tPt*math.sin(row.tPhi)
        taupz = row.tPt*math.sinh(row.tEta)
	tauP = math.sqrt(row.tPt*row.tPt + taupz*taupz)
        tauMass = row.tMass
	tauE = math.sqrt(tauMass*tauMass+tauP*tauP)

	metpx = row.pfMetEt*math.cos(row.pfMetPhi)
        metpy = row.pfMetEt*math.sin(row.pfMetPhi)
	metpz = METz(row)
	#metpz = -50
	metE = math.sqrt(row.pfMetEt*row.pfMetEt+metpz*metpz)

        mupx = row.mPt*math.cos(row.mPhi)
        mupy = row.mPt*math.sin(row.mPhi)
        mupz = row.mPt*math.sinh(row.mEta)
	muE = math.sqrt(row.mPt*row.mPt+mupz*mupz)

	Px = taupx+metpx+mupx
	Py = taupy+metpy+mupy
	Pz = taupz+metpz+mupz
	
	Psquare = Px*Px+Py*Py+Pz*Pz
	Esquare = tauE*tauE+metE*metE+muE*muE + 2*tauE*metE+2*tauE*muE+2*metE*muE
	higgsMass = math.sqrt(Esquare-Psquare)

	Pxtau = taupx+metpx
	Pytau = taupy + metpy
	Pztau = taupz+metpz
	Psquaretau = Pxtau*Pxtau+Pytau*Pytau+Pztau*Pztau
	Esquaretau = tauE*tauE+metE*metE+2*tauE*metE
	tauInvMass = math.sqrt(Esquaretau-Psquaretau) 
	#print "higgsMass: " + str(higgsMass)
	#print "tauInvMass: " + str(tauInvMass)
	return higgsMass
	

def METz(row):
	tauInvMass = 1.778
	tauMass = row.tMass
        taupx = row.tPt*math.cos(row.tPhi)
        taupy = row.tPt*math.sin(row.tPhi)
	taupz = row.tPt*math.sinh(row.tEta)
	taupt = row.tPt
        metpx = row.pfMetEt*math.cos(row.pfMetPhi)
        metpy = row.pfMetEt*math.sin(row.pfMetPhi)
	
	M = tauInvMass*tauInvMass
	T = math.sqrt(tauMass*tauMass+taupt*taupt+taupz*taupz)
	F = row.pfMetEt
	K = taupt*taupt+taupz*taupz+2*(taupx*metpx+taupy*metpy)
	Y = 2*taupz
 	
	#print "M: " + str(M)
	#print "T: " + str(T)
	#print "F: " + str(F)
	#print "K: " + str(K)
	#print "Y: " + str(Y)
	#print "Discriminant: " + str(-4*F*F*T*T*T*T+F*F*T*T*Y*Y + K*K*T*T + 2*K*M*T*T-2*K*T*T*T*T+M*M*T*T-2*M*T*T*T*T+T*T*T*T*T*T)
 	if (-4*F*F*T*T*T*T+F*F*T*T*Y*Y + K*K*T*T + 2*K*M*T*T-2*K*T*T*T*T+M*M*T*T-2*M*T*T*T*T+T*T*T*T*T*T) < 0:
		Z = (K*Y+M*Y-T*T*Y)/(4*T*T-Y*Y)
	else:

		Z = (2*math.sqrt(-4*F*F*T*T*T*T+F*F*T*T*Y*Y + K*K*T*T + 2*K*M*T*T-2*K*T*T*T*T+M*M*T*T-2*M*T*T*T*T+T*T*T*T*T*T)+K*Y+M*Y-T*T*Y)/(4*T*T-Y*Y)
		#print "Tau Eta for real MEz: " + str(row.tEta)
		#print "Tau Pt for real MEz: " + str(row.tPt)
		#print "Tau Phi for real MEz: " + str(row.tPhi)
		#print "MET Phi for real MEz: " + str(row.pfMetPhi)
		#print "MET for real MEz: " + str(row.pfMetEt)
	
	a = 4*(tauMass*tauMass+taupt*taupt)
	b = -4*(2*taupx*metpx*taupz + 2*taupy*metpy*taupz+taupz*(tauInvMass*tauInvMass-tauMass*tauMass))
	c = tauInvMass*tauInvMass*tauInvMass*tauInvMass+tauMass*tauMass*tauMass*tauMass+4*(taupx*taupx*metpx*metpx+taupy*taupy*metpy*metpy+2*taupx*taupy*metpx*metpy)-2*tauInvMass*tauInvMass*tauMass*tauMass+4*(taupx*metpx+taupy*metpy)*(tauInvMass*tauInvMass-tauMass*tauMass)-4*(tauMass*tauMass+taupx*taupx+taupy*taupy+taupz*taupz)*(metpx*metpx+metpy*metpy)
	#print "discriminant" + str(b*b-4*a*c)
	#print "a" + str(a)
	#print "b" + str(b)
	#print "c" + str(c)
        metpz = Z
	#print "metpz" + str(metpz)
	return metpz

def collMass(row):
        taupx = row.tPt*math.cos(row.tPhi)
        taupy = row.tPt*math.sin(row.tPhi)
        taupt=  row.tPt

        metpx = row.pfMetEt*math.cos(row.pfMetPhi)
        metpy = row.pfMetEt*math.sin(row.pfMetPhi)
        met = row.pfMetEt

        METproj= abs(metpx*taupx+metpy*taupy)/taupt

        xth=taupt/(taupt+METproj)
        den=math.sqrt(xth)

        mass=row.m_t_Mass/den

        return mass

def fullMT(met,mupt,taupt, metphi, muphi, tauphi):
	mux=mupt*math.cos(muphi)
	muy=mupt*math.sin(muphi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        taux=taupt*math.cos(tauphi)
        tauy=taupt*math.sin(tauphi)
	full_et=met+mupt+taupt # for muon and tau I am approximating pt~et (M<<P)
	full_x=metx+mux+taux
        full_y=mety+muy+tauy
	full_mt_2 = full_et*full_et-full_x*full_x-full_y*full_y
	full_mt=0
	if (full_mt_2>0):
		full_mt= math.sqrt(full_mt_2)
	return full_mt

def fullPT(met,mupt,taupt, metphi, muphi, tauphi):
        mux=mupt*math.cos(muphi)
        muy=mupt*math.sin(muphi)
        metx=met*math.cos(metphi)
        mety=met*math.sin(metphi)
        taux=taupt*math.cos(tauphi)
        tauy=taupt*math.sin(tauphi)
        full_x=metx+mux+taux
        full_y=mety+muy+tauy
        full_pt_2 = full_x*full_x+full_y*full_y
        full_pt=0
        if (full_pt_2>0):
                full_pt= math.sqrt(full_pt_2)
        return full_pt

################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################

# Determine MC-DATA corrections

# Make PU corrector from expected data PU distribution
# PU corrections .root files from pileupCalc.py
pu_distributions = glob.glob(os.path.join(
#    'inputs', os.environ['jobid'], 'data_TauPlusX*pu.root'))
	'inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight(
    'S6' if is7TeV else 'S10', *pu_distributions)

muon_pog_PFTight_2011 = MuonPOGCorrections.make_muon_pog_PFTight_2011()
#muon_pog_PFTight_2012 = MuonPOGCorrections.make_muon_pog_PFTight_2012()
muon_pog_PFTight_2012 = MuonPOGCorrections.make_muon_pog_PFTight_2012ABCD()

muon_pog_PFRelIsoDB02_2011 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2011()
#muon_pog_PFRelIsoDB02_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012()
muon_pog_PFRelIsoDB02_2012 = MuonPOGCorrections.make_muon_pog_PFRelIsoDB012_2012ABCD()

#muon_pog_IsoMu24eta2p1_2011 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2011() //  This does not exist,  yet :-)
muon_pog_IsoMu24eta2p1_2012 = MuonPOGCorrections.make_muon_pog_IsoMu24eta2p1_2012()



# Get object ID and trigger corrector functions
def mc_corrector_2011(row):
    if row.run > 2:
        return 1
    if isPU1signal == False:
    	pu = pu_corrector(row.nTruePU)
    else:
	pu = 1
    #pu = 1
    m1id = muon_pog_PFTight_2011(row.mPt, row.mEta)
    m1iso = muon_pog_PFRelIsoDB02_2011(row.mPt, row.mEta)
#    m_trg = H2TauCorrections.correct_mueg_mu_2011(row.mPt, row.mAbsEta)
 #   m_trg = muon_pog_IsoMu24eta2p1_2011(row.mPt, row.mAbsEta)     // Future: Figure out  how to fix this ones (see comment in FSA/T&P/MuonPOGCorrections
    return pu*m1id*m1iso*m_trg

def mc_corrector_2012(row):
    if row.run > 2:
        return 1
    #pu = 1
    if isPU1signal == False:
    	pu = pu_corrector(row.nTruePU)
    else:
	pu = 1
    m1id = muon_pog_PFTight_2012(row.mPt, row.mEta)
    m1iso = muon_pog_PFRelIsoDB02_2012(row.mPt, row.mEta)
#    m_trg = H2TauCorrections.correct_mueg_mu_2012(row.mPt, row.mAbsEta)
    m_trg = muon_pog_IsoMu24eta2p1_2012(row.mPt, row.mAbsEta)
    return pu*m1id*m1iso*m_trg

# Determine which set of corrections to use
mc_corrector = mc_corrector_2011
if not is7TeV:
    mc_corrector = mc_corrector_2012


# ApplyFakeRateMethod
def getFakeRateFactor(row,fakeName='Central'):
###INEFFICIENT APPLICATION OF FAKE RATE METHOD COMMENTED OUT#############################


    #if preselection == True:
    #mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    #mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
    #if preselection == False:
    #if row.tDecayMode == 0:
    #	mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan29_DecayMode0/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    #	mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan29_DecayMode0/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
    #if row.tDecayMode == 1:
    #	mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan30_DecayMode1/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    #	mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan30_DecayMode1/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
   # if row.tDecayMode ==10:
    #	mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan30_DecayMode10/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    #	mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan30_DecayMode10/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
    #else:
    #mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan29_NoBTaus/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    #mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan29_NoBTaus/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
    #mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan31_TightPlusLooseTauIso/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    #mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Jan31_TightPlusLooseTauIso/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
     #   mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/signal/isotrue/AnalyzeMuMuTauTight/"
     #   mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/signal/isofalse/AnalyzeMuMuTauTight/"
    #mumuidiso_dir = "MuMuIdIso/"
    #mumuid_dir = "MuMuId/"
    #data_file_str = "datammt_2012.root"
    #zjets_file_str = "Zjetsmmtvbf.root"
    #zjets_file_idiso = ROOT.TFile(mumuidiso_dir+zjets_file_str)
    #zjets_file_id = ROOT.TFile(mumuid_dir+zjets_file_str)
    #data_file_idiso = ROOT.TFile(mumuidiso_dir+data_file_str)
    #data_file_id = ROOT.TFile(mumuid_dir+data_file_str)
    #histoIdIso = data_file_idiso.Get("vbf/tJetPt").Clone()
    #histoId = data_file_id.Get("vbf/tJetPt").Clone()
#### try tPt test instead of tJetPt
    #histoIdIso = data_file_idiso.Get("vbf/tPt").Clone()
    #histoId = data_file_id.Get("vbf/tPt").Clone()


    #i = 0
    #bincount = -1
    #binning = array.array('d',[])
    #while (i <= 480):
    #	if i < 70:
     #   	binning.append(i)
      #          i = i+2
       # elif i< 100:
        #        binning.append(i)
         #       i = i+5
        #elif i<120:
        #	binning.append(i)
         #       i = i+10
        #else:
        #	binning.append(i)
         #       i = i+20
        #bincount = bincount+1
    #print bincount
    #print len(binning)
    #print binning
    #histoIdIsoRebin = histoIdIso.Rebin(bincount,"histoIdIsoRebin",binning)
    #histoIdRebin = histoId.Rebin(bincount,"histoIdRebin",binning)
    #histo_ftau = histoIdIsoRebin.Clone()
    #histo_ftau.Divide(histoIdRebin)
    
    #histoIdIso.Rebin(binwidth)
    #histoId.Rebin(binwidth)
    #histoIdIso = zjets_file_idiso.Get("vbf/tJetPt").Clone()
    #histoId = zjets_file_id.Get("vbf/tJetPt").Clone()
    #if histoIdRebin.GetBinContent(histoId.FindBin(row.tJetPt)) == 0:
#	fTauIso = 0.0
 #   else:
#	fTauIso = histo_ftau.GetBinContent(histo_ftau.FindBin(row.tJetPt))
 #   	if fakeName=='Down':
#		fTauIso = fTauIso - histo_ftau.GetBinError(histo_ftau.FindBin(row.tJetPt))
#		if fTauIso <= 0:
#			fTauIso = 0.0
#	elif fakeName=='Up':
 #               fTauIso = fTauIso + histo_ftau.GetBinError(histo_ftau.FindBin(row.tJetPt))
    #print "fTauIso: " + str(fTauIso)
    #print "tauJetPt: " + str(row.tJetPt)
    #print fTauIso
    #print row.tDecayMode
  #  if fTauIso >=1:
#	fTauIso = 0.0
    ##Scale Test
    #if row.tDecayMode==0:
    #    fTauIso = 0.53
    #elif row.tDecayMode==1 or row.tDecayMode==2:
    #    fTauIso = 0.48
    #elif row.tDecayMode==10:
    #    fTauIso = 0.46
    #fTauIso = 0.52 + (0.00002172)*row.tJetPt
    fTauIso = 0.52 -0.000901418*row.tJetPt
    fakeRateFactor = fTauIso/(1.0-fTauIso)
    #print fakeRateFactor
    return fakeRateFactor

class AnalyzeMuTauTightvbfNewpfMet(MegaBase):
    if treeSkim == True:
    	tree = 'New_Tree'
	print "tree set to New_tree"	
    else:
    	tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeMuTauTightvbfNewpfMet, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuTauTree.MuTauTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):

        names=["gg0","gg1","vbf","highMtgg0","highMtgg1","highMtvbf", "hightMtvbf","lowtMtvbf","antiisomuongg0","antiisomuongg1","antiisomuonvbf","highMtantiisomuonvbf","antiisotaugg0","antiisotaugg1","antiisotauvbf","antiisotauvbfdown","antiisotauvbfup","highMtssantiisotauvbf","highMtssantiisotauvbfdown","highMtssantiisotauvbfup","antiisotauhighMtvbf","highMtssantiisomuonvbf","ssgg0","ssgg1","highMtssgg0", "highMtssgg1","ssvbf","highMtssvbf", "ssantiisomuongg0", "ssantiisomuongg1", "ssantiisomuonvbf", "ssantiisomuonlowmMtvbf","sslowmMtvbf","lowmMtvbf","ttbarcontrolvbf","ztautaucontrolvbf","ztautaucontrolgg0","ztautaucontrolgg1","highMtztautaucontrolvbf","antiisomuonztautaucontrolvbf","highMtssztautaucontrolvbf","ssantiisomuonztautaucontrolvbf","ssztautaucontrolvbf","antiisotauztautaucontrolvbf","antiisotauztautaucontrolgg0","antiisotauztautaucontrolgg1"]
        namesize = len(names)
	for x in range(0,namesize):

            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "weight_nopu", "Event weight without PU", 100, 0, 5)
	    self.book(names[x], "EmbPtWeight", "embedded weight", 100, 0, 1)
	    self.book(names[x], "NUP", "NUP",11,0,11)
    
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)
    
            self.book(names[x], "mPt", "Muon  Pt", 100, 0, 100)
            self.book(names[x], "mEta", "Muon  eta", 100, -2.5, 2.5)
            self.book(names[x], "mMtToMVAMET", "Muon MT (MVA)", 200, 0, 200)
            self.book(names[x], "mMtToPFMET", "Muon MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "mCharge", "Muon Charge", 5, -2, 2)
            self.book(names[x], "tPt", "Tau  Pt", 200, 0, 200)
	    self.book(names[x], "tMass", "Tau Mass", 200,0,200)
            self.book(names[x], "tEta", "Tau  eta", 100, -2.5, 2.5)
            self.book(names[x], "tMtToMVAMET", "Tau MT (MVA)", 200, 0, 200)
            self.book(names[x], "tMtToPFMET", "Tau MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "tCharge", "Tau  Charge", 5, -2, 2)
	    self.book(names[x], "tJetPt", "Tau Jet Pt" , 500, 0 ,500)	    

            self.book(names[x], 'mPixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
            self.book(names[x], 'mJetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)
    	    
	    self.book(names[x],"fullMT_mva","fullMT_mva",500,0,500);
            self.book(names[x],"fullMT","fullMT",500,0,500);
	    self.book(names[x],"collMass","collMass",500,0,500);
	    self.book(names[x],"higgsMass","higgsMass",500,0,500);
            self.book(names[x],"fullPT_mva","fullPT_mva",500,0,500);
            self.book(names[x],"fullPT","fullPT",500,0,500);	    
    	    self.book(names[x], "LT", "ht", 400, 0, 400)
    
            self.book(names[x], "pfMetEt", "Type1 MET", 200, 0, 200)
            self.book(names[x], "mva_metEt", "MVA MET", 200, 0, 200)
    
    
            self.book(names[x], "m_t_Mass", "Muon + Tau Mass", 200, 0, 200)
            self.book(names[x], "m_t_Pt", "Muon + Tau Pt", 200, 0, 200)
            self.book(names[x], "m_t_DR", "Muon + Tau DR", 100, 0, 10)
            self.book(names[x], "m_t_DPhi", "Muon + Tau DPhi", 100, 0, 4)
            self.book(names[x], "m_t_SS", "Muon + Tau SS", 5, -2, 2)
            self.book(names[x], "m_t_ToMETDPhi_Ty1", "Muon Tau DPhi to MET", 100, 0, 4)
    
            # Vetoes
            self.book(names[x], 'bjetVeto', 'Number of b-jets', 5, -0.5, 4.5)
            self.book(names[x], 'bjetCSVVeto', 'Number of b-jets', 5, -0.5, 4.5)
            self.book(names[x], 'muVetoPt5IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
	    self.book(names[x], 'muVetoPt15IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
            #self.book(names[x], 'tauVetoPt20', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(names[x], 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
	    if isZtt:
		self.book(names[x], 'isZtautau', 'is Ztautau event', 2,-0.5,1.5)
    
            self.book(names[x], 'jetVeto20', 'Number of extra jets', 5, -0.5, 4.5)
            self.book(names[x], 'jetVeto30', 'Number of extra jets', 5, -0.5, 4.5)	
	    #Isolation
	    self.book(names[x], 'mRelPFIsoDB' ,'Muon Isolation', 100, 0.0,1.0)
   
 
            self.book(names[x], "mPhiMtPhi", "", 100, 0,4)
            self.book(names[x], "mPhiMETPhiMVA", "", 100, 0,4)
            self.book(names[x], "tPhiMETPhiMVA", "", 100, 0,4)
            self.book(names[x], "mPhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "tPhiMETPhiType1", "", 100, 0,4)
	    self.book(names[x], "tDecayMode", "" , 11,  -0.5 , 10.5)

### vbf ###
            self.book(names[x], "vbfJetVeto30", "central jet veto for vbf", 5, -0.5, 4.5)
	    self.book(names[x], "vbfJetVeto20", "", 5, -0.5, 4.5)
	    self.book(names[x], "vbfMVA", "", 100, 0,0.5)
	    self.book(names[x], "vbfMass", "", 200,0,1000.0)
	    self.book(names[x], "vbfDeta", "", 100, -0.5,10.0)
            self.book(names[x], "vbfj1eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfj2eta","",100,-2.5,2.5)
            self.book(names[x], "vbfj1pt","",100,0,400.0)
            self.book(names[x], "vbfj2pt","",100,0,400.0)
	    self.book(names[x], "vbfVispt","",100,0,200)
	    self.book(names[x], "vbfHrap","",100,0,5.0)
	    self.book(names[x], "vbfDijetrap","",100,0,5.0)
	    self.book(names[x], "vbfDphihj","",100,0,4)
            self.book(names[x], "vbfDphihjnomet","",100,0,4)


    def correction(self, row):
        return mc_corrector(row)
    
    def fakeRateMethod(self,row,fakeName='Central'):
	return getFakeRateFactor(row,fakeName)

    def fill_histos(self, row,name='gg', fakeRate = False, fakeName = 'Central'):
        histos = self.histograms
	if isEmbed == True:
		weight = row.EmbPtWeight*self.correction(row)
	else:
		weight = self.correction(row)
	if fakeRate == True:
		weight = weight * self.fakeRateMethod(row,fakeName)
        histos[name+'/weight'].Fill(weight)
        histos[name+'/weight_nopu'].Fill(self.correction(row))
	histos[name+'/EmbPtWeight'].Fill(row.EmbPtWeight)
	histos[name+'/NUP'].Fill(row.NUP, weight)
        histos[name+'/rho'].Fill(row.rho, weight)
        histos[name+'/nvtx'].Fill(row.nvtx, weight)
        histos[name+'/prescale'].Fill(row.doubleMuPrescale, weight)

        histos[name+'/mPt'].Fill(row.mPt, weight)
        histos[name+'/mEta'].Fill(row.mEta, weight)
	histos[name+'/mMtToMVAMET'].Fill(row.mMtToMVAMET,weight)
        histos[name+'/mMtToPFMET'].Fill(row.mMtToPFMET,weight)
        histos[name+'/mCharge'].Fill(row.mCharge, weight)
        histos[name+'/tPt'].Fill(row.tPt, weight)
        histos[name+'/tMass'].Fill(row.tMass, weight)
        histos[name+'/tEta'].Fill(row.tEta, weight)
        histos[name+'/tMtToMVAMET'].Fill(row.tMtToMVAMET,weight)
        histos[name+'/tMtToPFMET'].Fill(row.tMtToPFMET,weight)
        histos[name+'/tCharge'].Fill(row.tCharge, weight)
	histos[name+'/tJetPt'].Fill(row.tJetPt, weight)

	histos[name+'/LT'].Fill(row.LT,weight)

        histos[name+'/fullMT_mva'].Fill(fullMT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
        histos[name+'/fullMT'].Fill(fullMT(row.pfMetEt,row.mPt,row.tPt,row.pfMetPhi, row.mPhi, row.tPhi),weight)
	histos[name+'/collMass'].Fill(collMass(row), weight)
	histos[name+'/higgsMass'].Fill(MH(row),weight)
        histos[name+'/fullPT_mva'].Fill(fullPT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
        histos[name+'/fullPT'].Fill(fullPT(row.pfMetEt,row.mPt,row.tPt,row.pfMetPhi, row.mPhi, row.tPhi),weight) 


	histos[name+'/pfMetEt'].Fill(row.pfMetEt,weight)
	histos[name+'/mva_metEt'].Fill(row.mva_metEt,weight)

        histos[name+'/m_t_Mass'].Fill(row.m_t_Mass,weight)
        histos[name+'/m_t_Pt'].Fill(row.m_t_Pt,weight)
        histos[name+'/m_t_DR'].Fill(row.m_t_DR,weight)
        histos[name+'/m_t_DPhi'].Fill(row.m_t_DPhi,weight)
        histos[name+'/m_t_SS'].Fill(row.m_t_SS,weight)
	histos[name+'/m_t_ToMETDPhi_Ty1'].Fill(row.m_t_ToMETDPhi_Ty1,weight)

        histos[name+'/mPixHits'].Fill(row.mPixHits, weight)
        histos[name+'/mJetBtag'].Fill(row.mJetBtag, weight)

        histos[name+'/bjetVeto'].Fill(row.bjetVeto, weight)
        histos[name+'/bjetCSVVeto'].Fill(row.bjetCSVVeto, weight)
        histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
        histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
        #histos[name+'/tauVetoPt20'].Fill(row.tauVetoPt20, weight)
        histos[name+'/eVetoCicTightIso'].Fill(row.eVetoCicTightIso, weight)
	if isZtt:
		histos[name+'/isZtautau'].Fill(row.isZtautau, weight)
        histos[name+'/jetVeto20'].Fill(row.jetVeto20, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)

	histos[name+'/mRelPFIsoDB'].Fill(row.mRelPFIsoDB, weight)
        
	histos[name+'/mPhiMtPhi'].Fill(deltaPhi(row.mPhi,row.tPhi),weight)
        histos[name+'/mPhiMETPhiMVA'].Fill(deltaPhi(row.mPhi,row.mva_metPhi),weight)
        histos[name+'/tPhiMETPhiMVA'].Fill(deltaPhi(row.tPhi,row.mva_metPhi),weight)
        histos[name+'/mPhiMETPhiType1'].Fill(deltaPhi(row.mPhi,row.pfMetPhi),weight)
        histos[name+'/tPhiMETPhiType1'].Fill(deltaPhi(row.tPhi,row.pfMetPhi),weight)
	histos[name+'/tDecayMode'].Fill(row.tDecayMode, weight)
	histos[name+'/vbfJetVeto30'].Fill(row.vbfJetVeto30, weight)
     	histos[name+'/vbfJetVeto20'].Fill(row.vbfJetVeto20, weight)
        histos[name+'/vbfMVA'].Fill(row.vbfMVA, weight)
        histos[name+'/vbfMass'].Fill(row.vbfMass, weight)
        histos[name+'/vbfDeta'].Fill(row.vbfDeta, weight)
        histos[name+'/vbfj1eta'].Fill(row.vbfj1eta, weight)
        histos[name+'/vbfj2eta'].Fill(row.vbfj2eta, weight)
        histos[name+'/vbfj1pt'].Fill(row.vbfj1pt, weight)
        histos[name+'/vbfj2pt'].Fill(row.vbfj2pt, weight)
        histos[name+'/vbfVispt'].Fill(row.vbfVispt, weight)
        histos[name+'/vbfHrap'].Fill(row.vbfHrap, weight)
        histos[name+'/vbfDijetrap'].Fill(row.vbfDijetrap, weight)
        histos[name+'/vbfDphihj'].Fill(row.vbfDphihj, weight)
        histos[name+'/vbfDphihjnomet'].Fill(row.vbfDphihjnomet, weight)






    def presel(self, row):
	if not row.isoMu24eta2p1Pass:
            return False
        return True
    
    def twojets(self,row):
	if Twojets == True and row.vbfNJets<2:
		return False
	return True    	
    def kinematics(self, row):
        if row.mPt < 30:
            return False
        if abs(row.mEta) >= 2.1:
            return False
        if row.tPt < 30:
            return False
        if abs(row.tEta)>=2.3 :
            return False
        return True

    def lowMt(self, row):
	if preselection == True:
		return True
	elif optimized == False:
		#return True
                if row.tMtToPFMET > 35:
			return False
	else:
	###ORIGINAL FOR PREAPPROVAL###
		if row.tMtToPFMET > 35:
	###REOPTIMIZED FOR S/(SQRT(S+B))###
	#	if row.tMtToPFMET > 55:
	    		return False
	return True

    def lowmMt(self,row):
	if row.mMtToPFMET>20:
		return False
	return True

    def highMt(self,row):
	if row.mMtToPFMET<70:
	    return False
	return True
    
    def hightMt(self,row):
	if row.tMtToPFMET<70:
		return False
	return True

    def lowtMt(self,row):
        if row.tMtToPFMET>20:
                return False
        return True

    def mMtgg(self,row):
	if optimized == False and preselection == False:
		if row.mMtToPFMET < 30:
			return False	
	return True	
    def lowm_t_Mass(self,row):
	if row.m_t_Mass > 60:
		return False
	return True

    def lowm_t_DR(self,row):
	if row.m_t_DR > 2:
		return False
	return True
    def ggtight(self,row,jets):
       #if row.LT<75:
        #  return False
       if optimized == True and jets == 0:
          #if (deltaPhi(row.mPhi,row.pfMetPhi)) < 2.5:
	#	print "failing 1"
	 #  	return False
	  #if deltaPhi(row.tPhi,row.pfMetPhi) > 0.3:
	#	print "failing 2"
	#	return False
	  if deltaPhi(row.mPhi, row.tPhi) <2.7:
		return False
	  if row.mPt < 45:
		#print "failing 3"
		return False
	  if row.tPt < 35:
		#print "failing 4"
		return False
	  if row.tMtToPFMET > 50:
		#print "failing 5"
		return False
          if row.jetVeto30!=0:
		#print "failing 6"
		return False
       elif optimized == True and jets ==1:

          if row.jetVeto30!=1:
          	return False
          #if deltaPhi(row.tPhi,row.pfMetPhi) > 0.3:
         #       return False
          if row.mPt < 35:
                return False
          if row.tPt < 40:
                return False
          if row.tMtToPFMET > 35:
                return False


       else:
       		if row.mPt < 50:
          	 	return False
       		if (deltaPhi(row.tPhi,row.mva_metPhi)>0.5):
          		return False
		if jets==1:
			if row.jetVeto30!=1:
				return False
		if jets==0:
			if row.jetVeto30!=0:
				return False
       
#       if (deltaPhi(row.mPhi,row.tPhi)<2.):    
#          return False
       return True	


    def vbf(self,row):
	if(abs(row.vbfDeta)<3.5):
	    return False
	if vbfMassCut500 == True:
        	if row.vbfMass < 500:
	    		return False
	else:
		if row.vbfMass < 400:
			return False
	if row.jetVeto30 < 2:
	    return False
	if row.vbfJetVeto30 > 0:
	    return False
        return True

    def loosevbf(self,row):
	if(abs(row.vbfDeta)<2.0):
            return False
	if vbfMassCut500 == True:
        	if row.vbfMass < 200:
            		return False
        #if row.jetVeto20 < 2:
         #   return False
        if row.vbfJetVeto30 > 0:
            return False
        return True
    def optimizedvbf(self,row):
###ORIGINAL USED IN PREAPPROVAL####
	if row.mPt < 30:
		return False
	if row.tPt < 40:
		return False
	if (abs(row.vbfDeta)) < 3.5:
		return False
	if row.vbfMass < 550:
		return False
        if row.jetVeto30 < 2:
                return False
        if row.vbfJetVeto30 > 0:
                return False
        return True
    def intermediatevbf(self,row):
	if row.mPt < 30:
		return False
	if row.tPt < 30:
		return False
	if abs(row.vbfDeta) < 2.0:
		return False
	if row.vbfMass < 200:
		return False
        if row.jetVeto30 < 2:
                return False
        if row.vbfJetVeto30 > 0:
                return False
        return True
##OPTIMIZED AT BR=1%###
#        if row.mPt < 30:
#               return False
#        if row.tPt < 35:
#               return False
#        if (abs(row.vbfDeta)) < 3.2:
#               return False
#        if row.vbfMass < 700:
#               return False
#	if row.jetVeto30 < 2:
#		return False
#	if row.vbfJetVeto30 > 0:
#		return False
#OPTIMIZED AT S/SQRT(S+B)#####
#        if row.mPt < 30:
#               return False
#        if row.tPt < 30:
#               return False
#        if (abs(row.vbfDeta)) < 2.4:
#               return False
#        if row.vbfMass < 500:
#               return False
    def oppositesign(self,row):
	if row.mCharge*row.tCharge!=-1:
            return False
	return True

    def obj1_id(self, row):
        return bool(row.mPFIDTight)  and bool(abs(row.mDZ) < 0.2) 

    def obj2_id(self, row):
        if newTauID:
		return  row.tAntiElectronLoose and row.tAntiMuon2Tight and row.tDecayFindingNewDMs
	else:
                return  row.tAntiElectronLoose and row.tAntiMuonTight2 and row.tDecayFinding

    def vetos(self,row):
	if Twomu == False:# and newTauID:
		return  (bool (row.muVetoPt5IsoIdVtx<1) and bool (row.eVetoCicTightIso<1))# and bool (row.tauVetoPt20Loose3HitsVtx<1) )
	else:
		return (bool (row.muVetoPt15IsoIdVtx>0) and bool (row.eVetoCicTightIso<1))
    def Ztautauveto(self,row):
	if otherJets == False:
        	return (bool (row.isZtautau == 1 ))
	else:
		return (bool (row.isZtautau != 1))
    def NUPveto(self,row):
	return (bool (row.NUP == 5))
    def obj1_iso(self, row):
        return bool(row.mRelPFIsoDB <0.12)
    def obj2_iso(self, row, isgg=False):
	#if wjets==True and isgg==True:
	if newTauID:
		return row.tVTightIsoMVA3NewDMLT
	else:
        	#return  row.tLooseIso3Hits
		return row.tTightIso3Hits
	#else:
	#	return row.tTightIso3Hits

    def obj2_mediso(self, row):
	return row.tMediumIso3Hits

    def obj1_antiiso(self,row):
        return bool(row.mRelPFIsoDB >0.2) 

    def obj2_antiiso(self, row):
        #return not row.tLooseIso

	#anti-iso test
        if newTauID:
		return ((not row.tVTightIsoMVA3NewDMLT) and row.tMediumIsoMVA3NewDMLT)
	else:
		return ((not row.tTightIso3Hits) and row.tLooseIso3Hits)
 	
	#return not row.tTightIso3Hits
    def ttbarcontrol(self,row):
	if row.tMtToPFMET < 130:
		return False
	return True

    def zerojet(self,row):
	if row.jetVeto30 == 0:
		return True
    def onejet(self,row):
	if row.jetVeto30 ==1:
		return True
    def process(self):
        for row in self.tree:
	    if Isiso == False:
		obj1iso = True
		obj2iso = True
		obj2isogg = True
	    else:
		obj1iso = self.obj1_iso(row)
		obj2iso = self.obj2_iso(row)
		obj2isogg = self.obj2_iso(row,True)

	    if not self.presel(row):
		continue
            if not self.kinematics(row):
                continue
            if not self.obj1_id(row): 
                continue
            if not self.obj2_id(row):
                continue

	    if not self.vetos(row):
		continue 
	    if isZtt:	
            	if not self.Ztautauveto(row):
                	continue
	    if isZjetsM50:
		if not self.NUPveto(row):	    
			continue

	    if preselection == True:
	    	tightcutgg0=self.zerojet(row)
		tightcutgg1=self.onejet(row)
	    else: 
		tightcutgg0 = self.ggtight(row,0)
		tightcutgg1 = self.ggtight(row,1)		
	    if tightcutgg0:
		if self.mMtgg(row):
			if obj1iso and obj2isogg and self.oppositesign(row): 
	        	        self.fill_histos(row,'gg0')
                                if self.lowm_t_DR(row):
                                        self.fill_histos(row,'ztautaucontrolgg0')
                	if obj1iso and obj2isogg and not self.oppositesign(row):
                	        self.fill_histos(row,'ssgg0')
                	if self.obj1_antiiso(row) and obj2isogg and  self.oppositesign(row):
                        	self.fill_histos(row,'antiisomuongg0')
			if self.obj1_antiiso(row) and obj2isogg and not self.oppositesign(row):
                        	self.fill_histos(row,'ssantiisomuongg0')
			if obj1iso and self.obj2_antiiso(row) and self.oppositesign(row):
				self.fill_histos(row, 'antiisotaugg0',True)
                                if self.lowm_t_DR(row):
                                        self.fill_histos(row,'antiisotauztautaucontrolgg0',True)
		if self.highMt(row):
                        if obj1iso and obj2isogg and self.oppositesign(row):
                                self.fill_histos(row,'highMtgg0')
                        if obj1iso and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'highMtssgg0')
            if tightcutgg1:

                if self.mMtgg(row):
                        if obj1iso and obj2isogg and self.oppositesign(row):
                                self.fill_histos(row,'gg1')
                                if self.lowm_t_DR(row):
                                        self.fill_histos(row,'ztautaucontrolgg1')
                        if obj1iso and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'ssgg1')
                        if self.obj1_antiiso(row) and obj2isogg and  self.oppositesign(row):
                                self.fill_histos(row,'antiisomuongg1')
                        if self.obj1_antiiso(row) and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'ssantiisomuongg1')
                        if obj1iso and self.obj2_antiiso(row) and self.oppositesign(row):
                                self.fill_histos(row, 'antiisotaugg1',True)
                                if self.lowm_t_DR(row):
                                        self.fill_histos(row,'antiisotauztautaucontrolgg1',True)
                if self.highMt(row):
                        if obj1iso and obj2isogg and self.oppositesign(row):
                                self.fill_histos(row,'highMtgg1')
                        if obj1iso and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'highMtssgg1')	

	    if preselection == True:
	    	loosecutvbf = True
		tightcutvbf = True
	    else:
		loosecutvbf = self.loosevbf(row)
		if optimized == True:
			tightcutvbf = self.optimizedvbf(row)
		elif intermediate == True:
			tightcutvbf = self.intermediatevbf(row)
		else:
			tightcutvbf = self.vbf(row)

	    if loosecutvbf:
		if not self.twojets(row):
			continue
		if self.lowMt(row):
                        if self.obj1_antiiso(row) and self.obj2_mediso(row) and self.oppositesign(row):
                                self.fill_histos(row,'antiisomuonvbf')
				if self.lowm_t_DR(row):
					self.fill_histos(row,'antiisomuonztautaucontrolvbf') 
	    #print "about to apply tightcutvbf \n"
	    if tightcutvbf:	
		#print "passsed tightcutvbf"
		if not self.twojets(row):
                	continue
		if self.lowMt(row):
                	if obj1iso and obj2iso and self.oppositesign(row):
				#print "all cuts passed ready to fill \n"
                        	self.fill_histos(row,'vbf')
				#print "run: " + str(row.run)
				#print "evt: " + str(row.evt)
				if self.ttbarcontrol(row):
					self.fill_histos(row,'ttbarcontrolvbf')
				if self.lowm_t_DR(row):
					self.fill_histos(row,'ztautaucontrolvbf')
                	if obj1iso and obj2iso and not self.oppositesign(row):
                        	self.fill_histos(row,'ssvbf')
				if self.lowm_t_DR(row):
					self.fill_histos(row,'ssztautaucontrolvbf')
                        if self.obj1_antiiso(row) and obj2iso and not self.oppositesign(row):
                                self.fill_histos(row,'ssantiisomuonvbf')
				if self.lowm_t_DR(row):
					self.fill_histos(row,'ssantiisomuonztautaucontrolvbf')
	                if obj1iso and self.obj2_antiiso(row) and self.oppositesign(row):
                                self.fill_histos(row, 'antiisotauvbf',True,'Central')
				self.fill_histos(row, 'antiisotauvbfdown',True,'Down')
				self.fill_histos(row, 'antiisotauvbfup',True,'Up')
                                if self.lowm_t_DR(row):
                                        self.fill_histos(row,'antiisotauztautaucontrolvbf',True)

		if self.highMt(row):
                        if obj1iso and obj2iso and self.oppositesign(row):
                                self.fill_histos(row,'highMtvbf')
				if self.lowm_t_DR(row):
					self.fill_histos(row,'highMtztautaucontrolvbf')
			if self.obj1_antiiso(row) and obj2iso and self.oppositesign(row):
				self.fill_histos(row,'highMtantiisomuonvbf')
                        if obj1iso and obj2iso and not self.oppositesign(row):
                                self.fill_histos(row,'highMtssvbf')
				if self.lowm_t_DR(row):
					self.fill_histos(row,'highMtssztautaucontrolvbf')
			if self.obj1_antiiso(row) and obj2iso and not self.oppositesign(row):
				self.fill_histos(row,'highMtssantiisomuonvbf')
			if obj1iso and self.obj2_antiiso(row) and self.oppositesign(row):
				self.fill_histos(row, 'antiisotauhighMtvbf',True)
			if obj1iso and self.obj2_antiiso(row) and not self.oppositesign(row):
				self.fill_histos(row, 'highMtssantiisotauvbf',True)
				self.fill_histos(row, 'highMtssantiisotauvbfdown',True,'Down')
				self.fill_histos(row, 'highMtssantiisotauvbfup',True,'Up')

		if self.lowmMt(row):
			if obj1iso and obj2iso and self.oppositesign(row):
				self.fill_histos(row,'lowmMtvbf')
                        if obj1iso and obj2iso and not self.oppositesign(row):
                                self.fill_histos(row,'sslowmMtvbf')
                        if self.obj1_antiiso(row) and obj2iso and not self.oppositesign(row):
                                self.fill_histos(row,'ssantiisomuonlowmMtvbf')

		if self.hightMt(row):
			if obj1iso and obj2iso and self.oppositesign(row):
				self.fill_histos(row,'hightMtvbf')			
                if self.lowtMt(row):
                        if obj1iso and obj2iso and self.oppositesign(row):
                                self.fill_histos(row,'lowtMtvbf')



    def finish(self):
        self.write_histos()
