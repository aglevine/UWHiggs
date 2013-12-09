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
if preselection == True:
	optimized=False
Isiso = bool('true' in os.environ['iso'])
wjets = bool('true' in os.environ['wjets'])
Twomu = False
Twojets = True
#fakeRateDown=False
#fakeRateUp=False
#print "fakeRateDown" + str(fakeRateDown)
#print "fakeRateUp" + str(fakeRateUp)
if preselection == False:
	Twojets = True
vbfMassCut500 = True
#if vbfMassCut500 == True:
#	vbfMassCutstr = "500 GeV"
#else:
#	vbfMassCutstr = "400 GeV"
is7TeV = bool('7TeV' in os.environ['jobid'])
isPU1signal = bool ('false' in os.environ['PU'])
isEmbed = bool ('true' in os.environ['embed'])
isZtt = bool ('true' in os.environ['Ztt'])
print "Is wjets", wjets
print "Is 7TeV:", is7TeV
print "Is PU1 signal:", isPU1signal
print "Is Embed: ", isEmbed
print "Is Ztt: ", isZtt
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

def collMass_type1(row):
        taupx = row.tPt*math.cos(row.tPhi)
        taupy = row.tPt*math.sin(row.tPhi)
        taupt=  row.tPt

        metpx = row.type1_pfMetEt*math.cos(row.type1_pfMetPhi)
        metpy = row.type1_pfMetEt*math.sin(row.type1_pfMetPhi)
        met = row.type1_pfMetEt

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
    #if preselection == True:
    mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/preselection/isotrue/AnalyzeMuMuTauTight/"
    mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/preselection/isofalse/AnalyzeMuMuTauTight/"
    #if preselection == False:
     #   mumuidiso_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/signal/isotrue/AnalyzeMuMuTauTight/"
     #   mumuid_dir = "/afs/hep.wisc.edu/cms/aglevine/hlfv_5_3_9/src/UWHiggs/lfvmutau/MuMuTau_Sept30_loosevbf/mmt/signal/isofalse/AnalyzeMuMuTauTight/"
    #mumuidiso_dir = "MuMuIdIso/"
    #mumuid_dir = "MuMuId/"
    data_file_str = "datammt_2012.root"
    zjets_file_str = "Zjetsmmtvbf.root"
    zjets_file_idiso = ROOT.TFile(mumuidiso_dir+zjets_file_str)
    zjets_file_id = ROOT.TFile(mumuid_dir+zjets_file_str)
    data_file_idiso = ROOT.TFile(mumuidiso_dir+data_file_str)
    data_file_id = ROOT.TFile(mumuid_dir+data_file_str)
    histoIdIso = data_file_idiso.Get("vbf/tJetPt").Clone()
    histoId = data_file_id.Get("vbf/tJetPt").Clone()
    i = 0
    bincount = -1
    binning = array.array('d',[])
    while (i <= 480):
    	if i < 70:
        	binning.append(i)
                i = i+2
        elif i< 100:
                binning.append(i)
                i = i+5
        elif i<120:
        	binning.append(i)
                i = i+10
        else:
        	binning.append(i)
                i = i+20
        bincount = bincount+1
    #print bincount
    #print len(binning)
    #print binning
    histoIdIsoRebin = histoIdIso.Rebin(bincount,"histoIdIsoRebin",binning)
    histoIdRebin = histoId.Rebin(bincount,"histoIdRebin",binning)
    histo_ftau = histoIdIsoRebin.Clone()
    histo_ftau.Divide(histoIdRebin)
    
    #histoIdIso.Rebin(binwidth)
    #histoId.Rebin(binwidth)
    #histoIdIso = zjets_file_idiso.Get("vbf/tJetPt").Clone()
    #histoId = zjets_file_id.Get("vbf/tJetPt").Clone()
    if histoIdRebin.GetBinContent(histoId.FindBin(row.tJetPt)) == 0:
	fTauIso = 0
    else:
	fTauIso = histo_ftau.GetBinContent(histo_ftau.FindBin(row.tJetPt))
    	if fakeName=='Down':
		fTauIso = fTauIso - histo_ftau.GetBinError(histo_ftau.FindBin(row.tJetPt))
		if fTauIso < 0:
			fTauIso = 0
	elif fakeName=='Up':
                fTauIso = fTauIso + histo_ftau.GetBinError(histo_ftau.FindBin(row.tJetPt))
    #print "fTauIso: " + str(fTauIso)
    #print "tauJetPt: " + str(row.tJetPt)
    fakeRateFactor = fTauIso/(1-fTauIso)
    return fakeRateFactor

class AnalyzeMuTauTightvbfNew(MegaBase):
    tree = 'New_Tree'	
    #tree = 'mt/final/Ntuple'
    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeMuTauTightvbfNew, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuTauTree.MuTauTree(tree)
        self.out = outfile
        # Histograms for each category
        self.histograms = {}
        self.is7TeV = '7TeV' in os.environ['jobid']

    def begin(self):

        names=["gg0","gg1","vbf","highMtgg0","highMtgg1","highMtvbf", "hightMtvbf","lowtMtvbf","antiisomuongg0","antiisomuongg1","antiisomuonvbf","highMtantiisomuonvbf","antiisotaugg0","antiisotaugg1","antiisotauvbf","antiisotauvbfdown","antiisotauvbfup","highMtssantiisotauvbf","antiisotauhighMtvbf","highMtssantiisomuonvbf","ssgg0","ssgg1","highMtssgg0", "highMtssgg1","ssvbf","highMtssvbf", "ssantiisomuongg0", "ssantiisomuongg1", "ssantiisomuonvbf", "ssantiisomuonlowmMtvbf","sslowmMtvbf","lowmMtvbf","ttbarcontrolvbf","ztautaucontrolvbf","highMtztautaucontrolvbf","antiisomuonztautaucontrolvbf","highMtssztautaucontrolvbf","ssantiisomuonztautaucontrolvbf","ssztautaucontrolvbf"]
        namesize = len(names)
	for x in range(0,namesize):

            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "weight_nopu", "Event weight without PU", 100, 0, 5)
	    self.book(names[x], "EmbPtWeight", "embedded weight", 100, 0, 1)
    
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 31, -0.5, 30.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)
    
            self.book(names[x], "mPt", "Muon  Pt", 100, 0, 100)
            self.book(names[x], "mEta", "Muon  eta", 100, -2.5, 2.5)
            self.book(names[x], "mMtToMVAMET", "Muon MT (MVA)", 200, 0, 200)
            self.book(names[x], "mMtToPfMet_Ty1", "Muon MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "mCharge", "Muon Charge", 5, -2, 2)
            self.book(names[x], "tPt", "Tau  Pt", 200, 0, 200)
            self.book(names[x], "tEta", "Tau  eta", 100, -2.5, 2.5)
            self.book(names[x], "tMtToMVAMET", "Tau MT (MVA)", 200, 0, 200)
            self.book(names[x], "tMtToPfMet_Ty1", "Tau MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "tCharge", "Tau  Charge", 5, -2, 2)
	    self.book(names[x], "tJetPt", "Tau Jet Pt" , 500, 0 ,500)	    
		        
            self.book(names[x], 'mPixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
            self.book(names[x], 'mJetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)
    	    
	    self.book(names[x],"fullMT_mva","fullMT_mva",500,0,500);
            self.book(names[x],"fullMT_type1","fullMT_type1",500,0,500);
	    self.book(names[x],"collMass_type1","collMass_type1",500,0,500);
            self.book(names[x],"fullPT_mva","fullPT_mva",500,0,500);
            self.book(names[x],"fullPT_type1","fullPT_type1",500,0,500);	    
    	    self.book(names[x], "LT", "ht", 400, 0, 400)
    
            self.book(names[x], "type1_pfMetEt", "Type1 MET", 200, 0, 200)
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
            self.book(names[x], 'tauVetoPt20', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(names[x], 'eVetoCicTightIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
	    if isZtt:
		self.book(names[x], 'isZtautau', 'is Ztautua event', 2,-0.5,1.5)
    
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
		weight = row.EmbPtWeight
	else:
		weight = self.correction(row)
	if fakeRate == True:
		weight = weight * self.fakeRateMethod(row,fakeName)
        histos[name+'/weight'].Fill(weight)
        histos[name+'/weight_nopu'].Fill(self.correction(row))
	histos[name+'/EmbPtWeight'].Fill(row.EmbPtWeight)
        histos[name+'/rho'].Fill(row.rho, weight)
        histos[name+'/nvtx'].Fill(row.nvtx, weight)
        histos[name+'/prescale'].Fill(row.doubleMuPrescale, weight)

        histos[name+'/mPt'].Fill(row.mPt, weight)
        histos[name+'/mEta'].Fill(row.mEta, weight)
	histos[name+'/mMtToMVAMET'].Fill(row.mMtToMVAMET,weight)
        histos[name+'/mMtToPfMet_Ty1'].Fill(row.mMtToPfMet_Ty1,weight)
        histos[name+'/mCharge'].Fill(row.mCharge, weight)
        histos[name+'/tPt'].Fill(row.tPt, weight)
        histos[name+'/tEta'].Fill(row.tEta, weight)
        histos[name+'/tMtToMVAMET'].Fill(row.tMtToMVAMET,weight)
        histos[name+'/tMtToPfMet_Ty1'].Fill(row.tMtToPfMet_Ty1,weight)
        histos[name+'/tCharge'].Fill(row.tCharge, weight)
	histos[name+'/tJetPt'].Fill(row.tJetPt, weight)

	histos[name+'/LT'].Fill(row.LT,weight)

        histos[name+'/fullMT_mva'].Fill(fullMT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
        histos[name+'/fullMT_type1'].Fill(fullMT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight)
	histos[name+'/collMass_type1'].Fill(collMass_type1(row), weight)
        histos[name+'/fullPT_mva'].Fill(fullPT(row.mva_metEt,row.mPt,row.tPt,row.mva_metPhi, row.mPhi, row.tPhi),weight)
        histos[name+'/fullPT_type1'].Fill(fullPT(row.type1_pfMetEt,row.mPt,row.tPt,row.type1_pfMetPhi, row.mPhi, row.tPhi),weight) 


	histos[name+'/type1_pfMetEt'].Fill(row.type1_pfMetEt,weight)
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
        histos[name+'/tauVetoPt20'].Fill(row.tauVetoPt20, weight)
        histos[name+'/eVetoCicTightIso'].Fill(row.eVetoCicTightIso, weight)
	if isZtt:
		histos[name+'/isZtautau'].Fill(row.isZtautau, weight)
        histos[name+'/jetVeto20'].Fill(row.jetVeto20, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)

	histos[name+'/mRelPFIsoDB'].Fill(row.mRelPFIsoDB, weight)
        
	histos[name+'/mPhiMtPhi'].Fill(deltaPhi(row.mPhi,row.tPhi),weight)
        histos[name+'/mPhiMETPhiMVA'].Fill(deltaPhi(row.mPhi,row.mva_metPhi),weight)
        histos[name+'/tPhiMETPhiMVA'].Fill(deltaPhi(row.tPhi,row.mva_metPhi),weight)
        histos[name+'/mPhiMETPhiType1'].Fill(deltaPhi(row.mPhi,row.type1_pfMetPhi),weight)
        histos[name+'/tPhiMETPhiType1'].Fill(deltaPhi(row.tPhi,row.type1_pfMetPhi),weight)
	histos[name+'/tDecayMode'].Fill(row.tDecayMode, weight)
	histos[name+'/vbfJetVeto30'].Fill(row.vbfJetVeto30, weight)
     	histos[name+'/vbfJetVeto20'].Fill(row.vbfJetVeto20, weight)
        histos[name+'/vbfMVA'].Fill(row.vbfMVA, weight)
        histos[name+'/vbfMass'].Fill(row.vbfMass, weight)
        histos[name+'/vbfDeta'].Fill(row.vbfDeta, weight)
        histos[name+'/vbfj1eta'].Fill(row.vbfj1eta, weight)
        histos[name+'/vbfj2eta'].Fill(row.vbfj2eta, weight)
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
        if row.tPt<20 :
            return False
        if abs(row.tEta)>=2.3 :
            return False
        return True

    def lowMt(self, row):
	if optimized == False:
		if row.tMtToPfMet_Ty1>20:
			return False
	else:
		if row.tMtToPfMet_Ty1 > 30:
	    		return False
	return True

    def lowmMt(self,row):
	if row.mMtToPfMet_Ty1>20:
		return False
	return True

    def highMt(self,row):
	if row.mMtToPfMet_Ty1<70:
	    return False
	return True
    
    def hightMt(self,row):
	if row.tMtToPfMet_Ty1<70:
		return False
	return True
	
    def lowtMt(self,row):
        if row.tMtToPfMet_Ty1>20:
                return False
        return True

    def mMtgg(self,row):
	if optimized == False:
		if row.mMtToPfMet_Ty1 < 30:
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
    def ggtight(self,row, jets):
       if row.LT<75:
          return False
       if optimized == True and jets == 0:
          if (deltaPhi(row.mPhi,row.type1_pfMetPhi)) < 2.5:
		print "failing 1"
	   	return False
	  if deltaPhi(row.tPhi,row.type1_pfMetPhi) > 0.3:
		print "failing 2"
		return False
	  if row.mPt < 40:
		print "failing 3"
		return False
	  if row.tPt < 25:
		print "failing 4"
		return False
	  if row.tMtToPfMet_Ty1 > 10:
		print "failing 5"
		return False
          if row.jetVeto30!=0:
		print "failing 6"
		return False
       if optimized == True and jets ==1:
	  	
          if row.jetVeto30!=1:
          	return False
          if deltaPhi(row.tPhi,row.type1_pfMetPhi) > 0.3:
                return False
          if row.mPt < 40:
                return False
          if row.tPt < 25:
                return False
          if row.tMtToPfMet_Ty1 > 10:
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
	if row.mPt < 45:
		return False
	if row.tPt < 45:
		return False
	if (abs(row.vbfDeta)) < 3.5:
		return False
	if row.vbfMass < 600:
		return False
	if row.jetVeto30 < 2:
		return False
	if row.vbfJetVeto30 > 0:
		return False
	return True
    def oppositesign(self,row):
	if row.mCharge*row.tCharge!=-1:
            return False
	return True

    def obj1_id(self, row):
        return bool(row.mPFIDTight)  and bool(abs(row.mDZ) < 0.2) 

    def obj2_id(self, row):
	return  row.tAntiElectronLoose and row.tAntiMuonTight2 and row.tDecayFinding

    def vetos(self,row):
	if Twomu == False:
		return  (bool (row.muVetoPt5IsoIdVtx<1) and bool (row.eVetoCicTightIso<1))
	else:
		return (bool (row.muVetoPt15IsoIdVtx>0) and bool (row.eVetoCicTightIso<1))
    def Ztautauveto(self,row):
        if isZtt == True:
                return (bool (row.isZtautau == 1 ))
        else:
                return False


    def obj1_iso(self, row):
        return bool(row.mRelPFIsoDB <0.12)
    def obj2_iso(self, row, isgg=False):
	if wjets==True and isgg==True:
        	#return  row.tLooseIso3Hits
		return row.tTightIso3Hits
	else:
		return row.tTightIso3Hits

    def obj2_mediso(self, row):
	 return row.tMediumIso3Hits

    def obj1_antiiso(self, row):
        return bool(row.mRelPFIsoDB >0.2) 

    def obj2_antiiso(self, row):
        return  not row.tLooseIso
    def ttbarcontrol(self,row):
	if row.mMtToPfMet_Ty1 < 130:
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
	    
	    if preselection == True:
	    	tightcutgg0=self.zerojet(row)
		tightcutgg1=self.onejet(row)
	    else: 
		tightcutgg0 = self.ggtight(row,0)
		tightcutgg1 = self.ggtight(row,1)		
	    if tightcutgg0:
		if self.mMtgg(row):
		        print "passed mMtgg"
			if obj1iso and obj2isogg and self.oppositesign(row): 
				print "about to fill"
	        	        self.fill_histos(row,'gg0')
                	if obj1iso and obj2isogg and not self.oppositesign(row):
                	        self.fill_histos(row,'ssgg0')
                	if self.obj1_antiiso(row) and obj2isogg and  self.oppositesign(row):
                        	self.fill_histos(row,'antiisomuongg0')
			if self.obj1_antiiso(row) and obj2isogg and not self.oppositesign(row):
                        	self.fill_histos(row,'ssantiisomuongg0')
			if obj1iso and self.obj2_antiiso(row) and self.oppositesign(row):
				self.fill_histos(row, 'antiisotaugg0')
		if self.highMt(row):
                        if obj1iso and obj2isogg and self.oppositesign(row):
                                self.fill_histos(row,'highMtgg0')
                        if obj1iso and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'highMtssgg0')
            if tightcutgg1:

                if self.mMtgg(row):
                        if obj1iso and obj2isogg and self.oppositesign(row):
                                self.fill_histos(row,'gg1')
                        if obj1iso and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'ssgg1')
                        if self.obj1_antiiso(row) and obj2isogg and  self.oppositesign(row):
                                self.fill_histos(row,'antiisomuongg1')
                        if self.obj1_antiiso(row) and obj2isogg and not self.oppositesign(row):
                                self.fill_histos(row,'ssantiisomuongg1')
                        if obj1iso and self.obj2_antiiso(row) and self.oppositesign(row):
                                self.fill_histos(row, 'antiisotaugg1')
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
