'''

Run LFV H->MuMuTau analysis in the mu+tau channel.

Authors: Maria Cepeda, Aaron Levine, Evan K. Friis, UW

'''

import MuMuTauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
#import FinalStateAnalysis.TagAndProbe.H2TauCorrections as H2TauCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
import ROOT
import math

from math import sqrt, pi

data=bool ('true' in os.environ['isRealData'])
ZTauTau = bool('true' in os.environ['isZTauTau'])
ZeroJet = bool('true' in os.environ['isInclusive'])
systematic = os.environ['systematic']

def deltaPhi(phi1, phi2):
  PHI = abs(phi1-phi2)
  if PHI<=pi:
      return PHI
  else:
      return 2*pi-PHI

################################################################################
#### MC-DATA and PU corrections ################################################
################################################################################
pu_distributions = glob.glob(os.path.join(
#    'inputs', os.environ['jobid'], 'data_TauPlusX*pu.root'))
        'inputs', os.environ['jobid'], 'data_SingleMu*pu.root'))
pu_corrector = PileupWeight.PileupWeight('MC_Spring16', *pu_distributions)

muon_pog_PFTight_2016 = MuonPOGCorrections.make_muon_pog_PFTight_2016BCD()
muon_pog_TightIso_2016 = MuonPOGCorrections.make_muon_pog_TightIso_2016BCD()
muon_pog_IsoMu22oIsoTkMu22_2016 = MuonPOGCorrections.make_muon_pog_IsoMu22oIsoTkMu22_2016BCD()

def mc_corrector_2016(row):
  pu = pu_corrector(row.nTruePU)

  m1id = muon_pog_PFTight_2016(row.m1Pt,abs(row.m1Eta))
  m1iso = muon_pog_TightIso_2016('Tight',row.m1Pt,abs(row.m1Eta))
  m1_trg = muon_pog_IsoMu22oIsoTkMu22_2016(row.m1Pt,abs(row.m1Eta))
  m2id = muon_pog_PFTight_2016(row.m2Pt,abs(row.m2Eta))
  m2iso = muon_pog_TightIso_2016('Tight',row.m2Pt,abs(row.m2Eta))
  m2_trg = muon_pog_IsoMu22oIsoTkMu22_2016(row.m2Pt,abs(row.m2Eta))

  #print "pu"
  #print str(pu)
  return pu*m1id*m1iso*m1_trg*m2id*m2iso*m2_trg

mc_corrector = mc_corrector_2016

class AnalyzeLFVMuMuTau(MegaBase):
    tree = 'mmt/final/Ntuple'
    #tree = 'New_Tree'

    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeLFVMuMuTau, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = MuMuTauTree.MuMuTauTree(tree)
        self.out = outfile
        self.histograms = {}

    def begin(self):

        names=["preselectionSS", "preselectionDecay0","preselectionLooseIsoDecay0", "preselectionVLooseIsoDecay0","preselectionVTightIsoDecay0","preselectionMediumIsoDecay0", "preselection0JetDecay0", "preselection1JetDecay0", "preselection2JetDecay0","preselectionDecay1","preselectionLooseIsoDecay1", "preselectionVLooseIsoDecay1","preselectionVTightIsoDecay1","preselectionMediumIsoDecay1", "preselection0JetDecay1", "preselection1JetDecay1", "preselection2JetDecay1","preselectionDecay10","preselectionLooseIsoDecay10", "preselectionVLooseIsoDecay10","preselectionVTightIsoDecay10","preselectionMediumIsoDecay10", "preselection0JetDecay10", "preselection1JetDecay10", "preselection2JetDecay10"]
        namesize = len(names)
	for x in range(0,namesize):


            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "GenWeight", "Gen level weight", 200000 ,-1000000, 1000000)
            self.book(names[x], "genHTT", "genHTT", 1000 ,0,1000)
 
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 100, -0.5, 100.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)

            self.book(names[x], "m1Pt", "Muon  Pt", 300,0,300)
            self.book(names[x], "m1Eta", "Muon  eta", 100, -2.5, 2.5)
            self.book(names[x], "m1Charge", "Muon Charge", 5, -2, 2)
            self.book(names[x], "m2Pt", "Muon  Pt", 300,0,300)
            self.book(names[x], "m2Eta", "Muon  eta", 100, -2.5, 2.5)
            self.book(names[x], "m2Charge", "Muon Charge", 5, -2, 2)


            self.book(names[x], "tPt", "Tau  Pt", 300,0,300)
            self.book(names[x], "tEta", "Tau  eta", 100, -2.5, 2.5)
            self.book(names[x], "tMtToPfMet_type1", "Tau MT (PF Ty1)", 200, 0, 200)
            self.book(names[x], "tCharge", "Tau  Charge", 5, -2, 2)
	    self.book(names[x], "tJetPt", "Tau Jet Pt" , 500, 0 ,500)	    
            self.book(names[x], "tMass", "Tau  Mass", 1000, 0, 10)
            self.book(names[x], "tLeadTrackPt", "Tau  LeadTrackPt", 300,0,300)

		       

            #self.book(names[x], "tAgainstMuonLoose", "tAgainstMuonLoose", 2,-0.5,1.5)
            self.book(names[x], "tAgainstMuonLoose3", "tAgainstMuonLoose3", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonMedium", "tAgainstMuonMedium", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonTight", "tAgainstMuonTight", 2,-0.5,1.5)
            self.book(names[x], "tAgainstMuonTight3", "tAgainstMuonTight3", 2,-0.5,1.5)

            #self.book(names[x], "tAgainstMuonLooseMVA", "tAgainstMuonLooseMVA", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonMediumMVA", "tAgainstMuonMediumMVA", 2,-0.5,1.5)
            #self.book(names[x], "tAgainstMuonTightMVA", "tAgainstMuonTightMVA", 2,-0.5,1.5)

            self.book(names[x], "tDecayModeFinding", "tDecayModeFinding", 2,-0.5,1.5)
            self.book(names[x], "tDecayModeFindingNewDMs", "tDecayModeFindingNewDMs", 2,-0.5,1.5)
            self.book(names[x], "tDecayMode", "tDecayMode", 21,-0.5,20.5)


            self.book(names[x], "tByLooseCombinedIsolationDeltaBetaCorr3Hits", "tByLooseCombinedIsolationDeltaBetaCorr3Hits", 2,-0.5,1.5)
            self.book(names[x], "tByMediumCombinedIsolationDeltaBetaCorr3Hits", "tByMediumCombinedIsolationDeltaBetaCorr3Hits", 2,-0.5,1.5)
            self.book(names[x], "tByTightCombinedIsolationDeltaBetaCorr3Hits", "tByTightCombinedIsolationDeltaBetaCorr3Hits", 2,-0.5,1.5)

            self.book(names[x], 'm1PixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
            self.book(names[x], 'm1JetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)
            self.book(names[x], 'm2PixHits', 'Mu 1 pix hits', 10, -0.5, 9.5)
            self.book(names[x], 'm2JetBtag', 'Mu 1 JetBtag', 100, -5.5, 9.5)

    	  
    	    self.book(names[x], "LT", "ht", 400, 0, 400)
            self.book(names[x], "type1_pfMetEt", "Type1 MET", 200, 0, 200)
    
            self.book(names[x], "m1_t_Mass", "Muon + Tau Mass", 200, 0, 200)
            self.book(names[x], "m1_t_Pt", "Muon + Tau Pt", 200, 0, 200)
            self.book(names[x], "m1_t_DR", "Muon + Tau DR", 100, 0, 10)
            self.book(names[x], "m1_t_DPhi", "Muon + Tau DPhi", 100, 0, 4)
            self.book(names[x], "m1_t_SS", "Muon + Tau SS", 5, -2, 2)
            self.book(names[x], "m1_t_ToMETDPhi_Ty1", "Muon Tau DPhi to MET", 100, 0, 4)
            self.book(names[x], "m2_t_Mass", "Muon + Tau Mass", 200, 0, 200)
            self.book(names[x], "m2_t_Pt", "Muon + Tau Pt", 200, 0, 200)
            self.book(names[x], "m2_t_DR", "Muon + Tau DR", 100, 0, 10)
            self.book(names[x], "m2_t_DPhi", "Muon + Tau DPhi", 100, 0, 4)
            self.book(names[x], "m2_t_SS", "Muon + Tau SS", 5, -2, 2)
            self.book(names[x], "m2_t_ToMETDPhi_Ty1", "Muon Tau DPhi to MET", 100, 0, 4)
            self.book(names[x], "m1_m2_Mass", "Dimuon Mass", 200, 0, 200)
    
            # Vetoes
            self.book(names[x], 'muVetoPt5IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
	    self.book(names[x], 'muVetoPt15IsoIdVtx', 'Number of extra muons', 5, -0.5, 4.5)
            self.book(names[x], 'tauVetoPt20Loose3HitsVtx', 'Number of extra taus', 5, -0.5, 4.5)
            self.book(names[x], 'eVetoMVAIso', 'Number of extra CiC tight electrons', 5, -0.5, 4.5)
   
            #self.book(names[x], 'jetVeto30PUCleanedTight', 'Number of extra jets', 5, -0.5, 4.5)
            #self.book(names[x], 'jetVeto30PUCleanedLoose', 'Number of extra jets', 5, -0.5, 4.5)
            self.book(names[x], 'jetVeto30', 'Number of extra jets', 5, -0.5, 4.5)	
            self.book(names[x], 'jetVeto30ZTT', 'Number of extra jets', 5, -0.5, 4.5)
            self.book(names[x], 'jetVeto30Eta3', 'Number of extra jets within |eta| < 3', 5, -0.5, 4.5)
	    #Isolation
	    self.book(names[x], 'm1RelPFIsoDBDefault' ,'Muon Isolation', 100, 0.0,1.0)
            self.book(names[x], 'm2RelPFIsoDBDefault' ,'Muon Isolation', 100, 0.0,1.0)
   
 
            self.book(names[x], "m1PhiMtPhi", "", 100, 0,4)
            self.book(names[x], "m1PhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "m2PhiMtPhi", "", 100, 0,4)
            self.book(names[x], "m2PhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "tPhiMETPhiType1", "", 100, 0,4)

### vbf ###
            self.book(names[x], "vbfJetVeto30", "central jet veto for vbf", 5, -0.5, 4.5)
	    self.book(names[x], "vbfJetVeto20", "", 5, -0.5, 4.5)
	    self.book(names[x], "vbfMVA", "", 100, 0,0.5)
	    self.book(names[x], "vbfMass", "", 500,0,5000.0)
	    self.book(names[x], "vbfDeta", "", 100, -0.5,10.0)
            self.book(names[x], "vbfMassZTT", "", 500,0,5000.0)
            self.book(names[x], "vbfDetaZTT", "", 100, -0.5,10.0)
            self.book(names[x], "vbfj1eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfj2eta","",100,-2.5,2.5)
	    self.book(names[x], "vbfVispt","",100,0,200)
	    self.book(names[x], "vbfHrap","",100,0,5.0)
	    self.book(names[x], "vbfDijetrap","",100,0,5.0)
	    self.book(names[x], "vbfDphihj","",100,0,4)
            self.book(names[x], "vbfDphihjnomet","",100,0,4)
            self.book(names[x], "vbfNJets", "g", 5, -0.5, 4.5)
            #self.book(names[x], "vbfNJetsPULoose", "g", 5, -0.5, 4.5)
            #self.book(names[x], "vbfNJetsPUTight", "g", 5, -0.5, 4.5)

    def correction(self,row):
	return mc_corrector(row)
	
    def fakeRateMethod(self,row,isoName):
        return getFakeRateFactor(row,isoName)
	     
    def fill_histos(self, row,name='gg', fakeRate=False, isoName="old"):
        histos = self.histograms
        weight=1
        if (not(data)):
	   weight = row.GenWeight * self.correction(row) #apply gen and pu reweighting to MC
        if (fakeRate == True):
          weight=weight*self.fakeRateMethod(row,isoName) #apply fakerate method for given isolation definition


        histos[name+'/weight'].Fill(weight)
        histos[name+'/GenWeight'].Fill(row.GenWeight)
        histos[name+'/genHTT'].Fill(row.genHTT)
        histos[name+'/rho'].Fill(row.rho, weight)
        histos[name+'/nvtx'].Fill(row.nvtx, weight)
        histos[name+'/prescale'].Fill(row.doubleMuPrescale, weight)

        
        histos[name+'/m1Pt'].Fill(row.m1Pt, weight)
        histos[name+'/m1Eta'].Fill(row.m1Eta, weight)
        histos[name+'/m1Charge'].Fill(row.m1Charge, weight)
        histos[name+'/m2Pt'].Fill(row.m2Pt, weight)
        histos[name+'/m2Eta'].Fill(row.m2Eta, weight)
        histos[name+'/m2Charge'].Fill(row.m2Charge, weight)
        histos[name+'/tPt'].Fill(row.tPt, weight)
        histos[name+'/tEta'].Fill(row.tEta, weight)
        histos[name+'/tMtToPfMet_type1'].Fill(row.tMtToPfMet_type1,weight)
        histos[name+'/tCharge'].Fill(row.tCharge, weight)
	histos[name+'/tJetPt'].Fill(row.tJetPt, weight)

        histos[name+'/tMass'].Fill(row.tMass,weight)
        histos[name+'/tLeadTrackPt'].Fill(row.tLeadTrackPt,weight)
		       

        #histos[name+'/tAgainstMuonLoose'].Fill(row.tAgainstMuonLoose,weight)
        histos[name+'/tAgainstMuonLoose3'].Fill(row.tAgainstMuonLoose3,weight)
        #histos[name+'/tAgainstMuonMedium'].Fill(row.tAgainstMuonMedium,weight)
        #histos[name+'/tAgainstMuonTight'].Fill(row.tAgainstMuonTight,weight)
        histos[name+'/tAgainstMuonTight3'].Fill(row.tAgainstMuonTight3,weight)

        #histos[name+'/tAgainstMuonLooseMVA'].Fill(row.tAgainstMuonLooseMVA,weight)
        #histos[name+'/tAgainstMuonMediumMVA'].Fill(row.tAgainstMuonMediumMVA,weight)
        #histos[name+'/tAgainstMuonTightMVA'].Fill(row.tAgainstMuonTightMVA,weight)

        histos[name+'/tDecayModeFinding'].Fill(row.tDecayModeFinding,weight)
        histos[name+'/tDecayModeFindingNewDMs'].Fill(row.tDecayModeFindingNewDMs,weight)
        histos[name+'/tDecayMode'].Fill(row.tDecayMode,weight)

        histos[name+'/tByLooseCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByLooseCombinedIsolationDeltaBetaCorr3Hits,weight)
        histos[name+'/tByMediumCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByMediumCombinedIsolationDeltaBetaCorr3Hits,weight)
        histos[name+'/tByTightCombinedIsolationDeltaBetaCorr3Hits'].Fill(row.tByTightCombinedIsolationDeltaBetaCorr3Hits,weight)


	histos[name+'/LT'].Fill(row.LT,weight)


	histos[name+'/type1_pfMetEt'].Fill(row.type1_pfMetEt,weight)

        histos[name+'/m1_t_Mass'].Fill(row.m1_t_Mass,weight)
        histos[name+'/m1_t_Pt'].Fill(row.m1_t_Pt,weight)
        histos[name+'/m1_t_DR'].Fill(row.m1_t_DR,weight)
        histos[name+'/m1_t_DPhi'].Fill(row.m1_t_DPhi,weight)
        histos[name+'/m1_t_SS'].Fill(row.m1_t_SS,weight)
        histos[name+'/m2_t_Mass'].Fill(row.m2_t_Mass,weight)
        histos[name+'/m2_t_Pt'].Fill(row.m2_t_Pt,weight)
        histos[name+'/m2_t_DR'].Fill(row.m2_t_DR,weight)
        histos[name+'/m2_t_DPhi'].Fill(row.m2_t_DPhi,weight)
        histos[name+'/m2_t_SS'].Fill(row.m2_t_SS,weight)
        histos[name+'/m1_m2_Mass'].Fill(row.m1_m2_Mass,weight)
	#histos[name+'/m_t_ToMETDPhi_Ty1'].Fill(row.m_t_ToMETDPhi_Ty1,weight)

        histos[name+'/m1PixHits'].Fill(row.m1PixHits, weight)
        histos[name+'/m1JetBtag'].Fill(row.m1JetBtag, weight)
        histos[name+'/m2PixHits'].Fill(row.m2PixHits, weight)
        histos[name+'/m2JetBtag'].Fill(row.m2JetBtag, weight)

        histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
        histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
        histos[name+'/tauVetoPt20Loose3HitsVtx'].Fill(row.tauVetoPt20Loose3HitsVtx, weight)
        histos[name+'/eVetoMVAIso'].Fill(row.eVetoMVAIso, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)
	histos[name+'/jetVeto30ZTT'].Fill(row.jetVeto30ZTT,weight)
        histos[name+'/jetVeto30Eta3'].Fill(row.jetVeto30Eta3,weight)
        #histos[name+'/jetVeto30PUCleanedLoose'].Fill(row.jetVeto30PUCleanedLoose, weight)
        #histos[name+'/jetVeto30PUCleanedTight'].Fill(row.jetVeto30PUCleanedTight, weight)

	histos[name+'/m1RelPFIsoDBDefault'].Fill(row.m1RelPFIsoDBDefault, weight)
        histos[name+'/m2RelPFIsoDBDefault'].Fill(row.m2RelPFIsoDBDefault, weight)
        
	histos[name+'/m1PhiMtPhi'].Fill(deltaPhi(row.m1Phi,row.tPhi),weight)
        histos[name+'/m1PhiMETPhiType1'].Fill(deltaPhi(row.m1Phi,row.type1_pfMetPhi),weight)
        histos[name+'/m2PhiMtPhi'].Fill(deltaPhi(row.m2Phi,row.tPhi),weight)
        histos[name+'/m2PhiMETPhiType1'].Fill(deltaPhi(row.m2Phi,row.type1_pfMetPhi),weight)
        histos[name+'/tPhiMETPhiType1'].Fill(deltaPhi(row.tPhi,row.type1_pfMetPhi),weight)
	histos[name+'/tDecayMode'].Fill(row.tDecayMode, weight)
	histos[name+'/vbfJetVeto30'].Fill(row.vbfJetVeto30, weight)
     	#histos[name+'/vbfJetVeto20'].Fill(row.vbfJetVeto20, weight)
        #histos[name+'/vbfMVA'].Fill(row.vbfMVA, weight)
        histos[name+'/vbfMass'].Fill(row.vbfMass, weight)
        histos[name+'/vbfDeta'].Fill(row.vbfDeta, weight)
        histos[name+'/vbfMassZTT'].Fill(row.vbfMassZTT, weight)
        histos[name+'/vbfDetaZTT'].Fill(row.vbfDetaZTT, weight)
        #histos[name+'/vbfj1eta'].Fill(row.vbfj1eta, weight)
        #histos[name+'/vbfj2eta'].Fill(row.vbfj2eta, weight)
        #histos[name+'/vbfVispt'].Fill(row.vbfVispt, weight)
        #histos[name+'/vbfHrap'].Fill(row.vbfHrap, weight)
        #histos[name+'/vbfDijetrap'].Fill(row.vbfDijetrap, weight)
        #histos[name+'/vbfDphihj'].Fill(row.vbfDphihj, weight)
        #histos[name+'/vbfDphihjnomet'].Fill(row.vbfDphihjnomet, weight)
        histos[name+'/vbfNJets'].Fill(row.vbfNJets, weight)
        #histos[name+'/vbfNJetsPULoose'].Fill(row.vbfNJetsPULoose, weight)
        #histos[name+'/vbfNJetsPUTight'].Fill(row.vbfNJetsPUTight, weight)




    def presel(self, row):
        if not (row.singleIsoMu22Pass or row.singleIsoTkMu22Pass):
            return False
        return True

    def selectZtt(self,row):
        if (ZTauTau and not row.isZtautau):
            return False
        if (not ZTauTau and row.isZtautau):
            return False
        return True
 
    def selectZeroJet(self,row):
	if (ZeroJet and row.NUP != 5):
            return False
	return True

    def kinematics(self, row):
        if row.m1Pt < 30:
            return False
        if abs(row.m1Eta) >= 2.3:
            return False
        if row.m2Pt < 30:
            return False
        if abs(row.m2Eta) >= 2.3:
            return False
        if row.tPt<20 :
            return False
        if abs(row.tEta)>=2.3:
            return False
        return True

    def gg(self,row):
       if row.mPt < 25:    
           return False
       if deltaPhi(row.mPhi, row.tPhi) <2.7:
           return False
       if deltaPhi(row.tPhi,row.type1_pfMetPhi) > 3.0:
	   return False
       if row.tPt < 30:
           return False
       if row.tMtToPfMet_type1 > 75:
           return False
       if row.jetVeto30ZTT!=0:
           return False
       return True

    def boost(self,row):
          if row.jetVeto30ZTT!=1:
            return False
          if row.mPt < 25:
                return False
          if row.tPt < 30:
                return False
          if row.tMtToPfMet_type1 > 105:
                return False
          if deltaPhi(row.tPhi,row.type1_pfMetPhi) > 3.0:
                return False

          return True

    def vbfAntiTight(self,row):
        if row.tPt < 30:
                return False
        if row.mPt < 25:
		return False
        if row.tMtToPfMet_type1 > 75:
                return False
        if row.jetVeto30ZTT < 2:
            return False
	if(row.vbfNJets<2):
	    return False
	if(abs(row.vbfDetaZTT)>3.5):
	    return False
        if row.vbfMassZTT > 500:
	    return False
        if row.vbfJetVeto30 > 0:
            return False
        return True

    def vbf(self,row):
        if row.tPt < 30:
                return False
        if row.mPt < 25:
                return False
        if row.tMtToPfMet_type1 > 75:
                return False
        if row.jetVeto30ZTT < 2:
            return False
        if(row.vbfNJets<2):
            return False
        if(abs(row.vbfDetaZTT)<3.2):
            return False
        if row.vbfMassZTT < 500:
            return False
        if row.vbfJetVeto30 > 0:
            return False
        return True


    def oppositesign(self,row):
	if row.m1Charge*row.m2Charge!=-1:
            return False
	return True

    #def obj1_id(self, row):
    #    return bool(row.mPFIDTight)  
 
    def obj1_idICHEP(self,row):
	 if (row.m1IsGlobal and (row.m1NormalizedChi2 < 3) and (row.m1Chi2LocalPosition < 12) and (row.m1TrkKink < 20)):
		goodGlobal1=True
	 else:
		goodGlobal1=False
         if (row.m2IsGlobal and (row.m2NormalizedChi2 < 3) and (row.m2Chi2LocalPosition < 12) and (row.m2TrkKink < 20)):
                goodGlobal2=True
         else:
                goodGlobal2=False
    	 return ((row.m1PFIDLoose and row.m1ValidFraction > 0.49 and ((goodGlobal1 and row.m1SegmentCompatibility > 0.303) or row.m1SegmentCompatibility > 0.451)) and (row.m2PFIDLoose and row.m2ValidFraction > 0.49 and ((goodGlobal2 and row.m2SegmentCompatibility > 0.303) or row.m2SegmentCompatibility > 0.451)))

    def obj1_id(self,row):
         return (row.m1PixHits>0 and row.m2PixHits>0 and row.m1JetPFCISVBtag < 0.8 and row.m2JetPFCISVBtag < 0.8 and row.m1PVDZ < 0.2 and row.m2PVDZ < 0.2 and row.m1PFIDTight and row.m2PFIDTight)

    def m1m2Mass(self,row):
        if row.m1_m2_Mass < 70:
           return False
        if row.m1_m2_Mass > 110:
           return False
        return True

    def obj2_id(self, row):
	return  row.tAgainstElectronLooseMVA6 and row.tAgainstMuonTight3 and row.tDecayModeFinding and row.tPVDZ < 0.2 and row.tMuonIdIsoVtxOverlap == 0 and row.tElecOverlap == 0

    def vetos(self,row):
		return  ((row.eVetoZTTp001dxyzR0 == 0) and (row.dimuonVeto==1) and (row.muVetoZTTp001dxyzR0 < 3) and (row.tauVetoPt20Loose3HitsVtx<1) )

    #def obj1_iso(self, row):
    #    return bool(row.mRelPFIsoDBDefault <0.12)
   
    def obj1_iso(self,row):
         return bool(row.m1RelPFIsoDBDefault <0.25 and row.m2RelPFIsoDBDefault <0.25)

    def obj2_iso(self,row):
         return bool(row.tByTightIsolationMVArun2v1DBoldDMwLT)

    #def obj2_iso(self, row):
    #    return  row.tByTightCombinedIsolationDeltaBetaCorr3Hits

    def obj2_mediso(self, row):
	 return bool(row.tByMediumIsolationMVArun2v1DBoldDMwLT)

    def obj2_vtightiso(self,row):
         return bool(row.tByVTightIsolationMVArun2v1DBoldDMwLT)

    def obj2_looseiso(self,row):
         return bool(row.tByLooseIsolationMVArun2v1DBoldDMwLT)

    def obj2_vlooseiso(self,row):
         return bool(row.tByVLooseIsolationMVArun2v1DBoldDMwLT)

    def obj1_antiiso(self, row):
        return bool(row.m1RelPFIsoDBDefault >0.2 and row.m2RelPFIsoDBDefault >0.2) 

        return row.tByLooseCombinedIsolationDeltaBetaCorr3Hits

    def DecayMode0(self,row):
      return bool(row.tDecayMode==0)

    def DecayMode1(self,row):
      return bool(row.tDecayMode==1)

    def DecayMode10(self,row):
      return bool(row.tDecayMode==10)


    def process(self):
        event =0
        sel=False
        for row in self.tree:
            if event!=row.evt:   # This is just to ensure we get the (Mu,Tau) with the highest Pt
                event=row.evt    # In principle the code saves all the MU+Tau posibilities, if an event has several combinations
                sel = False      # it will save them all.
            if sel==True:
                continue

            if data and not self.presel(row): #only apply trigger selections for data
                continue
            if not self.selectZtt(row):
                continue

            if not self.kinematics(row): 
                continue
 
            if not self.obj1_iso(row):
                continue
            if not self.obj1_idICHEP(row):
                continue

            if not self.vetos (row):
                continue
            if not self.m1m2Mass (row):
                continue

            if not self.obj2_id (row):
                continue

            if self.obj2_iso(row) and not self.oppositesign(row):
              self.fill_histos(row,'preselectionSS',False)

            if self.DecayMode0(row):           
              if self.obj2_iso(row) and self.oppositesign(row):  
 
                self.fill_histos(row,'preselectionDecay0',False)
                if row.jetVeto30ZTT==0:
                  self.fill_histos(row,'preselection0JetDecay0',False)
                if row.jetVeto30ZTT==1:
                  self.fill_histos(row,'preselection1JetDecay0',False)
                if row.jetVeto30ZTT==2:
                  self.fill_histos(row,'preselection2JetDecay0',False)
 
              if self.obj2_looseiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionLooseIsoDecay0',False)
        
              if self.obj2_mediso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionMediumIsoDecay0',False)
 
              if self.obj2_vlooseiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionVLooseIsoDecay0',False)
 
              if self.obj2_vtightiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionVTightIsoDecay0',False)

            if self.DecayMode1(row):
              if self.obj2_iso(row) and self.oppositesign(row):

                self.fill_histos(row,'preselectionDecay1',False)
                if row.jetVeto30ZTT==0:
                  self.fill_histos(row,'preselection0JetDecay1',False)
                if row.jetVeto30ZTT==1:
                  self.fill_histos(row,'preselection1JetDecay1',False)
                if row.jetVeto30ZTT==2:
                  self.fill_histos(row,'preselection2JetDecay1',False)

              if self.obj2_looseiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionLooseIsoDecay1',False)
  
              if self.obj2_mediso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionMediumIsoDecay1',False)

              if self.obj2_vlooseiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionVLooseIsoDecay1',False)

              if self.obj2_vtightiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionVTightIsoDecay1',False)

            if self.DecayMode10(row):
              if self.obj2_iso(row) and self.oppositesign(row):

                self.fill_histos(row,'preselectionDecay10',False)
                if row.jetVeto30ZTT==0:
                  self.fill_histos(row,'preselection0JetDecay10',False)
                if row.jetVeto30ZTT==1:
                  self.fill_histos(row,'preselection1JetDecay10',False)
                if row.jetVeto30ZTT==2:
                  self.fill_histos(row,'preselection2JetDecay10',False)

              if self.obj2_looseiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionLooseIsoDecay10',False)
  
              if self.obj2_mediso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionMediumIsoDecay10',False)

              if self.obj2_vlooseiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionVLooseIsoDecay10',False)

              if self.obj2_vtightiso(row) and self.oppositesign(row):
                self.fill_histos(row,'preselectionVTightIsoDecay10',False)



            sel=True

    def finish(self):
        self.write_histos()
