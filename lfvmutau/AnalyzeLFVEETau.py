'''

Run LFV H->EETau analysis in the mu+tau channel.

Authors: Maria Cepeda, Aaron Levine, Evan K. Friis, UW

'''

import EETauTree
from FinalStateAnalysis.PlotTools.MegaBase import MegaBase
import glob
import os
import FinalStateAnalysis.TagAndProbe.EGammaPOGCorrections as EGammaPOGCorrections
import FinalStateAnalysis.TagAndProbe.HetauCorrection as HetauCorrections
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

egammaWP80_pog_2016 = EGammaPOGCorrections.make_egamma_pog_electronID_ICHEP2016('nontrigWP80')

def mc_corrector_2016(row):
  pu = pu_corrector(row.nTruePU)

  e1id = egammaWP80_pog_2016(abs(row.e1Eta),row.e1Pt)
  e2id = egammaWP80_pog_2016(abs(row.e2Eta),row.e2Pt)
  e1trg = HetauCorrections.single_ele_2016(row.e1Pt,abs(row.e1Eta))[1]
  e2trg = HetauCorrections.single_ele_2016(row.e2Pt,abs(row.e2Eta))[1]
  #print "pu"
  #print str(pu)
  return pu*e1id*e2id*e1trg*e2trg

mc_corrector = mc_corrector_2016

class AnalyzeLFVEETau(MegaBase):
    tree = 'eet/final/Ntuple'
    #tree = 'New_Tree'

    def __init__(self, tree, outfile, **kwargs):
        super(AnalyzeLFVEETau, self).__init__(tree, outfile, **kwargs)
        # Use the cython wrapper
        self.tree = EETauTree.EETauTree(tree)
        self.out = outfile
        self.histograms = {}

    def begin(self):

        names=["preselection","preselectionSS", "preselectionLooseIso", "preselectionVLooseIso","preselectionVTightIso","preselectionMediumIso", "preselection0Jet", "preselection1Jet", "preselection2Jet"]
        namesize = len(names)
	for x in range(0,namesize):


            self.book(names[x], "weight", "Event weight", 100, 0, 5)
            self.book(names[x], "GenWeight", "Gen level weight", 200000 ,-1000000, 1000000)
            self.book(names[x], "genHTT", "genHTT", 1000 ,0,1000)
 
            self.book(names[x], "rho", "Fastjet #rho", 100, 0, 25)
            self.book(names[x], "nvtx", "Number of vertices", 100, -0.5, 100.5)
            self.book(names[x], "prescale", "HLT prescale", 21, -0.5, 20.5)

            self.book(names[x], "e1Pt", "Electron  Pt", 300,0,300)
            self.book(names[x], "e1Eta", "Electron  eta", 100, -2.5, 2.5)
            self.book(names[x], "e1Charge", "Electron Charge", 5, -2, 2)
            self.book(names[x], "e2Pt", "Electron  Pt", 300,0,300)
            self.book(names[x], "e2Eta", "Electron  eta", 100, -2.5, 2.5)
            self.book(names[x], "e2Charge", "Electron Charge", 5, -2, 2)


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

    	  
    	    self.book(names[x], "LT", "ht", 400, 0, 400)
            self.book(names[x], "type1_pfMetEt", "Type1 MET", 200, 0, 200)
    
            self.book(names[x], "e1_t_Mass", "Electron + Tau Mass", 200, 0, 200)
            self.book(names[x], "e1_t_Pt", "Electron + Tau Pt", 200, 0, 200)
            self.book(names[x], "e1_t_DR", "Electron + Tau DR", 100, 0, 10)
            self.book(names[x], "e1_t_DPhi", "Electron + Tau DPhi", 100, 0, 4)
            self.book(names[x], "e1_t_SS", "Electron + Tau SS", 5, -2, 2)
            self.book(names[x], "e1_t_ToMETDPhi_Ty1", "Electron Tau DPhi to MET", 100, 0, 4)
            self.book(names[x], "e2_t_Mass", "Electron + Tau Mass", 200, 0, 200)
            self.book(names[x], "e2_t_Pt", "Electron + Tau Pt", 200, 0, 200)
            self.book(names[x], "e2_t_DR", "Electron + Tau DR", 100, 0, 10)
            self.book(names[x], "e2_t_DPhi", "Electron + Tau DPhi", 100, 0, 4)
            self.book(names[x], "e2_t_SS", "Electron + Tau SS", 5, -2, 2)
            self.book(names[x], "e2_t_ToMETDPhi_Ty1", "Electron Tau DPhi to MET", 100, 0, 4)
            self.book(names[x], "e1_e2_Mass", "Dielectron Mass", 200, 0, 200)
    
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
	    self.book(names[x], 'e1RelIso' ,'Electron Isolation', 100, 0.0,1.0)
            self.book(names[x], 'e2RelIso' ,'Electron Isolation', 100, 0.0,1.0)
   
 
            self.book(names[x], "e1PhiMtPhi", "", 100, 0,4)
            self.book(names[x], "e1PhiMETPhiType1", "", 100, 0,4)
            self.book(names[x], "e2PhiMtPhi", "", 100, 0,4)
            self.book(names[x], "e2PhiMETPhiType1", "", 100, 0,4)
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

        
        histos[name+'/e1Pt'].Fill(row.e1Pt, weight)
        histos[name+'/e1Eta'].Fill(row.e1Eta, weight)
        histos[name+'/e1Charge'].Fill(row.e1Charge, weight)
        histos[name+'/e2Pt'].Fill(row.e2Pt, weight)
        histos[name+'/e2Eta'].Fill(row.e2Eta, weight)
        histos[name+'/e2Charge'].Fill(row.e2Charge, weight)
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

        histos[name+'/e1_t_Mass'].Fill(row.e1_t_Mass,weight)
        histos[name+'/e1_t_Pt'].Fill(row.e1_t_Pt,weight)
        histos[name+'/e1_t_DR'].Fill(row.e1_t_DR,weight)
        histos[name+'/e1_t_DPhi'].Fill(row.e1_t_DPhi,weight)
        histos[name+'/e1_t_SS'].Fill(row.e1_t_SS,weight)
        histos[name+'/e2_t_Mass'].Fill(row.e2_t_Mass,weight)
        histos[name+'/e2_t_Pt'].Fill(row.e2_t_Pt,weight)
        histos[name+'/e2_t_DR'].Fill(row.e2_t_DR,weight)
        histos[name+'/e2_t_DPhi'].Fill(row.e2_t_DPhi,weight)
        histos[name+'/e2_t_SS'].Fill(row.e2_t_SS,weight)
        histos[name+'/e1_e2_Mass'].Fill(row.e1_e2_Mass,weight)
	#histos[name+'/m_t_ToMETDPhi_Ty1'].Fill(row.m_t_ToMETDPhi_Ty1,weight)


        histos[name+'/muVetoPt5IsoIdVtx'].Fill(row.muVetoPt5IsoIdVtx, weight)
        histos[name+'/muVetoPt15IsoIdVtx'].Fill(row.muVetoPt15IsoIdVtx, weight)
        histos[name+'/tauVetoPt20Loose3HitsVtx'].Fill(row.tauVetoPt20Loose3HitsVtx, weight)
        histos[name+'/eVetoMVAIso'].Fill(row.eVetoMVAIso, weight)
        histos[name+'/jetVeto30'].Fill(row.jetVeto30, weight)
	histos[name+'/jetVeto30ZTT'].Fill(row.jetVeto30ZTT,weight)
        histos[name+'/jetVeto30Eta3'].Fill(row.jetVeto30Eta3,weight)
        #histos[name+'/jetVeto30PUCleanedLoose'].Fill(row.jetVeto30PUCleanedLoose, weight)
        #histos[name+'/jetVeto30PUCleanedTight'].Fill(row.jetVeto30PUCleanedTight, weight)

	histos[name+'/e1RelIso'].Fill(row.e1RelIso, weight)
        histos[name+'/e2RelIso'].Fill(row.e2RelIso, weight)
        
	histos[name+'/e1PhiMtPhi'].Fill(deltaPhi(row.e1Phi,row.tPhi),weight)
        histos[name+'/e1PhiMETPhiType1'].Fill(deltaPhi(row.e1Phi,row.type1_pfMetPhi),weight)
        histos[name+'/e2PhiMtPhi'].Fill(deltaPhi(row.e2Phi,row.tPhi),weight)
        histos[name+'/e2PhiMETPhiType1'].Fill(deltaPhi(row.e2Phi,row.type1_pfMetPhi),weight)
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
        if not (row.singleE25eta2p1TightPass):
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
        if row.e1Pt < 24:
            return False
        if abs(row.e1Eta) >= 2.1 or (abs(row.e1Eta) > 1.4442 and abs(row.e1Eta) < 1.566):
            return False
        if row.e2Pt < 24:
            return False
        if abs(row.e2Eta) >= 2.1 or (abs(row.e2Eta) > 1.4442 and abs(row.e2Eta) < 1.566):
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
	if row.e1Charge*row.e2Charge!=-1:
            return False
	return True

    #def obj1_id(self, row):
    #    return bool(row.mPFIDTight)  
 

    def obj1_id(self,row):
         return (row.e1MVANonTrigWP80 and row.e2MVANonTrigWP80 and row.e1PassesConversionVeto and row.e2PassesConversionVeto and row.e1JetPFCISVBtag < 0.8 and row.e2JetPFCISVBtag < 0.8 and row.e1PVDZ < 0.2 and row.e2PVDZ < 0.2 and row.e1ChargeIdTight and row.e2ChargeIdTight)

    def e1e2Mass(self,row):
        if row.e1_e2_Mass < 70:
           return False
        if row.e1_e2_Mass > 110:
           return False
        return True

    def obj2_id(self, row):
	return  row.tAgainstElectronLooseMVA6 and row.tAgainstMuonTight3 and row.tDecayModeFinding and row.tPVDZ < 0.2 and row.tMuonIdIsoVtxOverlap == 0 and row.tElecOverlap == 0

    def vetos(self,row):
		return  ((row.eVetoZTTp001dxyzR0 < 3) and (row.dimuonVeto==0) and (row.muVetoZTTp001dxyzR0 ==0) and (row.tauVetoPt20Loose3HitsVtx<1) )

    #def obj1_iso(self, row):
    #    return bool(row.mRelPFIsoDBDefault <0.12)
   
    def obj1_iso(self,row):
         return bool(row.e1RelIso <0.1 and row.e2RelIso < 0.1)

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
        return bool(row.e1RelIso>0.1 and row.e2RelIso >0.1) 

        return row.tByLooseCombinedIsolationDeltaBetaCorr3Hits



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
            #print "passed data"
            if not self.selectZtt(row):
                continue
            #print "passed ztt"
            if not self.kinematics(row): 
                continue
            #print "passed kinematics"
            if not self.obj1_iso(row):
                continue
            #print "passed obj1iso"
            if not self.obj1_id(row):
                continue
            #print "passed obj1id"
            if not self.vetos (row):
                continue
            #print "passed vetos"
            if not self.e1e2Mass (row):
                continue

            if not self.obj2_id (row):
                continue

            if self.obj2_iso(row) and not self.oppositesign(row):
              self.fill_histos(row,'preselectionSS',False)


            if self.obj2_iso(row) and self.oppositesign(row):  
              #print "about to fill preselection"
              self.fill_histos(row,'preselection',False)
              if row.jetVeto30ZTT==0:
                self.fill_histos(row,'preselection0Jet',False)
              if row.jetVeto30ZTT==1:
                self.fill_histos(row,'preselection1Jet',False)
              if row.jetVeto30ZTT==2:
                self.fill_histos(row,'preselection2Jet',False)

            if self.obj2_looseiso(row) and self.oppositesign(row):
              self.fill_histos(row,'preselectionLooseIso',False)
      
            if self.obj2_mediso(row) and self.oppositesign(row):
              self.fill_histos(row,'preselectionMediumIso',False)

            if self.obj2_vlooseiso(row) and self.oppositesign(row):
              self.fill_histos(row,'preselectionVLooseIso',False)

            if self.obj2_vtightiso(row) and self.oppositesign(row):
              self.fill_histos(row,'preselectionVTightIso',False)



            sel=True

    def finish(self):
        self.write_histos()
