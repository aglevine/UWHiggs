#!/bin/bash

# Analyze the H->mutau channel

set -o nounset
set -o errexit

#export jobid=SubmitMuTauSingleMU
#export jobid=MuTauSingleMUV6
#export jobid=MuTauSingleMuJetReReco
#export jobid=MuTauSingleMuJetReRecoSkim
export otherJets=false
export ZjetsM50=false
export treeSkim=false
export jetcorrection=""
export jobid=SingleMu
#export jobid=MuTauSingleMuJetReRecoJes
#export jobid=dataJes
export PU=true
export iso=true
export selection="preselection"
export optimized=true
export wjets=false
export embed=false
export Ztt=false
export wjets=true
#rake mttightwjets
export wjets=false
export embed=true
#rake embedmttight
export embed=false
export Ztt=true
#rake DYmttight
export ZjetsM50=true
#rake mttightZjetsM50
export otherJets=true
#rake mttightotherM50
#rake mttightZjetsM50
export ZjetsM50=false
#rake mttightotherDY
#rake DYmttight
export Ztt=false
#rake dataonlymttightpfmet
#rake dataonlymttight
#rake mttight
#rake mttightRedoFakes
#rake mttightsingletop
#rake mttightsm
#rake mttightwjets
#rake mttightnewsignal
#export PU=false
#rake pu1mttight
export jobid=MuTauSingleMuJetReReco
#export jobid=MuTauSingleMuJetReRecoJes
#export jobid=dataJes
export otherJets=false
export ZjetsM50=false
export treeSkim=false
export jetcorrection=""
#export jobid=MuTauSingleMuJetReRecoJes
#export jobid=dataJes
export PU=true
export iso=true
export selection="signal"
export optimized=true
export wjets=false
export embed=false
export Ztt=false
export wjets=true
#rake mttightwjets
export wjets=false
export embed=true
#rake embedmttight
export embed=false
export Ztt=true
#rake DYmttight
export ZjetsM50=true
#rake mttightZjetsM50
export otherJets=true
#rake mttightotherM50
export ZjetsM50=false
#rake mttightotherDY
export Ztt=false
#rake mttightextra
rake dataonlymttightpfmet
#rake mttight
#rake mttightHWW 
#rake dataonlymttight
#rake mttightRedoFakes
#rake mttightsingletop
#rake mttightsm
#rake mttightwjets
#rake mttightnewsignal
#export PU=false
#rake pu1mttight

