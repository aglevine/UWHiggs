#!/bin/bash

set -o nounset
set -o errexit

#python MakeSysAnalyzers.py jesup
#python MakeSysAnalyzers.py jesdown
#python MakeSysAnalyzers.py uesup
#python MakeSysAnalyzers.py uesdown
#python MakeSysAnalyzers.py tesdown
#python MakeSysAnalyzers.py tesup

export systematic=none
export jobid=LFV_sep16_v2
#export jobid=SMHTT_aug16_v2
export isRealData=false
export isZTauTau=false
export isInclusive=false
#rake analyzeSpring2015Misc
#rake analyzeSpring2015WJets
#rake analyze2016mmtZJets
#rake analyze2016eetZJets
rake analyze2016mmmZJets
#rake analyze2016mmmDiBoson
#rake analyzeSpring2015DYAmcJets
#rake analyzeSpring2015TTBar
#rake analyze2016Higgs
export isZTauTau=true
#rake analyzeSpring2015ZTauTauJets
export isIncluse=false
#rake analyzeSpring2015ZTauTauAmcJets
export isZTauTau=false
#rake analyzeSpring2015DYZeroJets
export isInclusive=false
export isRealData=true
#rake analyzeLFVMuTauData
#rake analyze2016mmmDataSingleMuon
