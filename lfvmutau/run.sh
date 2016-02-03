#!/bin/bash

set -o nounset
set -o errexit

python MakeSysAnalyzers.py jesup
python MakeSysAnalyzers.py jesdown
python MakeSysAnalyzers.py uesup
python MakeSysAnalyzers.py uesdown

export systematic=none
export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetEtaFix_LFV
export isRealData=false
export isZTauTau=false
rake analyzeSpring2015Misc
rake analyzeSpring2015WJets
rake analyzeSpring2015ZJets
#rake analyzeSpring2015TTBar
export isZTauTau=true
rake analyzeSpring2015ZTauTauJets
export isZTauTau=false
export isRealData=true
rake analyzeLFVMuTauData


export systematic=jesup
export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetEtaFix_LFV
export isRealData=false
export isZTauTau=false
#rake analyzeSpring2015MiscJesUp
#rake analyzeSpring2015WJetsJesUp
#rake analyzeSpring2015ZJetsJesUp
#rake analyzeSpring2015TTBarJesUp
export isZTauTau=true
#rake analyzeSpring2015ZTauTauJetsJesUp
export isZTauTau=false
export isRealData=true
#rake analyzeLFVMuTauDataJesUp

export systematic=jesdown
export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetEtaFix_LFV
export isRealData=false
export isZTauTau=false
#rake analyzeSpring2015MiscJesDown
#rake analyzeSpring2015WJetsJesDown
#rake analyzeSpring2015ZJetsJesDown
#rake analyzeSpring2015TTBarJesDown
export isZTauTau=true
#rake analyzeSpring2015ZTauTauJetsJesDown
export isZTauTau=false
export isRealData=true
#rake analyzeLFVMuTauDataJesDown

export systematic=uesup
export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetEtaFix_LFV
export isRealData=false
export isZTauTau=false
#rake analyzeSpring2015MiscUesUp
#rake analyzeSpring2015WJetsUesUp
#rake analyzeSpring2015ZJetsUesUp
#rake analyzeSpring2015TTBarUesUp
export isZTauTau=true
#rake analyzeSpring2015ZTauTauJetsUesUp
export isZTauTau=false
export isRealData=true
#rake analyzeLFVMuTauDataUesUp

export systematic=uesdown
export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetEtaFix_LFV
export isRealData=false
export isZTauTau=false
#rake analyzeSpring2015MiscUesDown
#rake analyzeSpring2015WJetsUesDown
#rake analyzeSpring2015ZJetsUesDown
#rake analyzeSpring2015TTBarUesDown
export isZTauTau=true
#rake analyzeSpring2015ZTauTauJetsUesDown
export isZTauTau=false
export isRealData=true
#rake analyzeLFVMuTauDataUesDown
