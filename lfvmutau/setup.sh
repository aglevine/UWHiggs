#!/bin/bash

# Get the data
export datasrc=/hdfs/store/user/aglevine/
#export datasrc=/hdfs/store/user/cepeda/
#export datasrc=/hdfs/store/user/taroni/
#export jobid=MiniAODSIM-Spring15-25ns_LFV_V1_October10
#export jobid=MiniAODSIM-Spring15-25ns_LFV_MiniAODV2_Nov3
#export jobid=MiniAODSIMv2-Spring15-25ns_LFV_October13
#export jobid=MiniAODv2_2fb_v3
#export jobid=MiniAODSIM-Spring15-25ns_LFV_MiniAODV2_Dec2_LFV_NoHF_JetEta25
#export jobid=MiniAODSIM-Spring15-25ns_LFV_MiniAODV2_Dec2_Data_NoHF_JetEta
#export jobid=MiniAODSIM-Spring15-25ns_LFV_MiniAODV2_Dec2_LFV_NoHF_JetEta25_MissingHiggs
#export jobid=MiniAodV2For25ns_ExtraJets_LFV_Data
#export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetEtaFix_LFV
export jobid=MiniAodV2For25ns_ExtraJets_JesUes_JetFix_LFVData
#export jobid=MiniAODv2_2fb_v2
#export jobid=MiniAODSIM-Spring15-25ns_LFV_MiniAODV2_Data
#export jobid=MiniAODSIM-Spring15-25ns_LFV_Nov9_ZTTFakeStudy
#export jobid=MiniAODSIM-Spring15-25ns_LFV_Oct27
#export afile=`find $datasrc/$jobid | grep root | head -n 1`

## Build the cython wrappers
rake "make_wrapper[$afile, mt/final/Ntuple, MuTauTree]"

ls *pyx | sed "s|pyx|so|" | xargs rake 
#echo "finishing compilation" 
#bash compileTree.txt

rake "meta:getinputs[$jobid, $datasrc,mt/metaInfo, mt/summedWeights]"
rake "meta:getmeta[inputs/$jobid, mt/metaInfo, 13]"
