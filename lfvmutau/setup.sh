#!/bin/bash

# Get the data
#export datasrc=/hdfs/store/user/caillol/
export datasrc=/hdfs/store/user/taroni/
#export datasrc=/hdfs/store/user/ndev/
#export datasrc=/hdfs/store/user/cepeda/
#export datasrc=/hdfs/store/user/taroni/
#export jobid=SMHTT_aug16_v2
#export jobid=LFV_sep16_v2
export jobid=LFV_sep16_et
export afile=`find $datasrc/$jobid | grep root | head -n 1`

## Build the cython wrappers
rake "make_wrapper[$afile, eet/final/Ntuple, EETauTree]"

ls *pyx | sed "s|pyx|so|" | xargs rake 

rake "meta:getinputs[$jobid, $datasrc,eet/metaInfo, eet/summedWeights]"
rake "meta:getmeta[inputs/$jobid, eet/metaInfo, 13, eet/summedWeights]"
