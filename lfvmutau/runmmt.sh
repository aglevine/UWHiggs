#!/bin/bash

#script used for running AnalyzeMuMuTauTight
set -o nounset
set -o errexit

export Ztt=false
export embed=false
export PU=true
export iso=true
export jobid=mmt
export iso=true
export selection="preselection"
rake mmttight
export iso=false
rake mmttight
#export iso=true
#export selection="signal"
#rake mmttight
#export iso=false
#rake mmttight

