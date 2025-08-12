#!/bin/bash

root -l <<EOC
.L TestNewCorr.C
TChain *chain = new TChain("stage/ttree")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodMarch2023_CMSSW_10_6_30/hadd_WJetsToLNu_ntv2/prodMarch2023_CMSSW_10_6_30000*.root")
TestNewCorr t(chain)
t.Loop()
EOC

