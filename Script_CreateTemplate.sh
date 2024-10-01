#!/bin/bash

##### Template for Muons

root -l <<EOC
.L CreateTemplate.C
TChain *chain = new TChain("stage/ttree")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodJan2024_CMSSW_10_6_30/nt_mc_aod_eta2p1_eta3_wPU_ext1_v2.root")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodJan2024_CMSSW_10_6_30/ntuples_local_eta2p1_v2/nt_mc_aod_*.root")
chain->Add("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodJan2024_CMSSW_10_6_30/nt_mc_aod_eta2p1_eta3_wPU_ext2_v2.root")
CreateTemplate t(chain)
t.Loop()
EOC