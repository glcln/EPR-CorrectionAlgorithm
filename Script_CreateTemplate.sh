#!/bin/bash

##### Template for Muons

root -l <<EOC
.L CreateTemplate.C
TChain *chain = new TChain("stage/ttree")
chain->Add("ROOT_Files/nt_mc_aod_eta2p1_wPU_ext.root")
chain->Add("ROOT_Files/nt_mc_aod_eta2p1_eta3_wPU_ext.root")
CreateTemplate t(chain)
t.Loop()
EOC