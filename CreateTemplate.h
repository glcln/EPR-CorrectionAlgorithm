//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 28 16:47:08 2023 by ROOT version 6.26/10
// from TTree ttree/ttree
// found on file: nt_mc_aod.root
//////////////////////////////////////////////////////////

#ifndef CreateTemplate_h
#define CreateTemplate_h

#define nmax_gen 100
#define nmax_tr 1000
#define nmax_cl 10000
#define nmax_st 30000
#define nmax_simhit 1000
#define nmax_sigpix 30000
#define nmax_sumpix 30000
#define nmax_sigstrip 30000
#define nmax_sumstrip 30000
#define nmax_mu 200
#define nmax_hscp 200
#define nmax_jet 100


#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class CreateTemplate {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   UInt_t          event;
   Int_t           npv;
   Int_t           ngoodpv;
   Bool_t          goodIsFistPV;
   Float_t         pv_z;
   Float_t         pv_errz;
   Float_t         pv0_z;
   Float_t         pv0_errz;
   Bool_t          hlt_mu45;
   Bool_t          hlt_mu50;
   Bool_t          hlt_tkmu100;
   Bool_t          hlt_oldmu100;
   Bool_t          hlt_pfmet_mht;
   Bool_t          hlt_pfmet;
   Bool_t          hlt_tkmu50;
   Float_t         InstLumi;
   Float_t         gen_pv_z;
   Int_t           ngenpart;
   Int_t           gen_pdg[nmax_gen];   //[ngenpart]
   Float_t         gen_pt[nmax_gen];   //[ngenpart]
   Float_t         gen_eta[nmax_gen];   //[ngenpart]
   Float_t         gen_phi[nmax_gen];   //[ngenpart]
   Float_t         gen_mass[nmax_gen];   //[ngenpart]
   Bool_t          gen_isHardProcess[nmax_gen];   //[ngenpart]
   Int_t           gen_status[nmax_gen];   //[ngenpart]
   Int_t           gen_moth_pdg[nmax_gen];   //[ngenpart]
   Int_t           gen_ndaughter[nmax_gen];   //[ngenpart]
   Int_t           gen_daughter_pdg[nmax_gen];   //[ngenpart]
   Int_t           ntracks;
   Float_t         track_pt[nmax_tr];   //[ntracks]
   Float_t         track_pterr[nmax_tr];   //[ntracks]
   Float_t         track_p[nmax_tr];   //[ntracks]
   Float_t         track_eta[nmax_tr];   //[ntracks]
   Float_t         track_phi[nmax_tr];   //[ntracks]
   Float_t         track_charge[nmax_tr];   //[ntracks]
   Float_t         track_chi2[nmax_tr];   //[ntracks]
   Int_t           track_nvalidhits[nmax_tr];   //[ntracks]
   Int_t           track_npixhits[nmax_tr];   //[ntracks]
   Int_t           track_nonl1pixhits[nmax_tr];   //[ntracks]
   Int_t           track_missing[nmax_tr];   //[ntracks]
   Float_t         track_validfraction[nmax_tr];   //[ntracks]
   Float_t         track_validlast[nmax_tr];   //[ntracks]
   Int_t           track_qual[nmax_tr];   //[ntracks]
   Bool_t          track_qual2[nmax_tr];   //[ntracks]
   Float_t         track_dz[nmax_tr];   //[ntracks]
   Float_t         track_dxy[nmax_tr];   //[ntracks]
   Float_t         track_dz_0[nmax_tr];   //[ntracks]
   Float_t         track_dxy_0[nmax_tr];   //[ntracks]
   Float_t         track_dz_2[nmax_tr];   //[ntracks]
   Float_t         track_dxy_2[nmax_tr];   //[ntracks]
   Float_t         track_pvweight[nmax_tr];   //[ntracks]
   Float_t         track_pv0weight[nmax_tr];   //[ntracks]
   Int_t           track_index_hit[nmax_tr];   //[ntracks]
   Int_t           track_nhits[nmax_tr];   //[ntracks]
   Int_t           track_prescale[nmax_tr];   //[ntracks]
   Float_t         track_ih_ampl[nmax_tr];   //[ntracks]
   Float_t         track_ih_ampl_corr[nmax_tr];   //[ntracks]
   Float_t         track_ias_ampl[nmax_tr];   //[ntracks]
   Float_t         track_ias_ampl_corr[nmax_tr];   //[ntracks]
   Float_t         track_miniRelIso[nmax_tr];   //[ntracks]
   Float_t         track_TkRelIso[nmax_tr];   //[ntracks]
   Float_t         track_EoP[nmax_tr];   //[ntracks]
   Float_t         track_isoR005_sumChargedHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR005_sumNeutHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR005_sumPhotonPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR005_sumPUPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR01_sumChargedHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR01_sumNeutHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR01_sumPhotonPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR01_sumPUPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR03_sumChargedHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR03_sumNeutHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR03_sumPhotonPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR03_sumPUPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR05_sumChargedHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR05_sumNeutHadronPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR05_sumPhotonPt[nmax_tr];   //[ntracks]
   Float_t         track_isoR05_sumPUPt[nmax_tr];   //[ntracks]
   Float_t         track_probQ[nmax_tr];   //[ntracks]
   Float_t         track_probQNoL1[nmax_tr];   //[ntracks]
   Float_t         track_probXY[nmax_tr];   //[ntracks]
   Float_t         track_probXYNoL1[nmax_tr];   //[ntracks]
   Int_t           ndedxhits;
   UInt_t          dedx_detid[nmax_cl];   //[ndedxhits]
   Int_t           dedx_subdetid[nmax_cl];   //[ndedxhits]
   Int_t           dedx_modulgeom[nmax_cl];   //[ndedxhits]
   Int_t           dedx_layer[nmax_cl];   //[ndedxhits]
   Float_t         dedx_charge[nmax_cl];   //[ndedxhits]
   Float_t         dedx_pathlength[nmax_cl];   //[ndedxhits]
   Float_t         dedx_posx[nmax_cl];   //[ndedxhits]
   Float_t         dedx_posy[nmax_cl];   //[ndedxhits]
   Float_t         dedx_posz[nmax_cl];   //[ndedxhits]
   Bool_t          dedx_isstrip[nmax_cl];   //[ndedxhits]
   Bool_t          dedx_ispixel[nmax_cl];   //[ndedxhits]
   Bool_t          dedx_insideTkMod[nmax_cl];   //[ndedxhits]
   Float_t         dedx_probQ[nmax_cl];   //[ndedxhits]
   Float_t         dedx_probXY[nmax_cl];   //[ndedxhits]
   Float_t         dedx_probQNoL1[nmax_cl];   //[ndedxhits]
   Float_t         dedx_probXYNoL1[nmax_cl];   //[ndedxhits]
   Int_t           sclus_firstsclus[nmax_cl];   //[ndedxhits]
   Float_t         sclus_barycenter[nmax_cl];   //[ndedxhits]
   Float_t         sclus_charge[nmax_cl];   //[ndedxhits]
   Float_t         sclus_errorclus[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_ismerged[nmax_cl];   //[ndedxhits]
   Int_t           sclus_index_strip[nmax_cl];   //[ndedxhits]
   Int_t           sclus_nstrip[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_sat254[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_sat255[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_shape[nmax_cl];   //[ndedxhits]
   Int_t           sclus_index_strip_corr[nmax_cl];   //[ndedxhits]
   Int_t           sclus_nstrip_corr[nmax_cl];   //[ndedxhits]
   Float_t         sclus_charge_corr[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_clusclean[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_clusclean2[nmax_cl];   //[ndedxhits]
   Int_t           sclus_index_simhit[nmax_cl];   //[ndedxhits]
   Int_t           sclus_nsimhit[nmax_cl];   //[ndedxhits]
   Float_t         sclus_eloss[nmax_cl];   //[ndedxhits]
   Int_t           nstrips;
   Int_t           strip_ampl[nmax_st];   //[nstrips]
   Int_t           strip_channel[nmax_st];   //[nstrips]
   Int_t           nstrips_corr;
   Int_t           strip_ampl_corr[nmax_st];   //[nstrips_corr]
   Int_t           nsimhits;
   Int_t           simhit_pid[nmax_simhit];   //[nsimhits]
   Int_t           simhit_process[nmax_simhit];   //[nsimhits]
   Float_t         simhit_p[nmax_simhit];   //[nsimhits]
   Float_t         simhit_eloss[nmax_simhit];   //[nsimhits]
   Float_t         simhit_tof[nmax_simhit];   //[nsimhits]
   Float_t         simhit_segment[nmax_simhit];   //[nsimhits]
   Float_t         simhit_xentry[nmax_simhit];   //[nsimhits]
   Float_t         simhit_yentry[nmax_simhit];   //[nsimhits]
   Float_t         simhit_zentry[nmax_simhit];   //[nsimhits]
   Float_t         simhit_xexit[nmax_simhit];   //[nsimhits]
   Float_t         simhit_yexit[nmax_simhit];   //[nsimhits]
   Float_t         simhit_zexit[nmax_simhit];   //[nsimhits]
   Int_t           ndigi_sig_pix;
   UInt_t          sigdigi_pix_id[nmax_sigpix];   //[ndigi_sig_pix]
   Int_t           sigdigi_pix_adc[nmax_sigpix];   //[ndigi_sig_pix]
   Int_t           sigdigi_pix_channel[nmax_sigpix];   //[ndigi_sig_pix]
   Int_t           sigdigi_pix_row[nmax_sigpix];   //[ndigi_sig_pix]
   Int_t           sigdigi_pix_col[nmax_sigpix];   //[ndigi_sig_pix]
   Int_t           sigdigi_pix_layer[nmax_sigpix];   //[ndigi_sig_pix]
   Int_t           ndigi_sum_pix;
   UInt_t          sumdigi_pix_id[nmax_sumpix];   //[ndigi_sum_pix]
   Int_t           sumdigi_pix_adc[nmax_sumpix];   //[ndigi_sum_pix]
   Int_t           sumdigi_pix_channel[nmax_sumpix];   //[ndigi_sum_pix]
   Int_t           sumdigi_pix_row[nmax_sumpix];   //[ndigi_sum_pix]
   Int_t           sumdigi_pix_col[nmax_sumpix];   //[ndigi_sum_pix]
   Int_t           sumdigi_pix_layer[nmax_sumpix];   //[ndigi_sum_pix]
   Int_t           ndigi_sig_strip;
   UInt_t          sigdigi_strip_id[nmax_sigstrip];   //[ndigi_sig_strip]
   Int_t           sigdigi_strip_adc[nmax_sigstrip];   //[ndigi_sig_strip]
   Int_t           sigdigi_strip_channel[nmax_sigstrip];   //[ndigi_sig_strip]
   Int_t           sigdigi_strip_layer[nmax_sigstrip];   //[ndigi_sig_strip]
   Int_t           ndigi_sum_strip;
   UInt_t          sumdigi_strip_id[nmax_sumstrip];   //[ndigi_sum_strip]
   Int_t           sumdigi_strip_adc[nmax_sumstrip];   //[ndigi_sum_strip]
   Int_t           sumdigi_strip_channel[nmax_sumstrip];   //[ndigi_sum_strip]
   Int_t           sumdigi_strip_layer[nmax_sumstrip];   //[ndigi_sum_strip]
   Int_t           nmuons;
   Float_t         muon_pt[nmax_mu];   //[nmuons]
   Float_t         muon_ptSA[nmax_mu];   //[nmuons]
   Float_t         muon_ptIT[nmax_mu];   //[nmuons]
   Float_t         muon_p[nmax_mu];   //[nmuons]
   Float_t         muon_eta[nmax_mu];   //[nmuons]
   Float_t         muon_phi[nmax_mu];   //[nmuons]
   Bool_t          muon_isMatchesValid[nmax_mu];   //[nmuons]
   Bool_t          muon_isTrackerMuon[nmax_mu];   //[nmuons]
   Bool_t          muon_isGlobalMuon[nmax_mu];   //[nmuons]
   Bool_t          muon_isTightMuon[nmax_mu];   //[nmuons]
   Bool_t          muon_isMediumMuon[nmax_mu];   //[nmuons]
   Bool_t          muon_isLooseMuon[nmax_mu];   //[nmuons]
   Bool_t          muon_isHighPtMuon[nmax_mu];   //[nmuons]
   Float_t         muon_isoR04_sumChargedHadronPt[nmax_mu];   //[nmuons]
   Float_t         muon_isoR04_sumNeutHadronEt[nmax_mu];   //[nmuons]
   Float_t         muon_isoR04_sumPhotonEt[nmax_mu];   //[nmuons]
   Float_t         muon_isoR04_sumPUPt[nmax_mu];   //[nmuons]
   Float_t         muon_comb_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_comb_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_comb_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_comb_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_dt_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_dt_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_dt_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_dt_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_csc_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_csc_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_csc_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_csc_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_newcomb_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_newcomb_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_newcomb_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_newcomb_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_newdt_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_newdt_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_newdt_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_newdt_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_newcsc_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_newcsc_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_newcsc_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_newcsc_vertextime[nmax_mu];   //[nmuons]
   Int_t           nhscp;
   Int_t           hscp_type[nmax_hscp];   //[nhscp]
   Float_t         hscp_pt[nmax_hscp];   //[nhscp]
   Int_t           hscp_gen_id[nmax_hscp];   //[nhscp]
   Float_t         hscp_gen_dr[nmax_hscp];   //[nhscp]
   Int_t           hscp_track_idx[nmax_hscp];   //[nhscp]
   Int_t           hscp_muon_idx[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso0_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso0_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso0_hcal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso1_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso1_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso1_hcal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso2_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso2_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso2_hcal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso3_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso3_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso3_hcal[nmax_hscp];   //[nhscp]
   Float_t         pfmet;
   Float_t         calomet;
   Float_t         pfmht;
   Int_t           njets;
   Float_t         jet_pt[nmax_jet];   //[njets]
   Float_t         jet_eta[nmax_jet];   //[njets]
   Float_t         jet_phi[nmax_jet];   //[njets]
   Float_t         jet_E[nmax_jet];   //[njets]
   Float_t         jet_m[nmax_jet];   //[njets]
   Float_t         jet_id[nmax_jet];   //[njets]
   Float_t         jet_et[nmax_jet];   //[njets]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_event;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_ngoodpv;   //!
   TBranch        *b_goodIsFistPV;   //!
   TBranch        *b_pv_z;   //!
   TBranch        *b_pv_errz;   //!
   TBranch        *b_pv0_z;   //!
   TBranch        *b_pv0_errz;   //!
   TBranch        *b_hlt_mu45;   //!
   TBranch        *b_hlt_mu50;   //!
   TBranch        *b_hlt_tkmu100;   //!
   TBranch        *b_hlt_oldmu100;   //!
   TBranch        *b_hlt_pfmet_mht;   //!
   TBranch        *b_hlt_pfmet;   //!
   TBranch        *b_hlt_tkmu50;   //!
   TBranch        *b_InstLumi;   //!
   TBranch        *b_gen_pv_z;   //!
   TBranch        *b_ngenpart;   //!
   TBranch        *b_gen_pdg;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_mass;   //!
   TBranch        *b_gen_isHardProcess;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_moth_pdg;   //!
   TBranch        *b_gen_ndaughter;   //!
   TBranch        *b_gen_daughter_pdg;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_pterr;   //!
   TBranch        *b_track_p;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_chi2;   //!
   TBranch        *b_track_nvalidhits;   //!
   TBranch        *b_track_npixhits;   //!
   TBranch        *b_track_nonl1pixhits;   //!
   TBranch        *b_track_missing;   //!
   TBranch        *b_track_validfraction;   //!
   TBranch        *b_track_validlast;   //!
   TBranch        *b_track_qual;   //!
   TBranch        *b_track_qual2;   //!
   TBranch        *b_track_dz;   //!
   TBranch        *b_track_dxy;   //!
   TBranch        *b_track_dz_0;   //!
   TBranch        *b_track_dxy_0;   //!
   TBranch        *b_track_dz_2;   //!
   TBranch        *b_track_dxy_2;   //!
   TBranch        *b_track_pvweight;   //!
   TBranch        *b_track_pv0weight;   //!
   TBranch        *b_track_index_hit;   //!
   TBranch        *b_track_nhits;   //!
   TBranch        *b_track_prescale;   //!
   TBranch        *b_track_ih_ampl;   //!
   TBranch        *b_track_ih_ampl_corr;   //!
   TBranch        *b_track_ias_ampl;   //!
   TBranch        *b_track_ias_ampl_corr;   //!
   TBranch        *b_track_miniRelIso;   //!
   TBranch        *b_track_TkRelIso;   //!
   TBranch        *b_track_EoP;   //!
   TBranch        *b_track_isoR005_sumChargedHadronPt;   //!
   TBranch        *b_track_isoR005_sumNeutHadronPt;   //!
   TBranch        *b_track_isoR005_sumPhotonPt;   //!
   TBranch        *b_track_isoR005_sumPUPt;   //!
   TBranch        *b_track_isoR01_sumChargedHadronPt;   //!
   TBranch        *b_track_isoR01_sumNeutHadronPt;   //!
   TBranch        *b_track_isoR01_sumPhotonPt;   //!
   TBranch        *b_track_isoR01_sumPUPt;   //!
   TBranch        *b_track_isoR03_sumChargedHadronPt;   //!
   TBranch        *b_track_isoR03_sumNeutHadronPt;   //!
   TBranch        *b_track_isoR03_sumPhotonPt;   //!
   TBranch        *b_track_isoR03_sumPUPt;   //!
   TBranch        *b_track_isoR05_sumChargedHadronPt;   //!
   TBranch        *b_track_isoR05_sumNeutHadronPt;   //!
   TBranch        *b_track_isoR05_sumPhotonPt;   //!
   TBranch        *b_track_isoR05_sumPUPt;   //!
   TBranch        *b_track_probQ;   //!
   TBranch        *b_track_probQNoL1;   //!
   TBranch        *b_track_probXY;   //!
   TBranch        *b_track_probXYNoL1;   //!
   TBranch        *b_ndedxhits;   //!
   TBranch        *b_dedx_detid;   //!
   TBranch        *b_dedx_subdetid;   //!
   TBranch        *b_dedx_modulgeom;   //!
   TBranch        *b_dedx_layer;   //!
   TBranch        *b_dedx_charge;   //!
   TBranch        *b_dedx_pathlength;   //!
   TBranch        *b_dedx_posx;   //!
   TBranch        *b_dedx_posy;   //!
   TBranch        *b_dedx_posz;   //!
   TBranch        *b_dedx_isstrip;   //!
   TBranch        *b_dedx_ispixel;   //!
   TBranch        *b_dedx_insideTkMod;   //!
   TBranch        *b_dedx_probQ;   //!
   TBranch        *b_dedx_probXY;   //!
   TBranch        *b_dedx_probQNoL1;   //!
   TBranch        *b_dedx_probXYNoL1;   //!
   TBranch        *b_sclus_firstsclus;   //!
   TBranch        *b_sclus_barycenter;   //!
   TBranch        *b_sclus_charge;   //!
   TBranch        *b_sclus_errorclus;   //!
   TBranch        *b_sclus_ismerged;   //!
   TBranch        *b_sclus_index_strip;   //!
   TBranch        *b_sclus_nstrip;   //!
   TBranch        *b_sclus_sat254;   //!
   TBranch        *b_sclus_sat255;   //!
   TBranch        *b_sclus_shape;   //!
   TBranch        *b_sclus_index_strip_corr;   //!
   TBranch        *b_sclus_nstrip_corr;   //!
   TBranch        *b_sclus_charge_corr;   //!
   TBranch        *b_sclus_clusclean;   //!
   TBranch        *b_sclus_clusclean2;   //!
   TBranch        *b_sclus_index_simhit;   //!
   TBranch        *b_sclus_nsimhit;   //!
   TBranch        *b_sclus_eloss;   //!
   TBranch        *b_nstrips;   //!
   TBranch        *b_strip_ampl;   //!
   TBranch        *b_strip_channel;   //!
   TBranch        *b_nstrips_corr;   //!
   TBranch        *b_strip_ampl_corr;   //!
   TBranch        *b_nsimhits;   //!
   TBranch        *b_simhit_pid;   //!
   TBranch        *b_simhit_process;   //!
   TBranch        *b_simhit_p;   //!
   TBranch        *b_simhit_eloss;   //!
   TBranch        *b_simhit_tof;   //!
   TBranch        *b_simhit_segment;   //!
   TBranch        *b_simhit_xentry;   //!
   TBranch        *b_simhit_yentry;   //!
   TBranch        *b_simhit_zentry;   //!
   TBranch        *b_simhit_xexit;   //!
   TBranch        *b_simhit_yexit;   //!
   TBranch        *b_simhit_zexit;   //!
   TBranch        *b_ndigi_sig_pix;   //!
   TBranch        *b_sigdigi_pix_id;   //!
   TBranch        *b_sigdigi_pix_adc;   //!
   TBranch        *b_sigdigi_pix_channel;   //!
   TBranch        *b_sigdigi_pix_row;   //!
   TBranch        *b_sigdigi_pix_col;   //!
   TBranch        *b_sigdigi_pix_layer;   //!
   TBranch        *b_ndigi_sum_pix;   //!
   TBranch        *b_sumdigi_pix_id;   //!
   TBranch        *b_sumdigi_pix_adc;   //!
   TBranch        *b_sumdigi_pix_channel;   //!
   TBranch        *b_sumdigi_pix_row;   //!
   TBranch        *b_sumdigi_pix_col;   //!
   TBranch        *b_sumdigi_pix_layer;   //!
   TBranch        *b_ndigi_sig_strip;   //!
   TBranch        *b_sigdigi_strip_id;   //!
   TBranch        *b_sigdigi_strip_adc;   //!
   TBranch        *b_sigdigi_strip_channel;   //!
   TBranch        *b_sigdigi_strip_layer;   //!
   TBranch        *b_ndigi_sum_strip;   //!
   TBranch        *b_sumdigi_strip_id;   //!
   TBranch        *b_sumdigi_strip_adc;   //!
   TBranch        *b_sumdigi_strip_channel;   //!
   TBranch        *b_sumdigi_strip_layer;   //!
   TBranch        *b_nmuons;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_ptSA;   //!
   TBranch        *b_muon_ptIT;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_isMatchesValid;   //!
   TBranch        *b_muon_isTrackerMuon;   //!
   TBranch        *b_muon_isGlobalMuon;   //!
   TBranch        *b_muon_isTightMuon;   //!
   TBranch        *b_muon_isMediumMuon;   //!
   TBranch        *b_muon_isLooseMuon;   //!
   TBranch        *b_muon_isHighPtMuon;   //!
   TBranch        *b_muon_isoR04_sumChargedHadronPt;   //!
   TBranch        *b_muon_isoR04_sumNeutHadronEt;   //!
   TBranch        *b_muon_isoR04_sumPhotonEt;   //!
   TBranch        *b_muon_isoR04_sumPUPt;   //!
   TBranch        *b_muon_comb_inversebeta;   //!
   TBranch        *b_muon_comb_inversebetaerr;   //!
   TBranch        *b_muon_comb_tofndof;   //!
   TBranch        *b_muon_comb_vertextime;   //!
   TBranch        *b_muon_dt_inversebeta;   //!
   TBranch        *b_muon_dt_inversebetaerr;   //!
   TBranch        *b_muon_dt_tofndof;   //!
   TBranch        *b_muon_dt_vertextime;   //!
   TBranch        *b_muon_csc_inversebeta;   //!
   TBranch        *b_muon_csc_inversebetaerr;   //!
   TBranch        *b_muon_csc_tofndof;   //!
   TBranch        *b_muon_csc_vertextime;   //!
   TBranch        *b_muon_newcomb_inversebeta;   //!
   TBranch        *b_muon_newcomb_inversebetaerr;   //!
   TBranch        *b_muon_newcomb_tofndof;   //!
   TBranch        *b_muon_newcomb_vertextime;   //!
   TBranch        *b_muon_newdt_inversebeta;   //!
   TBranch        *b_muon_newdt_inversebetaerr;   //!
   TBranch        *b_muon_newdt_tofndof;   //!
   TBranch        *b_muon_newdt_vertextime;   //!
   TBranch        *b_muon_newcsc_inversebeta;   //!
   TBranch        *b_muon_newcsc_inversebetaerr;   //!
   TBranch        *b_muon_newcsc_tofndof;   //!
   TBranch        *b_muon_newcsc_vertextime;   //!
   TBranch        *b_nhscp;   //!
   TBranch        *b_hscp_type;   //!
   TBranch        *b_hscp_pt;   //!
   TBranch        *b_hscp_gen_id;   //!
   TBranch        *b_hscp_gen_dr;   //!
   TBranch        *b_hscp_track_idx;   //!
   TBranch        *b_hscp_muon_idx;   //!
   TBranch        *b_hscp_iso0_tk;   //!
   TBranch        *b_hscp_iso0_ecal;   //!
   TBranch        *b_hscp_iso0_hcal;   //!
   TBranch        *b_hscp_iso1_tk;   //!
   TBranch        *b_hscp_iso1_ecal;   //!
   TBranch        *b_hscp_iso1_hcal;   //!
   TBranch        *b_hscp_iso2_tk;   //!
   TBranch        *b_hscp_iso2_ecal;   //!
   TBranch        *b_hscp_iso2_hcal;   //!
   TBranch        *b_hscp_iso3_tk;   //!
   TBranch        *b_hscp_iso3_ecal;   //!
   TBranch        *b_hscp_iso3_hcal;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_calomet;   //!
   TBranch        *b_pfmht;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_E;   //!
   TBranch        *b_jet_m;   //!
   TBranch        *b_jet_id;   //!
   TBranch        *b_jet_et;   //!

   CreateTemplate(TTree *tree=0);
   virtual ~CreateTemplate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef CreateTemplate_cxx
CreateTemplate::CreateTemplate(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ROOT_Files/nt_mc_aod_gael_v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ROOT_Files/nt_mc_aod_gael_v2.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ROOT_Files/nt_mc_aod_gael_v2.root:/stage");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

CreateTemplate::~CreateTemplate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CreateTemplate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t CreateTemplate::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void CreateTemplate::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("ngoodpv", &ngoodpv, &b_ngoodpv);
   fChain->SetBranchAddress("goodIsFistPV", &goodIsFistPV, &b_goodIsFistPV);
   fChain->SetBranchAddress("pv_z", &pv_z, &b_pv_z);
   fChain->SetBranchAddress("pv_errz", &pv_errz, &b_pv_errz);
   fChain->SetBranchAddress("pv0_z", &pv0_z, &b_pv0_z);
   fChain->SetBranchAddress("pv0_errz", &pv0_errz, &b_pv0_errz);
   fChain->SetBranchAddress("hlt_mu45", &hlt_mu45, &b_hlt_mu45);
   fChain->SetBranchAddress("hlt_mu50", &hlt_mu50, &b_hlt_mu50);
   fChain->SetBranchAddress("hlt_tkmu100", &hlt_tkmu100, &b_hlt_tkmu100);
   fChain->SetBranchAddress("hlt_oldmu100", &hlt_oldmu100, &b_hlt_oldmu100);
   fChain->SetBranchAddress("hlt_pfmet_mht", &hlt_pfmet_mht, &b_hlt_pfmet_mht);
   fChain->SetBranchAddress("hlt_pfmet", &hlt_pfmet, &b_hlt_pfmet);
   fChain->SetBranchAddress("hlt_tkmu50", &hlt_tkmu50, &b_hlt_tkmu50);
   fChain->SetBranchAddress("InstLumi", &InstLumi, &b_InstLumi);
   fChain->SetBranchAddress("gen_pv_z", &gen_pv_z, &b_gen_pv_z);
   fChain->SetBranchAddress("ngenpart", &ngenpart, &b_ngenpart);
   fChain->SetBranchAddress("gen_pdg", gen_pdg, &b_gen_pdg);
   fChain->SetBranchAddress("gen_pt", gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_eta", gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_mass", gen_mass, &b_gen_mass);
   fChain->SetBranchAddress("gen_isHardProcess", gen_isHardProcess, &b_gen_isHardProcess);
   fChain->SetBranchAddress("gen_status", gen_status, &b_gen_status);
   fChain->SetBranchAddress("gen_moth_pdg", gen_moth_pdg, &b_gen_moth_pdg);
   fChain->SetBranchAddress("gen_ndaughter", gen_ndaughter, &b_gen_ndaughter);
   fChain->SetBranchAddress("gen_daughter_pdg", gen_daughter_pdg, &b_gen_daughter_pdg);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_pterr", track_pterr, &b_track_pterr);
   fChain->SetBranchAddress("track_p", track_p, &b_track_p);
   fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_charge", track_charge, &b_track_charge);
   fChain->SetBranchAddress("track_chi2", track_chi2, &b_track_chi2);
   fChain->SetBranchAddress("track_nvalidhits", track_nvalidhits, &b_track_nvalidhits);
   fChain->SetBranchAddress("track_npixhits", track_npixhits, &b_track_npixhits);
   fChain->SetBranchAddress("track_nonl1pixhits", track_nonl1pixhits, &b_track_nonl1pixhits);
   fChain->SetBranchAddress("track_missing", track_missing, &b_track_missing);
   fChain->SetBranchAddress("track_validfraction", track_validfraction, &b_track_validfraction);
   fChain->SetBranchAddress("track_validlast", track_validlast, &b_track_validlast);
   fChain->SetBranchAddress("track_qual", track_qual, &b_track_qual);
   fChain->SetBranchAddress("track_qual2", track_qual2, &b_track_qual2);
   fChain->SetBranchAddress("track_dz", track_dz, &b_track_dz);
   fChain->SetBranchAddress("track_dxy", track_dxy, &b_track_dxy);
   fChain->SetBranchAddress("track_dz_0", track_dz_0, &b_track_dz_0);
   fChain->SetBranchAddress("track_dxy_0", track_dxy_0, &b_track_dxy_0);
   fChain->SetBranchAddress("track_dz_2", track_dz_2, &b_track_dz_2);
   fChain->SetBranchAddress("track_dxy_2", track_dxy_2, &b_track_dxy_2);
   fChain->SetBranchAddress("track_pvweight", track_pvweight, &b_track_pvweight);
   fChain->SetBranchAddress("track_pv0weight", track_pv0weight, &b_track_pv0weight);
   fChain->SetBranchAddress("track_index_hit", track_index_hit, &b_track_index_hit);
   fChain->SetBranchAddress("track_nhits", track_nhits, &b_track_nhits);
   fChain->SetBranchAddress("track_prescale", track_prescale, &b_track_prescale);
   fChain->SetBranchAddress("track_ih_ampl", track_ih_ampl, &b_track_ih_ampl);
   fChain->SetBranchAddress("track_ih_ampl_corr", track_ih_ampl_corr, &b_track_ih_ampl_corr);
   fChain->SetBranchAddress("track_ias_ampl", track_ias_ampl, &b_track_ias_ampl);
   fChain->SetBranchAddress("track_ias_ampl_corr", track_ias_ampl_corr, &b_track_ias_ampl_corr);
   fChain->SetBranchAddress("track_miniRelIso", track_miniRelIso, &b_track_miniRelIso);
   fChain->SetBranchAddress("track_TkRelIso", track_TkRelIso, &b_track_TkRelIso);
   fChain->SetBranchAddress("track_EoP", track_EoP, &b_track_EoP);
   fChain->SetBranchAddress("track_isoR005_sumChargedHadronPt", track_isoR005_sumChargedHadronPt, &b_track_isoR005_sumChargedHadronPt);
   fChain->SetBranchAddress("track_isoR005_sumNeutHadronPt", track_isoR005_sumNeutHadronPt, &b_track_isoR005_sumNeutHadronPt);
   fChain->SetBranchAddress("track_isoR005_sumPhotonPt", track_isoR005_sumPhotonPt, &b_track_isoR005_sumPhotonPt);
   fChain->SetBranchAddress("track_isoR005_sumPUPt", track_isoR005_sumPUPt, &b_track_isoR005_sumPUPt);
   fChain->SetBranchAddress("track_isoR01_sumChargedHadronPt", track_isoR01_sumChargedHadronPt, &b_track_isoR01_sumChargedHadronPt);
   fChain->SetBranchAddress("track_isoR01_sumNeutHadronPt", track_isoR01_sumNeutHadronPt, &b_track_isoR01_sumNeutHadronPt);
   fChain->SetBranchAddress("track_isoR01_sumPhotonPt", track_isoR01_sumPhotonPt, &b_track_isoR01_sumPhotonPt);
   fChain->SetBranchAddress("track_isoR01_sumPUPt", track_isoR01_sumPUPt, &b_track_isoR01_sumPUPt);
   fChain->SetBranchAddress("track_isoR03_sumChargedHadronPt", track_isoR03_sumChargedHadronPt, &b_track_isoR03_sumChargedHadronPt);
   fChain->SetBranchAddress("track_isoR03_sumNeutHadronPt", track_isoR03_sumNeutHadronPt, &b_track_isoR03_sumNeutHadronPt);
   fChain->SetBranchAddress("track_isoR03_sumPhotonPt", track_isoR03_sumPhotonPt, &b_track_isoR03_sumPhotonPt);
   fChain->SetBranchAddress("track_isoR03_sumPUPt", track_isoR03_sumPUPt, &b_track_isoR03_sumPUPt);
   fChain->SetBranchAddress("track_isoR05_sumChargedHadronPt", track_isoR05_sumChargedHadronPt, &b_track_isoR05_sumChargedHadronPt);
   fChain->SetBranchAddress("track_isoR05_sumNeutHadronPt", track_isoR05_sumNeutHadronPt, &b_track_isoR05_sumNeutHadronPt);
   fChain->SetBranchAddress("track_isoR05_sumPhotonPt", track_isoR05_sumPhotonPt, &b_track_isoR05_sumPhotonPt);
   fChain->SetBranchAddress("track_isoR05_sumPUPt", track_isoR05_sumPUPt, &b_track_isoR05_sumPUPt);
   fChain->SetBranchAddress("track_probQ", track_probQ, &b_track_probQ);
   fChain->SetBranchAddress("track_probQNoL1", track_probQNoL1, &b_track_probQNoL1);
   fChain->SetBranchAddress("track_probXY", track_probXY, &b_track_probXY);
   fChain->SetBranchAddress("track_probXYNoL1", track_probXYNoL1, &b_track_probXYNoL1);
   fChain->SetBranchAddress("ndedxhits", &ndedxhits, &b_ndedxhits);
   fChain->SetBranchAddress("dedx_detid", dedx_detid, &b_dedx_detid);
   fChain->SetBranchAddress("dedx_subdetid", dedx_subdetid, &b_dedx_subdetid);
   fChain->SetBranchAddress("dedx_modulgeom", dedx_modulgeom, &b_dedx_modulgeom);
   fChain->SetBranchAddress("dedx_layer", dedx_layer, &b_dedx_layer);
   fChain->SetBranchAddress("dedx_charge", dedx_charge, &b_dedx_charge);
   fChain->SetBranchAddress("dedx_pathlength", dedx_pathlength, &b_dedx_pathlength);
   fChain->SetBranchAddress("dedx_posx", dedx_posx, &b_dedx_posx);
   fChain->SetBranchAddress("dedx_posy", dedx_posy, &b_dedx_posy);
   fChain->SetBranchAddress("dedx_posz", dedx_posz, &b_dedx_posz);
   fChain->SetBranchAddress("dedx_isstrip", dedx_isstrip, &b_dedx_isstrip);
   fChain->SetBranchAddress("dedx_ispixel", dedx_ispixel, &b_dedx_ispixel);
   fChain->SetBranchAddress("dedx_insideTkMod", dedx_insideTkMod, &b_dedx_insideTkMod);
   fChain->SetBranchAddress("dedx_probQ", dedx_probQ, &b_dedx_probQ);
   fChain->SetBranchAddress("dedx_probXY", dedx_probXY, &b_dedx_probXY);
   fChain->SetBranchAddress("dedx_probQNoL1", dedx_probQNoL1, &b_dedx_probQNoL1);
   fChain->SetBranchAddress("dedx_probXYNoL1", dedx_probXYNoL1, &b_dedx_probXYNoL1);
   fChain->SetBranchAddress("sclus_firstsclus", sclus_firstsclus, &b_sclus_firstsclus);
   fChain->SetBranchAddress("sclus_barycenter", sclus_barycenter, &b_sclus_barycenter);
   fChain->SetBranchAddress("sclus_charge", sclus_charge, &b_sclus_charge);
   fChain->SetBranchAddress("sclus_errorclus", sclus_errorclus, &b_sclus_errorclus);
   fChain->SetBranchAddress("sclus_ismerged", sclus_ismerged, &b_sclus_ismerged);
   fChain->SetBranchAddress("sclus_index_strip", sclus_index_strip, &b_sclus_index_strip);
   fChain->SetBranchAddress("sclus_nstrip", sclus_nstrip, &b_sclus_nstrip);
   fChain->SetBranchAddress("sclus_sat254", sclus_sat254, &b_sclus_sat254);
   fChain->SetBranchAddress("sclus_sat255", sclus_sat255, &b_sclus_sat255);
   fChain->SetBranchAddress("sclus_shape", sclus_shape, &b_sclus_shape);
   fChain->SetBranchAddress("sclus_index_strip_corr", sclus_index_strip_corr, &b_sclus_index_strip_corr);
   fChain->SetBranchAddress("sclus_nstrip_corr", sclus_nstrip_corr, &b_sclus_nstrip_corr);
   fChain->SetBranchAddress("sclus_charge_corr", sclus_charge_corr, &b_sclus_charge_corr);
   fChain->SetBranchAddress("sclus_clusclean", sclus_clusclean, &b_sclus_clusclean);
   fChain->SetBranchAddress("sclus_clusclean2", sclus_clusclean2, &b_sclus_clusclean2);
   fChain->SetBranchAddress("sclus_index_simhit", sclus_index_simhit, &b_sclus_index_simhit);
   fChain->SetBranchAddress("sclus_nsimhit", sclus_nsimhit, &b_sclus_nsimhit);
   fChain->SetBranchAddress("sclus_eloss", sclus_eloss, &b_sclus_eloss);
   fChain->SetBranchAddress("nstrips", &nstrips, &b_nstrips);
   fChain->SetBranchAddress("strip_ampl", strip_ampl, &b_strip_ampl);
   fChain->SetBranchAddress("strip_channel", strip_channel, &b_strip_channel);
   fChain->SetBranchAddress("nstrips_corr", &nstrips_corr, &b_nstrips_corr);
   fChain->SetBranchAddress("strip_ampl_corr", strip_ampl_corr, &b_strip_ampl_corr);
   fChain->SetBranchAddress("nsimhits", &nsimhits, &b_nsimhits);
   fChain->SetBranchAddress("simhit_pid", simhit_pid, &b_simhit_pid);
   fChain->SetBranchAddress("simhit_process", simhit_process, &b_simhit_process);
   fChain->SetBranchAddress("simhit_p", simhit_p, &b_simhit_p);
   fChain->SetBranchAddress("simhit_eloss", simhit_eloss, &b_simhit_eloss);
   fChain->SetBranchAddress("simhit_tof", simhit_tof, &b_simhit_tof);
   fChain->SetBranchAddress("simhit_segment", simhit_segment, &b_simhit_segment);
   fChain->SetBranchAddress("simhit_xentry", simhit_xentry, &b_simhit_xentry);
   fChain->SetBranchAddress("simhit_yentry", simhit_yentry, &b_simhit_yentry);
   fChain->SetBranchAddress("simhit_zentry", simhit_zentry, &b_simhit_zentry);
   fChain->SetBranchAddress("simhit_xexit", simhit_xexit, &b_simhit_xexit);
   fChain->SetBranchAddress("simhit_yexit", simhit_yexit, &b_simhit_yexit);
   fChain->SetBranchAddress("simhit_zexit", simhit_zexit, &b_simhit_zexit);
   fChain->SetBranchAddress("ndigi_sig_pix", &ndigi_sig_pix, &b_ndigi_sig_pix);
   fChain->SetBranchAddress("sigdigi_pix_id", sigdigi_pix_id, &b_sigdigi_pix_id);
   fChain->SetBranchAddress("sigdigi_pix_adc", sigdigi_pix_adc, &b_sigdigi_pix_adc);
   fChain->SetBranchAddress("sigdigi_pix_channel", sigdigi_pix_channel, &b_sigdigi_pix_channel);
   fChain->SetBranchAddress("sigdigi_pix_row", sigdigi_pix_row, &b_sigdigi_pix_row);
   fChain->SetBranchAddress("sigdigi_pix_col", sigdigi_pix_col, &b_sigdigi_pix_col);
   fChain->SetBranchAddress("sigdigi_pix_layer", sigdigi_pix_layer, &b_sigdigi_pix_layer);
   fChain->SetBranchAddress("ndigi_sum_pix", &ndigi_sum_pix, &b_ndigi_sum_pix);
   fChain->SetBranchAddress("sumdigi_pix_id", sumdigi_pix_id, &b_sumdigi_pix_id);
   fChain->SetBranchAddress("sumdigi_pix_adc", sumdigi_pix_adc, &b_sumdigi_pix_adc);
   fChain->SetBranchAddress("sumdigi_pix_channel", sumdigi_pix_channel, &b_sumdigi_pix_channel);
   fChain->SetBranchAddress("sumdigi_pix_row", sumdigi_pix_row, &b_sumdigi_pix_row);
   fChain->SetBranchAddress("sumdigi_pix_col", sumdigi_pix_col, &b_sumdigi_pix_col);
   fChain->SetBranchAddress("sumdigi_pix_layer", sumdigi_pix_layer, &b_sumdigi_pix_layer);
   fChain->SetBranchAddress("ndigi_sig_strip", &ndigi_sig_strip, &b_ndigi_sig_strip);
   fChain->SetBranchAddress("sigdigi_strip_id", sigdigi_strip_id, &b_sigdigi_strip_id);
   fChain->SetBranchAddress("sigdigi_strip_adc", sigdigi_strip_adc, &b_sigdigi_strip_adc);
   fChain->SetBranchAddress("sigdigi_strip_channel", sigdigi_strip_channel, &b_sigdigi_strip_channel);
   fChain->SetBranchAddress("sigdigi_strip_layer", sigdigi_strip_layer, &b_sigdigi_strip_layer);
   fChain->SetBranchAddress("ndigi_sum_strip", &ndigi_sum_strip, &b_ndigi_sum_strip);
   fChain->SetBranchAddress("sumdigi_strip_id", sumdigi_strip_id, &b_sumdigi_strip_id);
   fChain->SetBranchAddress("sumdigi_strip_adc", sumdigi_strip_adc, &b_sumdigi_strip_adc);
   fChain->SetBranchAddress("sumdigi_strip_channel", sumdigi_strip_channel, &b_sumdigi_strip_channel);
   fChain->SetBranchAddress("sumdigi_strip_layer", sumdigi_strip_layer, &b_sumdigi_strip_layer);
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_ptSA", muon_ptSA, &b_muon_ptSA);
   fChain->SetBranchAddress("muon_ptIT", muon_ptIT, &b_muon_ptIT);
   fChain->SetBranchAddress("muon_p", muon_p, &b_muon_p);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_isMatchesValid", muon_isMatchesValid, &b_muon_isMatchesValid);
   fChain->SetBranchAddress("muon_isTrackerMuon", muon_isTrackerMuon, &b_muon_isTrackerMuon);
   fChain->SetBranchAddress("muon_isGlobalMuon", muon_isGlobalMuon, &b_muon_isGlobalMuon);
   fChain->SetBranchAddress("muon_isTightMuon", muon_isTightMuon, &b_muon_isTightMuon);
   fChain->SetBranchAddress("muon_isMediumMuon", muon_isMediumMuon, &b_muon_isMediumMuon);
   fChain->SetBranchAddress("muon_isLooseMuon", muon_isLooseMuon, &b_muon_isLooseMuon);
   fChain->SetBranchAddress("muon_isHighPtMuon", muon_isHighPtMuon, &b_muon_isHighPtMuon);
   fChain->SetBranchAddress("muon_isoR04_sumChargedHadronPt", muon_isoR04_sumChargedHadronPt, &b_muon_isoR04_sumChargedHadronPt);
   fChain->SetBranchAddress("muon_isoR04_sumNeutHadronEt", muon_isoR04_sumNeutHadronEt, &b_muon_isoR04_sumNeutHadronEt);
   fChain->SetBranchAddress("muon_isoR04_sumPhotonEt", muon_isoR04_sumPhotonEt, &b_muon_isoR04_sumPhotonEt);
   fChain->SetBranchAddress("muon_isoR04_sumPUPt", muon_isoR04_sumPUPt, &b_muon_isoR04_sumPUPt);
   fChain->SetBranchAddress("muon_comb_inversebeta", muon_comb_inversebeta, &b_muon_comb_inversebeta);
   fChain->SetBranchAddress("muon_comb_inversebetaerr", muon_comb_inversebetaerr, &b_muon_comb_inversebetaerr);
   fChain->SetBranchAddress("muon_comb_tofndof", muon_comb_tofndof, &b_muon_comb_tofndof);
   fChain->SetBranchAddress("muon_comb_vertextime", muon_comb_vertextime, &b_muon_comb_vertextime);
   fChain->SetBranchAddress("muon_dt_inversebeta", muon_dt_inversebeta, &b_muon_dt_inversebeta);
   fChain->SetBranchAddress("muon_dt_inversebetaerr", muon_dt_inversebetaerr, &b_muon_dt_inversebetaerr);
   fChain->SetBranchAddress("muon_dt_tofndof", muon_dt_tofndof, &b_muon_dt_tofndof);
   fChain->SetBranchAddress("muon_dt_vertextime", muon_dt_vertextime, &b_muon_dt_vertextime);
   fChain->SetBranchAddress("muon_csc_inversebeta", muon_csc_inversebeta, &b_muon_csc_inversebeta);
   fChain->SetBranchAddress("muon_csc_inversebetaerr", muon_csc_inversebetaerr, &b_muon_csc_inversebetaerr);
   fChain->SetBranchAddress("muon_csc_tofndof", muon_csc_tofndof, &b_muon_csc_tofndof);
   fChain->SetBranchAddress("muon_csc_vertextime", muon_csc_vertextime, &b_muon_csc_vertextime);
   fChain->SetBranchAddress("muon_newcomb_inversebeta", muon_newcomb_inversebeta, &b_muon_newcomb_inversebeta);
   fChain->SetBranchAddress("muon_newcomb_inversebetaerr", muon_newcomb_inversebetaerr, &b_muon_newcomb_inversebetaerr);
   fChain->SetBranchAddress("muon_newcomb_tofndof", muon_newcomb_tofndof, &b_muon_newcomb_tofndof);
   fChain->SetBranchAddress("muon_newcomb_vertextime", muon_newcomb_vertextime, &b_muon_newcomb_vertextime);
   fChain->SetBranchAddress("muon_newdt_inversebeta", muon_newdt_inversebeta, &b_muon_newdt_inversebeta);
   fChain->SetBranchAddress("muon_newdt_inversebetaerr", muon_newdt_inversebetaerr, &b_muon_newdt_inversebetaerr);
   fChain->SetBranchAddress("muon_newdt_tofndof", muon_newdt_tofndof, &b_muon_newdt_tofndof);
   fChain->SetBranchAddress("muon_newdt_vertextime", muon_newdt_vertextime, &b_muon_newdt_vertextime);
   fChain->SetBranchAddress("muon_newcsc_inversebeta", muon_newcsc_inversebeta, &b_muon_newcsc_inversebeta);
   fChain->SetBranchAddress("muon_newcsc_inversebetaerr", muon_newcsc_inversebetaerr, &b_muon_newcsc_inversebetaerr);
   fChain->SetBranchAddress("muon_newcsc_tofndof", muon_newcsc_tofndof, &b_muon_newcsc_tofndof);
   fChain->SetBranchAddress("muon_newcsc_vertextime", muon_newcsc_vertextime, &b_muon_newcsc_vertextime);
   fChain->SetBranchAddress("nhscp", &nhscp, &b_nhscp);
   fChain->SetBranchAddress("hscp_type", hscp_type, &b_hscp_type);
   fChain->SetBranchAddress("hscp_pt", hscp_pt, &b_hscp_pt);
   fChain->SetBranchAddress("hscp_gen_id", hscp_gen_id, &b_hscp_gen_id);
   fChain->SetBranchAddress("hscp_gen_dr", hscp_gen_dr, &b_hscp_gen_dr);
   fChain->SetBranchAddress("hscp_track_idx", hscp_track_idx, &b_hscp_track_idx);
   fChain->SetBranchAddress("hscp_muon_idx", hscp_muon_idx, &b_hscp_muon_idx);
   fChain->SetBranchAddress("hscp_iso0_tk", hscp_iso0_tk, &b_hscp_iso0_tk);
   fChain->SetBranchAddress("hscp_iso0_ecal", hscp_iso0_ecal, &b_hscp_iso0_ecal);
   fChain->SetBranchAddress("hscp_iso0_hcal", hscp_iso0_hcal, &b_hscp_iso0_hcal);
   fChain->SetBranchAddress("hscp_iso1_tk", hscp_iso1_tk, &b_hscp_iso1_tk);
   fChain->SetBranchAddress("hscp_iso1_ecal", hscp_iso1_ecal, &b_hscp_iso1_ecal);
   fChain->SetBranchAddress("hscp_iso1_hcal", hscp_iso1_hcal, &b_hscp_iso1_hcal);
   fChain->SetBranchAddress("hscp_iso2_tk", hscp_iso2_tk, &b_hscp_iso2_tk);
   fChain->SetBranchAddress("hscp_iso2_ecal", hscp_iso2_ecal, &b_hscp_iso2_ecal);
   fChain->SetBranchAddress("hscp_iso2_hcal", hscp_iso2_hcal, &b_hscp_iso2_hcal);
   fChain->SetBranchAddress("hscp_iso3_tk", hscp_iso3_tk, &b_hscp_iso3_tk);
   fChain->SetBranchAddress("hscp_iso3_ecal", hscp_iso3_ecal, &b_hscp_iso3_ecal);
   fChain->SetBranchAddress("hscp_iso3_hcal", hscp_iso3_hcal, &b_hscp_iso3_hcal);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("calomet", &calomet, &b_calomet);
   fChain->SetBranchAddress("pfmht", &pfmht, &b_pfmht);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_E", jet_E, &b_jet_E);
   fChain->SetBranchAddress("jet_m", jet_m, &b_jet_m);
   fChain->SetBranchAddress("jet_id", jet_id, &b_jet_id);
   fChain->SetBranchAddress("jet_et", jet_et, &b_jet_et);
   Notify();
}

Bool_t CreateTemplate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void CreateTemplate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t CreateTemplate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif // #ifdef CreateTemplate_cxx
