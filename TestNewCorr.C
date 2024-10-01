#define TestNewCorr_cxx
#include "TestNewCorr.h"

#include "CorrFunctions.h"
#include <vector>
#include <stdio.h>
#include <TMath.h>

using namespace std;

void TestNewCorr::Loop()
{
    if (fChain == 0) return -1;
    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    TFile* SaveData = new TFile("TestOnSigdigi.root", "RECREATE");

    TH1F* QNewCorr_VS_Qtrue = new TH1F("QNewCorr_VS_Qtrue","QNewCorr_VS_Qtrue",150,0,2);
    QNewCorr_VS_Qtrue->GetXaxis()->SetTitle("Q^{new corr} / Q^{true}");

    cout << endl;
    cout << "Histograms saved in: " << SaveData->GetName() << endl;
    cout << endl;
    
    //nentries=10000;
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

            // TRACK LOOP
        for (int itr=0; itr<ntracks; itr++)
        {
            bool selection = true;

            if (abs(track_eta[itr]) >= 2.4) selection=false;
            if (track_npixhits[itr] < 2) selection=false;
            if (track_validfraction[itr] <= 0.8) selection=false;
            if (track_nvalidhits[itr] < 10) selection=false;
            if (!track_qual[itr]) selection=false;
            if (track_chi2[itr] >= 5) selection=false;
            if (abs(track_dz[itr]) >= 0.1) selection=false;
            if (abs(track_dxy[itr]) >= 0.02) selection=false;

            if (selection)
            {
                for (int idedx=track_index_hit[itr]; idedx<track_index_hit[itr]+track_nhits[itr]; idedx++)
                {
                    if (dedx_isstrip[idedx])
                    {
                        std::vector <int> Qsat;
                        std::vector <int> Qtrue;
                        int cpt_sat = 0;
                    
                        for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
                        {
                            if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
                            {
                                Qtrue.push_back(sigdigi_strip_adc[isigdigi]);

                                if (sigdigi_strip_adc[isigdigi]<=253) Qsat.push_back(sigdigi_strip_adc[isigdigi]);
                                else if (sigdigi_strip_adc[isigdigi]>253 && sigdigi_strip_adc[isigdigi]<=1023) Qsat.push_back(254);
                                else if (sigdigi_strip_adc[isigdigi]>1023) Qsat.push_back(255);

                                if (sigdigi_strip_adc[isigdigi]>253) cpt_sat++;
                            }
                        }

                        if (!Qsat.empty())
                        {     
                            int sumTrue = accumulate(Qtrue.begin(), Qtrue.end(), 0);

                            int Qcorr = ReturnCorr(Qsat, dedx_subdetid[idedx], dedx_detid[idedx]);

                            // Check if the correction is working (only if saturation) for:
                            //      - One single saturated strips
                            //      - Two consecutive saturated strips
                            if (cpt_sat > 0) QNewCorr_VS_Qtrue->Fill(1.*Qcorr/sumTrue);
                        }

                    } // end is_strip
                } // end dedx loop
            } // end selection on track
        } // end track loop

        if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100.0*jentry/nentries<<" %)"<<endl;

    } // end entry loop

    SaveData->cd();
    QNewCorr_VS_Qtrue->Write();
    SaveData->Close();
}