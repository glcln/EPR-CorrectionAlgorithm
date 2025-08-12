#define TestNewCorr_cxx
#include "TestNewCorr.h"

#include "CorrFunctions.h"
#include <vector>
#include <stdio.h>
#include <TMath.h>

using namespace std;

std::vector<int> SaturationCorrection(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat)
{
  std::vector<int> QII;

  // if clusters too small or too large --> no correction
  if(Q.size()<2 || Q.size()>8)
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    return QII;
  }
  
  if(way)
  {
    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
    
    if(*mQ>253)
    {
      // if more than one saturated strip --> no correction
      if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253) return Q;
      
      // if one saturated strip and adjacent strips not saturated and above a given threshold --> correction
      if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1))<40)
      {
        QII.push_back( (10*(*(mQ-1)) + 10*(*(mQ+1))) /2 );
        return QII;
      }
    }
    else return Q; // no saturation --> no correction
  }

  // do nothing else
  return Q;
}

std::vector<int> CrossTalkInv(const std::vector<int>& Q,
                              const float x1 = 0.10,
                              const float x2 = 0.04,
                              bool way = true,
                              float threshold = 20,
                              float thresholdSat = 25,
                              bool isClusterCleaning = false) {
  const unsigned N = Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N, 0);
  Double_t a = 1 - 2 * x1 - 2 * x2;
  TMatrix A(N, N);

  // EXCLUSION PART: ANOMALOUS SHAPES or SATURATION
  // if clusters too small or too large --> no correction
  if(Q.size()<2 || Q.size()>8)
  {
    for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
    return QII;
  }
  
  if(way)
  {
    vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
    
    if(*mQ>253)
    {
      // if more than one saturated strip --> no correction
      if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253) return Q;
      
      // if one saturated strip and adjacent strips not saturated and above a given threshold --> correction
      if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1))<40)
      {
        QII.push_back( (10*(*(mQ-1)) + 10*(*(mQ+1))) /2 );
        return QII;
      }
    }
  }
  //---

  for (unsigned int i = 0; i < N; i++) {
    A(i, i) = a;
    if (i < N - 1) {
      A(i + 1, i) = x1;
      A(i, i + 1) = x1;
    } else
      continue;
    if (i < N - 2) {
      A(i + 2, i) = x2;
      A(i, i + 2) = x2;
    }
  }

  if (N == 1)
    A(0, 0) = 1 / a;
  else
    A.InvertFast();

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = 0; j < N; j++) {
      QI[i] += A(i, j) * (float)Q[j];
    }
  }

  for (unsigned int i = 0; i < QI.size(); i++) {
    if (QI[i] < threshold)
      QI[i] = 0;
    QII.push_back((int)QI[i]);
  }

  return QII;
}

std::vector <int> CrossTalkInvInStrip(const std::vector<int>& Q,
                                      const int layer,
                                      bool IfCorrApplied = true,
                                      float threshold = 20, 
                                      const float x1 = 0.10,
                                      const float x2 = 0.04) {
 
        // EXCLUSION PART: ANOMALOUS SHAPES or SATURATION
    if (Q.empty()) return {};
    std::vector<int> QII;
    if(Q.size()<2 || Q.size()>8)
    {
        for (unsigned int i=0;i<Q.size();i++) QII.push_back((int) Q[i]);
        return QII;
    }

    if(IfCorrApplied)
    {
        int N_sat = 0;
        for (unsigned int i=0; i<Q.size(); i++) { if (Q[i]>=254) N_sat++;}

        if(N_sat>0)
        {
            bool AreSameCluster = true;
            QII = ReturnCorrVec(Q, layer, AreSameCluster);
            if (!AreSameCluster) return QII;    // saturation correction
        }
    }


      // CROSS-TALK INVERSION
    QII.clear();
    double a = 1 - 2*x1 - 2*x2;
    const unsigned N = Q.size();
    std::vector<float> QI(N, 0);
    TMatrix A(N, N);

    for (unsigned int i = 0; i < N; i++) {
        A(i, i) = a;
        if (i < N - 1) {
        A(i + 1, i) = x1;
        A(i, i + 1) = x1;
        }
        else continue;
        if (i < N - 2) {
        A(i + 2, i) = x2;
        A(i, i + 2) = x2;
        }
    }
  
    //    A =   a  x1 x2 0-------------0
    //          x1  a x1 x2 0          |
    //          x2 x1  a x1 x2 0       |
    //          0 x2 x1  a x1 x2 0     |
    //          |  0               0   |
    //          |     0              0 |
    //          |        0             0
    //          0----------0 x2 x1  a x1
    

    if (N == 1) A(0, 0) = 1/a;
    else A.InvertFast();

    for (unsigned int i = 0; i < N; i++) {
        for (unsigned int j = 0; j < N; j++) {
        QI[i] += A(i, j) * (float)Q[j];
        }
    }

    for (unsigned int i = 0; i < QI.size(); i++) {
        if (QI[i] < threshold) QI[i] = 0;
        QII.push_back((int)QI[i]);
    }

    return QII;
}

bool clusterCleaning(std::vector<int> ampls, int crosstalkInv = 0, uint8_t* exitCode = nullptr) {


  // ---------------- Count the number of maximas    --------------------------
  //----------------------------------------------------------------------------
  Int_t NofMax = 0;
  Int_t recur255 = 1;
  Int_t recur254 = 1;
  bool MaxOnStart = false;
  bool MaxInMiddle = false, MaxOnEnd = false;
  Int_t MaxPos = 0;
  // Start with a max
  if (ampls.size() != 1 &&
      ((ampls[0] > ampls[1]) ||
       (ampls.size() > 2 && ampls[0] == ampls[1] && ampls[1] > ampls[2] && ampls[0] != 254 && ampls[0] != 255) ||
       (ampls.size() == 2 && ampls[0] == ampls[1] && ampls[0] != 254 && ampls[0] != 255))) {
    NofMax = NofMax + 1;
    MaxOnStart = true;
  }

  // Max in the middle (strip between other ones)
  if (ampls.size() > 2) {
    for (unsigned int i = 1; i < ampls.size() - 1; i++) {
      if ((ampls[i] > ampls[i - 1] && ampls[i] > ampls[i + 1]) ||
          (ampls.size() > 3 && i > 0 && i < ampls.size() - 2 && ampls[i] == ampls[i + 1] && ampls[i] > ampls[i - 1] &&
           ampls[i] > ampls[i + 2] && ampls[i] != 254 && ampls[i] != 255)) {
        NofMax = NofMax + 1;
        MaxInMiddle = true;
        MaxPos = i;
      }
      if (ampls[i] == 255 && ampls[i] == ampls[i - 1]) {
        recur255 = recur255 + 1;
        MaxPos = i - (recur255 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
      if (ampls[i] == 254 && ampls[i] == ampls[i - 1]) {
        recur254 = recur254 + 1;
        MaxPos = i - (recur254 / 2);
        if (ampls[i] > ampls[i + 1]) {
          NofMax = NofMax + 1;
          MaxInMiddle = true;
        }
      }
    }
  }
  // Max at the end of the cluster
  if (ampls.size() > 1) {
    if (ampls[ampls.size() - 1] > ampls[ampls.size() - 2] ||
        (ampls.size() > 2 && ampls[ampls.size() - 1] == ampls[ampls.size() - 2] &&
         ampls[ampls.size() - 2] > ampls[ampls.size() - 3]) ||
        ampls[ampls.size() - 1] == 255) {
      NofMax = NofMax + 1;
      MaxOnEnd = true;
    }
  }
  // If only one strip is hit
  if (ampls.size() == 1) {
    NofMax = 1;
  }

  // --- SHAPE SELECTION
  //------------------------------------------------------------------------
  //
  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //
  //   bool shapetest=true;
  bool shapecdtn = false;
  if (exitCode)
    *exitCode = 255;

  if (crosstalkInv == 1) {
    if (NofMax == 1) {
      shapecdtn = true;
      if (exitCode)
        *exitCode = 0;
    }
    return shapecdtn;
  }

  //      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
  Float_t C_M = 0.0;
  Float_t C_D = 0.0;
  Float_t C_Mn = 10000;
  Float_t C_Dn = 10000;
  Float_t C_Mnn = 10000;
  Float_t C_Dnn = 10000;
  Int_t CDPos;
  Float_t coeff1 = 1.7;
  Float_t coeff2 = 2.0;
  Float_t coeffn = 0.10;
  Float_t coeffnn = 0.02;
  Float_t noise = 4.0;

  if (NofMax == 1) {
    if (MaxOnStart == true) {
      C_M = (Float_t)ampls[0];
      C_D = (Float_t)ampls[1];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[2];
        if (C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255)
          shapecdtn = true;
        else if (exitCode)
          *exitCode = 2;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[2];
        C_Dnn = (Float_t)ampls[3];
        if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
            C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        } else if (exitCode)
          *exitCode = 3;
      }
    }

    if (MaxOnEnd == true) {
      C_M = (Float_t)ampls[ampls.size() - 1];
      C_D = (Float_t)ampls[ampls.size() - 2];
      if (ampls.size() < 3)
        shapecdtn = true;
      else if (ampls.size() == 3) {
        C_Dn = (Float_t)ampls[0];
        if (C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255)
          shapecdtn = true;
        else if (exitCode)
          *exitCode = 4;
      } else if (ampls.size() > 3) {
        C_Dn = (Float_t)ampls[ampls.size() - 3];
        C_Dnn = (Float_t)ampls[ampls.size() - 4];
        if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
            C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
          shapecdtn = true;
        } else if (exitCode)
          *exitCode = 5;
      }
    }

    if (MaxInMiddle == true) {
      C_M = (Float_t)ampls[MaxPos];
      int LeftOfMaxPos = MaxPos - 1;
      if (LeftOfMaxPos <= 0)
        LeftOfMaxPos = 0;
      int RightOfMaxPos = MaxPos + 1;
      if (RightOfMaxPos >= (int)ampls.size())
        RightOfMaxPos = ampls.size() - 1;
      //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl;
      if (ampls[LeftOfMaxPos] < ampls[RightOfMaxPos]) {
        C_D = (Float_t)ampls[RightOfMaxPos];
        C_Mn = (Float_t)ampls[LeftOfMaxPos];
        CDPos = RightOfMaxPos;
      } else {
        C_D = (Float_t)ampls[LeftOfMaxPos];
        C_Mn = (Float_t)ampls[RightOfMaxPos];
        CDPos = LeftOfMaxPos;
      }
      if (C_Mn < coeff1 * coeffn * C_M + coeff2 * coeffnn * C_D + 2 * noise || C_M == 255) {
        if (ampls.size() == 3)
          shapecdtn = true;
        else if (ampls.size() > 3) {
          if (CDPos > MaxPos) {
            if (ampls.size() - CDPos - 1 == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 == 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = 0;
            }
            if (ampls.size() - CDPos - 1 > 1) {
              C_Dn = (Float_t)ampls[CDPos + 1];
              C_Dnn = (Float_t)ampls[CDPos + 2];
            }
            if (MaxPos >= 2) {
              C_Mnn = (Float_t)ampls[MaxPos - 2];
            } else if (MaxPos < 2)
              C_Mnn = 0;
          }
          if (CDPos < MaxPos) {
            if (CDPos == 0) {
              C_Dn = 0;
              C_Dnn = 0;
            }
            if (CDPos == 1) {
              C_Dn = (Float_t)ampls[0];
              C_Dnn = 0;
            }
            if (CDPos > 1) {
              C_Dn = (Float_t)ampls[CDPos - 1];
              C_Dnn = (Float_t)ampls[CDPos - 2];
            }
            if (ampls.size() - LeftOfMaxPos > 1 && MaxPos + 2 < (int)(ampls.size()) - 1) {
              C_Mnn = (Float_t)ampls[MaxPos + 2];
            } else
              C_Mnn = 0;
          }
          if ((C_Dn <= coeff1 * coeffn * C_D + coeff2 * coeffnn * C_M + 2 * noise || C_D == 255) &&
              C_Mnn <= coeff1 * coeffn * C_Mn + coeff2 * coeffnn * C_M + 2 * noise &&
              C_Dnn <= coeff1 * coeffn * C_Dn + coeff2 * coeffnn * C_D + 2 * noise) {
            shapecdtn = true;
          }
        }
      } else if (exitCode)
        *exitCode = 6;
    }
  } else if (NofMax > 1 && exitCode)
    *exitCode = 1;  // more than one maximum
  if (ampls.size() == 1) {
    shapecdtn = true;
  }
  if (shapecdtn && exitCode)
    *exitCode = 0;

  return shapecdtn;
}

void TestNewCorr::Loop()
{
    if (fChain == 0) return -1;
    Long64_t nentries = fChain->GetEntries();
    Long64_t nbytes = 0, nb = 0;

    TFile* ofile = new TFile("TestOnSigdigi.root", "RECREATE");

    TH1F* QNewCorr_VS_Qtrue = new TH1F("QNewCorr_VS_Qtrue","QNewCorr_VS_Qtrue",150,0,2);
    QNewCorr_VS_Qtrue->GetXaxis()->SetTitle("Q^{new corr} / Q^{true}");
    TH1F* QNewCorrVec_VS_Qtrue = new TH1F("QNewCorrVec_VS_Qtrue","QNewCorrVec_VS_Qtrue",150,0,2);
    QNewCorrVec_VS_Qtrue->GetXaxis()->SetTitle("Q^{new corr vec} / Q^{true}");
    TH1F* ChargeOverPathlength_NewCorr = new TH1F("ChargeOverPathlength_NewCorr","ChargeOverPathlength_NewCorr",400,0,50);
    ChargeOverPathlength_NewCorr->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    TH1F* ChargeOverPathlength_OldCorr = new TH1F("ChargeOverPathlength_OldCorr","ChargeOverPathlength_OldCorr",400,0,50);
    ChargeOverPathlength_OldCorr->GetXaxis()->SetTitle("dE/dx [MeV/cm]");

    cout << endl;
    cout << "Histograms saved in: " << ofile->GetName() << endl;
    cout << endl;
    //nentries=100000;
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

                        // Only on W+jets MC
                        /*for (int istrip=sclus_index_strip[idedx]; istrip<sclus_index_strip[idedx]+sclus_nstrip[idedx]; istrip++)
                        {
                          Qsat.push_back(strip_ampl[istrip]);
                        }*/
                    
                        // Only on SingleMuon MC 
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
                            std::vector <int> QcorrOld = SaturationCorrection(Qsat, 0.10, 0.04, true, 20, 25);
                            
                            int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);
                            bool AreSameCluster = false;              
                            std::vector <int> QcorrVec = ReturnCorrVec(Qsat, layer, AreSameCluster);
                            int Qcorr = ReturnCorr(Qsat, layer);
                            
                            // Check if the correction is working (only if saturation) for:
                            //      - One single saturated strips
                            //      - Two consecutive saturated strips
                            int sumTrue = 0; //accumulate(Qtrue.begin(), Qtrue.end(), 0);
                            if (cpt_sat > 0)
                            {
                                QNewCorr_VS_Qtrue->Fill(1.*Qcorr/sumTrue);
                                QNewCorrVec_VS_Qtrue->Fill(1.*accumulate(QcorrVec.begin(), QcorrVec.end(), 0)/sumTrue);
                            }

                            std::vector <int> Qinv_OldCorr = CrossTalkInv(Qsat, 0.10, 0.04, true);
                            std::vector <int> Qinv_NewCorr = CrossTalkInvInStrip(Qsat, layer, true, 20, 0.10, 0.04);

                            if (clusterCleaning(Qinv_OldCorr, 1)) ChargeOverPathlength_OldCorr->Fill(1.*accumulate(QcorrOld.begin(), QcorrOld.end(), 0)*3.61e-06*265/dedx_pathlength[idedx]);
                            if (clusterCleaning(Qinv_NewCorr, 1)) ChargeOverPathlength_NewCorr->Fill(1.*accumulate(QcorrVec.begin(), QcorrVec.end(), 0)*3.61e-06*265/dedx_pathlength[idedx]);
                        }

                    } // end is_strip
                } // end dedx loop
            } // end selection on track
        } // end track loop

        if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100.0*jentry/nentries<<" %)"<<endl;

    } // end entry loop

    ofile->cd();
    QNewCorr_VS_Qtrue->Write();
    QNewCorrVec_VS_Qtrue->Write();
    ChargeOverPathlength_NewCorr->Write();
    ChargeOverPathlength_OldCorr->Write();
    ofile->Close();
}