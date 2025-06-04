#define CreateTemplate_cxx
#include "CreateTemplate.h"

#include "CorrFunctions.h"

using namespace TMath;
using namespace std;

void CreateTemplate::Loop()
{
    if (fChain == 0) return -1;
    Long64_t nentries = fChain->GetEntries();
    Long64_t nb = 0;
    gROOT->SetBatch(kTRUE);

        // SETUP
    TFile* ofile = new TFile("Check_Template.root", "RECREATE");

    ofstream Template_FLFR("Template_correction/Template_FLFR.txt", std::ofstream::out);
    ofstream Template_CENTER("Template_correction/Template_CENTER.txt", std::ofstream::out);
    ofstream Template_FirstNeighbour("Template_correction/Template_LEFTRIGHT.txt", std::ofstream::out);
    
    // histograms to check the templates
    std::vector <TH1D*> FullLeft_template_TH1;
    std::vector <TH1D*> FullRight_template_TH1;
    std::vector <TH2D*> LeftRight_template_TH2;
    std::vector <TH1D*> Center_template_TH1;

    for (int i=0; i<=20; i++)
    {
        TString histName = Form("FL_%d", i);
        TH1D *histo = new TH1D(histName,histName,200,0,20);
        histo->GetXaxis()->SetTitle("Max/Vmax");
        FullLeft_template_TH1.push_back(histo);

        histName = Form("FR_%d", i);
        TH1D *histo2 = new TH1D(histName,histName,200,0,20);
        histo2->GetXaxis()->SetTitle("Max/Vmax");
        FullRight_template_TH1.push_back(histo2);

        histName = Form("LR_%d", i);
        TH2D *histo3 = new TH2D(histName,histName,100,0,10,200,0,20);
        histo3->GetXaxis()->SetTitle("Vmax/Vmin");
        histo3->GetYaxis()->SetTitle("Max/Vmax");
        LeftRight_template_TH2.push_back(histo3);

        histName = Form("Center_%d", i);
        TH1D *histo4 = new TH1D(histName,histName,200,0,20);
        histo4->GetXaxis()->SetTitle("Max/Vmax");
        Center_template_TH1.push_back(histo4);
    }
    
    // Template First Neighbour at a coarse cut
    ifstream Template_FirstNeighbour_cut("Template_correction/Template_LEFTRIGHT_firstcut.txt");
    std::string line;
    std::vector<double> template_a1;
    std::vector<double> template_a2;
    while (std::getline(Template_FirstNeighbour_cut, line))
    {
        std::istringstream iss(line);
        int i_layer=-1;
        double a1=-1, a2=-1;
        if (!(iss >> i_layer >> a1 >> a2 )) break; // error
        
        template_a1.push_back(a1);
        template_a2.push_back(a2);
    }
    Template_FirstNeighbour_cut.close();
    

    cout << endl;
    cout << "Histograms saved in: " << ofile->GetName() << endl;
    cout << endl;

    double DeltaR=0;
    double DeltaR_val;
    int index;
    //nentries = 20000;


        // ENTRY LOOP
    for (Long64_t jentry=0; jentry<nentries;jentry++)
    {
        nb = fChain->GetEntry(jentry);

            // MUON MATCHING
        int gen_index=0;
        for (int igen=0; igen<ngenpart; igen++)
        {
            if (gen_pdg[igen]==13 && gen_status[igen]==1) gen_index=igen;
        }

        DeltaR_val=100;
        index=-1;
        for (int itr=0; itr<ntracks; itr++)
        {
            DeltaR=Sqrt(Power(modulo2Pi(track_phi[itr]-gen_phi[gen_index]),2)+Power(track_eta[itr]-gen_eta[gen_index],2));

            if(DeltaR<DeltaR_val)
            {
                DeltaR_val=DeltaR;
                index=itr;
            }
        }


        if (track_nvalidhits[index]>8 && track_qual[index])
        {
            for (int idedx=track_index_hit[index]; idedx<track_index_hit[index]+track_nhits[index]; idedx++)
            {
                if (dedx_isstrip[idedx])
                {
                    std::vector<int> Q_sig;
                    bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
                    
                    int layer = FindLayer(dedx_subdetid[idedx], dedx_detid[idedx]);

                    for (int isigdigi=0;isigdigi<ndigi_sig_strip;isigdigi++)
                    {
                        if (sigdigi_strip_id[isigdigi]==dedx_detid[idedx] && sigdigi_strip_channel[isigdigi]>= strip_channel[sclus_index_strip[idedx]] && sigdigi_strip_channel[isigdigi]<= strip_channel[sclus_index_strip[idedx] + sclus_nstrip[idedx]-1])
                        {
                            Q_sig.push_back(sigdigi_strip_adc[isigdigi]);
                        }
                    }

                    if (!Q_sig.empty())
                    {
                        int i_max = find(Q_sig.begin(), Q_sig.end(), *max_element(Q_sig.begin(), Q_sig.end())) - Q_sig.begin();
                        ClusterShape(Q_sig, left, right, center, FullLeft, FullRight);


                            // FL AND FR TEMPLATES
                        if (FullLeft && ( ( Q_sig[i_max+1]>=16 && (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==4 || (layer>=14 && layer<=17)) )
                        || ( Q_sig[i_max+1]>=25 && (dedx_subdetid[idedx]==5 || layer>=18) ) ) )
                        {
                            FullLeft_template_TH1[layer]->Fill(1.0* Q_sig[i_max]/Q_sig[i_max+1]);  
                        }

                        if (FullRight && ( ( Q_sig[i_max-1]>=16 && (dedx_subdetid[idedx]==3 || dedx_subdetid[idedx]==4 || (layer>=14 && layer<=17)) )
                        || ( Q_sig[i_max-1]>=25 && (dedx_subdetid[idedx]==5 || layer>=18) ) ) )
                        {
                            FullRight_template_TH1[layer]->Fill(1.0* Q_sig[i_max]/Q_sig[i_max-1]);
                        }

                            // LEFT-RIGHT TEMPLATES
                        if (left || right)
                        {
                            if (Q_sig[i_max-1]>=Q_sig[i_max+1])
                            {
                                if (1.0*Q_sig[i_max] >= 0.85*(template_a1[layer-1]*Q_sig[i_max-1] + template_a2[layer-1]*Q_sig[i_max+1]) )
                                LeftRight_template_TH2[layer]->Fill(1.0* Q_sig[i_max-1]/Q_sig[i_max+1], 1.0* Q_sig[i_max]/Q_sig[i_max-1]);
                            }
                            else
                            {
                                if (1.0*Q_sig[i_max] >= 0.85*(template_a1[layer-1]*Q_sig[i_max+1] + template_a2[layer-1]*Q_sig[i_max-1]) ) 
                                LeftRight_template_TH2[layer]->Fill(1.0* Q_sig[i_max+1]/Q_sig[i_max-1], 1.0*Q_sig[i_max]/Q_sig[i_max+1]);
                            }
                        }

                            // CENTER TEMPLATES
                        if (center) Center_template_TH1[layer]->Fill( 2.*Q_sig[i_max] / (Q_sig[i_max+1] + Q_sig[i_max-1]));

                    } //end if index_sig
                    } //end if strip
            } //end hit loop
        } //end track loop

        if (jentry%20000==0) cout<<jentry<<" / "<<nentries<<" ("<<100.0* jentry/nentries<<" %)"<<endl;
    } // end entry loop

        // FL AND FR TEMPLATES
    for (int i=1;i<FullLeft_template_TH1.size();i++) Template_FLFR << i << " " << getMaxHisto_X(FullLeft_template_TH1[i]) << " " << getMaxHisto_X(FullRight_template_TH1[i]) << "\n";
    Template_FLFR.close();
    

        // LEFT-RIGHT TEMPLATES
    std::vector <float> p0;
    std::vector <float> p1;
    for (int i=0; i<LeftRight_template_TH2.size(); i++)
    {
        TF1 *fit = new TF1("fit","[0] + [1]/x",1,10);
        LeftRight_template_TH2[i]->Fit("fit","REMQ");

        p0.push_back(fit->GetParameter(0));
        p1.push_back(fit->GetParameter(1));
    }
    for (int i=1; i<p0.size(); i++) Template_FirstNeighbour << i << " " << p0[i] << " " << p1[i] << "\n";
    Template_FirstNeighbour.close();

        // CENTER TEMPLATES
    for (int i=1;i<Center_template_TH1.size();i++) Template_CENTER << i << " " << getMaxHisto_X(Center_template_TH1[i]) << "\n";
    Template_CENTER.close();


        // OUTPUT CHECK 
    ofile->cd();
    for (int i=1; i<LeftRight_template_TH2.size(); i++) LeftRight_template_TH2[i]->Write();
    for (int i=1; i<Center_template_TH1.size(); i++) Center_template_TH1[i]->Write();
    for (int i=0; i<FullLeft_template_TH1.size(); i++) {FullRight_template_TH1[i]->Write(); FullLeft_template_TH1[i]->Write();}
    ofile->Close();

    return;
}

