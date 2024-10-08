#define M_PI   3.14159265358979323846  /* pi */

    // GLOBAL FUNCTIONS
int FindLayer(int subdetid, int detid)
{
    if(subdetid==3)  // TIB
    {
        if(((detid>>14)&0x7)==1) return 1;
        else if(((detid>>14)&0x7)==2) return 2;
        else if(((detid>>14)&0x7)==3) return 3;
        else if(((detid>>14)&0x7)==4) return 4;
    }    
    else if(subdetid==5) // TOB
    {       
        if(((detid>>14)&0x7)==1) return 5;
        else if(((detid>>14)&0x7)==2) return 6;
        else if(((detid>>14)&0x7)==3) return 7;
        else if(((detid>>14)&0x7)==4) return 8;
        else if(((detid>>14)&0x7)==5) return 9;
        else if(((detid>>14)&0x7)==6) return 10;
    }       
    /*else if(subdetid==4)  //TID by disk 
    {    
        if(((detid>>11)&0x3)==1) return 11;
        else if(((detid>>11)&0x3)==2) return 12;
        else if(((detid>>11)&0x3)==3) return 13;
    }*/
    else if(subdetid==4)  //TID by ring 
    {    
        if(((detid>>9)&0x3)==1) return 11;
        else if(((detid>>9)&0x3)==2) return 12;
        else if(((detid>>9)&0x3)==3) return 13;
    }
    /*else if(subdetid==6) // TEC by wheel 
    {
        if(((detid>>14)&0xF)==1) return 14;
        else if(((detid>>14)&0xF)==2) return 15;
        else if(((detid>>14)&0xF)==3) return 16;
        else if(((detid>>14)&0xF)==4) return 17;
        else if(((detid>>14)&0xF)==5) return 18;
        else if(((detid>>14)&0xF)==6) return 19;
        else if(((detid>>14)&0xF)==7) return 20;
        else if(((detid>>14)&0xF)==8) return 21;
        else if(((detid>>14)&0xF)==9) return 22;
    }*/
    else if(subdetid==6) // TEC by ring 
    {
        if(((detid>>5)&0x7)==1) return 14;
        else if(((detid>>5)&0x7)==2) return 15;
        else if(((detid>>5)&0x7)==3) return 16;
        else if(((detid>>5)&0x7)==4) return 17;
        else if(((detid>>5)&0x7)==5) return 18;
        else if(((detid>>5)&0x7)==6) return 19;
        else if(((detid>>5)&0x7)==7) return 20;
        else return -1;
    }
    return -1;
}

Double_t modulo2Pi(Double_t nb)
{
   if (nb>0)
   {
      while(nb > M_PI) nb -= 2*M_PI;
   }
   else
   {
      while(nb < -M_PI) nb += 2*M_PI;
   }
   return nb;
}

double getMaxHisto_X(TH1D* h)
{
  return  h->GetMaximumBin() * (h->GetXaxis()->GetBinCenter(h->GetNbinsX()) - h->GetXaxis()->GetBinCenter(1)) * 1./h->GetNbinsX();
}



    // CORRECTION FUNCTIONS
int Correction_FL_FR(const std::vector <int>&  Q, int layer, std::string TemplateFile)
{
    int MaxCorr = 0;
    int N_sat = 0;
    int ThresholdSat = -1;

    if ((layer>=5 && layer<=10) || layer>=18) ThresholdSat = 25;
    else ThresholdSat = 16;

    // CORRECTION TEMPLATES
    std::string line;
    std::vector<double> template_FL;
    std::vector<double> template_FR;

    std::ifstream Template(TemplateFile);
    while (std::getline(Template, line))
    {
        std::istringstream iss(line);
        int i_layer=-1;
        double FL_a=-1, FR_a=-1;
        if (!(iss >> i_layer >> FL_a >> FR_a)) break; // error
            
        template_FL.push_back(FL_a);
        template_FR.push_back(FR_a);
    }
    Template.close();

    // NUMBER OF MAX
    for (unsigned int i=0; i<Q.size(); i++)
    {
        if (Q[i] >= 254) N_sat++;
    }

    int max_Q = *max_element(Q.begin(), Q.end());
    int i_max = find(Q.begin(), Q.end(), max_Q) - Q.begin();

    // NO CORRECTION IF TOO SMALL OR LARGE
    if (Q.size()<2 || Q.size()>8) return max_Q;

    if (N_sat == 1)
    {
        int MaxCorr = -1;
        // FULL LEFT
        if (max_Q==Q[0] && Q[i_max+1]>=ThresholdSat && Q[i_max+1]<254)
        {
            MaxCorr = template_FL[layer-1] * Q[i_max+1];

            if (max_Q == 255 && MaxCorr < 255) return 255;
            if (MaxCorr < 254) return 254;
            else return MaxCorr;
        }

        // FULL RIGHT
        if (max_Q==Q[Q.size()-1] && Q[i_max-1]>=ThresholdSat && Q[i_max-1]<254)
        {
            MaxCorr = template_FR[layer-1] * Q[i_max-1];

            if (max_Q == 255 && MaxCorr < 255) return 255;
            if (MaxCorr < 254) return 254;
            else return MaxCorr;
        }

        return max_Q;
    }

    // NO ONE SINGLE SATURATED STRIP --> no x-talk inversion
    return max_Q;
}

int Correction_LRC(const std::vector <int>&  Q, int layer, std::string TemplateFile, bool CENTER)
{
	// SETUP 
	int N_sat = 0;
	int Vmin = 0, Vmax = 0;

    for (unsigned int i=0; i<Q.size(); i++)
    {
        if (Q[i]>=254) N_sat++;
    }

	int max_Q = *max_element(Q.begin(), Q.end());
    int i_max = find(Q.begin(), Q.end(), max_Q) - Q.begin(); 

    if (Q[i_max-1] > Q[i_max+1])
    {
        Vmax = Q[i_max-1];
        Vmin = Q[i_max+1];
    }
    else
    {
        Vmax = Q[i_max+1];
        Vmin = Q[i_max-1];
    }

    // NO CORRECTION IF TOO SMALL OR LARGE
    if (Q.size()<2 || Q.size()>8) return max_Q;

    // CORRECTION TEMPLATES
    std::ifstream Template(TemplateFile);
    std::string line;
    std::vector<double> template_a1;
    std::vector<double> template_a2;
    
    if (CENTER)
    {
        while (std::getline(Template, line))
        {
            std::istringstream iss(line);
            int i_layer=-1;
            double a1=-1;
            if (!(iss >> i_layer >> a1)) break; // error
                
            template_a1.push_back(a1);
        }
        Template.close();

        int MaxCorr = 1.*(Vmax+Vmin)/2 * template_a1[layer-1];

        if (max_Q == 255 && MaxCorr < 255) return 255;
        if (MaxCorr < 254) return 254;
        else return MaxCorr;
    }

    while (std::getline(Template, line))
    {
        std::istringstream iss(line);
        int i_layer=-1;
        double a1=-1, a2=-1;
        if (!(iss >> i_layer >> a1 >> a2 )) break; // error
        
        template_a1.push_back(a1);
        template_a2.push_back(a2);
    }
    Template.close();

    int MaxCorr = Vmax * template_a1[layer-1] + Vmin * template_a2[layer-1];
    if (max_Q == 255 && MaxCorr < 255) return 255;
    if (MaxCorr < 254) return 254;
    else return MaxCorr;
}

int Correction_2strips(const std::vector <int>&  Q)
{
    int Qcorr = -1;
    int SumQ = accumulate(Q.begin(), Q.end(), 0);
    int index_maxLeft = -1;

        // NUMBER OF MAX + IF CONSECUTIVE OR NOT
    int N_sat = 0;
    bool two_cons_sat = false;

    if (Q[0]>=254) N_sat++;
    for (int i=1; i<Q.size(); i++)
    {
        if (Q[i]>=254) N_sat++;

        if (Q[i-1]>=254 && Q[i]>=254)
        {
            if (two_cons_sat) continue;
            two_cons_sat = true;
            index_maxLeft = i-1; 
        }
    }
    
        // NO CORRECTION if not 2 consecutive saturated strips, or if clusters too small or if max is at the edge
    if (N_sat!=2 || !two_cons_sat || Q.size() < 4
        || index_maxLeft==0 || index_maxLeft==Q.size()-2) return SumQ;
    
        // CORRECTION
    Qcorr = 1./0.0613 * (Q[index_maxLeft-1] + Q[index_maxLeft+2])/2;

    if (Qcorr > SumQ) return Qcorr;
    else return SumQ;
}

void ClusterShape(const std::vector <int>&  Q, bool &left, bool &right, bool &center, bool &FullLeft, bool &FullRight)
{
    // SETUP
    int max_Q = *max_element(Q.begin(), Q.end());
    int i_max = find(Q.begin(), Q.end(), max_Q) - Q.begin();

    int Nleft = -1, Nright = -1;
    if (Q.size()>=3) {Nleft = Q[i_max-1]; Nright = Q[i_max+1];}
    
    int N_sat = 0;
    for (int i=0; i<Q.size(); i++) { if (Q[i]>=254) N_sat++; }

    // SHAPE
    if (Q.size()>=3 && N_sat==1 && i_max>0 && i_max<Q.size()-1 && Nleft>1.1*Nright) left=true;
    if (Q.size()>=3 && N_sat==1 && i_max>0 && i_max<Q.size()-1 && Nleft<0.9*Nright) right=true;
    if (Q.size()>=3 && N_sat==1 && i_max>0 && i_max<Q.size()-1 && Nleft<=1.1*Nright && Nleft>=0.9*Nright) center=true;
    if (Q.size()>=2 && N_sat==1 && i_max==0) FullLeft=true;
    if (Q.size()>=2 && N_sat==1 && i_max==Q.size()-1) FullRight=true;

    return;
}

int ReturnCorr(const std::vector <int>& Q, int subdetid, int detid)
{
    // SETUP
    bool left=false, right=false, center=false, FullLeft=false, FullRight=false;
    int MaxCorr = -1;

    // SHAPE
    ClusterShape(Q, left, right, center, FullLeft, FullRight);
    

    // CORRECTION
        // 2 consecutive saturated strips
    MaxCorr = Correction_2strips(Q);
    if (MaxCorr > accumulate(Q.begin(), Q.end(), 0)) return MaxCorr;
    
    MaxCorr = *max_element(Q.begin(), Q.end());
        // 1 saturated strip
    if (center) MaxCorr = Correction_LRC(Q, FindLayer(subdetid, detid), "Template_correction/Template_CENTER.txt", true);
    else if (left || right) MaxCorr = Correction_LRC(Q, FindLayer(subdetid, detid), "Template_correction/Template_LEFTRIGHT.txt", false);
    else if (FullLeft || FullRight) MaxCorr = Correction_FL_FR(Q, FindLayer(subdetid, detid), "Template_correction/Template_FLFR.txt");
    

    // SUM CORRECTION
    int sum_Qcorr = 0;
    int i_max = find(Q.begin(), Q.end(), *max_element(Q.begin(), Q.end())) - Q.begin();
    for (int i=0; i<Q.size(); i++)
    {
        if (i != i_max) sum_Qcorr += Q[i];
    }
    sum_Qcorr += MaxCorr;

    return sum_Qcorr;
}

