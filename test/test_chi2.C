// # Copyright 2023  Adrián Irles (IFIC)

#include "../include/experimentalUnc.h"
#include "../include/analysis.h"
#include "../style/Style.C"
#include "../style/Labels.C"

void test_one_quark_one_energy(int quark_under_test=4, float energy_under_test=250)
{

    TString quarkst="c-quark";
    if(quark_under_test==5) quarkst="b-quark";
    TString energyst="250 GeV";
    if(energy_under_test==500) energyst="500 GeV";

    read_all_models(false);
    // stores the information for all  models (named in std::vector<TString> allmodels={"SMA","A1"};)
    // with 100% polarization options only
    // it stores the same info as the txt files
    // in a vector std::vector<model_struct_t> theory;
    // to access the info, you do: theory.at(i).model_values[ie].AFB_L[ifl]
    // indexes:
    //   i: as many as model strings
    //   ie: index for the energy: one index for every energy found in the file. If you want to know the value, you can acces to i: theory.at(i).model_values[ie].energy
    //       we can store up to 5 energy points
    //   ifl: for the flavour , 0 is empty, 1,2,3,4,5 are d,u,s,c,b


    // here I create the theoretical observables for the desired polarization scheme
    // I save a vector for all energies, all pol (2 in this case) and and the b and c flavours
    std::vector<observables_struct_t> observables=create_observables();

    // I calculate the chi2 matrix of tested models vs hypothesis models, using th2f
    std::vector<TH2F*> chi2histos=chi2_ndf_models(quark_under_test, energy_under_test, observables) ;
    
    //create final histograms
    TH2F * h_results[9];
    h_results[0] = new TH2F("AFB_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), %s #sqrt{s}=%s; ; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[1] = new TH2F("AFB_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[2] = new TH2F("AFB_comb",TString::Format("AFB_ILD_H20, %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[3] = new TH2F("R_L",TString::Format("R_ILD_H20 (-0.8,+0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[4] = new TH2F("R_R",TString::Format("R_ILD_H20 (+0.8,-0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[5] = new TH2F("R_comb",TString::Format("R_ILD_H20, %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[6] = new TH2F("AFB_&_R_L",TString::Format("AFB_&_R_ILD_H20 (-0.8,+0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[7] = new TH2F("AFB_&_R_R",TString::Format("AFB_&_R_ILD_H20 (+0.8,-0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[8] = new TH2F("AFB_&_R_comb",TString::Format("AFB_&_R_ILD_H20, %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);


   TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, full, %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),
        TString::Format("H20, P:(-0.8,+0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, full, %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),
        TString::Format("H20, P:(-0.8,+0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, full, %s #sqrt{s}=%s",quarkst.Data(),energyst.Data())};

    TString titles[]={"AFB_{q}","AFB_{q}","AFB_{q}","R_{q}","R_{q}","R_{q}","AFB_{q} and R_{q}","AFB_{q} and R_{q}","AFB_{q} and R_{q}"};

    for(int i=0; i<9; i++) {
     for(int j=0; j<theory.size(); j++) {
        h_results[i]->GetXaxis()->SetBinLabel(j+1,allmodels.at(j));
        h_results[i]->GetYaxis()->SetBinLabel(j+1,allmodels.at(j));
        h_results[i]->GetYaxis()->SetLabelSize(0.07);
        h_results[i]->GetXaxis()->SetLabelSize(0.07);
     }
    }  

   
    TCanvas *c1 = new TCanvas ("c1","c1",1200,800);
   // Int_t colors[] = {kRed+3, kRed+1, kMagenta+3, kBlue, kCyan+2, kTeal, kGreen}; // #colors >= #levels - 1
    //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Int_t colors[] = {kGray+1, kGray+2,kGray+3,kYellow-6, kSpring+5, kSpring-1, kSpring}; // #colors >= #levels - 1
    TString colors_st[] = {"<1#sigma","1-2#sigma","2-3#sigma","3-4#sigma","4-5#sigma","5-6#sigma","#geq 6#sigma"}; // #colors >= #levels - 1
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    
    c1->Divide(3,3);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    for(int i=0; i<9; i++) {
        c1->cd(i+1);
        //gPad->SetLogz();
        for(int j=0; j<theory.size(); j++) {
            for(int k=0; k<theory.size(); k++){
                if(j>=k) continue;
                float ndfv=1; // THIS IS ONLY VALID FOR THE CASE OF ONE QUARK AND ONE ENERGY !!

                float chi2v=0;
                if(i>5) {
                    chi2v=chi2histos.at(i-6)->GetBinContent(j+1,k+1)+chi2histos.at(i-3)->GetBinContent(j+1,k+1);
                    ndfv*=2;
                }
                else chi2v=chi2histos.at(i)->GetBinContent(j+1,k+1);

                float nsigma=0;
                if(TMath::Prob(chi2v/(ndfv),ndfv)>69.1/100.) nsigma=0.1;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<69.1/100. && TMath::Prob(chi2v/(ndfv),ndfv)>30.9/100.) nsigma=1;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<30.9/100. && TMath::Prob(chi2v/(ndfv),ndfv)>6.7/100.) nsigma=2;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<6.7/100. && TMath::Prob(chi2v/(ndfv),ndfv)>0.62/100.) nsigma=3;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<0.62/100. && TMath::Prob(chi2v/(ndfv),ndfv)>0.023/100.) nsigma=4;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<0.023/100.&& TMath::Prob(chi2v/(ndfv),ndfv)>0.00034/100) nsigma=5;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<0.00034/100.) nsigma=6;
                //if(TMath::Prob(chi2v/(ndfv),ndfv)<0.00034/100.) nsigma=6;
                h_results[i]->SetMarkerSize(2);
                h_results[i]->SetBinContent(j+1,k+1,nsigma);
            }
        }
        h_results[i]->SetTitle(titles[i]);
        h_results[i]->GetZaxis()->SetRangeUser(0.01,6.1);
        h_results[i]->Draw("col"    );
        QQBARLabel2(0.25,0.15, conditions[i],kBlack,0.055);

        TLegend *leg = new TLegend(0.7,0.3,0.85,0.72,"","blNDC");
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);
        TH2F* h_leg0[10];
        for(int ic=0; ic<7;ic++ ){
            h_leg0[ic]=new TH2F;
            h_leg0[ic]->SetFillStyle(1000);
            h_leg0[ic]->SetFillColor(colors[ic]);
            h_leg0[ic]->SetLineColor(colors[ic]);
            leg->AddEntry(h_leg0[ic],colors_st[ic],"f");
        }        
        // leg->AddEntry(bkg_zz_0,"ZZ","lf");
        
        // leg->SetFillStyle(1001);
        leg->SetFillColor(0);
        // leg->SetLineColor(1);
        leg->SetBorderSize(0);
        leg->Draw();
        //h_chi2[i]->Draw("col"    );
    }


}

// THIS IS TO BE CHECKED... 
// the combination is not well done!!!
void test_two_quarks_two_energies()
{

    TString quarkst="c&b-quark";
    TString energyst="250&500 GeV";

    read_all_models(false);
    // stores the information for all  models (named in std::vector<TString> allmodels={"SMA","A1"};)
    // with 100% polarization options only
    // it stores the same info as the txt files
    // in a vector std::vector<model_struct_t> theory;
    // to access the info, you do: theory.at(i).model_values[ie].AFB_L[ifl]
    // indexes:
    //   i: as many as model strings
    //   ie: index for the energy: one index for every energy found in the file. If you want to know the value, you can acces to i: theory.at(i).model_values[ie].energy
    //       we can store up to 5 energy points
    //   ifl: for the flavour , 0 is empty, 1,2,3,4,5 are d,u,s,c,b


    // here I create the thheoretical observables for the desired polarization scheme
    // I save a vector for all energies, all pol (2 in this case) and and the b and c flavours
    std::vector<observables_struct_t> observables=create_observables();

    // I calculate the chi2 matrix of tested models vs hypothesis models, using th2f
    std::vector<TH2F*> chi2histos_c=chi2_ndf_models(4, 250, observables) ;
    std::vector<TH2F*> chi2histos_b=chi2_ndf_models(5, 250, observables) ;
    std::vector<TH2F*> chi2histos_c500=chi2_ndf_models(4, 500, observables) ;
    std::vector<TH2F*> chi2histos_b500=chi2_ndf_models(5, 500, observables) ;

    //create final histograms
    TH2F * h_results[9];
    h_results[0] = new TH2F("AFB_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), %s #sqrt{s}=%s; ; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[1] = new TH2F("AFB_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[2] = new TH2F("AFB_comb",TString::Format("AFB_ILD_H20, %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[3] = new TH2F("R_L",TString::Format("R_ILD_H20 (-0.8,+0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[4] = new TH2F("R_R",TString::Format("R_ILD_H20 (+0.8,-0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[5] = new TH2F("R_comb",TString::Format("R_ILD_H20, %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[6] = new TH2F("AFB_&_R_L",TString::Format("AFB_&_R_ILD_H20 (-0.8,+0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[7] = new TH2F("AFB_&_R_R",TString::Format("AFB_&_R_ILD_H20 (+0.8,-0.3), %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);
    h_results[8] = new TH2F("AFB_&_R_comb",TString::Format("AFB_&_R_ILD_H20, %s #sqrt{s}=%s; ; ;",quarkst.Data(),energyst.Data()),9,-0.5,8.5,9,-0.5,8.5);


    TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, full, %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),
        TString::Format("H20, P:(-0.8,+0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, full, %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),
        TString::Format("H20, P:(-0.8,+0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), %s #sqrt{s}=%s",quarkst.Data(),energyst.Data()),TString::Format("H20, full, %s #sqrt{s}=%s",quarkst.Data(),energyst.Data())};

    TString titles[]={"AFB_{q}","AFB_{q}","AFB_{q}","R_{q}","R_{q}","R_{q}","AFB_{q} and R_{q}","AFB_{q} and R_{q}","AFB_{q} and R_{q}"};

    for(int i=0; i<9; i++) {
     for(int j=0; j<theory.size(); j++) {
        h_results[i]->GetXaxis()->SetBinLabel(j+1,allmodels.at(j));
        h_results[i]->GetYaxis()->SetBinLabel(j+1,allmodels.at(j));
        h_results[i]->GetYaxis()->SetLabelSize(0.07);
        h_results[i]->GetXaxis()->SetLabelSize(0.07);
     }
    }  

   
    TCanvas *c1 = new TCanvas ("c1","c1",1200,800);
   // Int_t colors[] = {kRed+3, kRed+1, kMagenta+3, kBlue, kCyan+2, kTeal, kGreen}; // #colors >= #levels - 1
    //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   // gStyle->SetPalette(kIsland);
    Int_t colors[] = {kGray+1, kGray+2,kGray+3,kYellow-6, kSpring+5, kSpring-1, kSpring}; // #colors >= #levels - 1
    TString colors_st[] = {"<1#sigma","1-2#sigma","2-3#sigma","3-4#sigma","4-5#sigma","5-6#sigma","#geq 6#sigma"}; // #colors >= #levels - 1
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    c1->Divide(3,3);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    for(int i=0; i<9; i++) {
        h_results[i]->SetTitle(titles[i]);
        c1->cd(i+1);
        //gPad->SetLogz();
        for(int j=0; j<theory.size(); j++) {
            for(int k=0; k<theory.size(); k++){
                if(j>=k) continue;
                float ndfv=4; // THIS IS ONLY VALID FOR THE CASE OF two QUARKs AND two ENERGies !!

                float chi2v=0;
                if(i>5 ) {
                    chi2v=chi2histos_c.at(i-6)->GetBinContent(j+1,k+1)+chi2histos_b.at(i-6)->GetBinContent(j+1,k+1)+chi2histos_c500.at(i-6)->GetBinContent(j+1,k+1)+chi2histos_b500.at(i-6)->GetBinContent(j+1,k+1)
                    +chi2histos_c.at(i-3)->GetBinContent(j+1,k+1)+chi2histos_b.at(i-3)->GetBinContent(j+1,k+1)+chi2histos_c500.at(i-3)->GetBinContent(j+1,k+1)+chi2histos_b500.at(i-3)->GetBinContent(j+1,k+1);                   
                    ndfv*=2;
                } else chi2v=chi2histos_c.at(i)->GetBinContent(j+1,k+1)+chi2histos_b.at(i)->GetBinContent(j+1,k+1)+chi2histos_c500.at(i)->GetBinContent(j+1,k+1)+chi2histos_b500.at(i)->GetBinContent(j+1,k+1);
                

                float nsigma=0;
                if(TMath::Prob(chi2v/(ndfv),ndfv)>69.1/100.) nsigma=0.1;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<69.1/100. && TMath::Prob(chi2v/(ndfv),ndfv)>30.9/100.) nsigma=1;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<30.9/100. && TMath::Prob(chi2v/(ndfv),ndfv)>6.7/100.) nsigma=2;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<6.7/100. && TMath::Prob(chi2v/(ndfv),ndfv)>0.62/100.) nsigma=3;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<0.62/100. && TMath::Prob(chi2v/(ndfv),ndfv)>0.023/100.) nsigma=4;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<0.023/100.&& TMath::Prob(chi2v/(ndfv),ndfv)>0.00034/100) nsigma=5;
                if(TMath::Prob(chi2v/(ndfv),ndfv)<0.00034/100.) nsigma=6;
                //if(TMath::Prob(chi2v/(ndfv),ndfv)<0.00034/100.) nsigma=6;
                h_results[i]->SetMarkerSize(2);
                h_results[i]->SetBinContent(j+1,k+1,nsigma);
            }
        }
        h_results[i]->GetZaxis()->SetRangeUser(0.01,6.1);
        h_results[i]->Draw("col");
        QQBARLabel2(0.25,0.15, conditions[i],kBlack,0.055);

        TLegend *leg = new TLegend(0.7,0.3,0.85,0.72,"","blNDC");
        leg->SetTextFont(42);
        leg->SetTextSize(0.04);
        TH2F* h_leg0[10];
        for(int ic=0; ic<7;ic++ ){
            h_leg0[ic]=new TH2F;
            h_leg0[ic]->SetFillStyle(1000);
            h_leg0[ic]->SetFillColor(colors[ic]);
            h_leg0[ic]->SetLineColor(colors[ic]);
            leg->AddEntry(h_leg0[ic],colors_st[ic],"f");
        }        
        // leg->AddEntry(bkg_zz_0,"ZZ","lf");
        
        // leg->SetFillStyle(1001);
        leg->SetFillColor(0);
        // leg->SetLineColor(1);
        leg->SetBorderSize(0);
        leg->Draw();

        //h_chi2[i]->Draw("col"    );
    }


}

void test_chi2(){

    test_one_quark_one_energy(5,500);
    //test_two_quarks_two_energies();

}

