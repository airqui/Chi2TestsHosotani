// # Copyright 2023  Adrián Irles (IFIC)

#include "../include/models_struct.h"
#include "../include/experimentalUnc.h"

void read_all_models(bool debug = true)
{

    //uses std::vector<TString> allmodels = {"SMB","A1","A2","A3","B1","B2","B3","B4"}; defined in models_struct
    //std::vector<model_struct_t> theory;

    for (int i = 0; i < allmodels.size(); i++)
    {
        model_struct_t newmodel = read_model(allmodels.at(i), i, debug);
        theory.push_back(newmodel);
    }

    if (debug == true)
    {

        for (int i = 0; i < theory.size(); i++)
        {
            cout << theory.at(i).modelname << " index:" << theory.at(i).id << endl;
            for (int ie = 0; ie < 5; ie++)
                for (int ifl = 1; ifl < 6; ifl++)
                    cout << theory.at(i).model_values[ie].energy << " " << theory.at(i).model_values[ie].cross_L[ifl] << " " << theory.at(i).model_values[ie].cross_R[ifl]
                         << " " << theory.at(i).model_values[ie].AFB_L[ifl] << " " << theory.at(i).model_values[ie].AFB_R[ifl] << endl;
        }
    }
}

void test2()
{

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


    // hereI create the thheoretical observables for the desired polarization scheme
    // I save a vector for all energies, all pol (2 in this case) and and the b and c flavours
    std::vector<observables_struct_t> observables;

    for (int imodel = 0; imodel < theory.size(); imodel++)
    {
        for (int ienergy = 0; ienergy < 5; ienergy++)
        {
            //if (theory.at(imodel).model_values[ienergy].energy != 250)
            //    continue;

            for (int iflav = 4; iflav < 6; iflav++)
            {
                observables_struct_t newobs;
                newobs.id = theory.at(imodel).id;
                newobs.modelname = theory.at(imodel).modelname;
                newobs.energy = theory.at(imodel).model_values[ienergy].energy;
                newobs.flav = iflav;
                newobs.Cross_L = ObsCross(imodel, ienergy, iflav, -0.8,0.3);
                newobs.Cross_R = ObsCross(imodel, ienergy, iflav, 0.8,-0.3);
                newobs.Cross_unpol = ObsCross(imodel, ienergy, iflav, 0,0);
                newobs.R_L = ObsR(imodel, ienergy, iflav, -0.8,0.3);
                newobs.AFB_L = ObsAFB(imodel, ienergy, iflav, -0.8,0.3);
                newobs.R_R = ObsR(imodel, ienergy, iflav, 0.8,-0.3);
                newobs.AFB_R = ObsAFB(imodel, ienergy, iflav,0.8,-0.3);
                newobs.R_unpol = ObsR(imodel, ienergy, iflav, 0,0);
                newobs.AFB_unpol = ObsAFB(imodel, ienergy, iflav, 0,0);
                observables.push_back(newobs);
            }
        }
    }
    //************************************************

    //EXAMPLE ANALYSIS -> Comparation of all models, only b-quark and only 250GeV
    int quark_under_test=4;//if -1, then all v.alues are used
    float energy_under_test=500;//if -1, then all values are used
    int rr=100000; //npseudo experiments

    //create final histograms
    TH2F * h_results[9];
    h_results[0] = new TH2F("AFB_L","AFB_L_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[1] = new TH2F("AFB_R","AFB_R_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[2] = new TH2F("AFB_L_R","AFB_L+R_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[3] = new TH2F("R_L","R_L_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[4] = new TH2F("R_R","R_R_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[5] = new TH2F("R_L_R","R_L+R_ILC; Reference Model ; Test Model",9,-0.5,8.5,9,-0.5,8.5);
    h_results[6] = new TH2F("AFB_and_R_L","AFB_and_R_L_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[7] = new TH2F("AFB_and_R_R","AFB_and_R_R_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);
    h_results[8] = new TH2F("AFB_and_R_L_R","AFB_and_R_L+R_ILC; Tested Model ; Hypothesis Model; N_{#sigma}>Z",9,-0.5,8.5,9,-0.5,8.5);

    for(int i=0; i<9; i++) {
     for(int j=0; j<theory.size(); j++) {
        h_results[i]->GetXaxis()->SetBinLabel(j+1,allmodels.at(j));
        h_results[i]->GetYaxis()->SetBinLabel(j+1,allmodels.at(j));
     }
    }  

    //auxiliary histograms, to record chi2 and the ndf2. It could be done also with matrixes or two dimensional float chi2[cases=9][models=9][models=9]
    TH2F * h_chi2[9];
    TH2F * h_ndf[9];

    for(int i=0; i<9; i++) {
        h_chi2[i]=new TH2F(TString::Format("h_chi2_%i",i),TString::Format("h_chi2_%i",i),9,-0.5,8.5,9,-0.5,8.5);
        h_ndf[i]=new TH2F(TString::Format("h_ndf_%i",i),TString::Format("h_ndf_%i",i),9,-0.5,8.5,9,-0.5,8.5);
    }


    float cross_SM_L;
    float cross_SM_R;
    float cross_SM_unpol;

    for(int iobs=0; iobs<observables.size(); iobs++) {//loop for reference model
        //this could be done with iterators, but this is good enough...

        if(observables.at(iobs).flav!=quark_under_test && quark_under_test>0) continue;
        if(observables.at(iobs).energy!=energy_under_test && energy_under_test>0) continue;

        int id1= observables.at(iobs).id;

        if(id1==0) {
            cross_SM_L=observables.at(iobs).Cross_L;
            cross_SM_R=observables.at(iobs).Cross_R;
            cross_SM_unpol=observables.at(iobs).Cross_unpol;
        }
   
        for(int jobs=0; jobs<observables.size(); jobs++) {

            int id2=observables.at(jobs).id;
            //if(id2>id1) continue;

            float chi2_AFB_L=0;
            float chi2_AFB_R=0;
            float chi2_AFB_unpol=0;

            float chi2_R_L=0;
            float chi2_R_R=0;
            float chi2_R_unpol=0;

            if(observables.at(jobs).flav!=quark_under_test && quark_under_test>0) continue;
            if(observables.at(jobs).energy!=energy_under_test && energy_under_test>0) continue;
            if(quark_under_test==-1 && observables.at(jobs).flav!=observables.at(iobs).flav) continue;
            if(energy_under_test==-1 && observables.at(jobs).energy!=observables.at(iobs).energy) continue;

            //read the uncertainties from experimentalUnc.h
            //float r_stat[2][4]={{0.18,0.27,0.12,0.23},{}};//{{250 c-left, c-right, b-left, b-right},{500, same }}
            int iquark=observables.at(iobs).flav;
            if(iquark==5) iquark=6;//index to read properly the uncertainties
            int index_energy=0;
            if(observables.at(iobs).energy==500) index_energy=1;

            float r_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*r_stat[index_energy][(iquark-4)]/100.;
            float r_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*r_stat[index_energy][(iquark-4)+1]/100.;
            float afb_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_stat[index_energy][(iquark-4)]/100.;
            float afb_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_stat[index_energy][(iquark-4)+1]/100.;
            float r_syst_L=r_syst[index_energy][(iquark-4)]/100.;
            float r_syst_R=r_syst[index_energy][(iquark-4)+1]/100.;
            float afb_syst_L=afb_syst[index_energy][(iquark-4)]/100.;
            float afb_syst_R=afb_syst[index_energy][(iquark-4)+1]/100.;

            float r_total_L=sqrt(r_stat_L*r_stat_L+r_syst_L*r_syst_L);
            float r_total_R=sqrt(r_stat_R*r_stat_R+r_syst_R*r_syst_R);
            float afb_total_L=sqrt(afb_stat_L*afb_stat_L+afb_syst_L*afb_syst_L);
            float afb_total_R=sqrt(afb_stat_R*afb_stat_R+afb_syst_R*afb_syst_R);

            // In reality, the statistical uncertainty should be changed for every model, we have calculated it for the SM, 
            // for a model X, sigma_stat= sigma_stat_SM * sqrt(cross_SM/cross_X)            
            chi2_AFB_L+=pow(observables.at(jobs).AFB_L-observables.at(iobs).AFB_L,2) / pow(observables.at(iobs).AFB_L*afb_total_L,2);
            chi2_AFB_R+=pow(observables.at(jobs).AFB_R-observables.at(iobs).AFB_R,2) / pow(observables.at(iobs).AFB_R*afb_total_R,2);

            //TRandom* r1= new TRandom1();
            //chi2_AFB_R=0;
           // for(int ir=0; ir<rr; ir++) {
           //     chi2_AFB_R+=pow(observables.at(iobs).AFB_R/rr-r1->Gaus(observables.at(jobs).AFB_R,observables.at(jobs).AFB_R*afb_total_R),2) / (observables.at(iobs).AFB_R/rr);//r1->Gaus(observables.at(jobs).AFB_R,observables.at(jobs).AFB_R*afb_total_R)/rr;//pow(observables.at(iobs).AFB_R*afb_total_R,2);
           // }
            cout<<observables.at(iobs).id<<" "<<observables.at(jobs).id<<observables.at(iobs).energy<<" "<<observables.at(jobs).energy<<" "<<observables.at(iobs).AFB_R<<" "<<observables.at(jobs).AFB_R<<" "<<chi2_AFB_R<<" "<<pow(observables.at(iobs).AFB_R-observables.at(jobs).AFB_R,2) /pow(observables.at(jobs).AFB_R*afb_total_R,2)<<" "<<pow(observables.at(iobs).AFB_R-observables.at(jobs).AFB_R,1) /pow(observables.at(jobs).AFB_R*afb_total_R,1)<<endl;
            chi2_AFB_unpol+=pow(observables.at(jobs).AFB_unpol-observables.at(iobs).AFB_unpol,2) / pow(observables.at(iobs).AFB_unpol*afb_total_L,2);

            chi2_R_L+=pow(observables.at(iobs).R_L-observables.at(jobs).R_L,2) / pow(observables.at(iobs).R_L*r_total_L,2);
            chi2_R_R+=pow(observables.at(iobs).R_R-observables.at(jobs).R_R,2) / pow(observables.at(iobs).R_R*r_total_R,2);
            chi2_R_unpol+=pow(observables.at(iobs).R_unpol-observables.at(jobs).R_unpol,2) / pow(observables.at(iobs).R_unpol*r_total_L,2);

            h_chi2[0]->SetBinContent(id1+1,id2+1,chi2_AFB_L);
            h_chi2[1]->SetBinContent(id1+1,id2+1,chi2_AFB_R);
            h_chi2[2]->SetBinContent(id1+1,id2+1,chi2_AFB_L+chi2_AFB_R);

            h_chi2[3]->SetBinContent(id1+1,id2+1,chi2_R_L);
            h_chi2[4]->SetBinContent(id1+1,id2+1,chi2_R_R);
            h_chi2[5]->SetBinContent(id1+1,id2+1,chi2_R_L+chi2_R_R);

            h_chi2[6]->SetBinContent(id1+1,id2+1,chi2_AFB_L+chi2_R_L);
            h_chi2[7]->SetBinContent(id1+1,id2+1,chi2_AFB_R+chi2_R_R);
            h_chi2[8]->SetBinContent(id1+1,id2+1,chi2_AFB_L+chi2_AFB_R+chi2_R_L+chi2_R_R);

        }
    }



    float degrees=1;
    if(quark_under_test==-1) degrees*=2;
    if(energy_under_test==-1) degrees*=2;

    for(int j=0; j<9; j++) {
        for(int k=0; k<9; k++){

            h_ndf[0]->Fill(j,k,degrees);
            h_ndf[1]->Fill(j,k,degrees);
            h_ndf[2]->Fill(j,k,2.*degrees);

            h_ndf[3]->Fill(j,k,degrees);
            h_ndf[4]->Fill(j,k,degrees);
            h_ndf[5]->Fill(j,k,2.*degrees);

            h_ndf[6]->Fill(j,k,2.*degrees);
            h_ndf[7]->Fill(j,k,2.*degrees);
            h_ndf[8]->Fill(j,k,4.*degrees);
        }
    }

    TCanvas *c1 = new TCanvas ("c1","c1",1200,800);
   // Int_t colors[] = {kRed+3, kRed+1, kMagenta+3, kBlue, kCyan+2, kTeal, kGreen}; // #colors >= #levels - 1
    //gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    gStyle->SetPalette(kIsland);
    c1->Divide(3,3);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    for(int i=0; i<9; i++) {
        c1->cd(i+1);
        //gPad->SetLogz();
        for(int j=0; j<9; j++) {
            for(int k=0; k<9; k++){
                if(j==k) continue;
                float chi2v=h_chi2[i]->GetBinContent(j+1,k+1);
                float ndfv=h_ndf[i]->GetBinContent(j+1,k+1);
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
        h_results[i]->Draw("col"    );
        //h_chi2[i]->Draw("col"    );
    }



}

