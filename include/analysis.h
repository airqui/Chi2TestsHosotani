// # Copyright 2023  Adrián Irles (IFIC)

#include "../include/models_struct.h"


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

 std::vector<observables_struct_t> create_observables() {


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

    return observables;

 }

 std::vector<TH2F*> chi2_ndf_models(int quark_under_test, float energy_under_test, std::vector<observables_struct_t> observables) {

    TH2F * h_chi2[8];

    for(int i=0; i<8; i++) {
        h_chi2[i]=new TH2F(TString::Format("h_chi2_%i_q%i_en%i",i,quark_under_test,int(energy_under_test)),TString::Format("h_chi2_%i_q%i_en%i",i,quark_under_test,int(energy_under_test)),31,-0.5,30.5,31,-0.5,30.5);
    }

    float cross_SM_L;
    float cross_SM_R;
    float cross_SM_unpol;

    for(int iobs=0; iobs<observables.size(); iobs++) {//loop for reference model
        //this could be done with iterators, but this is good enough...

        if(observables.at(iobs).flav!=quark_under_test) continue;
        if(observables.at(iobs).energy!=energy_under_test) continue;

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

            if(observables.at(jobs).flav!=quark_under_test) continue;
            if(observables.at(jobs).energy!=energy_under_test ) continue;

            int iquark=observables.at(iobs).flav;
            if(iquark==5) iquark=6;//index to read properly the uncertainties
            int index_energy=0;
            if(observables.at(iobs).energy==500) index_energy=1;

            //----------
            // polarized Beams
            float r_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*r_stat[index_energy][(iquark-4)]/100.; //renormalize of the stat. uncertainties to the cross section of the model
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

            //----------
            //we assume that the unpol syst unc is the same than lef handed
            //what about the unpol stat.. uncertainty? We should assume a larger luminosity, for the moment (13thApril) 
            //let's assume a factor 4000/900 in luminosity w.r.t the L scenarios: 900fb-1 for (-0.8,0.3) 250GeV
            
            float lum_factor=4000./900.;
            // for 500GeV, the unpolarised make no sense since this energy is not reachable by FCC
            float r_stat_unpol=sqrt(1./lum_factor)*sqrt(cross_SM_L/observables.at(jobs).Cross_unpol)*r_stat_L;
            float afb_stat_unpol=sqrt(1./lum_factor)*sqrt(cross_SM_L/observables.at(jobs).Cross_unpol)*afb_stat_L;

            float r_total_unpol=sqrt(r_stat_unpol*r_stat_unpol+r_syst_L*r_syst_L);
            float afb_total_unpol=sqrt(afb_stat_unpol*afb_stat_unpol+afb_syst_L*afb_syst_L);
            // chi2 formula, using uncertainty of the tested model           
            chi2_AFB_L+=pow(observables.at(jobs).AFB_L-observables.at(iobs).AFB_L,2) / pow(observables.at(iobs).AFB_L*afb_total_L,2);
            chi2_AFB_R+=pow(observables.at(jobs).AFB_R-observables.at(iobs).AFB_R,2) / pow(observables.at(iobs).AFB_R*afb_total_R,2);
            chi2_AFB_unpol+=pow(observables.at(jobs).AFB_unpol-observables.at(iobs).AFB_unpol,2) / pow(observables.at(iobs).AFB_unpol*afb_total_unpol,2);

            chi2_R_L+=pow(observables.at(iobs).R_L-observables.at(jobs).R_L,2) / pow(observables.at(iobs).R_L*r_total_L,2);
            chi2_R_R+=pow(observables.at(iobs).R_R-observables.at(jobs).R_R,2) / pow(observables.at(iobs).R_R*r_total_R,2);
            chi2_R_unpol+=pow(observables.at(iobs).R_unpol-observables.at(jobs).R_unpol,2) / pow(observables.at(iobs).R_unpol*r_total_unpol,2);


            h_chi2[0]->SetBinContent(id1+1,id2+1,chi2_AFB_L);
            h_chi2[1]->SetBinContent(id1+1,id2+1,chi2_AFB_R);
            h_chi2[2]->SetBinContent(id1+1,id2+1,chi2_AFB_L+chi2_AFB_R);
            h_chi2[3]->SetBinContent(id1+1,id2+1,chi2_AFB_unpol);

            h_chi2[4]->SetBinContent(id1+1,id2+1,chi2_R_L);
            h_chi2[5]->SetBinContent(id1+1,id2+1,chi2_R_R);
            h_chi2[6]->SetBinContent(id1+1,id2+1,chi2_R_L+chi2_R_R);
            h_chi2[7]->SetBinContent(id1+1,id2+1,chi2_R_unpol);

        }
    }

    for(int iobs=0; iobs<observables.size(); iobs++) {//loop for reference model
        for(int jobs=0; jobs<observables.size(); jobs++) {
            for(int i=0; i<8; i++) if(h_chi2[i]->GetBinContent(jobs+1,iobs+1) > h_chi2[i]->GetBinContent(iobs+1,jobs+1))  h_chi2[i]->SetBinContent(iobs+1,jobs+1, h_chi2[i]->GetBinContent(jobs+1,iobs+1));
        }
    }

    std::vector<TH2F*> result;
    for(int i=0; i<8; i++) {
        if(i==3  || i==7) continue; // I SKIP THE UNPOL CASE !!
        result.push_back(h_chi2[i]);
    }

    return result;

 }