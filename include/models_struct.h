// example about structures
#include "TFile.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TH1.h"
#include "TROOT.h"
#include "TRint.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TString.h"
#include "TRandom.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

std::vector<TString> allmodels = {"SMB","A1","A2","A3","B1","B2","B3","B4","B5"};

struct model_values_struct_t
{
    float energy = 0;
    // index 0 is empty, 1,2,3,4,5 are for d,u,s,c,b respectively
    float cross_L[6] = {0};
    float cross_R[6] = {0};
    float AFB_L[6] = {0};
    float AFB_R[6] = {0};
};

struct model_struct_t
{
    int id = -1;
    TString modelname = "";
    model_values_struct_t model_values[10]; // up to 10 energies
};

std::vector<model_struct_t> theory;

struct observables_struct_t
{
    int id=-1;
    TString modelname = "";
    float energy = 0;
    int flav = 0;
    float Cross_L=0;
    float R_L = 0;
    float AFB_L = 0;
    float Cross_R=0;
    float R_R = 0;
    float AFB_R = 0;
    float Cross_unpol=0;
    float R_unpol = 0;
    float AFB_unpol = 0;
};

float ObsCross(int imodel, int ienergy, int iflav, float pole, float polp)
{

    float tmpcross_L = theory.at(imodel).model_values[ienergy].cross_L[iflav];
    float tmpcross_R = theory.at(imodel).model_values[ienergy].cross_R[iflav];

    float peff = (pole - polp) / (1 - pole * polp);

    float tmpcross = 0.25 * (1 - pole * polp) * ((1 - peff) * tmpcross_L + (1 + peff) * tmpcross_R);
    float result = tmpcross;
    return result;
}

float ObsR(int imodel, int ienergy, int iflav, float pole, float polp)
{

    float tmpcrossHad_L = 0;
    float tmpcrossHad_R = 0;
    for (int i = 1; i < 6; i++)
    {
        tmpcrossHad_L += theory.at(imodel).model_values[ienergy].cross_L[i];
        tmpcrossHad_R += theory.at(imodel).model_values[ienergy].cross_R[i];
    }
    float tmpcross_L = theory.at(imodel).model_values[ienergy].cross_L[iflav];
    float tmpcross_R = theory.at(imodel).model_values[ienergy].cross_R[iflav];

    float peff = (pole - polp) / (1 - pole * polp);

    float tmpcrossHad = 0.25 * (1 - pole * polp) * ((1 - peff) * tmpcrossHad_L + (1 + peff) * tmpcrossHad_R);
    float tmpcross = 0.25 * (1 - pole * polp) * ((1 - peff) * tmpcross_L + (1 + peff) * tmpcross_R);
    float result = tmpcross / tmpcrossHad;
    return result;
}

float ObsAFB(int imodel, int ienergy, int iflav, float pole, float polp)
{

    float tmpcross_L = theory.at(imodel).model_values[ienergy].cross_L[iflav];
    float tmpcross_R = theory.at(imodel).model_values[ienergy].cross_R[iflav];
    float tmpAFB_L = theory.at(imodel).model_values[ienergy].AFB_L[iflav];
    float tmpAFB_R = theory.at(imodel).model_values[ienergy].AFB_R[iflav];

    float tmpcross_plus_L = tmpcross_L * (tmpAFB_L+1) / 2;
    float tmpcross_plus_R = tmpcross_R * (tmpAFB_R+1) / 2;

    float tmpcross_minus_L = tmpcross_L - tmpcross_plus_L;
    float tmpcross_minus_R = tmpcross_R - tmpcross_plus_R;

    float peff = (pole - polp) / (1 - pole * polp);

    float plus = 0.25 * (1 - pole * polp) * ((1 - peff) * tmpcross_plus_L + (1 + peff) * tmpcross_plus_R);
    float minus = 0.25 * (1 - pole * polp) * ((1 - peff) * tmpcross_minus_L + (1 + peff) * tmpcross_minus_R);

    float result = (plus - minus) / (plus + minus);
    return result;
}

// struct theory_values_struct_t
// {
//     float cross_b;
//     float cross_c;
//     float R_b;
//     float R_c;
//     float AFB_b;
//     float AFB_c;
// }

// struct theory_struct_t
// {
//     TString model;
//     float pole;
//     float polp;
//     float E;
//     theory_values_struct theory_values;
// }

model_struct_t read_model(TString st_model = "A1", int index = 0, bool debug = true)
{

    TString filename = "../models/DataModels310823_typoCorrection/Model_" + st_model + ".txt";
    std::ifstream reading_file(filename);
    if (!reading_file)
    {
        cout << " FILE DOES NOT EXIST" << endl;
    }
    model_struct_t newmodel;
    newmodel.id = index;
    newmodel.modelname = st_model;
    float firstenergy = 0;
    int index_energy = -1;
    while (reading_file)
    {
        int flav = 0;
        Float_t tmp_energy = 0, tmp_cross_L = 0, tmp_AFB_L = 0, tmp_cross_R = 0, tmp_AFB_R;
        reading_file >> tmp_energy >> flav >> tmp_cross_L >> tmp_AFB_L >> tmp_cross_R >> tmp_AFB_R;
        if (tmp_energy != firstenergy)
        {
            index_energy++;
            firstenergy = tmp_energy;
            newmodel.model_values[index_energy].energy = tmp_energy;
        }

        newmodel.model_values[index_energy].cross_L[flav] = tmp_cross_L;
        newmodel.model_values[index_energy].cross_R[flav] = tmp_cross_R;
        newmodel.model_values[index_energy].AFB_L[flav] = tmp_AFB_L;
        newmodel.model_values[index_energy].AFB_R[flav] = tmp_AFB_R;
    }

    if (debug == true)
    {

        cout << newmodel.modelname << " index:" << newmodel.id << endl;
        for (int ie = 0; ie < 10; ie++)
            for (int ifl = 1; ifl < 6; ifl++)
                cout << newmodel.model_values[ie].energy << " " << newmodel.model_values[ie].cross_L[ifl] << " " << newmodel.model_values[ie].cross_R[ifl]
                     << " " << newmodel.model_values[ie].AFB_L[ifl] << " " << newmodel.model_values[ie].AFB_R[ifl] << endl;
    }
    return newmodel;
}
