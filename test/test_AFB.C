// # Copyright 2023  Adrián Irles (IFIC)

#include "../include/experimentalUnc.h"
#include "../include/analysis.h"
#include "../style/Style.C"
#include "../style/Labels.C"


std::vector<TH2F*> DeltaAFB(TH2F* AFBc_ILC, TH2F* AFBb_ILC,TH2F* AFBc_CLIC, TH2F* AFBb_CLIC) {

    std::vector<TString> models = {"#font[52]{A_{1}}","#font[52]{A_{2}}","#font[52]{A_{3}}","#font[52]{B}","#font[52]{B_{L}}","#font[52]{B_{H}}","#font[52]{B_{#plus}}","#font[52]{B_{#minus}}"};
    std::vector<TString> quarks = {"#color[2]{c#bar{c}}, unpol.",
    "#color[2]{c#bar{c}}, Pol:(#minus0.8, 0.0)","#color[2]{c#bar{c}}, Pol:(#minus0.8,#plus0.3)","#color[2]{c#bar{c}}, Pol:(#plus0.8, 0.0)","#color[2]{c#bar{c}}, Pol:(#plus0.8,#minus0.3)",
    "#color[4]{b#bar{b}}, unpol.",
    "#color[4]{b#bar{b}}, Pol:(#minus0.8, 0.0)","#color[4]{b#bar{b}}, Pol:(#minus0.8,#plus0.3)","#color[4]{b#bar{b}}, Pol:(#plus0.8, 0.0)","#color[4]{b#bar{b}}, Pol:(#plus0.8,#minus0.3)"};

    std::vector<float> afb;
    afb.push_back(AFBc_ILC->GetBinContent(1,3));
    afb.push_back(AFBc_CLIC->GetBinContent(1,1));
    afb.push_back(AFBc_ILC->GetBinContent(1,1));
    afb.push_back(AFBc_CLIC->GetBinContent(1,2));
    afb.push_back(AFBc_ILC->GetBinContent(1,2));
    afb.push_back(AFBb_ILC->GetBinContent(1,3));
    afb.push_back(AFBb_CLIC->GetBinContent(1,1));
    afb.push_back(AFBb_ILC->GetBinContent(1,1));
    afb.push_back(AFBb_CLIC->GetBinContent(1,2));
    afb.push_back(AFBb_ILC->GetBinContent(1,2));



    TH2F *h_AFB=  new TH2F("h_AFB","h_AFB", 8, -0.5, 7.5, 10, -0.5, 9.5);
    TH2F *h_AFB2=  new TH2F("h_AFB2","h_AFB2", 8, -0.5, 7.5, 10, -0.5, 9.5);

    for(int i=0; i<models.size(); i++) {


        //j<2
        if(i<3) {
	  h_AFB->SetBinContent(i+1,1,fabs(AFBc_ILC->GetBinContent(i+3,3)-AFBc_ILC->GetBinContent(1,3)));//cunpol                                                                                    
	  h_AFB->SetBinContent(i+1,2,fabs(AFBc_CLIC->GetBinContent(i+3,1)-AFBc_CLIC->GetBinContent(1,1)));//cLCLIC                                                                                    
	  h_AFB->SetBinContent(i+1,3,fabs(AFBc_ILC->GetBinContent(i+3,1)-AFBc_ILC->GetBinContent(1,1)));//cL                                                                                          
	  h_AFB->SetBinContent(i+1,4,fabs(AFBc_CLIC->GetBinContent(i+3,2)-AFBc_CLIC->GetBinContent(1,2)));//cRCLIC                                                                                    
            h_AFB->SetBinContent(i+1,5,fabs(AFBc_ILC->GetBinContent(i+3,2)-AFBc_ILC->GetBinContent(1,2)));//cR                                                                                          
                                                                                                                                                                                                           
            h_AFB->SetBinContent(i+1,6,fabs(AFBb_ILC->GetBinContent(i+3,3)-AFBb_ILC->GetBinContent(1,3)));//cunpol                                                                                      
            h_AFB->SetBinContent(i+1,7,fabs(AFBb_CLIC->GetBinContent(i+3,1)-AFBb_CLIC->GetBinContent(1,1)));//cLCLIC                                                                                    
            h_AFB->SetBinContent(i+1,8,fabs(AFBb_ILC->GetBinContent(i+3,1)-AFBb_ILC->GetBinContent(1,1)));//cL                                                                                          
            h_AFB->SetBinContent(i+1,9,fabs(AFBb_CLIC->GetBinContent(i+3,2)-AFBb_CLIC->GetBinContent(1,2)));//cRCLIC                                                                                    
            h_AFB->SetBinContent(i+1,10,fabs(AFBb_ILC->GetBinContent(i+3,2)-AFBb_ILC->GetBinContent(1,2)));//cR  

        } else {
	  h_AFB->SetBinContent(i+1,1,fabs(AFBc_ILC->GetBinContent(i+3,3)-AFBc_ILC->GetBinContent(2,3)));//cunpol                                                                                      
	  h_AFB->SetBinContent(i+1,2,fabs(AFBc_CLIC->GetBinContent(i+3,1)-AFBc_CLIC->GetBinContent(2,1)));//cLCLIC                                                                                    
	  h_AFB->SetBinContent(i+1,3,fabs(AFBc_ILC->GetBinContent(i+3,1)-AFBc_ILC->GetBinContent(2,1)));//cL                                                                                          
	  h_AFB->SetBinContent(i+1,4,fabs(AFBc_CLIC->GetBinContent(i+3,2)-AFBc_CLIC->GetBinContent(2,2)));//cRCLIC                                                                                    
	  h_AFB->SetBinContent(i+1,5,fabs(AFBc_ILC->GetBinContent(i+3,2)-AFBc_ILC->GetBinContent(2,2)));//cR                                                                                          
	  
	  h_AFB->SetBinContent(i+1,6,fabs(AFBb_ILC->GetBinContent(i+3,3)-AFBb_ILC->GetBinContent(2,3)));//cunpol                                                                                      
	  h_AFB->SetBinContent(i+1,7,fabs(AFBb_CLIC->GetBinContent(i+3,1)-AFBb_CLIC->GetBinContent(2,1)));//cLCLIC                                                                                    
	  h_AFB->SetBinContent(i+1,8,fabs(AFBb_ILC->GetBinContent(i+3,1)-AFBb_ILC->GetBinContent(2,1)));//cL                                                                                          
	  h_AFB->SetBinContent(i+1,9,fabs(AFBb_CLIC->GetBinContent(i+3,2)-AFBb_CLIC->GetBinContent(2,2)));//cRCLIC                                                                                    
	  h_AFB->SetBinContent(i+1,10,fabs(AFBb_ILC->GetBinContent(i+3,2)-AFBb_ILC->GetBinContent(2,2)));//cR      
        }

    h_AFB->GetXaxis()->SetBinLabel(i+1,models.at(i));
    h_AFB->GetXaxis()->SetLabelOffset();


    }

    for(int i=0; i<quarks.size(); i++) {
      if(i<5) {
	if(afb.at(i)>0) h_AFB->GetYaxis()->SetBinLabel(i+1,TString::Format("#font[52]{%s} | A^{c}_{FB}= #plus%.3f",quarks.at(i).Data(),afb.at(i)));
      if(afb.at(i)<0) h_AFB->GetYaxis()->SetBinLabel(i+1,TString::Format("#font[52]{%s} | A^{c}_{FB}= #minus%.3f",quarks.at(i).Data(),fabs(afb.at(i))));
      } else {
	if(afb.at(i)>0) h_AFB->GetYaxis()->SetBinLabel(i+1,TString::Format("#font[52]{%s} | A^{b}_{FB}= #plus%.3f",quarks.at(i).Data(),afb.at(i)));
	if(afb.at(i)<0) h_AFB->GetYaxis()->SetBinLabel(i+1,TString::Format("#font[52]{%s} | A^{b}_{FB}= #minus%.3f",quarks.at(i).Data(),fabs(afb.at(i))));
      }

    }
    h_AFB->GetZaxis()->SetTitle("|A_{FB}^{X}-A_{FB}^{SM}|");//#Delta A_{FB}^{x} [%]");//frac{A_{FB}^{X}}{A_{FB}^{SM}}-1");
    h_AFB->GetZaxis()->SetTitleOffset(1.1);
    h_AFB->GetZaxis()->SetTitleSize(0.052);
    h_AFB->GetZaxis()->SetLabelOffset(0.001);

    h_AFB->GetZaxis()->CenterTitle(true);
    h_AFB->GetXaxis()->SetLabelSize(0.09);
    h_AFB->GetYaxis()->SetLabelSize(0.05);
    h_AFB->GetXaxis()->SetLabelOffset(0.02);

    h_AFB->GetZaxis()->SetRangeUser(0.0001,0.5);
    for(int i=0; i<models.size(); i++) { 
      for(int j=0; j<quarks.size(); j++) {
	h_AFB2->SetBinContent(i+1,j+1,h_AFB->GetBinContent(i+1,j+1));
        if(h_AFB->GetBinContent(i+1,j+1)<0.0001) {
	  h_AFB->SetBinContent(i+1,j+1,0.0001001);
	  h_AFB2->SetBinContent(i+1,j+1,0.);

	}
     }
    }


    std::vector<TH2F*> result;
    result.push_back(h_AFB);
    result.push_back(h_AFB2);

    return result;

}




void AFB2D(float energy_under_test=250)
{


    allmodels = {"SMA","SMB","A1","A2","A3","B1","B2","B3","B4","B5"};

    //Int_t colors[] = {kRed+2,kRed+4, kGreen+4,kGreen+2,kGreen+1,kGreen,kCyan+2}; // #colors >= #levels - 1
    //    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
   // Int_t colors[] = {kGray + 1, kGray + 2, kGray + 3, kYellow - 6, kSpring + 5, kSpring - 1, kSpring};                  // #colors >= #levels - 1
    //TString colors_st[] = {"<1#sigma", "1-2#sigma", "2-3#sigma", "3-4#sigma", "4-5#sigma", "5-6#sigma", "#geq 6#sigma"}; // #colors >= #levels - 1
    SetQQbarStyle();
    gStyle->SetPalette(kGreenPink);//BlackBody);//TemperatureMap);//ColorPrintableOnGrey);//BlueGreenYellow);//GreenPink);//kCool);//TemperatureMap);//kIsland
    gStyle->SetOptStat(0);    
    gStyle->SetOptTitle(1);
    gStyle->SetTitleStyle(0);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleY(0.999);
    gStyle->SetTitleX(0.53);
    gStyle->SetTitleFont(12);
    gStyle->SetMarkerSize(1.2);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.10);
    gStyle->SetPadLeftMargin(0.36);


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
    std::vector<observables_struct_t> observables_ILC = create_observables_pol(0.8, 0.3);
    std::vector<observables_struct_t> observables_CLIC = create_observables_pol(0.8, 0);

    // I calculate the chi2 matrix of tested models vs hypothesis models, using th2f
    TH2F*  AFBc_ILC=AFB(4,energy_under_test, observables_ILC);
    TH2F*  AFBb_ILC=AFB(5,energy_under_test, observables_ILC);
    TH2F*  AFBc_CLIC=AFB(4,energy_under_test, observables_CLIC);
    TH2F*  AFBb_CLIC=AFB(5,energy_under_test, observables_CLIC);
    
    std::vector<TH2F *> h_DeltaAFB=DeltaAFB(AFBc_ILC,AFBb_ILC,AFBc_CLIC,AFBb_CLIC);

    TCanvas *c1 = new TCanvas(TString::Format("c1_%.f",energy_under_test), TString::Format("c1_%.f",energy_under_test), 1100, 800);
    c1->cd();//Divide(3, 3);
    gPad->SetLogz();
    //if(energy_under_test>100) h_DeltaAFB.at(0)->SetTitle(TString::Format("e^{-}e^{+} collissions at %.f GeV",energy_under_test));
    //else h_DeltaAFB.at(0)->SetTitle("e^{-}e^{+} collissions at the Z-pole");
    h_DeltaAFB.at(0)->SetTitle("#font[62]{#cbar A^{SM}_{FB} #minus A^{X}_{FB}#cbar}");
    
    gStyle->SetPaintTextFormat("0.0e");
   
    h_DeltaAFB.at(1)->SetMarkerColor(kWhite);
    //h_DeltaAFB.at(0)->SetFillStyle(3002);
    // h_DeltaAFB.at(0)->SetFillColor(kGreen);
    h_DeltaAFB.at(0)->Draw("colz");
    //    h_DeltaAFB.at(1)->Draw("textsame");

    TH2F *h_AFB=  new TH2F("h_AFB","h_AFB", 8, -0.5, 7.5, 10, -0.5, 9.5);

    for (int i=1; i<=8; i++) {
      for (int j=1; j<=10; j++) {
	if(h_DeltaAFB.at(1)->GetBinContent(i,j)==0){
	  auto t = new TLatex(h_DeltaAFB.at(1)->GetXaxis()->GetBinCenter(i), h_DeltaAFB.at(1)->GetYaxis()->GetBinCenter(j), "<1#upoint10^{-4}");
	  t->SetTextAlign(22);
	  t->SetTextColor(kWhite);
	  t->SetTextSize(0.03);
	  t->Draw();
	}
	if(h_DeltaAFB.at(1)->GetBinContent(i,j)>=0.0001 && h_DeltaAFB.at(1)->GetBinContent(i,j)<0.001){
          auto t = new TLatex(h_DeltaAFB.at(1)->GetXaxis()->GetBinCenter(i), h_DeltaAFB.at(1)->GetYaxis()->GetBinCenter(j), TString::Format("%i#upoint10^{-4}", int(h_DeltaAFB.at(1)->GetBinContent(i,j)*10000)));
          t->SetTextAlign(22);
          t->SetTextColor(kWhite);
          t->SetTextSize(0.03);
	  t->Draw();
        }
	if(h_DeltaAFB.at(1)->GetBinContent(i,j)>=0.001 && h_DeltaAFB.at(1)->GetBinContent(i,j)<0.01){
	  auto t = new TLatex(h_DeltaAFB.at(1)->GetXaxis()->GetBinCenter(i), h_DeltaAFB.at(1)->GetYaxis()->GetBinCenter(j), TString::Format("%i#upoint10^{-3}", int(h_DeltaAFB.at(1)->GetBinContent(i,j)*1000)));
          t->SetTextAlign(22);
          t->SetTextColor(kWhite);
          t->SetTextSize(0.03);
          t->Draw();
        }
	
	if(h_DeltaAFB.at(1)->GetBinContent(i,j)>=0.01 && h_DeltaAFB.at(1)->GetBinContent(i,j)<0.1){
          auto t = new TLatex(h_DeltaAFB.at(1)->GetXaxis()->GetBinCenter(i), h_DeltaAFB.at(1)->GetYaxis()->GetBinCenter(j), TString::Format("%i#upoint10^{-2}", int(h_DeltaAFB.at(1)->GetBinContent(i,j)*100)));
          t->SetTextAlign(22);
          t->SetTextColor(kWhite);
          t->SetTextSize(0.03);
          t->Draw();
        }
	if(h_DeltaAFB.at(1)->GetBinContent(i,j)>=0.1 && h_DeltaAFB.at(1)->GetBinContent(i,j)<1){
	  auto t = new TLatex(h_DeltaAFB.at(1)->GetXaxis()->GetBinCenter(i), h_DeltaAFB.at(1)->GetYaxis()->GetBinCenter(j), TString::Format("%i#upoint10^{-1}", int(h_DeltaAFB.at(1)->GetBinContent(i,j)*10)));
          t->SetTextAlign(22);
          t->SetTextColor(kWhite);
          t->SetTextSize(0.03);
          t->Draw();
        }



      }
    }
    
    //    h_DeltaAFB.at(1)->Draw("textsame");
    //AFBc_ILC->SetFillStyle(3002);
    //AFBc_ILC->Draw("scat");
    TArrow *arr0 = new TArrow(-0.5, 4.5, 7.5, 4.5, 0.02, "");
    arr0->SetLineColor(kWhite);
    arr0->SetFillColor(kWhite);
    arr0->Draw();

    TArrow *arr1 = new TArrow(-5, 4.5, -0.5, 4.5, 0.02, "");
    arr1->SetLineColor(kBlack);
    arr1->SetFillColor(kBlack);
    arr1->Draw();

    //AFBc_ILC250->Draw("colz");

    if(energy_under_test>100)  QQBARLabel2(0.05,0.93,TString::Format("#sqrt{s}_{e^{-}e^{+}} =  %.f GeV",energy_under_test), kBlack,0.045,0,62);
    else QQBARLabel2(0.05,0.93, "#sqrt{s}_{e^{-}e^{+}} =  m_{Z}", kBlack, 0.045,0,62);

    c1->Print(TString::Format("SensTheory_%.fGeV.eps",energy_under_test));
}


void test_AFB()
{

  AFB2D(91.2);
  AFB2D(250);
	  
 // AFB2D(350);
    AFB2D(380);
    AFB2D(500);
    AFB2D(1000);
    

    // test_two_quarks_two_energies();
}
