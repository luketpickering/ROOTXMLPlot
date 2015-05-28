#include <utility>
#include <iostream>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVectorD.h"

#include "TRandom.h"

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"
#include "PlottingIO.hxx"
#include "PlottingSelections.hxx"

using namespace PlottingDefaults;
using namespace PlottingUtils;
using namespace DataSpecifics;
using namespace PlottingSelections;

namespace {
  TFile* OutputTCanviFile;
  std::string OutputPDFFileName;
  std::string OutputTCanvasFileName;
  std::string HistogramCacheFileName;
}

TH1* DummyRatioGuide(int color=kBlack,
 float xMin=PlottingSelections::DPhiTMin,
 float xMax=PlottingSelections::DPhiTMax){
  TH1* dummy = new TH1F("dummy","",1,xMin, xMax);
  dummy->Fill(1);
  dummy->SetLineColor(color);
  dummy->SetLineStyle(1);
  dummy->SetLineWidth(2);
  dummy->SetDirectory(0);

  return dummy;
}

TCanvas* SaveAndPrint(TCanvas* canv){
  std::cout << "Saving TCanvas: " << canv->GetName() << std::endl;
  OutputTCanviFile->WriteTObject(canv);
  canv->Print(OutputPDFFileName.c_str());
  return canv;
}

void SetSeriesStyle(TH1* hist, int color, int linestyle){
  hist->SetLineWidth(2);
  hist->SetLineColor(color);
  hist->SetLineStyle(linestyle);
  hist->SetMarkerStyle(0);
}

//**********************Start Plotting Functions**************************/

namespace TargetDependence {
TCanvas* TargetDeptDPTCanvas(bool IsNucEmit=false){
  std::string ModeString = (IsNucEmit?"NE":"NoNE");
  std::string CanvString = (IsNucEmit?
    "TDNuclearEmissionCanvas":"TDNoNuclearEmissionCanvas");
  std::string TitleString = (IsNucEmit?
    "Target Dependence Nuclear Emission":
    "Target Dependence No Nuclear Emission");

  TCanvas* canv = SplitCanvasFactory(CanvString.c_str());

  TH1* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][ModeString]);
  MakePDF(NEUT_C_Hist);

  TH1* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"][ModeString]);
  MakePDF(NEUT_O_Hist);

  TH1* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][ModeString]);
  MakePDF(NEUT_Pb_Hist);

  TH1* NEUT_Fe_Hist = NoSaveClone(Data::Histos["NEUT"]["Fe"][ModeString]);
  MakePDF(NEUT_Fe_Hist);

  TH1* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][ModeString]);
  MakePDF(GENIE_C_Hist);

  TH1* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"][ModeString]);
  MakePDF(GENIE_O_Hist);

  TH1* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][ModeString]);
  MakePDF(GENIE_Pb_Hist);

  TH1* GENIE_Fe_Hist = NoSaveClone(Data::Histos["GENIE"]["Fe"][ModeString]);
  MakePDF(GENIE_Fe_Hist);

  canv->cd(1);

  SetUserRangeToMaxYVal({NEUT_C_Hist,
                          NEUT_O_Hist,
                          NEUT_Pb_Hist,
                          NEUT_Fe_Hist,
                          GENIE_C_Hist,
                          GENIE_O_Hist,
                          GENIE_Pb_Hist,
                          GENIE_Fe_Hist});

  SetUpperPaneHistoDefaults(NEUT_C_Hist, DPhiT_PDF_Title);

  NEUT_C_Hist->Draw("EHIST");
  NEUT_O_Hist->Draw("EHIST SAME");
  NEUT_Pb_Hist->Draw("EHIST SAME");
  NEUT_Fe_Hist->Draw("EHIST SAME");
  GENIE_C_Hist->Draw("EHIST SAME");
  GENIE_O_Hist->Draw("EHIST SAME");
  GENIE_Pb_Hist->Draw("EHIST SAME");
  GENIE_Fe_Hist->Draw("EHIST SAME");

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(NEUT_C_Hist, ("C NEUT ("+GetIntegralAsString(NEUT_C_Hist)
      +")").c_str(), "l");
  leg->AddEntry(NEUT_O_Hist, ("O NEUT ("+GetIntegralAsString(NEUT_O_Hist)
      +")").c_str(), "l");
  leg->AddEntry(NEUT_Pb_Hist, ("Pb NEUT ("+GetIntegralAsString(NEUT_Pb_Hist)
      +")").c_str(), "l");
  leg->AddEntry(NEUT_Fe_Hist, ("Fe NEUT ("+GetIntegralAsString(NEUT_Fe_Hist)
      +")").c_str(), "l");
  leg->AddEntry(GENIE_C_Hist, ("C GENIE ("+GetIntegralAsString(GENIE_C_Hist)
      +")").c_str(), "l");
  leg->AddEntry(GENIE_O_Hist, ("O GENIE ("+GetIntegralAsString(GENIE_O_Hist)
      +")").c_str(), "l");
  leg->AddEntry(GENIE_Pb_Hist, ("Pb GENIE ("+GetIntegralAsString(GENIE_Pb_Hist)
      +")").c_str(), "l");
  leg->AddEntry(GENIE_Fe_Hist, ("Fe GENIE ("+GetIntegralAsString(GENIE_Fe_Hist)
      +")").c_str(), "l");

  leg->Draw();

  canv->cd(2);

  TH1* DRG = DummyRatioGuide(PlottingSelections::DPhiTMin_deg,
                              PlottingSelections::DPhiTMax_deg);

  TH1* PB_C_Hist_Ratio_NEUT = ScaleToRef(NEUT_Pb_Hist,NEUT_C_Hist).first;
  SetLowerPaneHistoDefaults(PB_C_Hist_Ratio_NEUT,"Ratio (Scaled To Denominator)",
   DPhiTTitle);
  TH1* O_C_Hist_Ratio_NEUT = ScaleToRef(NEUT_O_Hist, NEUT_C_Hist).first;
  TH1* Fe_C_Hist_Ratio_NEUT = ScaleToRef(NEUT_Fe_Hist,NEUT_C_Hist).first;

  TH1* PB_C_Hist_Ratio_GENIE = ScaleToRef(GENIE_Pb_Hist, GENIE_C_Hist).first;
  TH1* O_C_Hist_Ratio_GENIE = ScaleToRef(GENIE_O_Hist, GENIE_C_Hist).first;
  TH1* Fe_C_Hist_Ratio_GENIE = ScaleToRef(GENIE_Fe_Hist, GENIE_C_Hist).first;



 SetUserRangeToMaxYVal({PB_C_Hist_Ratio_NEUT,
                        O_C_Hist_Ratio_NEUT,
                        Fe_C_Hist_Ratio_NEUT,
                        PB_C_Hist_Ratio_GENIE,
                        O_C_Hist_Ratio_GENIE,
                        Fe_C_Hist_Ratio_GENIE});

  PB_C_Hist_Ratio_NEUT->Draw("EHIST");
  O_C_Hist_Ratio_NEUT->Draw("EHIST SAME");
  Fe_C_Hist_Ratio_NEUT->Draw("EHIST SAME");

  PB_C_Hist_Ratio_GENIE->Draw("EHIST SAME");
  O_C_Hist_Ratio_GENIE->Draw("EHIST SAME");
  Fe_C_Hist_Ratio_GENIE->Draw("EHIST SAME");

  DRG->Draw("HIST SAME");

  leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(PB_C_Hist_Ratio_NEUT, "Pb/C (NEUT)","l");
  leg->AddEntry(PB_C_Hist_Ratio_GENIE, "Pb/C (GENIE)","l");
  leg->AddEntry(Fe_C_Hist_Ratio_GENIE, "Fe/C (GENIE)","l");
  leg->AddEntry(O_C_Hist_Ratio_NEUT, "O/C (NEUT)","l");
  leg->AddEntry(O_C_Hist_Ratio_GENIE, "O/C (GENIE)","l");
  leg->AddEntry(Fe_C_Hist_Ratio_GENIE, "Fe/C (GENIE)","l");

  leg->Draw();

  canv->cd(0);

  TitleFactory(TitleString.c_str())->Draw();

  return canv;
}

TCanvas* PionProduction(std::string Generator, std::string Target,
  std::string STVar="", bool IncludeCCPiProdInclusive=false){

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar==""){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(STVar=="Pt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  }


  TCanvas* canv = CanvasFactory(Generator+Target+"PionProduction");

  TH1* CohPiProd = NoSaveClone(Data::Histos[Generator][Target]["CohPiProd"+STVar]);
  TH1* NonCoh = NoSaveClone(Data::Histos[Generator][Target]["NonCoh"+STVar]);
  TH1* CCQEAll = NoSaveClone(Data::Histos[Generator][Target]["CCQE"+STVar]);

  AddCountToTHTitle({std::make_pair("CCCoh #pi^{+} Production",CohPiProd),
                      std::make_pair("CC non-Coh #pi^{+} Production",NonCoh),
                      std::make_pair("(CCQE Proton MisId)",CCQEAll)});

  CohPiProd->SetLineColor(kPink);
  NonCoh->SetLineColor(kRed-2);
  CCQEAll->SetLineColor(kBlue-2);

  CohPiProd->SetLineStyle(1);
  NonCoh->SetLineStyle(1);
  CCQEAll->SetLineStyle(1);

  if(IncludeCCPiProdInclusive){
    TH1* CCPiProdInc = NoSaveClone(Data::Histos[Generator][Target]["CCPiProdInc"+STVar]);
    std::initializer_list<TH1*> hists = {CohPiProd,
                                          NonCoh,
                                          CCQEAll,
                                          CCPiProdInc};
    CCPiProdInc->SetLineColor(kOrange-2);
    CCPiProdInc->SetLineStyle(1);

    AddCountToTHTitle({std::make_pair("CC #pi^{+} Production",CCPiProdInc)});
    MakePDF(hists);
    SetUserRangeToMaxYVal(hists);
    SetOnePaneHistDefaults(CohPiProd, YTitleText.c_str(), XTitleText.c_str());
    DrawTHs(hists);

    TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*3.0 - Pad0XMargin),
      (1.0 - P0LegendHeight*0.5 - Pad0YMargin),
      (1.0 - Pad0XMargin),
      (1.0 - Pad0YMargin));

    AddToLegend(hists,leg);
    leg->Draw();
  } else {
    std::initializer_list<TH1*> hists = {CohPiProd,
                                          NonCoh,
                                          CCQEAll};
    MakePDF(hists);
    SetUserRangeToMaxYVal(hists);
    SetOnePaneHistDefaults(CohPiProd, YTitleText.c_str(), XTitleText.c_str());
    DrawTHs(hists);

    TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*3.0 - Pad0XMargin),
      (1.0 - P0LegendHeight*0.5 - Pad0YMargin),
      (1.0 - Pad0XMargin),
      (1.0 - Pad0YMargin));

    AddToLegend(hists,leg);
    leg->Draw();
  }

  TitleFactory(Generator + " Target: " + Target
    + " Coherent Pion Production " + XTitleText)->Draw();

  return canv;
}

TCanvas* CCQEAllDistributions(){
  TCanvas* canv = CanvasFactory("CCQEAllCanvas",
    false);

  TH1* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["Inc"]);
  MakePDF(NEUT_C_Hist);

  // TH1* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"]["Inc"]);
  // MakePDF(NEUT_O_Hist);

  // TH1* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"]["Inc"]);
  // MakePDF(NEUT_Pb_Hist);

  TH1* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"]["Inc"]);
  MakePDF(GENIE_C_Hist);

  // TH1* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"]["Inc"]);
  // MakePDF(GENIE_O_Hist);

  // TH1* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"]["Inc"]);
  // MakePDF(GENIE_Pb_Hist);

 SetUserRangeToMaxYVal({NEUT_C_Hist,
                        // NEUT_O_Hist,
                        // NEUT_Pb_Hist,
                        GENIE_C_Hist});
                        // GENIE_O_Hist,
                        // GENIE_Pb_Hist});

NEUT_C_Hist->SetLineWidth(2);
NEUT_C_Hist->SetLineStyle(1);
// NEUT_O_Hist->SetLineWidth(2);
// NEUT_O_Hist->SetLineStyle(1);
// NEUT_Pb_Hist->SetLineWidth(2);
// NEUT_Pb_Hist->SetLineStyle(1);
GENIE_C_Hist->SetLineWidth(2);
GENIE_C_Hist->SetLineStyle(2);
// GENIE_O_Hist->SetLineWidth(2);
// GENIE_O_Hist->SetLineStyle(2);
// GENIE_Pb_Hist->SetLineWidth(2);
// GENIE_Pb_Hist->SetLineStyle(2);

  NEUT_C_Hist->Draw("EHIST");
  // NEUT_O_Hist->Draw("EHIST SAME");
  // NEUT_Pb_Hist->Draw("EHIST SAME");
  GENIE_C_Hist->Draw("EHIST SAME");
  // GENIE_O_Hist->Draw("EHIST SAME");
  // GENIE_Pb_Hist->Draw("EHIST SAME");

  TitleFactory("bladeblablablabladebla")->Draw();
  return canv;

}

TCanvas* SpeciesCanvas(char const * species){

  std::string CanvString = species;

  TCanvas* canv = CanvasFactory(CanvString.c_str(), false);

  TH1* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][species]);
  MakePDF(NEUT_C_Hist);

  TH1* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"][species]);
  MakePDF(NEUT_O_Hist);

  TH1* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][species]);
  MakePDF(NEUT_Pb_Hist);

  TH1* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][species]);
  MakePDF(GENIE_C_Hist);

  TH1* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"][species]);
  MakePDF(GENIE_O_Hist);

  TH1* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][species]);
  MakePDF(GENIE_Pb_Hist);

  // NEUT_C_Hist->GetYaxis()->SetRangeUser(1E-2,1.05);

  NEUT_C_Hist->GetXaxis()->SetTitle(species);
  NEUT_C_Hist->GetYaxis()->SetTitle(("#frac{dN}{Nd("
    +std::string(species)+")}").c_str());
  NEUT_C_Hist->GetXaxis()->SetTitleOffset(XTitleOffset);
  NEUT_C_Hist->GetXaxis()->SetLabelOffset(XLabelOffset);

  canv->SetLogy(true);

  NEUT_Pb_Hist->SetLineWidth(4);
  NEUT_O_Hist->SetLineWidth(4);
  NEUT_C_Hist->SetLineWidth(4);
  NEUT_C_Hist->SetLineWidth(4);
  GENIE_Pb_Hist->SetLineWidth(4);
  GENIE_O_Hist->SetLineWidth(4);
  GENIE_C_Hist->SetLineWidth(4);
  GENIE_C_Hist->SetLineWidth(4);

  NEUT_C_Hist->Draw("EHIST");
  NEUT_O_Hist->Draw("EHIST SAME");
  NEUT_Pb_Hist->SetLineWidth(2);
  NEUT_Pb_Hist->Draw("EHIST SAME");
  GENIE_C_Hist->Draw("EHIST SAME");
  GENIE_O_Hist->Draw("EHIST SAME");
  GENIE_Pb_Hist->SetLineWidth(2);
  GENIE_Pb_Hist->Draw("EHIST SAME");

  TitleFactory((std::string(species) + " CCQE").c_str())->Draw();

  TLegend* leg = LegendFactory(
    (1.0 - P0LegendWidth - Pad0XMargin),
    (1.0 - (P0LegendHeight*1.25) - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin)
    );

  leg->AddEntry(NEUT_C_Hist,"C (NEUT)","l");
  leg->AddEntry(NEUT_O_Hist,"O (NEUT)","l");
  leg->AddEntry(NEUT_Pb_Hist,"Pb (NEUT)","l");
  leg->AddEntry(GENIE_C_Hist,"C (GENIE)","l");
  leg->AddEntry(GENIE_O_Hist,"O (GENIE)","l");
  leg->AddEntry(GENIE_Pb_Hist,"Pb (GENIE)","l");

  leg->Draw();

  return canv;
}

}

namespace GeneratorDependence {

TCanvas* EmissionFractions(std::string generator, bool FractionNumerIsEmit){

  std::string CanvString = generator + (FractionNumerIsEmit?
    "NuclearEmission":"NoNuclearEmission");

  TCanvas* canv = CanvasFactory((generator+CanvString+"FractionsCanv").c_str());

  int linestyle = 1;
  (void)linestyle;

//Carbon************************************************************************
  TH1* C_Numer = NoSaveClone(FractionNumerIsEmit ?
      Data::Histos[generator]["C"]["NE"] : Data::Histos[generator]["C"]["NoNE"]);
  TH1* C_Denom =  NoSaveClone(Data::Histos[generator]["C"]["Inc"]);

  TH1* C_Ratio = ScaleToRef(C_Numer,C_Denom).first;

  C_Ratio->GetYaxis()->SetRangeUser(0,1.05);
  C_Ratio->GetYaxis()->SetTitle( (FractionNumerIsEmit ?
    "Nuclear Emission Fraction" : "No Nuclear Emission Fraction") );

  SetOnePaneHistDefaults(C_Ratio);
  C_Ratio->Draw("EHIST");

//Oxygen************************************************************************

  TH1* O_Numer = NoSaveClone(FractionNumerIsEmit ?
      Data::Histos[generator]["O"]["NE"] : Data::Histos[generator]["O"]["NoNE"]);
  // BinWidthNormalise(O_Numer);
  TH1* O_Denom = NoSaveClone(Data::Histos[generator]["O"]["Inc"]);
  // BinWidthNormalise(O_Denom);

  TH1* O_Ratio = ScaleToRef(O_Numer,O_Denom).first;

  O_Ratio->GetYaxis()->SetRangeUser(0,1.05);
  O_Ratio->GetYaxis()->SetTitle( (FractionNumerIsEmit ?
    "Nuclear Emission Fraction" : "No Nuclear Emission Fraction") );

  SetOnePaneHistDefaults(O_Ratio);
  O_Ratio->Draw("EHIST SAME");

//Lead**************************************************************************

  TH1* Pb_Numer = NoSaveClone(FractionNumerIsEmit ?
      Data::Histos[generator]["Pb"]["NE"] : Data::Histos[generator]["Pb"]["NoNE"]);
  // BinWidthNormalise(Pb_Numer);
  TH1* Pb_Denom = NoSaveClone(Data::Histos[generator]["Pb"]["Inc"]);
  // BinWidthNormalise(Pb_Denom);

  TH1* Pb_Ratio = ScaleToRef(Pb_Numer,Pb_Denom).first;

  Pb_Ratio->GetYaxis()->SetRangeUser(0,1.05);
  Pb_Ratio->GetYaxis()->SetTitle( (FractionNumerIsEmit ?
    "Nuclear Emission Fraction" : "No Nuclear Emission Fraction") );

  SetOnePaneHistDefaults(Pb_Ratio);
  Pb_Ratio->Draw("EHIST SAME");

  TitleFactory( (std::string((FractionNumerIsEmit ?
      "Nuclear Emission Fraction" : "No Nuclear Emission Fraction"))
      + " CCQE Fraction " + generator).c_str())->Draw();

  TLegend* leg = LegendFactory(
    (1.0 - P0LegendWidth - Pad0XMargin),
    FractionNumerIsEmit ? (Pad0YMargin) : (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    FractionNumerIsEmit ? (Pad0YMargin + P0LegendHeight) : (1.0 - Pad0YMargin)
    );

  leg->AddEntry(C_Ratio,"Carbon","l");
  leg->AddEntry(O_Ratio,"Oxygen","l");
  leg->AddEntry(Pb_Ratio,"Lead","l");

  leg->Draw();

  return canv;
}

TCanvas* GeneratorDeptDPTCanvas(bool IsNucEmit=false){
  std::string ModeString = (IsNucEmit?"NE":"NoNE");
  std::string CanvString = (IsNucEmit?
    "GDNuclearEmissionCanvas":"GDNoNuclearEmissionCanvas");
  std::string TitleString = (IsNucEmit?
    "Generator Dependence Nuclear Emission":
    "Generator Dependence No Nuclear Emission");

  TCanvas* canv = SplitCanvasFactory(CanvString.c_str());

  TH1* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][ModeString]);
  MakePDF(NEUT_C_Hist);

  TH1* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"][ModeString]);
  MakePDF(NEUT_O_Hist);

  TH1* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][ModeString]);
  MakePDF(NEUT_Pb_Hist);

  TH1* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][ModeString]);
  MakePDF(GENIE_C_Hist);

  TH1* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"][ModeString]);
  MakePDF(GENIE_O_Hist);

  TH1* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][ModeString]);
  MakePDF(GENIE_Pb_Hist);

  canv->cd(1);

  SetUserRangeToMaxYVal({NEUT_C_Hist,
                          NEUT_O_Hist,
                          NEUT_Pb_Hist,
                          GENIE_C_Hist,
                          GENIE_O_Hist,
                          GENIE_Pb_Hist});

  SetUpperPaneHistoDefaults(NEUT_C_Hist, DPhiT_PDF_Title);

  NEUT_C_Hist->Draw("EHIST");
  NEUT_O_Hist->Draw("EHIST SAME");
  NEUT_Pb_Hist->Draw("EHIST SAME");
  GENIE_C_Hist->Draw("EHIST SAME");
  GENIE_O_Hist->Draw("EHIST SAME");
  GENIE_Pb_Hist->Draw("EHIST SAME");

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(NEUT_C_Hist, "C NEUT", "l");
  leg->AddEntry(NEUT_O_Hist, "O NEUT", "l");
  leg->AddEntry(NEUT_Pb_Hist, "Pb NEUT", "l");
  leg->AddEntry(GENIE_C_Hist, "C GENIE", "l");
  leg->AddEntry(GENIE_O_Hist, "O GENIE", "l");
  leg->AddEntry(GENIE_Pb_Hist, "Pb GENIE", "l");

  leg->Draw();

  canv->cd(2);

  TH1* Pb_GENIE_NEUT_Ratio = ScaleToRef(GENIE_Pb_Hist, NEUT_Pb_Hist).first;
  TH1* O_GENIE_NEUT_Ratio = ScaleToRef(GENIE_O_Hist, NEUT_O_Hist).first;
  TH1* C_GENIE_NEUT_Ratio = ScaleToRef(GENIE_C_Hist, NEUT_C_Hist).first;

  SetUserRangeToMaxYVal({Pb_GENIE_NEUT_Ratio,
                          O_GENIE_NEUT_Ratio,
                          C_GENIE_NEUT_Ratio});

  SetLowerPaneHistoDefaults(Pb_GENIE_NEUT_Ratio, "Ratio (Scaled To Denominator)",
   DPhiTTitle);

  SetSeriesStyle(Pb_GENIE_NEUT_Ratio,kGreen,1);
  Pb_GENIE_NEUT_Ratio->Draw("EHIST");

  SetSeriesStyle(O_GENIE_NEUT_Ratio, kRed,1);
  O_GENIE_NEUT_Ratio->Draw("EHIST SAME");

  SetSeriesStyle(C_GENIE_NEUT_Ratio, kBlue,1);
  C_GENIE_NEUT_Ratio->Draw("EHIST SAME");

  TH1* DRG = DummyRatioGuide(PlottingSelections::DPhiTMin_deg,
                              PlottingSelections::DPhiTMax_deg);
  DRG->Draw("HIST SAME");

  leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(Pb_GENIE_NEUT_Ratio, "GENIE/NEUT (Pb)","l");
  leg->AddEntry(O_GENIE_NEUT_Ratio, "GENIE/NEUT (O)","l");
  leg->AddEntry(C_GENIE_NEUT_Ratio, "GENIE/NEUT (C)","l");

  leg->Draw();

  canv->cd(0);

  TitleFactory(TitleString.c_str())->Draw();

  return canv;
}

TCanvas* ProtonMomentumCanvas(char const * Target){

  TCanvas* canv = CanvasFactory((std::string(Target)
    + "ProtonMomentumCanv").c_str(), false);

  canv->SetLogy(true);

  TH1* NEUT_T_First = NoSaveClone(Data::Histos["NEUT"][Target]["FirstPMom"]);
  MakePDF(NEUT_T_First); SetSeriesStyle(NEUT_T_First,kRed,1);

  TH1* NEUT_T_Highest = NoSaveClone(Data::Histos["NEUT"][Target]["HighestPMom"]);
  MakePDF(NEUT_T_Highest); SetSeriesStyle(NEUT_T_Highest,kBlue,1);

  TH1* GENIE_T_First = NoSaveClone(Data::Histos["GENIE"][Target]["FirstPMom"]);
  MakePDF(GENIE_T_First); SetSeriesStyle(GENIE_T_First,kRed,2);

  TH1* GENIE_T_Highest = NoSaveClone(Data::Histos["GENIE"][Target]["HighestPMom"]);
  MakePDF(GENIE_T_Highest); SetSeriesStyle(GENIE_T_Highest,kBlue,2);

  NEUT_T_First->GetXaxis()->SetTitle("Proton Momentum [MeV/C]");
  NEUT_T_First->GetXaxis()->SetTitleOffset(XTitleOffset);
  NEUT_T_First->GetXaxis()->SetLabelOffset(XLabelOffset);

  SetUserRangeToMaxYVal({NEUT_T_First,
                          NEUT_T_Highest,
                          GENIE_T_First,
                          GENIE_T_Highest}, 1E-10);

  NEUT_T_First->Draw("EHIST");
  NEUT_T_Highest->Draw("EHIST SAME");
  GENIE_T_First->Draw("EHIST SAME");
  GENIE_T_Highest->Draw("EHIST SAME");

  TLegend* leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));

  leg->AddEntry(NEUT_T_First, (std::string(Target) + " First NEUT").c_str(), "l");
  leg->AddEntry(NEUT_T_Highest, (std::string(Target) + " Highest NEUT").c_str(), "l");
  leg->AddEntry(GENIE_T_First, (std::string(Target) + " First GENIE").c_str(), "l");
  leg->AddEntry(GENIE_T_Highest, (std::string(Target) + " Highest GENIE").c_str(), "l");

  leg->Draw();

  TitleFactory((std::string(Target) + " First and Highest Proton Momenta").c_str())->Draw();

  return canv;
}

}

namespace ModeDependence {
  TCanvas* EOverNE(char const * Target, std::string generator){

  TCanvas* canv = SplitCanvasFactory(
    (Target + std::string("NEoverNoNEcanv") + generator).c_str());

  TH1* T_NoNE_Hist = NoSaveClone(Data::Histos[generator][Target]["NoNE"]);
  MakePDF(T_NoNE_Hist); SetSeriesStyle(T_NoNE_Hist,kGreen-3,1);

  TH1* T_NE_Hist = NoSaveClone(Data::Histos[generator][Target]["NE"]);
  MakePDF(T_NE_Hist); SetSeriesStyle(T_NE_Hist,kCyan,1);

  canv->cd(1);

  std::vector<double> YVals;

  YVals.push_back(GetMaximumBinContents(T_NoNE_Hist));
  YVals.push_back(GetMaximumBinContents(T_NE_Hist));

  double maxY = (*std::max_element(YVals.begin(),YVals.end()));
  YVals.clear();

  T_NoNE_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  T_NE_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);

  SetUpperPaneHistoDefaults(T_NoNE_Hist, DPhiT_PDF_Title);

  T_NoNE_Hist->Draw("EHIST");
  T_NE_Hist->Draw("EHIST SAME");

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight*0.3 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(T_NoNE_Hist, (std::string("No Emit ") + Target + " "
    + generator).c_str(), "l");
  leg->AddEntry(T_NE_Hist, (std::string("Emit ") + Target + " "
    + generator).c_str(), "l");

  leg->Draw();

  canv->cd(2);

  TH1* T_NE_NoNE_Ratio = ScaleToRef(T_NE_Hist, T_NoNE_Hist).first;
  SetSeriesStyle(T_NE_NoNE_Ratio,kGreen-3,1);

  SetLowerPaneHistoDefaults(T_NE_NoNE_Ratio, "Ratio (Scaled To Denominator)",
    DPhiTTitle);

  T_NE_NoNE_Ratio->Draw("EHIST");

  TH1* DRG = DummyRatioGuide(PlottingSelections::DPhiTMin_deg,
                              PlottingSelections::DPhiTMax_deg);
  DRG->Draw("HIST SAME");

  leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight*0.25 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(T_NE_NoNE_Ratio, (std::string("Emit/No Emit ") + Target + " "
    + generator).c_str(),"l");

  leg->Draw();

  canv->cd(0);

  TitleFactory((std::string("Emission dependence ") + Target
    + " " + generator).c_str())->Draw();

  return canv;
  }
}

// TCanvas* LinearInset(){
//   TCanvas* canv = CanvasFactory("LinearInset");

//   TH1* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"]["NoNE"]);
//   MakePDF(GENIE_O_Hist);

//   TH1* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"]["NoNE"]);
//   MakePDF(GENIE_C_Hist);

//   TH1* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"]["NoNE"]);
//   MakePDF(GENIE_Pb_Hist);

//   TH1* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"]["NoNE"]);
//   MakePDF(NEUT_O_Hist);

//   TH1* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["NoNE"]);
//   MakePDF(NEUT_C_Hist);

//   TH1* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"]["NoNE"]);
//   MakePDF(NEUT_Pb_Hist);

//   SetOnePaneHistDefaults(GENIE_Pb_Hist);
//   GENIE_Pb_Hist->GetYaxis()->SetTitleOffset(YTitleOffset*1.1);
//   GENIE_Pb_Hist->GetYaxis()->SetTitle(DPhiT_PDF_Title);

//   GENIE_Pb_Hist->Draw("EHIST");
//   GENIE_C_Hist->Draw("EHIST SAME");
//   GENIE_O_Hist->Draw("EHIST SAME");
//   NEUT_Pb_Hist->Draw("EHIST SAME");
//   NEUT_C_Hist->Draw("EHIST SAME");
//   NEUT_O_Hist->Draw("EHIST SAME");

//   TLegend* leg = LegendFactory(
//     (1.0 - P0LegendWidth - (Pad0XMargin*1.05)),
//     (Pad0YMargin*1.2),
//     (1.0 - (Pad0XMargin*1.05)),
//     ((Pad0YMargin*1.2) + (P0LegendHeight*0.7) ));

//   leg->AddEntry(GENIE_O_Hist,"O (GENIE)","l");
//   leg->AddEntry(GENIE_C_Hist,"C (GENIE)","l");
//   leg->AddEntry(GENIE_Pb_Hist,"Pb (GENIE)","l");
//   leg->AddEntry(NEUT_O_Hist,"O (NEUT)","l");
//   leg->AddEntry(NEUT_C_Hist,"C (NEUT)","l");
//   leg->AddEntry(NEUT_Pb_Hist,"Pb (NEUT)","l");

//   double InsetWidth = 0.6;
//   double InsetHeight = 0.55;
//   // leg->Draw();

//   TPad *subpad = new TPad("subpad","",
//     (1.0 - InsetWidth - (Pad0XMargin*1.1)),
//     (1.0 - InsetHeight - (Pad0YMargin*1.2)),
//     (1.0 - (Pad0XMargin*1.1)),
//     (1.0 - (Pad0YMargin*1.2)));
//   subpad->Draw();
//   subpad->cd();
//   subpad->SetRightMargin(0.05);


//   TH1* GENIE_O_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["O"]["NoNEZoom"]);
//   MakePDF(GENIE_O_Hist_Zoom);

//   TH1* GENIE_C_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["C"]["NoNEZoom"]);
//   MakePDF(GENIE_C_Hist_Zoom);

//   TH1* GENIE_Pb_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["Pb"]["NoNEZoom"]);
//   MakePDF(GENIE_Pb_Hist_Zoom);

//   TH1* NEUT_O_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["O"]["NoNEZoom"]);
//   MakePDF(NEUT_O_Hist_Zoom);

//   TH1* NEUT_C_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["C"]["NoNEZoom"]);
//   MakePDF(NEUT_C_Hist_Zoom);

//   TH1* NEUT_Pb_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["Pb"]["NoNEZoom"]);
//   MakePDF(NEUT_Pb_Hist_Zoom);

//   SetOnePaneHistDefaults(GENIE_Pb_Hist_Zoom);
//   GENIE_Pb_Hist_Zoom->GetYaxis()->SetTitleOffset(YTitleOffset*1.1);
//   GENIE_Pb_Hist_Zoom->GetYaxis()->SetTitle(DPhiT_PDF_Title);
//   GENIE_Pb_Hist_Zoom->GetYaxis()->SetRangeUser(0,
//     GetMaximumBinContents(GENIE_Pb_Hist_Zoom)*1.1);

//   GENIE_Pb_Hist_Zoom->Draw("EHIST");
//   GENIE_C_Hist_Zoom->Draw("EHIST SAME");
//   GENIE_O_Hist_Zoom->Draw("EHIST SAME");
//   NEUT_Pb_Hist_Zoom->Draw("EHIST SAME");
//   NEUT_C_Hist_Zoom->Draw("EHIST SAME");
//   NEUT_O_Hist_Zoom->Draw("EHIST SAME");

//   canv->cd();
//   leg->Draw();

//   return canv;
// }

TCanvas* TargetDependenceRatio(std::string mode){
  TCanvas* canv = SplitCanvasFactory("GaussianExampleCanvas");
  canv->cd(1);
  TH1* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][mode]);
  TH1* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][mode]);
  TH1* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][mode]);
  TH1* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][mode]);
  TH1* NuWro_C_Hist = NoSaveClone(Data::Histos["NuWro"]["C"][mode]);
  TH1* NuWro_Pb_Hist = NoSaveClone(Data::Histos["NuWro"]["Pb"][mode]);
  TH1* GiBUU_C_Hist = NoSaveClone(Data::Histos["GiBUU"]["C"][mode]);
  TH1* GiBUU_Pb_Hist = NoSaveClone(Data::Histos["GiBUU"]["Pb"][mode]);

  AddCountToTHTitle({std::make_pair("GENIE C",GENIE_C_Hist),
        std::make_pair("GENIE Pb",GENIE_Pb_Hist),
        std::make_pair("NEUT C",NEUT_C_Hist),
        std::make_pair("NEUT Pb",NEUT_Pb_Hist),
        std::make_pair("NuWro C",NuWro_C_Hist),
        std::make_pair("NuWro Pb",NuWro_Pb_Hist),
        std::make_pair("GiBUU C",GiBUU_C_Hist),
        std::make_pair("GiBUU Pb",GiBUU_Pb_Hist)});

  MakePDF(GENIE_C_Hist);
  MakePDF(GENIE_Pb_Hist);
  MakePDF(NEUT_C_Hist);
  MakePDF(NEUT_Pb_Hist);
  MakePDF(NuWro_C_Hist);
  MakePDF(NuWro_Pb_Hist);
  MakePDF(GiBUU_C_Hist);
  MakePDF(GiBUU_Pb_Hist);

  GENIE_C_Hist->SetLineColor(kBlue);
  GENIE_Pb_Hist->SetLineColor(kBlack);
  NEUT_C_Hist->SetLineColor(kBlue);
  NEUT_Pb_Hist->SetLineColor(kBlack);

  NuWro_C_Hist->SetLineStyle(1);
  NuWro_Pb_Hist->SetLineStyle(1);
  GiBUU_C_Hist->SetLineStyle(2);
  GiBUU_Pb_Hist->SetLineStyle(2);

  NuWro_C_Hist->SetLineColor(kRed);
  NuWro_Pb_Hist->SetLineColor(kGreen-3);
  GiBUU_C_Hist->SetLineColor(kRed);
  GiBUU_Pb_Hist->SetLineColor(kGreen-3);

  SetUpperPaneHistoDefaults(GENIE_Pb_Hist, DPhiT_PDF_Title, DPhiTTitle);

  SetUserRangeToMaxYVal({GENIE_Pb_Hist,
                        GENIE_C_Hist,
                        NEUT_C_Hist,
                        NEUT_Pb_Hist,
                        NuWro_C_Hist,
                        NuWro_Pb_Hist,
                        GiBUU_C_Hist,
                        GiBUU_Pb_Hist});
  DrawTHs({GENIE_Pb_Hist,
                        GENIE_C_Hist,
                        NEUT_C_Hist,
                        NEUT_Pb_Hist,
                        NuWro_C_Hist,
                        NuWro_Pb_Hist,
                        GiBUU_C_Hist,
                        GiBUU_Pb_Hist});

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth*2.3 - Pad1XMargin),
    (1.0 - P1LegendHeight*1.5 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));
  AddToLegend({GENIE_Pb_Hist,
                        GENIE_C_Hist,
                        NEUT_C_Hist,
                        NEUT_Pb_Hist,
                        NuWro_C_Hist,
                        NuWro_Pb_Hist,
                        GiBUU_C_Hist,
                        GiBUU_Pb_Hist},leg);

  leg->Draw();

  canv->cd(2);
  TH1* GENIE_Ratio = ScaleToRef(GENIE_Pb_Hist,GENIE_C_Hist).first;
  TH1* NEUT_Ratio = ScaleToRef(NEUT_Pb_Hist,NEUT_C_Hist).first;
  TH1* NuWro_Ratio = ScaleToRef(NuWro_Pb_Hist,NuWro_C_Hist).first;
  TH1* GiBUU_Ratio = ScaleToRef(GiBUU_Pb_Hist,GiBUU_C_Hist).first;

  GENIE_Ratio->SetTitle("GENIE Pb/C");
  NEUT_Ratio->SetTitle("NEUT Pb/C");
  NuWro_Ratio->SetTitle("NuWro Pb/C");
  GiBUU_Ratio->SetTitle("GiBUU Pb/C");

  TH1* DRG = DummyRatioGuide(kBlack,0,180);

  GENIE_Ratio->SetLineStyle(1);
  NEUT_Ratio->SetLineStyle(1);
  NuWro_Ratio->SetLineStyle(1);
  GiBUU_Ratio->SetLineStyle(1);

  GENIE_Ratio->SetLineColor(kRed+2);
  NEUT_Ratio->SetLineColor(kBlue-2);
  NuWro_Ratio->SetLineColor(kGreen-3);
  GiBUU_Ratio->SetLineColor(kOrange-2);

  TLegend* leg2 = LegendFactory( (Pad1XMargin),
  (1.0 - P1LegendHeight - Pad1YMargin),
  (Pad1XMargin + P1LegendWidth*1.1),
  (1.0 - Pad1YMargin));
  AddToLegend({GENIE_Ratio,NEUT_Ratio,NuWro_Ratio,GiBUU_Ratio},leg2);
  SetLowerPaneHistoDefaults(GENIE_Ratio, "Ratio (Scaled To Denominator)",
    DPhiTTitle);
  DRG->SetLineColor(kBlack);
  SetUserRangeToMaxYVal({GENIE_Ratio,NEUT_Ratio,NuWro_Ratio,GiBUU_Ratio});
  DrawTHs({GENIE_Ratio,NEUT_Ratio,NuWro_Ratio,GiBUU_Ratio});
  DRG->Draw("HIST SAME");
  leg2->Draw();
  return canv;
}

TCanvas* FSIModeCompare(std::string Generator, std::string Target,
  std::string STVar="", bool ShowMecAndOthers=false){

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar==""){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(STVar=="Pt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  }

  TCanvas* canv = CanvasFactory(Generator+Target+"FSIMODE");

  TH1* CCQEAll_norescat =
    NoSaveClone(Data::Histos[Generator][Target+"_norescat"]["CCQE"+STVar]);
  TH1* CCQEAll = NoSaveClone(Data::Histos[Generator][Target]["CCQE"+STVar]);
  TH1* CCQEAll_NE = NoSaveClone(Data::Histos[Generator][Target]["CCQENE"+STVar]);

  if(ShowMecAndOthers){
    TH1* CCNonQELike =
      NoSaveClone(Data::Histos[Generator][Target]["CC_not_QELike"+STVar]);
    TH1* MEC = NoSaveClone(Data::Histos[Generator][Target]["MEC"+STVar]);
    TH1* CCQELike = NoSaveClone(Data::Histos[Generator][Target]["CCQELike"+STVar]);


      std::initializer_list<TH1*> hists = {CCQEAll_norescat,
                    CCQEAll,
                    CCQELike,
                    CCQEAll_NE,
                    MEC,
                    CCNonQELike};

      AddCountToTHTitle({std::make_pair("CCQE (Generator HT Off)",CCQEAll_norescat),
                          std::make_pair("CCQE",CCQEAll),
                          std::make_pair("CCQE With Nuclear Emission",CCQEAll_NE),
                          std::make_pair("CCQE+MEC",CCQELike),
                          std::make_pair("MEC",MEC),
                          std::make_pair("CC non-QE Not MEC",CCNonQELike)});
      CCQEAll_NE->SetLineColor(kBlue-2);
      CCQELike->SetLineColor(kGreen-3);
      MEC->SetLineColor(kMagenta-2);
      CCNonQELike->SetLineColor(kOrange-2);
      MakePDF(hists);
      SetUserRangeToMaxYVal(hists);


      CCQEAll_norescat->SetLineStyle(1);
      CCQEAll->SetLineStyle(1);
      CCQEAll_NE->SetLineStyle(1);
      CCQELike->SetLineStyle(1);
      MEC->SetLineStyle(1);
      CCNonQELike->SetLineStyle(1);
      SetOnePaneHistDefaults(CCQEAll_norescat,YTitleText.c_str(),
      XTitleText.c_str());

      DrawTHs(hists);

      TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*3.0 - Pad0XMargin),
        (1.0 - P0LegendHeight - Pad0YMargin),
        (1.0 - Pad0XMargin),
        (1.0 - Pad0YMargin));

      AddToLegend(hists,leg);
      leg->Draw();
  } else {

    std::initializer_list<TH1*> hists = {CCQEAll_norescat,
                  CCQEAll,
                  CCQEAll_NE};

    AddCountToTHTitle({std::make_pair("CCQE (Generator HT Off)",CCQEAll_norescat),
                        std::make_pair("CCQE",CCQEAll),
                        std::make_pair("CCQE With Nuclear Emission",CCQEAll_NE)});
    CCQEAll_NE->SetLineColor(kBlue-2);
    MakePDF(hists);
    SetUserRangeToMaxYVal(hists);


    CCQEAll_norescat->SetLineStyle(1);
    CCQEAll->SetLineStyle(1);
    CCQEAll_NE->SetLineStyle(1);
    SetOnePaneHistDefaults(CCQEAll_norescat,YTitleText.c_str(),
    XTitleText.c_str());

    DrawTHs(hists);

    TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*3.0 - Pad0XMargin),
      (1.0 - P0LegendHeight*0.5 - Pad0YMargin),
      (1.0 - Pad0XMargin),
      (1.0 - Pad0YMargin));

    AddToLegend(hists,leg);
    leg->Draw();
  }
  TitleFactory((std::string(Generator+ " " + Target
    + " FSI Mode Shapes")).c_str())->Draw();
  return canv;
}

TCanvas* ArgueTheCaseForNE_neut(){

  TCanvas* canv = CanvasFactory("ArgueForNE_NEUT");

  TH1* T_INC_norescat_Hist =
    NoSaveClone(Data::Histos["NEUT"]["C_norescat"]["Inc"]);
  MakePDF(T_INC_norescat_Hist);

  TH1* T_INC_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["Inc"]);
  MakePDF(T_INC_Hist);

  TH1* T_NE_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["NE"]);
  MakePDF(T_NE_Hist);

  TH1* T_MEC_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["MEC"]);
  MakePDF(T_MEC_Hist);

  TH1* T_CCQELike_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["CCQELike"]);
  MakePDF(T_CCQELike_Hist);

  SetUserRangeToMaxYVal({T_INC_norescat_Hist,
                          T_INC_Hist,
                          T_NE_Hist,T_MEC_Hist,
T_CCQELike_Hist});

  T_INC_norescat_Hist->SetLineColor(kBlue+2);
  T_INC_Hist->SetLineColor(kRed-3);
  T_NE_Hist->SetLineColor(kGreen-3);
  T_MEC_Hist->SetLineColor(kOrange-2);
  T_CCQELike_Hist->SetLineColor(kMagenta-2);

  T_INC_norescat_Hist->Draw("EHIST");
  T_INC_Hist->Draw("EHIST SAME");
  T_NE_Hist->Draw("EHIST SAME");
  T_MEC_Hist->Draw("EHIST SAME");
  T_CCQELike_Hist->Draw("EHIST SAME");

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight*0.3 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(T_INC_norescat_Hist,"C (CCQE No HT)","l");
  leg->AddEntry(T_INC_Hist,"C (CCQE default HT)","l");
  leg->AddEntry(T_NE_Hist,"C (CCQE Nuclear Emission)","l");
  leg->AddEntry(T_MEC_Hist,"C (MEC)","l");
  leg->AddEntry(T_CCQELike_Hist,"C (CCQE-Like)","l");

  leg->Draw();

  canv->cd(0);

  TitleFactory("Why Nuclear Emission NEUT?")->Draw();

  return canv;
}

TCanvas* ArgueTheCaseForNE_genie(){

  TCanvas* canv = CanvasFactory("ArgueForNE_GENIE");

  TH1* T_INC_norescat_Hist =
    NoSaveClone(Data::Histos["GENIE"]["C_norescat"]["Inc"]);
  MakePDF(T_INC_norescat_Hist);

  TH1* T_INC_Hist = NoSaveClone(Data::Histos["GENIE"]["C"]["Inc"]);
  MakePDF(T_INC_Hist);

  TH1* T_NE_Hist = NoSaveClone(Data::Histos["GENIE"]["C"]["NE"]);
  MakePDF(T_NE_Hist);

  TH1* T_hN_Hist = NoSaveClone(Data::Histos["GENIE"]["C_hN"]["Inc"]);
  MakePDF(T_hN_Hist);


  SetUserRangeToMaxYVal({T_INC_norescat_Hist,
                          T_INC_Hist,
                          T_NE_Hist,
                          T_hN_Hist});

  T_INC_norescat_Hist->SetLineColor(kBlue+2);
  T_INC_Hist->SetLineColor(kRed-3);
  T_NE_Hist->SetLineColor(kGreen-3);
  T_hN_Hist->SetLineColor(kOrange-2);

  T_INC_norescat_Hist->Draw("EHIST");
  T_INC_Hist->Draw("EHIST SAME");
  T_NE_Hist->Draw("EHIST SAME");
  T_hN_Hist->Draw("EHIST SAME");

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight*0.3 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(T_INC_norescat_Hist,"C (CCQE No HT)","l");
  leg->AddEntry(T_INC_Hist,"C (CCQE default HT)","l");
  leg->AddEntry(T_hN_Hist,"C (CCQE hN Cascade HT)","l");
  leg->AddEntry(T_NE_Hist,"C (CCQE Nuclear Emission)","l");

  leg->Draw();

  canv->cd(0);

  TitleFactory("Why Nuclear Emission GENIE")->Draw();

  return canv;
}

TCanvas* DPHITNUCEMITO(std::string generator){
  TH1* GEN_C_CCQEAll = NoSaveClone(Data::Histos[generator]["O"]["CCQE"]);
  TH1* GEN_C_CCQEAllNE = NoSaveClone(Data::Histos[generator]["O"]["CCQENE"]);
  TH1* GEN_C_CCQEAllNoNE = NoSaveClone(Data::Histos[generator]["O"]["CCQENoNE"]);

  TH1* GEN_C_CCQEAll_zoom =
    NoSaveClone(Data::Histos[generator]["O"]["CCQEZoom"]);
  TH1* GEN_C_CCQEAllNE_zoom =
    NoSaveClone(Data::Histos[generator]["O"]["NEZoom"]);
  TH1* GEN_C_CCQEAllNoNE_zoom =
    NoSaveClone(Data::Histos[generator]["O"]["CCQENoNEZoom"]);

  std::initializer_list<TH1*> BigHists = {GEN_C_CCQEAll,
    GEN_C_CCQEAllNoNE,GEN_C_CCQEAllNE};

  std::initializer_list<TH1*> InsetHists = {GEN_C_CCQEAll_zoom,
    GEN_C_CCQEAllNoNE_zoom,GEN_C_CCQEAllNE_zoom};

  AddCountToTHTitle({std::make_pair("CCQE All",GEN_C_CCQEAll),
        std::make_pair("CCQE 1p1mu",GEN_C_CCQEAllNoNE),
        std::make_pair("CCQE not 1p1mu",GEN_C_CCQEAllNE)});

  MakePDF(BigHists);
  MakePDF(InsetHists);

  GEN_C_CCQEAll->SetLineStyle(1);
  GEN_C_CCQEAllNE->SetLineStyle(1);
  GEN_C_CCQEAllNoNE->SetLineStyle(1);
  GEN_C_CCQEAll_zoom->SetLineStyle(1);
  GEN_C_CCQEAllNE_zoom->SetLineStyle(1);
  GEN_C_CCQEAllNoNE_zoom->SetLineStyle(1);

  // GEN_C_CCQEAll->SetLineColor(1);
  GEN_C_CCQEAllNE->SetLineColor(kMagenta-2);
  GEN_C_CCQEAllNoNE->SetLineColor(kBlack);
  // GEN_C_CCQEAll_zoom->SetLineColor(1);
  GEN_C_CCQEAllNE_zoom->SetLineColor(kMagenta-2);
  GEN_C_CCQEAllNoNE_zoom->SetLineColor(kBlack);

  SetOnePaneHistDefaults(GEN_C_CCQEAll, DPhiT_PDF_Title, DPhiTTitle);
  SetUserRangeToMaxYVal(BigHists);
  SetUserRangeToMaxYVal(InsetHists);

  std::vector<TH1*> BigHistsV = BigHists;
  std::vector<TH1*> InsetHistsV = InsetHists;

  return GetInsetCanvas(BigHistsV,InsetHistsV);
}

TCanvas* DPHITNUCEMIT(std::string generator){
  TH1* GEN_C_CCQEAll = NoSaveClone(Data::Histos[generator]["C"]["CCQE"]);
  TH1* GEN_C_CCQEAllNE = NoSaveClone(Data::Histos[generator]["C"]["CCQENE"]);
  TH1* GEN_C_CCQEAllNoNE = NoSaveClone(Data::Histos[generator]["C"]["CCQENoNE"]);

  TH1* GEN_C_CCQEAll_zoom =
    NoSaveClone(Data::Histos[generator]["C"]["CCQEZoom"]);
  TH1* GEN_C_CCQEAllNE_zoom =
    NoSaveClone(Data::Histos[generator]["C"]["NEZoom"]);
  TH1* GEN_C_CCQEAllNoNE_zoom =
    NoSaveClone(Data::Histos[generator]["C"]["CCQENoNEZoom"]);

  std::initializer_list<TH1*> BigHists = {GEN_C_CCQEAll,
    GEN_C_CCQEAllNoNE,GEN_C_CCQEAllNE};

  std::initializer_list<TH1*> InsetHists = {GEN_C_CCQEAll_zoom,
    GEN_C_CCQEAllNoNE_zoom,GEN_C_CCQEAllNE_zoom};

  AddCountToTHTitle({std::make_pair("CCQE All",GEN_C_CCQEAll),
        std::make_pair("CCQE 1p1mu",GEN_C_CCQEAllNoNE),
        std::make_pair("CCQE not 1p1mu",GEN_C_CCQEAllNE)});

  MakePDF(BigHists);
  MakePDF(InsetHists);

  GEN_C_CCQEAll->SetLineStyle(1);
  GEN_C_CCQEAllNE->SetLineStyle(1);
  GEN_C_CCQEAllNoNE->SetLineStyle(1);
  GEN_C_CCQEAll_zoom->SetLineStyle(1);
  GEN_C_CCQEAllNE_zoom->SetLineStyle(1);
  GEN_C_CCQEAllNoNE_zoom->SetLineStyle(1);

  // GEN_C_CCQEAll->SetLineColor(1);
  GEN_C_CCQEAllNE->SetLineColor(kMagenta-2);
  GEN_C_CCQEAllNoNE->SetLineColor(kBlack);
  // GEN_C_CCQEAll_zoom->SetLineColor(1);
  GEN_C_CCQEAllNE_zoom->SetLineColor(kMagenta-2);
  GEN_C_CCQEAllNoNE_zoom->SetLineColor(kBlack);

  SetOnePaneHistDefaults(GEN_C_CCQEAll, DPhiT_PDF_Title, DPhiTTitle);
  SetUserRangeToMaxYVal(BigHists);
  SetUserRangeToMaxYVal(InsetHists);

  std::vector<TH1*> BigHistsV = BigHists;
  std::vector<TH1*> InsetHistsV = InsetHists;

  return GetInsetCanvas(BigHistsV,InsetHistsV);
}

TCanvas* DPHITCCQEGENCOMPNORESCAT(){

  TH1* NEUT_C_CCQEAll = NoSaveClone(Data::Histos["NEUT"]["C"]["CCQE"]);
  TH1* NEUT_C_CCQEAll_norescat =
    NoSaveClone(Data::Histos["NEUT"]["C_norescat"]["CCQE"]);
  TH1* GENIE_C_CCQEAll = NoSaveClone(Data::Histos["GENIE"]["C"]["CCQE"]);
  TH1* GENIE_C_CCQEAll_norescat =
    NoSaveClone(Data::Histos["GENIE"]["C_norescat"]["CCQE"]);
  TH1* GENIE_C_CCQEAll_hN =
    NoSaveClone(Data::Histos["GENIE"]["C_hN"]["CCQE"]);
  TH1* NuWro_C_CCQEAll = NoSaveClone(Data::Histos["NuWro"]["C"]["CCQE"]);
  TH1* NuWro_C_CCQEAll_norescat =
    NoSaveClone(Data::Histos["NuWro"]["C_norescat"]["CCQE"]);

  TH1* NEUT_C_CCQEAll_zoom = NoSaveClone(Data::Histos["NEUT"]["C"]["CCQEZoom"]);
  TH1* NEUT_C_CCQEAll_norescat_zoom =
    NoSaveClone(Data::Histos["NEUT"]["C_norescat"]["CCQEZoom"]);
  TH1* GENIE_C_CCQEAll_zoom = NoSaveClone(Data::Histos["GENIE"]["C"]["CCQEZoom"]);
  TH1* GENIE_C_CCQEAll_norescat_zoom =
    NoSaveClone(Data::Histos["GENIE"]["C_norescat"]["CCQEZoom"]);
  TH1* GENIE_C_CCQEAll_hN_zoom =
    NoSaveClone(Data::Histos["GENIE"]["C_hN"]["CCQEZoom"]);
  TH1* NuWro_C_CCQEAll_zoom = NoSaveClone(Data::Histos["NuWro"]["C"]["CCQEZoom"]);
  TH1* NuWro_C_CCQEAll_norescat_zoom =
    NoSaveClone(Data::Histos["NuWro"]["C_norescat"]["CCQEZoom"]);

  std::initializer_list<TH1*> BigHists = {NEUT_C_CCQEAll,GENIE_C_CCQEAll,
    NuWro_C_CCQEAll,GENIE_C_CCQEAll_hN,
    NEUT_C_CCQEAll_norescat,
    NuWro_C_CCQEAll_norescat,GENIE_C_CCQEAll_norescat};

  std::initializer_list<TH1*> InsetHists = {NEUT_C_CCQEAll_zoom,GENIE_C_CCQEAll_zoom,
    NuWro_C_CCQEAll_zoom,GENIE_C_CCQEAll_hN_zoom,
    NEUT_C_CCQEAll_norescat_zoom,
    NuWro_C_CCQEAll_norescat_zoom,GENIE_C_CCQEAll_norescat_zoom};

  AddCountToTHTitle({std::make_pair("NEUT",NEUT_C_CCQEAll),
        std::make_pair("GENIE (hA)",GENIE_C_CCQEAll),
        std::make_pair("NuWro",NuWro_C_CCQEAll),
        std::make_pair("GENIE (hN Cascade Model)",GENIE_C_CCQEAll_hN),
        std::make_pair("NEUT (NUCRESCAT = 0)",NEUT_C_CCQEAll_norescat),
        std::make_pair("GENIE (Hadronic Transport = false)",GENIE_C_CCQEAll_norescat),
        std::make_pair("NuWro (Kaskada_on = 0)",NuWro_C_CCQEAll_norescat) });

  MakePDF(BigHists);
  MakePDF(InsetHists);

  NEUT_C_CCQEAll->SetLineStyle(1);
  NEUT_C_CCQEAll_norescat->SetLineStyle(2);
  GENIE_C_CCQEAll->SetLineStyle(1);
  GENIE_C_CCQEAll_norescat->SetLineStyle(2);
  GENIE_C_CCQEAll_hN->SetLineStyle(1);
  NuWro_C_CCQEAll->SetLineStyle(1);
  NuWro_C_CCQEAll_norescat->SetLineStyle(2);

  NEUT_C_CCQEAll->SetLineColor(kRed+2);
  NEUT_C_CCQEAll_norescat->SetLineColor(kRed+2);
  GENIE_C_CCQEAll->SetLineColor(kBlue-2);
  GENIE_C_CCQEAll_norescat->SetLineColor(kBlue-2);
  GENIE_C_CCQEAll_hN->SetLineColor(kCyan-2);
  NuWro_C_CCQEAll->SetLineColor(kGreen-3);
  NuWro_C_CCQEAll_norescat->SetLineColor(kGreen-3);

  NEUT_C_CCQEAll_zoom->SetLineStyle(1);
  NEUT_C_CCQEAll_norescat_zoom->SetLineStyle(2);
  GENIE_C_CCQEAll_zoom->SetLineStyle(1);
  GENIE_C_CCQEAll_norescat_zoom->SetLineStyle(2);
  GENIE_C_CCQEAll_hN_zoom->SetLineStyle(1);
  NuWro_C_CCQEAll_zoom->SetLineStyle(1);
  NuWro_C_CCQEAll_norescat_zoom->SetLineStyle(2);

  NEUT_C_CCQEAll_zoom->SetLineColor(kRed+2);
  NEUT_C_CCQEAll_norescat_zoom->SetLineColor(kRed+2);
  GENIE_C_CCQEAll_zoom->SetLineColor(kBlue-2);
  GENIE_C_CCQEAll_norescat_zoom->SetLineColor(kBlue-2);
  GENIE_C_CCQEAll_hN_zoom->SetLineColor(kCyan-2);
  NuWro_C_CCQEAll_zoom->SetLineColor(kGreen-3);
  NuWro_C_CCQEAll_norescat_zoom->SetLineColor(kGreen-3);



  SetOnePaneHistDefaults(NEUT_C_CCQEAll, DPhiT_PDF_Title, DPhiTTitle);
  SetUserRangeToMaxYVal(BigHists);
  SetUserRangeToMaxYVal(InsetHists);

  std::vector<TH1*> BigHistsV = BigHists;
  std::vector<TH1*> InsetHistsV = InsetHists;

  return GetInsetCanvas(BigHistsV,InsetHistsV,
    0.12,0.12,
    0.55,0.5,
    0.45,0.3,
    false,true,
    0.07,0.005,
    0.02,0.135,
    0,0);

}

TCanvas* DPHITCCQEGENCOMP(){
  TH1* NEUT_C_CCQEAll = NoSaveClone(Data::Histos["NEUT"]["C"]["CCQE"]);
  TH1* GENIE_C_CCQEAll = NoSaveClone(Data::Histos["GENIE"]["C"]["CCQE"]);
  TH1* NuWro_C_CCQEAll = NoSaveClone(Data::Histos["NuWro"]["C"]["CCQE"]);
  TH1* GiBUU_C_CCQEAll = NoSaveClone(Data::Histos["GiBUU"]["C"]["CCQE"]);

  TCanvas* canv = CanvasFactory("DPHITCCQEGENCOMP");

  std::initializer_list<TH1*> hists = {NEUT_C_CCQEAll,GENIE_C_CCQEAll,
    NuWro_C_CCQEAll,GiBUU_C_CCQEAll};

  AddCountToTHTitle({std::make_pair("NEUT",NEUT_C_CCQEAll),
                      std::make_pair("GENIE",GENIE_C_CCQEAll),
                      std::make_pair("NuWro",NuWro_C_CCQEAll),
                      std::make_pair("GiBUU",GiBUU_C_CCQEAll)});

  MakePDF(hists);

  // NEUT_C_CCQEAll->SetLineColor();
  GENIE_C_CCQEAll->SetLineColor(kBlue-2);
  NuWro_C_CCQEAll->SetLineColor(kGreen-3);
  GiBUU_C_CCQEAll->SetLineColor(kOrange-2);

  NEUT_C_CCQEAll->SetLineStyle(1);
  GENIE_C_CCQEAll->SetLineStyle(1);
  NuWro_C_CCQEAll->SetLineStyle(1);
  GiBUU_C_CCQEAll->SetLineStyle(1);

  SetOnePaneHistDefaults(NEUT_C_CCQEAll, DPhiT_PDF_Title, DPhiTTitle);
  SetUserRangeToMaxYVal(hists);
  DrawTHs(hists);

  TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*1.2 - Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));

  AddToLegend(hists,leg);
  leg->Draw();
  return canv;

}

TCanvas* DPHITMODEC(){
  TH1* NEUT_C_CCQEAll = NoSaveClone(Data::Histos["NEUT"]["C"]["CCQE"]);
  TH1* NEUT_C_MEC = NoSaveClone(Data::Histos["NEUT"]["C"]["MEC"]);
  TH1* NEUT_C_RES = NoSaveClone(Data::Histos["NEUT"]["C"]["CCRes"]);
  TH1* NEUT_C_CCOTHER = NoSaveClone(Data::Histos["NEUT"]["C"]["CCOther"]);

  TCanvas* canv = CanvasFactory("DeltaPhiTModeNeutC");

  std::initializer_list<TH1*> hists = {NEUT_C_CCQEAll,
                NEUT_C_MEC,
                NEUT_C_RES,
                /*NEUT_C_CCOTHER*/};
  AddCountToTHTitle({std::make_pair("CCQE",NEUT_C_CCQEAll),
                      std::make_pair("MEC",NEUT_C_MEC),
                      std::make_pair("CCRes",NEUT_C_RES),
                      std::make_pair("CCOther",NEUT_C_CCOTHER)});
  MakePDF(hists);

  NEUT_C_MEC->SetLineColor(kBlue-2);
  NEUT_C_RES->SetLineColor(kGreen-3);
  SetOnePaneHistDefaults(NEUT_C_CCQEAll, DPhiT_PDF_Title, DPhiTTitle);
  DrawTHs(hists);

  TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*1.2 - Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));

  AddToLegend(hists,leg);
  leg->Draw();
  return canv;
}

TCanvas* GaussianExampleCanvas(){
  TCanvas* canv = SplitCanvasFactory("GaussianExampleCanvas");
  TRandom gaus;

  TH1* HThin = new TH1F("HThin","",50,0,180);
  TH1* HFat = new TH1F("HFat","",50,0,180);

  HThin->SetDirectory(0);
  HFat->SetDirectory(0);

  for(long i = 0; i < 1E6L; ++i){
    HThin->Fill(gaus.Gaus(0.0,80.0));
  }
  MakePDF(HThin);

  for(long i = 0; i < 1E6L; ++i){
    HFat->Fill(gaus.Gaus(0.0,140.0));
  }
  MakePDF(HFat);

  canv->cd(1);

  SetUpperPaneHistoDefaults(HFat, DPhiT_PDF_Title);

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight*0.3 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  SetSeriesStyle(HThin,kRed-3,1);
  SetSeriesStyle(HFat,kBlue-2,1);

  leg->AddEntry(HThin,"Thinner","l");
  leg->AddEntry(HFat,"Wider","l");

  SetUserRangeToMaxYVal({HFat,
    HThin});
  HFat->Draw("EHIST");
  HThin->Draw("EHIST SAME");

  leg->Draw();

  canv->cd(2);

  TH1* Ratio = NoSaveClone(HFat);

  Ratio->Divide(HThin);

  SetLowerPaneHistoDefaults(Ratio, "Ratio (Scaled To Denominator)", DPhiTTitle);

  Ratio->GetYaxis()->SetRangeUser(0,
    Ratio->GetBinContent(Ratio->GetNbinsX())*1.05);

  Ratio->Draw("EHIST");

  leg = LegendFactory( (Pad1XMargin),
    (1.0 - P1LegendHeight*0.3 - Pad1YMargin),
    (Pad1XMargin + P1LegendWidth),
    (1.0 - Pad1YMargin));

  leg->AddEntry(Ratio,"Fat/Thin","l");

  leg->Draw();

  return canv;
}

int MakeSomePlotzz() {

  TH1F::SetDefaultSumw2();

  TCanvas* dummyCanv = new TCanvas("dummyCanv","");
  dummyCanv->Update();dummyCanv->Print((std::string(OutputPDFFileName)+"[").c_str());

  if(!PlottingIO::InitialiseHistogramCache(HistogramCacheFileName.c_str())){
    return 1;
  }
  std::cout <<
"******************************************************************************"
"**\n\t\tDone Filling Histos.\n**"
"******************************************************************************"
<< std::endl;
  MakeT2KStyle()->cd();

  OutputTCanviFile = new TFile(OutputTCanvasFileName.c_str(),"RECREATE");

  // SaveAndPrint(DPHITMODEC());
  // SaveAndPrint(GaussianExampleCanvas());
  // SaveAndPrint(DPHITCCQEGENCOMP());
  // SaveAndPrint(DPHITCCQEGENCOMPNORESCAT());
  // SaveAndPrint(DPHITNUCEMIT("NEUT"));
  // SaveAndPrint(DPHITNUCEMIT("GENIE"));
  // SaveAndPrint(DPHITNUCEMIT("NuWro"));
  // SaveAndPrint(DPHITNUCEMIT("GiBUU"));
  // SaveAndPrint(DPHITNUCEMITO("NEUT"));
  // SaveAndPrint(DPHITNUCEMITO("GENIE"));
  // SaveAndPrint(TargetDependenceRatio("CCQE"));

  // SaveAndPrint(ArgueTheCaseForNE_neut());
  // SaveAndPrint(ArgueTheCaseForNE_genie());
  // SaveAndPrint(ArgueTheCaseForNE_nuwro());

  // SaveAndPrint(GeneratorDependence::EmissionFractions("NEUT", true));
  // SaveAndPrint(GeneratorDependence::EmissionFractions("GENIE", true));
  // SaveAndPrint(GeneratorDependence::EmissionFractions("NEUT", false));
  // SaveAndPrint(GeneratorDependence::EmissionFractions("GENIE", false));
  // // SaveAndPrint(TargetDependence::TargetDeptDPTCanvas(true));
  // // SaveAndPrint(TargetDependence::TargetDeptDPTCanvas(false));
  // // // SaveAndPrint(GeneratorDependence::GeneratorDeptDPTCanvas(true));
  // // // SaveAndPrint(GeneratorDependence::GeneratorDeptDPTCanvas(false));

  // SaveAndPrint(FSIModeCompare("NEUT","C","", true));
  // SaveAndPrint(FSIModeCompare("NEUT","C","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NEUT","C","Pt"));

  // SaveAndPrint(FSIModeCompare("GENIE","C",""));
  // SaveAndPrint(FSIModeCompare("GENIE","C","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GENIE","C","Pt"));

  // SaveAndPrint(FSIModeCompare("NEUT","Fe","", true));
  // SaveAndPrint(FSIModeCompare("NEUT","Fe","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NEUT","Fe","Pt"));

  // SaveAndPrint(FSIModeCompare("GENIE","Fe",""));
  // SaveAndPrint(FSIModeCompare("GENIE","Fe","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GENIE","Fe","Pt"));

  // SaveAndPrint(FSIModeCompare("NEUT","Pb","", true));
  // SaveAndPrint(FSIModeCompare("NEUT","Pb","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NEUT","Pb","Pt"));

  // SaveAndPrint(FSIModeCompare("GENIE","Pb",""));
  // SaveAndPrint(FSIModeCompare("GENIE","Pb","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GENIE","Pb","Pt"));

  // SaveAndPrint(FSIModeCompare("NuWro","C","", true));
  // SaveAndPrint(FSIModeCompare("NuWro","C","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NuWro","C","Pt"));

  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNProtons"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNGammma"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNNeutron"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNPiPlus"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNPiZero"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNPiMinus"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("ccqeNPrimaryParticles"));

  // SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("C"));
  // SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("O"));
  // SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("Pb"));
  // SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("Fe"));

  // SaveAndPrint(TargetDependenceInset("NE"));
  // SaveAndPrint(TargetDependenceInset("NoNE"));
  // SaveAndPrint(IronFSIModeCompare());

  // SaveAndPrint(TargetDependence::CCQEAllDistributions());

  // SaveAndPrint(TargetDependence::PionProduction("NEUT","C","",true));
  // SaveAndPrint(TargetDependence::PionProduction("NEUT","C","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("NEUT","C","Pt"));

  // SaveAndPrint(TargetDependence::PionProduction("GENIE","C","",true));
  // SaveAndPrint(TargetDependence::PionProduction("GENIE","C","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("GENIE","C","Pt"));

  // SaveAndPrint(TargetDependence::PionProduction("NEUT","Fe","",true));
  // SaveAndPrint(TargetDependence::PionProduction("NEUT","Fe","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("NEUT","Fe","Pt"));

  // SaveAndPrint(TargetDependence::PionProduction("GENIE","Fe","",true));
  // SaveAndPrint(TargetDependence::PionProduction("GENIE","Fe","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("GENIE","Fe","Pt"));

  // SaveAndPrint(TargetDependence::PionProduction("NEUT","Pb","",true));
  // SaveAndPrint(TargetDependence::PionProduction("NEUT","Pb","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("NEUT","Pb","Pt"));

  // SaveAndPrint(TargetDependence::PionProduction("GENIE","Pb","",true));
  // SaveAndPrint(TargetDependence::PionProduction("GENIE","Pb","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("GENIE","Pb","Pt"));

  // SaveAndPrint(TargetDependence::PionProduction("NuWro","C","",true));
  // SaveAndPrint(TargetDependence::PionProduction("NuWro","C","AlphaT"));
  // SaveAndPrint(TargetDependence::PionProduction("NuWro","C","Pt"));


  // SaveAndPrint(ModeDependence::EOverNE("O","GENIE"));
  dummyCanv->Print((std::string(OutputPDFFileName)+"]").c_str());
  return 0;
}


void SayUsage(const char * runcmd){
  std::cout << "Run like:\n\t" << runcmd
  << " [<OutputPDFFileName>=TransversityPlots.pdf]"
  << " [<OutputTCanvasFileName>=TCanvi.root]"
  " [<HistogramCacheFileName>=CachedHistograms.root]"
  << std::endl;
}

int main(int argc, char const * argv[]){

  std::string arg1;
  OutputPDFFileName = "TransversityPlots.pdf";
  OutputTCanvasFileName = "TCanvi.root";
  HistogramCacheFileName = "CachedHistograms.root";

  if(argc >= 2){
    arg1 = argv[1];
    if((arg1 == "-h")||(arg1 == "-?")||(arg1 == "--help")){
      SayUsage(argv[0]);
      return 0;
    }
  }
  if(argc >= 3){
    OutputTCanvasFileName = argv[2];
  }
  if(argc >= 4){
    HistogramCacheFileName = argv[3];
  }

  return MakeSomePlotzz();
}
