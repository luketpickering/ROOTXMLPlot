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
#include "TProfile.h"
#include "TRandom.h"

#include "CLITools.hxx"

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"
#include "PlottingIO.hxx"
#include "PlottingSelections.hxx"
#include "PlottingAutoPlots.hxx"

#include "NeutrinoTools_Xianguo.hxx"
#include "style_Xianguo.hxx"

using namespace PlottingDefaults;
using namespace PlottingUtils;
using namespace DataSpecifics;
using namespace PlottingSelections;
using namespace PlottingIO;

namespace {
  std::string SamplesXML = "";
  std::string SelectionsXML = "";
  std::string PlotsXML = "";

  bool MakeAutoPlots = true;
  bool MakeSpecialPlots = false;
  bool PreCache = false;
}

TH1* DummyRatioGuide(int color=kGray,
 float xMin=PlottingSelections::DPhiTMin,
 float xMax=PlottingSelections::DPhiTMax){
  TH1* dummy = new TH1F("dummy","",1,xMin, xMax);
  dummy->Fill(1);
  dummy->SetLineColor(color);
  dummy->SetLineStyle(kDotted);
  dummy->SetLineWidth(2);
  dummy->SetDirectory(0);

  return dummy;
}

void SetSeriesStyle(TH1* hist, int color, int linestyle){
  hist->SetLineWidth(2);
  hist->SetLineColor(color);
  hist->SetLineStyle(linestyle);
  hist->SetMarkerStyle(0);
}

//**********************Start Plotting Functions**************************/

namespace TargetDependence {
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

TCanvas* SpeciesCanvas(std::string Generator, std::string species,
  bool UseNoRescat=false){

  std::string CanvString = Generator+species+(UseNoRescat?"norescat":"");

  TCanvas* canv = CanvasFactory(CanvString.c_str(), false);

  TH1* GEN_C_Hist;
  TH1* GEN_Pb_Hist;

  std::initializer_list<TH1*> hists;

  if(UseNoRescat){
    GEN_C_Hist = NoSaveClone(FillHistogram1D(Generator,"C_norescat",species));
    GEN_Pb_Hist = NoSaveClone(FillHistogram1D(Generator,"Pb_norescat",species));

    hists = {GEN_C_Hist, GEN_Pb_Hist};

    AddCountToTHTitle({std::make_pair("Carbon",GEN_C_Hist),
                    std::make_pair("Lead",GEN_Pb_Hist)});
  } else {
    GEN_C_Hist = NoSaveClone(FillHistogram1D(Generator,"C",species));
    GEN_Pb_Hist = NoSaveClone(FillHistogram1D(Generator,"Pb",species));
    TH1* GEN_O_Hist = NoSaveClone(FillHistogram1D(Generator,"O",species));

    hists = {GEN_C_Hist, GEN_O_Hist, GEN_Pb_Hist};

    AddCountToTHTitle({std::make_pair("Carbon",GEN_C_Hist),
                    std::make_pair("Lead",GEN_Pb_Hist),
                    std::make_pair("Oxygen",GEN_O_Hist)});
  }

  SetOnePaneHistDefaults(*hists.begin());

  MakePDF(hists);
  SetStyles(hists);
  SetUserRangeToMaxYVal(hists,1E-3);

  GEN_C_Hist->GetXaxis()->SetTitle(species.c_str());
  GEN_C_Hist->GetYaxis()->SetTitle(("#frac{dN}{Nd("
    +std::string(species)+")}").c_str());
  GEN_C_Hist->GetXaxis()->SetTitleOffset(XTitleOffset);
  GEN_C_Hist->GetXaxis()->SetLabelOffset(XLabelOffset);

  canv->SetLogy(true);

  DrawTHs(hists);

  TitleFactory(( Generator + " " + species + " CCQE"
    + (UseNoRescat?" NoRescat":"")).c_str())->Draw();

  TLegend* leg = LegendFactory(
    (1.0 - P0LegendWidth - Pad0XMargin),
    (1.0 - (P0LegendHeight*1.25) - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin)
    );
  AddToLegend(hists,leg);

  leg->Draw();

  return canv;
}

}

namespace GeneratorDependence {

TCanvas* EmissionFractionsSep(std::string Generator, std::string Target, std::string var=""){

  std::string CanvString = Generator + Target + "NuclearEmission";

  TCanvas* canv = CanvasFactory(CanvString.c_str());

  TH1* C_Numer = NoSaveClone(FillHistogram1D(Generator,Target,"CCQENE"+var));
  TH1* C_Denom = NoSaveClone(FillHistogram1D(Generator,Target,"CCQE"+var));

  SetUserRangeToMaxYVal({C_Numer,C_Denom});
  AddCountToTHTitle({
    std::make_pair("CCQE With Nuclear Emission",C_Numer),
    std::make_pair("CCQE All",C_Denom)});

  TLegend* leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));

  C_Numer->SetMarkerStyle(0);
  C_Numer->SetLineWidth(2);
  C_Numer->SetLineColor(kRed);
  C_Denom->SetMarkerStyle(0);
  C_Denom->SetLineWidth(2);
  C_Denom->SetLineColor(kBlue);

  AddToLegend({C_Numer,C_Denom},leg);
  SetOnePaneHistDefaults(C_Numer,"Count",(var=="DPt")?DPTTitle:MuPtTitle);

  DrawTHs({C_Numer,C_Denom});

  leg->Draw();

  TitleFactory(CanvString)->Draw();

  return canv;
}

TCanvas* EmissionFractionsRatio(std::string Generator, std::string Target, std::string var=""){

  std::string XTitleText = "";
  if(var==""){
    XTitleText = DPhiTTitle;
  } else if(var=="AlphaT"){
    XTitleText = DAlphaTTitle;
  } else if(var=="DPt"){
    XTitleText = DPTTitle;
  }  else if(var=="MuPt"){
    XTitleText = MuPtTitle;
  }else {
    return new TCanvas("broken");
  }

  std::string CanvString = Generator + Target + "NuclearEmissionRatio";

  TCanvas* canv = CanvasFactory(CanvString.c_str());

  TH1* C_Numer = NoSaveClone(FillHistogram1D(Generator,Target,"CCQENE"+var));
  TH1* C_Denom = NoSaveClone(FillHistogram1D(Generator,Target,"CCQE"+var));

  TH1* C_Ratio = NoSaveClone(C_Numer);
  C_Ratio->Divide(C_Numer, C_Denom,1,1,"B");

  C_Ratio->SetMarkerStyle(0);
  C_Ratio->SetLineWidth(2);
  C_Ratio->SetLineColor(kOrange-2);

  SetOnePaneHistDefaults(C_Ratio,"Emission Fraction",XTitleText.c_str());

  DrawTHs({C_Ratio});

  TitleFactory(CanvString)->Draw();

  return canv;
}

TCanvas* EmissionFractionsGenComp(std::string Target, std::string var="",
  bool regen=false){

  std::string XTitleText = "";
  if(var==""){
    XTitleText = DPhiTTitle;
  } else if(var=="AlphaT"){
    XTitleText = DAlphaTTitle;
  } else if(var=="DPt"){
    XTitleText = DPTTitle;
  }  else if(var=="MuPt"){
    XTitleText = MuPtTitle;
  }else {
    return new TCanvas("broken");
  }

  std::string CanvString = Target + "NuclearEmissionRatio";

  TCanvas* canv = CanvasFactory(CanvString.c_str());

  TH1* NEUT_C_Numer = NoSaveClone(FillHistogram1D("NEUT",Target,"CCQENE"+var,regen));
  TH1* NEUT_C_Denom = NoSaveClone(FillHistogram1D("NEUT",Target,"CCQE"+var,regen));

  TH1* NEUT_C_Ratio = NoSaveClone(NEUT_C_Numer);
  NEUT_C_Ratio->Divide(NEUT_C_Numer, NEUT_C_Denom,1,1,"B");

  TH1* GENIE_C_Numer = NoSaveClone(FillHistogram1D("GENIE",Target,"CCQENE"+var,regen));
  TH1* GENIE_C_Denom = NoSaveClone(FillHistogram1D("GENIE",Target,"CCQE"+var,regen));

  TH1* GENIE_C_Ratio = NoSaveClone(GENIE_C_Numer);
  GENIE_C_Ratio->Divide(GENIE_C_Numer, GENIE_C_Denom,1,1,"B");

  TH1* NuWro_C_Numer = NoSaveClone(FillHistogram1D("NuWro",Target,"CCQENE"+var,regen));
  TH1* NuWro_C_Denom = NoSaveClone(FillHistogram1D("NuWro",Target,"CCQE"+var,regen));

  TH1* NuWro_C_Ratio = NoSaveClone(NuWro_C_Numer);
  NuWro_C_Ratio->Divide(NuWro_C_Numer, NuWro_C_Denom,1,1,"B");

  TH1* GiBUU_C_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQENE"+var,regen));
  TH1* GiBUU_C_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQE"+var,regen));

  TH1* GiBUU_C_Ratio = NoSaveClone(GiBUU_C_Numer);
  GiBUU_C_Ratio->Divide(GiBUU_C_Numer, GiBUU_C_Denom,1,1,"B");

  SetStyles({NEUT_C_Ratio,
            GENIE_C_Ratio,
            NuWro_C_Ratio,
            GiBUU_C_Ratio});

  NEUT_C_Ratio->SetTitle("NEUT");
  GENIE_C_Ratio->SetTitle("GENIE");
  NuWro_C_Ratio->SetTitle("NuWro");
  GiBUU_C_Ratio->SetTitle("GiBUU");

  SetOnePaneHistDefaults(NEUT_C_Ratio,"Emission Fraction",
    XTitleText.c_str());

  SetUserRangeToMaxYVal({NEUT_C_Ratio,
            GENIE_C_Ratio,
            NuWro_C_Ratio,
            GiBUU_C_Ratio});
  SetUserRangeToMaxXValWithContent({NEUT_C_Ratio,
            GENIE_C_Ratio,
            NuWro_C_Ratio,
            GiBUU_C_Ratio});

  DrawTHs({NEUT_C_Ratio,
            GENIE_C_Ratio,
            NuWro_C_Ratio,
            GiBUU_C_Ratio});

  TLegend* leg = LegendFactory( (Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (Pad0XMargin + P0LegendWidth),
    (1.0 - Pad0YMargin));

  AddToLegend({NEUT_C_Ratio,
            GENIE_C_Ratio,
            NuWro_C_Ratio,
            GiBUU_C_Ratio},leg);

  TitleFactory(CanvString)->Draw();
  leg->Draw();

  return canv;
}

TCanvas* EmissionFractionsGenComp_Xianguo(std::string Target, std::string var="",
  bool regen=false, float xmin=0, float xmax=0xDEADBEEF){

  std::string XTitleText = "";
  if(var==""){
    XTitleText = DPhiTTitle;
  } else if(var=="AlphaT"){
    XTitleText = DAlphaTTitle;
  } else if(var=="DPt"){
    XTitleText = DPTTitle;
  }  else if(var=="MuPt"){
    XTitleText = MuPtTitle;
  }else {
    return new TCanvas("broken");
  }

  std::string CanvString = Target + "NuclearEmissionRatio";

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas(CanvString.c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);


  TH1* NEUT_C_Numer = NoSaveClone(FillHistogram1D("NEUT",Target,"CCQENE"+var,regen));
  TH1* NEUT_C_Denom = NoSaveClone(FillHistogram1D("NEUT",Target,"CCQE"+var,regen));

  TH1* NEUT_C_Ratio = NoSaveClone(NEUT_C_Numer);
  NEUT_C_Ratio->Divide(NEUT_C_Numer, NEUT_C_Denom,1,1,"B");

  TH1* GENIE_C_Numer = NoSaveClone(FillHistogram1D("GENIE",Target,"CCQENE"+var,regen));
  TH1* GENIE_C_Denom = NoSaveClone(FillHistogram1D("GENIE",Target,"CCQE"+var,regen));

  TH1* GENIE_C_Ratio = NoSaveClone(GENIE_C_Numer);
  GENIE_C_Ratio->Divide(GENIE_C_Numer, GENIE_C_Denom,1,1,"B");

  TH1* NuWro_C_Numer = NoSaveClone(FillHistogram1D("NuWro",Target,"CCQENE"+var,regen));
  TH1* NuWro_C_Denom = NoSaveClone(FillHistogram1D("NuWro",Target,"CCQE"+var,regen));

  TH1* NuWro_C_Ratio = NoSaveClone(NuWro_C_Numer);
  NuWro_C_Ratio->Divide(NuWro_C_Numer, NuWro_C_Denom,1,1,"B");

  // TH1* GiBUU_C_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQENE"+var,regen));
  // TH1* GiBUU_C_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQE"+var,regen));

  // TH1* GiBUU_C_Ratio = NoSaveClone(GiBUU_C_Numer);
  // GiBUU_C_Ratio->Divide(GiBUU_C_Numer, GiBUU_C_Denom,1,1,"B");

  TH1* NEUT_C_2GeV_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_2GeV","CCQENE"+var,regen));
  TH1* NEUT_C_2GeV_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_2GeV","CCQE"+var,regen));

  TH1* NEUT_C_2GeV_Ratio = NoSaveClone(NEUT_C_2GeV_Numer);
  NEUT_C_2GeV_Ratio->Divide(NEUT_C_2GeV_Numer, NEUT_C_2GeV_Denom,1,1,"B");

  // TH1* GENIE_C_2GeV_Numer = NoSaveClone(FillHistogram1D("GENIE",Target+"_2GeV","CCQENE"+var,regen));
  // TH1* GENIE_C_2GeV_Denom = NoSaveClone(FillHistogram1D("GENIE",Target+"_2GeV","CCQE"+var,regen));

  // TH1* GENIE_C_2GeV_Ratio = NoSaveClone(GENIE_C_2GeV_Numer);
  // GENIE_C_2GeV_Ratio->Divide(GENIE_C_2GeV_Numer, GENIE_C_2GeV_Denom,1,1,"B");

  TH1* NuWro_C_2GeV_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_2GeV","CCQENE"+var,regen));
  TH1* NuWro_C_2GeV_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_2GeV","CCQE"+var,regen));

  TH1* NuWro_C_2GeV_Ratio = NoSaveClone(NuWro_C_2GeV_Numer);
  NuWro_C_2GeV_Ratio->Divide(NuWro_C_2GeV_Numer, NuWro_C_2GeV_Denom,1,1,"B");

  // TH1* GiBUU_C_2GeV_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQENE"+var,regen));
  // TH1* GiBUU_C_2GeV_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQE"+var,regen));

  // TH1* GiBUU_C_2GeV_Ratio = NoSaveClone(GiBUU_C_2GeV_Numer);
  // GiBUU_C_2GeV_Ratio->Divide(GiBUU_C_2GeV_Numer, GiBUU_C_2GeV_Denom,1,1,"B");


  TH1* NEUT_C_8GeV_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_8GeV","CCQENE"+var,regen));
  TH1* NEUT_C_8GeV_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_8GeV","CCQE"+var,regen));

  TH1* NEUT_C_8GeV_Ratio = NoSaveClone(NEUT_C_8GeV_Numer);
  NEUT_C_8GeV_Ratio->Divide(NEUT_C_8GeV_Numer, NEUT_C_8GeV_Denom,1,1,"B");

  // TH1* GENIE_C_8GeV_Numer = NoSaveClone(FillHistogram1D("GENIE",Target+"_8GeV","CCQENE"+var,regen));
  // TH1* GENIE_C_8GeV_Denom = NoSaveClone(FillHistogram1D("GENIE",Target+"_8GeV","CCQE"+var,regen));

  // TH1* GENIE_C_8GeV_Ratio = NoSaveClone(GENIE_C_8GeV_Numer);
  // GENIE_C_8GeV_Ratio->Divide(GENIE_C_8GeV_Numer, GENIE_C_8GeV_Denom,1,1,"B");

  TH1* NuWro_C_8GeV_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_8GeV","CCQENE"+var,regen));
  TH1* NuWro_C_8GeV_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_8GeV","CCQE"+var,regen));

  TH1* NuWro_C_8GeV_Ratio = NoSaveClone(NuWro_C_8GeV_Numer);
  NuWro_C_8GeV_Ratio->Divide(NuWro_C_8GeV_Numer, NuWro_C_8GeV_Denom,1,1,"B");

  // TH1* GiBUU_C_8GeV_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQENE"+var,regen));
  // TH1* GiBUU_C_8GeV_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQE"+var,regen));

  // TH1* GiBUU_C_8GeV_Ratio = NoSaveClone(GiBUU_C_8GeV_Numer);
  // GiBUU_C_8GeV_Ratio->Divide(GiBUU_C_8GeV_Numer, GiBUU_C_8GeV_Denom,1,1,"B");

  std::initializer_list<TH1*> Hists = {
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    // GiBUU_C_Ratio
    NEUT_C_2GeV_Ratio,
    // GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    // GiBUU_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio,
    // GENIE_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
    // GiBUU_C_8GeV_Ratio
  };

  std::vector<TH1*> Hists_v = {
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    // GiBUU_C_Ratio
    NEUT_C_2GeV_Ratio, 0,
    // GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    // GiBUU_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio, 0,
    // GENIE_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
    // GiBUU_C_8GeV_Ratio
  };

  NEUT_C_Ratio->SetTitle("Neut");
  GENIE_C_Ratio->SetTitle("GENIE");
  NuWro_C_Ratio->SetTitle("NuWro");

  NEUT_C_2GeV_Ratio->SetTitle("Neut");
  // GENIE_C_2GeV_Ratio->SetTitle("GENIE");
  NuWro_C_2GeV_Ratio->SetTitle("NuWro");
  // GiBUU_C_Ratio->SetTitle("GiBUU");

  NEUT_C_2GeV_Ratio->GetXaxis()->SetTitle(XTitleText.c_str());
  NEUT_C_2GeV_Ratio->GetYaxis()->SetTitle("#it{P}_{NE}");

  SetUserRangeToMaxYVal(Hists);
  if(xmax!=0xDEADBEEF){
    NEUT_C_2GeV_Ratio->GetXaxis()->SetRangeUser(xmin,xmax);
  } else {
    SetUserRangeToMaxXValWithContent(Hists);
  }

  for(auto& histo : Hists){
    style::ResetStyle((TH1*)(histo));
  }

  Hists_v[0]->SetLineColor(kBlack);
  Hists_v[1]->SetLineColor(kBlack);
  Hists_v[2]->SetLineColor(kBlack);

  Hists_v[0]->SetLineStyle(kSolid);
  Hists_v[1]->SetLineStyle(kDashed);
  Hists_v[2]->SetLineStyle(kDotted);

  // Hists_v[3]->SetLineWidth(2);
  // Hists_v[4]->SetLineWidth(2);
  // Hists_v[5]->SetLineWidth(2);

  Hists_v[3]->SetLineColor(kRed);
  // Hists_v[4]->SetLineColor(kRed);
  Hists_v[5]->SetLineColor(kRed);

  Hists_v[3]->SetLineStyle(kSolid);
  // Hists_v[4]->SetLineStyle(kDashed);
  Hists_v[5]->SetLineStyle(kDotted);

  Hists_v[6]->SetLineColor(kGreen);
  // Hists_v[7]->SetLineColor(kGreen);
  Hists_v[8]->SetLineColor(kGreen);

  Hists_v[6]->SetLineStyle(kSolid);
  // Hists_v[7]->SetLineStyle(kDashed);
  Hists_v[8]->SetLineStyle(kDotted);

  DrawTHs({
    NEUT_C_2GeV_Ratio,
    // GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    // GiBUU_C_2GeV_Ratio
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    // GiBUU_C_Ratio
    NEUT_C_8GeV_Ratio,
    // GENIE_C_Rato,
    NuWro_C_8GeV_Ratio,
    // GiBUU_C_Ratio
  }, "HIST");

  TLegend* leg, *leg2;

  leg = LegendFactory( (0.12 + P0LegendWidth*0.1),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*1.1),
    (1.0 - Pad0YMargin*0.2));

  leg2 = LegendFactory( (0.12 + P0LegendWidth*1.2),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*2.2),
    (1.0 - Pad0YMargin*0.2));


  leg->SetHeader("0.6 GeV");
  leg2->SetHeader("2.0 GeV");

  style::ResetStyle(leg);
  style::ResetStyle(leg2);

  AddToLegend({
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    },leg);
  AddToLegend({
    NEUT_C_2GeV_Ratio,
    // GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    },leg2);

  leg->Draw();
  leg2->Draw();

  return canv;
}

TCanvas* EmissionFractionsGenComp_Xianguo_thresh(
  std::string Target, std::string selBase, std::string lt, std::string ht, std::string var="",
  bool regen=false, float xmin=0, float xmax=0xDEADBEEF, std::string genieVer="2_9_0"){

  std::string XTitleText = "";
  if(var==""){
    XTitleText = DPhiTTitle;
  } else if(var=="AlphaT"){
    XTitleText = DAlphaTTitle;
  } else if(var=="DPt"){
    XTitleText = DPTTitle;
  }  else if(var=="MuPt"){
    XTitleText = MuPtTitle;
  }else {
    return new TCanvas("broken");
  }

  std::string CanvString = Target + "NuclearEmissionRatio";

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas(CanvString.c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);


  TH1* NEUT_C_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NEUT_C_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NEUT_C_Ratio = NoSaveClone(NEUT_C_Numer);
  NEUT_C_Ratio->Divide(NEUT_C_Numer, NEUT_C_Denom,1,1,"B");

  TH1* GENIE_C_Numer = NoSaveClone(FillHistogram1D("GENIE_"+genieVer,Target+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* GENIE_C_Denom = NoSaveClone(FillHistogram1D("GENIE_"+genieVer,Target+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* GENIE_C_Ratio = NoSaveClone(GENIE_C_Numer);
  GENIE_C_Ratio->Divide(GENIE_C_Numer, GENIE_C_Denom,1,1,"B");

  TH1* NuWro_C_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NuWro_C_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NuWro_C_Ratio = NoSaveClone(NuWro_C_Numer);
  NuWro_C_Ratio->Divide(NuWro_C_Numer, NuWro_C_Denom,1,1,"B");

  TH1* GiBUU_C_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* GiBUU_C_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* GiBUU_C_Ratio = NoSaveClone(GiBUU_C_Numer);
  GiBUU_C_Ratio->Divide(GiBUU_C_Numer, GiBUU_C_Denom,1,1,"B");

  TH1* NEUT_C_2GeV_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NEUT_C_2GeV_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NEUT_C_2GeV_Ratio = NoSaveClone(NEUT_C_2GeV_Numer);
  NEUT_C_2GeV_Ratio->Divide(NEUT_C_2GeV_Numer, NEUT_C_2GeV_Denom,1,1,"B");

  TH1* GENIE_C_2GeV_Numer = NoSaveClone(FillHistogram1D("GENIE_"+genieVer,Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* GENIE_C_2GeV_Denom = NoSaveClone(FillHistogram1D("GENIE_"+genieVer,Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* GENIE_C_2GeV_Ratio = NoSaveClone(GENIE_C_2GeV_Numer);
  GENIE_C_2GeV_Ratio->Divide(GENIE_C_2GeV_Numer, GENIE_C_2GeV_Denom,1,1,"B");

  TH1* NuWro_C_2GeV_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NuWro_C_2GeV_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NuWro_C_2GeV_Ratio = NoSaveClone(NuWro_C_2GeV_Numer);
  NuWro_C_2GeV_Ratio->Divide(NuWro_C_2GeV_Numer, NuWro_C_2GeV_Denom,1,1,"B");

  TH1* GiBUU_C_2GeV_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* GiBUU_C_2GeV_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* GiBUU_C_2GeV_Ratio = NoSaveClone(GiBUU_C_2GeV_Numer);
  GiBUU_C_2GeV_Ratio->Divide(GiBUU_C_2GeV_Numer, GiBUU_C_2GeV_Denom,1,1,"B");

  TH1* NEUT_C_8GeV_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NEUT_C_8GeV_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NEUT_C_8GeV_Ratio = NoSaveClone(NEUT_C_8GeV_Numer);
  NEUT_C_8GeV_Ratio->Divide(NEUT_C_8GeV_Numer, NEUT_C_8GeV_Denom,1,1,"B");

  TH1* GENIE_C_8GeV_Numer = NoSaveClone(FillHistogram1D("GENIE_"+genieVer,Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* GENIE_C_8GeV_Denom = NoSaveClone(FillHistogram1D("GENIE_"+genieVer,Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* GENIE_C_8GeV_Ratio = NoSaveClone(GENIE_C_8GeV_Numer);
  GENIE_C_8GeV_Ratio->Divide(GENIE_C_8GeV_Numer, GENIE_C_8GeV_Denom,1,1,"B");

  TH1* NuWro_C_8GeV_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NuWro_C_8GeV_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NuWro_C_8GeV_Ratio = NoSaveClone(NuWro_C_8GeV_Numer);
  NuWro_C_8GeV_Ratio->Divide(NuWro_C_8GeV_Numer, NuWro_C_8GeV_Denom,1,1,"B");

  TH1* GiBUU_C_8GeV_Numer = NoSaveClone(FillHistogram1D("GiBUU",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* GiBUU_C_8GeV_Denom = NoSaveClone(FillHistogram1D("GiBUU",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* GiBUU_C_8GeV_Ratio = NoSaveClone(GiBUU_C_8GeV_Numer);
  GiBUU_C_8GeV_Ratio->Divide(GiBUU_C_8GeV_Numer, GiBUU_C_8GeV_Denom,1,1,"B");

  std::initializer_list<TH1*> Hists = {
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    GiBUU_C_Ratio,
    NEUT_C_2GeV_Ratio,
    GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    GiBUU_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio,
    GENIE_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
    GiBUU_C_8GeV_Ratio
  };

  std::vector<TH1*> Hists_v = {
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    GiBUU_C_Ratio,
    NEUT_C_2GeV_Ratio,
    GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    GiBUU_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio,
    GENIE_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
    GiBUU_C_8GeV_Ratio
  };

  NEUT_C_8GeV_Ratio->GetXaxis()->SetTitle(XTitleText.c_str());
  NEUT_C_8GeV_Ratio->GetYaxis()->SetTitle("#it{P}_{NE}");

  SetUserRangeToMaxYVal(Hists);
  if(xmax!=0xDEADBEEF){
    NEUT_C_8GeV_Ratio->GetXaxis()->SetRangeUser(xmin,xmax);
  } else {
    SetUserRangeToMaxXValWithContent(Hists);
  }

  for(auto& histo : Hists){
    style::ResetStyle((TH1*)(histo));
  }

  Hists_v[0]->SetLineWidth(2);
  Hists_v[1]->SetLineWidth(2);
  Hists_v[2]->SetLineWidth(2);
  Hists_v[3]->SetLineWidth(2);

  Hists_v[0]->SetLineColor(kBlack);
  Hists_v[1]->SetLineColor(kBlue);
  Hists_v[2]->SetLineColor(kRed);
  Hists_v[3]->SetLineColor(kMagenta);

  Hists_v[0]->SetLineStyle(kSolid);
  Hists_v[1]->SetLineStyle(kSolid);
  Hists_v[2]->SetLineStyle(kSolid);
  Hists_v[3]->SetLineStyle(kSolid);

  Hists_v[4]->SetLineWidth(2);
  Hists_v[5]->SetLineWidth(2);
  Hists_v[6]->SetLineWidth(2);
  Hists_v[7]->SetLineWidth(2);

  Hists_v[4]->SetLineColor(kBlack);
  Hists_v[5]->SetLineColor(kBlue);
  Hists_v[6]->SetLineColor(kRed);
  Hists_v[7]->SetLineColor(kMagenta);

  Hists_v[4]->SetLineStyle(kDashed);
  Hists_v[5]->SetLineStyle(kDashed);
  Hists_v[6]->SetLineStyle(kDashed);
  Hists_v[7]->SetLineStyle(kDashed);

  Hists_v[8]->SetLineWidth(2);
  Hists_v[9]->SetLineWidth(2);
  Hists_v[10]->SetLineWidth(2);
  Hists_v[11]->SetLineWidth(2);

  Hists_v[8]->SetLineColor(kBlack);
  Hists_v[9]->SetLineColor(kBlue);
  Hists_v[10]->SetLineColor(kRed);
  Hists_v[11]->SetLineColor(kMagenta);

  Hists_v[8]->SetLineStyle(kDotted);
  Hists_v[9]->SetLineStyle(kDotted);
  Hists_v[10]->SetLineStyle(kDotted);
  Hists_v[11]->SetLineStyle(kDotted);

  DrawTHs({
    NEUT_C_8GeV_Ratio,
    GENIE_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
    GiBUU_C_8GeV_Ratio,
    NEUT_C_2GeV_Ratio,
    GENIE_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    GiBUU_C_2GeV_Ratio,
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    GiBUU_C_Ratio,
  }, "HIST");

  TLegend* leg, *leg2;

  leg = LegendFactory( (0.12 + P0LegendWidth*0.1),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*1.1),
    (1.0 - Pad0YMargin*0.2));

  leg2 = LegendFactory( (0.12 + P0LegendWidth*1.2),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*2.2),
    (1.0 - Pad0YMargin*0.2));

  style::ResetStyle(leg);
  style::ResetStyle(leg2);

  NEUT_C_Ratio->SetTitle("Neut");
  GENIE_C_Ratio->SetTitle("GENIE");
  NuWro_C_Ratio->SetTitle("NuWro");
  GiBUU_C_Ratio->SetTitle("GiBUU");
  AddToLegend({
    NEUT_C_Ratio,
    GENIE_C_Ratio,
    NuWro_C_Ratio,
    GiBUU_C_Ratio
    },leg);


  NEUT_C_Ratio->SetTitle("0.6 GeV");
  NEUT_C_2GeV_Ratio->SetTitle("2.0 GeV");
  NEUT_C_8GeV_Ratio->SetTitle("8.0 GeV");

  AddToLegend({
    NEUT_C_Ratio,
    NEUT_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio,
    },leg2);

  leg->Draw();
  leg2->Draw();

  return canv;
}

TCanvas* EmissionFractionsGenComp_Xianguo_thresh_CC0Pi(
  std::string Target, std::string selBase, std::string lt, std::string ht, std::string var="",
  bool regen=false, float xmin=0, float xmax=0xDEADBEEF){

  std::string XTitleText = "";
  if(var==""){
    XTitleText = DPhiTTitle;
  } else if(var=="AlphaT"){
    XTitleText = DAlphaTTitle;
  } else if(var=="DPt"){
    XTitleText = DPTTitle;
  }  else if(var=="MuPt"){
    XTitleText = MuPtTitle;
  }else {
    return new TCanvas("broken");
  }

  std::string CanvString = Target + "NuclearEmissionRatio_CC0Pi";

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas(CanvString.c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);


  TH1* NEUT_C_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NEUT_C_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NEUT_C_Ratio = NoSaveClone(NEUT_C_Numer);
  NEUT_C_Ratio->Divide(NEUT_C_Numer, NEUT_C_Denom,1,1,"B");

  TH1* NuWro_C_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NuWro_C_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NuWro_C_Ratio = NoSaveClone(NuWro_C_Numer);
  NuWro_C_Ratio->Divide(NuWro_C_Numer, NuWro_C_Denom,1,1,"B");

  TH1* NEUT_C_2GeV_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NEUT_C_2GeV_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NEUT_C_2GeV_Ratio = NoSaveClone(NEUT_C_2GeV_Numer);
  NEUT_C_2GeV_Ratio->Divide(NEUT_C_2GeV_Numer, NEUT_C_2GeV_Denom,1,1,"B");

  TH1* NuWro_C_2GeV_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NuWro_C_2GeV_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_2GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NuWro_C_2GeV_Ratio = NoSaveClone(NuWro_C_2GeV_Numer);
  NuWro_C_2GeV_Ratio->Divide(NuWro_C_2GeV_Numer, NuWro_C_2GeV_Denom,1,1,"B");

  TH1* NEUT_C_8GeV_Numer = NoSaveClone(FillHistogram1D("NEUT",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NEUT_C_8GeV_Denom = NoSaveClone(FillHistogram1D("NEUT",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NEUT_C_8GeV_Ratio = NoSaveClone(NEUT_C_8GeV_Numer);
  NEUT_C_8GeV_Ratio->Divide(NEUT_C_8GeV_Numer, NEUT_C_8GeV_Denom,1,1,"B");

  TH1* NuWro_C_8GeV_Numer = NoSaveClone(FillHistogram1D("NuWro",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+"NE"+var,regen));
  TH1* NuWro_C_8GeV_Denom = NoSaveClone(FillHistogram1D("NuWro",Target+"_8GeV"+"_lt"+lt+"_ht"+ht,selBase+""+var,regen));

  TH1* NuWro_C_8GeV_Ratio = NoSaveClone(NuWro_C_8GeV_Numer);
  NuWro_C_8GeV_Ratio->Divide(NuWro_C_8GeV_Numer, NuWro_C_8GeV_Denom,1,1,"B");


  std::initializer_list<TH1*> Hists = {
    NEUT_C_Ratio,
    NuWro_C_Ratio,
    NEUT_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
  };

  std::vector<TH1*> Hists_v = {
    NEUT_C_Ratio,
    NuWro_C_Ratio,
    NEUT_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    NEUT_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio
  };

  NEUT_C_Ratio->SetTitle("Neut");
  NuWro_C_Ratio->SetTitle("NuWro");

  NEUT_C_8GeV_Ratio->SetTitle("Neut");
  NuWro_C_8GeV_Ratio->SetTitle("NuWro");

  NEUT_C_8GeV_Ratio->GetXaxis()->SetTitle(XTitleText.c_str());
  NEUT_C_8GeV_Ratio->GetYaxis()->SetTitle("#it{P}_{NE}");

  SetUserRangeToMaxYVal(Hists);
  if(xmax!=0xDEADBEEF){
    NEUT_C_2GeV_Ratio->GetXaxis()->SetRangeUser(xmin,xmax);
  } else {
    SetUserRangeToMaxXValWithContent(Hists);
  }

  for(auto& histo : Hists){
    style::ResetStyle((TH1*)(histo));
  }

  Hists_v[0]->SetLineWidth(2);
  Hists_v[1]->SetLineWidth(2);

  Hists_v[0]->SetLineColor(kBlack);
  Hists_v[1]->SetLineColor(kRed);

  Hists_v[0]->SetLineStyle(kSolid);
  Hists_v[1]->SetLineStyle(kSolid);

  Hists_v[2]->SetLineWidth(2);
  Hists_v[3]->SetLineWidth(2);

  Hists_v[2]->SetLineColor(kBlack);
  Hists_v[3]->SetLineColor(kRed);

  Hists_v[2]->SetLineStyle(kDashed);
  Hists_v[3]->SetLineStyle(kDashed);


  DrawTHs({
    NEUT_C_8GeV_Ratio,
    NuWro_C_8GeV_Ratio,
    NEUT_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    NEUT_C_Ratio,
    NuWro_C_Ratio,
  }, "HIST");

  TLegend* leg, *leg2;

  leg = LegendFactory( (0.12 + P0LegendWidth*0.1),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*1.1),
    (1.0 - Pad0YMargin*0.2));

  leg2 = LegendFactory( (0.12 + P0LegendWidth*1.2),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*2.2),
    (1.0 - Pad0YMargin*0.2));


  leg->SetHeader("0.6 GeV");
  leg2->SetHeader("2.0 GeV");

  style::ResetStyle(leg);
  style::ResetStyle(leg2);

  AddToLegend({
    NEUT_C_Ratio,
    NuWro_C_Ratio,
    },leg);
  AddToLegend({
    NEUT_C_2GeV_Ratio,
    NuWro_C_2GeV_Ratio,
    },leg2);

  leg->Draw();
  leg2->Draw();

  return canv;
}

TCanvas* GenComp_STV_Xianguo(std::string Target, std::string STVar,
  bool regen=false){

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas((STVar+Target+"vsmupt").c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar==""){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(STVar=="DPt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  } else {
    return new TCanvas("broken");
  }

  TH1* CCQEAll_NEUT =
    NoSaveClone(FillHistogram1D("NEUT",Target+"_norescat",
      "CCQE"+STVar,regen));
  TH1* CCQEAll_GENIE =
    NoSaveClone(FillHistogram1D("GENIE",Target+"_norescat",
      "CCQE"+STVar,regen));
  TH1* CCQEAll_NuWro =
    NoSaveClone(FillHistogram1D("NuWro",Target+"_norescat",
      "CCQE"+STVar,regen));

  TH1* CCQEAll_NEUT_FSI =
    NoSaveClone(FillHistogram1D("NEUT",Target,
      "CCQE"+STVar,regen));
  TH1* CCQEAll_GENIE_FSI =
    NoSaveClone(FillHistogram1D("GENIE",Target,
      "CCQE"+STVar,regen));
  TH1* CCQEAll_NuWro_FSI =
    NoSaveClone(FillHistogram1D("NuWro",Target,
      "CCQE"+STVar,regen));

  std::initializer_list<TH1*> hists = {
    CCQEAll_NEUT_FSI,
    CCQEAll_GENIE_FSI,
    CCQEAll_NuWro_FSI,
    CCQEAll_NEUT,
    CCQEAll_GENIE,
    CCQEAll_NuWro,
  };

  CCQEAll_NEUT->SetTitle("Neut (RFG)");
  CCQEAll_GENIE->SetTitle("GENIE (Bodek-Ritchie)");
  CCQEAll_NuWro->SetTitle("NuWro (LFG)");

  CCQEAll_NEUT_FSI->SetTitle("Neut (RFG)");
  CCQEAll_GENIE_FSI->SetTitle("GENIE (Bodek-Ritchie)");
  CCQEAll_NuWro_FSI->SetTitle("NuWro (LFG)");

  MakePDF(hists);

  (*hists.begin())->GetXaxis()->SetTitle(XTitleText.c_str());
  (*hists.begin())->GetYaxis()->SetTitle(YTitleText.c_str());

  for(auto& histo : hists){
    style::ResetStyle((TH1*)(histo));
  }

  std::vector<TH1*> Hists_v = hists;

  Hists_v[3]->SetLineColor(kBlack);
  Hists_v[4]->SetLineColor(kBlack);
  Hists_v[5]->SetLineColor(kBlack);

  Hists_v[3]->SetLineStyle(kSolid);
  Hists_v[4]->SetLineStyle(kDashed);
  Hists_v[5]->SetLineStyle(kDotted);

  // Hists_v[0]->SetLineWidth(2);
  // Hists_v[1]->SetLineWidth(2);
  // Hists_v[2]->SetLineWidth(2);

  Hists_v[0]->SetLineColor(kGreen);
  Hists_v[1]->SetLineColor(kGreen);
  Hists_v[2]->SetLineColor(kGreen);

  Hists_v[0]->SetLineStyle(kSolid);
  Hists_v[1]->SetLineStyle(kDashed);
  Hists_v[2]->SetLineStyle(kDotted);

  SetUserRangeToMaxYVal(hists);

  DrawTHs(hists, "HIST");

  TLegend* leg, *leg2;
  if(STVar==""){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*2.8),
    (1.0 - P0LegendHeight*2.0 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else if(STVar=="AlphaT"){

    leg = LegendFactory( (0.12),
      (1.0 - P0LegendHeight - Pad0YMargin*0.2),
      (0.12 + 0.35),
      (1.0 - Pad0YMargin*0.2));

    leg2 = LegendFactory( (0.12 + 0.36),
      (1.0 - P0LegendHeight - Pad0YMargin*0.2),
      (0.12 + 0.36 + 0.35),
      (1.0 - Pad0YMargin*0.2));

  } else if(STVar=="DPt"){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*2.0),
    (1.0 - P0LegendHeight*1.3 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  }

  leg->SetHeader("No FSI");
  leg2->SetHeader("With FSI");

  style::ResetStyle(leg);
  style::ResetStyle(leg2);

  AddToLegend({
    CCQEAll_NEUT,
    CCQEAll_GENIE,
    CCQEAll_NuWro,
  },leg);

  AddToLegend({
    CCQEAll_NEUT_FSI,
    CCQEAll_GENIE_FSI,
    CCQEAll_NuWro_FSI,
  },leg2);

  leg->Draw();
  leg2->Draw();

  return canv;
}

TCanvas* STVarMuPTCorrelations_NoNE_Xianguo(std::string STVar, std::string Target,
  bool regen=false, bool DoProfile=true){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas((STVar+Target+"vsmupt_NoNE").c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);


  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    YTitleText =  DoProfile?"#LT#delta#phi_{T}#GT (deg)":
      "{#delta#phi_{T}}_{MLE} (deg)";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else if(STVar=="AlphaT"){
    YTitleText = DoProfile?"#LT#delta#alpha_{T}#GT (deg)":
      "{#delta#alpha_{T}}_{MLE} (deg)";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else if(STVar=="DPT"){
    YTitleText =  DoProfile?"#LT#delta#it{p}_{T}#GT (MeV/#it{c})":
      "{#deltap_{T}}_{MLE} (MeV/#it{c})";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else {
    return new TCanvas("broken");
  }

  TH2* NEUT_C_vsmupt =
    NoSaveClone(FillHistogram2D("NEUT",Target,STVar+"HMPvsMuonPt",regen));
  TH2* GENIE_C_vsmupt =
    NoSaveClone(FillHistogram2D("GENIE",Target,STVar+"HMPvsMuonPt",regen));
  TH2* NuWro_C_vsmupt =
    NoSaveClone(FillHistogram2D("NuWro",Target,STVar+"HMPvsMuonPt",regen));
  // TH2* GiBUU_C_vsmupt =
    // NoSaveClone(FillHistogram2D("GiBUU",Target,STVar+"HMPvsMuonPt",regen));

  TH2* NEUT_C_NoNE_vsmupt =
    NoSaveClone(FillHistogram2D("NEUT",Target,STVar+"HMPvsMuonPt"+"_NoNE",regen));
  TH2* GENIE_C_NoNE_vsmupt =
    NoSaveClone(FillHistogram2D("GENIE",Target,STVar+"HMPvsMuonPt"+"_NoNE",regen));
  TH2* NuWro_C_NoNE_vsmupt =
    NoSaveClone(FillHistogram2D("NuWro",Target,STVar+"HMPvsMuonPt"+"_NoNE",regen));
  // TH2* GiBUU_C_NoNE_vsmupt =
    // NoSaveClone(FillHistogram2D("GiBUU",Target+"_2GeV",STVar+"HMPvsMuonPt",regen));

  TH1* NEUT_C_vsmupt_p;
  TH1* GENIE_C_vsmupt_p;
  TH1* NuWro_C_vsmupt_p;
  // TH1* GiBUU_C_vsmupt_p;

  TH1* NEUT_C_NoNE_vsmupt_p;
  TH1* GENIE_C_NoNE_vsmupt_p;
  TH1* NuWro_C_NoNE_vsmupt_p;
  // TH1* GiBUU_C_NoNE_vsmupt_p;

  if(DoProfile){
    NEUT_C_vsmupt_p = NEUT_C_vsmupt->ProfileX();
    GENIE_C_vsmupt_p = GENIE_C_vsmupt->ProfileX();
    NuWro_C_vsmupt_p = NuWro_C_vsmupt->ProfileX();
    // GiBUU_C_vsmupt_p = GiBUU_C_vsmupt->ProfileX();
    NEUT_C_NoNE_vsmupt_p = NEUT_C_NoNE_vsmupt->ProfileX();
    GENIE_C_NoNE_vsmupt_p = GENIE_C_NoNE_vsmupt->ProfileX();
    NuWro_C_NoNE_vsmupt_p = NuWro_C_NoNE_vsmupt->ProfileX();
    // GiBUU_C_2GeVvsmupt_p = GiBUU_C_2GeVvsmupt->ProfileX();
  } else {
    NEUT_C_vsmupt_p = MakePeakCorrXTH1(NEUT_C_vsmupt);
    GENIE_C_vsmupt_p = MakePeakCorrXTH1(GENIE_C_vsmupt);
    NuWro_C_vsmupt_p = MakePeakCorrXTH1(NuWro_C_vsmupt);
    // GiBUU_C_vsmupt_p = MakePeakCorrXTH1(GiBUU_C_vsmupt);
    NEUT_C_NoNE_vsmupt_p = MakePeakCorrXTH1(NEUT_C_NoNE_vsmupt);
    GENIE_C_NoNE_vsmupt_p = MakePeakCorrXTH1(GENIE_C_NoNE_vsmupt);
    NuWro_C_NoNE_vsmupt_p = MakePeakCorrXTH1(NuWro_C_NoNE_vsmupt);
    // GiBUU_C_NoNE_vsmupt_p = MakePeakCorrXTH1(GiBUU_C_NoNE_vsmupt);
  }

  NEUT_C_vsmupt_p->SetMarkerStyle(0);
  GENIE_C_vsmupt_p->SetMarkerStyle(0);
  NuWro_C_vsmupt_p->SetMarkerStyle(0);
  // GiBUU_C_vsmupt_p->SetMarkerStyle(0);

  NEUT_C_NoNE_vsmupt_p->SetMarkerStyle(0);
  GENIE_C_NoNE_vsmupt_p->SetMarkerStyle(0);
  NuWro_C_NoNE_vsmupt_p->SetMarkerStyle(0);
  // GiBUU_C_vsmupt_p->SetMarkerStyle(0);

  std::vector<TH1*> Hists_v = {
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p,
    // GiBUU_C_vsmupt_p,
    NEUT_C_NoNE_vsmupt_p,
    GENIE_C_NoNE_vsmupt_p,
    NuWro_C_NoNE_vsmupt_p
  };

  std::initializer_list<TH1*> Hists = {
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p,
    // GiBUU_C_vsmupt_p,
    NEUT_C_NoNE_vsmupt_p,
    GENIE_C_NoNE_vsmupt_p,
    NuWro_C_NoNE_vsmupt_p
  };

  NEUT_C_vsmupt_p->SetTitle("Neut");
  GENIE_C_vsmupt_p->SetTitle("GENIE");
  NuWro_C_vsmupt_p->SetTitle("NuWro");

  NEUT_C_NoNE_vsmupt_p->SetTitle("Neut");
  GENIE_C_NoNE_vsmupt_p->SetTitle("GENIE");
  NuWro_C_NoNE_vsmupt_p->SetTitle("NuWro");

  Hists_v[0]->SetLineColor(kBlack);
  Hists_v[1]->SetLineColor(kBlack);
  Hists_v[2]->SetLineColor(kBlack);

  Hists_v[0]->SetLineStyle(kSolid);
  Hists_v[1]->SetLineStyle(kDashed);
  Hists_v[2]->SetLineStyle(kDotted);

  // Hists_v[3]->SetLineWidth(2);
  // Hists_v[4]->SetLineWidth(2);
  // Hists_v[5]->SetLineWidth(2);

  Hists_v[3]->SetLineColor(kBlue);
  Hists_v[4]->SetLineColor(kBlue);
  Hists_v[5]->SetLineColor(kBlue);

  Hists_v[3]->SetLineStyle(kSolid);
  Hists_v[4]->SetLineStyle(kDashed);
  Hists_v[5]->SetLineStyle(kDotted);

  NEUT_C_NoNE_vsmupt_p->GetXaxis()->SetTitle(XTitleText.c_str());
  NEUT_C_NoNE_vsmupt_p->GetYaxis()->SetTitle(YTitleText.c_str());

  SetUserRangeToMaxYVal(Hists);
  SetUserRangeToMaxXValWithContent(Hists);

  for(auto& histo : Hists){
    style::ResetStyle((TH1*)(histo));
  }
  DrawTHs({
    NEUT_C_NoNE_vsmupt_p,
    GENIE_C_NoNE_vsmupt_p,
    NuWro_C_NoNE_vsmupt_p,
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p
  }, "HIST");

  TLegend* leg, *leg2;

  if(STVar=="DeltaPhiT"){

    leg = LegendFactory( (1.0 - P0LegendWidth*2.1),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (1.0 - P0LegendWidth*1.1),
      (1.0 - Pad0YMargin*0.2));

    leg2 = LegendFactory( (1.0 - P0LegendWidth*1.05),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (1.0 - P0LegendWidth*0.05),
      (1.0 - Pad0YMargin*0.2));

  } else if(STVar=="AlphaT"){
    leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
      (Pad0YMargin),
      (1.0 - Pad0XMargin),
      (Pad0YMargin + P0LegendHeight));
  } else if(STVar=="DPT"){
    leg = LegendFactory( (1.0 - P0LegendWidth*2.1),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (1.0 - P0LegendWidth*1.1),
      (1.0 - Pad0YMargin*0.2));

    leg2 = LegendFactory( (1.0 - P0LegendWidth*1.05),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (1.0 - P0LegendWidth*0.05),
      (1.0 - Pad0YMargin*0.2));
  }

  leg->SetHeader("CCQE All");
  leg2->SetHeader("CCQE No N.E.");

  style::ResetStyle(leg);
  style::ResetStyle(leg2);

  AddToLegend({
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p,
    },leg);
  AddToLegend({
    NEUT_C_NoNE_vsmupt_p,
    GENIE_C_NoNE_vsmupt_p,
    NuWro_C_NoNE_vsmupt_p,
    },leg2);

  leg->Draw();
  leg2->Draw();

  return canv;

}

TCanvas* STVarMuPTCorrelations_Single_Xianguo(std::string Gen, std::string STVar, std::string Target,
  bool regen=false, bool DoProfile=true){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas((STVar+Gen+Target+"vsmupt_NoNE").c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);


  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    YTitleText =  DoProfile?"#LT#delta#phi_{T}#GT (deg)":
      "{#delta#phi_{T}}_{MLE} (deg)";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else if(STVar=="AlphaT"){
    YTitleText = DoProfile?"#LT#delta#alpha_{T}#GT (deg)":
      "{#delta#alpha_{T}}_{MLE} (deg)";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else if(STVar=="DPTHMPvsMuonPt"){
    YTitleText =  DoProfile?"#LT#delta#it{p}_{T}#GT (MeV/#it{c})":
      "{#deltap_{T}}_{MLE} (MeV/#it{c})";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else {
    return new TCanvas("broken");
  }

  TH2* Tar_vsmupt =
    NoSaveClone(FillHistogram2D(Gen,Target,STVar,regen));

  TH2* Tar_NoRescat_vsmupt =
    NoSaveClone(FillHistogram2D(Gen,Target+"_norescat",STVar,regen));

  TH2* Tar_NoNE_vsmupt =
    NoSaveClone(FillHistogram2D(Gen,Target,STVar+"_NoNE",regen));

  TH1* Tar_vsmupt_p, *Tar_NoRescat_vsmupt_p, *Tar_NoNE_vsmupt_p;

  if(DoProfile){
    Tar_vsmupt_p = Tar_vsmupt->ProfileX();
    Tar_NoRescat_vsmupt_p = Tar_NoRescat_vsmupt->ProfileX();
    Tar_NoNE_vsmupt_p = Tar_NoNE_vsmupt->ProfileX();
  } else {
    Tar_vsmupt_p = MakePeakCorrXTH1(Tar_vsmupt);
  }

  Tar_vsmupt_p->SetMarkerStyle(0);
  Tar_NoRescat_vsmupt_p->SetMarkerStyle(0);
  Tar_NoNE_vsmupt_p->SetMarkerStyle(0);

  std::vector<TH1*> Hists_v = {
    Tar_vsmupt_p,
    Tar_NoRescat_vsmupt_p,
    Tar_NoNE_vsmupt_p,
  };

  std::initializer_list<TH1*> Hists = {
    Tar_vsmupt_p,
    Tar_NoRescat_vsmupt_p,
    Tar_NoNE_vsmupt_p,
  };

  Tar_vsmupt_p->SetTitle("CCQE All");
  Tar_NoRescat_vsmupt_p->SetTitle("HT Off");
  Tar_NoNE_vsmupt_p->SetTitle("No N.E.");

  Hists_v[0]->SetLineColor(kBlack);
  Hists_v[0]->SetLineStyle(kSolid);

  Hists_v[1]->SetLineColor(kBlack);
  Hists_v[1]->SetLineStyle(kDotted);

  Hists_v[2]->SetLineColor(kBlack);
  Hists_v[2]->SetLineStyle(kDashed);

  Tar_vsmupt_p->GetXaxis()->SetTitle(XTitleText.c_str());
  Tar_vsmupt_p->GetYaxis()->SetTitle(YTitleText.c_str());

  SetUserRangeToMaxYVal(Hists);
  SetUserRangeToMaxXValWithContent(Hists);

  for(auto& histo : Hists){
    style::ResetStyle((TH1*)(histo));
  }

  DrawTHs({
    Tar_vsmupt_p,
    Tar_NoRescat_vsmupt_p,
    Tar_NoNE_vsmupt_p,
  }, "HIST");

  TLegend* leg;

  if(STVar=="DeltaPhiT"){

    leg = LegendFactory( (1.0 - P0LegendWidth*2.1),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (1.0 - P0LegendWidth*1.1),
      (1.0 - Pad0YMargin*0.2));

  } else if(STVar=="AlphaT"){
    leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
      (Pad0YMargin),
      (1.0 - Pad0XMargin),
      (Pad0YMargin + P0LegendHeight));

  } else if(STVar=="DPTHMPvsMuonPt"){
    leg = LegendFactory( (1.0 - P0LegendWidth*2.1),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (1.0 - P0LegendWidth*1.1),
      (1.0 - Pad0YMargin*0.2));
  }

  style::ResetStyle(leg);

  leg->SetHeader("NuWro C-Target 0.6 GeV");

  AddToLegend({
    Tar_vsmupt_p,
    Tar_NoRescat_vsmupt_p,
    Tar_NoNE_vsmupt_p,
  },leg);
  leg->Draw();

  return canv;

}

TCanvas* STVarMuPTCorrelations_EDep_Xianguo(std::string STVar, std::string Target,
   bool regen=false, bool DoProfile=true){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  TCanvas* canv = new TCanvas((STVar+Target+"vsmupt_NoNE").c_str(),"");
  style::PadSetup(canv);
  canv->SetRightMargin(0.03);
  canv->SetTopMargin(0.02);
  canv->SetLeftMargin(0.12);
  canv->SetBottomMargin(0.15);


  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    YTitleText =  DoProfile?"#LT#delta#phi_{T}#GT (deg)":
      "{#delta#phi_{T}}_{MLE} (deg)";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else if(STVar=="AlphaT"){
    YTitleText = DoProfile?"#LT#delta#alpha_{T}#GT (deg)":
      "{#delta#alpha_{T}}_{MLE} (deg)";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else if(STVar=="DPT"){
    YTitleText =  DoProfile?"#LT#delta#it{p}_{T}#GT (MeV/#it{c})":
      "{#deltap_{T}}_{MLE} (MeV/#it{c})";
    XTitleText = "#it{p}^{#mu}_{T} (MeV/#it{c})";
  } else {
    return new TCanvas("broken");
  }

  TH2* NEUT_C_vsmupt =
    NoSaveClone(FillHistogram2D("NEUT",Target,STVar+"HMPvsMuonPt",regen));
  TH2* GENIE_C_vsmupt =
    NoSaveClone(FillHistogram2D("GENIE",Target,STVar+"HMPvsMuonPt",regen));
  TH2* NuWro_C_vsmupt =
    NoSaveClone(FillHistogram2D("NuWro",Target,STVar+"HMPvsMuonPt",regen));
  // TH2* GiBUU_C_vsmupt =
    // NoSaveClone(FillHistogram2D("GiBUU",Target,STVar+"HMPvsMuonPt",regen));

  TH2* NEUT_C_2GeV_vsmupt =
    NoSaveClone(FillHistogram2D("NEUT",Target+"_2GeV",STVar+"HMPvsMuonPt",regen));
  // TH2* GENIE_C_2GeV_vsmupt =
  //   NoSaveClone(FillHistogram2D("GENIE",Target+"_2GeV",STVar+"HMPvsMuonPt",regen));
  TH2* NuWro_C_2GeV_vsmupt =
    NoSaveClone(FillHistogram2D("NuWro",Target+"_2GeV",STVar+"HMPvsMuonPt",regen));
  // TH2* GiBUU_C_2GeV_vsmupt =
    // NoSaveClone(FillHistogram2D("GiBUU",Target+"_2GeV",STVar+"HMPvsMuonPt",regen));

  TH2* NEUT_C_8GeV_vsmupt =
    NoSaveClone(FillHistogram2D("NEUT",Target+"_8GeV",STVar+"HMPvsMuonPt",regen));
  // TH2* GENIE_C_8GeV_vsmupt =
  //   NoSaveClone(FillHistogram2D("GENIE",Target+"_8GeV",STVar+"HMPvsMuonPt",regen));
  TH2* NuWro_C_8GeV_vsmupt =
    NoSaveClone(FillHistogram2D("NuWro",Target+"_8GeV",STVar+"HMPvsMuonPt",regen));

  TH1* NEUT_C_vsmupt_p;
  TH1* GENIE_C_vsmupt_p;
  TH1* NuWro_C_vsmupt_p;
  // TH1* GiBUU_C_vsmupt_p;

  TH1* NEUT_C_2GeV_vsmupt_p;
  // TH1* GENIE_C_2GeV_vsmupt_p;
  TH1* NuWro_C_2GeV_vsmupt_p;
  // TH1* GiBUU_C_2GeV_vsmupt_p;

  TH1* NEUT_C_8GeV_vsmupt_p;
  // TH1* GENIE_C_8GeV_vsmupt_p;
  TH1* NuWro_C_8GeV_vsmupt_p;
  // TH1* GiBUU_C_8GeV_vsmupt_p;

  if(DoProfile){
    NEUT_C_vsmupt_p = NEUT_C_vsmupt->ProfileX();
    GENIE_C_vsmupt_p = GENIE_C_vsmupt->ProfileX();
    NuWro_C_vsmupt_p = NuWro_C_vsmupt->ProfileX();
    // GiBUU_C_vsmupt_p = GiBUU_C_vsmupt->ProfileX();
    NEUT_C_2GeV_vsmupt_p = NEUT_C_2GeV_vsmupt->ProfileX();
    // GENIE_C_2GeV_vsmupt_p = GENIE_C_2GeV_vsmupt->ProfileX();
    NuWro_C_2GeV_vsmupt_p = NuWro_C_2GeV_vsmupt->ProfileX();
    // GiBUU_C_2GeVvsmupt_p = GiBUU_C_2GeVvsmupt->ProfileX();

    NEUT_C_8GeV_vsmupt_p = NEUT_C_8GeV_vsmupt->ProfileX();
    // GENIE_C_8GeV_vsmupt_p = GENIE_C_8GeV_vsmupt->ProfileX();
    NuWro_C_8GeV_vsmupt_p = NuWro_C_8GeV_vsmupt->ProfileX();
    // GiBUU_C_8GeVvsmupt_p = GiBUU_C_8GeVvsmupt->ProfileX();

  } else {
    NEUT_C_vsmupt_p = MakePeakCorrXTH1(NEUT_C_vsmupt);
    GENIE_C_vsmupt_p = MakePeakCorrXTH1(GENIE_C_vsmupt);
    NuWro_C_vsmupt_p = MakePeakCorrXTH1(NuWro_C_vsmupt);
    // GiBUU_C_vsmupt_p = MakePeakCorrXTH1(GiBUU_C_vsmupt);
    NEUT_C_2GeV_vsmupt_p = MakePeakCorrXTH1(NEUT_C_2GeV_vsmupt);
    // GENIE_C_2GeV_vsmupt_p = MakePeakCorrXTH1(GENIE_C_2GeV_vsmupt);
    NuWro_C_2GeV_vsmupt_p = MakePeakCorrXTH1(NuWro_C_2GeV_vsmupt);
    // GiBUU_C_2GeV_vsmupt_p = MakePeakCorrXTH1(GiBUU_C_2GeV_vsmupt);
  }

  NEUT_C_vsmupt_p->SetMarkerStyle(0);
  GENIE_C_vsmupt_p->SetMarkerStyle(0);
  NuWro_C_vsmupt_p->SetMarkerStyle(0);
  // GiBUU_C_vsmupt_p->SetMarkerStyle(0);

  NEUT_C_2GeV_vsmupt_p->SetMarkerStyle(0);
  // GENIE_C_2GeV_vsmupt_p->SetMarkerStyle(0);
  NuWro_C_2GeV_vsmupt_p->SetMarkerStyle(0);
  // GiBUU_C_vsmupt_p->SetMarkerStyle(0);

  NEUT_C_8GeV_vsmupt_p->SetMarkerStyle(0);
  // GENIE_C_8GeV_vsmupt_p->SetMarkerStyle(0);
  NuWro_C_8GeV_vsmupt_p->SetMarkerStyle(0);
  // GiBUU_C_vsmupt_p->SetMarkerStyle(0);

  std::vector<TH1*> Hists_v = {
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p,
    // GiBUU_C_vsmupt_p,
    NEUT_C_2GeV_vsmupt_p, 0,
    // GENIE_C_2GeV_vsmupt_p,
    NuWro_C_2GeV_vsmupt_p,
    NEUT_C_8GeV_vsmupt_p, 0,
    // GENIE_C_8GeV_vsmupt_p,
    NuWro_C_8GeV_vsmupt_p
  };

  std::initializer_list<TH1*> Hists = {
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p,
    // GiBUU_C_vsmupt_p,
    NEUT_C_2GeV_vsmupt_p,
    // GENIE_C_2GeV_vsmupt_p,
    NuWro_C_2GeV_vsmupt_p,
    NEUT_C_8GeV_vsmupt_p,
    NuWro_C_8GeV_vsmupt_p,
    // GENIE_C_8GeV_vsmupt_p,

  };

  NEUT_C_vsmupt_p->SetTitle("Neut");
  GENIE_C_vsmupt_p->SetTitle("GENIE");
  NuWro_C_vsmupt_p->SetTitle("NuWro");
  NEUT_C_2GeV_vsmupt_p->SetTitle("Neut");
  // GENIE_C_2GeV_vsmupt_p->SetTitle("GENIE");
  NuWro_C_2GeV_vsmupt_p->SetTitle("NuWro");
  NEUT_C_8GeV_vsmupt_p->SetTitle("Neut");
  // GENIE_C_8GeV_vsmupt_p->SetTitle("GENIE");
  NuWro_C_8GeV_vsmupt_p->SetTitle("NuWro");

  Hists_v[0]->SetLineColor(kBlack);
  Hists_v[1]->SetLineColor(kBlack);
  Hists_v[2]->SetLineColor(kBlack);

  Hists_v[0]->SetLineStyle(kSolid);
  Hists_v[1]->SetLineStyle(kDashed);
  Hists_v[2]->SetLineStyle(kDotted);

  // Hists_v[3]->SetLineWidth(2);
  // Hists_v[4]->SetLineWidth(2);
  // Hists_v[5]->SetLineWidth(2);

  Hists_v[3]->SetLineColor(kRed);
  // Hists_v[4]->SetLineColor(kRed);
  Hists_v[5]->SetLineColor(kRed);

  Hists_v[3]->SetLineStyle(kSolid);
  // Hists_v[4]->SetLineStyle(kDashed);
  Hists_v[5]->SetLineStyle(kDotted);

  Hists_v[6]->SetLineColor(kGreen);
  // Hists_v[7]->SetLineColor(kGreen);
  Hists_v[8]->SetLineColor(kGreen);

  Hists_v[6]->SetLineStyle(kSolid);
  // Hists_v[7]->SetLineStyle(kDashed);
  Hists_v[8]->SetLineStyle(kDotted);

  NEUT_C_2GeV_vsmupt_p->GetXaxis()->SetTitle(XTitleText.c_str());
  NEUT_C_2GeV_vsmupt_p->GetYaxis()->SetTitle(YTitleText.c_str());

  SetUserRangeToMaxYVal(Hists);
  SetUserRangeToMaxXValWithContent(Hists);

  for(auto& histo : Hists){
    style::ResetStyle((TH1*)(histo));
  }
  DrawTHs({
    NEUT_C_8GeV_vsmupt_p,
    // GENIE_C_8GeV_vsmupt_p,
    NuWro_C_8GeV_vsmupt_p,
    NEUT_C_2GeV_vsmupt_p,
    // GENIE_C_2GeV_vsmupt_p,
    NuWro_C_2GeV_vsmupt_p,
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p
  }, "HIST");

  TLegend* leg, *leg2;

  if(STVar=="DeltaPhiT"){

    leg = LegendFactory( (0.12 + P0LegendWidth*0.1),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (0.12 + P0LegendWidth*1.1),
      (1.0 - Pad0YMargin*0.2));

    leg2 = LegendFactory( (0.12 + P0LegendWidth*1.2),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (0.12 + P0LegendWidth*2.2),
      (1.0 - Pad0YMargin*0.2));

  } else if(STVar=="AlphaT"){

    leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
      (Pad0YMargin),
      (1.0 - Pad0XMargin),
      (Pad0YMargin + P0LegendHeight));

  } else if(STVar=="DPT"){
    leg = LegendFactory( (0.12 + P0LegendWidth*0.1),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (0.12 + P0LegendWidth*1.1),
      (1.0 - Pad0YMargin*0.2));

    leg2 = LegendFactory( (0.12 + P0LegendWidth*1.2),
      (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
      (0.12 + P0LegendWidth*2.2),
      (1.0 - Pad0YMargin*0.2));
  }

  leg->SetHeader("0.6 GeV");
  leg2->SetHeader("2.0 GeV");

  style::ResetStyle(leg);
  style::ResetStyle(leg2);

  AddToLegend({
    NEUT_C_vsmupt_p,
    GENIE_C_vsmupt_p,
    NuWro_C_vsmupt_p,
    },leg);
  AddToLegend({
    NEUT_C_2GeV_vsmupt_p,
    // GENIE_C_2GeV_vsmupt_p,
    NuWro_C_2GeV_vsmupt_p,
    },leg2);

  leg->Draw();
  leg2->Draw();

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

TCanvas* TargetDependenceRatio(std::string mode, std::string other){
  TCanvas* canv = CanvasFactory("tdep"+mode+other);

  std::string XTitleText = "";
  if(mode=="CCQE"){
    XTitleText = DPhiTTitle;
  } else if(mode=="CCQEAlphaT"){
    XTitleText = DAlphaTTitle;
  } else if(mode=="CCQEDPt"){
    XTitleText = DPTTitle;
  }  else if(mode=="CCQEMuPt"){
    XTitleText = MuPtTitle;
  }else {
    return new TCanvas("broken");
  }

  TH1* GENIE_C_Hist = NoSaveClone(FillHistogram1D("GENIE","C",mode));
  TH1* GENIE_Pb_Hist = NoSaveClone(FillHistogram1D("GENIE",other,mode));
  TH1* NEUT_C_Hist = NoSaveClone(FillHistogram1D("NEUT","C_SFOff",mode));
  TH1* NEUT_Pb_Hist = NoSaveClone(FillHistogram1D("NEUT",other,mode));
  TH1* NuWro_C_Hist = NoSaveClone(FillHistogram1D("NuWro","C",mode));
  TH1* NuWro_Pb_Hist = NoSaveClone(FillHistogram1D("NuWro",other,mode));
  TH1* GiBUU_C_Hist = NoSaveClone(FillHistogram1D("GiBUU","C",mode));
  TH1* GiBUU_Pb_Hist = NoSaveClone(FillHistogram1D("GiBUU",other,mode));

  MakePDF({
    NEUT_C_Hist,
    NEUT_Pb_Hist,
    GENIE_C_Hist,
    GENIE_Pb_Hist,
    NuWro_C_Hist,
    NuWro_Pb_Hist,
    GiBUU_C_Hist,
    GiBUU_Pb_Hist}, true);

  TH1* GENIE_Ratio = NoSaveClone(GENIE_Pb_Hist);
  GENIE_Ratio->Divide(GENIE_Pb_Hist, GENIE_C_Hist,1,1,"B");
  TH1* NEUT_Ratio = NoSaveClone(NEUT_Pb_Hist);
  NEUT_Ratio->Divide(NEUT_Pb_Hist, NEUT_C_Hist,1,1,"B");
  TH1* NuWro_Ratio = NoSaveClone(NuWro_Pb_Hist);
  NuWro_Ratio->Divide(NuWro_Pb_Hist, NuWro_C_Hist,1,1,"B");
  TH1* GiBUU_Ratio = NoSaveClone(GiBUU_Pb_Hist);
  GiBUU_Ratio->Divide(GiBUU_Pb_Hist, GiBUU_C_Hist,1,1,"B");

  SetOnePaneHistDefaults(NEUT_Ratio, "Ratio", XTitleText.c_str());

  std::initializer_list<TH1*> Ratio_Hists = {NEUT_Ratio,
                                            GENIE_Ratio,
                                            NuWro_Ratio,
                                            GiBUU_Ratio};

  GENIE_Ratio->SetTitle(("GENIE (hN INC) " + other + "/C").c_str());
  NEUT_Ratio->SetTitle(("NEUT " + other + "/C (RFG)").c_str());
  NuWro_Ratio->SetTitle(("NuWro (LFG) " + other + "/C").c_str());
  GiBUU_Ratio->SetTitle(("GiBUU (LFG) " + other + "/C").c_str());

  TH1* DRG = DummyRatioGuide(kBlack,0,180);

  SetStyles(Ratio_Hists);

  TLegend* leg2 = LegendFactory( (Pad0XMargin),
  (1.0 - P0LegendHeight - Pad0YMargin),
  (Pad1XMargin + P0LegendWidth*1.1),
  (1.0 - Pad0YMargin));

  AddToLegend(Ratio_Hists,leg2);
  SetOnePaneHistDefaults(GENIE_Ratio, "Ratio",
    DPhiTTitle);
  SetUserRangeToMaxYVal(Ratio_Hists);
  DrawTHs(Ratio_Hists);
  DRG->Draw("HIST SAME");
  leg2->Draw();

  TitleFactory("Target Dependence "+mode+other)->Draw();
  return canv;
}

TCanvas* GenVarComparison(std::string mode, std::string other, bool regen=false){
  TCanvas* canv = CanvasFactory("tdep"+mode+other);

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(mode=="CCQE" || mode == "CC0Pi1P"){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(mode=="CCQEAlphaT" || mode == "CC0Pi1PAlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(mode=="CCQEDPt" || mode == "CC0Pi1PDPt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  } else if(mode=="CCQEStruckNucleonMomentum"){
    XTitleText = StruckNucleonMomTitle;
    YTitleText = StruckNucleonMom_PDF_Title;
  } else {
    return new TCanvas("broken");
  }
  TLegend* leg2;
  if(mode=="CCQE"){
    leg2 = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else if(mode=="CCQEAlphaT"){
    leg2 = LegendFactory( (Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (Pad0XMargin + P0LegendWidth*1.1),
    (1.0 - Pad0YMargin));
  } else if(mode=="CCQEDPt"){
    leg2 = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else{
    leg2 = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  }

  TH1* GENIE_Hist = NoSaveClone(FillHistogram1D("GENIE",other,mode, regen));
  TH1* NEUT_Hist = NoSaveClone(FillHistogram1D("NEUT",other,mode, regen));
  TH1* NuWro_Hist = NoSaveClone(FillHistogram1D("NuWro",other,mode, regen));
  TH1* GENIE_hA_Hist = NoSaveClone(FillHistogram1D("GENIE",other+"_hA",mode, regen));
  TH1* GiBUU_Hist = NoSaveClone(FillHistogram1D("GiBUU",other,mode, regen));

  std::initializer_list<TH1*> Hists = {NEUT_Hist,
            GENIE_Hist,
            NuWro_Hist,
            GENIE_hA_Hist,
            GiBUU_Hist};

  MakePDF(Hists);

  SetStyles(Hists);

  SetOnePaneHistDefaults((*Hists.begin()), YTitleText.c_str(), XTitleText.c_str());

  GENIE_Hist->SetTitle(("GENIE (hN)"));
  NEUT_Hist->SetTitle(("NEUT"));
  NuWro_Hist->SetTitle(("NuWro"));
  GENIE_hA_Hist->SetTitle("GENIE (hA)");
  GiBUU_Hist->SetTitle(("GiBUU"));

  AddToLegend(Hists,leg2);


  SetUserRangeToMaxYVal(Hists);
  DrawTHs(Hists);
  leg2->Draw();

  TitleFactory(std::string("Target Dependence ")+mode+other)->Draw();
  return canv;
}

TCanvas* CArPbTargetComparison(std::string Generator, std::string STVar,
  bool regen=false){

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar==""){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(STVar=="DPt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  } else if(STVar=="StruckNucleonMomentum"){
    XTitleText = StruckNucleonMomTitle;
    YTitleText = StruckNucleonMom_PDF_Title;
  } else {
    return new TCanvas("broken");
  }

  TCanvas* canv = CanvasFactory("TDep_CCQE"+STVar+Generator);

  TH1* C_Hist = NoSaveClone(FillHistogram1D(Generator,"C"+std::string(Generator=="NEUT"?"_SFOff":""),"CCQE"+STVar,regen));
  // TH1* C_Silly = NoSaveClone(FillHistogram1D(Generator,"C_silly","CCQE"+STVar));
  TH1* Ar_Hist = NoSaveClone(FillHistogram1D(Generator,"Ar","CCQE"+STVar,regen));
  TH1* Pb_Hist = NoSaveClone(FillHistogram1D(Generator,"Pb","CCQE"+STVar,regen));

  AddCountToTHTitle({std::make_pair("C",C_Hist),
    // std::make_pair("C (Eb 100MeV, kf 300MeV)",C_Silly),
    std::make_pair("Ar",Ar_Hist),
    std::make_pair("Pb",Pb_Hist)});

  MakePDF({C_Hist,
    Ar_Hist,//C_Silly,
    Pb_Hist});

  SetStyles({C_Hist,
    Ar_Hist,//C_Silly,
    Pb_Hist});

  SetOnePaneHistDefaults(C_Hist, YTitleText.c_str(), XTitleText.c_str());

  SetUserRangeToMaxYVal({C_Hist,//C_Silly,
      Ar_Hist,
      Pb_Hist});
  DrawTHs({C_Hist,//C_Silly,
      Ar_Hist,
      Pb_Hist});

  TLegend* leg2;
  if(STVar==""){
    leg2 = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else if(STVar=="AlphaT"){
    leg2 = LegendFactory( (Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (Pad0XMargin + P0LegendWidth*1.1),
    (1.0 - Pad0YMargin));
  } else if(STVar=="DPt"){
    leg2 = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else{
    leg2 = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  }

  AddToLegend({C_Hist,//C_Silly,
      Ar_Hist,
      Pb_Hist},leg2);

  leg2->Draw();
  TitleFactory(Generator + " Target Dependence "
    + (STVar.length()? STVar : std::string("#delta#phi_{T}") ) )->Draw();
  return canv;
}

TCanvas* DrawDistrib(std::string Generator, std::string Target,
  std::string STVar="", std::string selection="", std::string Xtitle="",
  std::string Ytitle="", bool doPDF=false, bool redraw=false){

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar==""){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(STVar=="DPt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  } else {
    XTitleText = Xtitle;
    YTitleText = Ytitle;
  }

  TCanvas* canv = CanvasFactory(Generator+Target+"FSIMODE"+STVar);
  TH1* CCQEAll = NoSaveClone(FillHistogram1D(Generator,Target, selection,redraw));

  if(doPDF){
    MakePDF(CCQEAll);
  }

  CCQEAll->SetLineStyle(1);
  CCQEAll->SetLineWidth(2);
  CCQEAll->SetMarkerStyle(0);

  SetOnePaneHistDefaults(CCQEAll,YTitleText.c_str(),
    XTitleText.c_str());

  DrawTHs({CCQEAll});

  TitleFactory((std::string(Generator+ " " + Target
    + " " + selection + " " + STVar)).c_str())->Draw();
  return canv;
}

TCanvas* SVTHistCompare(std::initializer_list<TH1*> hists,
  std::string Var="", std::string canvname="SVTHistCompare",
  std::string Title=""){

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(Var=="CCQE"){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(Var=="CCQEAlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(Var=="CCQEDPt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  } else if(Var=="CCQEStruckNucleonMomentum"){
    XTitleText = StruckNucleonMomTitle;
    YTitleText = StruckNucleonMom_PDF_Title;
  } else if(Var=="CCQEStruckNucleon_BE"){
    XTitleText = StruckNucleonBETitle;
    YTitleText = StruckNucleonBE_PDF_Title;
  } else {
    return new TCanvas("broken");
  }

  TCanvas* canv = CanvasFactory(canvname);

  SetOnePaneHistDefaults((*hists.begin()),YTitleText.c_str(), XTitleText.c_str());
  MakePDF(hists);
  SetUserRangeToMaxYVal(hists);
  SetStyles(hists,1);
  DrawTHs(hists);

  TLegend* leg;
  if(Var=="CCQE"){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*2.8),
    (1.0 - P0LegendHeight*2.0 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else if(Var=="CCQEAlphaT"){
    leg = LegendFactory( (Pad0XMargin),
    (1.0 - P0LegendHeight*1.5 - Pad0YMargin),
    (Pad0XMargin + P0LegendWidth*2.1),
    (1.0 - Pad0YMargin));
  } else if(Var=="CCQEDPt"){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*2.0),
    (1.0 - P0LegendHeight*1.3 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else if(Var == "CCQEStruckNucleonMomentum" || Var == "CCQEStruckNucleon_BE" ){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*1.6),
    (1.0 - P0LegendHeight*1.3 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  }

  AddToLegend(hists,leg);
  leg->Draw();

  TitleFactory(Title.c_str())->Draw();
  return canv;
}

TCanvas* SVTFMModelCompare(std::string Var="", bool Kaskada_on=false){
  std::string pfx = Kaskada_on?"":"_norescat";
  // C.LFG
  // C.E_SF
  // C.E_SF_M1
  // C.FG
  // C.B_R
  // C_RPA1


  //Want
  //-LFG
  //-SF+MEC
  //RFG+RPA+MEC
  //B_R

  bool RegenHisto = false;//(Var=="CCQEStruckNucleon_BE");

  TH1* CCQE_NUWRO_C_LFG =
    NoSaveClone(FillHistogram1D("NuWro","C.LFG"+pfx,Var,RegenHisto));
  TH1* CCQE_NUWRO_C_E_SF_M1 =
    NoSaveClone(FillHistogram1D("NuWro","C.E_SF_M1"+pfx,Var,RegenHisto));
  TH1* CCQE_NUWRO_C_FG =
    NoSaveClone(FillHistogram1D("NuWro","C.FG"+pfx,Var,RegenHisto));
  TH1* CCQE_NUWRO_C_B_R =
    NoSaveClone(FillHistogram1D("NuWro","C.B_R"+pfx,Var,RegenHisto));
  TH1* CCQE_NUWRO_C_RPA1 =
    NoSaveClone(FillHistogram1D("NuWro","C_RPA1"+pfx,Var,RegenHisto));

  TH1* NEUT_C_SFOff =
    NoSaveClone(FillHistogram1D("NEUT","C_SFOff"+pfx,Var,RegenHisto));
  TH1* NEUT_C_DEFAULT =
    NoSaveClone(FillHistogram1D("NEUT","C"+pfx,Var,RegenHisto));

  TH1* NEUT_TE_dipole =
    NoSaveClone(FillHistogram1D("NEUT","C_TE_smithmoniz"+pfx,Var,RegenHisto));

  if((Var=="CCQEStruckNucleon_BE")){
    TH1* CCQE_NUWRO_C_RFG_MEC =
      NoSaveClone(FillHistogram1D("NuWro","C.FG"+pfx,"CCQE_MEC_StruckNucleon_BE",RegenHisto));

  AddCountToTHTitle({std::make_pair("NuWro LFG",CCQE_NUWRO_C_LFG),
                      std::make_pair("NuWro Eff SF M1",CCQE_NUWRO_C_E_SF_M1),
                      std::make_pair("NuWro RFG",CCQE_NUWRO_C_FG),
                      std::make_pair("NuWro RFG+MEC",CCQE_NUWRO_C_RFG_MEC),
                      std::make_pair("NuWro Bodek-Ritchie",CCQE_NUWRO_C_B_R),
                      std::make_pair("LFG + RPA M1",CCQE_NUWRO_C_RPA1),
                      std::make_pair("NEUT SF Off",NEUT_C_SFOff),
                      std::make_pair("NEUT Dipole+TE",NEUT_TE_dipole),
                      std::make_pair("NEUT ESF",NEUT_C_DEFAULT)});

  return SVTHistCompare({CCQE_NUWRO_C_LFG,
                        CCQE_NUWRO_C_E_SF_M1,
                        CCQE_NUWRO_C_FG,
                        CCQE_NUWRO_C_RFG_MEC,
                        CCQE_NUWRO_C_B_R,
                        CCQE_NUWRO_C_RPA1,
                        NEUT_C_SFOff,
                        NEUT_TE_dipole,
                        NEUT_C_DEFAULT}, Var, Var+"SVTFMModelCompare_NuWro",
                        "NuWro Nuclear Model Compare "+ Var + pfx);
  } else {
    AddCountToTHTitle({std::make_pair("NuWro LFG",CCQE_NUWRO_C_LFG),
                        std::make_pair("NuWro Eff SF M1",CCQE_NUWRO_C_E_SF_M1),
                        std::make_pair("NuWro RFG",CCQE_NUWRO_C_FG),
                        std::make_pair("NuWro Bodek-Ritchie",CCQE_NUWRO_C_B_R),
                        std::make_pair("NuWro LFG + RPA M1",CCQE_NUWRO_C_RPA1),
                        std::make_pair("NEUT SF Off",NEUT_C_SFOff),
                        std::make_pair("NEUT Dipole+TE",NEUT_TE_dipole),
                        std::make_pair("NEUT ESF",NEUT_C_DEFAULT)});

    return SVTHistCompare({CCQE_NUWRO_C_LFG,
                          CCQE_NUWRO_C_E_SF_M1,
                          CCQE_NUWRO_C_FG,
                          CCQE_NUWRO_C_B_R,
                          CCQE_NUWRO_C_RPA1,
                          NEUT_C_SFOff,
                          NEUT_TE_dipole,
                          NEUT_C_DEFAULT}, Var, Var+"SVTFMModelCompare_NuWro",
                          "Nuclear Model Compare "+ Var + pfx);
  }
}

TCanvas* SVTFancyMModelCompare(std::string Var="", bool Kaskada_on=false){
  std::string pfx = Kaskada_on?"":"_norescat";

  //Want
  //-LFG
  //-SF+MEC
  //RFG+RPA+MEC
  //B_R

  bool RegenHisto = false;//(Var=="CCQEStruckNucleon_BE");

  TH1* CCQE_NUWRO_C_LFG =
    NoSaveClone(FillHistogram1D("NuWro","C.LFG"+pfx,Var,RegenHisto));
  TH1* CCQE_NUWRO_C_E_SF_M1 =
    NoSaveClone(FillHistogram1D("NuWro","C.E_SF_M1"+pfx,Var+"_MEC",RegenHisto));
  TH1* CCQE_NUWRO_C_FGRPA =
    NoSaveClone(FillHistogram1D("NuWro","C.RFGRPA"+pfx,Var+"_MEC",RegenHisto));
  TH1* CCQE_NUWRO_C_B_R =
    NoSaveClone(FillHistogram1D("NuWro","C.B_R"+pfx,Var,RegenHisto));

  TH1* NEUT_C_DEFAULT =
    NoSaveClone(FillHistogram1D("NEUT","C"+pfx,Var+"_MEC",RegenHisto));

  TH1* NEUT_C_RFGRPAMEC =
    NoSaveClone(FillHistogram1D("NEUT","C_RFGRPA"+pfx,Var+"_MEC",RegenHisto));

  TH1* GENIE_C_DEFAULT =
    NoSaveClone(FillHistogram1D("GENIE","C"+pfx,Var,RegenHisto));

  // TH1* GiBUU_C_DEFAULT =
    // NoSaveClone(FillHistogram1D("GENIE","C"+pfx,Var,RegenHisto));


  AddCountToTHTitle({std::make_pair("NuWro LFG",CCQE_NUWRO_C_LFG),
                      std::make_pair("NuWro SF+MEC",CCQE_NUWRO_C_E_SF_M1),
                      std::make_pair("NuWro RFG+RPA+MEC",CCQE_NUWRO_C_FGRPA),
                      std::make_pair("NuWro Bodek-Ritchie",CCQE_NUWRO_C_B_R),
                      std::make_pair("NEUT SF+MEC",NEUT_C_DEFAULT),
                      std::make_pair("NEUT RFG+RPA+MEC",NEUT_C_RFGRPAMEC),
                      std::make_pair("GENIE Bodek-Ritchie (hN INC)",GENIE_C_DEFAULT),
                    });

  return SVTHistCompare({CCQE_NUWRO_C_LFG,
                        CCQE_NUWRO_C_E_SF_M1,
                        CCQE_NUWRO_C_FGRPA,
                        CCQE_NUWRO_C_B_R,
                        NEUT_C_DEFAULT,
                        NEUT_C_RFGRPAMEC,
                        GENIE_C_DEFAULT,
                      }, Var, Var+"SVTFMModelCompare_NuWro_Fancy",
                        "Nuclear Model Compare "+ Var + pfx);
}

TCanvas* NuWroModel_DPTVsDAlphaT(bool Kaskada_on=false){
  std::string pfx = Kaskada_on?"":"_norescat";
  TCanvas* mycanv = CanvasFactory("blabla");
  mycanv->Divide(2,2,0.011,0.011);

  TH2* NuWro_C_DPTALPHA = NoSaveClone(FillHistogram2D("NuWro","C"+pfx,"DPTHMPvsDAlphaT"));
  TH2* NuWro_C_B_R_DPTALPHA = NoSaveClone(FillHistogram2D("NuWro","C.B_R"+pfx,"DPTHMPvsDAlphaT"));
  TH2* NuWro_C_E_SF_DPTALPHA = NoSaveClone(FillHistogram2D("NuWro","C.E_SF_M1"+pfx,"DPTHMPvsDAlphaT"));
  TH2* NuWro_C_FG_DPTALPHA = NoSaveClone(FillHistogram2D("NuWro","C.FG"+pfx,"DPTHMPvsDAlphaT"));

  SetOnePaneHistDefaults(NuWro_C_DPTALPHA,DPTTitle,DAlphaTTitle);
  SetOnePaneHistDefaults(NuWro_C_B_R_DPTALPHA,DPTTitle,DAlphaTTitle);
  SetOnePaneHistDefaults(NuWro_C_E_SF_DPTALPHA,DPTTitle,DAlphaTTitle);
  SetOnePaneHistDefaults(NuWro_C_FG_DPTALPHA,DPTTitle,DAlphaTTitle);

  NuWro_C_DPTALPHA->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);
  NuWro_C_B_R_DPTALPHA->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);
  NuWro_C_E_SF_DPTALPHA->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);
  NuWro_C_FG_DPTALPHA->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);

  NuWro_C_DPTALPHA->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);
  NuWro_C_B_R_DPTALPHA->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);
  NuWro_C_E_SF_DPTALPHA->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);
  NuWro_C_FG_DPTALPHA->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);

  mycanv->cd(1);
  NuWro_C_DPTALPHA->Draw("colz");
  TitleFactory("#scale[2]{#color[2]{NuWro_C_DPTALPHA}}")->Draw();
  mycanv->cd(2);
  NuWro_C_B_R_DPTALPHA->Draw("colz");
  TitleFactory("#scale[2]{#color[2]{NuWro_C_B_R_DPTALPHA}}")->Draw();
  mycanv->cd(3);
  NuWro_C_E_SF_DPTALPHA->Draw("colz");
  TitleFactory("#scale[2]{#color[2]{NuWro_C_E_SF_DPTALPHA}}")->Draw();
  mycanv->cd(4);
  NuWro_C_FG_DPTALPHA->Draw("colz");
  TitleFactory("#scale[2]{#color[2]{NuWro_C_FG_DPTALPHA}}")->Draw();

  mycanv->cd();

  return mycanv;
}

TCanvas* DPhiTVsDAlphaT(std::string gen, std::string tar){
  TCanvas* mycanv = CanvasFactory("DPhiTVsDAlphaT" + gen + tar);

  TH2* DPhiVsDalphaT_Hist = NoSaveClone(FillHistogram2D(gen,tar,"DPhiTHMPvsDAlphaT"));

  SetOnePaneHistDefaults(DPhiVsDalphaT_Hist,DPhiTTitle,DAlphaTTitle);

  DPhiVsDalphaT_Hist->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);

  mycanv->SetLogz(true);

  DPhiVsDalphaT_Hist->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);
  DPhiVsDalphaT_Hist->Draw("colz");
  TitleFactory("DPhiTVsDAlphaT" + gen + tar)->Draw();

  return mycanv;
}

TCanvas* ERecDistribComp(std::string Generator, std::string Target,
  double minx=500, double maxx=700){

  TCanvas* canv = CanvasFactory(Generator+Target+"ERecDistribComp");

  TH1* ERec_CCQE =
    FillHistogram1D(Generator,Target,"CCQENuERec");
  TH1* ERec_MuP =
    FillHistogram1D(Generator,Target,"CCQENuMuPERec");
  TH1* ETrue =
    FillHistogram1D(Generator,Target,"CCQENuETrue");

  SetOnePaneHistDefaults(ERec_CCQE,"Count","E_{#nu} [MeV]");

  ERec_CCQE->SetTitle("CCQE Formula E_{rec}");
  ERec_MuP->SetTitle("E_{rec} = E_{#mu} + E_{p} - m_n + E_b");
  ETrue->SetTitle("E_{true}");

  SetStyles({ERec_CCQE,
    ERec_MuP,
    ETrue});

  ERec_CCQE->GetXaxis()->SetRangeUser(minx,maxx);

  SetUserRangeToMaxYVal({ERec_CCQE,
    ERec_MuP});

  TLegend* leg = LegendFactory( (Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (Pad0XMargin + P0LegendWidth*2.2),
    (1.0 - Pad0YMargin));
  AddToLegend({ERec_CCQE,
    ERec_MuP,
    ETrue}, leg);

  DrawTHs({ERec_CCQE,
    ERec_MuP,
    ETrue});

  leg->Draw();
  TitleFactory(Generator + " " + Target + " ERec Comparison" )->Draw();
  return canv;
}

TCanvas* pHMP_pMu_FS_2D(std::string gen, std::string tar){
// CCQE_fsmmtm
  TCanvas* mycanv = CanvasFactory("pHMP_pMu_FS_2D" + gen + tar);

  TH2* pHMP_pMu_FS_2D_hist = NoSaveClone(FillHistogram2D(gen,tar,"CCQE_fsmmtm"));

  SetOnePaneHistDefaults(pHMP_pMu_FS_2D_hist,"Proton Momenta [MeV/C]","Muon Momenta [MeV/C]");

  pHMP_pMu_FS_2D_hist->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);

  pHMP_pMu_FS_2D_hist->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);
  pHMP_pMu_FS_2D_hist->Draw("colz");
  TitleFactory("pHMP_pMu_FS_2D" + gen + tar)->Draw();

  return mycanv;
}

TCanvas* pHMP_pfm_FS_2D(std::string gen, std::string tar){
// HMPMomVsStruckNucleonMom
  TCanvas* mycanv = CanvasFactory("pHMP_pfm_FS_2D" + gen + tar);

  TH2* pHMP_pfm_FS_2D_hist = NoSaveClone(FillHistogram2D(gen,tar,
    "HMPMomVsStruckNucleonMom"));

  SetOnePaneHistDefaults(pHMP_pfm_FS_2D_hist,"Proton Momenta [MeV/C]",
    "Nucleon Initial Momenta [MeV/C]");

  pHMP_pfm_FS_2D_hist->GetYaxis()->SetTitleOffset(YTitleOffset*0.6);

  pHMP_pfm_FS_2D_hist->GetXaxis()->SetTitleOffset(XTitleOffset*0.75);
  pHMP_pfm_FS_2D_hist->Draw("colz");
  TitleFactory("pHMP_pfm_FS_2D" + gen + tar)->Draw();

  return mycanv;
}

TCanvas* ERec_pmu_pp(std::string en_suff=""){
  TCanvas* canv = CanvasFactory("ERec_pmu_pp"+en_suff);

  TH1* NEUT_pmu = FillHistogram1D("NEUT","C"+en_suff,"CCQE_mmtm_mu");
  TH1* NEUT_pp = FillHistogram1D("NEUT","C"+en_suff,"CCQE_mmtm_p");

  TH1* GENIE_pmu = FillHistogram1D("GENIE","C"+en_suff,"CCQE_mmtm_mu");
  TH1* GENIE_pp = FillHistogram1D("GENIE","C"+en_suff,"CCQE_mmtm_p");

  TH1* NuWro_pmu = FillHistogram1D("NuWro","C"+en_suff,"CCQE_mmtm_mu");
  TH1* NuWro_pp = FillHistogram1D("NuWro","C"+en_suff,"CCQE_mmtm_p");

  // TH1* NEUT_pmu = FillHistogram1D("NEUT","C"+en_suff,"CCQE_mmtm_mu");
  // TH1* NEUT_pmu = FillHistogram1D("NEUT","C"+en_suff,"CCQE_mmtm_mu");

  std::initializer_list<TH1*> Hists = {
            NEUT_pmu,
            GENIE_pmu,
            NuWro_pmu,
            NEUT_pp,
            GENIE_pp,
            NuWro_pp};

  SetOnePaneHistDefaults(NEUT_pmu,"Count","Momentum [MeV/C]");

  SetUserRangeToMaxYVal(Hists);

  SetStyles(Hists);

  AddCountToTHTitle({std::make_pair("NEUT p_{#mu}",NEUT_pmu),
    std::make_pair("NEUT p_{p}",NEUT_pp),
    std::make_pair("GENIE p_{#mu}",GENIE_pmu),
    std::make_pair("GENIE p_{p}",GENIE_pp),
    std::make_pair("NuWro p_{#mu}",NuWro_pmu),
    std::make_pair("NuWro p_{p}",NuWro_pp)});

  TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*2.0 -Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));

  AddToLegend(Hists,leg);

  DrawTHs(Hists);
  leg->Draw();

  TitleFactory("ERec_pmu_pp"+en_suff)->Draw();

  return canv;
}

TCanvas* ERecRatioComp(std::string Generator, std::string Target){

  TCanvas* canv = CanvasFactory(Generator+Target+"ERecRatioComp");

  TH1* ERec_CCQE =
    FillHistogram1D(Generator,Target,"CCQENuERec_CCQENuETrue");
  TH1* ERec_MuP =
    FillHistogram1D(Generator,Target,"CCQENuMuPERec_CCQENuETrue");

  SetOnePaneHistDefaults(ERec_CCQE,"Count","(E_Rec/E_True)");


  ERec_CCQE->SetTitle("CCQE Formula E_{rec}/E_{true}");
  ERec_MuP->SetTitle("E_{rec}/E_{true} = (E_{#mu} + E_{p} - m_n + E_b)/E_{true}");

  SetStyles({ERec_CCQE,
    ERec_MuP});

  ERec_CCQE->GetXaxis()->SetRangeUser(0,2.5);

  SetUserRangeToMaxYVal({ERec_CCQE,
    ERec_MuP});

  TLegend* leg = LegendFactory( (1.0 - P0LegendWidth*2.0 -Pad0XMargin),
    (1.0 - P0LegendHeight - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  leg->AddEntry(ERec_CCQE,ERec_CCQE->GetTitle(),"l");

  std::stringstream ss("");
  ss << "Mean: " << ERec_CCQE->GetMean() << ", #sigma: "
    << ERec_CCQE->GetStdDev();

  leg->AddEntry((TObject*)0, ss.str().c_str(), "");
  leg->AddEntry(ERec_MuP,ERec_MuP->GetTitle(),"l");
  ss.str("");
  ss << "Mean: " << ERec_MuP->GetMean() << ", #sigma: "
    << ERec_MuP->GetStdDev();
  leg->AddEntry((TObject*)0, ss.str().c_str(), "");


  DrawTHs({ERec_CCQE,
    ERec_MuP});

  leg->Draw();
  TitleFactory(Generator + " " + Target + " ERec Ratio Comparison" )->Draw();
  return canv;
}

TCanvas* ERecSigmaVsETrueCorr(std::string Generator, std::string Target){

  TCanvas* canv = CanvasFactory(Generator+Target+"ERecSigmaVsETrueCorr");

  TH1* ERec_CCQE_600MeV =
    FillHistogram1D(Generator,Target,"CCQENuERec_CCQENuETrue");
  TH1* ERec_MuP_600MeV =
    FillHistogram1D(Generator,Target,"CCQENuMuPERec_CCQENuETrue");

  TH1* ERec_CCQE_1GeV =
    FillHistogram1D(Generator,Target+"_1GeV","CCQENuERec_CCQENuETrue");
  TH1* ERec_MuP_1GeV =
    FillHistogram1D(Generator,Target+"_1GeV","CCQENuMuPERec_CCQENuETrue");

  TH1* ERec_CCQE_2GeV =
    FillHistogram1D(Generator,Target+"_2GeV","CCQENuERec_CCQENuETrue");
  TH1* ERec_MuP_2GeV =
    FillHistogram1D(Generator,Target+"_2GeV","CCQENuMuPERec_CCQENuETrue");

  TH1* ERec_CCQE_4GeV =
    FillHistogram1D(Generator,Target+"_4GeV","CCQENuERec_CCQENuETrue");
  TH1* ERec_MuP_4GeV =
    FillHistogram1D(Generator,Target+"_4GeV","CCQENuMuPERec_CCQENuETrue");

  TH1* Correlations_CCQE = new TH1F("ERecCorrelations_CCQE","CCQE Formula E_{rec}",20,0,5000);
  TH1* Correlations_MuP = new TH1F("ERecCorrelations_MuP","E_{rec} = E_{#mu} + E_{p} - m_n + E_b",20,0,5000);

  Correlations_CCQE->SetDirectory(0);
  Correlations_MuP->SetDirectory(0);

  for (auto const &pair : {std::make_pair(600,ERec_CCQE_600MeV),
                            std::make_pair(1000,ERec_CCQE_1GeV),
                            std::make_pair(2000,ERec_CCQE_2GeV),
                            std::make_pair(4000,ERec_CCQE_4GeV)}){
    int bin =
      Correlations_CCQE->FindBin(pair.first);
    Correlations_CCQE->SetBinContent(bin,pair.second->GetMean());
    Correlations_CCQE->SetBinError(bin,pair.second->GetStdDev());
    std::cout << "CCQE - ETrue [" << pair.first << "] - <ERec/ETrue> ["
      << pair.second->GetMean()
      << "], Sigma(ERec/ETrue) [" << pair.second->GetStdDev() << "], bin == "
      << bin << std::endl;
  }

  for (auto const &pair : {std::make_pair(600,ERec_MuP_600MeV),
                            std::make_pair(1000,ERec_MuP_1GeV),
                            std::make_pair(2000,ERec_MuP_2GeV),
                            std::make_pair(4000,ERec_MuP_4GeV)}){
    int bin =
      Correlations_MuP->FindBin(pair.first);
    Correlations_MuP->SetBinContent(bin,pair.second->GetMean());
    Correlations_MuP->SetBinError(bin,pair.second->GetStdDev());
    std::cout << "Mu P - ETrue [" << pair.first << "] - <ERec/ETrue> ["
      << pair.second->GetMean()
      << "], Sigma(ERec/ETrue) [" << pair.second->GetStdDev() << "], bin == "
      << bin << std::endl;
  }

  SetOnePaneHistDefaults(Correlations_CCQE,"#LT(E_Rec/E_True)#GT","E_True [MeV]");

  SetUserRangeToMaxYVal({Correlations_CCQE,
    Correlations_MuP},0.8);

  SetStyles({Correlations_CCQE,
    Correlations_MuP});

  DrawTHs({Correlations_CCQE,
    Correlations_MuP});

  TLegend* leg = LegendFactory( (Pad0XMargin),
    (Pad0YMargin),
    (Pad0XMargin + P0LegendWidth*1.5),
    (Pad0YMargin + P0LegendHeight) );

  AddToLegend({Correlations_CCQE,
    Correlations_MuP},leg);
  leg->Draw();

  TitleFactory(Generator + " " + Target + " ERec Method Comparison" )->Draw();

  return canv;
}

TCanvas* FSIModeCompare(std::string Generator, std::string Target,
  std::string STVar="", bool DoPDF=PlottingDefaults::DoPDF, bool regen=false){

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar==""){
    XTitleText = DPhiTTitle;
    YTitleText = DPhiT_PDF_Title;
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = DAlphaT_PDF_Title;
  } else if(STVar=="DPt"){
    XTitleText = DPTTitle;
    YTitleText = DPT_PDF_Title;
  } else {
    return new TCanvas("broken");
  }

  if(!DoPDF){
    YTitleText = "Count";
  }

  TCanvas* canv = CanvasFactory(Generator+Target+"FSIMODE"+STVar);

  TH1* CCQEAll_norescat =
    NoSaveClone(FillHistogram1D(Generator,Target+"_norescat","CCQE"+STVar,regen));
  // TH1* CCQEAll_NoNE_norescat =
  //   NoSaveClone(FillHistogram1D(Generator,Target+"_norescat","CCQENoNE"+STVar,regen));
  TH1* CCQEAll = NoSaveClone(FillHistogram1D(Generator,Target,"CCQE"+STVar,regen));
  TH1* CCQEAll_NE = NoSaveClone(FillHistogram1D(Generator,Target,"CCQENE"+STVar,regen));
  TH1* CCQEAll_NoNE = NoSaveClone(FillHistogram1D(Generator,Target,"CCQENoNE"+STVar,regen));

  std::initializer_list<TH1*> hists;
  if(false && Generator == "GENIE" && Target == "C"){
    TH1* CCQEAll_hA = NoSaveClone(FillHistogram1D(Generator,"C_hA","CCQE"+STVar,regen));

    hists = {CCQEAll_norescat,
              CCQEAll,
              CCQEAll_NE,
              CCQEAll_NoNE,
              // CCQEAll_NoNE_norescat,
              CCQEAll_hA
            };
    AddCountToTHTitle({std::make_pair("CCQE (hA Model)",CCQEAll_hA)});
    CCQEAll_hA->SetMarkerStyle(0);
    CCQEAll_hA->SetLineWidth(2);
    CCQEAll_hA->SetLineColor(kCyan-2);
  } else {
    hists = {CCQEAll_norescat,
              CCQEAll,
              CCQEAll_NE,
              CCQEAll_NoNE,
              // CCQEAll_NoNE_norescat
            };
  }

  AddCountToTHTitle({std::make_pair("CCQE (Generator HT Off)",CCQEAll_norescat),
              std::make_pair("CCQE",CCQEAll),
              std::make_pair("CCQE With Nuclear Emission",CCQEAll_NE),
              std::make_pair("CCQE With No Nuclear Emission",CCQEAll_NoNE),
              // std::make_pair("CCQE (Generator HT Off) No NE",CCQEAll_NoNE_norescat)
            });

  if(DoPDF){
    MakePDF(hists);
  }
  SetUserRangeToMaxYVal(hists);

  SetStyles(hists);

  CCQEAll->SetLineColor(kBlue-3);
  CCQEAll_norescat->SetLineColor(kRed+2);
  CCQEAll_NE->SetLineColor(kMagenta-2);
  CCQEAll_NoNE->SetLineColor(kBlack);

  SetOnePaneHistDefaults(CCQEAll_norescat,YTitleText.c_str(),
  XTitleText.c_str());

  DrawTHs(hists);

  TLegend* leg;
  if(STVar==""){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*2.8),
    (1.0 - P0LegendHeight*2.0 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  } else if(STVar=="AlphaT"){
    leg = LegendFactory( (Pad0XMargin),
    (1.0 - P0LegendHeight*1.5 - Pad0YMargin),
    (Pad0XMargin + P0LegendWidth*2.1),
    (1.0 - Pad0YMargin));
  } else if(STVar=="DPt"){
    leg = LegendFactory( (1.0 - Pad0XMargin - P0LegendWidth*2.0),
    (1.0 - P0LegendHeight*1.3 - Pad0YMargin),
    (1.0 - Pad0XMargin),
    (1.0 - Pad0YMargin));
  }

  AddToLegend(hists,leg);
  leg->Draw();

  TitleFactory((std::string(Generator+ " " + Target
    + " FSI Mode Shapes")).c_str())->Draw();
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

TCanvas* NuclearEmission(std::string generator){
  TH1* GEN_C_CCQEAll = NoSaveClone(FillHistogram1D(generator,"C","CCQE"));
  TH1* GEN_C_CCQEAllNE = NoSaveClone(FillHistogram1D(generator,"C","CCQENE"));
  TH1* GEN_C_CCQEAllNoNE = NoSaveClone(FillHistogram1D(generator,"C","CCQENoNE"));

  TH1* GEN_C_CCQEAll_zoom =
    NoSaveClone(FillHistogram1D(generator,"C","CCQEZoom"));
  TH1* GEN_C_CCQEAllNE_zoom =
    NoSaveClone(FillHistogram1D(generator,"C","NEZoom"));
  TH1* GEN_C_CCQEAllNoNE_zoom =
    NoSaveClone(FillHistogram1D(generator,"C","CCQENoNEZoom"));

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

TCanvas* DPHITCCQEGENCOMPNORESCAT(bool usecop=false, std::string Vers=""){

  std::string var = usecop? "CCQECoplanarity" : "CCQE";

  // TH1* NEUT_C_CCQEAll = NoSaveClone(FillHistogram1D("NEUT","C",var));
  // TH1* NEUT_C_CCQEAll_norescat =
  //   NoSaveClone(FillHistogram1D("NEUT","C_norescat",var));
  TH1* GENIE_C_CCQEAll = NoSaveClone(FillHistogram1D("GENIE"+Vers,"C_ha",var));
  TH1* GENIE_C_CCQEAll_norescat =
    NoSaveClone(FillHistogram1D("GENIE"+Vers,"C_norescat",var));
  TH1* GENIE_C_CCQEAll_hN =
    NoSaveClone(FillHistogram1D("GENIE"+Vers,"C_hn",var));
  // TH1* NuWro_C_CCQEAll = NoSaveClone(FillHistogram1D("NuWro","C",var));
  // TH1* NuWro_C_CCQEAll_norescat =
  //   NoSaveClone(FillHistogram1D("NuWro","C_norescat",var));

  // TH1* GiBUU_C_CCQEAll = NoSaveClone(FillHistogram1D("GiBUU","C",var));
  // TH1* GiBUU_C_CCQEAll_norescat =
  //   NoSaveClone(FillHistogram1D("GiBUU","C_norescat",var));

  // TH1* NEUT_C_CCQEAll_zoom = NoSaveClone(FillHistogram1D("NEUT","C",var+"Zoom"));
  // TH1* NEUT_C_CCQEAll_norescat_zoom =
  //   NoSaveClone(FillHistogram1D("NEUT","C_norescat",var+"Zoom"));
  TH1* GENIE_C_CCQEAll_zoom = NoSaveClone(FillHistogram1D("GENIE"+Vers,"C_ha",var+"Zoom"));
  TH1* GENIE_C_CCQEAll_norescat_zoom =
    NoSaveClone(FillHistogram1D("GENIE"+Vers,"C_norescat",var+"Zoom"));
  TH1* GENIE_C_CCQEAll_hN_zoom =
    NoSaveClone(FillHistogram1D("GENIE"+Vers,"C_hn",var+"Zoom"));
  // TH1* NuWro_C_CCQEAll_zoom = NoSaveClone(FillHistogram1D("NuWro","C",var+"Zoom"));
  // TH1* NuWro_C_CCQEAll_norescat_zoom =
  //   NoSaveClone(FillHistogram1D("NuWro","C_norescat",var+"Zoom"));

  // TH1* GiBUU_C_CCQEAll_zoom = NoSaveClone(FillHistogram1D("GiBUU","C",var+"Zoom"));
  // TH1* GiBUU_C_CCQEAll_norescat_zoom =
  //   NoSaveClone(FillHistogram1D("GiBUU","C_norescat",var+"Zoom"));


  GENIE_C_CCQEAll->SetStats(false);
  GENIE_C_CCQEAll_norescat->SetStats(false);
  GENIE_C_CCQEAll_hN->SetStats(false);
  GENIE_C_CCQEAll_zoom->SetStats(false);
  GENIE_C_CCQEAll_norescat_zoom->SetStats(false);
  GENIE_C_CCQEAll_hN_zoom->SetStats(false);

  std::initializer_list<TH1*> BigHists = {
    // NEUT_C_CCQEAll,
    GENIE_C_CCQEAll,
    // NuWro_C_CCQEAll,
    GENIE_C_CCQEAll_hN,
    // GiBUU_C_CCQEAll,
    // NEUT_C_CCQEAll_norescat,
    // NuWro_C_CCQEAll_norescat,
    GENIE_C_CCQEAll_norescat,
    // GiBUU_C_CCQEAll_norescat,
  };

  std::initializer_list<TH1*> InsetHists = {
    // NEUT_C_CCQEAll_zoom,
    GENIE_C_CCQEAll_zoom,
    // NuWro_C_CCQEAll_zoom,
    GENIE_C_CCQEAll_hN_zoom,
    // GiBUU_C_CCQEAll_zoom,
    // NEUT_C_CCQEAll_norescat_zoom,
    // NuWro_C_CCQEAll_norescat_zoom,
    GENIE_C_CCQEAll_norescat_zoom,
    // GiBUU_C_CCQEAll_norescat_zoom,
  };

  AddCountToTHTitle({
    // std::make_pair("NEUT",NEUT_C_CCQEAll),
    std::make_pair("GENIE (hA)",GENIE_C_CCQEAll),
    // std::make_pair("NuWro",NuWro_C_CCQEAll),
    std::make_pair("GENIE (hN Cascade Model)",GENIE_C_CCQEAll_hN),
    // std::make_pair("GiBUU",GiBUU_C_CCQEAll),
    // std::make_pair("NEUT (NUCRESCAT = 0)",NEUT_C_CCQEAll_norescat),
    std::make_pair("GENIE (Hadronic Transport = false)",GENIE_C_CCQEAll_norescat),
    // std::make_pair("NuWro (Kaskada_on = 0)",NuWro_C_CCQEAll_norescat),
    // std::make_pair("GiBUU (numTimeSteps=0)",GiBUU_C_CCQEAll_norescat)
  });

  MakePDF(BigHists);
  MakePDF(InsetHists);

  SetStyles(BigHists);
  SetStyles(InsetHists);

  std::string YTitle=DPhiT_PDF_Title, XTitle=DPhiTTitle;
  if(usecop){
    YTitle = CopA_PDF_Title;
    XTitle = CopATitle;
  }

  SetOnePaneHistDefaults((*BigHists.begin()), YTitle.c_str(), XTitle.c_str());
  (*BigHists.begin())->GetYaxis()->SetTitleOffset(
    (*BigHists.begin())->GetYaxis()->GetTitleOffset()*0.6);
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

TCanvas* MuPTGenComp(std::string Target){
// CCQEMuon_Pt

  TCanvas* canv = CanvasFactory(Target+"mupt");

  TH1* NEUT_C_pt =
    NoSaveClone(FillHistogram1D("NEUT",Target,"CCQEMuon_Pt"));
  TH1* GENIE_C_pt =
    NoSaveClone(FillHistogram1D("GENIE",Target,"CCQEMuon_Pt"));
  TH1* NuWro_C_pt =
    NoSaveClone(FillHistogram1D("NuWro",Target,"CCQEMuon_Pt"));
  TH1* GiBUU_C_pt =
    NoSaveClone(FillHistogram1D("GiBUU",Target,"CCQEMuon_Pt"));

  NEUT_C_pt->SetTitle("NEUT");
  GENIE_C_pt->SetTitle("GENIE");
  NuWro_C_pt->SetTitle("NuWro");
  GiBUU_C_pt->SetTitle("GiBUU");

  std::initializer_list<TH1*> hists = {NEUT_C_pt,GENIE_C_pt,
    NuWro_C_pt,GiBUU_C_pt};

  MakePDF(hists);
  SetStyles(hists);

  SetOnePaneHistDefaults(NEUT_C_pt, "P.D.F", "#mu_{P_{T}} [MeV]");
  SetUserRangeToMaxYVal({NEUT_C_pt,
                        GENIE_C_pt,
                        NuWro_C_pt,
                        GiBUU_C_pt});
  DrawTHs(hists);

  TLegend* leg;

  leg = LegendFactory( (Pad0XMargin),
    (Pad0YMargin),
    (Pad0XMargin + P0LegendWidth),
    (Pad0YMargin + P0LegendHeight));


  AddToLegend(hists,leg);
  leg->Draw();

  TitleFactory("#mu_{P_{T}} from " + Target + " Target" )->Draw();

  return canv;

}

TCanvas* STVarMuPT2D(std::string Generator, std::string STVar, std::string Target, bool NoNe=false){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  TCanvas* canv = CanvasFactory(STVar+Target+"vsmupt");

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    YTitleText = DPhiTTitle;
    XTitleText = "#mu_{p_{T}} [MeV]";
  } else if(STVar=="AlphaT"){
    YTitleText = DAlphaTTitle;
    XTitleText = "#mu_{p_{T}} [MeV]";
  } else if(STVar=="DPT"){
    YTitleText = DPTTitle;
    XTitleText = "#mu_{p_{T}} [MeV]";
  } else {
    return new TCanvas("broken");
  }

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,STVar+"HMPvsMuonPt"+(NoNe?"_NoNE":"")));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  TitleFactory(Generator + " " + STVar + " " + Target)->Draw();


  C_vsmupt->Draw("colz");
  return canv;
}

TCanvas* STVarMuPT2D_slicenorm(std::string Generator,
  std::string STVar, std::string Target, bool lastbool, bool NoNe){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  TCanvas* canv = CanvasFactory(STVar+Target+"vsmupt");

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    YTitleText = DPhiTTitle;
    XTitleText = "p^{#mu}_{T} [MeV/c]";
  } else if(STVar=="AlphaT"){
    YTitleText = DAlphaTTitle;
    XTitleText = "p^{#mu}_{T} [MeV/c]";
  } else if(STVar=="DPT"){
    YTitleText = DPTTitle;
    XTitleText = "p^{#mu}_{T} [MeV/c]";
  } else {
    return new TCanvas("broken");
  }

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,STVar+"HMPvsMuonPt"+(NoNe?"_NoNE":"")));

  TH2* C_vsmupt_norm = NeutrinoTools::NormalHist(C_vsmupt,20,lastbool);

  SetOnePaneHistDefaults(C_vsmupt_norm,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt_norm->Draw("colz");
  SetUserRangeToMaxXValWithContent({C_vsmupt_norm});

  TitleFactory(Generator + " " + STVar + " " + Target +(NoNe?"_NoNE":"") + "Y slice normalised")->Draw();
  return canv;
}

TCanvas* STVarStruckNucleonZMom2D(std::string Generator, std::string STVar, std::string Target){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  TCanvas* canv = CanvasFactory(STVar+Target+"vsStruckNucleonZMomentum");

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    YTitleText = DPhiTTitle;
    XTitleText = "Nucleon Initial Z Momentum [MeV]";
  } else if(STVar=="AlphaT"){
    YTitleText = DAlphaTTitle;
    XTitleText = "Nucleon Initial Z Momentum [MeV]";
  } else if(STVar=="DPT"){
    YTitleText = DPTTitle;
    XTitleText = "Nucleon Initial Z Momentum [MeV]";
  } else {
    return new TCanvas("broken");
  }

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,STVar+"HMPvsStruckNucleonZMomentum"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + STVar + " " + Target)->Draw();

  return canv;
}

TCanvas* STVarStruckNucleonZMomProj(std::string Generator, std::string STVar, std::string Target){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  TCanvas* canv = CanvasFactory(STVar+Target+"vsStruckNucleonZMomentum");

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT"){
    XTitleText = DPhiTTitle;
    YTitleText = "Struck Nucleon #LTZ Momentum#GT [MeV]";
  } else if(STVar=="AlphaT"){
    XTitleText = DAlphaTTitle;
    YTitleText = "Struck Nucleon #LTZ Momentum#GT [MeV]";
  } else if(STVar=="DPT"){
    XTitleText = DPTTitle;
    YTitleText = "Struck Nucleon #LTZ Momentum #GT[MeV]";
  } else {
    return new TCanvas("broken");
  }

  TH2* Tar_vsStruckNucZMom =
    NoSaveClone(FillHistogram2D(Generator,Target,STVar+"HMPvsStruckNucleonZMomentum"));

  TProfile* Tar_vsStruckNucZMom_p = Tar_vsStruckNucZMom->ProfileY();

  SetOnePaneHistDefaults(Tar_vsStruckNucZMom_p,YTitleText.c_str(),XTitleText.c_str());

  Tar_vsStruckNucZMom_p->Draw("EHIST");

  TitleFactory(Generator + " " + STVar + " " + Target + " Correlation")->Draw();

  return canv;
}

//OmegaHMPvsDAlphaT

TCanvas* OmegaHMPvsDAlphaT(){
  std::string Generator="NEUT";
  std::string Target="C";
  TCanvas* canv = CanvasFactory(Target+"OmegaHMPvsDAlphaT");

  std::string YTitleText = "#omega = E^{lab}_{#nu} - E^{lab}_{#nu} [MeV]";
  std::string XTitleText = DAlphaTTitle;

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"OmegaHMPvsDAlphaT"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " OmegaHMPvsDAlphaT")->Draw();

  return canv;
}

//NProtonsVsCCQEDphiT

TCanvas* NProtonsVsCCQE(std::string Generator, std::string Target,
  std::string var, bool regen=false){
  // std::string Generator="NEUT";
  // std::string Target="C";
  TCanvas* canv = CanvasFactory(Target+Generator+"NProtonsVsCCQE"+var);

  std::string YTitleText = "NProtons";
  std::string XTitleText = var;

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"NProtonsVsCCQE"+var, regen));

  TH2* C_vsmupt_norm = NeutrinoTools::NormalHist(C_vsmupt,20,true);

  SetOnePaneHistDefaults(C_vsmupt_norm,YTitleText.c_str(),XTitleText.c_str());
  canv->SetLogz(true);
  C_vsmupt_norm->Draw("colz");

  TitleFactory(Generator + " " + Target + " NProtonsVsCCQE"+var)->Draw();

  return canv;
}


//LabMomVsStruckNucleonMom

TCanvas* LabMomVsStruckNucleonMom(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory(Target+"LabMomVsStruckNucleonMom");

  std::string YTitleText = "p^{lab}_{#mu} + p^{lab}_{p} - p^{lab}_{#nu} [MeV/C]";
  std::string XTitleText = "Nucleon Initial p^{lab} [MeV/C]";

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"LabMomVsStruckNucleonMom"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " LabMomVsStruckNucleonMom")->Draw();

  return canv;
}

TCanvas* DPTVsStruckNucleonPt(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory(Target+"vsStruckNucleonPt");

  std::string YTitleText = "p^{lab}_{#mu} + p^{lab}_{p} - p^{lab}_{#nu} [MeV/C]";
  std::string XTitleText = "Nucleon Initial p^{lab}_{T} [MeV]";

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"DPTHMPvsStruckNucleonPt"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " DPTHMPvsStruckNucleonPt")->Draw();

  return canv;
}

TCanvas* DPTFirstPVsStruckNucleonPt(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory(Target+"vsStruckNucleonPt_firstp");

  std::string YTitleText = "p^{lab}_{#mu} + p^{lab}_{p} - p^{lab}_{#nu} [MeV/C]";
  std::string XTitleText = "Nucleon Initial p^{lab}_{T} [MeV]";

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"DPTFirstPvsStruckNucleonPt"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " DPTHMPvsStruckNucleonPt_firstp")->Draw();

  return canv;
}

TCanvas* DPTHMPvsStruckNucleonPtNE(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory(Target+"vsStruckNucleonPtNE");

  std::string YTitleText = "p^{lab}_{#mu} + p^{lab}_{p} - p^{lab}_{#nu} [MeV/C]";
  std::string XTitleText = "Nucleon Initial p^{lab}_{T} [MeV]";

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"DPTHMPvsStruckNucleonPtNE"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " DPTHMPvsStruckNucleonPtNE")->Draw();

  return canv;
}


TCanvas* DPTHMPvsStruckNucleonPtNE_excludediagonal(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory(Target+"DPTHMPvsStruckNucleonPtNE_excludediagonal");

  std::string YTitleText = "p^{lab}_{#mu} + p^{lab}_{p} - p^{lab}_{#nu} [MeV/C]";
  std::string XTitleText = "Nucleon Initial p^{lab}_{T} [MeV]";

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"DPTHMPvsStruckNucleonPtNE_excludediagonal"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " DPTHMPvsStruckNucleonPtNE_excludediagonal")->Draw();

  return canv;
}


TCanvas* DPTNoNEVsStruckNucleonPt(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory(Target+"vsStruckNucleonPt");

  std::string YTitleText = DPTTitle;
  std::string XTitleText = "Nucleon Initial P_{T} [MeV]";

  TH2* C_vsmupt =
    NoSaveClone(FillHistogram2D(Generator,Target,"DPTHMPvsStruckNucleonPtNoNE"));

  SetOnePaneHistDefaults(C_vsmupt,YTitleText.c_str(),XTitleText.c_str());

  C_vsmupt->Draw("colz");

  TitleFactory(Generator + " " + Target + " DPTHMPvsStruckNucleonPtNoNE")->Draw();

  return canv;
}

TCanvas* STVarMuPTCorrelations(std::string STVar, std::string Target,
  bool DoProfile=true, bool regen=false){
// DeltaPhiTHMPvsMuonPt
// AlphaTHMPvsMuonPt
// DeltaPTHMPvsMuonPt

  TCanvas* canv = CanvasFactory(STVar+Target+"vsmupt");

  std::string XTitleText = "";
  std::string YTitleText = "";
  if(STVar=="DeltaPhiT" || STVar=="ECalTop2Prong_sm_DeltaPhiT"){
    YTitleText =  DoProfile?"#LT#delta#phi_{T}#GT [deg]":
      "{#delta#phi_{T}}_{MLE} [deg]";
    XTitleText = "p^{#mu}_{T} [MeV/c]";
  } else if(STVar=="AlphaT"){
    YTitleText = DoProfile?"#LT#delta#alpha_{T}#GT [deg]":
      "{#delta#alpha_{T}}_{MLE} [deg]";
    XTitleText = "p^{#mu}_{T} [MeV/c]";
  } else if(STVar=="DPT"){
    YTitleText =  DoProfile?"#LT#deltap_{T}#GT [MeV/c]":
      "{#deltap_{T}}_{MLE} [MeV/c]";
    XTitleText = "p^{#mu}_{T} [MeV/c]";
  } else {
    return new TCanvas("broken");
  }

  TH2* NEUT_C_vsmupt =
    NoSaveClone(FillHistogram2D("NEUT",Target,STVar+"HMPvsMuonPt",regen));
  TH2* GENIE_C_vsmupt =
    NoSaveClone(FillHistogram2D("GENIE",Target,STVar+"HMPvsMuonPt",regen));
  TH2* NuWro_C_vsmupt =
    NoSaveClone(FillHistogram2D("NuWro",Target,STVar+"HMPvsMuonPt",regen));
  TH2* GiBUU_C_vsmupt =
    NoSaveClone(FillHistogram2D("GiBUU",Target,STVar+"HMPvsMuonPt",regen));

  TH1* NEUT_C_vsmupt_p;
  TH1* GENIE_C_vsmupt_p;
  TH1* NuWro_C_vsmupt_p;
  TH1* GiBUU_C_vsmupt_p;

  if(DoProfile){
    NEUT_C_vsmupt_p = NEUT_C_vsmupt->ProfileX();
    GENIE_C_vsmupt_p = GENIE_C_vsmupt->ProfileX();
    NuWro_C_vsmupt_p = NuWro_C_vsmupt->ProfileX();
    GiBUU_C_vsmupt_p = GiBUU_C_vsmupt->ProfileX();
  } else {
    NEUT_C_vsmupt_p = MakePeakCorrXTH1(NEUT_C_vsmupt);
    GENIE_C_vsmupt_p = MakePeakCorrXTH1(GENIE_C_vsmupt);
    NuWro_C_vsmupt_p = MakePeakCorrXTH1(NuWro_C_vsmupt);
    GiBUU_C_vsmupt_p = MakePeakCorrXTH1(GiBUU_C_vsmupt);
  }

  NEUT_C_vsmupt_p->SetTitle("NEUT");
  GENIE_C_vsmupt_p->SetTitle("GENIE");
  NuWro_C_vsmupt_p->SetTitle("NuWro");
  GiBUU_C_vsmupt_p->SetTitle("GiBUU");

  SetStyles({NEUT_C_vsmupt_p,
            GENIE_C_vsmupt_p,
            NuWro_C_vsmupt_p,
            GiBUU_C_vsmupt_p});


  SetOnePaneHistDefaults(NEUT_C_vsmupt_p, YTitleText.c_str(), XTitleText.c_str());
  SetUserRangeToMaxYVal({NEUT_C_vsmupt_p,
                          GENIE_C_vsmupt_p,
                          NuWro_C_vsmupt_p,
                          GiBUU_C_vsmupt_p});

  if(STVar=="ECalTop2Prong_sm_DeltaPhiT"){
    NEUT_C_vsmupt_p->GetYaxis()->SetRangeUser(0,180);
    GENIE_C_vsmupt_p->GetYaxis()->SetRangeUser(0,180);
    NuWro_C_vsmupt_p->GetYaxis()->SetRangeUser(0,180);
    GiBUU_C_vsmupt_p->GetYaxis()->SetRangeUser(0,180);
  } else {
    SetUserRangeToMaxXValWithContent({NEUT_C_vsmupt_p,
                                      GENIE_C_vsmupt_p,
                                      NuWro_C_vsmupt_p,
                                      GiBUU_C_vsmupt_p});
  }

  NEUT_C_vsmupt_p->Draw("EHIST");
  GENIE_C_vsmupt_p->Draw("EHIST SAME");
  NuWro_C_vsmupt_p->Draw("EHIST SAME");
  GiBUU_C_vsmupt_p->Draw("EHIST SAME");

  TLegend* leg;

  if(STVar=="DeltaPhiT" || STVar=="ECalTop2Prong_sm_DeltaPhiT"){
    leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
      (1.0 - Pad0YMargin - P0LegendHeight),
      (1.0 - Pad0XMargin),
      (1.0 - Pad0YMargin));
  } else if(STVar=="AlphaT"){
    leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
      (Pad0YMargin),
      (1.0 - Pad0XMargin),
      (Pad0YMargin + P0LegendHeight));
  } else if(STVar=="DPT"){
    leg = LegendFactory( (1.0 - P0LegendWidth - Pad0XMargin),
      (1.0 - Pad0YMargin - P0LegendHeight),
      (1.0 - Pad0XMargin),
      (1.0 - Pad0YMargin));
  }

  AddToLegend({NEUT_C_vsmupt_p,
                GENIE_C_vsmupt_p,
                NuWro_C_vsmupt_p,
                GiBUU_C_vsmupt_p},leg);
  leg->Draw();

  TitleFactory(STVar + " " + Target + " vs #mu_{P_{T}}" )->Draw();

  return canv;

}

TCanvas* StruckNucleonZMomentum(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory("StruckNucleonZMomentum" + Generator + Target);
  TH1* Gen_C_CCQEAll = NoSaveClone(FillHistogram1D(Generator,Target,
    "CCQEStruckNucleonZMomentum"));

  SetOnePaneHistDefaults(Gen_C_CCQEAll, "Count", "Struck Nucleon Z Momentum [MeV]");
  SetStyles({Gen_C_CCQEAll});
  Gen_C_CCQEAll->Draw("EHIST");

  TitleFactory(Generator + " " + Target + " Struck Nucleon Z Momentum")->Draw();

  return canv;
}

TCanvas* StruckNucleonMomentum(std::string Generator, std::string Target){
  TCanvas* canv = CanvasFactory("StruckNucleonMomentum" + Generator + Target);
  TH1* Gen_C_CCQEAll = NoSaveClone(FillHistogram1D(Generator,Target,
    "CCQEStruckNucleonMomentum"));

  SetOnePaneHistDefaults(Gen_C_CCQEAll, "Count", "Struck Nucleon Initial Momentum [MeV]");
  SetStyles({Gen_C_CCQEAll});
  Gen_C_CCQEAll->Draw("EHIST");

  TitleFactory(Generator + " " + Target + " Struck Nucleon  Momentum")->Draw();

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



  // MakeT2KStyle()->cd();

  // SaveAndPrint(TargetDependenceRatio("CCQE","Pb"));
  // SaveAndPrint(TargetDependenceRatio("CCQEAlphaT","Pb"));
  // SaveAndPrint(TargetDependenceRatio("ECalTop2Prong_CCQE","Pb"));
  // SaveAndPrint(GenVarComparison("ECalTop2Prong_CCQE","Pb"));
  // SaveAndPrint(GenVarComparison("CCQE","C"));
  // SaveAndPrint(GenVarComparison("CCQEDPt","C"));
  // SaveAndPrint(GenVarComparison("CCQEAlphaT","C"));
  // SaveAndPrint(TargetDependenceRatio("ECalTop2Prong","Pb"));
  // SaveAndPrint(TargetDependenceRatio("CCQE","Ar"));

  // SaveAndPrint(GenVarComparison("CC0Pi1P","C"));
  // SaveAndPrint(GenVarComparison("CC0Pi1PDPt","C"));
  // SaveAndPrint(GenVarComparison("CC0Pi1PAlphaT","C"));

  // SaveAndPrint(DrawDistrib("NEUT","Pb","","ECalTop2Prong_CCQE"));
  // SaveAndPrint(DrawDistrib("GENIE","Pb","","ECalTop2Prong_CCQE"));
  // SaveAndPrint(DrawDistrib("NuWro","Pb","","ECalTop2Prong_CCQE"));
  // SaveAndPrint(DrawDistrib("GiBUU","Pb","","ECalTop2Prong_CCQE"));

  // SaveAndPrint(DrawDistrib("NEUT","C_ND5","","CCQE"));
  // SaveAndPrint(DrawDistrib("NEUT","C","","CCQE"));

  // SaveAndPrint(FSIModeCompare("NEUT","C",""));
  // SaveAndPrint(FSIModeCompare("GENIE","C","",true));
  // SaveAndPrint(FSIModeCompare("NuWro","C","",false));
  // SaveAndPrint(FSIModeCompare("GiBUU","C",""));

  // SaveAndPrint(FSIModeCompare("NEUT","Pb",""));
  // SaveAndPrint(FSIModeCompare("GENIE","Pb",""));
  // SaveAndPrint(FSIModeCompare("NuWro","Pb",""));
  // SaveAndPrint(FSIModeCompare("GiBUU","Pb",""));

  // SaveAndPrint(CArPbTargetComparison("NEUT",""));
  // SaveAndPrint(CArPbTargetComparison("GENIE",""));
  // SaveAndPrint(CArPbTargetComparison("NuWro",""));
  // SaveAndPrint(CArPbTargetComparison("GiBUU",""));

  // SaveAndPrint(FSIModeCompare("NEUT","C","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GENIE","C","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NuWro","C","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NuWro","C.B_R","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GiBUU","C","AlphaT"));

  // SaveAndPrint(FSIModeCompare("NEUT","Pb","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GENIE","Pb","AlphaT"));
  // SaveAndPrint(FSIModeCompare("NuWro","Pb","AlphaT"));
  // SaveAndPrint(FSIModeCompare("GiBUU","Pb","AlphaT"));

  // SaveAndPrint(CArPbTargetComparison("NEUT","AlphaT"));
  // SaveAndPrint(CArPbTargetComparison("GENIE","AlphaT"));
  // SaveAndPrint(CArPbTargetComparison("NuWro","AlphaT"));
  // SaveAndPrint(CArPbTargetComparison("GiBUU","AlphaT"));

  // SaveAndPrint(FSIModeCompare("NEUT","C","DPt"));
  // SaveAndPrint(FSIModeCompare("GENIE","C","DPt",true));
  // SaveAndPrint(FSIModeCompare("NuWro","C","DPt"));
  // SaveAndPrint(FSIModeCompare("NuWro","C","DPt",false));
  // SaveAndPrint(FSIModeCompare("GiBUU","C","DPt"));

  // SaveAndPrint(FSIModeCompare("NEUT","Pb","DPt"));
  // SaveAndPrint(FSIModeCompare("GENIE","Pb","DPt"));
  // SaveAndPrint(FSIModeCompare("NuWro","Pb","DPt"));
  // SaveAndPrint(FSIModeCompare("GiBUU","Pb","DPt"));

  // SaveAndPrint(CArPbTargetComparison("NEUT","DPt"));
  // SaveAndPrint(CArPbTargetComparison("GENIE","DPt"));
  // SaveAndPrint(CArPbTargetComparison("NuWro","DPt"));
  // SaveAndPrint(CArPbTargetComparison("NuWro","StruckNucleonMomentum"));
  // SaveAndPrint(CArPbTargetComparison("GiBUU","DPt"));

  // SaveAndPrint(STVarMuPTCorrelations("DeltaPhiT","C"));
  // SaveAndPrint(STVarMuPTCorrelations("DeltaPhiT","C_norescat"));
  // SaveAndPrint(STVarMuPTCorrelations("AlphaT","C"));
  // SaveAndPrint(STVarMuPTCorrelations("AlphaT","C_norescat"));
  // SaveAndPrint(STVarMuPTCorrelations("DPT","C", true));
  // SaveAndPrint(STVarMuPTCorrelations("DPT","C_norescat", true));
  // SaveAndPrint(STVarMuPTCorrelations("DeltaPhiT","Pb", true));
  // SaveAndPrint(STVarMuPTCorrelations("DeltaPhiT","C", true));
  // SaveAndPrint(STVarMuPTCorrelations("DPT","C",false));
  // SaveAndPrint(STVarMuPTCorrelations("DPT","C_norescat",false));

  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_NoNE_Xianguo("DeltaPhiT","C"),"mean_deltaphi_t__vs_mupt__NoNE.eps");
  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_NoNE_Xianguo("DPT","C", false),"mean_dp_t__vs_mupt__NoNE.eps");

  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_Single_Xianguo("NuWro","DPTHMPvsMuonPt","C",false),"mean_dp_t__vs_mupt__NuWro_C_single.eps");

  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_NoNE_Xianguo("DeltaPhiT","C.10mil",false),"mean_deltaphi_t__vs_mupt__NoNE.10mil.eps");
  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_NoNE_Xianguo("DPT","C.10mil",false),"mean_dp_t__vs_mupt__NoNE.10mil.eps");

  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_EDep_Xianguo("DeltaPhiT","C",false),"mean_deltaphi_t__vs_mupt__2GeV.eps");
  // SaveAndPrint(GeneratorDependence::STVarMuPTCorrelations_EDep_Xianguo("DPT","C", false),"mean_dp_t__vs_mupt__2GeV.eps");
  // SaveAndPrint(GeneratorDependence::GenComp_STV_Xianguo("C","AlphaT", false),"dalpha_t__gencomp.eps");

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo("C","DPt",false, 0,600),"emission_fraction_dp_t.eps");
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo("C","MuPt",false, 0,600),"emission_fraction_mu_p_t.eps");

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh_CC0Pi("C","CC0Pi","5","400","DPt",false, 0,600),"emission_fraction_CC0Pi_dp_t_lt5_ht400.eps");
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh_CC0Pi("C","CC0Pi","5","400","MuPt",false, 0,600),"emission_fraction_CC0Pi_mu_p_t_lt5_ht400.eps");

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh("C","CCQE", "5","400","DPt",false, 0,600),"emission_fraction_CCQE_dp_t_lt5_ht400.eps");
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh("C","CCQE", "5","400","MuPt",false, 0,600),"emission_fraction_CCQE_mu_p_t_lt5_ht400.eps");

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh("C","CCQE_lphi_", "5","400","DPt",false, 0,600),"emission_fraction_CCQE_lphi_dp_t_lt5_ht400.eps");
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh("C","CCQE_hphi_", "5","400","DPt",false, 0,600),"emission_fraction_CCQE_hphi_dp_t_lt5_ht400.eps");

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh("C","CCQE", "5","400","DPt",false, 0,600),"emission_fraction_CCQEdp_t_lt5_ht400_genie_2_9_0.eps");
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp_Xianguo_thresh("C","CCQE", "5","400","DPt",false, 0,600,"2_8_6"),"emission_fraction_CCQEdp_t_lt5_ht400_genie_2_8_6.eps");

  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C",1,false));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C_norescat",1,false));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C",1,false));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C_ND5",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","AlphaT","C_ND5",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DeltaPhiT","C_ND5",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C_1GeV",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C_2GeV",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C_4GeV",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","DPT","C_norescat",1));

  SaveAndPrint(DPHITCCQEGENCOMPNORESCAT(false,"_2_9_0"));
  SaveAndPrint(DPHITCCQEGENCOMPNORESCAT(false,"_2_8_6"));

  // SaveAndPrint(STVarMuPTCorrelations("ECalTop2Prong_sm_DeltaPhiT","Pb"));

  // SaveAndPrint(TargetDependence::SpeciesCanvas("NEUT","ccqeNProtons"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("GENIE","ccqeNProtons"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("NuWro","ccqeNProtons"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("GiBUU","ccqeNProtons"));

  // SaveAndPrint(TargetDependence::SpeciesCanvas("NEUT","ccqeNeNGamma"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("GENIE","ccqeNeNGamma"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("NuWro","ccqeNeNGamma"));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("GiBUU","ccqeNeNGamma"));

  // SaveAndPrint(TargetDependence::SpeciesCanvas("NEUT","ccqeNeNGamma", true));
  // SaveAndPrint(TargetDependence::SpeciesCanvas("GiBUU","ccqeNeNGamma", true));

  // SaveAndPrint(MuPTGenComp("C"));
  // SaveAndPrint(MuPTGenComp("Ar"));
  // SaveAndPrint(MuPTGenComp("Pb"));
  // SaveAndPrint(MuPTGenComp("O"));

  // SaveAndPrint(GeneratorDependence::EmissionFractionsSep("NEUT","C"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsSep("GENIE","C"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsSep("NuWro","C"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsSep("GiBUU","C"));

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("C","DPt"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("C","MuPt"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("C","AlphaT"));

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("Pb","DPt"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("Pb","MuPt"));

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("Ar","DPt"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("Pb","DPt"));

  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("C","MuPt"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("Ar","MuPt"));
  // SaveAndPrint(GeneratorDependence::EmissionFractionsGenComp("Pb","MuPt"));

  // SaveAndPrint(STVarMuPT2D("NEUT","AlphaT","C"));
  // SaveAndPrint(STVarMuPT2D("GENIE","AlphaT","C"));

  // SaveAndPrint(STVarMuPT2D_slicenorm("NEUT","AlphaT","C",1));
  // SaveAndPrint(STVarMuPT2D_slicenorm("GENIE","AlphaT","C",1));

  // SaveAndPrint(STVarStruckNucleonZMom2D("NEUT","AlphaT","C"));
  // SaveAndPrint(STVarStruckNucleonZMom2D("GENIE","AlphaT","C"));
  // SaveAndPrint(STVarStruckNucleonZMom2D("NuWro","AlphaT","C"));
  // SaveAndPrint(STVarStruckNucleonZMom2D("NEUT","AlphaT","C_norescat"));
  // SaveAndPrint(STVarStruckNucleonZMom2D("GENIE","AlphaT","C_norescat"));
  // SaveAndPrint(STVarStruckNucleonZMom2D("NuWro","AlphaT","C_norescat"));
  // SaveAndPrint(STVarStruckNucleonZMomProj("NEUT","AlphaT","C"));

  // SaveAndPrint(SVTFMModelCompare("CCQE"));
  // SaveAndPrint(SVTFMModelCompare("CCQEDPt"));
  // SaveAndPrint(SVTFMModelCompare("CCQEAlphaT"));
  // SaveAndPrint(SVTFMModelCompare("CCQEStruckNucleonMomentum"));
  // SaveAndPrint(SVTFMModelCompare("CCQEStruckNucleon_BE"));

  // SaveAndPrint(SVTFMModelCompare("CCQE",true));
  // SaveAndPrint(SVTFMModelCompare("CCQEDPt",true));
  // SaveAndPrint(SVTFMModelCompare("CCQEAlphaT",true));
  // SaveAndPrint(SVTFMModelCompare("CCQEStruckNucleonMomentum",true));
  // SaveAndPrint(SVTFMModelCompare("CCQEStruckNucleon_BE",true));

  // SaveAndPrint(SVTFancyMModelCompare("CCQE",true));
  // SaveAndPrint(SVTFancyMModelCompare("CCQEDPt",true));
  // SaveAndPrint(SVTFancyMModelCompare("CCQEAlphaT",true));
  // SaveAndPrint(SVTFancyMModelCompare("CCQEStruckNucleonMomentum",true));

  // SaveAndPrint(SVTFancyMModelCompare("CCQE"));
  // SaveAndPrint(SVTFancyMModelCompare("CCQEDPt"));
  // SaveAndPrint(SVTFancyMModelCompare("CCQEAlphaT"));
  // SaveAndPrint(SVTFancyMModelCompare("CCQEStruckNucleonMomentum"));

  // SaveAndPrint(ERecDistribComp("NEUT","C"));
  // SaveAndPrint(ERecDistribComp("NEUT","C_1GeV",900,1100));
  // SaveAndPrint(ERecDistribComp("NEUT","C_2GeV",1850,2100));
  // SaveAndPrint(ERecDistribComp("NEUT","C_4GeV",3850,4100));
  // SaveAndPrint(ERecRatioComp("NEUT","C"));
  // SaveAndPrint(ERecRatioComp("NEUT","C_1GeV"));
  // SaveAndPrint(ERecRatioComp("NEUT","C_2GeV"));
  // SaveAndPrint(ERecRatioComp("NEUT","C_4GeV"));

  // SaveAndPrint(ERecDistribComp("GENIE","C"));
  // SaveAndPrint(ERecDistribComp("GENIE","C_1GeV",900,1100));
  // SaveAndPrint(ERecDistribComp("GENIE","C_2GeV",1850,2100));
  // SaveAndPrint(ERecDistribComp("GENIE","C_4GeV",3850,4100));
  // SaveAndPrint(ERecRatioComp("GENIE","C"));
  // SaveAndPrint(ERecRatioComp("GENIE","C_1GeV"));
  // SaveAndPrint(ERecRatioComp("GENIE","C_2GeV"));
  // SaveAndPrint(ERecRatioComp("GENIE","C_4GeV"));

  // SaveAndPrint(ERecDistribComp("NuWro","C"));
  // SaveAndPrint(ERecDistribComp("NuWro","C_1GeV",900,1100));
  // SaveAndPrint(ERecDistribComp("NuWro","C_2GeV",1850,2100));
  // SaveAndPrint(ERecDistribComp("NuWro","C_4GeV",3850,4100));
  // SaveAndPrint(ERecRatioComp("NuWro","C"));
  // SaveAndPrint(ERecRatioComp("NuWro","C_1GeV"));
  // SaveAndPrint(ERecRatioComp("NuWro","C_2GeV"));
  // SaveAndPrint(ERecRatioComp("NuWro","C_4GeV"));

  // SaveAndPrint(ERec_pmu_pp());
  // SaveAndPrint(ERec_pmu_pp("_1GeV"));
  // SaveAndPrint(ERec_pmu_pp("_2GeV"));
  // SaveAndPrint(ERec_pmu_pp("_4GeV"));

  // SaveAndPrint(pHMP_pMu_FS_2D("NEUT","C"));
  // SaveAndPrint(pHMP_pMu_FS_2D("NEUT","C_norescat"));

  // SaveAndPrint(pHMP_pfm_FS_2D("NEUT","C"));
  // SaveAndPrint(pHMP_pfm_FS_2D("NEUT","C_norescat"));

  // SaveAndPrint(ERecSigmaVsETrueCorr("NEUT","C"));

  // SaveAndPrint(DPhiTVsDAlphaT("NEUT","C"));
  // SaveAndPrint(DPhiTVsDAlphaT("NEUT","C_norescat"));
  // SaveAndPrint(DPhiTVsDAlphaT("GENIE","C"));
  // SaveAndPrint(DPhiTVsDAlphaT("GENIE","C_norescat"));

  // SaveAndPrint(DPhiTVsDAlphaT("NuWro","C"));
  // SaveAndPrint(DPhiTVsDAlphaT("NuWro","C_norescat"));

  // SaveAndPrint(DPhiTVsDAlphaT("GiBUU","C"));
  // SaveAndPrint(DPhiTVsDAlphaT("GiBUU","C_norescat"));

  // SaveAndPrint(DPTVsStruckNucleonPt("NEUT","C"));
  // SaveAndPrint(DPTVsStruckNucleonPt("GENIE","C"));
  // SaveAndPrint(DPTVsStruckNucleonPt("NuWro","C"));
  // SaveAndPrint(DPTVsStruckNucleonPt("NEUT","C_norescat"));
  // SaveAndPrint(DrawDistrib("NEUT","C_norescat","other","HMProtonStdHepPos", "Count",
  //   "StdHepPos"));

  // SaveAndPrint(DPTFirstPVsStruckNucleonPt("NEUT","C_norescat"));

  // SaveAndPrint(DPTHMPvsStruckNucleonPtNE("NEUT","C_norescat"));
  // SaveAndPrint(DrawDistrib("NEUT","C_norescat","other","HMProtonStdHepPosNE", "Count",
  //   "StdHepPos"));
  // SaveAndPrint(DrawDistrib("NEUT","C_norescat","other","CCQEHMProtonMom_NE",
  //   "FS Proton Momentum [MeV/C]", "Count"));


  // SaveAndPrint(DPTHMPvsStruckNucleonPtNE_excludediagonal("NEUT","C_norescat"));
  // SaveAndPrint(DrawDistrib("NEUT","C_norescat","other",
  //   "HMProtonStdHepPos_excludediagonal",
  //   "StdHepPos", "Count"));
  // SaveAndPrint(DrawDistrib("NEUT","C_norescat","other",
  //   "CCQEHMProtonMom_excludediagonal",
  //   "FS Proton Momentum [MeV/C]", "Count"));

  // SaveAndPrint(DPTNoNEVsStruckNucleonPt("NEUT","C_norescat"));
  // SaveAndPrint(DPTVsStruckNucleonPt("GENIE","C_norescat"));
  // SaveAndPrint(DPTVsStruckNucleonPt("NuWro","C_norescat"));

  // SaveAndPrint(NuWroModel_DPTVsDAlphaT());
  // SaveAndPrint(NuWroModel_DPTVsDAlphaT(true));

  // SaveAndPrint(OmegaHMPvsDAlphaT());

  // SaveAndPrint(NProtonsVsCCQE("NEUT","C","DphiT"));
  // SaveAndPrint(NProtonsVsCCQE("GENIE","C","DphiT"));
  // SaveAndPrint(NProtonsVsCCQE("NuWro","C","DphiT"));
  // SaveAndPrint(NProtonsVsCCQE("GiBUU","C","DphiT"));

  // SaveAndPrint(NProtonsVsCCQE("NEUT","C","MuPt"));
  // SaveAndPrint(NProtonsVsCCQE("GENIE","C","MuPt"));
  // SaveAndPrint(NProtonsVsCCQE("NuWro","C","MuPt"));
  // SaveAndPrint(NProtonsVsCCQE("GiBUU","C","MuPt"));

  // SaveAndPrint(NProtonsVsCCQE("NEUT","C","DPt",true));
  // SaveAndPrint(NProtonsVsCCQE("GENIE","C","DPt",true));
  // SaveAndPrint(NProtonsVsCCQE("NuWro","C","DPt",true));
  // SaveAndPrint(NProtonsVsCCQE("GiBUU","C","DPt",true));

  // SaveAndPrint(LabMomVsStruckNucleonMom("NEUT","C_norescat"));
  // SaveAndPrint(LabMomVsStruckNucleonMom("GENIE","C_norescat"));


  // SaveAndPrint(DrawDistrib("NEUT","C","","LabMom", "","",true, true));
  // SaveAndPrint(DrawDistrib("NEUT","C_norescat","","LabMom", "","",true, true));
  return 0;
}


namespace {

void SetOpts(){
  CLIArgs::OptSpec.emplace_back("-h","--help", false,
    [&] (std::string const &opt) -> bool {
      CLIArgs::SayRunLike();
      exit(0);
    });

  CLIArgs::OptSpec.emplace_back("-Sa", "--Samples-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading samples from file : " << opt << std::endl;
      SamplesXML = opt;
      return true;
    }, true,[](){},"<Samples.xml>");

  CLIArgs::OptSpec.emplace_back("-Se", "--Selections-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading selections from file : " << opt << std::endl;
      SelectionsXML = opt;
      return true;
    }, false,[](){SelectionsXML = SamplesXML;},
    "<Selections.xml {Default: Value of -Sa}>");

  CLIArgs::OptSpec.emplace_back("-P", "--Plots-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading plots from file : " << opt << std::endl;
      PlotsXML = opt;
      MakeAutoPlots = true;
      return true;
    }, false,[](){PlotsXML = SelectionsXML;},
    "<Plots.xml {Default: Value of -Se}>");

  CLIArgs::OptSpec.emplace_back("-o", "--output-pdf", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Writing special plots to : " << opt << std::endl;
      OutputPDFFileName = opt;
      MakeSpecialPlots = true;
      return true;
    }, false,[](){MakeSpecialPlots = false; OutputPDFFileName="output-plots.pdf";},
    "<output-plots.pdf {Default:output-plots.pdf}>");

    CLIArgs::OptSpec.emplace_back("-C", "--Canvas-File", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Writing TCanvases to: " << opt << std::endl;
      OutputTCanvasFileName = opt;
      return true;
    }, false,[](){OutputTCanvasFileName = "TCanvi.root";},
    "<TFileName.root {Default: TCanvi.root}>");

    CLIArgs::OptSpec.emplace_back("-x", "--Cache-File", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Using " << opt << " as histogram Cache." <<  std::endl;
      HistogramCacheFileName = opt;
      return true;
    }, false,[](){HistogramCacheFileName = "CachedHistograms.root";},
    "<TFileName.root {Default: CachedHistograms.root}>");

    CLIArgs::OptSpec.emplace_back("-H", "--PreCache", false,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Pre-caching all the histograms " <<  std::endl;
      PreCache = true;
      return true;
    }, false,[](){PreCache = false;},
    "{Default: false} N.B. This may take an obscenely long time.");
  }
}

int main(int argc, char const * argv[]){
  SetOpts();

  CLIArgs::AddArguments(argc,argv);
  if(!CLIArgs::GetOpts()){
    CLIArgs::SayRunLike();
    return 1;
  }

  (void)Verbosity;
  if(!MakeAutoPlots && !MakeSpecialPlots && !PreCache){
    std::cout << "[INFO]: I haven't actually been asked to do anything."
      << std::endl;
    return 0;
  }

  TH1F::SetDefaultSumw2();
  if(!PlottingIO::InitialiseHistogramCache(HistogramCacheFileName.c_str(),
    SamplesXML, SelectionsXML, PreCache)){
    return 1;
  }

  InitOutputCanvasFile();

  if(MakeAutoPlots){
    try {
      std::cout << "[INFO]: Making XML configured plots." << std::endl;
      AutoPlots::MakeAutoPlotsXML(PlotsXML);
    } catch (EInvalidPlot e) {
      std::cout << "Caught Invalid Plot exception. Leaving." << std::endl;
    }
  }

  if(MakeSpecialPlots){
    try {
      std::cout << "[INFO]: Making special plots." << std::endl;
      MakeSomePlotzz();
    } catch (EInvalidPlot e) {
      std::cout << "Caught Invalid Plot exception. Leaving." << std::endl;
    }
  }

  CloseOutputCanvasFile();

}
