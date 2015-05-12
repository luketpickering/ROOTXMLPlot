#include <utility>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TCut.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TVectorD.h"

#include "TRandom.h"

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"
#include "PlottingIO.hxx"

using namespace PlottingDefaults;
using namespace PlottingUtils;

namespace {

  TFile* OutputFile;
  char const * PDFPrintName;

}

TH1F* DummyRatioGuide(int color=kBlack,
 float xMin=DataSpecifics::DPhiTMin, float xMax=DataSpecifics::DPhiTMax){
  TH1F* dummy = new TH1F("dummy","",1,xMin, xMax);
  dummy->Fill(1);
  dummy->SetLineColor(color);
  dummy->SetLineWidth(2);
  dummy->SetDirectory(0);

  return dummy;
}

TCanvas* SaveAndPrint(TCanvas* canv){
  std::cout << "Saving TCanvas: " << canv->GetName() << std::endl;
  OutputFile->WriteTObject(canv);
  canv->Print(PDFPrintName);
  return canv;
}

void PrintSampleStats(){
  if(Data::Trees.size() != DataSpecifics::Generators.size()){
    std::cerr << "Don't Appear to have the original trees in memory. "
      "Not printing Stats." << std::endl;
    return;
  } else {
    std::cout << "**************\nSample stats:\n**************" << std::endl;
  }
  for(std::vector<PlottingTypes::Generator>::iterator gen =
    DataSpecifics::Generators.begin();
    gen != DataSpecifics::Generators.end(); ++gen){

    for(std::vector<PlottingTypes::Target>::iterator tar = gen->Targets.begin();
      tar != gen->Targets.end(); ++tar){

      TTree &gen_nuc_tree = (*Data::Trees[gen->Name][tar->Name]);

      int Num_CCQE_Inc =
        gen_nuc_tree.GetEntries(DataSpecifics::Selections1D["Inc"].Cut);
      int Num_CCQE_NoNucEmit =
        gen_nuc_tree.GetEntries(DataSpecifics::Selections1D["NoNE"].Cut);
      int Num_CCQE_NucEmit =
        gen_nuc_tree.GetEntries(DataSpecifics::Selections1D["NE"].Cut);
      int Num_CCQE_ZeroMomProt =
        gen_nuc_tree.GetEntries(DataSpecifics::MuonValidSelect &&
          DataSpecifics::HighestMomProtonValidSelect &&
          (!DataSpecifics::HighestMomProtonNonZeroSelect));

      std::cout << "\t" << tar->Name << ":" << std::endl;
      std::cout << "\t\tNum_CCQE_Inc: "
        << Num_CCQE_Inc << std::endl;
      std::cout << "\t\tNum_CCQE_NoNucEmit: "
        << Num_CCQE_NoNucEmit << std::endl;
      std::cout << "\t\tNum_CCQE_NucEmit: "
        << Num_CCQE_NucEmit << std::endl;
      std::cout << "\t\tEmission Fraction: "
        << (Num_CCQE_NucEmit/float(Num_CCQE_Inc)) << std::endl;
      std::cout << "\t\t(No Emission Fraction): "
        << (Num_CCQE_NoNucEmit/float(Num_CCQE_Inc)) << std::endl;
      std::cout << "\t\tNum_CCQE_ZeroMomProt: "
        << Num_CCQE_ZeroMomProt << std::endl;

        std::cout << "\n" << std::endl;
    }
    std::cout << "\n\n" << std::endl;
  }
}

void SetSeriesStyle(TH1F* hist, int color, int linestyle){
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

  TH1F* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][ModeString]);
  MakePDF(NEUT_C_Hist);

  TH1F* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"][ModeString]);
  MakePDF(NEUT_O_Hist);

  TH1F* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][ModeString]);
  MakePDF(NEUT_Pb_Hist);

  TH1F* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][ModeString]);
  MakePDF(GENIE_C_Hist);

  TH1F* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"][ModeString]);
  MakePDF(GENIE_O_Hist);

  TH1F* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][ModeString]);
  MakePDF(GENIE_Pb_Hist);

  canv->cd(1);

  std::vector<double> YVals;

  YVals.push_back(GetMaximumBinContents(NEUT_C_Hist));
  YVals.push_back(GetMaximumBinContents(NEUT_O_Hist));
  YVals.push_back(GetMaximumBinContents(NEUT_Pb_Hist));
  YVals.push_back(GetMaximumBinContents(GENIE_C_Hist));
  YVals.push_back(GetMaximumBinContents(GENIE_O_Hist));
  YVals.push_back(GetMaximumBinContents(GENIE_Pb_Hist));

  double maxY = (*max_element(YVals.begin(),YVals.end()));
  YVals.clear();

  NEUT_C_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  NEUT_O_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  NEUT_Pb_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_C_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_O_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_Pb_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);

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

  TH1F* DRG = DummyRatioGuide();

  TH1F* PB_C_Hist_Ratio_NEUT = ScaleToRef(NEUT_Pb_Hist,NEUT_C_Hist).first;
  SetLowerPaneHistoDefaults(PB_C_Hist_Ratio_NEUT,"Ratio (Scaled To Denominator)",
   DPhiTTitle);
  PB_C_Hist_Ratio_NEUT->Draw("EHIST");

  TH1F* O_C_Hist_Ratio_NEUT = ScaleToRef(NEUT_O_Hist, NEUT_C_Hist).first;
  O_C_Hist_Ratio_NEUT->Draw("EHIST SAME");

  TH1F* PB_C_Hist_Ratio_GENIE = ScaleToRef(GENIE_Pb_Hist, GENIE_C_Hist).first;
  PB_C_Hist_Ratio_GENIE->Draw("EHIST SAME");

  TH1F* O_C_Hist_Ratio_GENIE = ScaleToRef(GENIE_O_Hist, GENIE_C_Hist).first;
  O_C_Hist_Ratio_GENIE->Draw("EHIST SAME");

  YVals.push_back(GetMaximumBinContents(PB_C_Hist_Ratio_NEUT));
  YVals.push_back(GetMaximumBinContents(PB_C_Hist_Ratio_NEUT));
  YVals.push_back(GetMaximumBinContents(O_C_Hist_Ratio_NEUT));
  YVals.push_back(GetMaximumBinContents(PB_C_Hist_Ratio_GENIE));
  YVals.push_back(GetMaximumBinContents(O_C_Hist_Ratio_GENIE));

  maxY = (*max_element(YVals.begin(),YVals.end()));
  YVals.clear();

  PB_C_Hist_Ratio_NEUT->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  PB_C_Hist_Ratio_NEUT->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  O_C_Hist_Ratio_NEUT->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  PB_C_Hist_Ratio_GENIE->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  O_C_Hist_Ratio_GENIE->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);

  DRG->Draw("HIST SAME");

  leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  leg->AddEntry(PB_C_Hist_Ratio_NEUT, "Pb/C (NEUT)","l");
  leg->AddEntry(PB_C_Hist_Ratio_GENIE, "Pb/C (GENIE)","l");
  leg->AddEntry(O_C_Hist_Ratio_NEUT, "O/C (NEUT)","l");
  leg->AddEntry(O_C_Hist_Ratio_GENIE, "O/C (GENIE)","l");

  leg->Draw();

  canv->cd(0);

  TitleFactory(TitleString.c_str())->Draw();

  return canv;
}

TCanvas* SpeciesCanvas(char const * species){

  std::string CanvString = species;

  TCanvas* canv = CanvasFactory(CanvString.c_str(), false);

  std::cout << "Trying to clone: Data::Histos[\"NEUT\"][\"C\"]["
    << species << "] == " << Data::Histos["NEUT"]["C"][species] << std::endl;

  TH1F* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][species]);
  MakePDF(NEUT_C_Hist);

  TH1F* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"][species]);
  MakePDF(NEUT_O_Hist);

  TH1F* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][species]);
  MakePDF(NEUT_Pb_Hist);

  TH1F* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][species]);
  MakePDF(GENIE_C_Hist);

  TH1F* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"][species]);
  MakePDF(GENIE_O_Hist);

  TH1F* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][species]);
  MakePDF(GENIE_Pb_Hist);

  // NEUT_C_Hist->GetYaxis()->SetRangeUser(1E-2,1.05);

  NEUT_C_Hist->GetXaxis()->SetTitle(species);
  NEUT_C_Hist->GetYaxis()->SetTitle(("#frac{dN}{Nd("+std::string(species)+")}").c_str());
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

//Carbon************************************************************************
  TH1F* C_Numer = NoSaveClone(FractionNumerIsEmit ?
      Data::Histos[generator]["C"]["NE"] : Data::Histos[generator]["C"]["NoNE"]);
  TH1F* C_Denom =  NoSaveClone(Data::Histos[generator]["C"]["Inc"]);

  TH1F* C_Ratio = ScaleToRef(C_Numer,C_Denom).first;

  C_Ratio->GetYaxis()->SetRangeUser(0,1.05);
  C_Ratio->GetYaxis()->SetTitle( (FractionNumerIsEmit ?
    "Nuclear Emission Fraction" : "No Nuclear Emission Fraction") );

  SetOnePaneHistDefaults(C_Ratio);
  C_Ratio->Draw("EHIST");

//Oxygen************************************************************************

  TH1F* O_Numer = NoSaveClone(FractionNumerIsEmit ?
      Data::Histos[generator]["O"]["NE"] : Data::Histos[generator]["O"]["NoNE"]);
  // BinWidthNormalise(O_Numer);
  TH1F* O_Denom = NoSaveClone(Data::Histos[generator]["O"]["Inc"]);
  // BinWidthNormalise(O_Denom);

  TH1F* O_Ratio = ScaleToRef(O_Numer,O_Denom).first;

  O_Ratio->GetYaxis()->SetRangeUser(0,1.05);
  O_Ratio->GetYaxis()->SetTitle( (FractionNumerIsEmit ?
    "Nuclear Emission Fraction" : "No Nuclear Emission Fraction") );

  SetOnePaneHistDefaults(O_Ratio);
  O_Ratio->Draw("EHIST SAME");

//Lead**************************************************************************

  TH1F* Pb_Numer = NoSaveClone(FractionNumerIsEmit ?
      Data::Histos[generator]["Pb"]["NE"] : Data::Histos[generator]["Pb"]["NoNE"]);
  // BinWidthNormalise(Pb_Numer);
  TH1F* Pb_Denom = NoSaveClone(Data::Histos[generator]["Pb"]["Inc"]);
  // BinWidthNormalise(Pb_Denom);

  TH1F* Pb_Ratio = ScaleToRef(Pb_Numer,Pb_Denom).first;

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

  TH1F* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"][ModeString]);
  MakePDF(NEUT_C_Hist);

  TH1F* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"][ModeString]);
  MakePDF(NEUT_O_Hist);

  TH1F* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"][ModeString]);
  MakePDF(NEUT_Pb_Hist);

  TH1F* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"][ModeString]);
  MakePDF(GENIE_C_Hist);

  TH1F* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"][ModeString]);
  MakePDF(GENIE_O_Hist);

  TH1F* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"][ModeString]);
  MakePDF(GENIE_Pb_Hist);

  canv->cd(1);

  std::vector<double> YVals;

  YVals.push_back(GetMaximumBinContents(NEUT_C_Hist));
  YVals.push_back(GetMaximumBinContents(NEUT_O_Hist));
  YVals.push_back(GetMaximumBinContents(NEUT_Pb_Hist));
  YVals.push_back(GetMaximumBinContents(GENIE_C_Hist));
  YVals.push_back(GetMaximumBinContents(GENIE_O_Hist));
  YVals.push_back(GetMaximumBinContents(GENIE_Pb_Hist));

  double maxY = (*max_element(YVals.begin(),YVals.end()));
  YVals.clear();

  NEUT_C_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  NEUT_O_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  NEUT_Pb_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_C_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_O_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_Pb_Hist->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);

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

  TH1F* Pb_GENIE_NEUT_Ratio = ScaleToRef(GENIE_Pb_Hist, NEUT_Pb_Hist).first;
  TH1F* O_GENIE_NEUT_Ratio = ScaleToRef(GENIE_O_Hist, NEUT_O_Hist).first;
  TH1F* C_GENIE_NEUT_Ratio = ScaleToRef(GENIE_C_Hist, NEUT_C_Hist).first;

  YVals.push_back(GetMaximumBinContents(Pb_GENIE_NEUT_Ratio));
  YVals.push_back(GetMaximumBinContents(O_GENIE_NEUT_Ratio));
  YVals.push_back(GetMaximumBinContents(C_GENIE_NEUT_Ratio));

  maxY = (*max_element(YVals.begin(),YVals.end()));
  YVals.clear();

  Pb_GENIE_NEUT_Ratio->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  O_GENIE_NEUT_Ratio->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  C_GENIE_NEUT_Ratio->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);

  SetLowerPaneHistoDefaults(Pb_GENIE_NEUT_Ratio, "Ratio (Scaled To Denominator)",
   DPhiTTitle);

  SetSeriesStyle(Pb_GENIE_NEUT_Ratio,kGreen,1);
  Pb_GENIE_NEUT_Ratio->Draw("EHIST");

  SetSeriesStyle(O_GENIE_NEUT_Ratio, kRed,1);
  O_GENIE_NEUT_Ratio->Draw("EHIST SAME");

  SetSeriesStyle(C_GENIE_NEUT_Ratio, kBlue,1);
  C_GENIE_NEUT_Ratio->Draw("EHIST SAME");

  TH1F* DRG = DummyRatioGuide();
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

  TCanvas* canv = CanvasFactory((std::string(Target) + "ProtonMomentumCanv").c_str(), false);

  TH1F* NEUT_T_First = NoSaveClone(Data::Histos["NEUT"][Target]["FirstPMom"]);
  MakePDF(NEUT_T_First); SetSeriesStyle(NEUT_T_First,kRed,1);

  TH1F* NEUT_T_Highest = NoSaveClone(Data::Histos["NEUT"][Target]["HighestPMom"]);
  MakePDF(NEUT_T_Highest); SetSeriesStyle(NEUT_T_Highest,kBlue,1);

  TH1F* GENIE_T_First = NoSaveClone(Data::Histos["GENIE"][Target]["FirstPMom"]);
  MakePDF(GENIE_T_First); SetSeriesStyle(GENIE_T_First,kRed,2);

  TH1F* GENIE_T_Highest = NoSaveClone(Data::Histos["GENIE"][Target]["HighestPMom"]);
  MakePDF(GENIE_T_Highest); SetSeriesStyle(GENIE_T_Highest,kBlue,2);

  NEUT_T_First->GetXaxis()->SetTitle("Proton Momentum [MeV/C]");
  NEUT_T_First->GetXaxis()->SetTitleOffset(XTitleOffset);
  NEUT_T_First->GetXaxis()->SetLabelOffset(XLabelOffset);

  std::vector<double> YVals;

  YVals.push_back(GetMaximumBinContents(NEUT_T_First));
  YVals.push_back(GetMaximumBinContents(NEUT_T_Highest));
  YVals.push_back(GetMaximumBinContents(GENIE_T_First));
  YVals.push_back(GetMaximumBinContents(GENIE_T_Highest));


  double maxY = (*max_element(YVals.begin(),YVals.end()));
  YVals.clear();

  NEUT_T_First->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  NEUT_T_Highest->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_T_First->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);
  GENIE_T_Highest->GetYaxis()->SetRangeUser(MinYVal,maxY*1.2);


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

  TH1F* T_NoNE_Hist = NoSaveClone(Data::Histos[generator][Target]["NoNE"]);
  MakePDF(T_NoNE_Hist); SetSeriesStyle(T_NoNE_Hist,kGreen-3,1);

  TH1F* T_NE_Hist = NoSaveClone(Data::Histos[generator][Target]["NE"]);
  MakePDF(T_NE_Hist); SetSeriesStyle(T_NE_Hist,kCyan,1);

  canv->cd(1);

  std::vector<double> YVals;

  YVals.push_back(GetMaximumBinContents(T_NoNE_Hist));
  YVals.push_back(GetMaximumBinContents(T_NE_Hist));

  double maxY = (*max_element(YVals.begin(),YVals.end()));
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

  TH1F* T_NE_NoNE_Ratio = ScaleToRef(T_NE_Hist, T_NoNE_Hist).first;
  SetSeriesStyle(T_NE_NoNE_Ratio,kGreen-3,1);

  SetLowerPaneHistoDefaults(T_NE_NoNE_Ratio, "Ratio (Scaled To Denominator)",
    DPhiTTitle);

  T_NE_NoNE_Ratio->Draw("EHIST");

  TH1F* DRG = DummyRatioGuide();
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

TCanvas* LinearInset(){
  TCanvas* canv = CanvasFactory("LinearInset");

  TH1F* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"]["NoNE"]);
  MakePDF(GENIE_O_Hist);

  TH1F* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"]["NoNE"]);
  MakePDF(GENIE_C_Hist);

  TH1F* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"]["NoNE"]);
  MakePDF(GENIE_Pb_Hist);

  TH1F* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"]["NoNE"]);
  MakePDF(NEUT_O_Hist);

  TH1F* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["NoNE"]);
  MakePDF(NEUT_C_Hist);

  TH1F* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"]["NoNE"]);
  MakePDF(NEUT_Pb_Hist);

  SetOnePaneHistDefaults(GENIE_Pb_Hist);
  GENIE_Pb_Hist->GetYaxis()->SetTitleOffset(YTitleOffset*1.1);
  GENIE_Pb_Hist->GetYaxis()->SetTitle(DPhiT_PDF_Title);

  GENIE_Pb_Hist->Draw("EHIST");
  GENIE_C_Hist->Draw("EHIST SAME");
  GENIE_O_Hist->Draw("EHIST SAME");
  NEUT_Pb_Hist->Draw("EHIST SAME");
  NEUT_C_Hist->Draw("EHIST SAME");
  NEUT_O_Hist->Draw("EHIST SAME");

  TLegend* leg = LegendFactory(
    (1.0 - P0LegendWidth - (Pad0XMargin*1.05)),
    (Pad0YMargin*1.2),
    (1.0 - (Pad0XMargin*1.05)),
    ((Pad0YMargin*1.2) + (P0LegendHeight*0.7) ));

  leg->AddEntry(GENIE_O_Hist,"O (GENIE)","l");
  leg->AddEntry(GENIE_C_Hist,"C (GENIE)","l");
  leg->AddEntry(GENIE_Pb_Hist,"Pb (GENIE)","l");
  leg->AddEntry(NEUT_O_Hist,"O (NEUT)","l");
  leg->AddEntry(NEUT_C_Hist,"C (NEUT)","l");
  leg->AddEntry(NEUT_Pb_Hist,"Pb (NEUT)","l");

  double InsetWidth = 0.6;
  double InsetHeight = 0.55;
  // leg->Draw();

  TPad *subpad = new TPad("subpad","",
    (1.0 - InsetWidth - (Pad0XMargin*1.1)),
    (1.0 - InsetHeight - (Pad0YMargin*1.2)),
    (1.0 - (Pad0XMargin*1.1)),
    (1.0 - (Pad0YMargin*1.2)));
  subpad->Draw();
  subpad->cd();
  subpad->SetRightMargin(0.05);


  TH1F* GENIE_O_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["O"]["NoNEZoom"]);
  MakePDF(GENIE_O_Hist_Zoom);

  TH1F* GENIE_C_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["C"]["NoNEZoom"]);
  MakePDF(GENIE_C_Hist_Zoom);

  TH1F* GENIE_Pb_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["Pb"]["NoNEZoom"]);
  MakePDF(GENIE_Pb_Hist_Zoom);

  TH1F* NEUT_O_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["O"]["NoNEZoom"]);
  MakePDF(NEUT_O_Hist_Zoom);

  TH1F* NEUT_C_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["C"]["NoNEZoom"]);
  MakePDF(NEUT_C_Hist_Zoom);

  TH1F* NEUT_Pb_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["Pb"]["NoNEZoom"]);
  MakePDF(NEUT_Pb_Hist_Zoom);

  SetOnePaneHistDefaults(GENIE_Pb_Hist_Zoom);
  GENIE_Pb_Hist_Zoom->GetYaxis()->SetTitleOffset(YTitleOffset*1.1);
  GENIE_Pb_Hist_Zoom->GetYaxis()->SetTitle(DPhiT_PDF_Title);
  GENIE_Pb_Hist_Zoom->GetYaxis()->SetRangeUser(0,
    GetMaximumBinContents(GENIE_Pb_Hist_Zoom)*1.1);

  GENIE_Pb_Hist_Zoom->Draw("EHIST");
  GENIE_C_Hist_Zoom->Draw("EHIST SAME");
  GENIE_O_Hist_Zoom->Draw("EHIST SAME");
  NEUT_Pb_Hist_Zoom->Draw("EHIST SAME");
  NEUT_C_Hist_Zoom->Draw("EHIST SAME");
  NEUT_O_Hist_Zoom->Draw("EHIST SAME");

  canv->cd();
  leg->Draw();

  return canv;
}

TCanvas* LinearInsetTest(){

  TH1F* GENIE_O_Hist = NoSaveClone(Data::Histos["GENIE"]["O"]["NoNE"]);
  MakePDF(GENIE_O_Hist);

  TH1F* GENIE_C_Hist = NoSaveClone(Data::Histos["GENIE"]["C"]["NoNE"]);
  MakePDF(GENIE_C_Hist);

  TH1F* GENIE_Pb_Hist = NoSaveClone(Data::Histos["GENIE"]["Pb"]["NoNE"]);
  MakePDF(GENIE_Pb_Hist);

  TH1F* NEUT_O_Hist = NoSaveClone(Data::Histos["NEUT"]["O"]["NoNE"]);
  MakePDF(NEUT_O_Hist);

  TH1F* NEUT_C_Hist = NoSaveClone(Data::Histos["NEUT"]["C"]["NoNE"]);
  MakePDF(NEUT_C_Hist);

  TH1F* NEUT_Pb_Hist = NoSaveClone(Data::Histos["NEUT"]["Pb"]["NoNE"]);
  MakePDF(NEUT_Pb_Hist);

  SetOnePaneHistDefaults(GENIE_Pb_Hist);
  GENIE_Pb_Hist->GetYaxis()->SetTitleOffset(YTitleOffset*1.1);
  GENIE_Pb_Hist->GetYaxis()->SetTitle(DPhiT_PDF_Title);

  GENIE_O_Hist->SetTitle("O (GENIE)");
  GENIE_C_Hist->SetTitle("C (GENIE)");
  GENIE_Pb_Hist->SetTitle("Pb (GENIE)");
  NEUT_O_Hist->SetTitle("O (NEUT)");
  NEUT_C_Hist->SetTitle("C (NEUT)");
  NEUT_Pb_Hist->SetTitle("Pb (NEUT)");

  TH1F* GENIE_O_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["O"]["NoNEZoom"]);
  MakePDF(GENIE_O_Hist_Zoom);

  TH1F* GENIE_C_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["C"]["NoNEZoom"]);
  MakePDF(GENIE_C_Hist_Zoom);

  TH1F* GENIE_Pb_Hist_Zoom = NoSaveClone(Data::Histos["GENIE"]["Pb"]["NoNEZoom"]);
  MakePDF(GENIE_Pb_Hist_Zoom);

  TH1F* NEUT_O_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["O"]["NoNEZoom"]);
  MakePDF(NEUT_O_Hist_Zoom);

  TH1F* NEUT_C_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["C"]["NoNEZoom"]);
  MakePDF(NEUT_C_Hist_Zoom);

  TH1F* NEUT_Pb_Hist_Zoom = NoSaveClone(Data::Histos["NEUT"]["Pb"]["NoNEZoom"]);
  MakePDF(NEUT_Pb_Hist_Zoom);

  SetOnePaneHistDefaults(GENIE_Pb_Hist_Zoom);
  GENIE_Pb_Hist_Zoom->GetYaxis()->SetTitleOffset(YTitleOffset*1.1);
  GENIE_Pb_Hist_Zoom->GetYaxis()->SetTitle(DPhiT_PDF_Title);
  GENIE_Pb_Hist_Zoom->GetYaxis()->SetRangeUser(0,
    GetMaximumBinContents(GENIE_Pb_Hist_Zoom)*1.1);

  std::vector<TH1F*> BigPlots;
  BigPlots.push_back(GENIE_Pb_Hist);
  BigPlots.push_back(GENIE_C_Hist);
  BigPlots.push_back(GENIE_O_Hist);
  BigPlots.push_back(NEUT_O_Hist);
  BigPlots.push_back(NEUT_C_Hist);
  BigPlots.push_back(NEUT_Pb_Hist);
  std::vector<TH1F*> InsetPlots;
  InsetPlots.push_back(GENIE_Pb_Hist_Zoom);
  InsetPlots.push_back(GENIE_C_Hist_Zoom);
  InsetPlots.push_back(GENIE_O_Hist_Zoom);
  InsetPlots.push_back(NEUT_O_Hist_Zoom);
  InsetPlots.push_back(NEUT_C_Hist_Zoom);
  InsetPlots.push_back(NEUT_Pb_Hist_Zoom);

  return GetInsetCanvas(BigPlots,InsetPlots);
}

TCanvas* GaussianExampleCanvas(){
  TCanvas* canv = SplitCanvasFactory("GaussianExampleCanvas");
  TRandom gaus;

  TH1F* HThin = new TH1F("HThin","",100,0,30);
  TH1F* HFat = new TH1F("HFat","",100,0,30);

  HThin->SetDirectory(0);
  HFat->SetDirectory(0);

  for(long i = 0; i < 1E6L; ++i){
    HThin->Fill(gaus.Gaus(0.0,5.0));
  }
  MakePDF(HThin);

  for(long i = 0; i < 1E6L; ++i){
    HFat->Fill(gaus.Gaus(0.0,12.0));
  }
  MakePDF(HFat);

  canv->cd(1);

  SetUpperPaneHistoDefaults(HFat, "#frac{dN}{NdArb}}");

  TLegend* leg = LegendFactory( (1.0 - P1LegendWidth - Pad1XMargin),
    (1.0 - P1LegendHeight*0.3 - Pad1YMargin),
    (1.0 - Pad1XMargin),
    (1.0 - Pad1YMargin));

  SetSeriesStyle(HThin,kRed-3,1);
  SetSeriesStyle(HFat,kBlue-2,1);

  leg->AddEntry(HThin,"Thinner","l");
  leg->AddEntry(HFat,"Wider","l");

  HThin->GetYaxis()->SetRangeUser(0,HThin->GetBinContent(1)*1.05);
  HFat->GetYaxis()->SetRangeUser(0,HThin->GetBinContent(1)*1.05);

  HFat->Draw("EHIST");
  HThin->Draw("EHIST SAME");

  leg->Draw();

  canv->cd(2);

  TH1F* Ratio = NoSaveClone(HFat);

  Ratio->Divide(HThin);

  SetLowerPaneHistoDefaults(Ratio, "Ratio (Scaled To Denominator)", "Arb");

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

void MakeSomePlotzz(char const * outputName,
  bool LoadHists = false,
  char const * InputFileName = "") {


  PDFPrintName = outputName;
  TH1F::SetDefaultSumw2();

  TCanvas* dummyCanv = new TCanvas("dummyCanv","");
  dummyCanv->Update();dummyCanv->Print((std::string(outputName)+"[").c_str());

  if(!PlottingIO::LoadHistograms(LoadHists, InputFileName)){
    return;
  }

  OutputFile = new TFile("TR_AP_Generator_O_Pb_C.root","RECREATE");

  PrintSampleStats();

  SaveAndPrint(GeneratorDependence::EmissionFractions("NEUT", true));
  SaveAndPrint(GeneratorDependence::EmissionFractions("GENIE", true));
  SaveAndPrint(GeneratorDependence::EmissionFractions("NEUT", false));
  SaveAndPrint(GeneratorDependence::EmissionFractions("GENIE", false));
  SaveAndPrint(TargetDependence::TargetDeptDPTCanvas(true));
  SaveAndPrint(TargetDependence::TargetDeptDPTCanvas(false));
  SaveAndPrint(GeneratorDependence::GeneratorDeptDPTCanvas(true));
  SaveAndPrint(GeneratorDependence::GeneratorDeptDPTCanvas(false));

  SaveAndPrint(TargetDependence::SpeciesCanvas("NProt"));
  SaveAndPrint(TargetDependence::SpeciesCanvas("NGam"));
  SaveAndPrint(TargetDependence::SpeciesCanvas("NNeutron"));
  SaveAndPrint(TargetDependence::SpeciesCanvas("NPiPlus"));
  SaveAndPrint(TargetDependence::SpeciesCanvas("NPiZero"));
  SaveAndPrint(TargetDependence::SpeciesCanvas("NPiMin"));

  SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("C"));
  SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("O"));
  SaveAndPrint(GeneratorDependence::ProtonMomentumCanvas("Pb"));

  SaveAndPrint(ModeDependence::EOverNE("O","GENIE"));
  SaveAndPrint(LinearInset());
  SaveAndPrint(LinearInsetTest());
  SaveAndPrint(GaussianExampleCanvas());

//   /**************MuonPt 2D Plots******************/

//   TH2F* C_CCQE_NoNE_DPT_MuonPT_nt = new TH2F("C_CCQE_NoNE_DPT_MuonPT_nt","",
//     NXBINS,MINXBIN,PIPIPI,100,0,500);
//   BinLogX(C_CCQE_NoNE_DPT_MuonPT_nt->GetXaxis());
//   trees["NEUT"]["C"]->Draw("TransversityInfo.ModMuonPt:TransversityInfo.Delta"
//     "PhiT_HighestMomProton >> C_CCQE_NoNE_DPT_MuonPT_nt",
//     (BaseSelect+CCQECleanSelect).c_str());

//   TH2F* C_CCQE_NE_DPT_MuonPT_gn = new TH2F("C_CCQE_NE_DPT_MuonPT_gn","",
//     NXBINS,MINXBIN,PIPIPI,100,0,0.5);
//   BinLogX(C_CCQE_NE_DPT_MuonPT_gn->GetXaxis());
//   trees["GENIE"]["C"]->Draw("TransversityInfo.ModMuonPt:TransversityInfo.Delta"
//     "PhiT_HighestMomProton >> C_CCQE_NE_DPT_MuonPT_gn",
//     (BaseSelect+CCQECleanSelect).c_str());

//   c1->cd(0);
//   c1->GetPad(0)->SetRightMargin(0.2);
//   C_CCQE_NE_DPT_MuonPT_nt->GetYaxis()->SetTitle("Muon P_{t} [MeV/c]");
//   C_CCQE_NE_DPT_MuonPT_nt->GetXaxis()->SetTitle("#delta#phi_{T} [rad]");
//   C_CCQE_NE_DPT_MuonPT_nt->Draw("colz");
//   TLatex title_bla2(0.1,0.95,"CCQE No Nuclear Emission NEUT");
//   title_bla2.SetTextSize(0.04);
//   title_bla2.SetNDC();
//   title_bla2.Draw();
//   c1->Update(); c1->Print(outputName.c_str());

//   C_CCQE_NE_DPT_MuonPT_gn->GetYaxis()->SetTitle("Muon P_{t} [GeV/c]");
//   C_CCQE_NE_DPT_MuonPT_gn->GetXaxis()->SetTitle("#delta#phi_{T} [rad]");
//   C_CCQE_NE_DPT_MuonPT_gn->Draw("colz");
//   TLatex title_bla3(0.1,0.95,"CCQE No Nuclear Emission GENIE");
//   title_bla3.SetTextSize(0.04);
//   title_bla3.SetNDC();
//   title_bla3.Draw();
//   c1->Update(); c1->Print(outputName.c_str());

// std::cout << "blablalbal22222222" << std::endl;


// //   /**************MuonPt 2D Plots******************/

//   TH2F* C_CCQE_NE_DPT_MuonPT_nt = new TH2F("C_CCQE_NE_DPT_MuonPT_nt","",
//     NXBINS,MINXBIN,PIPIPI,100,0,500);
//   BinLogX(C_CCQE_E_DPT_MuonPT_nt->GetXaxis());
//   trees["NEUT"]["C"]->Draw("TransversityInfo.ModMuonPt:TransversityInfo.Delta"
//     "PhiT_HighestMomProton >> C_CCQE_E_DPT_MuonPT_nt",
//     (BaseSelect+CCQEExplodeSelect).c_str());

//   TH2F* C_CCQE_E_DPT_MuonPT_gn = new TH2F("C_CCQE_E_DPT_MuonPT_gn","",
//     NXBINS,MINXBIN,PIPIPI,100,0,0.5);
//   BinLogX(C_CCQE_E_DPT_MuonPT_gn->GetXaxis());
//   trees["GENIE"]["C"]->Draw("TransversityInfo.ModMuonPt:TransversityInfo.Delta"
//     "PhiT_HighestMomProton >> C_CCQE_E_DPT_MuonPT_gn",
//     (BaseSelect+CCQEExplodeSelect).c_str());

//   c1->cd(0);

//   C_CCQE_E_DPT_MuonPT_nt->GetYaxis()->SetTitle("Muon P_{t} [MeV/c]");
//   C_CCQE_E_DPT_MuonPT_nt->GetXaxis()->SetTitle("#delta#phi_{T} [rad]");

//   C_CCQE_E_DPT_MuonPT_nt->Draw("colz");

//   TLatex title_bla1(0.1,0.95,"CCQE Nuclear Emission NEUT");
//   title_bla1.SetTextSize(0.04);
//   title_bla1.SetNDC();
//   title_bla1.Draw();

//   c1->Update(); c1->Print(outputName.c_str());

//   C_CCQE_E_DPT_MuonPT_gn->GetYaxis()->SetTitle("Muon P_{t} [GeV/c]");
//   C_CCQE_E_DPT_MuonPT_gn->GetXaxis()->SetTitle("#delta#phi_{T} [rad]");

//   C_CCQE_E_DPT_MuonPT_gn->Draw("colz");

//   TLatex title_bla(0.1,0.95,"CCQE Nuclear Emission GENIE");
//   title_bla.SetTextSize(0.04);
//   title_bla.SetNDC();
//   title_bla.Draw();

//   c1->Update(); c1->Print(outputName.c_str());

  dummyCanv->Print((std::string(outputName)+"]").c_str());

  OutputFile->Write();
}
