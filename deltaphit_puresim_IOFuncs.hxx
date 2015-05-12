#ifndef PLOTTINGIO_HXX_SEEN
#define PLOTTINGIO_HXX_SEEN
#include <map>
#include <string>
#include <vector>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCut.h"

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"

namespace Data {
  std::map< std::string,
            std::map< std::string,
                      TFile*>
          > Files;
  std::map< std::string,
            std::map< std::string,
                      TTree*>
          > Trees;
  std::map< std::string,
            std::map< std::string,
                      std::map< std::string,
                                TH1F*
                              >
                    >
          > Histos;
}

namespace DataSpecifics {

  bool const DoLogx = PlottingDefaults::DoLogx;

  int const NProtonsMax = 50;
  int const NGammaMax = 5;
  int const NNeutronsMax = 50;
  int const NPiPlusMax = 5;
  int const NPiZeroMax = 5;
  int const NPiMinusMax = 5;

  int const ProtonMomBins = 100;
  double const ProtonMomMax = 1000;

  int const DPhiTBins = 100;
  double const DPhiTMin = ((!DoLogx)?0.0:1E-2);
  double const DPhiTMax = 3.15;

  int const DPhiTZoomxBins = 50;
  double const DPhiTZoomxMin = ((!DoLogx)?0.0:1E-5);
  double const DPhiTZoomxMax = 1E-1;

  int const C_Color = kRed;
  int const O_Color = kBlue;
  int const Pb_Color = kOrange;

  int const NEUT_Line = 1;
  int const GENIE_Line = 2;

  char const * VarTreeName = "TransversitudenessPureSim";
  char const * HistoOutputFileName =
    (DoLogx?"PureSimPlots_logx.root":"PureSimPlots_nologx.root");

  std::vector<PlottingTypes::Generator> Generators;
  std::map<std::string, PlottingTypes::Selection1D> Selections1D;

  TCut MuonValidSelect;
  TCut HighestMomProtonValidSelect;
  TCut HighestMomProtonNonZeroSelect;

  TCut BaseDPhiTSelect;
  TCut CCQESelect;

  TCut NucEmitSelect;
  TCut NoNucEmitSelect;

  TCut NPartsSelect;
  TCut ProtonMomentumSelect;

  void InitTCuts(){
    MuonValidSelect = "MuonPDG==13";
    HighestMomProtonValidSelect = "HighestMomProtonPDG==2212";
    HighestMomProtonNonZeroSelect = "ProtonMom_HighestMomProton>1E-8";

    BaseDPhiTSelect = MuonValidSelect && HighestMomProtonValidSelect &&
      HighestMomProtonNonZeroSelect;

    CCQESelect = "NeutReacCode==1";
    NoNucEmitSelect = (CCQESelect && "NProtons==1") && "NPrimTraj==2";
    NucEmitSelect = (CCQESelect && "NProtons!=1");
    NPartsSelect = MuonValidSelect && CCQESelect;
    ProtonMomentumSelect = MuonValidSelect && HighestMomProtonValidSelect &&
      CCQESelect;
  }
  ///You edit what things are loaded, and thus plottable, here.
  void SetSpecifics(){
    InitTCuts();
    //Add NEUT Targets
    Generators.push_back(PlottingTypes::Generator("NEUT", NEUT_Line));
    Generators.back().Targets.push_back(PlottingTypes::Target("C", C_Color,
      "../neut/C/nt_puresim_C_transversityVars.root"));
    Generators.back().Targets.push_back(PlottingTypes::Target("O", O_Color,
      "../neut/O/nt_puresim_O_transversityVars.root"));
    Generators.back().Targets.push_back(PlottingTypes::Target("Pb", Pb_Color,
      "../neut/Pb/nt_puresim_Pb_transversityVars.root"));

    //Add GENIE Targets
    Generators.push_back(PlottingTypes::Generator("GENIE", GENIE_Line));
    Generators.back().Targets.push_back(PlottingTypes::Target("C", C_Color,
      "../genie/C/gn_puresim_C_transversityVars.root"));
    Generators.back().Targets.push_back(PlottingTypes::Target("O", O_Color,
      "../genie/O/gn_puresim_O_transversityVars.root"));
    Generators.back().Targets.push_back(PlottingTypes::Target("Pb", Pb_Color,
      "../genie/Pb/gn_puresim_Pb_transversityVars.root"));

    //Add generic selections
    ///Delta Phi T Distribution Selections
    Selections1D["Inc"] =
      (PlottingTypes::Selection1D("Inclusive","Inc",
      "TransversityInfo.DeltaPhiT_HighestMomProton",
      DPhiTBins, DPhiTMin, DPhiTMax, BaseDPhiTSelect && CCQESelect, DoLogx));
    Selections1D["NE"] =
      (PlottingTypes::Selection1D("Nuclear Emission","NE",
      "TransversityInfo.DeltaPhiT_HighestMomProton",
      DPhiTBins, DPhiTMin, DPhiTMax, BaseDPhiTSelect && NucEmitSelect, DoLogx));
    Selections1D["NoNE"] =
      (PlottingTypes::Selection1D("No Nuclear Emission","NoNE",
      "TransversityInfo.DeltaPhiT_HighestMomProton",
      DPhiTBins, DPhiTMin, DPhiTMax, BaseDPhiTSelect && NoNucEmitSelect, DoLogx));
    Selections1D["NEZoom"] =
      (PlottingTypes::Selection1D("Nuclear Emission Zoom","NEZoom",
      "TransversityInfo.DeltaPhiT_HighestMomProton",
      DPhiTZoomxBins, DPhiTZoomxMin, DPhiTZoomxMax,
      BaseDPhiTSelect&&NucEmitSelect, DoLogx));
    Selections1D["NoNEZoom"] =
      (PlottingTypes::Selection1D("No Nuclear Emission Zoom","NoNEZoom",
      "TransversityInfo.DeltaPhiT_HighestMomProton",
      DPhiTZoomxBins, DPhiTZoomxMin, DPhiTZoomxMax,
      BaseDPhiTSelect && NoNucEmitSelect, DoLogx));

    ///Emission Species Distributions
    Selections1D["NProt"] =
      (PlottingTypes::Selection1D("NProtons","NProt",
      "TransversityInfo.NProtons",
      NProtonsMax, 0, NProtonsMax, NPartsSelect));
    Selections1D["NGam"] =
      (PlottingTypes::Selection1D("NGamma","NGam",
      "TransversityInfo.NGamma",
      NGammaMax, 0, NGammaMax, NPartsSelect));
    Selections1D["NNeutron"] =
      (PlottingTypes::Selection1D("NNeutrons","NNeutron",
      "TransversityInfo.NNeutrons",
      NNeutronsMax, 0, NNeutronsMax, NPartsSelect));
    Selections1D["NPiPlus"] =
      (PlottingTypes::Selection1D("NPiPlus","NPiPlus",
      "TransversityInfo.NPiPlus",
      NPiPlusMax, 0, NPiPlusMax, NPartsSelect));
    Selections1D["NPiZero"] =
      (PlottingTypes::Selection1D("NPiZero","NPiZero",
      "TransversityInfo.NPiZero",
      NPiZeroMax, 0, NPiZeroMax, NPartsSelect));
    Selections1D["NPiMin"] =
      (PlottingTypes::Selection1D("NPiMinus","NPiMin",
      "TransversityInfo.NPiMinus",
      NPiMinusMax, 0, NPiMinusMax, NPartsSelect));

    ///Proton Momenta
    Selections1D["FirstPMom"] =
      (PlottingTypes::Selection1D("FirstPMom","FirstPMom",
      "TransversityInfo.ProtonMom_FirstProton",
      ProtonMomBins, 0, ProtonMomMax, ProtonMomentumSelect));
    Selections1D["HighestPMom"] =
      (PlottingTypes::Selection1D("HighestPMom","HighestPMom",
      "TransversityInfo.ProtonMom_HighestMomProton",
      ProtonMomBins, 0, ProtonMomMax, ProtonMomentumSelect));
  }
}

namespace PlottingIO {

  bool LoadFilesAndTrees(){
    for(std::vector<PlottingTypes::Generator>::iterator gen = DataSpecifics::Generators.begin();
      gen != DataSpecifics::Generators.end(); ++gen){

      for(std::vector<PlottingTypes::Target>::iterator tar = gen->Targets.begin();
        tar != gen->Targets.end(); ++tar){

        if(!(Data::Files[gen->Name][tar->Name] =
            new TFile(tar->FileLocation.c_str()))){
          std::cerr << "Failed to load file: " << tar->FileLocation
            << " for Generator: " << gen->Name << ", Target: " << tar->Name
            << std::endl;
            return false;
        }
        if( !(Data::Trees[gen->Name][tar->Name] =
              dynamic_cast<TTree*>(
                Data::Files[gen->Name][tar->Name]->Get(DataSpecifics::VarTreeName)) ) ){
          std::cerr << "Failed to load TTree: for Generator: " << gen->Name
            << ", Target: " << tar->Name << std::endl;
            return false;
        }

      }
    }
    return true;
  }
  bool FillHistograms(){
    TFile* OutputFile = new TFile(DataSpecifics::HistoOutputFileName,"CREATE");
    if(!OutputFile->IsOpen()){
      std::cerr << "Cannot open " << DataSpecifics::HistoOutputFileName
        << " for writing." << std::endl;
      return false;
    }
    for(std::vector<PlottingTypes::Generator>::iterator gen = DataSpecifics::Generators.begin();
      gen != DataSpecifics::Generators.end(); ++gen){

      std::cout << "Filling Generator: " << gen->Name << std::endl;
      for(std::vector<PlottingTypes::Target>::iterator tar = gen->Targets.begin();
        tar != gen->Targets.end(); ++tar){
        std::cout << "\tFilling Target: " << tar->Name << std::endl;

        for(std::map<std::string, PlottingTypes::Selection1D>::iterator sel_it = DataSpecifics::Selections1D.begin();
          sel_it != DataSpecifics::Selections1D.end(); ++sel_it){

          PlottingTypes::Selection1D const &sel = sel_it->second;
          std::cout << "\t\tSel: " << sel.Name
            << "\n\t\t  Draw: \"" << sel.DrawString
            << "\"\n\t\t  Cut: \"" << sel.Cut
            <<"\"\n\t\t  Bins: [ " << sel.XBinLow << " .. " << sel.NBins << " .. "
            << sel.XBinUpper << " ]" << std::endl;

          std::string HistName = PlottingTypes::MakeHistoName(*gen, *tar, sel);
          ///Make the TH1F
          Data::Histos[gen->Name][tar->Name][sel.Tag] =
            new TH1F(HistName.c_str(),"",
              sel.NBins,sel.XBinLow,sel.XBinUpper);

          if(sel.UseLogX){
            PlottingUtils::BinLogX(Data::Histos[gen->Name][tar->Name][sel.Tag]->GetXaxis());
          }
          gen->SetSeriesStyle(
            tar->SetSeriesStyle(
              Data::Histos[gen->Name][tar->Name][sel.Tag]));
          ///Fill it
          Data::Trees[gen->Name][tar->Name]->Draw(
            (sel.DrawString + " >> " + HistName).c_str(),
            sel.Cut);
          std::cout << "\t\t  Contained: "
            << Data::Histos[gen->Name][tar->Name][sel.Tag]->Integral()
            << std::endl;
        }
        std::cout << "********************************" << std::endl;
      }
      std::cout << "-------------------------------" << std::endl;
    }
    OutputFile->Write();
    return true;
  }

  bool LoadHistogramsFromFile(char const * InputFileName){
    TFile* InputFile = new TFile(InputFileName);
    for(std::vector<PlottingTypes::Generator>::iterator gen = DataSpecifics::Generators.begin();
      gen != DataSpecifics::Generators.end(); ++gen){
      std::cout << "Filling Generator: " << gen->Name << std::endl;

      for(std::vector<PlottingTypes::Target>::iterator tar = gen->Targets.begin();
        tar != gen->Targets.end(); ++tar){
        std::cout << "\tFilling Target: " << tar->Name << std::endl;

        for(std::map<std::string, PlottingTypes::Selection1D>::iterator sel_it = DataSpecifics::Selections1D.begin();
          sel_it != DataSpecifics::Selections1D.end(); ++sel_it){

          PlottingTypes::Selection1D const &sel = sel_it->second;

          if(!(Data::Histos[gen->Name][tar->Name][sel.Tag] =
            dynamic_cast<TH1F*>(InputFile->Get(
              PlottingTypes::MakeHistoName(*gen,*tar,sel).c_str())))){
            std::cerr << "Couldn't find Histo: "
              << PlottingTypes::MakeHistoName(*gen,*tar,sel)
              << " in file: " << InputFileName << std::endl;
            return false;
          }
        }
      }
    }
    return true;
  }

  bool LoadHistograms(bool LoadFromFile=false, char const* InputFileName=""){
    DataSpecifics::SetSpecifics();

    if(LoadFromFile){
      if(!LoadHistogramsFromFile(InputFileName)){
        return false;
      }
    } else {
      if(!LoadFilesAndTrees()){
        return false;
      }
      if(!FillHistograms()){
        return false;
      }
    }
    return true;
  }
}
#endif
