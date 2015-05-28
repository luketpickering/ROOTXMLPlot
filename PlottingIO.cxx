#include <algorithm>
#include "TH2F.h"
#include "TROOT.h"

#include "PlottingIO.hxx"
#include "PlottingSelections.hxx"

using namespace PlottingUtils;
using namespace PlottingTypes;
using namespace PlottingSelections;
using namespace DataSpecifics;

namespace {
  TFile* ogf = 0;
  int FileOpsVerb = 2;
  void SaveGFile(){
    ogf = gROOT->GetFile();
  }
  void RevertGFile(){
    if(ogf && (ogf != gROOT->GetFile())){
      if(FileOpsVerb){
        std::cout << "Moving back to file: " << ogf->GetName() << std::endl;
      }
      ogf->cd();
      ogf = 0;
    }
  }
  bool MkMv(TFile* file, std::string dir){
    bool NotMk = false;
    if(file->GetDirectory(dir.c_str())){
      file->cd(dir.c_str());
      return true;
    }
    if(FileOpsVerb>1){
      std::cout << "Couldn't GetDirectory " << dir << " in file: "
        << file->GetName() << std::endl;
    }
    if(file->mkdir(dir.c_str())){
      if(FileOpsVerb){
        std::cout << "In File: " << file->GetName() << " made dir: " << dir
          << std::endl;
      }
    } else {
      if(FileOpsVerb){
        std::cout << "Attempted to make another version of " << dir
          << " in " << file->GetName() << std::endl;
      }
      NotMk = true;
    }
    if(!file->cd(dir.c_str())){
      if(FileOpsVerb){
        std::cout << "Failed to cd into: " << dir
          << std::endl;
      }
    }
    return NotMk;
  }
}

//**********************************Globals*************************************
namespace Data {
  TFile* HistoCacheFile = 0;
  bool HaveTrees = false;
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
                                TH1*
                              >
                    >
          > Histos;

  std::vector<PlottingTypes::Generator> Generators;

  PlottingTypes::Generator const * FindGen(std::string name){
    if(std::find_if(Generators.begin(), Generators.end(),
        [&](Generator const & el){return (el.Name == name);}) !=
          Generators.end()){
      return &(*std::find_if(Generators.begin(), Generators.end(),
        [&](Generator const & el){return (el.Name == name);}));
    }
    std::cout << "[WARN]: Couldn't find Generator named: " << name
      << std::endl;
    return 0;
  }

  PlottingTypes::Target const * FindTar(PlottingTypes::Generator const & Gen,
    std::string name){
    if(std::find_if(Gen.Targets.begin(),Gen.Targets.end(),
        [&](Target const & el){return (el.Name == name);}) !=
        Gen.Targets.end()){
      return &(*std::find_if(Gen.Targets.begin(),Gen.Targets.end(),
          [&](Target const & el){return (el.Name == name);}));
    }
    std::cout << "[WARN]: Couldn't find Target named: " << name
      << std::endl;
    return 0;
  }

///Add Generators and samples here.
void DefineGenerators(){
  Generators.emplace_back("NEUT", NEUT_Line);
  auto &NEUTTargets = Generators.back().Targets;
  //***************************NEUT Targets***********************************
  NEUTTargets.emplace_back("C_4GeV", C_Color,
    "../neut/C_4GeV/neut_C_4GeV_puresim_transversityVars.root");
  NEUTTargets.emplace_back("C", C_Color,
    "../neut/C/neut_C_puresim_transversityVars.root");
  NEUTTargets.emplace_back("C_geniemaqe", C_hN_Color,
    "../neut/C_geniemaqe/neut_C_geniemaqe_puresim_transversityVars.root");
  NEUTTargets.emplace_back("C_nobe", C_Color,
    "../neut/C_nobe/neut_C_nobe_puresim_transversityVars.root");
  NEUTTargets.emplace_back("C_nofm", C_Color,
    "../neut/C_nofm/neut_C_nofm_puresim_transversityVars.root");
  NEUTTargets.emplace_back("C_nonucrescat", C_norescat_Color,
    "../neut/C_nonucrescat/neut_C_nonucrescat_puresim_transversityVars.root");
  NEUTTargets.emplace_back("Fe", Fe_norescat_Color,
    "../neut/Fe/neut_Fe_puresim_transversityVars.root");
  NEUTTargets.emplace_back("O", O_Color,
    "../neut/O/neut_O_puresim_transversityVars.root");
  NEUTTargets.emplace_back("Pb", Pb_Color,
    "../neut/Pb/neut_Pb_puresim_transversityVars.root");
  NEUTTargets.emplace_back("Pb_nonucrescat", Fe_nofm_Color,
    "../neut/Pb_nonucrescat/neut_Pb_nonucrescat_puresim_transversityVars.root");

  Generators.emplace_back("GENIE", GENIE_Line);
  auto &GENIETargets = Generators.back().Targets;
  // ***************************GENIE Targets**********************************
  GENIETargets.emplace_back("C", C_Color,
    "../genie/C/genie_C_puresim_transversityVars.root");
  GENIETargets.emplace_back("C", C_Color,
    "../genie/C/genie_C_puresim_transversityVars.root");
  GENIETargets.emplace_back("C_noht", C_Color,
    "../genie/C_noht/genie_C_noht_puresim_transversityVars.root");
  GENIETargets.emplace_back("Fe", Fe_Color,
    "../genie/Fe/genie_Fe_puresim_transversityVars.root");
  GENIETargets.emplace_back("O", O_Color,
    "../genie/O/genie_O_puresim_transversityVars.root");
  GENIETargets.emplace_back("Pb", Pb_Color,
    "../genie/Pb/genie_Pb_puresim_transversityVars.root");
  GENIETargets.emplace_back("Pb_noht", Pb_Color,
    "../genie/Pb_noht/genie_Pb_noht_puresim_transversityVars.root");

  Generators.emplace_back("NuWro", NuWro_Line);
  auto &NuWroTargets = Generators.back().Targets;
  //***************************NuWro Targets**********************************
  NuWroTargets.emplace_back("C_nobe", C_Color,
    "../nuwro/C_nobe/nuwro_C_nobe_puresim_transversityVars.root");
  NuWroTargets.emplace_back("C_nofm", C_Color,
    "../nuwro/C_nofm/nuwro_C_nofm_puresim_transversityVars.root");
  NuWroTargets.emplace_back("C_nokasada", C_Color,
    "../nuwro/C_nokasada/nuwro_C_nokasada_puresim_transversityVars.root");
  NuWroTargets.emplace_back("C", C_Color,
    "../nuwro/C/nuwro_C_puresim_transversityVars.root");
  NuWroTargets.emplace_back("Fe", Fe_Color,
    "../nuwro/Fe/nuwro_Fe_puresim_transversityVars.root");
  NuWroTargets.emplace_back("O", O_Color,
    "../nuwro/O/nuwro_O_puresim_transversityVars.root");
  NuWroTargets.emplace_back("Pb_nokaskada", Pb_Color,
    "../nuwro/Pb_nokaskada/nuwro_Pb_nokaskada_puresim_transversityVars.root");
  NuWroTargets.emplace_back("Pb", Pb_Color,
    "../nuwro/Pb/nuwro_Pb_puresim_transversityVars.root");

  Generators.emplace_back("GiBUU", NuWro_Line);
  auto &GiBUUTargets = Generators.back().Targets;
  //***************************GIBUU Targets**********************************
  GiBUUTargets.emplace_back("C", C_Color,
    "../gibuu/C/gibuu_C_puresim_transversityVars.root");
  GiBUUTargets.emplace_back("C_mec2", C_Color,
    "../gibuu/C_mecround2/gibuu_C_mecround2_puresim_transversityVars.root");
  GiBUUTargets.emplace_back("Pb", Pb_Color,
    "../gibuu/Pb/gibuu_Pb_puresim_transversityVars.root");
  GiBUUTargets.emplace_back("Pb_moretimesteps", Pb_Color,
    "../gibuu/Pb_moretimesteps/gibuu_Pb_moretimesteps_puresim_transversityVars.root");
  GiBUUTargets.emplace_back("Pb_notimesteps", Pb_Color,
    "../gibuu/Pb_notimesteps/gibuu_Pb_notimesteps_puresim_transversityVars.root");
  GiBUUTargets.emplace_back("Fe", Fe_Color,
    "../gibuu/Fe/gibuu_Fe_puresim_transversityVars.root");
  GiBUUTargets.emplace_back("O", O_Color,
    "../gibuu/O/gibuu_O_puresim_transversityVars.root");
}

}

namespace PlottingIO {

  ///Loops through defined generator and target files and loads the TTrees into
  ///memory
  bool LoadFilesAndTrees(){
    for(auto const &Gen : Data::Generators){
      for(auto const &Tar : Gen.Targets){

        if(!(Data::Files[Gen.Name][Tar.Name] =
            new TFile(Tar.FileLocation.c_str()))){
          std::cerr << "Failed to load file: " << Tar.FileLocation
            << " for Generator: " << Gen.Name << ", Target: " << Tar.Name
            << std::endl;
            return false;
        }
        if( !(Data::Trees[Gen.Name][Tar.Name] =
              dynamic_cast<TTree*>(
                Data::Files[Gen.Name][Tar.Name]\
                ->Get(VarTreeName)) ) ){
          std::cerr << "Failed to load TTree: for Generator: " << Gen.Name
            << ", Target: " << Tar.Name << std::endl;
            return false;
        }

      }
    }
    return true;
  }

  std::string AddOffSetToDrawString(std::string draw, double off){
    std::stringstream ss;
    ss << "(" + draw + "+";
    ss << off << ")";
    return ss.str();
  }

  ///Creates a new histogram and fills with the specified selection from the
  ///tree specified by Generator and Target. 1D Selections
  TH1* FillHistogram(Generator const & Gen, Target const & Tar,
    Selection1D const & Sel){
    std::string HistName = MakeHistoName(Gen, Tar, Sel);
    std::cout << "\t\tSel: " << Sel.Name
      << "\n\t\t  Draw: \"" << Sel.DrawString
      << "\"\n\t\t  HistName: \"" << HistName
      << "\"\n\t\t  Cut: \""
        << ConfineWidthStringWithIndent("\t\t\t",Sel.Cut.GetTitle(),50)
      <<"\"\n\t\t  Bins: [ " << Sel.XBinLow << " .. " << Sel.NBins << " .. "
      << Sel.XBinUpper << " ]" << std::endl;

    std::string PlotDrawString = Sel.DrawString;

    SaveGFile();
    //Stick stuff in sensible places
    MkMv(Data::HistoCacheFile,(Gen.Name + "/" + Tar.Name));

    ///Make the TH1F
    Double_t offset = 0;
    if(Sel.DoPerUseXOffset){ //If we want each Target to be slightly
      //offset from each other
      size_t TarNum = std::distance(
        Gen.Targets.begin(),
        std::find(Gen.Targets.begin(), Gen.Targets.end(),Tar));
      if(TarNum == Gen.Targets.size()){
        std::cout << "[ERROR] Couldn't find this Target descriptor: "
          << Gen.Name << std::endl;
      }
      offset = ((Sel.XBinUpper - Sel.XBinLow)/Sel.NBins)*0.005*TarNum;
      PlotDrawString = AddOffSetToDrawString(Sel.DrawString, offset);
    }

    try{
      Data::Histos[Gen.Name][Tar.Name][Sel.Tag] =
        new TH1F(HistName.c_str(),"",
          Sel.NBins,Sel.XBinLow+offset,Sel.XBinUpper+offset);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        RevertGFile(); //Return to previous gfile
        return 0;
    }

    if(Sel.UseLogX){
      PlottingUtils::BinLogX(Data::Histos[Gen.Name][Tar.Name][Sel.Tag]\
        ->GetXaxis());
    }
    ///Fill it
    Data::Trees[Gen.Name][Tar.Name]->Draw(
      (PlotDrawString + " >> " + HistName).c_str(),
      Sel.Cut,"",DrawNoMore);
    std::cout << "\t\t\"" << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetName()
      << "\" Contained: [Uf:"
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetBinContent(0) << "]"
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->Integral()
      << "[Of:" << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetBinContent(
          Sel.NBins+1)
      << "]"
      << std::endl << std::endl;
    RevertGFile();
    return Data::Histos[Gen.Name][Tar.Name][Sel.Tag];
  }

  ///Creates a new histogram and fills with the specified selection from the
  ///tree specified by Generator and Target. 2D Selections
  TH1* FillHistogram(Generator const & Gen, Target const & Tar,
    Selection2D const & Sel){
    std::string HistName = MakeHistoName(Gen, Tar, Sel);
    std::cout << "\t\tSel: " << Sel.Name
      << "\n\t\t  Draw: \"" << Sel.DrawString
      << "\"\n\t\t  HistName: \"" << HistName
      << "\"\n\t\t  Cut: \""
        << ConfineWidthStringWithIndent("\t\t\t",Sel.Cut.GetTitle(),50)
      <<"\"\n\t\t  XBins: [ " << Sel.XBinLow << " .. " << Sel.NXBins << " .. "
      << Sel.XBinUpper << " ]"
      <<"\"\n\t\t  YBins: [ " << Sel.YBinLow << " .. " << Sel.NYBins << " .. "
      << Sel.YBinUpper << " ]" << std::endl;

    SaveGFile();
    //Stick stuff in sensible places
    MkMv(Data::HistoCacheFile,(Gen.Name + "/" + Tar.Name));

    ///Make the TH1F
    try {
      Data::Histos[Gen.Name][Tar.Name][Sel.Tag] =
        new TH2F(HistName.c_str(),"",
          Sel.NXBins,Sel.XBinLow,Sel.XBinUpper,
          Sel.NYBins,Sel.YBinLow,Sel.YBinUpper);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        RevertGFile(); // Return to previous gfile
        return 0;
    }

    if(Sel.UseLogX){
      PlottingUtils::BinLogX(Data::Histos[Gen.Name][Tar.Name][Sel.Tag]\
        ->GetXaxis());
    }
    if(Sel.UseLogY){
      PlottingUtils::BinLogX(Data::Histos[Gen.Name][Tar.Name][Sel.Tag]\
        ->GetYaxis());
    }
    ///Fill it
    Data::Trees[Gen.Name][Tar.Name]->Draw(
      (Sel.DrawString + " >> " + HistName).c_str(), Sel.Cut,"",DrawNoMore);
    std::cout << "\t\t\""
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetName()
      << "\" Contained: "
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->Integral()
      << std::endl << std::endl;
    RevertGFile(); // Return to previous gfile
    return Data::Histos[Gen.Name][Tar.Name][Sel.Tag];
  }

  //Wrapper for 1D FillHistogram
  TH1* FillHistogram(std::tuple<PlottingTypes::Generator const &,
    PlottingTypes::Target const &, PlottingTypes::Selection1D const &>
      HistoDescriptor){
    return FillHistogram(std::get<0>(HistoDescriptor),
      std::get<1>(HistoDescriptor),
      std::get<2>(HistoDescriptor));
  }

  //Wrapper for 2D FillHistogram
  TH1* FillHistogram(std::tuple<PlottingTypes::Generator const &,
    PlottingTypes::Target const &, PlottingTypes::Selection2D const &>
      HistoDescriptor){
    return FillHistogram(std::get<0>(HistoDescriptor),
      std::get<1>(HistoDescriptor),
      std::get<2>(HistoDescriptor));
  }

  ///Wrapper for 1D FillHistogram which takes generator, target and selection
  ///names rather than instances/pointers
  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  TH1* FillHistogram1D(std::string GenName, std::string TarName,
    std::string SelName){
    PlottingTypes::Generator const * Gen = Data::FindGen(GenName);
    if(!Gen){ return 0; }
    PlottingTypes::Target const * Tar = Data::FindTar(*Gen,TarName);
    PlottingTypes::Selection1D const * Sel1D = FindSel1D(SelName);
    if(!Tar || !Sel1D){ return 0; }
    SaveGFile();
    if(MkMv(Data::HistoCacheFile,(Gen->Name + "/" + Tar->Name))){
      TDirectory *dir =
        Data::HistoCacheFile->GetDirectory(
          (Gen->Name + "/" + Tar->Name).c_str());
      if((Data::Histos[Gen->Name][Tar->Name][Sel1D->Tag] =
              dynamic_cast<TH1*>(dir->Get(
                MakeHistoName(*Gen,*Tar,*Sel1D).c_str())))){
        RevertGFile(); // Return to previous gfile
        return Data::Histos[Gen->Name][Tar->Name][Sel1D->Tag];
      }
    }
    RevertGFile(); // Return to previous gfile
    return FillHistogram(std::make_tuple(*Gen,*Tar,*Sel1D));
  }

  ///Wrapper for 1D FillHistogram which takes generator, target and selection
  ///names rather than instances/pointers.
  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  TH1* FillHistogram2D(std::string GenName, std::string TarName,
    std::string SelName){
    PlottingTypes::Generator const * Gen = Data::FindGen(GenName);
    if(!Gen){ return 0; }
    PlottingTypes::Target const * Tar = Data::FindTar(*Gen,TarName);
    PlottingTypes::Selection2D const * Sel2D = FindSel2D(SelName);
    if(!Tar || !Sel2D){ return 0; }
    SaveGFile();
    if(MkMv(Data::HistoCacheFile,(Gen->Name + "/" + Tar->Name))){
      TDirectory *dir =
        Data::HistoCacheFile->GetDirectory(
          (Gen->Name + "/" + Tar->Name).c_str());
      if((Data::Histos[Gen->Name][Tar->Name][Sel2D->Tag] =
              dynamic_cast<TH2*>(dir->Get(
                MakeHistoName(*Gen,*Tar,*Sel2D).c_str())))){
        RevertGFile(); // Return to previous gfile
        return Data::Histos[Gen->Name][Tar->Name][Sel2D->Tag];
      }
    }
    RevertGFile(); // Return to previous gfile
    return FillHistogram(std::make_tuple(*Gen,*Tar,*Sel2D));
  }

  ///Systematically loads all defined histograms. This may take a while
  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  bool LoadHistogramsFromFile(char const * HistogramCacheFileName){
    Data::HistoCacheFile = new TFile(HistogramCacheFileName, "UPDATE");
    Data::HistoCacheFile->cd();
    if(!Data::HistoCacheFile->IsOpen()){
      std::cout << "Couldnt open the histogram cache file: "
        << HistogramCacheFileName << std::endl;
      return false;
    }
    for(auto const &Gen : Data::Generators){
      std::cout << "Filling Generator: " << Gen.Name << std::endl;
      for(auto const &Tar : Gen.Targets){
        std::cout << "\tFilling Target: " << Tar.Name << std::endl;
        for(auto const &selpair1D : Selections1D){
          Selection1D const &Sel1D = selpair1D.second;

          SaveGFile();
          TDirectory *dir;
          if((dir = Data::HistoCacheFile->GetDirectory(
            (Gen.Name + "/" + Tar.Name).c_str()))){

            if((Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag] =
                dynamic_cast<TH1*>(dir->Get(
                  MakeHistoName(Gen,Tar,Sel1D).c_str())))){
              RevertGFile();
              continue; //If we have found it, bail.
            }
          }

          std::cout << "[NOTE]: Couldn't find Histo: "
            << MakeHistoName(Gen,Tar,Sel1D)
            << " in file: " << HistogramCacheFileName << std::endl;
          if(Data::HaveTrees){
            if(!FillHistogram(Gen,Tar,Sel1D)){
              Data::HistoCacheFile->Write();
              Data::HistoCacheFile->Close();
              RevertGFile();
              return false; //Probably caught a std::bad_alloc, it's time to leave
            }
            Gen.SetSeriesStyle(
              Tar.SetSeriesStyle(
                Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag]));
            Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag]->Write();
          } else {
            std::cout << "[WARN]: We don't seem to have the data files, "
              "cannot load this histogram." << std::endl;
          }
          RevertGFile();
        }

        for(auto const &selpair2D : Selections2D){
          Selection2D const &Sel2D = selpair2D.second;
          SaveGFile();
          TDirectory *dir;
          if((dir = Data::HistoCacheFile->GetDirectory(
            (Gen.Name + "/" + Tar.Name).c_str()))){
            if((Data::Histos[Gen.Name][Tar.Name][Sel2D.Tag] =
              dynamic_cast<TH2*>(dir->Get(
                MakeHistoName(Gen,Tar,Sel2D).c_str())))){
              RevertGFile();
              continue;
            }
          }
          std::cerr << "Couldn't find Histo: "
            << MakeHistoName(Gen,Tar,Sel2D)
            << " in file: " << HistogramCacheFileName << std::endl;
          if(Data::HaveTrees){
            if(!FillHistogram(Gen,Tar,Sel2D)){
              Data::HistoCacheFile->Write();
              Data::HistoCacheFile->Close();
              RevertGFile();
              return false;
            }
            Data::Histos[Gen.Name][Tar.Name][Sel2D.Tag]->Write();
          } else {
            std::cout << "We don't seem to have the data files, cannot load"
            " this histogram." << std::endl;
          }
          RevertGFile();
        }
      }
    }
    Data::HistoCacheFile->Write();
    return true;
  }

  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  bool LoadSpecificHistograms(char const * HistogramCacheFileName){

    Data::HistoCacheFile = new TFile(HistogramCacheFileName, "UPDATE");
    Data::HistoCacheFile->cd();
    if(!Data::HistoCacheFile->IsOpen()){
      std::cout << "Couldnt open the histogram cache file: "
        << HistogramCacheFileName << std::endl;
      return false;
    } else {
      std::cout << "Opened " << HistogramCacheFileName
        << " to cache histograms in." << std::endl;
    }
  //ECalLike
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_CCQE"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_CCRes"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_CCMEC"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECal2ProngHMTrackableMom"));

    MakePDFClone(FillHistogram1D("NuWro","Pb","CCQE"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_CCQE"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_CCRes"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_CCMEC"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECal2ProngHMTrackableMom"));

  //ECalLike-Smeared
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_CCQE_sm"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_CCRes_sm"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_CCMEC_sm"));
    MakePDFClone(FillHistogram1D("NuWro","C","ECalTop2Prong_sm"));

    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_CCQE_sm"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_CCRes_sm"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_CCMEC_sm"));
    MakePDFClone(FillHistogram1D("NuWro","Pb","ECalTop2Prong_sm"));
  //STVars 1D FSIOff
    MakePDFClone(FillHistogram1D("NuWro","C_nokasada","CCQE"));
    MakePDFClone(FillHistogram1D("NuWro","C_nokasada","CCQEDPt"));
    MakePDFClone(FillHistogram1D("NuWro","C_nokasada","CCQEAlphaT"));
    MakePDFClone(FillHistogram1D("NuWro","C_nokasada","ECal2ProngHMTrackableMom"));

  //STVars 2D FSIOff
    MakePDFClone(FillHistogram2D("NuWro","C_nokasada","DeltaPhiTHMPvsMuonPt"));
    MakePDFClone(FillHistogram2D("NuWro","C_nokasada","DeltaAlphaTHMPvsMuonPt"));
    MakePDFClone(FillHistogram2D("NuWro","C_nokasada","DeltaPTHMPvsMuonPt"));

  //STVars 1D FSIOn
    MakePDFClone(FillHistogram1D("NuWro","C","CCQE"));
    MakePDFClone(FillHistogram1D("NuWro","C","CCQEDPt"));
    MakePDFClone(FillHistogram1D("NuWro","C","CCQEAlphaT"));
  //STVars 2D FSIOn
    MakePDFClone(FillHistogram2D("NuWro","C","DeltaPhiTHMPvsMuonPt"));
    MakePDFClone(FillHistogram2D("NuWro","C","DeltaAlphaTHMPvsMuonPt"));
    MakePDFClone(FillHistogram2D("NuWro","C","DeltaPTHMPvsMuonPt"));
    Data::HistoCacheFile->Write();
    return true;
  }

  ///Function which loads the files and draws any pre-requested histograms.
  bool InitialiseHistogramCache(char const* HistogramCacheFileName,
    bool PreLoadAll){

    (void)Verbosity;
    Data::DefineGenerators();
    PlottingSelections::InitSelections();
    Data::HaveTrees = LoadFilesAndTrees();
    if(PreLoadAll){
      return LoadHistogramsFromFile(HistogramCacheFileName);
    }
    return LoadSpecificHistograms(HistogramCacheFileName);
  }
}
