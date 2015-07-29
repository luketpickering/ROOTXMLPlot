#include <algorithm>
#include <set>

#include "TH2F.h"
#include "TROOT.h"
#include "TXMLEngine.h"

#include "PureGenUtils.hxx"

#include "PlottingIO.hxx"
#include "PlottingSelections.hxx"

using namespace PlottingUtils;
using namespace PlottingTypes;
using namespace PlottingSelections;
using namespace DataSpecifics;

namespace {
  TFile* ogf = 0;
  int FileOpsVerb = 0;

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
      if(FileOpsVerb){
        std::cout << "Moving to " << dir << " in file: "
          << file->GetName() << std::endl;
      }
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

  std::string MakeHistoName(Generator const& gen, Sample const& tar,
    Selection const& sel){
    return (gen.Name + "_" + tar.Name + "_" + sel.Name);
  }
}

void SetIOVerbosity(int verb){
  FileOpsVerb = verb;
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
    throw EInvalidGenerator();
    return 0;
  }

  PlottingTypes::Sample const * FindSample(PlottingTypes::Generator const & Gen,
    std::string name){
    if(std::find_if(Gen.Samples.begin(),Gen.Samples.end(),
        [&](Sample const & el){return (el.Name == name);}) !=
        Gen.Samples.end()){
      return &(*std::find_if(Gen.Samples.begin(),Gen.Samples.end(),
          [&](Sample const & el){return (el.Name == name);}));
    }
    std::cout << "[WARN]: Couldn't find Sample named: " << name
      << " for generator " << Gen.Name << std::endl;
    throw EInvalidSample();
    return 0;
  }
}

namespace PlottingIO {

  TFile* OutputTCanviFile = 0;
  int PDFsWritten = 0;
  std::string OutputPDFFileName = "";
  std::string OutputTCanvasFileName = "";
  std::string HistogramCacheFileName = "";

  void InitOutputCanvasFile(){
    if(!OutputTCanviFile && OutputTCanvasFileName.length()){
      OutputTCanviFile = new TFile(OutputTCanvasFileName.c_str(),"RECREATE");
    }
  }

  void CloseOutputCanvasFile(){
    if(PDFsWritten){
      TCanvas dummyCanv("dummyCanv","");
      dummyCanv.Print((std::string(OutputPDFFileName)+"]").c_str());
    }
    OutputTCanviFile->Write();
    OutputTCanviFile->Close();
    delete OutputTCanviFile;
    OutputTCanviFile = 0;
  }

  TCanvas* SaveAndPrint(TCanvas* canv, std::string const &outputname){
    if(canv == nullptr){
      std::cout << "[ERROR]: canv is a nullptr." << std::endl;
      return canv;
    }
    std::cout << "Saving TCanvas: " << canv->GetName() << std::endl;
    if(OutputTCanviFile){
      OutputTCanviFile->WriteTObject(canv);
    } else {
      std::cout << "[WARN]: TCanvas output file does not exist." << std::endl;
    }

    if(!outputname.length()){
      if(!PDFsWritten){
        TCanvas dummyCanv("dummyCanv","");
        dummyCanv.Update();dummyCanv.Print(
          (std::string(OutputPDFFileName)+"[").c_str());
        std::cout << "[INFO]: Written the dummy canvas to " << OutputPDFFileName
          << std::endl;
        PDFsWritten++;
      }
      canv->Print(OutputPDFFileName.c_str());
    } else {
      canv->Print(outputname.c_str());
    }
    return canv;
  }

  PlottingTypes::Sample AddSampleDescriptionXML(
    TXMLEngine &xmlengine, XMLDocPointer_t sampleNode){

    std::string const &sampleName =
      IOUtils::GetAttrValue(xmlengine, sampleNode,"Name");
    std::string const &sampleLocation =
      IOUtils::GetChildNodeContent(xmlengine, sampleNode, "FlatTreeLocation");

    bool valid, Invalidate;
    Invalidate = PGUtils::str2bool(IOUtils::GetAttrValue(xmlengine, sampleNode,
       "InvalidateCache"), valid);
    if(!valid){
      Invalidate = true;
    }

    return PlottingTypes::Sample(sampleName,sampleLocation, Invalidate);
  }

  void AddGeneratorDescriptionXML(TXMLEngine &xmlengine,
    XMLDocPointer_t genDescriptNode){

    std::string genName;

    for(XMLAttrPointer_t nodeAttr_ptr = xmlengine.GetFirstAttr(genDescriptNode);
        (nodeAttr_ptr != NULL);
        nodeAttr_ptr = xmlengine.GetNextAttr(nodeAttr_ptr)){
      if(std::string(xmlengine.GetAttrName(nodeAttr_ptr)) != "Name"){
        std::cout << "[WARN]: Unexpected attribute on "
          << xmlengine.GetNodeName(genDescriptNode) << " node: "
          << xmlengine.GetAttrName(nodeAttr_ptr) << std::endl;
        continue;
      }
      genName = xmlengine.GetAttrValue(nodeAttr_ptr);
      break;
    }

    if(genName.length() == 0){
      std::cout << "[WARN]: Found un-Named generator node." << std::endl;
      return;
    }

    Data::Generators.emplace_back(genName.c_str());
    auto &Samples = Data::Generators.back().Samples;

    for(XMLNodePointer_t sampleNode =
        xmlengine.GetChild(genDescriptNode);
      (sampleNode != NULL);
      sampleNode = xmlengine.GetNext(sampleNode)){
      PlottingTypes::Sample const &sample =
        AddSampleDescriptionXML(xmlengine,sampleNode);

      if(!sample.IsValid()){
        std::cout << "[WARN]: Found invalid sample: " << sample << std::endl;
        continue;
      }

      std::cout << "[INFO]: Adding Sample: " << sample << std::endl;

      Samples.emplace_back(sample);
    }
  }

  bool DefineGeneratorsXML(std::string configXML){
    TXMLEngine* xmlengine = new TXMLEngine;
    xmlengine->SetSkipComments();

    XMLNodePointer_t confSecNode =
      IOUtils::GetNamedChildElementOfDocumentRoot(*xmlengine,
        configXML, "Samples");
    if(!confSecNode){
      delete xmlengine; return false;
    }

    for(XMLNodePointer_t genDescriptNode = xmlengine->GetChild(confSecNode);
      (genDescriptNode != NULL);
      genDescriptNode = xmlengine->GetNext(genDescriptNode)){
      AddGeneratorDescriptionXML(*xmlengine,genDescriptNode);
    }

    delete xmlengine;
    return true;
  }

  ///Loops through defined generator and Sample files and loads the TTrees into
  ///memory
  bool LoadFilesAndTrees(){
    for(auto const &Gen : Data::Generators){
      for(auto const &Sm : Gen.Samples){

        if(!(Data::Files[Gen.Name][Sm.Name] =
            new TFile(Sm.FileLocation.c_str()))){
          std::cerr << "[ERROR]: Failed to load file: " << Sm.FileLocation
            << " for Generator: " << Gen.Name << ", Sample: " << Sm.Name
            << std::endl;
            return false;
        }
        if( !(Data::Trees[Gen.Name][Sm.Name] =
              dynamic_cast<TTree*>(
                Data::Files[Gen.Name][Sm.Name]\
                ->Get(VarTreeName)) ) ){
          std::cerr << "[ERROR]: Failed to load TTree(" << VarTreeName
            << "): for Generator: " << Gen.Name
            << ", Sample: " << Sm.Name << std::endl;
            return false;
        }
      }
    }
    return true;
  }

  std::string AddOffSetToDrawString(std::string draw, double off){
    std::stringstream ss;
    ss << "(" + draw + "+";
    ss << off*1.1 << ")";
    return ss.str();
  }

  std::set<std::string> DrawnTH1s;

  ///Creates a new histogram and fills with the specified selection from the
  ///tree specified by Generator and Sample. 1D Selections
  TH1* FillHistogram(Generator const & Gen, Sample const & Sm,
    Selection1D const & Sel){
    std::string HistName = MakeHistoName(Gen, Sm, Sel);

    //Short circuit to stop invalidate cache requiring multiple redraws
    if(DrawnTH1s.count(HistName)){
      std::cout << "[INFO]: Short circuit returning " << HistName << std::endl;
      return Data::Histos[Gen.Name][Sm.Name][Sel.Name];
    }

    std::cout << "\t\"HistName: \"" << HistName << std::endl;
    std::cout << "\t" << Sel << std::endl;

    std::string PlotDrawString = Sel.DrawString;

    SaveGFile();
    //Stick stuff in sensible places
    MkMv(Data::HistoCacheFile,(Gen.Name + "/" + Sm.Name));

    ///Make the TH1F
    Double_t offset = 0;
    if(Sel.DoPerUseXOffset){ //If we want each Sample to be slightly
      //offset from each other
      size_t TarNum = std::distance(
        Gen.Samples.begin(),
        std::find(Gen.Samples.begin(), Gen.Samples.end(),Sm));
      if(TarNum == Gen.Samples.size()){
        std::cout << "[ERROR] Couldn't find this Sample descriptor: "
          << Gen.Name << std::endl;
      }
      offset = ((Sel.XBinUpper - Sel.XBinLow)/Sel.NBins)*Sel.PerUseXOffset*TarNum;
      PlotDrawString = AddOffSetToDrawString(Sel.DrawString, offset);
    }

    try{
      Data::Histos[Gen.Name][Sm.Name][Sel.Name] =
        new TH1F(HistName.c_str(),"",
          Sel.NBins,Sel.XBinLow+offset,Sel.XBinUpper+offset);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        RevertGFile(); //Return to previous gfile
        return 0;
    }

    ///Fill it
    Data::Trees[Gen.Name][Sm.Name]->Draw(
      (PlotDrawString + " >> " + HistName).c_str(),
      Sel.Cut,"",DrawNoMore);

    std::cout << "\t\t\"" << Data::Histos[Gen.Name][Sm.Name][Sel.Name]->GetName()
      << "\" Contained: [Uf:"
      << Data::Histos[Gen.Name][Sm.Name][Sel.Name]->GetBinContent(0) << "]"
      << Data::Histos[Gen.Name][Sm.Name][Sel.Name]->Integral()
      << "[Of:" << Data::Histos[Gen.Name][Sm.Name][Sel.Name]->GetBinContent(
          Sel.NBins+1)
      << "]"
      << std::endl << std::endl;
    Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
    DrawnTH1s.insert(HistName);
    RevertGFile();
    return Data::Histos[Gen.Name][Sm.Name][Sel.Name];
  }

  ///Creates a new histogram and fills with the specified selection from the
  ///tree specified by Generator and Sample. 2D Selections
  TH2* FillHistogram(Generator const & Gen, Sample const & Sm,
    Selection2D const & Sel){
    std::string HistName = MakeHistoName(Gen, Sm, Sel);
    std::cout << "\t\"HistName: \"" << HistName << std::endl;
    std::cout << "\t" << Sel << std::endl;

    SaveGFile();
    //Stick stuff in sensible places
    MkMv(Data::HistoCacheFile,(Gen.Name + "/" + Sm.Name));

    ///Make the TH1F
    try {
      Data::Histos[Gen.Name][Sm.Name][Sel.Name] =
        new TH2F(HistName.c_str(),"",
          Sel.NXBins,Sel.XBinLow,Sel.XBinUpper,
          Sel.NYBins,Sel.YBinLow,Sel.YBinUpper);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        RevertGFile(); // Return to previous gfile
        return 0;
    }

    ///Fill it
    Data::Trees[Gen.Name][Sm.Name]->Draw(
      (Sel.DrawString + " >> " + HistName).c_str(), Sel.Cut,"",DrawNoMore);
    std::cout << "\t\t\""
      << Data::Histos[Gen.Name][Sm.Name][Sel.Name]->GetName()
      << "\" Contained: "
      << Data::Histos[Gen.Name][Sm.Name][Sel.Name]->Integral()
      << std::endl << std::endl;
    Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
    RevertGFile(); // Return to previous gfile
    return static_cast<TH2*>(Data::Histos[Gen.Name][Sm.Name][Sel.Name]);
  }

  //Wrapper for 1D FillHistogram
  TH1* FillHistogram(std::tuple<PlottingTypes::Generator const &,
    PlottingTypes::Sample const &, PlottingTypes::Selection1D const &>
      HistoDescriptor){
    return FillHistogram(std::get<0>(HistoDescriptor),
      std::get<1>(HistoDescriptor),
      std::get<2>(HistoDescriptor));
  }

  //Wrapper for 2D FillHistogram
  TH2* FillHistogram(std::tuple<PlottingTypes::Generator const &,
    PlottingTypes::Sample const &, PlottingTypes::Selection2D const &>
      HistoDescriptor){
    return FillHistogram(std::get<0>(HistoDescriptor),
      std::get<1>(HistoDescriptor),
      std::get<2>(HistoDescriptor));
  }

  ///Wrapper for 1D FillHistogram which takes generator, Sample and selection
  ///names rather than instances/pointers
  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  TH1* FillHistogram1D(std::string GenName, std::string SmName,
    std::string SelName, bool Redraw){
    PlottingTypes::Generator const * Gen = Data::FindGen(GenName);
    if(!Gen){ return 0; }
    PlottingTypes::Sample const * Sm = Data::FindSample(*Gen,SmName);
    PlottingTypes::Selection1D const * Sel1D = FindSel1D(SelName);
    if(!Sm || !Sel1D){ return 0; }
    Redraw = Redraw || Sm->InvalidateCache || Sel1D->InvalidateCache;
    SaveGFile();
    if(!Redraw && MkMv(Data::HistoCacheFile,(Gen->Name + "/" + Sm->Name))){
      TDirectory *dir =
        Data::HistoCacheFile->GetDirectory(
          (Gen->Name + "/" + Sm->Name).c_str());
      if((Data::Histos[Gen->Name][Sm->Name][Sel1D->Name] =
              dynamic_cast<TH1*>(dir->Get(
                MakeHistoName(*Gen,*Sm,*Sel1D).c_str())))){
        RevertGFile(); // Return to previous gfile
        return Data::Histos[Gen->Name][Sm->Name][Sel1D->Name];
      } else {
        if(FileOpsVerb > 1){
          std::cout << "Couldn't find histo(" << MakeHistoName(*Gen,*Sm,*Sel1D)
            << ") in directory..." << std::endl;
          dir->ls();
        }
      }
    }
    RevertGFile(); // Return to previous gfile
    return FillHistogram(std::make_tuple(*Gen,*Sm,*Sel1D));
  }

  ///Wrapper for 1D FillHistogram which takes generator, Sample and selection
  ///names rather than instances/pointers.
  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  TH2* FillHistogram2D(std::string GenName, std::string SmName,
    std::string SelName, bool Redraw){
    PlottingTypes::Generator const * Gen = Data::FindGen(GenName);
    if(!Gen){ return 0; }
    PlottingTypes::Sample const * Sm = Data::FindSample(*Gen,SmName);
    PlottingTypes::Selection2D const * Sel2D = FindSel2D(SelName);
    if(!Sm || !Sel2D){ return 0; }
    Redraw = Redraw || Sm->InvalidateCache || Sel2D->InvalidateCache;
    SaveGFile();
    if(!Redraw &&  MkMv(Data::HistoCacheFile,(Gen->Name + "/" + Sm->Name))){
      TDirectory *dir =
        Data::HistoCacheFile->GetDirectory(
          (Gen->Name + "/" + Sm->Name).c_str());
      if((Data::Histos[Gen->Name][Sm->Name][Sel2D->Name] =
              dynamic_cast<TH2*>(dir->Get(
                MakeHistoName(*Gen,*Sm,*Sel2D).c_str())))){
        RevertGFile(); // Return to previous gfile
        return static_cast<TH2*>(Data::Histos[Gen->Name][Sm->Name][Sel2D->Name]);
      }
    }
    RevertGFile(); // Return to previous gfile
    return FillHistogram(std::make_tuple(*Gen,*Sm,*Sel2D));
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
      for(auto const &Sm : Gen.Samples){
        std::cout << "\tFilling Sample: " << Sm.Name << std::endl;
        for(auto const &selpair1D : Selections1D){
          Selection1D const &Sel1D = selpair1D.second;

          SaveGFile();
          TDirectory *dir;
          if((dir = Data::HistoCacheFile->GetDirectory(
            (Gen.Name + "/" + Sm.Name).c_str()))){

            if((Data::Histos[Gen.Name][Sm.Name][Sel1D.Name] =
                dynamic_cast<TH1*>(dir->Get(
                  MakeHistoName(Gen,Sm,Sel1D).c_str())))){
              RevertGFile();
              continue; //If we have found it, bail.
            }
          }

          std::cout << "[NOTE]: Couldn't find Histo: "
            << MakeHistoName(Gen,Sm,Sel1D)
            << " in file: " << HistogramCacheFileName << std::endl;
          if(Data::HaveTrees){
            if(!FillHistogram(Gen,Sm,Sel1D)){
              Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
              Data::HistoCacheFile->Close();
              RevertGFile();
              return false; //Probably caught a std::bad_alloc, it's time to leave
            }
            Data::Histos[Gen.Name][Sm.Name][Sel1D.Name]->Write();
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
            (Gen.Name + "/" + Sm.Name).c_str()))){
            if((Data::Histos[Gen.Name][Sm.Name][Sel2D.Name] =
              dynamic_cast<TH2*>(dir->Get(
                MakeHistoName(Gen,Sm,Sel2D).c_str())))){
              RevertGFile();
              continue;
            }
          }
          std::cerr << "Couldn't find Histo: "
            << MakeHistoName(Gen,Sm,Sel2D)
            << " in file: " << HistogramCacheFileName << std::endl;
          if(Data::HaveTrees){
            if(!FillHistogram(Gen,Sm,Sel2D)){
              Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
              Data::HistoCacheFile->Close();
              RevertGFile();
              return false;
            }
            Data::Histos[Gen.Name][Sm.Name][Sel2D.Name]->Write();
          } else {
            std::cout << "We don't seem to have the data files, cannot load"
            " this histogram." << std::endl;
          }
          RevertGFile();
        }
      }
    }
    Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
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
    Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
    return true;
  }

  ///Function which loads the files and draws any pre-requested histograms.
  bool InitialiseHistogramCache(char const* HistogramCacheFileName,
    std::string const & XMLConfFile, std::string const & SelectionsXMLFile,
    bool PreLoadAll){

    DefineGeneratorsXML(XMLConfFile);

    (void)Verbosity;
    // Data::DefineGenerators();
    PlottingSelections::InitSelectionsXML(SelectionsXMLFile);
    Data::HaveTrees = LoadFilesAndTrees();

    // return false;

    if(!Data::HaveTrees){ return false; }

    if(PreLoadAll){
      return LoadHistogramsFromFile(HistogramCacheFileName);
    }

    return LoadSpecificHistograms(HistogramCacheFileName);
  }
}
