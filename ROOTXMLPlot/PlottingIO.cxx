#include <algorithm>
#include <set>

#include "TH2F.h"
#include "TROOT.h"

#include "PureGenUtils.hxx"

#include "PlottingIO.hxx"
#include "PlottingSelections.hxx"

using namespace PlottingUtils;
using namespace PlottingTypes;
using namespace PlottingSelections;
using namespace DataSpecifics;

namespace DataSpecifics {
  std::string VarTreeName = "TransversitudenessPureSim";
  long long MaxEntriesToDraw = 1000000000; // ROOT6 defines TTree::kMaxEntries;
}

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

  std::string MakeHistoName(SampleGroup const& gen, Sample const& tar,
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

  std::vector<PlottingTypes::SampleGroup> SampleGroups;

  PlottingTypes::SampleGroup const * FindGen(std::string name){
    if(std::find_if(SampleGroups.begin(), SampleGroups.end(),
        [&](SampleGroup const & el){return (el.Name == name);}) !=
          SampleGroups.end()){
      return &(*std::find_if(SampleGroups.begin(), SampleGroups.end(),
        [&](SampleGroup const & el){return (el.Name == name);}));
    }
    std::cout << "[WARN]: Couldn't find SampleGroup named: " << name
      << std::endl;
    throw EInvalidSampleGroup();
    return nullptr;
  }

  PlottingTypes::Sample const * FindSample(
    PlottingTypes::SampleGroup const & Grp,
    std::string name){

    if(std::find_if(Grp.Samples.begin(),Grp.Samples.end(),
        [&](Sample const & el){return (el.Name == name);}) !=
        Grp.Samples.end()){
      return &(*std::find_if(Grp.Samples.begin(),Grp.Samples.end(),
          [&](Sample const & el){return (el.Name == name);}));
    }
    std::cout << "[WARN]: Couldn't find Sample named: " << name
      << " for generator " << Grp.Name << std::endl;
    throw EInvalidSample();
    return nullptr;
  }
}

namespace ExternalData {
  std::map<std::string, TFile*> Files;
  std::vector<PlottingTypes::ExternalFile> ExternalDescriptors;
}

namespace PlottingIO {

  TFile* OutputTCanviFile = nullptr;
  int PDFsWritten = 0;
  std::string OutputPDFFileName = "";
  std::string OutputTCanvasFileName = "";
  std::string HistogramCacheFileName = "";

  void SetTreeName(std::string newName){
    DataSpecifics::VarTreeName = newName;
  }

  void SetMaxEntriesToDraw(long long nentries){
    DataSpecifics::MaxEntriesToDraw = nentries;
  }

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
    OutputTCanviFile = nullptr;
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

  bool AddExternalDataXML(std::string configXML){

    TXMLEngine& xmlengine = XMLDoc::Get("Samples");

    XMLNodePointer_t confSecNode =
      IOUtils::GetNamedChildElementOfDocumentRoot(xmlengine, "ExternalFiles");
    if(!confSecNode){
      return false;
    }

    for(XMLAttrPointer_t fileNode = xmlengine.GetChild(confSecNode);
        (fileNode != nullptr);
        fileNode = xmlengine.GetNext(fileNode)){
      if(std::string(xmlengine.GetNodeName(fileNode)) != "TFile"){
        std::cout << "[WARN]" << "Found an unexpected node named \""
          << xmlengine.GetNodeName(fileNode) << "\" in the ExternalFiles "
          "element. Skipping" << std::endl;
        continue;
      }

      std::string const &externalName =
        IOUtils::GetAttrValue(xmlengine, fileNode,"Name");
      std::string const &externalLocation =
        IOUtils::GetAttrValue(xmlengine, fileNode, "Location");

      if(!externalLocation.length() || ! externalName.length()){
        std::cout << "[WARN] Found dud External File descriptor { "
        << externalName
        << ", " << externalLocation << " }." << std::endl;
        continue;
      }
      std::cout << "[INFO] Found External File Descriptor: { " << externalName
        << ", " << externalLocation << " }." << std::endl;

      ExternalData::ExternalDescriptors.emplace_back(externalName,
        externalLocation);
    }
    return bool(ExternalData::ExternalDescriptors.size());
  }

  void LoadExternalDataTFiles(){
    for(auto const & ef : ExternalData::ExternalDescriptors){
      TFile* extFile = TFile::Open(ef.Location.c_str());
      if(!extFile || !extFile->IsOpen()){
        std::cerr << "[ERROR] Couldn't open TFile: " << ef.Location << std::endl;
        throw EInvalidTFile();
      }
      ExternalData::Files[ef.Name] = extFile;
    }
  }

  TH1* GetExternalTH1(std::string Name, std::string TH1Name){
    if(!ExternalData::Files.count(Name)){
      std::cerr << "[ERROR] Don't know about External TFile: " << Name
        << std::endl;
      throw EInvalidTFile();
    }
    TH1* rtn = dynamic_cast<TH1*>(ExternalData::Files[Name]->Get(TH1Name.c_str()));
    if(!rtn){
      std::cerr << "[ERROR] External TFile: " << Name
        << " didn't contain TH1 " << TH1Name << std::endl;
      throw EInvalidPlot();
    }
    return rtn;
  }

  PlottingTypes::Sample AddSampleDescriptionXML(XMLDocPointer_t sampleNode){

    TXMLEngine& xmlengine = XMLDoc::Get("Samples");

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

  void AddSampleGroupDescriptionXML(XMLDocPointer_t smplGrpDescriptNode){

    TXMLEngine& xmlengine = XMLDoc::Get("Samples");

    std::string SampleGroupName = XMLUtils::GetAttrValue(xmlengine,
      smplGrpDescriptNode, "Name");

    if(SampleGroupName.length() == 0){
      std::cout << "[WARN]: Found un-Named SampleGroup node." << std::endl;
      return;
    }

    Data::SampleGroups.emplace_back(SampleGroupName.c_str());
    auto &Samples = Data::SampleGroups.back().Samples;

    for(XMLNodePointer_t sampleNode =
        xmlengine.GetChild(smplGrpDescriptNode);
      (sampleNode != nullptr);
      sampleNode = xmlengine.GetNext(sampleNode)){
      PlottingTypes::Sample const &sample =
        AddSampleDescriptionXML(sampleNode);

      if(!sample.IsValid()){
        std::cout << "[WARN]: Found invalid sample: " << sample << std::endl;
        continue;
      }
      std::cout << "[INFO]: Adding Sample: " << sample << std::endl;
      Samples.emplace_back(sample);
    }
  }

  bool DefineSampleGroupsXML(){

    TXMLEngine& xmlengine = XMLDoc::Get("Samples");

    XMLNodePointer_t confSecNode =
      IOUtils::GetNamedChildElementOfDocumentRoot(xmlengine, "Samples");
    if(!confSecNode){
      return false;
    }

    for(XMLNodePointer_t smplGrpDescriptNode = xmlengine.GetChild(confSecNode);
      (smplGrpDescriptNode != nullptr);
      smplGrpDescriptNode = xmlengine.GetNext(smplGrpDescriptNode)){
      if(std::string(xmlengine.GetNodeName(smplGrpDescriptNode)) ==
          "SampleGroup"){
        AddSampleGroupDescriptionXML(xmlengine,smplGrpDescriptNode);
      }
    }
    return true;
  }

  ///Loops through defined generator and Sample files and loads the TTrees into
  ///memory
  bool LoadFilesAndTrees(){
    for(auto const &Grp : Data::SampleGroups){
      for(auto const &Sm : Grp.Samples){

        if(!(Data::Files[Grp.Name][Sm.Name] =
            new TFile(Sm.FileLocation.c_str()))){
          std::cerr << "[ERROR]: Failed to load file: " << Sm.FileLocation
            << " for SampleGroup: " << Grp.Name << ", Sample: " << Sm.Name
            << std::endl;
            return false;
        }
        if( !(Data::Trees[Grp.Name][Sm.Name] =
              dynamic_cast<TTree*>(
                Data::Files[Grp.Name][Sm.Name]\
                ->Get(VarTreeName.c_str())) ) ){
          std::cerr << "[ERROR]: Failed to load TTree(" << VarTreeName
            << "): for SampleGroup: " << Grp.Name
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

  ///Creates a new histogram and fills with the specified selection from the
  ///tree specified by SampleGroup and Sample. 1D Selections
  TH1* FillHistogram(SampleGroup const & Grp, Sample const & Sm,
    Selection1D const & Sel){
    std::string HistName = MakeHistoName(Grp, Sm, Sel);

    //Short circuit to stop invalidate cache requiring multiple redraws
    static std::set<std::string> DrawnTH1s;
    if(DrawnTH1s.count(HistName)){
      std::cout << "[INFO]: Short circuit returning " << HistName << std::endl;
      return Data::Histos[Grp.Name][Sm.Name][Sel.Name];
    }

    std::cout << "\t\"HistName: \"" << HistName << std::endl;
    std::cout << "\t" << Sel << std::endl;

    std::string PlotDrawString = Sel.DrawString;

    SaveGFile();
    //Stick stuff in sensible places
    MkMv(Data::HistoCacheFile,(Grp.Name + "/" + Sm.Name));

    ///Make the TH1F
    Double_t offset = 0;
    if(Sel.DoPerUseXOffset){ //If we want each Sample to be slightly
      //offset from each other
      size_t TarNum = std::distance(
        Grp.Samples.begin(),
        std::find(Grp.Samples.begin(), Grp.Samples.end(),Sm));
      if(TarNum == Grp.Samples.size()){
        std::cout << "[ERROR] Couldn't find this Sample descriptor: "
          << Grp.Name << std::endl;
      }
      offset = ((Sel.XBinUpper - Sel.XBinLow)/Sel.NBins)*Sel.PerUseXOffset*TarNum;
      PlotDrawString = AddOffSetToDrawString(Sel.DrawString, offset);
    }

    try{
      Data::Histos[Grp.Name][Sm.Name][Sel.Name] =
        new TH1F(HistName.c_str(),"",
          Sel.NBins,Sel.XBinLow+offset,Sel.XBinUpper+offset);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        RevertGFile(); //Return to previous gfile
        return nullptr;
    }

    ///Fill it
    Data::Trees[Grp.Name][Sm.Name]->Draw(
      (PlotDrawString + " >> " + HistName).c_str(),
      Sel.Cut,"",DataSpecifics::MaxEntriesToDraw);

    std::cout << "\t\t\"" << Data::Histos[Grp.Name][Sm.Name][Sel.Name]->GetName()
      << "\" Contained: [Uf:"
      << Data::Histos[Grp.Name][Sm.Name][Sel.Name]->GetBinContent(0) << "]"
      << Data::Histos[Grp.Name][Sm.Name][Sel.Name]->Integral()
      << "[Of:" << Data::Histos[Grp.Name][Sm.Name][Sel.Name]->GetBinContent(
          Sel.NBins+1)
      << "]"
      << std::endl << std::endl;
    Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
    DrawnTH1s.insert(HistName);
    RevertGFile();
    return Data::Histos[Grp.Name][Sm.Name][Sel.Name];
  }

  ///Creates a new histogram and fills with the specified selection from the
  ///tree specified by SampleGroup and Sample. 2D Selections
  TH2* FillHistogram(SampleGroup const & Grp, Sample const & Sm,
    Selection2D const & Sel){
    std::string HistName = MakeHistoName(Grp, Sm, Sel);

    //Short circuit to stop invalidate cache requiring multiple redraws
    static std::set<std::string> DrawnTH2s;
    if(DrawnTH2s.count(HistName)){
      std::cout << "[INFO]: Short circuit returning " << HistName << std::endl;
      return static_cast<TH2*>(Data::Histos[Grp.Name][Sm.Name][Sel.Name]);
    }

    std::cout << "\t\"HistName: \"" << HistName << std::endl;
    std::cout << "\t" << Sel << std::endl;

    SaveGFile();
    //Stick stuff in sensible places
    MkMv(Data::HistoCacheFile,(Grp.Name + "/" + Sm.Name));

    ///Make the TH1F
    try {
      Data::Histos[Grp.Name][Sm.Name][Sel.Name] =
        new TH2F(HistName.c_str(),"",
          Sel.NXBins,Sel.XBinLow,Sel.XBinUpper,
          Sel.NYBins,Sel.YBinLow,Sel.YBinUpper);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        RevertGFile(); // Return to previous gfile
        return nullptr;
    }

    ///Fill it
    Data::Trees[Grp.Name][Sm.Name]->Draw(
      (Sel.DrawString + " >> " + HistName).c_str(), Sel.Cut,"",
      DataSpecifics::MaxEntriesToDraw);
    std::cout << "\t\t\""
      << Data::Histos[Grp.Name][Sm.Name][Sel.Name]->GetName()
      << "\" Contained: "
      << Data::Histos[Grp.Name][Sm.Name][Sel.Name]->Integral()
      << std::endl << std::endl;
    Data::HistoCacheFile->Write(0,TObject::kWriteDelete);
    DrawnTH2s.insert(HistName);
    RevertGFile(); // Return to previous gfile
    return static_cast<TH2*>(Data::Histos[Grp.Name][Sm.Name][Sel.Name]);
  }

  //Wrapper for 1D FillHistogram
  TH1* FillHistogram(std::tuple<PlottingTypes::SampleGroup const &,
    PlottingTypes::Sample const &, PlottingTypes::Selection1D const &>
      HistoDescriptor){
    return FillHistogram(std::get<0>(HistoDescriptor),
      std::get<1>(HistoDescriptor),
      std::get<2>(HistoDescriptor));
  }

  //Wrapper for 2D FillHistogram
  TH2* FillHistogram(std::tuple<PlottingTypes::SampleGroup const &,
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
  TH1* FillHistogram1D(std::string SampleGroupName, std::string SmName,
    std::string SelName, bool Redraw){
    PlottingTypes::SampleGroup const * Grp = Data::FindGen(SampleGroupName);
    if(!Grp){ return nullptr; }
    PlottingTypes::Sample const * Sm = Data::FindSample(*Grp,SmName);
    PlottingTypes::Selection1D const * Sel1D = FindSel1D(SelName);
    if(!Sm || !Sel1D){ return nullptr; }
    Redraw = Redraw || Sm->InvalidateCache || Sel1D->InvalidateCache;
    SaveGFile();
    if(!Redraw && MkMv(Data::HistoCacheFile,(Grp->Name + "/" + Sm->Name))){
      TDirectory *dir =
        Data::HistoCacheFile->GetDirectory(
          (Grp->Name + "/" + Sm->Name).c_str());
      if((Data::Histos[Grp->Name][Sm->Name][Sel1D->Name] =
              dynamic_cast<TH1*>(dir->Get(
                MakeHistoName(*Grp,*Sm,*Sel1D).c_str())))){
        RevertGFile(); // Return to previous gfile
        return Data::Histos[Grp->Name][Sm->Name][Sel1D->Name];
      } else {
        if(FileOpsVerb > 1){
          std::cout << "Couldn't find histo(" << MakeHistoName(*Grp,*Sm,*Sel1D)
            << ") in directory..." << std::endl;
          dir->ls();
        }
      }
    }
    RevertGFile(); // Return to previous gfile
    return FillHistogram(std::make_tuple(*Grp,*Sm,*Sel1D));
  }

  ///Wrapper for 1D FillHistogram which takes generator, Sample and selection
  ///names rather than instances/pointers.
  ///If a histogram exists in the Cache file it will be loaded from there
  ///rather than redrawn.
  TH2* FillHistogram2D(std::string SampleGroupName, std::string SmName,
    std::string SelName, bool Redraw){
    PlottingTypes::SampleGroup const * Grp = Data::FindGen(SampleGroupName);
    if(!Grp){ return nullptr; }
    PlottingTypes::Sample const * Sm = Data::FindSample(*Grp,SmName);
    PlottingTypes::Selection2D const * Sel2D = FindSel2D(SelName);
    if(!Sm || !Sel2D){ return nullptr; }
    Redraw = Redraw || Sm->InvalidateCache || Sel2D->InvalidateCache;
    SaveGFile();
    if(!Redraw &&  MkMv(Data::HistoCacheFile,(Grp->Name + "/" + Sm->Name))){
      TDirectory *dir =
        Data::HistoCacheFile->GetDirectory(
          (Grp->Name + "/" + Sm->Name).c_str());
      if((Data::Histos[Grp->Name][Sm->Name][Sel2D->Name] =
              dynamic_cast<TH2*>(dir->Get(
                MakeHistoName(*Grp,*Sm,*Sel2D).c_str())))){
        RevertGFile(); // Return to previous gfile
        return static_cast<TH2*>(Data::Histos[Grp->Name][Sm->Name][Sel2D->Name]);
      }
    }
    RevertGFile(); // Return to previous gfile
    return FillHistogram(std::make_tuple(*Grp,*Sm,*Sel2D));
  }

  ///Function which loads the files and draws any pre-requested histograms.
  bool InitialiseHistogramCache(char const* HistogramCacheFileName){

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

  bool InitialisedInputData(std::string const & XMLConfFile){

    (void)Verbosity;
    XMLDoc::LoadDoc("Samples", XMLConfFile);

    try {
      DefineSampleGroupsXML();
      if(AddExternalDataXML()){
        LoadExternalDataTFiles();
      }
    } catch (...){
      XMLDoc::Destroy("Samples");
      return false;
    }

    XMLDoc::Destroy("Samples");

    return LoadFilesAndTrees();
  }
}
