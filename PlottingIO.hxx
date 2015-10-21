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
  extern std::map< std::string,
            std::map< std::string,
                      TFile*>
          > Files;
  extern std::map< std::string,
            std::map< std::string,
                      TTree*>
          > Trees;
  extern std::map< std::string,
            std::map< std::string,
                      std::map< std::string,
                                TH1*
                              >
                    >
          > Histos;
}

namespace ExternalData {
  extern std::map<std::string, TFile*> Files;
}

namespace DataSpecifics {
  extern std::string VarTreeName;
}

namespace Data {
  PlottingTypes::Generator const * FindGen(std::string name);
  PlottingTypes::Sample const * FindTar(PlottingTypes::Generator const & Gen,
    std::string name);
}

namespace PlottingIO {

  extern std::string OutputPDFFileName;
  extern std::string OutputTCanvasFileName;
  extern std::string HistogramCacheFileName;

  void SetTreeName(std::string newName);
  void SetMaxEntriesToDraw(long long nentries);

  void InitOutputCanvasFile();
  void CloseOutputCanvasFile();
  TCanvas* SaveAndPrint(TCanvas* canv, std::string const &outputname="");

  void SetIOVerbosity(int verb);

  TH1* GetExternalTH1(std::string GenName, std::string TarName);
  TH1* FillHistogram1D(std::string GenName, std::string TarName,
    std::string SelName, bool Redraw=false);
  TH2* FillHistogram2D(std::string GenName, std::string TarName,
    std::string SelName, bool Redraw=false);

  bool InitialiseHistogramCache(char const* HistogramCacheFileName,
    std::string const & XMLConfFile, std::string const & SelectionsXMLFile,
    bool PreLoadAll=false);
}
#endif
