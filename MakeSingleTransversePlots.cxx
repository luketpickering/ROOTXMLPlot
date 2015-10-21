#include <utility>
#include <string>
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

#include "PlottingIO.hxx"
#include "PlottingAutoPlots.hxx"

using namespace PlottingIO;

namespace {
  std::string SamplesXML = "";
  std::string SelectionsXML = "";
  std::string PlotsXML = "";

  bool MakeAutoPlots = true;
  bool MakeSpecialPlots = false;
  bool PreCache = false;
}

namespace {

void SetOpts(){

  CLIArgs::AddOpt("-Sa", "--Samples-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading samples from file : " << opt << std::endl;
      SamplesXML = opt;
      return true;
    }, true,[](){},"<Samples.xml>");

  CLIArgs::AddOpt("-Se", "--Selections-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading selections from file : " << opt << std::endl;
      SelectionsXML = opt;
      return true;
    }, false,[](){SelectionsXML = SamplesXML;},
    "<Selections.xml {Default: Value of -Sa}>");

  CLIArgs::AddOpt("-P", "--Plots-file", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading plots from file : " << opt << std::endl;
      PlotsXML = opt;
      MakeAutoPlots = true;
      return true;
    }, false,[](){PlotsXML = SelectionsXML;},
    "<Plots.xml {Default: Value of -Se}>");

  CLIArgs::AddOpt("-o", "--output-pdf", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Writing special plots to : " << opt << std::endl;
      OutputPDFFileName = opt;
      MakeSpecialPlots = true;
      return true;
    }, false,[](){MakeSpecialPlots = false; OutputPDFFileName="output-plots.pdf";},
    "<output-plots.pdf {Default:output-plots.pdf}>");

    CLIArgs::AddOpt("-C", "--Canvas-File", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Writing TCanvases to: " << opt << std::endl;
      OutputTCanvasFileName = opt;
      return true;
    }, false,[](){OutputTCanvasFileName = "TCanvi.root";},
    "<TFileName.root {Default: TCanvi.root}>");

    CLIArgs::AddOpt("-x", "--Cache-File", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Using " << opt << " as histogram Cache." <<  std::endl;
      HistogramCacheFileName = opt;
      return true;
    }, false,[](){HistogramCacheFileName = "CachedHistograms.root";},
    "<TFileName.root {Default: CachedHistograms.root}>");

    CLIArgs::AddOpt("-H", "--PreCache", false,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Pre-caching all the histograms " <<  std::endl;
      PreCache = true;
      return true;
    }, false,[](){PreCache = false;},
    "{Default: false} N.B. This may take an obscenely long time.");

  CLIArgs::AddOpt("-T", "--Tree-Name", true,
    [&] (std::string const &opt) -> bool {
      std::cout << "\t--Reading from TTree Named : " << opt << std::endl;
      PlottingIO::SetTreeName(opt);
      return true;
    }, false,[](){},"<TransversitudenessPureSim>");

  CLIArgs::AddOpt("-n", "--nentries-to-draw", true,
    [&] (std::string const &opt) -> bool {
      long long hldr;
      try{
        hldr = std::stol(opt);
      } catch(...){
        std::cerr << "Couldn't parse: " << opt << " to a long." << std::endl;
        return false;
      }
      std::cout << "\t--Histograms drawn will contain a maximum of: " << opt
        <<  " entries." << std::endl;
      PlottingIO::SetMaxEntriesToDraw(hldr);
      return true;
    }, false,[](){},"<TransversitudenessPureSim>");
  }
}

int main(int argc, char const * argv[]){
  SetOpts();

  CLIArgs::AddArguments(argc,argv);
  if(!CLIArgs::HandleArgs()){
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

  CloseOutputCanvasFile();

}
