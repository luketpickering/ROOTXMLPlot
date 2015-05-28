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

namespace DataSpecifics {

  constexpr int C_Color = kRed+2;
  constexpr int C_hN_Color = kAzure-7;
  constexpr int C_norescat_Color = kPink+9;
  constexpr int C_2_9_Color = kBlack;
  constexpr int O_Color = kBlue-3;
  constexpr int Pb_Color = kOrange;
  constexpr int Fe_Color = kMagenta+2;
  constexpr int Fe_norescat_Color = kCyan;
  constexpr int Fe_nofm_Color = kGreen+2;

  constexpr int NEUT_Line = 1;
  constexpr int GENIE_Line = 2;
  constexpr int NuWro_Line = 3;

  constexpr char const * VarTreeName = "TransversitudenessPureSim";

  constexpr Long64_t DrawNoMore = 1000000;

}

namespace Data {
  PlottingTypes::Generator const * FindGen(std::string name);
  PlottingTypes::Target const * FindTar(PlottingTypes::Generator const & Gen,
    std::string name);
}

namespace PlottingIO {

  TH1* FillHistogram1D(std::string GenName, std::string TarName,
    std::string SelName);
  TH1* FillHistogram2D(std::string GenName, std::string TarName,
    std::string SelName);

  bool InitialiseHistogramCache(char const* HistogramCacheFileName,
    bool PreLoadAll=false);
}
#endif
