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

  constexpr bool DoLogx = PlottingDefaults::DoLogx;

  constexpr int ProtonMomBins = 100;
  constexpr double ProtonMomMax = 1000;

  constexpr int DPhiTBins = 100;
  constexpr double DPhiTMin = ((!DoLogx)?0.0:1E-2);
  constexpr double DPhiTMax = 3.15;

  constexpr int DPhiTZoomxBins = 50;
  constexpr double DPhiTZoomxMin = ((!DoLogx)?0.0:1E-5);
  constexpr double DPhiTZoomxMax = 1E-1;

  constexpr int DPhiTBins_deg = 100;
  constexpr double DPhiTMin_deg = ((!DoLogx)?0.0:1E-2);
  constexpr double DPhiTMax_deg = 180.0;

  constexpr int DPhiTZoomxBins_deg = 50;
  constexpr double DPhiTZoomxMin_deg = ((!DoLogx)?0.0:1E-5);
  constexpr double DPhiTZoomxMax_deg = 6.0;

  constexpr int DAlphaTBins_deg = 100;
  constexpr double DAlphaTMin_deg = ((!DoLogx)?0.0:1E-2);
  constexpr double DAlphaTMax_deg = 180.0;

  constexpr int DPTBins = 100;
  constexpr double DPTMin = ((!DoLogx)?0.0:1E-2);
  constexpr double DPTMax = 700; //MeV

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

  extern std::vector<PlottingTypes::Generator> Generators;
  extern std::map<std::string, PlottingTypes::Selection1D> Selections1D;

  extern std::string MuonValidSelect;
  extern std::string FirstProtonValidSelect;
  extern std::string FirstProtonNonZeroSelect;
  extern std::string CCQERCSelect;
  extern std::string CCQEDPhiTSelect;
  extern std::string CCQENoCCQENucEmitSelect;
  extern std::string CCQENucEmitSelect;
  extern std::string CCQESelect;
  extern std::string MECSelect;
  extern std::string CCQELikeSelect;
  extern std::string CC_not_QELikeSelect;
  extern std::string CCQEProtonMomentumSelect;

  extern std::string FirstPiPlusValidSelect;
  extern std::string FirstPiPlusNonZeroSelect;

  extern std::string ResPPiPlusProdSelect;
  extern std::string ResNPiPlusProdSelect;
  extern std::string CohPiPlusProdSelect;
  extern std::string OtherPiPlusSelect;
  extern std::string CCSelect;

  void InitTCuts();
  ///You edit what things are loaded, and thus plottable, here.
  void SetSpecifics();
}

namespace PlottingIO {

  bool LoadFilesAndTrees();
  bool FillHistogram(PlottingTypes::Generator const & gen,
    PlottingTypes::Target const & tar,
    PlottingTypes::Selection1D const & sel);
  bool FillHistogram(PlottingTypes::Generator const & Gen,
    PlottingTypes::Target const & Tar,
    PlottingTypes::Selection2D const & Sel);
  bool LoadHistogramsFromFile(char const * InputFileName);
  bool LoadHistograms(char const* InputFileName="");
}
#endif
