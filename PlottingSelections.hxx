#ifndef __PLOTTTINGSELECTIONS_HXX_SEEN__
#define __PLOTTTINGSELECTIONS_HXX_SEEN__
#include <map>
#include <string>

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"

namespace PlottingSelections {

  extern std::map<std::string, PlottingTypes::Selection1D> Selections1D;
  extern std::map<std::string, PlottingTypes::Selection2D> Selections2D;

  constexpr bool DoLogx = PlottingDefaults::DoLogx;

  constexpr int ProtonMomBins = 100;
  constexpr double ProtonMomMax = 1000;

  constexpr int DPhiTBins = 100;
  constexpr double DPhiTMin = ((!DoLogx)?0.0:1E-2);
  constexpr double DPhiTMax = 3.15;

  constexpr int DPhiTZoomxBins = 50;
  constexpr double DPhiTZoomxMin = ((!DoLogx)?0.0:1E-5);
  constexpr double DPhiTZoomxMax = 1E-1;

  constexpr int DPhiTBins_deg = 500;
  constexpr double DPhiTMin_deg = ((!DoLogx)?0.0:1E-2);
  constexpr double DPhiTMax_deg = 180.0;

  constexpr int DPhiTZoomxBins_deg = 50;
  constexpr double DPhiTZoomxMin_deg = ((!DoLogx)?0.0:1E-5);
  constexpr double DPhiTZoomxMax_deg = 6.0;

  constexpr int DAlphaTBins_deg = 500;
  constexpr double DAlphaTMin_deg = ((!DoLogx)?0.0:1E-2);
  constexpr double DAlphaTMax_deg = 180.0;

  constexpr int DPTBins = 500;
  constexpr double DPTMin = ((!DoLogx)?0.0:1E-2);
  constexpr double DPTMax = 700; //MeV

  void InitSelections();
  void InitSelectionsXML(std::string const &confFile);

  PlottingTypes::Selection1D const * FindSel1D(std::string name);
  PlottingTypes::Selection2D const * FindSel2D(std::string name);
}
#endif
