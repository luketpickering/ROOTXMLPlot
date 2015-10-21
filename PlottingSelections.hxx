#ifndef __PLOTTTINGSELECTIONS_HXX_SEEN__
#define __PLOTTTINGSELECTIONS_HXX_SEEN__
#include <map>
#include <string>

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"

namespace PlottingSelections {

  extern std::map<std::string, PlottingTypes::Selection1D> Selections1D;
  extern std::map<std::string, PlottingTypes::Selection2D> Selections2D;

  void InitSelectionsXML(std::string const &confFile);

  PlottingTypes::Selection1D const * FindSel1D(std::string name);
  PlottingTypes::Selection2D const * FindSel2D(std::string name);
}
#endif
