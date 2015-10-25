#ifndef PLOTTINGAUTOPLOTS_HXX_SEEN
#define PLOTTINGAUTOPLOTS_HXX_SEEN

#include <iostream>
#include <string>

#include "TXMLEngine.h"

#include "PlottingTypes.hxx"
#include "PlottingUtils.hxx"

namespace AutoPlots {

extern std::vector<PlottingTypes::PlotDescriptor> AutoPlots;

bool MakeAutoPlotsXML(std::string const &fileName);

}

#endif
