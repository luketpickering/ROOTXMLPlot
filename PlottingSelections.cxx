#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "TH1.h"
#include "TH2.h"

#include "PureGenUtils.hxx"

#include "PlottingSelections.hxx"

using namespace PlottingUtils;
using namespace PlottingTypes;
using namespace PlottingDefaults;

namespace PlottingSelections {

  std::map<std::string, PlottingTypes::Selection1D> Selections1D;
  std::map<std::string, PlottingTypes::Selection2D> Selections2D;
  std::map<std::string,std::string> PreSelections;

// //**********************************Globals*************************************

  PlottingTypes::Selection1D const * FindSel1D(std::string name){
    if(Selections1D.count(name)){
      return &Selections1D[name];
    }
    std::cout << "[WARN]: Couldn't find Selection1D named: " << name
      << std::endl;
    return 0;
  }

  PlottingTypes::Selection2D const * FindSel2D(std::string name){
    if(Selections2D.count(name)){
      return &Selections2D[name];
    }
    std::cout << "[WARN]: Couldn't find Selection2D named: " << name
      << std::endl;
    return 0;
  }


void AddPreSelectionDescriptionXML(TXMLEngine &xmlengine,
  XMLDocPointer_t PreSelsNode){

  for(XMLNodePointer_t SelsDescriptNode = xmlengine.GetChild(PreSelsNode);
    (SelsDescriptNode != NULL);
    SelsDescriptNode = xmlengine.GetNext(SelsDescriptNode)){

    if(std::string(xmlengine.GetNodeName(SelsDescriptNode)) != "PreSelection"){
      std::cout << "[WARN]: Found unexpected node: "
        << xmlengine.GetNodeName(SelsDescriptNode)
        << " in the PreSelection element." << std::endl;
      continue;
    }
    std::string const &selName =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"Name");

    std::string const &selectionString =
      IOUtils::GetNodeContent(xmlengine, SelsDescriptNode);

    std::cout << "[INFO]: Found Preselection { Name: " << selName
    << ", Selection: " << selectionString << "}" << std::endl;

    if(!selName.length() || !selectionString.length()){
      return;
    }

    PreSelections[selName] = selectionString;
  }
}

void Add1DSelection(TXMLEngine &xmlengine,
  XMLDocPointer_t SelsDescriptNode){

    std::string const &selName =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"Name");

    if(selName.length() == 0){
      std::cout << "[ERROR]: Selection had no Name." << std::endl;
      return;
    }

    std::string const &selPreSelection =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"PreSelection");

    std::string const &selNBins =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"NBins");
    std::string const &selXBinLow =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinLow");
    std::string const &selXBinUpper =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinUpper");

    int NBins;
    double XBinLow;
    double XBinUpper;

    try{
      NBins = std::stoi(selNBins);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << selNBins << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinLow = std::stod(selXBinLow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << selXBinLow << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinUpper = std::stod(selXBinUpper);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << selXBinUpper << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    bool valid, Invalidate;
    Invalidate = PGUtils::str2bool(
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,
       "InvalidateCache"), valid);
    if(!valid){
      Invalidate = true;
    }

    std::string selectionString = "";
    if(selPreSelection.length() && PreSelections.count(selPreSelection)){
      selectionString = CombineWeightStrings(PreSelections[selPreSelection],
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString"));
    } else {
      selectionString =
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString");
    }

    std::string const &drawString =
      IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "DrawString");

    if(!drawString.length() || !selectionString.length()){
      return;
    }

    Selections1D[selName] = Selection1D(
      selName,
      drawString,
      selectionString.c_str(),
      NBins,
      XBinLow,
      XBinUpper,
      Invalidate);

    std::cout << "[INFO]: Added Selection: " << Selections1D[selName]
      << std::endl;

}
void Add2DSelection(TXMLEngine &xmlengine,
  XMLDocPointer_t SelsDescriptNode){

    std::string const &selName =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"Name");

    if(selName.length() == 0){
      std::cout << "[ERROR]: Selection had no Name." << std::endl;
      return;
    }

    std::string const &selPreSelection =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"PreSelection");

    std::string const &selNXBins =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"NXBins");
    std::string const &selXBinLow =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinLow");
    std::string const &selXBinUpper =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinUpper");

    std::string const &selNYBins =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"NYBins");
    std::string const &selYBinLow =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"YBinLow");
    std::string const &selYBinUpper =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"YBinUpper");

    int NXBins;
    double XBinLow;
    double XBinUpper;

    int NYBins;
    double YBinLow;
    double YBinUpper;

    try{
      NXBins = std::stoi(selNXBins);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse NXBins: " << selNXBins << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinLow = std::stod(selXBinLow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse XBinLow: " << selXBinLow << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinUpper = std::stod(selXBinUpper);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse XBinUpper: " << selXBinUpper << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      NYBins = std::stoi(selNYBins);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse NYBins: " << selNYBins << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      YBinLow = std::stod(selYBinLow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse YBinLow: " << selYBinLow << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      YBinUpper = std::stod(selYBinUpper);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse YBinUpper: " << selYBinUpper << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    bool valid, Invalidate;
    Invalidate = PGUtils::str2bool(
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,
       "InvalidateCache"), valid);

    if(!valid){
      Invalidate = true;
    }

    std::string selectionString = "";
    if(selPreSelection.length() && PreSelections.count(selPreSelection)){
      selectionString = CombineWeightStrings(PreSelections[selPreSelection],
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString"));
    } else {
      selectionString =
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString");
    }

    std::string const &drawString =
      IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "DrawString");



    if(!drawString.length() || !selectionString.length()){
      return;
    }

    Selections2D[selName] = Selection2D(
      selName,
      drawString,
      selectionString.c_str(),
      NXBins,
      XBinLow,
      XBinUpper,
      NYBins,
      YBinLow,
      YBinUpper,
      Invalidate
      );

    std::cout << "[INFO]: Added Selection2D: " << Selections1D[selName]
      << std::endl;

}

void AddPlotSelectionDescriptionXML(TXMLEngine &xmlengine,
  XMLDocPointer_t PlotSelsNode){

  for(XMLNodePointer_t SelsDescriptNode = xmlengine.GetChild(PlotSelsNode);
    (SelsDescriptNode != NULL);
    SelsDescriptNode = xmlengine.GetNext(SelsDescriptNode)){

    if(std::string(xmlengine.GetNodeName(SelsDescriptNode)) == "Selection1D"){
      Add1DSelection(xmlengine, SelsDescriptNode);
    } else if(std::string(xmlengine.GetNodeName(SelsDescriptNode)) == "Selection2D"){
      Add2DSelection(xmlengine, SelsDescriptNode);
    } else {
      std::cout << "[WARN]: Found unexpected node: "
        << xmlengine.GetNodeName(SelsDescriptNode)
        << " in the Selection element." << std::endl;
      continue;
    }
  }
}

void InitSelectionsXML(std::string const &confFile){
  (void)Verbosity;

  TXMLEngine* xmlengine = new TXMLEngine;
  xmlengine->SetSkipComments();

  XMLNodePointer_t confSecNode =
    IOUtils::GetNamedChildElementOfDocumentRoot(*xmlengine,
      confFile, "Selections");
  if(!confSecNode){
    delete xmlengine; return;
  }

  for(XMLNodePointer_t PreSelsNode = xmlengine->GetChild(confSecNode);
    (PreSelsNode != NULL);
    PreSelsNode = xmlengine->GetNext(PreSelsNode)){
    if(std::string(xmlengine->GetNodeName(PreSelsNode)) != "PreSelections"){
      continue;
    }
    std::cout << "Found PreSelections element." << std::endl;
    AddPreSelectionDescriptionXML(*xmlengine,PreSelsNode);
    break;
  }

  for(XMLNodePointer_t PlotSelsNode = xmlengine->GetChild(confSecNode);
    (PlotSelsNode != NULL);
    PlotSelsNode = xmlengine->GetNext(PlotSelsNode)){
    if(std::string(xmlengine->GetNodeName(PlotSelsNode)) != "PlotSelections"){
      continue;
    }
    std::cout << "Found PlotSelections element." << std::endl;
    AddPlotSelectionDescriptionXML(*xmlengine,PlotSelsNode);
    break;

  }
  std::cout << "Done reading selections" << std::endl;
  delete xmlengine;
}
}
