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
using namespace PlottingDefaults;

namespace PlottingSelections {

  std::map<std::string, SelectionDescriptors::Selection1D> Selections1D;
  std::map<std::string, SelectionDescriptors::Selection2D> Selections2D;
  std::map<std::string,std::string> PreSelections;

// //**********************************Globals*************************************

  SelectionDescriptors::Selection1D const * FindSel1D(std::string name){
    if(Selections1D.count(name)){
      return &Selections1D[name];
    }
    std::cout << "[WARN]: Couldn't find Selection1D named: " << name
      << std::endl;
    return 0;
  }

  SelectionDescriptors::Selection2D const * FindSel2D(std::string name){
    if(Selections2D.count(name)){
      return &Selections2D[name];
    }
    std::cout << "[WARN]: Couldn't find Selection2D named: " << name
      << std::endl;
    return 0;
  }


void AddPreSelectionDescriptionXML(XMLDocPointer_t PreSelsNode){

  TXMLEngine &xmlengine = XMLDoc::Get("Selections");

  for(auto SelsDescriptNode : XMLUtils::GetNamedChildNodes(xmlengine,
    PreSelsNode,"PreSelection")){

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

void AddSelectionXML(XMLDocPointer_t PlotSelsNode){

  TXMLEngine &xmlengine = XMLDoc::Get("Selections");

  for(auto SelsDescriptNode : XMLUtils::GetNamedChildNodes(xmlengine,
    PlotSelsNode,"Selection1D")){

    Selection1D sel1D(xmlengine,SelsDescriptNode);
    if(!!sel1D){
      Selections1D[sel1D.Name] = sel1D;
    }
  }

  for(auto SelsDescriptNode : XMLUtils::GetNamedChildNodes(xmlengine,
    PlotSelsNode,"Selection2D")){

    Selection2D sel2D(xmlengine,SelsDescriptNode);
    if(!!sel2D){
      Selections2D[sel2D.Name] = sel2D;
    }
  }
}

void InitSelectionsXML(std::string const &confFile){
  (void)Verbosity;

  XMLDoc::LoadDoc("Selections",confFile);
  TXMLEngine &xmlengine = XMLDoc::Get("Selections");

  XMLNodePointer_t confSecNode =
    IOUtils::GetNamedChildElementOfDocumentRoot(xmlengine,
      confFile, "SelectionDescriptors");
  if(!confSecNode){
    XMLDoc::Destroy("Selections");
    return;
  }

  XMLNodePointer_t PreSelsNode = GetFirstNamedChildNode(xmlengine,confSecNode,
    "PreSelections");

  if(PreSelsNode != NULL){
    AddPreSelectionDescriptionXML(PreSelsNode);
  }
  AddSelectionXML(confSecNode);

  std::cout << "Done reading selections" << std::endl;
  XMLDoc::Destroy("Selections");
}
}
