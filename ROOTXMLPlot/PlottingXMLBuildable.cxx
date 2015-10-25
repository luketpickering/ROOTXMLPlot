
namespace SelectionDescriptors {

  Selection1D::Selection1D() : Selection(), IXMLBuildable() {
    NBins = 0;
    XBinLow = 0;
    XBinUpper = 0;
    DoPerUseXOffset = false;
    PerUseXOffset = 0;
  }

  Selection1D(TXMLEngine& xmlengine, XMLNode_Pointer_t const & node){
    IsValid = Build(xmlengine,node);
  }

  bool Build(TXMLEngine& xmlengine, XMLNode_Pointer_t const & node){

   Name = IOUtils::GetAttrValue(xmlengine, node,"Name");

    if(!Name.length()){
      std::cout << "[ERROR]: Selection had no Name." << std::endl;
      return false;
    }

    std::vector<std::string> PreSelsUsed =
      PGUtils::SplitStringByDelim(
        IOUtils::GetAttrValue(xmlengine, node,"PreSelections"), " ");

    try{
      NXBins = GetIntAttrValue(xmlengine,node,"NXBins", true);
      XBinLow = GetDoubleAttrValue(xmlengine,node,"XBinLow", true);
      XBinUpper = GetDoubleAttrValue(xmlengine,node,"XBinUpper", true);
    } catch (...){
      return false;
    }

    bool valid;
    InvalidateCache =
      GetBoolAttrValue(xmlengine, node, "InvalidateCache", valid);
    if(!valid){
      InvalidateCache = true;
    }

    std::string SelectionString = "";

    for(auto const & preSel : PreSelsUsed) {
      if(preSel.length() && PreSelections.count(preSel)){
        SelectionString = CombineWeightStrings(PreSelections[preSel],
          SelectionString);
      }
    }
    SelectionString =
        IOUtils::GetChildNodeContent(xmlengine, node, "SelString");

    DrawCommand =
      IOUtils::GetChildNodeContent(xmlengine, node, "DrawString");

    if((!DrawCommand.length()) || (!SelectionString.length())){
      return false;
    }
    Cut = SelectionString.c_str();
    return true;
  }

  Selection2D::Selection2D() : Selection(), IXMLBuildable() {
    NXBins = 0;
    XBinLow = 0;
    XBinUpper = 0;
    NYBins = 0;
    YBinLow = 0;
    YBinUpper = 0;
  }

  Selection2D(TXMLEngine& xmlengine, XMLNode_Pointer_t const & node){
    IsValid = Build(xmlengine,node);
  }

  bool Build(TXMLEngine& xmlengine, XMLNode_Pointer_t const & node){

    //We know that Selection1D looks for a subset of this information so
    //use that as a proxy to get it
    Selection1D sel1D(xmlengine,node);
    if(!sel1D){
      return false;
    }

    try{
      NYBins = GetIntAttrValue(xmlengine,node,"NYBins", true);
      YBinLow = GetDoubleAttrValue(xmlengine,node,"YBinLow", true);
      YBinUpper = GetDoubleAttrValue(xmlengine,node,"YBinUpper", true);
    } catch (...){
      return false;
    }

    Name = sel1D.Name;
    NXBins = sel1D.NXBins;
    XBinLow = sel1D.XBinLow;
    XBinUpper = sel1D.XBinUpper;
    DrawCommand = sel1D.DrawCommand;
    Cut = sel1D.Cut;
    InvalidateCache = sel1D.InvalidateCache;

    return true;
  }
}
