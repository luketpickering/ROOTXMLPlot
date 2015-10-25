#ifndef PLOTTINGXMLBUILDABLE_HXX_SEEN
#define PLOTTINGXMLBUILDABLE_HXX_SEEN

class IXMLBuildable {
  bool IsValid;
  IXMLBuildable(){IsValid = false;}
public:
  bool Build(TXMLEngine&, XMLNode_Pointer_t const &) = 0;

};

bool operator!(IXMLBuildable const & xb){ return !xb.IsValid; }

namespace SelectionDescriptors {

  class Selection1D : public PlottingTypes::Selection, public IXMLBuildable {
  public:
    int NXBins;
    double XBinLow;
    double XBinUpper;
    bool DoPerUseXOffset;
    double PerUseXOffset;

    Selection1D() : Selection() {
      NBins = 0;
      XBinLow = 0;
      XBinUpper = 0;
      DoPerUseXOffset = false;
      PerUseXOffset = 0;
    }

    Selection1D(TXMLEngine&, XMLNode_Pointer_t const &);
    bool Build(TXMLEngine&, XMLNode_Pointer_t const &);

  };

  struct Selection2D : public Selection, public IXMLBuildable {

    int NXBins;
    double XBinLow;
    double XBinUpper;

    int NYBins;
    double YBinLow;
    double YBinUpper;

    bool Normalise;

    Selection2D() : Selection() {
      NXBins = 0;
      XBinLow = 0;
      XBinUpper = 0;

      NYBins = 0;
      YBinLow = 0;
      YBinUpper = 0;

    }

    Selection2D(TXMLEngine&, XMLNode_Pointer_t const &);
    bool Build(TXMLEngine&, XMLNode_Pointer_t const &);
  }
}

class Canvas : public IXMLBuildable {

};

class Pad : public IXMLBuildable {

};



#endif
