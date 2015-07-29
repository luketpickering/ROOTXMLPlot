#include <ctype.h>
#include <algorithm>
#include <functional>

#include "TMath.h"

#include "PlottingUtils.hxx"

namespace {
  int Verbosity = 0;
}

namespace PlottingUtils {

//**********START*****************TH1 Utils***********************************/

///Adjusts input TAxis to to log bins. Left edge is either the current or if
///less than epsilon = 1E-12, it is set to epsilon.
bool BinLogX(TAxis* axis){
  static double const epsilon = 1E-12;
  int const bins = axis->GetNbins();

  double const from = axis->GetXmin();
  double const to = axis->GetXmax();
  if (from < epsilon) return false;
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  const Double_t factor = TMath::Power(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
  return true;
}

///Scales scaled bin-by-bin to ref, errors on both are scaled by the bin
///value of ref as well.
std::pair<TH1*,TH1*>
ScaleToRef(const TH1& scaled, const TH1& ref){
  if(ref.GetNbinsX() != scaled.GetNbinsX()){
    std::cerr << "Error: bin num not consistent: Ref("
      << ref.GetNbinsX() << "), scaled(" << scaled.GetNbinsX()
      << ")." << std::endl;
    return std::pair<TH1*,TH1*>(0,0);
  }
  TH1* refClone = NoSaveClone(&ref);
  TH1* scaledClone = NoSaveClone(&scaled);

  for(int i = 1; i < ref.GetNbinsX() + 1; ++i){
    float rbc = ref.GetBinContent(i);
    refClone->SetBinContent(i,1.0);
    refClone->SetBinError(i,ref.GetBinError(i)/rbc);
    scaledClone->SetBinContent(i,scaled.GetBinContent(i)/rbc);
    scaledClone->SetBinError(i,scaled.GetBinError(i)/rbc);
  }
  return std::pair<TH1*,TH1*>(scaledClone,refClone);
}

///Ease of use interface for ScaleToRef(const TH1& scaled, const TH1& ref)
std::pair<TH1*,TH1*>
ScaleToRef(const TH1* scaled, const TH1* ref){
  return ScaleToRef(*scaled,*ref);
}

///Divides the contents and error of each bin by with width of the
//bin
void BinWidthNormalise(TH1* th){
  for(int i = 1; i < th->GetNbinsX() + 1; ++i){
    float bc = th->GetBinContent(i);
    bc /= th->GetBinWidth(i);
    th->SetBinContent(i,bc);

    float be = th->GetBinError(i);
    be /= th->GetBinWidth(i);
    th->SetBinError(i,be);
  }
}

///Divides the contents and error of each bin by the total populations,
///including the under and overflow bins.
void AreaNormalise(TH1* th){
  double integral = th->Integral(0, th->GetNbinsX() + 1);

  for(int i = 1; i < th->GetNbinsX() + 1; ++i){
    float bc = th->GetBinContent(i);
    bc /= integral;
    th->SetBinContent(i,bc);

    float be = th->GetBinError(i);
    be /= integral;
    th->SetBinError(i,be);
  }
}

///Combines the functions of BinWidthNormalise and AreaNormalise and adds
///some debugging output
void MakePDF(TH1* th, bool reallyMakePDF){
  if(!reallyMakePDF){return;}
  double integral = th->Integral(0, th->GetNbinsX() + 1);

  std::cout << th->GetName() << "(" << integral << ")" << std::endl;

  for(int i = 1; i < th->GetNbinsX() + 1; ++i){
    float bc = th->GetBinContent(i);
    float pbc = bc;
    bc /= (th->GetBinWidth(i)*integral);
    th->SetBinContent(i,bc);

    float be = th->GetBinError(i);
    float pbe = be;
    be /= (th->GetBinWidth(i)*integral);
    th->SetBinError(i,be);

    if(Verbosity > 1){
      std::cout << "\t" << "(" << th->GetBinLowEdge(i) << ")[ "
        << pbc << " +- " << pbe << "  -> " << bc << " +- " << be << "]("
        << th->GetBinLowEdge(i+1) << ") bw ("
        << (th->GetBinLowEdge(i+1) - th->GetBinLowEdge(i))
        << ":" << th->GetBinWidth(i) << ")" << std::endl;
    }
  }
}

TH1* MakePDFClone(TH1 const * const &inp){
  if(!inp){
    std::cout << "[WARN]: Zeroed pointer passed to MakePDFClone" << std::endl;
    return 0;
  }
  TH1* clone = static_cast<TH1*>(inp->Clone());
  clone->SetName((std::string(clone->GetName())+"_pdfc").c_str());
  MakePDF(clone);
  return clone;
}

//***********END******************TH1 Utils***********************************/

//**********START***************Factory Utils**********************************/

///Creates a standard canvas with two sub TPads.
TCanvas* SplitCanvasFactory(std::string const &CanvName){
  TCanvas* canv = new TCanvas(CanvName.c_str(),"");
  canv->SetTopMargin(0);
  canv->SetBottomMargin(0);
  canv->SetLeftMargin(0);
  canv->SetRightMargin(0);
  TVirtualPad* P0 = canv->GetPad(0);

  P0->Divide(1,2,0,0);

  TVirtualPad* P1 = P0->GetPad(1);

  P1->SetTopMargin(PlottingDefaults::Pad1YMargin);
  P1->SetBottomMargin(PlottingDefaults::Pad1YMargin);

  P1->SetLeftMargin(PlottingDefaults::Pad1XMargin);
  P1->SetRightMargin(PlottingDefaults::Pad1XMargin);
  P1->SetLogx(PlottingDefaults::DoLogx);

  TVirtualPad* P2 = P0->GetPad(2);

  P2->SetTopMargin(PlottingDefaults::Pad1YMargin);
  P2->SetBottomMargin(PlottingDefaults::Pad1YMargin);

  P2->SetLeftMargin(PlottingDefaults::Pad1XMargin);
  P2->SetRightMargin(PlottingDefaults::Pad1XMargin);
  P2->SetLogx(PlottingDefaults::DoLogx);

  return canv;
}

TCanvas* CanvasFactory(std::string const &CanvName,
  bool logx, bool tight){
  TCanvas* canv = new TCanvas(CanvName.c_str(),"");
  canv->SetTopMargin(0);
  canv->SetBottomMargin(0);
  canv->SetLeftMargin(0);
  canv->SetRightMargin(0);
  TVirtualPad* P0 = canv->GetPad(0);

  if(!tight){
    P0->SetTopMargin(PlottingDefaults::Pad0YMargin);
    P0->SetBottomMargin(PlottingDefaults::Pad0YMargin);
    P0->SetLeftMargin(PlottingDefaults::Pad0XMargin);
    P0->SetRightMargin(PlottingDefaults::Pad0XMargin);
  } else {
    P0->SetRightMargin(0.03);
    P0->SetTopMargin(0.02);
    P0->SetLeftMargin(0.12);
    P0->SetBottomMargin(0.15);
  }
  P0->SetBorderMode(-1);
  P0->SetBorderSize(0);
  P0->SetLogx(logx);

  return canv;
}

TLatex* TitleFactory(std::string const &TitleText,
  double x0, double y0, float TextSize){
  TLatex* title = new TLatex(x0,y0, TitleText.c_str());
  title->SetTextSize(TextSize);
  title->SetNDC();
  return title;
}

///Where legends are made
TLegend* LegendFactory(float lx0,float ly0,float lx1,float ly1){
  TLegend* leg = new TLegend(lx0,ly0,lx1,ly1);
  leg->SetFillColor(kWhite);
  return leg;
}

//***********END****************Factory Utils**********************************/

//**********START**********Pretty Defaults Functions***************************/

void SetUpperPaneHistoDefaults(TH1* hist, char const * ytitle,
  char const * xtitle){
  hist->GetYaxis()->SetTitleOffset(PlottingDefaults::TopPaneYTitleOffset);
  hist->GetYaxis()->SetLabelOffset(PlottingDefaults::TopPaneYLabelOffset);
  hist->GetXaxis()->SetTitleOffset(PlottingDefaults::TopPaneXTitleOffset);
  hist->GetXaxis()->SetLabelOffset(PlottingDefaults::TopPaneXLabelOffset);
  hist->GetYaxis()->SetTitle(ytitle);
  hist->GetXaxis()->SetTitle(xtitle);
}

void SetLowerPaneHistoDefaults(TH1* hist, char const * ytitle,
  char const * xtitle){
  hist->GetYaxis()->SetTitleOffset(PlottingDefaults::BottomPaneYTitleOffset);
  hist->GetYaxis()->SetLabelOffset(PlottingDefaults::BottomPaneYLabelOffset);
  hist->GetXaxis()->SetTitleOffset(PlottingDefaults::BottomPaneXTitleOffset);
  hist->GetXaxis()->SetLabelOffset(PlottingDefaults::BottomPaneXLabelOffset);
  hist->GetYaxis()->SetTitle(ytitle);
  hist->GetXaxis()->SetTitle(xtitle);
}

void SetOnePaneHistDefaults(TH1* hist,
  char const * Ytitle,
  char const * Xtitle){
  hist->GetYaxis()->SetTitle(Ytitle);
  hist->GetXaxis()->SetTitle(Xtitle);
  hist->GetYaxis()->SetTitleOffset(PlottingDefaults::YTitleOffset);
  hist->GetYaxis()->SetLabelOffset(PlottingDefaults::YLabelOffset);
  hist->GetXaxis()->SetTitleOffset(PlottingDefaults::XTitleOffset);
  hist->GetXaxis()->SetLabelOffset(PlottingDefaults::XLabelOffset);

  hist->GetYaxis()->SetLabelSize(PlottingDefaults::AxisLabelSize);
  hist->GetXaxis()->SetLabelSize(PlottingDefaults::AxisLabelSize);
  hist->GetYaxis()->SetTitleSize(PlottingDefaults::AxisTitleSize);
  hist->GetXaxis()->SetTitleSize(PlottingDefaults::AxisTitleSize);

}

//***********END***********Pretty Defaults Functions***************************/

// Notes on Usage:
//   The TH1*s wil be Added to the legend and plotted as is, you need to set up
//   any styling, titles, and axes ranges beforehand.
//   (myth->GetYaxis()->SetUserRange(low,high).

//   All of the margin and width values are in units of drawing area fraction, x or
//   y.

//   Tips: The default Inset R and T Margins will just keep the inset off the top
//   and right axes tics of a T2K plot. If you set these to 0.0 then the Inset
//   will appear flush with the main plot.

//   Currently the inset can only be anchored to the top right. I will add final
//   positioning later.
TCanvas* GetInsetCanvas(
  std::vector<TH1*> &MainPlotSeries,
  std::vector<TH1*> &InsetSeries,
  float MainXMargin, float MainYMargin,
  float InsetWidth, float InsetHeight,
  float LegendWidth, float LegendHeight,
  bool LegendOnLeft, bool LegendOnTop,
  float InsetLMargin, float InsetRMargin,
  float InsetTMargin, float InsetBMargin,
  float LegendClosestXMargin, float LegendClosestYMargin){

  static std::stringstream ss("");
  static int i = 0;
  ss << (i++);
  TCanvas* canv = new TCanvas(("CallMeWhatYouWill"+ss.str()).c_str(),"");

  canv->SetTopMargin(0);
  canv->SetBottomMargin(0);
  canv->SetLeftMargin(0);
  canv->SetRightMargin(0);
  canv->SetBorderMode(-1);
  canv->SetBorderSize(0);
  TVirtualPad* P0 = canv->GetPad(0);

  P0->SetTopMargin(MainYMargin);
  P0->SetBottomMargin(MainYMargin);
  P0->SetLeftMargin(MainXMargin);
  P0->SetRightMargin(MainXMargin);
  P0->SetBorderMode(-1);
  P0->SetBorderSize(0);

  float LegX0,LegX1,LegY0,LegY1;

  if(LegendOnLeft){
    LegX0 = MainXMargin + LegendClosestXMargin;
    LegX1 = MainXMargin + LegendWidth + LegendClosestXMargin;
  } else {
    LegX0 = 1.0 - MainXMargin - LegendWidth - LegendClosestXMargin;
    LegX1 = 1.0 - MainXMargin - + LegendClosestXMargin;
  }
  if(LegendOnTop){
    LegY0 = 1.0 - MainYMargin - LegendHeight - LegendClosestYMargin;
    LegY1 = 1.0 - MainYMargin - LegendClosestYMargin;
  } else {
    LegY0 = MainYMargin + LegendClosestYMargin;
    LegY1 = MainYMargin + LegendHeight + LegendClosestYMargin;
  }

  TLegend* leg = new TLegend(LegX0,LegY0,LegX1,LegY1);
  leg->SetFillColor(kWhite);

  DrawTHs(MainPlotSeries);
  AddToLegend(MainPlotSeries, leg);

  // for(std::vector<TH1*>::iterator main_plots_it = MainPlotSeries.begin();
  //   main_plots_it != MainPlotSeries.end(); ++main_plots_it){

  //   (*main_plots_it)->Draw((main_plots_it==MainPlotSeries.begin())?"EHIST":
  //     "EHIST SAME");
  //   leg->AddEntry(*main_plots_it,(*main_plots_it)->GetTitle(), "l");
  // }

  TPad *insetPad = new TPad("inset","",
    (1.0 - MainXMargin - InsetWidth - InsetLMargin - InsetRMargin),
    (1.0 - MainYMargin - InsetHeight - InsetTMargin - InsetBMargin),
    (1.0 - MainXMargin - InsetRMargin),
    (1.0 - MainYMargin - InsetTMargin));
  insetPad->SetRightMargin(InsetRMargin);
  insetPad->SetLeftMargin(InsetLMargin);
  insetPad->SetTopMargin(InsetTMargin);
  insetPad->SetBottomMargin(InsetBMargin);

  insetPad->Draw();
  insetPad->cd();

  for(std::vector<TH1*>::iterator inset_plots_it = InsetSeries.begin();
    inset_plots_it != InsetSeries.end(); ++inset_plots_it){

    (*inset_plots_it)->Draw((inset_plots_it==InsetSeries.begin())?"EHIST":
      "EHIST SAME");
  }
  P0->cd();
  leg->Draw();

  return canv;
}

}

namespace IOUtils {
  XMLDocPointer_t GetNamedChildElementOfDocumentRoot(TXMLEngine &xmlengine,
      std::string const &FileName, std::string const & ElementName){

    XMLDocPointer_t xmldoc = xmlengine.ParseFile(FileName.c_str());

    if(!xmldoc){
      std::cout << "In File: " << FileName << std::endl;
      std::cout << "[ERROR]: TXMLEngine could not parse file." << std::endl;
      return XMLDocPointer_t();
    }

    XMLNodePointer_t rootNode = xmlengine.DocGetRootElement(xmldoc);

    for(XMLNodePointer_t confSecNode = xmlengine.GetChild(rootNode);
      (confSecNode != NULL);
      confSecNode = xmlengine.GetNext(confSecNode)){
      if(std::string(xmlengine.GetNodeName(confSecNode)) != ElementName){
        continue;
      }
      return confSecNode;
    }
    std::cout << "[WARN]: In File: " << FileName << " could not find a root "
      "element child named: " << ElementName << "." << std::endl;
    return XMLDocPointer_t();
  }

  std::string GetAttrValue(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &attrName, bool expectExclusive){

    for(XMLAttrPointer_t nodeAttr_ptr = xmlengine.GetFirstAttr(node);
        (nodeAttr_ptr != NULL);
        nodeAttr_ptr = xmlengine.GetNextAttr(nodeAttr_ptr)){
      if(std::string(xmlengine.GetAttrName(nodeAttr_ptr)) != attrName){
        if(expectExclusive){
          std::cout << "[WARN]: Unexpected attribute on "
            << xmlengine.GetNodeName(node) << " node: "
            << xmlengine.GetAttrName(nodeAttr_ptr) << std::endl;
        }
        continue;
      }
      return xmlengine.GetAttrValue(nodeAttr_ptr);
    }
    return std::string();
  }

  std::string GetNodeContent(TXMLEngine &xmlengine, XMLDocPointer_t node){
    std::string elementString;
    try {
      if(!xmlengine.GetNodeContent(node)){
        throw std::logic_error("Empty Element.");
      }
      elementString = xmlengine.GetNodeContent(node);
    } catch(std::logic_error& le){
      std::cout << "[WARN]: Found empty element: "
        << xmlengine.GetNodeName(node) << "." << std::endl;
      return "";
    }
    return elementString;
  }

  std::string GetChildNodeContent(TXMLEngine &xmlengine, XMLDocPointer_t node,
    std::string const &ChildNodeName){

    for(XMLNodePointer_t childNode_ptr = xmlengine.GetChild(node);
      (childNode_ptr != NULL);
      childNode_ptr = xmlengine.GetNext(childNode_ptr)){

      if(std::string(xmlengine.GetNodeName(childNode_ptr)) ==
          ChildNodeName){
        return GetNodeContent(xmlengine, childNode_ptr);
      }
    }
    return std::string();
  }

  EColor ParseRootColor(std::string const &colorString){
    if(!colorString.length()){ return kBlack; }
    // enum EColor { kWhite =0,   kBlack =1,   kGray=920,
    //         kRed   =632, kGreen =416, kBlue=600, kYellow=400, kMagenta=616, kCyan=432,
    //         kOrange=800, kSpring=820, kTeal=840, kAzure =860, kViolet =880, kPink=900 };

    std::string intString = colorString;
    EColor BaseColor = kWhite;
    int offset = 0;

    if(colorString[0] == 'k'){
      size_t eoc = colorString.find_first_of(" +");

      std::string cSubStr = colorString.substr(0,eoc);

      size_t sois = colorString.substr(eoc+1).find_first_not_of(" +");
      intString = colorString.substr(eoc+1+sois);

      if((eoc == std::string::npos) || (sois == std::string::npos)){
        intString = "0";
      }

      std::cout << "[INFO]: Found Color String: \"" << cSubStr
        << "\" (" << eoc << "), OffSet String \"" << intString << "\" ("
        << sois << ")" << std::endl;

      std::transform(cSubStr.begin(), cSubStr.end(), cSubStr.begin(), tolower);

      if(cSubStr == "kwhite"){
        BaseColor = kWhite;
      } else if (cSubStr == "kblack"){
        BaseColor = kBlack;
      } else if (cSubStr == "kgray"){
        BaseColor = kGray;
      } else if (cSubStr == "kred"){
        BaseColor = kRed;
      } else if (cSubStr == "kgreen"){
        BaseColor = kGreen;
      } else if (cSubStr == "kblue"){
        BaseColor = kBlue;
      } else if (cSubStr == "kyellow"){
        BaseColor = kYellow;
      } else if (cSubStr == "kmagenta"){
        BaseColor = kMagenta;
      } else if (cSubStr == "kcyan"){
        BaseColor = kCyan;
      } else if (cSubStr == "korange"){
        BaseColor = kOrange;
      } else if (cSubStr == "kspring"){
        BaseColor = kSpring;
      } else if (cSubStr == "kteal"){
        BaseColor = kTeal;
      } else if (cSubStr == "kazure"){
        BaseColor = kAzure;
      } else if (cSubStr == "kviolet"){
        BaseColor = kViolet;
      } else if (cSubStr == "kpink"){
        BaseColor = kPink;
      } else {
        std::cout << "[ERROR]: Couldn't parse EColor string: " << cSubStr
          << std::endl;
      }
    }

    offset = std::stoi(intString);
    return EColor(BaseColor + offset);
  }

}

std::ostream & operator<<(std::ostream& os, PlottingTypes::Sample const & sm){
  os << "{";
  int added = 0;
  if(sm.Name.length()){
    os << "Name: " << sm.Name;
    added++;
  }
  if(sm.Flux.length()){
    os << (added?", ":"") << "Flux: " << sm.Flux;
    added++;
  }
  if(sm.FileLocation.length()){
    os << (added?", ":"") << "FileLocation: " << sm.FileLocation;
    added++;
  }
  if(sm.Description.length()){
    os << (added?", ":"") << "Description: " << sm.Description;
    added++;
  }
  if(sm.InvalidateCache){
    os << (added?", ":"") << "InvalidateCache: true";
    added++;
  }
  os << " }";

  return os;
}

std::ostream& operator<<(std::ostream& os,
  PlottingTypes::Selection1D sel){
    os << "{ Name : " << sel.Name
    << "\n\tDraw: \"" << sel.DrawString
    << "\"\n\tCut: \""
    << PlottingUtils::ConfineWidthStringWithIndent("\t\t\t",sel.Cut.GetTitle(),70)
    <<"\"\n\tBins: [ " << sel.XBinLow << " .. " << sel.NBins << " .. "
    << sel.XBinUpper << " ] }";
  return os;
}

std::ostream& operator<<(std::ostream& os,
  PlottingTypes::Selection2D sel){
    os << "{ Name : " << sel.Name
    << "\n\tDraw: \"" << sel.DrawString
    << "\"\n\tCut: \""
    << PlottingUtils::ConfineWidthStringWithIndent("\t\t\t",sel.Cut.GetTitle(),70)
    <<"\"\n\tBins: [ " << sel.XBinLow << " .. " << sel.NXBins << " .. "
    << sel.XBinUpper << " ]"
    <<"\"\n\tBins: [ " << sel.YBinLow << " .. " << sel.NYBins << " .. "
    << sel.YBinUpper << " ] }";
  return os;
}

std::ostream& operator<<(std::ostream& os,
  PlottingTypes::SeriesDescriptor::SeriesType st){
  switch(st){
    case PlottingTypes::SeriesDescriptor::kUnSet:{
      return os << "N/A";
    }
    case PlottingTypes::SeriesDescriptor::kTH1:{
      return os << "1D Histogram";
    }
    case PlottingTypes::SeriesDescriptor::kTH1RatioNumer:{
      return os << "Ratio of 1D Histograms---Numerator";
    }
    case PlottingTypes::SeriesDescriptor::kTH1RatioDenom:{
      return os << "Ratio of 1D Histograms---Denominator";
    }
    case PlottingTypes::SeriesDescriptor::kTH2Profile:{
      return os << "Profile of 2D Histogram";
    }
    case PlottingTypes::SeriesDescriptor::kTH2:{
      return os << "2D Histogram";
    }
    default:{
      return os << "N/A";
    }
  }
}

std::ostream& operator<<(std::ostream& os, PlottingTypes::SeriesDescriptor sd){
  os << "{ Gen: " << sd.Generator << ", Sample: " << sd.Sample
    << ", Sel: " << sd.Selection << ", Type : " << sd.Type
    << ", LineStyle : {C : " << sd.LineColor
    << " , W : " <<  sd.LineWidth << " , S : " << sd.LineStyle << " } }";
  return os;
}

std::ostream& operator<<(std::ostream& os,
  PlottingTypes::PlotDescriptor pd){

  os << "{ Name : " << pd.Name << ", Title : " << pd.Title << ", Print: "
    << pd.PrintName << " Series : { \n\t";
  for(size_t i = 0; i < pd.Series.size(); ++i){
    os << pd.Series[i] << ((i == (pd.Series.size()-1))?"}":",\n\t");
  }
  os << ", Logs: X:" << (pd.Logs&1) << ", Y:" <<  (pd.Logs&2) << ", Z: "
    << (pd.Logs&4) << "}";
  return os;
}
