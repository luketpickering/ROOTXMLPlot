#ifndef PLOTTINGUTILS_HXX_SEEN
#define PLOTTINGUTILS_HXX_SEEN

#include <initializer_list>
#include <vector>
#include <sstream>
#include <iomanip>

#include "TAxis.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"

namespace {
  int Verbosity = 0;
}

namespace PlottingDefaults {
  constexpr double TopPaneXLabelOffset = 50;
  constexpr double TopPaneXTitleOffset = 50;
  constexpr double TopPaneYLabelOffset = 0.01;
  constexpr double TopPaneYTitleOffset = 0.45;

  constexpr double BottomPaneXLabelOffset = 0.005;
  constexpr double BottomPaneXTitleOffset = 0.7;
  constexpr double BottomPaneYLabelOffset = 0.001;
  constexpr double BottomPaneYTitleOffset = 0.45;

  constexpr double YLabelOffset = 0.004;
  constexpr double YTitleOffset = 1.5;
  constexpr double XLabelOffset = 0.005;
  constexpr double XTitleOffset = 1.0;

  constexpr double Pad0YMargin = 0.135;
  constexpr double Pad0XMargin = 0.135;

  constexpr double Pad1YMargin = 0.12;
  constexpr double Pad1XMargin = 0.12;

  constexpr double P0LegendWidth = 0.2;
  constexpr double P0LegendHeight = 0.2;

  constexpr double P1LegendWidth = 0.2;
  constexpr double P1LegendHeight = 0.36;

  constexpr double TitleX0 = 0.3;
  constexpr double TitleY0 = 0.95;
  constexpr double TitleTextSize = 0.04;

  constexpr bool DoLogx = false;
  constexpr bool DoLogy = false;

  constexpr Float_t MinYVal = ((!DoLogy)?0.0:1E-5);
  constexpr Float_t YStretchFactor = 1.2;

  constexpr char const * DPhiTTitle = "#delta#phi_{T} [degrees]";
  constexpr char const * DPhiT_PDF_Title = "#frac{dN}{Nd#delta#phi_{T}}";

  constexpr char const * DAlphaTTitle = "#delta#alpha_{T} [degrees]";
  constexpr char const * DAlphaT_PDF_Title = "#frac{dN}{Nd#delta#alpha_{T}}";

  constexpr char const * DPTTitle = "#delta p_{T} [MeV/C]";
  constexpr char const * DPT_PDF_Title = "#frac{dN}{Nd#delta p_{T}}";

}

namespace PlottingUtils {

//**********START*****************TH1 Utils***********************************/

///Adjusts input TAxis to to log bins. Left edge is either the current or if
///less than epsilon = 1E-12, it is set to epsilon.
bool BinLogX(TAxis* axis);

///Makes a clone of a T and sets it to be not saved.
///It is obviously, therefore, your job to delete it
template<typename T>
T* NoSaveClone(T const * const &inp){
  T* clone = static_cast<T*>(inp->Clone());
  clone->SetName((std::string(clone->GetName())+"_c").c_str());
  clone->SetDirectory(0);
  return clone;
}


///Scales scaled bin-by-bin to ref, errors on both are scaled by the bin
///value of ref as well.
std::pair<TH1*,TH1*>
ScaleToRef(const TH1& scaled, const TH1& ref);

///Ease of use interface for ScaleToRef(const TH1& scaled, const TH1& ref)
std::pair<TH1*,TH1*>
ScaleToRef(const TH1* scaled, const TH1* ref);

///Divides the contents and error of each bin by with width of the
//bin
void BinWidthNormalise(TH1* th);

///Divides the contents and error of each bin by the total populations,
///including the under and overflow bins.
void AreaNormalise(TH1* th);

///Combines the functions of BinWidthNormalise and AreaNormalise and adds
///some debugging output
void MakePDF(TH1* th);

inline double GetMaximumBinContents(TH1 const * inp){
  return inp->GetBinContent(inp->GetMaximumBin());
}
//***********END******************TH1 Utils***********************************/

//**********START***************Factory Utils**********************************/

///Creates a standard canvas with two sub TPads.
TCanvas* SplitCanvasFactory(std::string const &CanvName="SplitCanvas");

TCanvas* CanvasFactory(std::string const &CanvName="Canvas",
  bool logx=PlottingDefaults::DoLogx);

TLatex* TitleFactory(std::string const &TitleText="TitleBeHere",
  double x0=PlottingDefaults::TitleX0,
  double y0=PlottingDefaults::TitleY0,
  float TextSize=PlottingDefaults::TitleTextSize);

///Where legends are made
TLegend* LegendFactory(float lx0,float ly0,float lx1,float ly1);

//***********END****************Factory Utils**********************************/

//**********START**********Pretty Defaults Functions***************************/

void SetUpperPaneHistoDefaults(TH1* hist, char const * ytitle="",
  char const * xtitle="");

void SetLowerPaneHistoDefaults(TH1* hist, char const * ytitle="",
  char const * xtitle="");

void SetOnePaneHistDefaults(TH1* hist,
  char const * Ytitle=PlottingDefaults::DPhiT_PDF_Title,
  char const * Xtitle=PlottingDefaults::DPhiTTitle);

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
  float MainXMargin=0.15, float MainYMargin=0.15,
  float InsetWidth=0.4, float InsetHeight=0.3,
  float LegendWidth=0.25, float LegendHeight=0.25,
  bool LegendOnLeft=false, bool LegendOnTop=false,
  float InsetLMargin=0.135, float InsetRMargin=0.018,
  float InsetTMargin=0.02, float InsetBMargin=0.135,
  float LegendClosestXMargin=0.018, float LegendClosestYMargin=0.05);

inline void MakePDF(std::initializer_list<TH1 *> const histos){
  for(auto const & histo : histos){
    MakePDF(histo);
  }
}

inline Float_t GetLargestYVal(std::initializer_list<TH1 *> const histos){
  Float_t max = 0;
  for(auto const & histo : histos){
    Float_t esto = GetMaximumBinContents(histo) ;
    max = esto > max ? esto : max;
  }
  return max;
}

inline void SetUserRangeToMaxYVal(std::initializer_list<TH1*> histos,
  Float_t minval=0xDEADDEAD){
  Float_t ymin = (minval==0xDEADDEAD)?PlottingDefaults::MinYVal:minval;
  Float_t max = GetLargestYVal(histos);
  for(auto & histo : histos){
    histo->GetYaxis()->SetRangeUser(ymin,
      max*PlottingDefaults::YStretchFactor);
  }
}

//Splits an input string over multiple lines of a given width. Does't try
//anything clever so will happily split a word in half.
inline std::string ConfineWidthStringWithIndent(
  const char* indent, std::string inp,
  int width=50,
  const char* firstindent=""){

  std::stringstream ss("");

  ss << firstindent;
  int sslen = 0;
  while(true){
    ss << inp.substr(sslen,width);
    if((sslen+width) < inp.length()){
      ss << "\n" << indent;
      sslen += width;
    } else {
      break;
    }
  }
  return ss.str();
}

inline std::string GetIntegralAsString(TH1 const * hist){
  std::stringstream ss("");
  ss << std::setprecision(2) << hist->Integral();
  return ss.str();
}

inline void AddCountToTHTitle(
  std::initializer_list< std::pair< char const*, TH1*> >histos){
  for(auto const &namepair : histos){
    namepair.second->SetTitle((std::string(namepair.first) + " ("
      + GetIntegralAsString(namepair.second) + ")").c_str());
  }
}

inline void DrawTHs(std::initializer_list<TH1*> histos,
  std::string drawOpt="EHIST"){
  for(std::initializer_list<TH1*>::iterator th_it = histos.begin();
    th_it != histos.end(); ++th_it){
    (*th_it)->Draw((drawOpt+(th_it!=histos.begin()?" SAME":"")).c_str());
  }
}

inline void AddToLegend(std::initializer_list<TH1*> histos, TLegend* leg,
  char const * drawOpt="l"){
  for(auto const &hist : histos){
    leg->AddEntry(hist,hist->GetTitle(),drawOpt);
  }
}

inline TStyle* MakeT2KStyle(){
  TStyle *t2kStyle = new TStyle("T2K","T2K approved plots style");

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20,26);

  // use large Times-Roman fonts
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.15);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.07,"x");
  t2kStyle->SetTitleSize(0.07,"x");
  t2kStyle->SetLabelSize(0.07,"y");
  t2kStyle->SetTitleSize(0.07,"y");
  t2kStyle->SetLabelSize(0.07,"z");
  t2kStyle->SetTitleSize(0.07,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.1);
  t2kStyle->SetTitleFont(132,"pad");

  // use bold lines and markers
  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  // Add a greyscale palette for 2D plots
  int ncol=50;
  double dcol = 1./float(ncol);
  double gray = 1;
  TColor **theCols = new TColor*[ncol];
  for(int i = 0; i < ncol; ++i) {
    theCols[i] = new TColor(999-i,0.0,0.7,0.7);
  }
  for(int j = 0; j < ncol; j++) {
    theCols[j]->SetRGB(gray,gray,gray);
    gray -= dcol;
  }
  int ColJul[100];
  for(int i=0; i<100; i++) {
    ColJul[i]=999-i;
  }
  t2kStyle->SetPalette(ncol,ColJul);

  // Define a nicer color palette (red->blue)
  // Uncomment these lines for a color palette (default is B&W)
  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  t2kStyle->SetNumberContours(NCont);
  return t2kStyle;
}

}

#endif
