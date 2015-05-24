#include <iostream>

#include "TMath.h"

#include "PlottingUtils.hxx"

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
void MakePDF(TH1* th){
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
  bool logx){
  TCanvas* canv = new TCanvas(CanvName.c_str(),"");
  canv->SetTopMargin(0);
  canv->SetBottomMargin(0);
  canv->SetLeftMargin(0);
  canv->SetRightMargin(0);
  TVirtualPad* P0 = canv->GetPad(0);

  P0->SetTopMargin(PlottingDefaults::Pad0YMargin);
  P0->SetBottomMargin(PlottingDefaults::Pad0YMargin);
  P0->SetLeftMargin(PlottingDefaults::Pad0XMargin);
  P0->SetRightMargin(PlottingDefaults::Pad0XMargin);
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

  for(std::vector<TH1*>::iterator main_plots_it = MainPlotSeries.begin();
    main_plots_it != MainPlotSeries.end(); ++main_plots_it){

    (*main_plots_it)->Draw((main_plots_it==MainPlotSeries.begin())?"EHIST":
      "EHIST SAME");
    leg->AddEntry(*main_plots_it,(*main_plots_it)->GetTitle(), "l");
  }

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
