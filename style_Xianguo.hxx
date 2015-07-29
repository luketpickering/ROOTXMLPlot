#ifndef STYLE_H
#define STYLE_H

#include "TVirtualPad.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH1.h"
#include "THStack.h"
#include "TPaveText.h"

class style
{
 public:
  static TGraphAsymmErrors* ScaleGraph(const TGraphAsymmErrors * gin, const  TGraphAsymmErrors * gref);
  static void GraphMinMaxXY(const TGraphAsymmErrors * gr, Double_t & xmin, Double_t & xmax, Double_t &ymin, Double_t &ymax );
  static void PadSetup(TPad *currentPad, const Double_t currentLeft=0.12, const Double_t currentTop=0.09, const Double_t currentRight=0.13, const Double_t currentBottom=0.14);
  static void ToNaturalScale(TGraph* gr);
  static void ToNaturalScale(TGraphErrors* gr);
  static void ToNaturalScale(TGraphAsymmErrors* gr);
  static void ToNaturalScale(TAxis *ax);
  static void ToNaturalScale(TH1 *hh);
  static void AxisStyle(TAxis *ax, Bool_t kcen=0);
  static void ResetStyle(TLegend *obj, Double_t mar=-999, const Double_t ff=0.8);
  static void ResetStyle(TLatex *obj, const Double_t ff=1);
  static void ResetStyle(TPaveText *obj, const Double_t ff=1);
  static void ResetStyle(TGraph *obj);
  static void ResetStyle(TF1 * obj, Bool_t kcen=0);
  static void ResetStyle(TH1 * obj, TVirtualPad* cpad=0x0, Bool_t kcen=0);
  static void ResetStyle(THStack * obj,  Bool_t kcen=0);
  static void SetGlobalStyle(const Int_t lStat=0, Bool_t kcolor = 1);
  static void SetColor();

  static void BinLogX(TAxis *axis);
  static void UpdateLogX(TH1 *obj);
  static void UpdateLogX(TGraphAsymmErrors *obj);
  static void UpdateLogX(TGraphErrors *obj);

  static void DividegPad(int nx, int ny, double l, double r, double t, double b);
  static void DividePad(TPad * pin, int nx, int ny, double l, double r, double t, double b, TList * lout, Int_t &npad, const Int_t lastny);

  static TGraphAsymmErrors * GetInverse(const TGraphAsymmErrors * graw);
  static TH1D * Func2Hist(TF1 * ff, const Bool_t klogx);

  static Double_t fgkTextSize;
  static Double_t fgkTitleSize;
  static Double_t fgkMarkerSize;
  static Double_t fgkLineWidth;
  static Int_t fgkTextFont;
  static Double_t fgkLabelOffset;
  static Double_t fgkXTitleOffset;
  static Double_t fgkYTitleOffset;
  static Double_t fgkTickLength;

 private:
};



#endif
