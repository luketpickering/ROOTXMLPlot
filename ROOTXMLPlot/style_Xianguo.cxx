#include <math.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
using namespace std;

#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Math/Minimizer.h"

//#include "TASImage.h"
#include "TAxis.h"
#include "TColor.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TDatabasePDG.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include "TDirectory.h"
#include "TEventList.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TGeoManager.h"
#include "TGeoGlobalMagField.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TGraphPolar.h"
#include "TGrid.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLinearFitter.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TPaletteAxis.h"
#include "TPaveText.h"
#include "TPolyMarker.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TUUID.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TVirtualPad.h"

#include "style_Xianguo.hxx"

#define EPSILON 1e-12

Double_t style::fgkTextSize = 0.05;
Double_t style::fgkTitleSize = 0.05;
Double_t style::fgkMarkerSize = 1;
Double_t style::fgkLineWidth = 2;
Int_t style::fgkTextFont = 42;
Double_t style::fgkLabelOffset = 0.01;
Double_t style::fgkXTitleOffset = 1.25;//1.1;//1.25;
Double_t style::fgkYTitleOffset = 1.1;//1.2;
Double_t style::fgkTickLength = 0.02;

TGraphAsymmErrors* style::ScaleGraph(const TGraphAsymmErrors * gin, const  TGraphAsymmErrors * gref)
{
  TGraphAsymmErrors *gout = new TGraphAsymmErrors;
  //for each ref, search in, save only when both exist
  for(Int_t iref=0; iref<gref->GetN(); iref++){
    const Double_t xx = gref->GetX()[iref];

    Int_t targetid = -999;
    for(Int_t iin = 0; iin<gin->GetN(); iin++){
      if( fabs(gin->GetX()[iin]-xx)<EPSILON){
        targetid = iin;
        break;
      }
    }
    if(targetid<0){
      continue;
    }

    const Double_t yref = gref->GetY()[iref];

    //only for yref > 0
    if(yref<EPSILON){
      continue;
    }

    const Double_t xin = gin->GetX()[targetid];
    const Double_t exl = gin->GetEXlow()[targetid];
    const Double_t exh = gin->GetEXhigh()[targetid];

    const Double_t yin = gin->GetY()[targetid];
    const Double_t eyl = gin->GetEYlow()[targetid];
    const Double_t eyh = gin->GetEYhigh()[targetid];

    const Int_t iout = gout->GetN();

    gout->SetPoint(iout, xin, yin/yref);
    gout->SetPointError(iout, exl, exh, eyl/yref, eyh/yref);
  }

  gout->SetName(Form("%sDivBy%s", gin->GetName(), gref->GetName()));
  return gout;
}

void style::GraphMinMaxXY(const TGraphAsymmErrors * gr, Double_t & xmin, Double_t & xmax, Double_t &ymin, Double_t &ymax )
{
  for(Int_t ii=0; ii<gr->GetN(); ii++){
    const Double_t yy = gr->GetY()[ii];
    const Double_t xx = gr->GetX()[ii];
    const Double_t ey = gr->GetEYlow()[ii]*2;
    const Double_t ex = gr->GetEXlow()[ii]*2;

    if(yy-ey<ymin)
      ymin = yy-ey;

    if(yy+ey>ymax)
      ymax = yy+ey;

    if(xx-ex<xmin)
      xmin = xx-ex;

    if(xx+ex>xmax)
      xmax = xx+ex;
  }
}

void style::PadSetup(TPad *currentPad, const Double_t currentLeft, const Double_t currentTop, const Double_t currentRight, const Double_t currentBottom)
{
  currentPad->SetTicks(1,1);
  currentPad->SetLeftMargin(currentLeft);
  currentPad->SetTopMargin(currentTop);
  currentPad->SetRightMargin(currentRight);
  currentPad->SetBottomMargin(currentBottom);

  currentPad->SetFillColor(0);//this is the desired one!!!
}

void style::ToNaturalScale(TGraph* gr)
{
  for(Int_t ip=0; ip<gr->GetN(); ip++){
    Double_t xx=-999, yy=-999;
    gr->GetPoint(ip,xx,yy);
    gr->SetPoint(ip,TMath::Power(10,xx), yy);
  }
}

void style::ToNaturalScale(TGraphErrors* gr)
{
  for(Int_t ip=0; ip<gr->GetN(); ip++){
    Double_t xx=-999, yy=-999;
    gr->GetPoint(ip,xx,yy);
    //printf("%d %e %e\n",ip,xx,yy);
    gr->SetPoint(ip,TMath::Power(10,xx), yy);
  }
}

void style::ToNaturalScale(TGraphAsymmErrors* gr)
{
  for(Int_t ip=0; ip<gr->GetN(); ip++){
    const Double_t xx = gr->GetX()[ip];
    const Double_t exl = gr->GetEXlow()[ip];
    const Double_t exh = gr->GetEXhigh()[ip];

    const Double_t yy = gr->GetY()[ip];
    const Double_t eyl = gr->GetEYlow()[ip];
    const Double_t eyh = gr->GetEYhigh()[ip];

    //printf("%d %e %e\n",ip,xx,yy);
    gr->SetPoint(ip,TMath::Power(10,xx), yy);
    gr->SetPointError(ip, TMath::Power(10, xx) - TMath::Power(10, xx-exl), TMath::Power(10, xx+exh)-TMath::Power(10, xx), eyl, eyh);
  }
}

void style::ToNaturalScale(TAxis *ax)
{
  TAxis* oldx = (TAxis*)ax->Clone("oldx");
  ax->SetLimits(TMath::Power(10,oldx->GetXmin()), TMath::Power(10,oldx->GetXmax()));
  const Int_t nb = oldx->GetNbins();
  Double_t *bins = new Double_t[nb+1];
  bins[0]=TMath::Power(10,oldx->GetXmin());
  for(Int_t ii=1; ii<=nb; ii++){
    bins[ii]=TMath::Power(10,oldx->GetBinUpEdge(ii));
  }
  ax->Set(nb, bins);

  delete oldx;
  delete bins;
}

void style::ToNaturalScale(TH1 *hh)
{
  ToNaturalScale(hh->GetXaxis());
}


void style::AxisStyle(TAxis *ax, Bool_t kcen)
{
  ax->SetTickLength(fgkTickLength);

  ax->SetLabelFont(fgkTextFont);
  ax->SetLabelSize(fgkTextSize);
  ax->SetLabelOffset(fgkLabelOffset);

  ax->SetTitleFont(fgkTextFont);
  //printf("title %f text %f \n", fgkTitleSize, fgkTextSize);
  ax->SetTitleSize(fgkTitleSize);

  ax->CenterTitle(kcen);
}


void style::ResetStyle(TLegend *obj, Double_t mar, const Double_t ff)
{
  if(mar>0){
    obj->SetMargin(mar);
  }
  obj->SetFillStyle(-1);
  obj->SetBorderSize(-1);
  obj->SetTextFont(fgkTextFont);
  obj->SetTextSize(fgkTextSize*ff);
}

void style::ResetStyle(TLatex *obj, const Double_t ff)
{
  obj->SetNDC();
  obj->SetTextFont(fgkTextFont);
  obj->SetTextSize(fgkTextSize*ff);
}

void style::ResetStyle(TPaveText *obj, const Double_t ff)
{
  obj->SetFillStyle(-1);
  obj->SetBorderSize(-1);
  obj->SetTextFont(fgkTextFont);
  obj->SetTextSize(fgkTextSize*ff);
  obj->SetTextAlign(11);
}

void style::ResetStyle(TGraph *obj)
{
  obj->SetMarkerSize(fgkMarkerSize);
  obj->SetLineWidth(fgkLineWidth);

  AxisStyle(obj->GetXaxis());
  AxisStyle(obj->GetYaxis());

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
}

void style::ResetStyle(TF1 * obj, Bool_t kcen)
{
  obj->SetMarkerSize(fgkMarkerSize);

  AxisStyle(obj->GetXaxis(), kcen);
  AxisStyle(obj->GetYaxis(), kcen);

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
}


void style::ResetStyle(TH1 * obj, TVirtualPad* cpad, Bool_t kcen)
{
  obj->SetMarkerSize(fgkMarkerSize);

  AxisStyle(obj->GetXaxis(), kcen);
  AxisStyle(obj->GetYaxis(), kcen);

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);

  if(cpad){
    TPaletteAxis *palette = (TPaletteAxis*)obj->GetListOfFunctions()->FindObject("palette");
    if(!palette){
      printf("ResetStyle no palette!!\n");
      obj->GetListOfFunctions()->Print();
    }
    else{
      palette->SetX1NDC(1-cpad->GetRightMargin()+0.005);
      palette->SetX2NDC(1-cpad->GetRightMargin()/3*2);
      palette->SetY1NDC(cpad->GetBottomMargin());
      palette->SetY2NDC(1-cpad->GetTopMargin());
      palette->SetLabelFont(fgkTextFont);
      palette->SetLabelSize(fgkTextSize);
      palette->SetLabelOffset(fgkLabelOffset);
    }
  }

}


void style::ResetStyle(THStack * obj, Bool_t kcen)
{
  AxisStyle(obj->GetXaxis(), kcen);
  AxisStyle(obj->GetYaxis(), kcen);

  obj->GetXaxis()->SetTitleOffset(fgkXTitleOffset);
  obj->GetYaxis()->SetTitleOffset(fgkYTitleOffset);
}


void style::SetGlobalStyle(const Int_t lStat, Bool_t kcolor)
{
  // Set gStyle
  // From plain

  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(-1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(-1);
  gStyle->SetLegendBorderSize(-1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(fgkTextFont);
  gStyle->SetStatFont(fgkTextFont);
  gStyle->SetStatFontSize(fgkTextSize);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(fgkTickLength,"xy");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(fgkTextSize,"xyz");
  gStyle->SetLabelFont(fgkTextFont,"xyz");
  gStyle->SetLabelOffset(fgkLabelOffset,"xyz");
  gStyle->SetTitleFont(fgkTextFont,"xyz");
  gStyle->SetTitleFont(fgkTextFont,"");
  gStyle->SetTitleFontSize(fgkTitleSize);
  gStyle->SetTitleOffset(fgkXTitleOffset,"x");
  gStyle->SetTitleOffset(fgkYTitleOffset,"y");
  gStyle->SetTitleOffset(1.0,"z");
  gStyle->SetTitleSize(fgkTitleSize,"xyz");
  gStyle->SetTitleSize(fgkTitleSize,"");
  gStyle->SetMarkerSize(fgkMarkerSize);
  gStyle->SetPalette(1,0);
  if (lStat){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }

  TGaxis::SetMaxDigits(3);
  gStyle->SetTitleBorderSize(-1);

  if(kcolor){
    SetColor();
  }

  gROOT->ForceStyle();
}


void style::SetColor()
{
  gStyle->SetHistFillColor(0);
  //gStyle->SetFillColor(0);//it conflicts with color palette
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetTitleFillColor(0);

  //---

  const Int_t nRGBs = 5;
  const Int_t nCont = 100;
  Double_t stops[nRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[nRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[nRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[nRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  const Int_t cgc = TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nCont);

  const Int_t nc = nCont;
  Int_t colors[nc];
  gStyle->SetNumberContours(nc);
  for(Int_t ii=0; ii<nc; ii++){
    colors[ii]=cgc+ii;
  }
  gStyle->SetPalette(nc, colors);
}

void style::BinLogX(TAxis* axis)
{
  //void XGLUtils::BinLogX(TAxis *axis)
  //
  // Method for the correct logarithmic binning of histograms
  // copied and modified from AliTPCcalibBase

  const Int_t bins = axis->GetNbins();

  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  if (from<EPSILON) return;
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  const Double_t factor = TMath::Power(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;

}

void style::UpdateLogX(TH1 *obj)
{
  TAxis * axis = obj->GetXaxis();
  double xmin = axis->GetXmin();
  double xmax = axis->GetXmax();
  axis->SetLimits(TMath::Power(10,xmin), TMath::Power(10,xmax));
  BinLogX(axis);
}

void style::UpdateLogX(TGraphAsymmErrors *obj)
{
  const Int_t np = obj->GetN();
  for(Int_t ii=0; ii<np; ii++){
    Double_t xx =-999, yy=-999;
    obj->GetPoint(ii,xx,yy);
    Double_t ex0 = obj->GetErrorXlow(ii);
    Double_t ex1 = obj->GetErrorXhigh(ii);
    const Double_t ey0 = obj->GetErrorYlow(ii);
    const Double_t ey1 = obj->GetErrorYhigh(ii);

    ex0 = TMath::Power(10, xx) - TMath::Power(10, xx-ex0);
    ex1 = TMath::Power(10, xx+ex1) - TMath::Power(10, xx);
    xx = TMath::Power(10, xx);

    obj->SetPoint(ii, xx, yy);
    obj->SetPointError(ii, ex0, ex1, ey0, ey1);
  }
}


void style::UpdateLogX(TGraphErrors *obj)
{
  const Int_t np = obj->GetN();
  for(Int_t ii=0; ii<np; ii++){
    Double_t xx =-999, yy=-999;
    obj->GetPoint(ii,xx,yy);
    Double_t ex = obj->GetErrorX(ii);
    const Double_t ey = obj->GetErrorY(ii);

    ex = TMath::Power(10, xx) - TMath::Power(10, xx-ex);
    xx = TMath::Power(10, xx);

    obj->SetPoint(ii, xx, yy);
    obj->SetPointError(ii, ex, ey);
  }
}


void style::DividegPad(int nx, int ny, double l, double r, double t, double b)
{
  gStyle->SetPadTopMargin(0.);
  gStyle->SetPadBottomMargin(0.);
  gStyle->SetPadLeftMargin(0.);
  gStyle->SetPadRightMargin(0.);

   int ix, iy, n=0;
   double x1, x2, y1, y2;
   double dx = ((1-r)*(1-l))/((1-r)*(1-l)*(nx-2)-r+2-l);
   double dl = dx/(1-l);
   double dy = ((1-t)*(1-b))/((1-t)*(1-b)*(ny-2)-t+2-b);
   double db = dy/(1-b);
   char *name  = new char [strlen(gPad->GetName())+6];

   y1 = 0;
   y2 = db;
   for (iy=0; iy<ny; iy++) {
      x1 = 0;
      x2 = dl;
      for (ix=0;ix<nx;ix++) {
         if (x1 > x2) continue;
         n++;
         sprintf(name,"%s_%d",gPad->GetName(),n);
         TPad *pad = new TPad(name,name,x1,y1,x2,y2,0);
         if (ix==0)    pad->SetLeftMargin(l);
         if (ix==nx-1) pad->SetRightMargin(r);
         if (iy==ny-1) pad->SetTopMargin(t);
         if (iy==0)    pad->SetBottomMargin(b);
         x1 = x2;
         if (ix==nx-2) x2 = 1;
         else          x2 = x1+dx;
         pad->SetNumber(n);
         pad->Draw();
         pad->SetTicks();
      }
      y1 = y2;
      if (iy==ny-2) y2 = 1;
      else          y2 = y1+dy;
   }
}


void style::DividePad(TPad *pin, int nx, int ny, double l, double r, double t, double b, TList * lout, Int_t &npad, const Int_t lastny)
{
  pin->SetTopMargin(0.);
  pin->SetBottomMargin(0.);
  pin->SetLeftMargin(0.);
  pin->SetRightMargin(0.);

   int ix, iy, n=0;
   double x1, x2, y1, y2;
   double dx = nx>=2? ( ((1-r)*(1-l))/((1-r)*(1-l)*(nx-2)-r+2-l) ) : 0;
   double dl = nx>=2? dx/(1-l) : 1;
   double dy = ((1-t)*(1-b))/((1-t)*(1-b)*(ny-2)-t+2-b);
   double db = dy/(1-b);
   char *name  = new char [strlen(pin->GetName())+6];

   //printf("test dx %f dl %f dy %f db %f\n", dx, dl, dy, db);

   y1 = 0;
   y2 = db;
   for (iy=0; iy<lastny+1; iy++) {
     //printf("test1 iy %d\n", iy);
      x1 = 0;
      x2 = dl;
      for (ix=0;ix<nx;ix++) {
         if (x1 > x2) continue;
         //printf("test2 ix %d\n", ix);

         n++;
         sprintf(name,"%s_%d",pin->GetName(),n);
         TPad *pad = new TPad(name,name,x1,y1,x2,y2,0);lout->Add(pad);
         PadSetup(pad);

  pad->SetTopMargin(0.);
  pad->SetBottomMargin(0.);
  pad->SetLeftMargin(0.);
  pad->SetRightMargin(0.);

         if (ix==0)    pad->SetLeftMargin(l);
         if (ix==nx-1) pad->SetRightMargin(r);
         if (iy==lastny) pad->SetTopMargin(t);
         if (iy==0)    pad->SetBottomMargin(b);
         x1 = x2;
         if (ix==nx-2) x2 = 1;
         else          x2 = x1+dx;
         pad->SetNumber(npad);npad++;
         pad->Draw();
         pad->SetTicks();
      }
      y1 = y2;
      if (iy==lastny-1) y2 = 1;
      else          y2 = y1+dy;
   }
}

TGraphAsymmErrors * style::GetInverse(const TGraphAsymmErrors * graw)
{
  TGraphAsymmErrors * ginv = (TGraphAsymmErrors*) graw->Clone(Form("%sInv", graw->GetName()));

  const Int_t np = graw->GetN();

  Int_t ipoint = 0;
  for(Int_t ii=0; ii<np; ii++){
    const Double_t xx = graw->GetX()[ii];
    const Double_t exl = graw->GetEXlow()[ii];
    const Double_t exh = graw->GetEXhigh()[ii];

    const Double_t tmpyy = graw->GetY()[ii];
    const Double_t tmpeyl = graw->GetEYlow()[ii];
    const Double_t tmpeyh = graw->GetEYhigh()[ii];

    if(tmpyy<EPSILON)
      continue;

    const Double_t yy = 1/tmpyy;
    const Double_t eyl = yy * tmpeyl/tmpyy;
    const Double_t eyh = yy * tmpeyh/tmpyy;

    ginv->SetPoint(ipoint, xx, yy);
    ginv->SetPointError(ipoint, exl, exh, eyl, eyh);
    ipoint++;
  }

  ginv->Set(ipoint);

  //have to reset title, at least, after Set
  ginv->GetYaxis()->SetTitle(Form("inverse of %s",graw->GetYaxis()->GetTitle()));
  ginv->GetXaxis()->SetTitle(graw->GetXaxis()->GetTitle());

  return ginv;
}


TH1D * style::Func2Hist(TF1 * ff, const Bool_t klogx)
{
  TH1D * hh=0x0;
  TString hname(Form("hist%s%d", ff->GetName(), gRandom->Integer(1000)));
  if(gDirectory->Get(hname)){
    hname = Form("hist%s%d", ff->GetName(), gRandom->Integer(1000));
  }

  hh=new TH1D(hname,"",500, ff->GetXmin(), ff->GetXmax());
  if(klogx){
    BinLogX(hh->GetXaxis());
  }

  hh->SetLineColor(ff->GetLineColor());
  hh->SetLineStyle(ff->GetLineStyle());

  for(Int_t ii=1; ii<=hh->GetNbinsX(); ii++){
    const Double_t xi = hh->GetBinCenter(ii);
    const Double_t yi = ff->Eval(xi);
    hh->SetBinContent(ii, yi);
  }

  hh->GetXaxis()->SetTitle(ff->GetXaxis()->GetTitle());
  hh->GetYaxis()->SetTitle(ff->GetYaxis()->GetTitle());
  return hh;
}
