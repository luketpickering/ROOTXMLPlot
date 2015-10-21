#ifndef PLOTTINGTYPES_HXX_SEEN
#define PLOTTINGTYPES_HXX_SEEN

#include <tuple>

#include "TObject.h"
#include "TH1.h"
#include "TCut.h"

namespace PlottingTypes {

  struct ExternalFile {
    std::string Name;
    std::string Location;

    ExternalFile(std::string name, std::string location){
      Name = name;
      Location = location;
    }
  };

  struct Sample {
    std::string Name;
    std::string Flux;
    std::string FileLocation;
    std::string Description;
    bool InvalidateCache;

    Sample(std::string const name="", std::string const filelocation="",
      bool invalidatecache=true,
      std::string const flux="", std::string const description=""){
      Name = name;
      FileLocation = filelocation;
      Flux = flux;
      Description = description;
      InvalidateCache = invalidatecache;
    }

    bool IsValid() const {
      return ((Name.length() > 0) && (FileLocation.length() > 0));
    }

    bool operator==(Sample const &other) const {
      return ((Name==other.Name)&&(FileLocation==other.FileLocation));
    }

  };

  struct Generator {
    std::string Name;
    std::vector<Sample> Samples;

    Generator(char const* name=""){
      Name = name;
    }

    bool operator==(Generator const &other) const {
      return (Name==other.Name);
    }
  };

  struct Selection {
    std::string Name;
    std::string DrawString;
    TCut Cut;
    bool InvalidateCache;

    Selection(){
      Name = "";
      DrawString = "";
      Cut = "";
      InvalidateCache = true;
    }
    Selection(std::string const &name,
              std::string const &drawString,
              TCut const &cut,
              bool const &invalidatecache=true){
      Name = name;
      DrawString = drawString;
      Cut = cut;
      InvalidateCache = invalidatecache;
    }
  };

  struct Selection1D : public Selection {

    int NBins;
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

    Selection1D(std::string name,
                std::string drawString,
                TCut const &cut,
                int nBins,
                double xBinLow,
                double xBinUpper,
                bool invalidatecache=true,
                bool doPerUseXOffset=false,
                double perUseXOffset=0.005) :
      Selection(name, drawString, cut, invalidatecache) {

      NBins = nBins;
      XBinLow = xBinLow;
      XBinUpper = xBinUpper;
      DoPerUseXOffset = doPerUseXOffset;
      PerUseXOffset = perUseXOffset;
    }
  };

  struct Selection2D : public Selection {

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

    Selection2D(std::string name,
                std::string drawString,
                TCut const &cut,
                int nxBins,
                double xBinLow,
                double xBinUpper,
                int nyBins,
                double yBinLow,
                double yBinUpper,
                bool invalidatecache=true) :
      Selection(name, drawString, cut, invalidatecache) {

      NXBins = nxBins;
      XBinLow = xBinLow;
      XBinUpper = xBinUpper;

      NYBins = nyBins;
      YBinLow = yBinLow;
      YBinUpper = yBinUpper;
    }
  };

  struct SeriesDescriptor {

    enum SeriesType {
      kUnSet,
      kTH1,
      kTH1RatioNumer,
      kTH1RatioDenom,
      kTH2Profile,
      kNormTH2,
      kTH2,
      kTH1QuadRatio1,
      kTH1QuadRatio2,
      kTH1QuadRatio3,
      kTH1QuadRatio4,
      kTH1External
    };

    std::string Generator;
    std::string Sample;
    std::string Selection;
    std::string LegendTitle;
    std::string DrawOpt;
    EColor LineColor;
    Int_t LineWidth;
    Int_t LineStyle;
    EColor MarkerColor;
    Int_t MarkerSize;
    Int_t MarkerStyle;
    Int_t FillStyle;
    EColor FillColor;
    SeriesType Type;

    Float_t Scale;


    SeriesDescriptor(
      std::string Generator,
      std::string Sample,
      std::string Selection,
      std::string LegendTitle,
      EColor LineColor,
      Int_t LineWidth,
      Int_t LineStyle,
      EColor MarkerColor,
      Int_t MarkerSize,
      Int_t MarkerStyle,
      EColor FillColor,
      Int_t FillStyle,
      SeriesType Type,
      Float_t Scale,
      std::string DrawOpt="HIST"){

      this->Generator = Generator;
      this->Sample = Sample;
      this->Selection = Selection;
      this->LineColor = LineColor;
      this->LineWidth = LineWidth;
      this->LineStyle = LineStyle;
      this->MarkerColor = MarkerColor;
      this->MarkerSize = MarkerSize;
      this->MarkerStyle = MarkerStyle;
      this->FillColor = FillColor;
      this->FillStyle = FillStyle;
      this->LegendTitle = LegendTitle;
      this->Type = Type;
      this->DrawOpt = DrawOpt;
      this->Scale = Scale;
    }
  };

  struct LegendDescriptor {
    Int_t NColumns;
    Float_t X1,X2,Y1,Y2;
    std::string Title;
    LegendDescriptor(){
      NColumns = 0;
      X1 = 0xdeadbeef;
      X2 = 0;
      Y1 = 0;
      Y2 = 0;
      Title = "";
    }
    LegendDescriptor(
      Int_t NColumns,
      Float_t X1, Float_t X2, Float_t Y1, Float_t Y2,
      std::string Title){
      this->NColumns = NColumns;
      this->X1 = X1;
      this->X2 = X2;
      this->Y1 = Y1;
      this->Y2 = Y2;
      this->Title = Title;
    }
  };

  struct LineDescriptor {
    Float_t X1,X2,Y1,Y2;
    Float_t X3Arrow;
    Int_t Style, Width;
    EColor Color;

    LineDescriptor(){
      X1 = 0xdeadbeef;
      X2 = 0;
      Y1 = 0;
      Y2 = 0;
      X3Arrow = 0xdeadbeef;
      Style = 0;
      Width = 0;
      Color = kWhite;
    }
    LineDescriptor(Float_t X1, Float_t X2, Float_t Y1, Float_t Y2,
      Float_t X3Arrow, Int_t Style, Int_t Width, EColor Color){
      this->X1 = X1;
      this->X2 = X2;
      this->Y1 = Y1;
      this->Y2 = Y2;
      this->X3Arrow = X3Arrow;
      this->Style = Style;
      this->Width = Width;
      this->Color = Color;
    }
    std::tuple<Float_t, Float_t, Float_t, Float_t> GetArrowCoords() const {
      Float_t X_Mid = (X1+X2)/2.0;
      Float_t Y_Mid = (Y1+Y2)/2.0;
      return std::make_tuple(X_Mid,Y_Mid, X3Arrow,Y_Mid);
    }
  };

  struct PlotDescriptor {
    std::string Name;
    std::string Title;
    std::string XAxisTitle;
    std::string YAxisTitle;
    std::string PrintName;
    bool DoPDF;
    bool AddNormToTitle;
    Int_t Logs;
    bool HaveXYMargins;
    Float_t X1,X2,Y1,Y2;
    Float_t YStretch;

    PlotDescriptor(std::string const &name,
      std::string const &title,
      std::string const &xAxisTitle,
      std::string const &yAxisTitle,
      std::string const &printName,
      Int_t logs = 0, bool doPDF=false,
      bool addNormToTitle=false,
      Float_t X1=0xdeadbeef, Float_t X2=0xdeadbeef,
      Float_t Y1=0xdeadbeef, Float_t Y2=0xdeadbeef,
      Float_t yStretch=0xDEADDEAD){

      Name = name;
      Title = title;
      PrintName = printName;
      XAxisTitle = xAxisTitle;
      YAxisTitle = yAxisTitle;
      Logs = logs;
      DoPDF = doPDF;
      AddNormToTitle = addNormToTitle;
      YStretch = yStretch;

      HaveXYMargins = false;
      if(X1 != 0xdeadbeef){
        this->X1 = X1;
        this->X2 = X2;
        this->Y1 = Y1;
        this->Y2 = Y2;
        HaveXYMargins = true;
      }
    }

    std::vector<SeriesDescriptor> Series;
    std::vector<LineDescriptor> Lines;
    LegendDescriptor Legend;
  };
}

#endif
