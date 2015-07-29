#ifndef PLOTTINGTYPES_HXX_SEEN
#define PLOTTINGTYPES_HXX_SEEN

#include "TObject.h"
#include "TH1.h"
#include "TCut.h"

namespace PlottingTypes {

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
      kTH2
    };

    std::string Generator;
    std::string Sample;
    std::string Selection;
    std::string LegendTitle;
    Int_t LineColor;
    Int_t LineWidth;
    Int_t LineStyle;
    SeriesType Type;


    SeriesDescriptor(
      std::string Generator,
      std::string Sample,
      std::string Selection,
      std::string LegendTitle,
      Int_t LineColor,
      Int_t LineWidth,
      Int_t LineStyle,
      SeriesType Type){

      this->Generator = Generator;
      this->Sample = Sample;
      this->Selection = Selection;
      this->LineColor = LineColor;
      this->LineWidth = LineWidth;
      this->LineStyle = LineStyle;
      this->LegendTitle = LegendTitle;
      this->Type = Type;
    }
  };


  struct PlotDescriptor {
    std::string Name;
    std::string Title;
    std::string XAxisTitle;
    std::string YAxisTitle;
    std::string PrintName;
    Int_t Logs;

    PlotDescriptor(std::string const &name, std::string const &title,
      std::string const &xAxisTitle, std::string const &yAxisTitle,
      std::string const &printName,
      Int_t logs = 0){
      Name = name;
      Title = title;
      PrintName = printName;
      XAxisTitle = xAxisTitle;
      YAxisTitle = yAxisTitle;
      this->Logs = logs;
    }
    std::vector<SeriesDescriptor> Series;

    ///11 top left, 10 top right, 01 bottom left, 00 bottom right
    Int_t LegendLocation;
  };
}

#endif
