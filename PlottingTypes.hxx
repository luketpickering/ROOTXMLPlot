#ifndef PLOTTINGTYPES_HXX_SEEN
#define PLOTTINGTYPES_HXX_SEEN

#include "TCut.h"

namespace PlottingTypes {

  struct CutString {
    private:
      size_t len(){return cutstr.length();}
    public:
    std::string cutstr;

    CutString(){cutstr="";}
    CutString(std::string const &str){cutstr = "(" + str + ")";}
    CutString(char const * str){cutstr = str; cutstr = "(" + cutstr + ")";}
    CutString(CutString const & other){
      cutstr = other.cutstr;
    }
    CutString operator &=(CutString const & other){
      if(!len()){
        return CutString(other);
      }
      cutstr = "(" + cutstr + "&&" + other.cutstr + ")";
      return (*this);
    }
    CutString operator &&(CutString const & other){
      return (*this) &= other;
    }
    CutString operator |=(CutString const & other){
      if(!len()){
        return CutString(other);
      }
      cutstr = "(" + cutstr + "||" + other.cutstr + ")";
      return (*this);
    }
    CutString operator ||(CutString const & other){
       return (*this) |= other;
    }
    CutString operator !(){
      if(!len()){
        return CutString();
      }
      cutstr = "(!" + cutstr + ")";
      return (*this);
    }
  };

  CutString operator &&(std::string const &str, CutString const & Cut){
    return CutString(str) && Cut;
  }
  CutString operator ||(std::string const &str, CutString const & Cut){
    return CutString(str) || Cut;
  }

  struct Target {
    std::string Name;
    std::string FileLocation;
    int LineColor;

    Target(char const* name="", int linecolor=kBlack,
      char const* filelocation=""){
      Name = name;
      LineColor = linecolor;
      FileLocation = filelocation;
    }
    TH1* SetSeriesStyle(TH1* series) const {
      series->SetLineColor(LineColor);
      return series;
    }
    bool operator==(Target const &other) const {
      return ((Name==other.Name)&&
              (LineColor == other.LineColor));
    }

  };

  struct Generator {
    std::string Name;
    int LineStyle;
    int LineWidth;
    int MarkerStyle;

    std::vector<Target> Targets;

    Generator(char const* name="", int linestyle=1,
      int linewidth=2, int markerstyle=0){
      Name = name;
      LineStyle = linestyle;
      LineWidth = linewidth;
      MarkerStyle = markerstyle;
    }

    TH1* SetSeriesStyle(TH1* series) const {
      series->SetLineWidth(LineWidth);
      series->SetLineStyle(LineStyle);
      series->SetMarkerStyle(MarkerStyle);
      return series;
    }

    bool operator==(Generator const &other) const {
      return ((Name==other.Name)&&
              (LineStyle == other.LineStyle)&&
              (LineWidth == other.LineWidth)&&
              (MarkerStyle == other.MarkerStyle)&&
              (Targets.size() == other.Targets.size()));
    }
  };

  struct Selection {
    std::string Name;
    std::string Tag;
    std::string DrawString;
    TCut Cut;
    bool UseLogX;
  };

  struct Selection1D : public Selection {
    int NBins;
    double XBinLow;
    double XBinUpper;
    bool DoPerUseXOffset;
    Selection1D(){
      Name = "";
      Tag = "";
      DrawString = "";
      NBins = 0;
      XBinLow = 0;
      XBinUpper = 0;
      Cut = "";
      UseLogX = false;
      DoPerUseXOffset = false;
    }
    Selection1D(std::string name,
                std::string tag,
                std::string drawString,
                int nBins,
                double xBinLow,
                double xBinUpper,
                TCut const &cut,
                bool useLogX=false,
                bool doPerUseXOffset=false){
      Name = name;
      Tag = tag;
      DrawString = drawString;
      NBins = nBins;
      XBinLow = xBinLow;
      XBinUpper = xBinUpper;
      Cut = cut;
      UseLogX = useLogX;
      DoPerUseXOffset = doPerUseXOffset;
    }
  };
  struct Selection2D : public Selection {
    int NXBins;
    double XBinLow;
    double XBinUpper;
    int NYBins;
    double YBinLow;
    double YBinUpper;
    bool UseLogY;
    Selection2D(){
      Name = "";
      Tag = "";
      DrawString = "";
      NXBins = 0;
      XBinLow = 0;
      XBinUpper = 0;
      NYBins = 0;
      YBinLow = 0;
      YBinUpper = 0;
      Cut = "";
      UseLogX = false;
      UseLogY = false;
    }
    Selection2D(std::string name,
                std::string tag,
                std::string drawString,
                int nxBins,
                double xBinLow,
                double xBinUpper,
                int nyBins,
                double yBinLow,
                double yBinUpper,
                TCut const &cut,
                bool useLogX=false,
                bool useLogY=false){
      Name = name;
      Tag = tag;
      DrawString = drawString;
      NXBins = nxBins;
      XBinLow = xBinLow;
      XBinUpper = xBinUpper;
      NYBins = nyBins;
      YBinLow = yBinLow;
      YBinUpper = yBinUpper;
      Cut = cut;
      UseLogX = useLogX;
      UseLogY = useLogY;
    }
  };

  std::string MakeHistoName(Generator const& gen, Target const& tar,
    Selection const& sel){
    return (gen.Name + "_" + tar.Name + "_" + sel.Tag
      + (sel.UseLogX?"_logx":""));
  }
}
#endif
