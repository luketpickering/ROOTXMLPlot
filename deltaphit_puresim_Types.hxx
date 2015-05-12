#ifndef PLOTTINGTYPES_HXX_SEEN
#define PLOTTINGTYPES_HXX_SEEN

namespace PlottingTypes {

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
    TH1F* SetSeriesStyle(TH1F* series){
      series->SetLineColor(LineColor);
      return series;
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

    TH1F* SetSeriesStyle(TH1F* series){
      series->SetLineWidth(LineWidth);
      series->SetLineStyle(LineStyle);
      series->SetMarkerStyle(MarkerStyle);
      return series;
    }

  };

  struct Selection1D {
    std::string Name;
    std::string Tag;
    std::string DrawString;
    int NBins;
    double XBinLow;
    double XBinUpper;
    TCut Cut;
    bool UseLogX;
    Selection1D(){
      Name = "";
      Tag = "";
      DrawString = "";
      NBins = 0;
      XBinLow = 0;
      XBinUpper = 0;
      Cut = "";
      UseLogX = false;
    }
    Selection1D(std::string name,
                std::string tag,
                std::string drawString,
                int nBins,
                double xBinLow,
                double xBinUpper,
                TCut cut,
                bool useLogX=false){
      Name = name;
      Tag = tag;
      DrawString = drawString;
      NBins = nBins;
      XBinLow = xBinLow;
      XBinUpper = xBinUpper;
      Cut = cut;
      UseLogX = useLogX;
    }
  };

  std::string MakeHistoName(Generator const& gen, Target const& tar,
    Selection1D const& sel){
    return (gen.Name + "_" + tar.Name + "_" + sel.Tag
      + (sel.UseLogX?"_logx":""));
  }
}
#endif
