#include <iterator>
#include <tuple>

#include "TH1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TArrow.h"
#include "TProfile.h"
#include "TGaxis.h"

#include "PureGenUtils.hxx"

#include "PlottingAutoPlots.hxx"
#include "PlottingSelections.hxx"
#include "style_Xianguo.hxx"
#include "PlottingIO.hxx"
#include "NeutrinoTools_Xianguo.hxx"

using namespace PlottingTypes;
using namespace PlottingDefaults;
using namespace PlottingUtils;

namespace AutoPlots {
std::vector<PlotDescriptor> AutoPlots;

std::tuple<Float_t,Float_t,Float_t,Float_t>
  GetX1X2Y1Y2Attributes(TXMLEngine &xmlengine,
  XMLDocPointer_t node){

  Float_t X1,X2,Y1,Y2;

  std::string const &nodeX1 =
    IOUtils::GetAttrValue(xmlengine, node, "X1");
  std::string const &nodeX2 =
    IOUtils::GetAttrValue(xmlengine, node, "X2");
  std::string const &nodeY1 =
    IOUtils::GetAttrValue(xmlengine, node, "Y1");
  std::string const &nodeY2 =
    IOUtils::GetAttrValue(xmlengine, node, "Y2");

  try{
    X1 = std::stof(nodeX1);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[ERROR]: Couldn't parse X1: \"" << nodeX1 << "\"\n\t"
      << ia.what() << std::endl;
    X1 = 0xdeadbeef;
  }
  try{
    X2 = std::stof(nodeX2);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[ERROR]: Couldn't parse X2: \"" << nodeX2 << "\"\n\t"
      << ia.what() << std::endl;
    X2 = 0xdeadbeef;
  }

  try{
    Y1 = std::stof(nodeY1);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[ERROR]: Couldn't parse Y1: \"" << nodeY1 << "\"\n\t"
      << ia.what() << std::endl;
    Y1 = 0xdeadbeef;
  }

  try{
    Y2 = std::stof(nodeY2);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[ERROR]: Couldn't parse Y2: \"" << nodeY2 << "\"\n\t"
      << ia.what() << std::endl;
    Y2 = 0xdeadbeef;
  }
  std::cout << "[INFO]: Read " << X1 << ", " << X2 << ", " << Y1 << ", " << Y2
    << " from: "  << nodeX1 << ", " << nodeX2 << ", " << nodeY1 << ", "
    << nodeY2 << std::endl;
  return std::make_tuple(X1,X2,Y1,Y2);

}

std::tuple<Int_t,Int_t,EColor>
  GetROOTLineAttributes(TXMLEngine &xmlengine,
  XMLDocPointer_t node){

  Int_t Style = 0, Width = 0;
  EColor Color = kBlack;

  std::string const &nodeLineStyle =
    IOUtils::GetAttrValue(xmlengine, node, "LineStyle");
  std::string const &nodeLineColor =
    IOUtils::GetAttrValue(xmlengine, node, "LineColor");
  std::string const &nodeLineWidth =
    IOUtils::GetAttrValue(xmlengine, node, "LineWidth");

  if(nodeLineStyle.length()){
    try{
      Style = std::stoi(nodeLineStyle);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << nodeLineStyle << "\"\n\t"
        << ia.what() << std::endl;
    }
  }
  if(nodeLineColor.length()){
    try{
      Color = IOUtils::ParseRootColor(nodeLineColor);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << nodeLineColor << "\"\n\t"
        << ia.what() << std::endl;
    }
  }

  if(nodeLineWidth.length()){
    try{
      Width = std::stoi(nodeLineWidth);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << nodeLineWidth << "\"\n\t"
        << ia.what() << std::endl;
    }
  }

  return std::make_tuple(Style,Width,Color);

}

std::tuple<Int_t,Int_t,EColor,Int_t,Int_t,EColor,Int_t,EColor>
  GetROOTTH1StyleAttributes(TXMLEngine &xmlengine,
  XMLDocPointer_t node){

  Int_t LineStyle = 0, LineWidth = 0;
  EColor LineColor = kBlack;

  std::tie(LineStyle,LineWidth,LineColor) =
    GetROOTLineAttributes(xmlengine, node);

  Int_t MarkerStyle = 0, MarkerSize = 0;
  EColor MarkerColor = kWhite;

  std::string const &nodeMarkerStyle =
    IOUtils::GetAttrValue(xmlengine, node, "MarkerStyle");
  std::string const &nodeMarkerColor =
    IOUtils::GetAttrValue(xmlengine, node, "MarkerColor");
  std::string const &nodeMarkerSize =
    IOUtils::GetAttrValue(xmlengine, node, "MarkerSize");

  if(nodeMarkerStyle.length()){
    try{
      MarkerStyle = std::stoi(nodeMarkerStyle);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << nodeMarkerStyle << "\"\n\t"
        << ia.what() << std::endl;
    }
  }
  if(nodeMarkerColor.length()){
    try{
      MarkerColor = IOUtils::ParseRootColor(nodeMarkerColor);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << nodeMarkerColor << "\"\n\t"
        << ia.what() << std::endl;
    }
  }

  if(nodeMarkerSize.length()){
    try{
      MarkerSize = std::stoi(nodeMarkerSize);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << nodeMarkerSize << "\"\n\t"
        << ia.what() << std::endl;
    }
  }

  Int_t FillStyle = 0;
  EColor FillColor = kWhite;

  std::string const &seriesFillStyle =
    IOUtils::GetAttrValue(xmlengine, node, "FillStyle");
  std::string const &seriesFillColor =
    IOUtils::GetAttrValue(xmlengine, node, "FillColor");

  if(seriesFillStyle.length()){
    try{
      FillStyle = std::stoi(seriesFillStyle);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << seriesFillStyle << "\"\n\t"
        << ia.what() << std::endl;
    }
  }
  if(seriesFillColor.length()){
    try{
      FillColor = IOUtils::ParseRootColor(seriesFillColor);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: \"" << seriesFillColor << "\"\n\t"
        << ia.what() << std::endl;
    }
  }

  return std::make_tuple(LineStyle, LineWidth, LineColor, MarkerStyle,
    MarkerSize, MarkerColor, FillStyle, FillColor);

}

LineDescriptor GetLineDescriptorXML(TXMLEngine &xmlengine,
  XMLDocPointer_t lineDescriptNode){

  Float_t X1,X2,Y1,Y2;
  std::tie(X1,X2,Y1,Y2) = GetX1X2Y1Y2Attributes(xmlengine,lineDescriptNode);

  Int_t lineStyle, lineWidth;
  EColor lineColor;
  std::tie(lineStyle,lineWidth,lineColor) =
    GetROOTLineAttributes(xmlengine, lineDescriptNode);

  if( (X1 == 0xdeadbeef) ||
      (X2 == 0xdeadbeef) ||
      (Y1 == 0xdeadbeef) ||
      (Y2 == 0xdeadbeef) ){
    return LineDescriptor();
  }

  std::string const &lineX3Arrow =
    IOUtils::GetAttrValue(xmlengine, lineDescriptNode, "X3Arrow");

  Float_t X3Arrow = 0xdeadbeef;

  if(lineX3Arrow.length()){
    try{
      X3Arrow = std::stof(lineX3Arrow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse X3Arrow: \"" << lineX3Arrow << "\"\n\t"
        << ia.what() << std::endl;
      return LineDescriptor();
    }
  }
  return LineDescriptor(X1, X2, Y1, Y2, X3Arrow, lineStyle, lineWidth,
    lineColor);
}

LegendDescriptor GetLegendDescriptorXML(TXMLEngine &xmlengine,
  XMLDocPointer_t legendDescriptNode){

  Float_t X1,X2,Y1,Y2;
  std::tie(X1,X2,Y1,Y2) = GetX1X2Y1Y2Attributes(xmlengine,legendDescriptNode);

  if( (X1 == 0xdeadbeef) ||
      (X2 == 0xdeadbeef) ||
      (Y1 == 0xdeadbeef) ||
      (Y2 == 0xdeadbeef) ){
    return LegendDescriptor();
  }

  std::string const &legendNColumns =
    IOUtils::GetAttrValue(xmlengine, legendDescriptNode, "NColumns");
  std::string const &legendTitle =
    IOUtils::GetAttrValue(xmlengine, legendDescriptNode, "Title");

  Int_t NColumns = 1;

  try{
    NColumns = std::stof(legendNColumns);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[ERROR]: Couldn't parse NColumns: \"" << legendNColumns
      << "\"\n\t" << ia.what() << std::endl;
    NColumns = 1;
  }

  return LegendDescriptor(NColumns, X1, X2, Y1, Y2, legendTitle);
}

void AddCanvasDescriptionXML(TXMLEngine &xmlengine,
  XMLDocPointer_t canvDescriptNode){

  std::string const &canvName =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "Name");
  std::string const &canvTitle =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "Title");
  std::string const &canvPrintName =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "PrintName");
  std::string const &canvXAxisTitle =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "XTitle");
  std::string const &canvYAxisTitle =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "YTitle");
  std::string const &canvYStretch =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "YStretch");


  std::string const &canvlogs =
    IOUtils::GetAttrValue(xmlengine, canvDescriptNode, "Logs");
  Int_t Logs = 0;
  try{
    Logs = std::stoi(canvlogs);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[Warn]: Couldn't parse Logs: \"" << canvlogs << "\"\n\t"
      << ia.what() << std::endl;
  }

  if(!canvName.length()){
    std::cout << "[WARN]: Canvas descriptor had no name."
      << std::endl;
    return;
  }

  Float_t YStretch = 0XDEADDEAD;
  try{
    YStretch = std::stof(canvYStretch);
  } catch (const std::invalid_argument& ia) {
    std::cout << "[ERROR]: Couldn't parse YStretch: \"" << canvYStretch << "\"\n\t"
      << ia.what() << std::endl;
    YStretch = 0XDEADDEAD;
  }

  bool doPDF = PGUtils::str2bool(IOUtils::GetAttrValue(xmlengine,
    canvDescriptNode, "DoPDF"));
  bool AddNorm = PGUtils::str2bool(IOUtils::GetAttrValue(xmlengine,
    canvDescriptNode, "AddNorm"));

  Float_t X1,X2,Y1,Y2;
  std::tie(X1,X2,Y1,Y2) = GetX1X2Y1Y2Attributes(xmlengine,canvDescriptNode);

  PlotDescriptor pd(canvName, canvTitle,canvXAxisTitle, canvYAxisTitle,
    canvPrintName, Logs, doPDF, AddNorm, X1,X2,Y1,Y2, YStretch);

  bool ExpectingDenom = false;
  size_t FoundNQuadRatios = 0;


  for(XMLNodePointer_t SeriesDescriptor = xmlengine.GetChild(canvDescriptNode);
    (SeriesDescriptor != NULL);
    SeriesDescriptor = xmlengine.GetNext(SeriesDescriptor)){

    SeriesDescriptor::SeriesType type = SeriesDescriptor::kUnSet;

    if(std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Series1D"){
      type = SeriesDescriptor::kTH1;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Ratio1D"){
      type = ExpectingDenom? SeriesDescriptor::kTH1RatioDenom:SeriesDescriptor::kTH1RatioNumer;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "QuadRatio"){

      switch(FoundNQuadRatios){
        case 0:{
          type = SeriesDescriptor::kTH1QuadRatio1;
          break;
        }
        case 1:{
          type = SeriesDescriptor::kTH1QuadRatio2;
          break;
        }
        case 2:{
          type = SeriesDescriptor::kTH1QuadRatio3;
          break;
        }
        case 3:{
          type = SeriesDescriptor::kTH1QuadRatio4;
          break;
        }
        default:{
          std::cout << "[ERROR]: Didn't expect this..." << std::endl;
          break;
        }
      }
      FoundNQuadRatios++;

    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Profile2D"){
      type = SeriesDescriptor::kTH2Profile;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "External1D"){
      type = SeriesDescriptor::kTH1External;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Series2D"){
      if(pd.Series.size()){
        std::cout << "[ERROR]: I was not expecting a Series2D element unless it"
         " was the first and only." << std::endl;
        return;
      }
      if(PGUtils::str2bool(IOUtils::GetAttrValue(xmlengine, SeriesDescriptor,
          "Normalise"))){
        type = SeriesDescriptor::kNormTH2;
      } else {
        type = SeriesDescriptor::kTH2;
      }
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Line"){
      LineDescriptor ld = GetLineDescriptorXML(xmlengine,SeriesDescriptor);
      if(ld.X1 != 0xdeadbeef){
        pd.Lines.push_back(ld);
      }
      continue;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Legend"){
      LegendDescriptor legd =
        GetLegendDescriptorXML(xmlengine,SeriesDescriptor);
      if(legd.X1 != 0xdeadbeef){
        pd.Legend = legd;
      }
      continue;
    } else {
      continue;
    }

    std::string seriesGenerator;
    std::string seriesSample;
    std::string seriesSelection;
    std::string const &seriesLegendTitle =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "LegendTitle");

    if(type != SeriesDescriptor::kTH1External){

      seriesGenerator =
        IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Generator");
      seriesSample =
        IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Sample");
      seriesSelection =
        IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Selection");

      if( (!seriesGenerator.length()) || (!seriesSample.length()) ||
          (!seriesSelection.length()) ){
        std::cout << "[WARN]: series descriptor had no Generator, Sample, "
          "or Selection attached to it.." << std::endl;
        continue;
      }

      //Ratio checks
      if(ExpectingDenom){
        if(type != SeriesDescriptor::kTH1RatioDenom){
          std::cout << "[WARN]: Failed to find second <Ratio1D /> element. You "
            "need to specify a numerator series and a denominator series."
            << std::endl;
          return;
        } else {
          //Either we found the numerator or we found the denominator
          ExpectingDenom = !ExpectingDenom;
        }

        PlottingTypes::Selection1D const & NumerSel =
          PlottingSelections::Selections1D[pd.Series.back().Sample];

        PlottingTypes::Selection1D const & DenomSel =
          PlottingSelections::Selections1D[seriesSample];

        if( (NumerSel.NBins != DenomSel.NBins) ||
            (NumerSel.XBinLow != DenomSel.XBinLow) ||
            (NumerSel.XBinUpper != DenomSel.XBinUpper) ){

          std::cout << "[WARN]: Numerator [" << NumerSel.XBinLow << ".."
            << NumerSel.NBins << ".." << NumerSel.XBinUpper << "] and Denominator"
            " [" << DenomSel.XBinLow << ".."  << DenomSel.NBins << ".."
            << DenomSel.XBinUpper << "] bins are incompatible." << std::endl;
          return;
        }
      }

      if(FoundNQuadRatios){
        if( (FoundNQuadRatios != 4) &&
            (SeriesDescriptor::kTH1QuadRatio1 != type ) &&
            (SeriesDescriptor::kTH1QuadRatio2 != type ) &&
            (SeriesDescriptor::kTH1QuadRatio3 != type ) ){

            std::cout << "[ERROR]: Was expecting another QuadRatio element (have "
              << FoundNQuadRatios << ") but got " << type
              << std::endl;

            throw 0;
        } else if (FoundNQuadRatios == 4){
          FoundNQuadRatios = 0;
        }
      }
    } else { //It is an external. For the moment, hack it in.
      seriesGenerator =
        IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "ExternalFile");
      seriesSample =
        IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "TH1Name");

      if( (!seriesGenerator.length()) || (!seriesSample.length()) ){
        std::cout << "[WARN]: External series descriptor had no ExternalFile "
        "or TH1Name attached to it.." << std::endl;
        continue;
      }
    }

    Int_t MarkerStyle = 0, MarkerSize = 0;
    EColor MarkerColor = kBlack;
    Int_t LineStyle = 0, LineWidth = 0;
    EColor LineColor = kBlack;
    Int_t FillStyle = 0;
    EColor FillColor = kWhite;

    std::tie(LineStyle, LineWidth, LineColor, MarkerStyle, MarkerSize,
      MarkerColor, FillStyle, FillColor) =
    GetROOTTH1StyleAttributes(xmlengine, SeriesDescriptor);

    std::string const &seriesDrawOpts =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "DrawOption");

    std::string const &seriesScale =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Scale");

    Float_t scale = 1.0;
    try{
      scale = std::stof(seriesScale);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse scale: \"" << seriesScale << "\"\n\t"
        << ia.what() << std::endl;
      scale = 1.0;
    }

    pd.Series.emplace_back(
      seriesGenerator,
      seriesSample,
      seriesSelection,
      seriesLegendTitle,
      LineColor,
      LineWidth,
      LineStyle,
      MarkerColor,
      MarkerSize,
      MarkerStyle,
      FillColor,
      FillStyle,
      type,
      scale,
      seriesDrawOpts.length()?seriesDrawOpts:"HIST");

    //Only currently handle TH2s on their own.
    if(type == SeriesDescriptor::kTH2){
      break;
    }
  }

  if(!pd.Series.size()){
    std::cout << "[WARN]: Canvas element: " << pd.Name << " contained no valid"
    " series." << std::endl;
    return;
  }

  std::cout << "[INFO]: " << pd << std::endl;

  AutoPlots.push_back(pd);
}


TCanvas* AutoPlots____RollOut(PlotDescriptor const &plot){

  style::SetGlobalStyle();

  TCanvas* canv = new TCanvas(plot.Name.c_str(),"",600,400);

  style::fgkYTitleOffset = 1.1;
  TGaxis::SetMaxDigits(3);

  style::PadSetup(canv);

  if(plot.HaveXYMargins){
    canv->SetTopMargin(plot.Y2);
    canv->SetBottomMargin(plot.Y1);
    canv->SetLeftMargin(plot.X1);
    canv->SetRightMargin(plot.X2);
  }

  canv->SetLogx(plot.Logs&1);
  canv->SetLogy(plot.Logs&2);
  canv->SetLogz(plot.Logs&4);

  TLegend* leg;
  if(plot.Legend.X1 == 0xdeadbeef){
    leg = LegendFactory();
    std::cout << "[INFO]: Making Legend." << std::endl;
  } else {
    leg = LegendFactory(plot.Legend.X1,plot.Legend.Y1,
      plot.Legend.X2,plot.Legend.Y2);
    leg->SetNColumns(plot.Legend.NColumns);
    std::cout << "[INFO]: Making Legend {" << plot.Legend.X1 << ", "
      << plot.Legend.Y1 << ", " << plot.Legend.X2 << ", " << plot.Legend.Y2
      << "}" << std::endl;
  }
  if(plot.Legend.Title.length()){
    leg->SetHeader(plot.Legend.Title.c_str());
  }
  style::ResetStyle(leg);

  int ctr = 0;
  std::vector<std::pair<TH1*,std::string>> seriesTH1s;
  for(std::vector<SeriesDescriptor>::const_iterator s_it = plot.Series.begin();
    s_it != plot.Series.end(); std::advance(s_it,1)){

    auto const &series = *s_it;
    TH1* plotTH1 = nullptr;
    switch(series.Type){
      case SeriesDescriptor::kTH1:{
        plotTH1 = NoSaveClone(
          PlottingIO::FillHistogram1D(
            series.Generator,series.Sample,series.Selection));
        if(plot.AddNormToTitle){
          AddCountToTHTitle(series.LegendTitle.c_str(),plotTH1);
        } else if(series.LegendTitle.length()){
          plotTH1->SetTitle(series.LegendTitle.c_str());
        }
        break;
      }
      case SeriesDescriptor::kTH1External:{
        plotTH1 = NoSaveClone(
          PlottingIO::GetExternalTH1(
            series.Generator,series.Sample));
        if(plot.AddNormToTitle){
          AddCountToTHTitle(series.LegendTitle.c_str(),plotTH1);
        } else if(series.LegendTitle.length()){
          plotTH1->SetTitle(series.LegendTitle.c_str());
        }
        break;
      }
      case SeriesDescriptor::kTH1RatioNumer:{
        std::advance(s_it,1);
        auto const & denom_series = *s_it;
        TH1* NumerTH1 = NoSaveClone(
          PlottingIO::FillHistogram1D(series.Generator,series.Sample,
            series.Selection));
        TH1* DenomTH1 = NoSaveClone(
          PlottingIO::FillHistogram1D(denom_series.Generator,
            denom_series.Sample, denom_series.Selection));

        plotTH1 = NoSaveClone(NumerTH1);
        plotTH1->Divide(NumerTH1, DenomTH1,1,1,"B");

        if(plot.AddNormToTitle){
          plotTH1->SetTitle((series.LegendTitle + "("
            + GetIntegralAsString(NumerTH1) + "/"
            + GetIntegralAsString(DenomTH1) + ")").c_str());
        } else if(series.LegendTitle.length()){
          plotTH1->SetTitle(series.LegendTitle.c_str());
        }
        break;
      }
      case SeriesDescriptor::kTH1QuadRatio1:{
        TH1* QRat1 = NoSaveClone(
          PlottingIO::FillHistogram1D(series.Generator,series.Sample,
            series.Selection));
        std::advance(s_it,1);

        auto const &next_series1 = *s_it;
        if(next_series1.Type != SeriesDescriptor::kTH1QuadRatio2){
          throw 1;
        }
        TH1* QRat2 = NoSaveClone(
          PlottingIO::FillHistogram1D(next_series1.Generator,
            next_series1.Sample,
            next_series1.Selection));
        std::advance(s_it,1);

        auto const &next_series2 = *s_it;
        if(next_series2.Type != SeriesDescriptor::kTH1QuadRatio3){
          throw 1;
        }
        TH1* QRat3 = NoSaveClone(
          PlottingIO::FillHistogram1D(next_series2.Generator,
            next_series2.Sample,
            next_series2.Selection));
        std::advance(s_it,1);

        auto const &next_series3 = *s_it;
        if(next_series3.Type != SeriesDescriptor::kTH1QuadRatio4){
          throw 1;
        }
        TH1* QRat4 = NoSaveClone(
          PlottingIO::FillHistogram1D(next_series3.Generator,
            next_series3.Sample,
            next_series3.Selection));

        TH1* NumerRat = NoSaveClone(QRat1);
        NumerRat->Divide(QRat1,QRat2,1,1,"B");

        TH1* DenomRat = NoSaveClone(QRat3);
        DenomRat->Divide(QRat3,QRat4,1,1,"B");

        plotTH1 = NoSaveClone(NumerRat);
        plotTH1->Divide(NumerRat, DenomRat,1,1,"B");

        if(plot.AddNormToTitle){
          plotTH1->SetTitle((series.LegendTitle + "("
            + GetIntegralAsString(QRat1) + "/"
            + GetIntegralAsString(QRat2) + ") / ("
            + GetIntegralAsString(QRat3) + "/"
            + GetIntegralAsString(QRat4) + ")").c_str());
        } else if(series.LegendTitle.length()){
          plotTH1->SetTitle(series.LegendTitle.c_str());
        }
        break;
      }
      case SeriesDescriptor::kTH1QuadRatio2:
      case SeriesDescriptor::kTH1QuadRatio3:
      case SeriesDescriptor::kTH1QuadRatio4:
      case SeriesDescriptor::kTH1RatioDenom:{
        std::cout << "[ERROR]: I came across a denominator that should have "
          "been hoovered up and skipped." << std::endl;
        return nullptr;
      }
      case SeriesDescriptor::kTH2Profile:{
        static std::map<std::string, uint> UsedProfileNames;
        TH2* Plot2D =
          NoSaveClone(
            PlottingIO::FillHistogram2D(series.Generator,series.Sample,
            series.Selection));

        std::string pfxName = std::string(Plot2D->GetName())+"_pfx";
        if(UsedProfileNames.count(pfxName)){
          UsedProfileNames[pfxName]++;
          std::stringstream ss;
          ss << UsedProfileNames[pfxName];
          pfxName += ss.str();
        } else {
          UsedProfileNames[pfxName] = 0;
        }

        plotTH1 = Plot2D->ProfileX(pfxName.c_str(),1,-1,"");
        break;
      }
      case SeriesDescriptor::kNormTH2:
      case SeriesDescriptor::kTH2:{
        plotTH1 = NoSaveClone(
          PlottingIO::FillHistogram2D(series.Generator,series.Sample,
            series.Selection));
        if(ctr){
          std::cout << "[ERROR]: Expected Series2D to be plotted alone but it"
           " was not the first plot in the Canvas1D element." << std::endl;
           return nullptr;
        }

        if(series.Type == SeriesDescriptor::kNormTH2){
          plotTH1 = NeutrinoTools::NormalHist(static_cast<TH2*>(plotTH1),0,true);
        }
        //Leave early for TH2
        style::ResetStyle(plotTH1);
        plotTH1->GetXaxis()->SetTitle(plot.XAxisTitle.c_str());
        plotTH1->GetYaxis()->SetTitle(plot.YAxisTitle.c_str());
        plotTH1->Draw("colz");
        TitleFactory(plot.Title)->Draw();

        return canv;
      }
      case SeriesDescriptor::kUnSet:
      default:{
        std::cout << "[ERROR]: I came across an unset or unhandled series type."
          << std::endl;
        return nullptr;
      }
    }

    style::ResetStyle(plotTH1);

    plotTH1->SetLineWidth(series.LineWidth);
    plotTH1->SetLineColor(series.LineColor);
    plotTH1->SetLineStyle(series.LineStyle);

    plotTH1->SetMarkerSize(series.MarkerSize);
    plotTH1->SetMarkerColor(series.MarkerColor);
    plotTH1->SetMarkerStyle(series.MarkerStyle);

    plotTH1->SetFillStyle(series.FillStyle);
    plotTH1->SetFillColor(series.FillColor);

    if(plot.DoPDF){
      MakePDF(plotTH1);
    }

    plotTH1->Scale(series.Scale);

    if(series.LegendTitle.length()){
      AddToLegend({plotTH1},leg);
    }

    if(!ctr){
      plotTH1->GetXaxis()->SetTitle(plot.XAxisTitle.c_str());
      plotTH1->GetYaxis()->SetTitle(plot.YAxisTitle.c_str());
    }


    std::cout << "[INFO]: Adding " << plotTH1 << ", Drawn with: "
      << series << " To draw list." << std::endl;

    //Set up some defaults if the drawing options haven't been changed.
    std::string drawOpt = series.DrawOpt;
    if(drawOpt == "HIST"){
      drawOpt += ( std::string(plotTH1->GetMarkerStyle()?"P":"")
                          + (plotTH1->GetLineStyle()?"H":"")
            + ((plotTH1->GetLineStyle()||plotTH1->GetMarkerStyle())?"E":""));
    }
    seriesTH1s.emplace_back(
      std::make_pair( plotTH1, drawOpt) );
    ctr++;
  }

  SetUserRangeToMaxYVal(seriesTH1s, canv->GetLogy()?0.9:0, plot.YStretch);

  DrawTHs(seriesTH1s);

  if(leg->GetListOfPrimitives()->GetSize()){
    leg->Draw();
  } else {
    delete leg;
  }

  if(plot.Title.length()){
    TitleFactory(plot.Title)->Draw();
  }

  for(auto & line : plot.Lines){
    TLine* ln = new TLine(line.X1, line.Y1, line.X2, line.Y2);
    ln->SetLineColor(line.Color);
    ln->SetLineWidth(line.Width);
    ln->SetLineStyle(line.Style);
    ln->Draw();
    std::cout << "[INFO]: Drawing Line: {" << line.X1 << ","
        << line.Y1 << "," << line.X2 << ","
        << line.Y2 << "}" << std::endl;

    if(line.X3Arrow != 0xdeadbeef){
      Float_t X1,X2,Y1,Y2;
      std::tie(X1,Y1,X2,Y2) = line.GetArrowCoords();
      TArrow* ar = new TArrow(X1,Y1,X2,Y2,0.015,"|>");
      ar->SetLineColor(line.Color);
      ar->SetFillColor(line.Color);
      ar->SetLineWidth(line.Width);
      ar->SetLineStyle(line.Style);
      ar->Draw();

      std::cout << "[INFO]: Drawing Arrow: {" << X1 << ","
        << Y1 << "," << X2 << "," << Y2 << " }" << std::endl;
    }
  }

  return canv;
}

bool MakeAutoPlotsXML(std::string const &fileName){
  TXMLEngine* xmlengine = new TXMLEngine;
  xmlengine->SetSkipComments();

  XMLNodePointer_t confSecNode =
    IOUtils::GetNamedChildElementOfDocumentRoot(*xmlengine, fileName, "Plots");
  if(!confSecNode){
    delete xmlengine; return false;
  }

  for(XMLNodePointer_t canvDescriptNode = xmlengine->GetChild(confSecNode);
    (canvDescriptNode != NULL);
    canvDescriptNode = xmlengine->GetNext(canvDescriptNode)){

    if(std::string(xmlengine->GetNodeName(canvDescriptNode)) != "Canvas1D"){
      continue;
    }

    AddCanvasDescriptionXML(*xmlengine,canvDescriptNode);
  }

  delete xmlengine;

  if(!AutoPlots.size()){
    std::cout << "[WARN]: Couldn't find any valid configured plots."
    << std::endl;
    return true;
  }

  std::cout << "[INFO]: Found " << AutoPlots.size() << " XML configured plots."
  << std::endl;

  //Start Making the plots

  for(auto const &plot : AutoPlots){
    PlottingIO::SaveAndPrint(AutoPlots____RollOut(plot),plot.PrintName);
  }

  return true;
}

}
