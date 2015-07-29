#include <iterator>

#include "TH1.h"
#include "TCanvas.h"

#include "PlottingAutoPlots.hxx"
#include "PlottingSelections.hxx"
#include "style_Xianguo.hxx"
#include "PlottingIO.hxx"


using namespace PlottingTypes;
using namespace PlottingDefaults;
using namespace PlottingUtils;

namespace AutoPlots {
std::vector<PlotDescriptor> AutoPlots;

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

  PlotDescriptor pd(canvName, canvTitle,canvXAxisTitle, canvYAxisTitle,
    canvPrintName, Logs);

  bool ExpectingDenom = false;

  for(XMLNodePointer_t SeriesDescriptor = xmlengine.GetChild(canvDescriptNode);
    (SeriesDescriptor != NULL);
    SeriesDescriptor = xmlengine.GetNext(SeriesDescriptor)){

    SeriesDescriptor::SeriesType type = SeriesDescriptor::kUnSet;

    if(std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Series1D"){
      type = SeriesDescriptor::kTH1;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Ratio1D"){
      type = ExpectingDenom? SeriesDescriptor::kTH1RatioDenom:SeriesDescriptor::kTH1RatioNumer;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Profile2D"){
      type = SeriesDescriptor::kTH2Profile;
    } else if (std::string(xmlengine.GetNodeName(SeriesDescriptor)) == "Series2D"){
      if(pd.Series.size()){
        std::cout << "[ERROR]: I was not expecting a Series2D element unless it"
         " was the first and only." << std::endl;
        return;
      }
      type = SeriesDescriptor::kTH2;
    } else {
      continue;
    }

    std::string const &seriesGenerator =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Generator");
    std::string const &seriesSample =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Sample");
    std::string const &seriesSelection =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "Selection");
    std::string const &seriesLegendTitle =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "LegendTitle");

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

    std::string const &seriesLineStyle =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "LineStyle");
    std::string const &seriesLineColor =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "LineColor");
    std::string const &seriesLineWidth =
      IOUtils::GetAttrValue(xmlengine, SeriesDescriptor, "LineWidth");

    Int_t lineStyle = 0, lineWidth = 0, lineColor = 0;

    if(seriesLineStyle.length()){
      try{
        lineStyle = std::stoi(seriesLineStyle);
      } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR]: Couldn't parse: \"" << seriesLineStyle << "\"\n\t"
          << ia.what() << std::endl;
        continue;
      }
    }
    if(seriesLineColor.length()){
      try{
        lineColor = IOUtils::ParseRootColor(seriesLineColor);
      } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR]: Couldn't parse: \"" << seriesLineColor << "\"\n\t"
          << ia.what() << std::endl;
        continue;
      }
    }

    if(seriesLineWidth.length()){
      try{
        lineWidth = std::stoi(seriesLineWidth);
      } catch (const std::invalid_argument& ia) {
        std::cout << "[ERROR]: Couldn't parse: \"" << seriesLineWidth << "\"\n\t"
          << ia.what() << std::endl;
        continue;
      }
    }

    pd.Series.emplace_back(
      seriesGenerator,
      seriesSample,
      seriesSelection,
      seriesLegendTitle,
      lineColor,
      lineWidth,
      lineStyle,
      type);

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
  TCanvas* canv = CanvasFactory(plot.Name.c_str());
  style::PadSetup(canv);

  canv->SetLogx(plot.Logs&1);
  canv->SetLogy(plot.Logs&2);
  canv->SetLogz(plot.Logs&4);

  TLegend* leg = LegendFactory( (0.12 + P0LegendWidth*0.1),
    (1.0 - Pad0YMargin*0.2 - P0LegendHeight),
    (0.12 + P0LegendWidth*1.1),
    (1.0 - Pad0YMargin*0.2));
  // style::ResetStyle(leg);

  int ctr = 0;
  std::vector<TH1*> seriesTH1s;
  for(std::vector<SeriesDescriptor>::const_iterator s_it = plot.Series.begin();
    s_it != plot.Series.end(); std::advance(s_it,1)){

    auto const &series = *s_it;
    TH1* plotTH1 = nullptr;
    switch(series.Type){
      case SeriesDescriptor::kTH1:{
        plotTH1 = NoSaveClone(
          PlottingIO::FillHistogram1D(
            series.Generator,series.Sample,series.Selection));
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
        break;
      }
      case SeriesDescriptor::kTH1RatioDenom:{
        std::cout << "[ERROR]: I came across a denominator that should have "
          "been hoovered up and skipped." << std::endl;
        return nullptr;
      }
      case SeriesDescriptor::kTH2Profile:{
        TH2* Plot2D =
          NoSaveClone(
            PlottingIO::FillHistogram2D(series.Generator,series.Sample,
            series.Selection));
        plotTH1 = (TH1*)Plot2D->ProfileX();
      }
      case SeriesDescriptor::kTH2:{
        plotTH1 = NoSaveClone(
          PlottingIO::FillHistogram2D(series.Generator,series.Sample,
            series.Selection));
        if(ctr){
          std::cout << "[ERROR]: Expected Series2D to be plotted alone but it"
           " was not the first plot in the Canvas1D element." << std::endl;
           return nullptr;
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

    if(series.LegendTitle.length()){
      leg->AddEntry(plotTH1, series.LegendTitle.c_str(), "l");
    }

    if(!ctr){
      plotTH1->GetXaxis()->SetTitle(plot.XAxisTitle.c_str());
      plotTH1->GetYaxis()->SetTitle(plot.YAxisTitle.c_str());
    }

    seriesTH1s.push_back(plotTH1);
    ctr++;
  }

  SetUserRangeToMaxYVal(seriesTH1s, canv->GetLogy()?0.9:0);

  DrawTHs(seriesTH1s,"HIST");

  if(leg->GetListOfPrimitives()->GetSize()){
    leg->Draw();
  } else {
    delete leg;
  }

  TitleFactory(plot.Title)->Draw();

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
