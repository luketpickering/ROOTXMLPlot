#include <algorithm>
#include "PlottingIO.hxx"
#include "TH2F.h"

using namespace PlottingUtils;
using namespace PlottingTypes;

//**********************************Globals*************************************
namespace Data {
  bool HaveTrees = false;
  std::map< std::string,
            std::map< std::string,
                      TFile*>
          > Files;
  std::map< std::string,
            std::map< std::string,
                      TTree*>
          > Trees;
  std::map< std::string,
            std::map< std::string,
                      std::map< std::string,
                                TH1*
                              >
                    >
          > Histos;
}

namespace DataSpecifics {

  std::vector<PlottingTypes::Generator> Generators;
  std::map<std::string, PlottingTypes::Selection1D> Selections1D;
  std::map<std::string, PlottingTypes::Selection2D> Selections2D;

//**********************************CCQEFSI*************************************
  std::string MuonValidSelect;
  std::string FirstProtonValidSelect;
  std::string FirstProtonNonZeroSelect;
  std::string CCQERCSelect;
  std::string CCQEDPhiTSelect;
  std::string CCQENoCCQENucEmitSelect;
  std::string CCQENucEmitSelect;
  std::string CCQESelect;
  std::string MECSelect;
  std::string CCQELikeSelect;
  std::string CC_not_QELikeSelect;
  std::string CCQEProtonMomentumSelect;

  std::string FirstPiPlusValidSelect;
  std::string FirstPiPlusNonZeroSelect;

  std::string ResPPiPlusProdSelect;
  std::string ResNPiPlusProdSelect;
  std::string CohPiPlusProdSelect;
  std::string OtherPiPlusSelect;
  std::string CCSelect;

  std::string CCResSelect;
  std::string CCOtherSelect;
//**********************************Globals*************************************

  void InitTCuts(){

    static std::string const MomentumCut_MeV = "10";

    MuonValidSelect = "(MuonProtonTransversity.MuonPDG==13)";

    FirstProtonValidSelect = "(MuonProtonTransversity.FirstProtonPDG==2212)";
    FirstProtonNonZeroSelect =
      "(MuonProtonTransversity.FirstProtonMomentum_MeV > "
        + MomentumCut_MeV +")";
    CCQERCSelect = "(MuonProtonTransversity.NeutConventionReactionCode==1)";

    CCQEDPhiTSelect = FirstProtonNonZeroSelect;

    CCQENoCCQENucEmitSelect = MuonValidSelect + "&&" + FirstProtonValidSelect +
      "&&(MuonProtonTransversity.NFinalStateParticles==2)&&" + CCQERCSelect;

    CCQENucEmitSelect = "(!("+CCQENoCCQENucEmitSelect+"))&&"+CCQERCSelect;

    CCQESelect =  MuonValidSelect+"&&"+FirstProtonValidSelect+"&&"+CCQERCSelect;

    MECSelect = "(MuonProtonTransversity.NeutConventionReactionCode==2)&&"
      + MuonValidSelect + "&&" + FirstProtonValidSelect;
    CCQELikeSelect = "(" + CCQERCSelect
        + "||(MuonProtonTransversity.NeutConventionReactionCode==2))&&"
        + MuonValidSelect + "&&" + FirstProtonValidSelect;
    CC_not_QELikeSelect = "((MuonProtonTransversity.NeutConventionReactionCode"
      ">2)&&(MuonProtonTransversity.NeutConventionReactionCode<30))&&"
       + MuonValidSelect + "&&" + FirstProtonValidSelect;

    ResPPiPlusProdSelect =
      "(PionProductionTransversity.NeutConventionReactionCode==11)";
    ResNPiPlusProdSelect =
      "(PionProductionTransversity.NeutConventionReactionCode==13)";
    CohPiPlusProdSelect =
      "(PionProductionTransversity.NeutConventionReactionCode==16)";
    CCSelect = "(PionProductionTransversity.NeutConventionReactionCode<30)";

    OtherPiPlusSelect =
      "((PionProductionTransversity.NeutConventionReactionCode!=11)"
      "&&(PionProductionTransversity.NeutConventionReactionCode==13)"
      "&&(PionProductionTransversity.NeutConventionReactionCode!=16))";

    FirstPiPlusValidSelect = "(PionProductionTransversity.FirstPiPlusPDG==211)";
    FirstPiPlusNonZeroSelect =
      "(PionProductionTransversity.FirstPiPlusMomentum_MeV > " + MomentumCut_MeV
        + ")";

    CCQEProtonMomentumSelect = MuonValidSelect
      + "&&" + FirstProtonValidSelect
      + "&&" + CCQERCSelect;

    CCResSelect = MuonValidSelect + "&&"
      + FirstProtonValidSelect
      + "&&((MuonProtonTransversity.NeutConventionReactionCode>=11)"
        "&&(MuonProtonTransversity.NeutConventionReactionCode<=13))";
    CCOtherSelect = MuonValidSelect + "&&"
      + FirstProtonValidSelect
      + "&&((MuonProtonTransversity.NeutConventionReactionCode>=16)"
          "&&(MuonProtonTransversity.NeutConventionReactionCode<=30))";

  }

  void SetProtonMomentaSelections(){
    Selections1D["CCQEMuonMom"] =
      (Selection1D("CCQEMuonMom","CCQEMuonMom",
      "PionProductionTransversity.MuonMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
        CCQESelect.c_str(),
        DoLogx));

    Selections1D["CCQEHMProtonMom"] =
      (Selection1D("CCQEHMProtonMom","CCQEHMProtonMom",
      "MuonProtonTransversity.HMProtonMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
        CCQESelect.c_str(),
        DoLogx));

    Selections1D["CCQEMuon_Pt"] =
      (Selection1D("CCQEMuon_Pt","CCQEMuon_Pt",
      "MuonProtonTransversity.MuonPt_MeV.Mag()",
      ProtonMomBins, 0, ProtonMomMax,
        CCQESelect.c_str(),
        DoLogx));

    Selections1D["CCQEHMProton_Pt"] =
      (Selection1D("CCQEHMProton_Pt","CCQEHMProton_Pt",
      "MuonProtonTransversity.HMProtonPt_MeV.Mag()",
      ProtonMomBins, 0, ProtonMomMax,
        CCQESelect.c_str(),
        DoLogx));

    Selections1D["CCQEHMProton_PHI__"] =
      (Selection1D("CCQEHMProton_PHI__","CCQEHMProton_PHI__",
      "MuonProtonTransversity.HMProtonDirection.Phi()",
      100, -1.0*M_PI, M_PI,
      (CCQESelect + "&&" + FirstProtonNonZeroSelect).c_str(),
      DoLogx));
    Selections1D["CCQEHMProton_Phi__"] =
      (Selection1D("CCQEHMProton_Pt__","CCQEHMProton_Pt__",
      "MuonProtonTransversity.HMProtonDirection.CosTheta()",
      100, -1, 1,
      (CCQESelect + "&&" + FirstProtonNonZeroSelect).c_str(),
      DoLogx));
  }

  void SetModeSelections(){
    Selections1D["MEC"] =
      (Selection1D("MEC","MEC",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + MECSelect).c_str(), DoLogx));

    Selections1D["MECZoom"] =
      (Selection1D("MECZoom","MECZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + MECSelect).c_str(), DoLogx));

    // Selections1D["CCQELike"] =
    //   (Selection1D("CCQELike","CCQELike",
    //   "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
    //   DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
    //   (CCQEDPhiTSelect + "&&" + CCQELikeSelect).c_str(), DoLogx));

    // Selections1D["CCQELikeZoom"] =
    //   (Selection1D("CCQELikeZoom","CCQELikeZoom",
    //   "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
    //   DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
    //   (CCQEDPhiTSelect + "&&" + CCQELikeSelect).c_str(), DoLogx));

    // Selections1D["CC_not_QELike"] =
    //   (Selection1D("CC_not_QELike","CC_not_QELike",
    //   "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
    //   DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
    //   (CCQEDPhiTSelect + "&&" + CC_not_QELikeSelect).c_str(), DoLogx));

    // Selections1D["CC_not_QELikeZoom"] =
    //   (Selection1D("CC_not_QELikeZoom","CC_not_QELikeZoom",
    //   "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
    //   DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
    //   (CCQEDPhiTSelect + "&&" + CCQELikeSelect).c_str(), DoLogx));

    Selections1D["CCRes"] =
      (Selection1D("CCRes","CCRes",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCResSelect).c_str(), DoLogx));

    Selections1D["CCResZoom"] =
      (Selection1D("CCResZoom","CCResZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + CCResSelect).c_str(), DoLogx));

    Selections1D["CCOther"] =
      (Selection1D("CCOther","CCOther",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCOtherSelect).c_str(), DoLogx));

    Selections1D["CCOtherZoom"] =
      (Selection1D("CCOtherZoom","CCOtherZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + CCOtherSelect).c_str(), DoLogx));
  }

  void TargetMassSelection(){

    static std::string const selccqe =
      "&&(MuonProtonTransversity.NeutConventionReactionCode==1)";
    static std::string const selccqeNE =
      "&&(MuonProtonTransversity.NeutConventionReactionCode==1)"
      "&&((MuonProtonTransversity.FirstProtonPDG!=2212)||"
        "(MuonProtonTransversity.NFinalStateParticles!=2))";
    static std::string const selmec =
      "&&(MuonProtonTransversity.NeutConventionReactionCode==2)";
    static std::string const selres =
      "&&(MuonProtonTransversity.NeutConventionReactionCode>=11)"
      "&&(MuonProtonTransversity.NeutConventionReactionCode<=13)";
    static std::string const selcoh =
      "&&(MuonProtonTransversity.NeutConventionReactionCode==16)";
    static std::string const selccother =
      "&&(MuonProtonTransversity.NeutConventionReactionCode>=17)"
      "&&(MuonProtonTransversity.NeutConventionReactionCode<=23)";
    static std::string const seldis =
      "&&(MuonProtonTransversity.NeutConventionReactionCode==26)";


    Selections1D["ReconTargetMassselccqe"] =
      Selection1D("ReconTargetMass","ReconTargetMassselccqe",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + selccqe).c_str(),DoLogx);
    Selections1D["ReconTargetMassselccqeNE"] =
      Selection1D("ReconTargetMassNE","ReconTargetMassselccqeNE",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + selccqeNE).c_str(),DoLogx);
    Selections1D["ReconTargetMassselmec"] =
      Selection1D("ReconTargetMass","ReconTargetMassselmec",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + selmec).c_str(),DoLogx);
    // Selections1D["ReconTargetMassselres"] =
    //   Selection1D("ReconTargetMass","ReconTargetMassselres",
    //     "MuonProtonTransversity.ReconTargetMass",
    //     100,0,2000,
    //     (MuonValidSelect
    //       + "&&" + FirstProtonValidSelect
    //       + "&&" + FirstProtonNonZeroSelect
    //       + selres).c_str(),DoLogx);
    // Selections1D["ReconTargetMassselcoh"] =
    //   Selection1D("ReconTargetMass","ReconTargetMassselcoh",
    //     "MuonProtonTransversity.ReconTargetMass",
    //     100,0,2000,
    //     (MuonValidSelect
    //       + "&&" + FirstProtonValidSelect
    //       + "&&" + FirstProtonNonZeroSelect
    //       + selcoh).c_str(),DoLogx);
    // Selections1D["ReconTargetMassselccother"] =
    //   Selection1D("ReconTargetMass","ReconTargetMassselccother",
    //     "MuonProtonTransversity.ReconTargetMass",
    //     100,0,2000,
    //     (MuonValidSelect
    //       + "&&" + FirstProtonValidSelect
    //       + "&&" + FirstProtonNonZeroSelect
    //       + selccother).c_str(),DoLogx);
    // Selections1D["ReconTargetMassseldis"] =
    //   Selection1D("ReconTargetMass","ReconTargetMassseldis",
    //     "MuonProtonTransversity.ReconTargetMass",
    //     100,0,2000,
    //     (MuonValidSelect
    //       + "&&" + FirstProtonValidSelect
    //       + "&&" + FirstProtonNonZeroSelect
    //       + seldis).c_str(),DoLogx);
  }

  void SetCCQESelections(){

    Selections1D["CCQE"] =
      (Selection1D("CCQE","CCQE",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

    Selections1D["CCQEZoom"] =
      (Selection1D("CCQE All Zoom","CCQEZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

    Selections1D["CCQENE"] =
      (Selection1D("CCQE Nuclear Emission","CCQENE",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

    Selections1D["NEZoom"] =
      (Selection1D("Nuclear Emission Zoom","NEZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

    Selections1D["CCQENoNE"] =
      (Selection1D("CCQE No Nuclear Emission","CCQENoNE",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

    Selections1D["CCQENoNEZoom"] =
      (Selection1D("CCQE No Nuclear Emission Zoom","CCQENoNEZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//DeltaPT
//     Selections1D["CCQEPt"] =
//       (Selection1D("CCQE Pt","CCQEPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENEPt"] =
//       (Selection1D("CCQE Nuclear Emission Pt","CCQENEPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQENoNEPt"] =
//       (Selection1D("CCQE No Nuclear Emission Pt","CCQENoNEPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

// //AlphaT
//     Selections1D["CCQEAlphaT"] =
//       (Selection1D("CCQE AlphaT","CCQEAlphaT",
//       "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENEAlphaT"] =
//       (Selection1D("CCQE Nuclear Emission AlphaT","CCQENEAlphaT",
//       "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQENoNEAlphaT"] =
//       (Selection1D("No Nuclear Emission","CCQENoNEAlphaT",
//       "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//     Selections2D["DeltaPhiTHMPvsHMPpt"] =
//       (Selection2D("DeltaPhiTHMPvsHMPpt","DeltaPhiTHMPvsHMPpt",
//         "MuonProtonTransversity.HMProtonPt_MeV.Mag():"
//         "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//         DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//         100,0,1000,(CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DeltaPhiTHMPvsMuonPt"] =
//       (Selection2D("DeltaPhiTHMPvsMuonPt","DeltaPhiTHMPvsMuonPt",
//         "MuonProtonTransversity.MuonPt_MeV.Mag():"
//         "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//         DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//         100,0,1000,(CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

  }

  void FinalStateSpeciesSelections(std::string drawObject,
    std::string selectString,
    std::string prefix=""){

    constexpr int NProtonsMax = 50;
    constexpr int NGammaMax = 5;
    constexpr int NNeutronsMax = 50;
    constexpr int NPiPlusMax = 5;
    constexpr int NPiZeroMax = 5;
    constexpr int NPiMinusMax = 5;
    constexpr int NFinalStateParticles = 100;

    Selections1D[prefix + "NProtons"] =
      (Selection1D("NProtons", prefix + "NProtons",
      drawObject + ".NProtons",
      NProtonsMax, 0, NProtonsMax, selectString.c_str(), DoLogx, true));

    Selections1D[prefix + "NGammma"] =
      (Selection1D("NGamma", prefix + "NGammma",
      drawObject + ".NGamma",
      NGammaMax, 0, NGammaMax, selectString.c_str(), DoLogx, true));

    Selections1D[prefix + "NNeutron"] =
      (Selection1D("NNeutrons", prefix + "NNeutron",
      drawObject + ".NNeutrons",
      NNeutronsMax, 0, NNeutronsMax, selectString.c_str(), DoLogx, true));

    Selections1D[prefix + "NPiPlus"] =
      (Selection1D("NPiPlus", prefix + "NPiPlus",
      drawObject + ".NPiPlus",
      NPiPlusMax, 0, NPiPlusMax, selectString.c_str(), DoLogx, true));

    Selections1D[prefix + "NPiZero"] =
      (Selection1D("NPiZero", prefix + "NPiZero",
      drawObject + ".NPiZero",
      NPiZeroMax, 0, NPiZeroMax, selectString.c_str(), DoLogx, true));

    Selections1D[prefix + "NPiMinus"] =
      (Selection1D("NPiMinus", prefix + "NPiMinus",
      drawObject + ".NPiMinus",
      NPiMinusMax, 0, NPiMinusMax, selectString.c_str(), DoLogx, true));

    Selections1D[prefix + "NPrimaryParticles"] =
      (Selection1D("NFinalStateParticles", prefix + "NPrimaryParticles",
      drawObject + ".NFinalStateParticles",
      NFinalStateParticles, 0, NFinalStateParticles, selectString.c_str(),
      DoLogx, true));
  }

  void DetailsPiProdSelection(){
    Selections1D["Res11"] =
      (Selection1D("Res11","Res11",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["Res11Zoom"] =
      (Selection1D("Res11Zoom","Res11Zoom",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["Res13"] =
      (Selection1D("Res13","Res13",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["Res13Zoom"] =
      (Selection1D("Res13Zoom","Res13Zoom",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["PiProdOther"] =
      (Selection1D("PiProdOther","PiProdOther",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (OtherPiPlusSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["PiProdOtherZoom"] =
      (Selection1D("PiProdOtherZoom","PiProdOtherZoom",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (OtherPiPlusSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//COH Momenta START
    Selections1D["CohMuonMom"] =
      (Selection1D("CohMuonMom","CohMuonMom",
      "PionProductionTransversity.MuonMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
        "(PionProductionTransversity.NeutConventionReactionCode==16)",
        DoLogx));

    Selections1D["CohHMPiPlusMom"] =
      (Selection1D("CohHMPiPlusMom","CohHMPiPlusMom",
      "PionProductionTransversity.HMPiPlusMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
        "(PionProductionTransversity.NeutConventionReactionCode==16)",
        DoLogx));

    Selections1D["CohMuon_Pt"] =
      (Selection1D("CohMuon_Pt","CohMuon_Pt",
      "PionProductionTransversity.MuonPt_MeV.Mag()",
      ProtonMomBins, 0, ProtonMomMax,
        "(PionProductionTransversity.NeutConventionReactionCode==16)",
        DoLogx));

    Selections1D["CohHMPiPlus_Pt"] =
      (Selection1D("CohHMPiPlus_Pt","CohHMPiPlus_Pt",
      "PionProductionTransversity.HMPiPlusPt_MeV.Mag()",
      ProtonMomBins, 0, ProtonMomMax,
        "(PionProductionTransversity.NeutConventionReactionCode==16)",
        DoLogx));
  }

  void BasicPiProdSelections(){
    Selections1D["CCPiProdInc"] =
      (Selection1D("CCPiProdInc","CCPiProdInc",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["CCPiProdIncZoom"] =
      (Selection1D("CCPiProdIncZoom","CCPiProdIncZoom",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["CohPiProd"] =
      (Selection1D("CohPiProd","CohPiProd",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["CohPiProdZoom"] =
      (Selection1D("Coherent Pion Production Zoom","CohPiProdZoom",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["NonCoh"] =
      (Selection1D("NonCoh","NonCoh",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
        + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str()
      , DoLogx));

    Selections1D["NonCohZoom"] =
      (Selection1D("NonCohZoom","NonCohZoom",
      "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
        + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str(),
      DoLogx));

    Selections1D["CohPiProdPt"] =
      (Selection1D("Coh Pt","CohPiProdPt",
      "PionProductionTransversity.DeltaPT_HMPiPlus_MeV.Mag()",
      DPTBins, DPTMin, DPTMax,
      (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["CohPiProdAlphaT"] =
      (Selection1D("Coh AlphaT","CohPiProdAlphaT",
      "PionProductionTransversity.DeltaAlphaT_HMPiPlus_deg",
      DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
      (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
        + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

    Selections1D["NonCohPt"] =
      (Selection1D("NonCoh Pt","NonCohPt",
      "PionProductionTransversity.DeltaPT_HMPiPlus_MeV.Mag()",
      DPTBins, DPTMin, DPTMax,
      ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
        + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str(),
      DoLogx));

    Selections1D["NonCohAlphaT"] =
      (Selection1D("NonCoh AlphaT","NonCohAlphaT",
      "PionProductionTransversity.DeltaAlphaT_HMPiPlus_deg",
      DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
      ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
        + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str(),
      DoLogx));

  }

  ///You edit what things are loaded, and thus plottable, here.
  void SetSpecifics(){
    InitTCuts();
    //Add NEUT Targets
    Generators.emplace_back("NEUT", NEUT_Line);

    Generators.back().Targets.emplace_back("C", C_Color,
      "../neut/C/neut_C_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("C_norescat", C_norescat_Color,
      "../neut/C_norescat/neut_C_norescat_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("C_nofm", C_hN_Color,
    //   "../neut/C_nofm/neut_C_nofm_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("O", O_Color,
      "../neut/O/neut_O_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("Pb", Pb_Color,
      "../neut/Pb/neut_Pb_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("Fe", Fe_Color,
    //   "../neut/Fe/neut_Fe_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("Fe_norsct", Fe_norescat_Color,
    //   "../neut/Fe_norescat/neut_Fe_norescat_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("Fe_nofm", Fe_nofm_Color,
      // "../neut/Fe_nofm/neut_Fe_nofm_puresim_transversityVars.root");

    //Add GENIE Targets
    Generators.emplace_back("GENIE", GENIE_Line);
    Generators.back().Targets.emplace_back("C", C_Color,
      "../genie/C/genie_C_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("C_norescat", C_norescat_Color,
      "../genie/C_norescat/genie_C_norescat_puresim_transversityVars.root");
     Generators.back().Targets.emplace_back("C_hN", C_hN_Color,
      "../genie/C_hN/genie_C_hN_puresim_transversityVars.root");
     // Generators.back().Targets.emplace_back("C_2_9", C_2_9_Color,
     //  "../genie/C_2_9/genie_C_2_9_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("O", O_Color,
      "../genie/O/genie_O_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("Pb", Pb_Color,
      "../genie/Pb/genie_Pb_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("Fe", Fe_Color,
    //   "../genie/Fe/genie_Fe_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("Fe_norsct", Fe_norescat_Color,
    //   "../genie/Fe_norescat/genie_Fe_norescat_puresim_transversityVars.root");

    Generators.emplace_back("NuWro", NuWro_Line);
    Generators.back().Targets.emplace_back("C", C_Color,
      "../nuwro/C/nuwro_C_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("C_norescat", C_norescat_Color,
      "../nuwro/C_norescat/nuwro_C_norescat_puresim_transversityVars.root");
    // Generators.back().Targets.emplace_back("C_nofm", C_nofm_Color,
    //   "../nuwro/C_nofm/nuwro_C_nofm_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("Pb", Pb_Color,
      "../nuwro/Pb/nuwro_Pb_puresim_transversityVars.root");

    Generators.emplace_back("GiBUU", NuWro_Line);
    Generators.back().Targets.emplace_back("C", C_Color,
      "../gibuu/C/gibuu_C_puresim_transversityVars.root");
    Generators.back().Targets.emplace_back("Pb", Pb_Color,
      "../gibuu/Pb/gibuu_Pb_puresim_transversityVars.root");

    SetCCQESelections();
    SetModeSelections();
    FinalStateSpeciesSelections("MuonProtonTransversity",
        CCQESelect, "ccqe");
    SetProtonMomentaSelections();
    // BasicPiProdSelections();
    TargetMassSelection();
    // DetailsPiProdSelection();

    ///Proton Momenta
    Selections1D["FirstPMom"] =
      (Selection1D("FirstPMom","FirstPMom",
      "MuonProtonTransversity.FirstProtonMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
      CCQEProtonMomentumSelect.c_str(), DoLogx));

    Selections1D["HighestPMom"] =
      (Selection1D("HighestPMom","HighestPMom",
      "MuonProtonTransversity.HMProtonMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
      CCQEProtonMomentumSelect.c_str(), DoLogx));

  }
}

namespace PlottingIO {

  bool LoadFilesAndTrees(){
    for(auto const &Gen : DataSpecifics::Generators){
      for(auto const &Tar : Gen.Targets){

        if(!(Data::Files[Gen.Name][Tar.Name] =
            new TFile(Tar.FileLocation.c_str()))){
          std::cerr << "Failed to load file: " << Tar.FileLocation
            << " for Generator: " << Gen.Name << ", Target: " << Tar.Name
            << std::endl;
            return false;
        }
        if( !(Data::Trees[Gen.Name][Tar.Name] =
              dynamic_cast<TTree*>(
                Data::Files[Gen.Name][Tar.Name]\
                ->Get(DataSpecifics::VarTreeName)) ) ){
          std::cerr << "Failed to load TTree: for Generator: " << Gen.Name
            << ", Target: " << Tar.Name << std::endl;
            return false;
        }

      }
    }
    return true;
  }
  std::string AddOffSetToDrawString(std::string draw, double off){
    std::stringstream ss;
    ss << "(" + draw + "+";
    ss << off << ")";
    return ss.str();
  }
  bool FillHistogram(Generator const & Gen, Target const & Tar,
    Selection1D const & Sel){
    std::string HistName = MakeHistoName(Gen, Tar, Sel);

    std::cout << "\t\tSel: " << Sel.Name
      << "\n\t\t  Draw: \"" << Sel.DrawString
      << "\"\n\t\t  HistName: \"" << HistName
      << "\"\n\t\t  Cut: \""
        << ConfineWidthStringWithIndent("\t\t\t",Sel.Cut.GetTitle(),50)
      <<"\"\n\t\t  Bins: [ " << Sel.XBinLow << " .. " << Sel.NBins << " .. "
      << Sel.XBinUpper << " ]" << std::endl;

    std::string PlotDrawString = Sel.DrawString;

    ///Make the TH1F
    if(Sel.DoPerUseXOffset){ //If we want each Generator to be slightly
      //offset from each other
      size_t TarNum = std::distance(
        Gen.Targets.begin(),
        std::find(Gen.Targets.begin(),
          Gen.Targets.end(),Tar));
      if(TarNum == Gen.Targets.size()){
        std::cout << "[ERROR] Couldn't find this Target descriptor: "
          << Gen.Name << std::endl;
      }
      Double_t offset = ((Sel.XBinUpper - Sel.XBinLow)/Sel.NBins)*0.005*TarNum;
      PlotDrawString = AddOffSetToDrawString(Sel.DrawString, offset);
      try{
        Data::Histos[Gen.Name][Tar.Name][Sel.Tag] =
          new TH1F(HistName.c_str(),"",
            Sel.NBins,Sel.XBinLow+offset,Sel.XBinUpper+offset);
      } catch(std::bad_alloc& ba){
        std::cerr << "Caught bad alloc, exiting grace-ish-fully."
          << std::endl;
          return false;
      }
    } else { //Easy
      try{
        Data::Histos[Gen.Name][Tar.Name][Sel.Tag] =
          new TH1F(HistName.c_str(),"",
            Sel.NBins,Sel.XBinLow,Sel.XBinUpper);
      } catch(std::bad_alloc& ba){
        std::cerr << "Caught bad alloc, exiting grace-ish-fully."
          << std::endl;
          return false;
      }
    }

    if(Sel.UseLogX){
      PlottingUtils::BinLogX(Data::Histos[Gen.Name][Tar.Name][Sel.Tag]\
        ->GetXaxis());
    }
    ///Fill it
    Data::Trees[Gen.Name][Tar.Name]->Draw(
      (PlotDrawString + " >> " + HistName).c_str(),
      Sel.Cut);
    std::cout << "\t\t\"" << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetName()
      << "\" Contained: [Uf:"
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetBinContent(0) << "]"
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->Integral()
      << "[Of:" << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetBinContent(
          Sel.NBins+1)
      << "]"
      << std::endl << std::endl;
    return true;
  }

  bool FillHistogram(Generator const & Gen, Target const & Tar,
    Selection2D const & Sel){
    std::string HistName = MakeHistoName(Gen, Tar, Sel);

    std::cout << "\t\tSel: " << Sel.Name
      << "\n\t\t  Draw: \"" << Sel.DrawString
      << "\"\n\t\t  HistName: \"" << HistName
      << "\"\n\t\t  Cut: \""
        << ConfineWidthStringWithIndent("\t\t\t",Sel.Cut.GetTitle(),50)
      <<"\"\n\t\t  XBins: [ " << Sel.XBinLow << " .. " << Sel.NXBins << " .. "
      << Sel.XBinUpper << " ]"
      <<"\"\n\t\t  YBins: [ " << Sel.YBinLow << " .. " << Sel.NYBins << " .. "
      << Sel.YBinUpper << " ]" << std::endl;

    ///Make the TH1F

    try {
      Data::Histos[Gen.Name][Tar.Name][Sel.Tag] =
        new TH2F(HistName.c_str(),"",
          Sel.NXBins,Sel.XBinLow,Sel.XBinUpper,
          Sel.NYBins,Sel.YBinLow,Sel.YBinUpper);
    } catch(std::bad_alloc& ba){
      std::cerr << "Caught bad alloc, exiting grace-ish-fully."
        << std::endl;
        return false;
    }

    if(Sel.UseLogX){
      PlottingUtils::BinLogX(Data::Histos[Gen.Name][Tar.Name][Sel.Tag]\
        ->GetXaxis());
    }
    if(Sel.UseLogY){
      PlottingUtils::BinLogX(Data::Histos[Gen.Name][Tar.Name][Sel.Tag]\
        ->GetYaxis());
    }
    ///Fill it
    Data::Trees[Gen.Name][Tar.Name]->Draw(
      (Sel.DrawString + " >> " + HistName).c_str(), Sel.Cut);
    std::cout << "\t\t\""
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->GetName()
      << "\" Contained: "
      << Data::Histos[Gen.Name][Tar.Name][Sel.Tag]->Integral()
      << std::endl << std::endl;
    return true;
  }

  bool LoadHistogramsFromFile(char const * HistogramCacheFileName){
    TFile* InputFile = new TFile(HistogramCacheFileName, "UPDATE");
    InputFile->cd();
    if(!InputFile->IsOpen()){
      std::cout << "Couldnt open the histogram cache file: "
        << HistogramCacheFileName << std::endl;
      return false;
    }
    for(auto const &Gen : DataSpecifics::Generators){
      std::cout << "Filling Generator: " << Gen.Name << std::endl;
      for(auto const &Tar : Gen.Targets){
        std::cout << "\tFilling Target: " << Tar.Name << std::endl;

        for(auto const &selpair1D : DataSpecifics::Selections1D){
          Selection1D const &Sel1D = selpair1D.second;

          if(!(Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag] =
            dynamic_cast<TH1F*>(InputFile->Get(
              MakeHistoName(Gen,Tar,Sel1D).c_str())))){
            std::cout << "[NOTE]: Couldn't find Histo: "
              << MakeHistoName(Gen,Tar,Sel1D)
              << " in file: " << HistogramCacheFileName << std::endl;
            if(Data::HaveTrees){
              if(!FillHistogram(Gen,Tar,Sel1D)){
                InputFile->Write();
                InputFile->Close();
                return false;
              }
              Gen.SetSeriesStyle(
                Tar.SetSeriesStyle(
                  Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag]));
              Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag]->Write();
            } else {
              std::cout << "[WARN]: We don't seem to have the data files, "
                "cannot load this histogram." << std::endl;
            }
          } else {
            std::cout << "\t\t" << Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag]\
            ->GetName() << std::endl;
            Gen.SetSeriesStyle(
              Tar.SetSeriesStyle(
                Data::Histos[Gen.Name][Tar.Name][Sel1D.Tag]));
          }
        }
        for(auto const &selpair2D : DataSpecifics::Selections2D){
          Selection2D const &Sel2D = selpair2D.second;
          if(!(Data::Histos[Gen.Name][Tar.Name][Sel2D.Tag] =
            dynamic_cast<TH2F*>(InputFile->Get(
              MakeHistoName(Gen,Tar,Sel2D).c_str())))){
            std::cerr << "Couldn't find Histo: "
              << MakeHistoName(Gen,Tar,Sel2D)
              << " in file: " << HistogramCacheFileName << std::endl;
            if(Data::HaveTrees){
              if(!FillHistogram(Gen,Tar,Sel2D)){
                InputFile->Write();
                InputFile->Close();
                return false;
              }
              Data::Histos[Gen.Name][Tar.Name][Sel2D.Tag]->Write();
            } else {
              std::cout << "We don't seem to have the data files, cannot load"
              " this histogram." << std::endl;
            }
          } else {
            std::cout << "\t\t" << Data::Histos[Gen.Name][Tar.Name][Sel2D.Tag]\
            ->GetName() << std::endl;
          }
        }
      }
    }
    InputFile->Write();
    return true;
  }

  bool LoadHistograms(char const* HistogramCacheFileName){
    DataSpecifics::SetSpecifics();
    Data::HaveTrees = LoadFilesAndTrees();
    return LoadHistogramsFromFile(HistogramCacheFileName);
  }
}
