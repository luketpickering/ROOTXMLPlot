#include <map>
#include <string>
#include <iostream>
#include <algorithm>

#include "TH1.h"
#include "TH2.h"

#include "PlottingSelections.hxx"

using namespace PlottingUtils;
using namespace PlottingTypes;
using namespace PlottingDefaults;

namespace PlottingSelections {

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

  PlottingTypes::Selection1D const * FindSel1D(std::string name){
    if(Selections1D.count(name)){
      return &Selections1D[name];
    }
    std::cout << "[WARN]: Couldn't find Selection1D named: " << name
      << std::endl;
    return 0;
  }

  PlottingTypes::Selection2D const * FindSel2D(std::string name){
    if(Selections2D.count(name)){
      return &Selections2D[name];
    }
    std::cout << "[WARN]: Couldn't find Selection2D named: " << name
      << std::endl;
    return 0;
  }

  void InitCutStrings(){

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

    Selections1D["CCQELike"] =
      (Selection1D("CCQELike","CCQELike",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQELikeSelect).c_str(), DoLogx));

    Selections1D["CCQELikeZoom"] =
      (Selection1D("CCQELikeZoom","CCQELikeZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQELikeSelect).c_str(), DoLogx));

    Selections1D["CC_not_QELike"] =
      (Selection1D("CC_not_QELike","CC_not_QELike",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (CCQEDPhiTSelect + "&&" + CC_not_QELikeSelect).c_str(), DoLogx));

    Selections1D["CC_not_QELikeZoom"] =
      (Selection1D("CC_not_QELikeZoom","CC_not_QELikeZoom",
      "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
      DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQELikeSelect).c_str(), DoLogx));

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

  void SetTopologyDefinedECal(){

    std::string GotMuon = "(MuonProtonTransversity.MuonPDG==13)";
    std::string OneAboveThreshChargedTrack =
    "((MuonProtonTransversity.NAboveThresholdProtons"
        "+ MuonProtonTransversity.NAboveThresholdChargedPions) == 1)";

    Selections1D["ECalTop2Prong_CCQE"] =
      (Selection1D("ECalTop2Prong_CCQE","ECalTop2Prong_CCQE",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==1)&&"
        + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong_CCRes"] =
      (Selection1D("ECalTop2Prong_CCRes","ECalTop2Prong_CCRes",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode>=11)&&"
      "(MuonProtonTransversity.NeutConventionReactionCode<=13)&&"
        + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong_CCMEC"] =
      (Selection1D("ECalTop2Prong_CCMEC","ECalTop2Prong_CCMEC",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==2)&&"
        + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong"] =
      (Selection1D("ECalTop2Prong","ECalTop2Prong",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&" + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECal2ProngHMTrackableMom"] =
      (Selection1D("ECal2ProngHMTrackableMom","ECal2ProngHMTrackableMom",
      "MuonProtonTransversity.HMTrackableMomentum_MeV",
      ProtonMomBins, 0, ProtonMomMax,
      (GotMuon + "&&" + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong_CCQE_sm"] =
      (Selection1D("ECalTop2Prong_CCQE_sm","ECalTop2Prong_CCQE_sm",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==1)&&"
        + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong_CCRes_sm"] =
      (Selection1D("ECalTop2Prong_CCRes_sm","ECalTop2Prong_CCRes_sm",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode>=11)&&"
      "(MuonProtonTransversity.NeutConventionReactionCode<=13)&&"
        + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong_CCMEC_sm"] =
      (Selection1D("ECalTop2Prong_CCMEC_sm","ECalTop2Prong_CCMEC_sm",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==2)&&"
        + OneAboveThreshChargedTrack).c_str(), DoLogx));

    Selections1D["ECalTop2Prong_sm"] =
      (Selection1D("ECalTop2Prong_sm","ECalTop2Prong_sm",
      "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
      DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
      (GotMuon + "&&" + OneAboveThreshChargedTrack).c_str(), DoLogx));
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

    Selections1D["ReconTargetMassselres"] =
      Selection1D("ReconTargetMass","ReconTargetMassselres",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + selres).c_str(),DoLogx);
    Selections1D["ReconTargetMassselcoh"] =
      Selection1D("ReconTargetMass","ReconTargetMassselcoh",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + selcoh).c_str(),DoLogx);
    Selections1D["ReconTargetMassselccother"] =
      Selection1D("ReconTargetMass","ReconTargetMassselccother",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + selccother).c_str(),DoLogx);
    Selections1D["ReconTargetMassseldis"] =
      Selection1D("ReconTargetMass","ReconTargetMassseldis",
        "MuonProtonTransversity.ReconTargetMass",
        100,0,2000,
        (MuonValidSelect
          + "&&" + FirstProtonValidSelect
          + "&&" + FirstProtonNonZeroSelect
          + seldis).c_str(),DoLogx);
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
    Selections1D["CCQEDPt"] =
      (Selection1D("CCQE DPt","CCQEDPt",
      "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
      DPTBins, DPTMin, DPTMax,
      (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

    Selections1D["CCQENEDPt"] =
      (Selection1D("CCQE Nuclear Emission DPt","CCQENEDPt",
      "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
      DPTBins, DPTMin, DPTMax,
      (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

    Selections1D["CCQENoNEDPt"] =
      (Selection1D("CCQE No Nuclear Emission DPt","CCQENoNEDPt",
      "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
      DPTBins, DPTMin, DPTMax,
      (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//AlphaT
    Selections1D["CCQEAlphaT"] =
      (Selection1D("CCQE AlphaT","CCQEAlphaT",
      "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
      DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
      (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

    Selections1D["CCQENEAlphaT"] =
      (Selection1D("CCQE Nuclear Emission AlphaT","CCQENEAlphaT",
      "MuonProtonTransversity.",
      DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

    Selections1D["CCQENoNEAlphaT"] =
      (Selection1D("No Nuclear Emission","CCQENoNEAlphaT",
      "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
      DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
      (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//2D Distributions
    Selections2D["DeltaPhiTHMPvsMuonPt"] =
      (Selection2D("DeltaPhiTHMPvsMuonPt","DeltaPhiTHMPvsMuonPt",
        "MuonProtonTransversity.MuonPt_MeV.Mag():"
        "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
        DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
        100,0,1000,(CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

    Selections2D["DeltaAlphaTHMPvsMuonPt"] =
      (Selection2D("DeltaAlphaTHMPvsMuonPt","DeltaAlphaTHMPvsMuonPt",
        "MuonProtonTransversity.MuonPt_MeV.Mag():"
        "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
        DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
        100,0,1000,(CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

    Selections2D["DeltaPTHMPvsMuonPt"] =
      (Selection2D("DeltaPTHMPvsMuonPt","DeltaPTHMPvsMuonPt",
        "MuonProtonTransversity.MuonPt_MeV.Mag():"
        "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
        DPTBins, DPTMin, DPTMax,
        100,0,1000,(CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

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

  void InitSelections(){
    (void)Verbosity;
    InitCutStrings();
    SetProtonMomentaSelections();
    SetModeSelections();
    SetTopologyDefinedECal();
    TargetMassSelection();
    SetCCQESelections();
    FinalStateSpeciesSelections("MuonProtonTransversity",CCQESelect,"ccqe");
    DetailsPiProdSelection();
    BasicPiProdSelections();

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
