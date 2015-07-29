#include <map>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdexcept>

#include "TH1.h"
#include "TH2.h"

#include "PureGenUtils.hxx"

#include "PlottingSelections.hxx"

using namespace PlottingUtils;
using namespace PlottingTypes;
using namespace PlottingDefaults;

namespace PlottingSelections {

  std::map<std::string, PlottingTypes::Selection1D> Selections1D;
  std::map<std::string, PlottingTypes::Selection2D> Selections2D;

// //**********************************CCQEFSI*************************************
  std::string MuonValidSelect;
  std::string FirstProtonValidSelect;
  std::string FirstProtonNonZeroSelect;
  std::string CCQERCSelect;
  std::string CCQEDPhiTSelect;
  std::string CCQENoCCQENucEmitSelect;
  std::string CCQENucEmitSelect;
  std::string CCQESelect;
  std::string MECSelect;
  std::string CCQE_MECSelect;
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

  std::string CC0Pi1P_Select_thresh;
  std::string CC0Pi1P_Select;
// //**********************************Globals*************************************

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
    CCQEDPhiTSelect = "(MuonProtonTransversity.HMProtonMomentum_MeV > "
        + MomentumCut_MeV +")";

    CCQENoCCQENucEmitSelect = MuonValidSelect + "&&" + FirstProtonValidSelect +
      "&&(MuonProtonTransversity.NFinalStateParticles==2)&&" + CCQERCSelect;

    CCQENucEmitSelect = "(!("+CCQENoCCQENucEmitSelect+"))&&"+CCQERCSelect;

    CCQENucEmitSelect = MuonValidSelect + "&&" + FirstProtonValidSelect +
      "&&(MuonProtonTransversity.NFinalStateParticles!=2)&&" + CCQERCSelect;

    CCQESelect = MuonValidSelect+"&&"+FirstProtonValidSelect+"&&"+CCQERCSelect;

    MECSelect = "(MuonProtonTransversity.NeutConventionReactionCode==2)&&"
      + MuonValidSelect + "&&" + FirstProtonValidSelect + "&&" + FirstProtonNonZeroSelect;

    CC0Pi1P_Select_thresh = "((MuonProtonTransversity.NAboveThresholdProtons>0)&&(MuonProtonTransversity.NAboveThresholdPiPlus==0)&&(MuonProtonTransversity.NAboveThresholdPiMinus==0))"
      "&&" + MuonValidSelect + "&&" + FirstProtonValidSelect;
    CC0Pi1P_Select = "((MuonProtonTransversity.NPiPlus==0)&&(MuonProtonTransversity.NPiMinus==0))"
      "&&" + MuonValidSelect + "&&" + FirstProtonValidSelect + "&&" + FirstProtonNonZeroSelect;

    CCQE_MECSelect = "(" + CCQERCSelect
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

//   void SetProtonMomentaSelections(){
//     Selections1D["CCQEMuonMom"] =
//       (Selection1D("CCQEMuonMom","CCQEMuonMom",
//       "PionProductionTransversity.MuonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//         CCQESelect.c_str(),
//         DoLogx));

//     Selections1D["HMProtonStdHepPos"] =
//       (Selection1D("HMProtonStdHepPos","HMProtonStdHepPos",
//       "MuonProtonTransversity.HMProtonStdHepPosition",
//       15, 0, 15,
//       (CCQESelect + "&&" + FirstProtonNonZeroSelect).c_str(), DoLogx));

//     Selections1D["HMProtonStdHepPosNE"] =
//       (Selection1D("HMProtonStdHepPosNE","HMProtonStdHepPosNE",
//       "MuonProtonTransversity.HMProtonStdHepPosition",
//       15, 0, 15,
//       (CCQESelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["HMProtonStdHepPos_excludediagonal"] =
//       (Selection1D("HMProtonStdHepPos_excludediagonal","HMProtonStdHepPos_excludediagonal",
//       "MuonProtonTransversity.HMProtonStdHepPosition",
//       20, 0, 20,
//      (CCQESelect + "&&" + CCQENucEmitSelect +
//       "&&(abs(MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()-"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag())>10)").c_str(),
//      DoLogx));

//     Selections1D["CCQEHMProtonMom"] =
//       (Selection1D("CCQEHMProtonMom","CCQEHMProtonMom",
//       "MuonProtonTransversity.HMProtonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//         CCQESelect.c_str(),
//         DoLogx));

//     Selections1D["CCQEHMProtonMom_NE"] =
//       (Selection1D("CCQEHMProtonMom_NE","CCQEHMProtonMom",
//       "MuonProtonTransversity.HMProtonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//       (CCQESelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQEHMProtonMom_excludediagonal"] =
//       (Selection1D("CCQEHMProtonMom_excludediagonal","CCQEHMProtonMom_excludediagonal",
//       "MuonProtonTransversity.HMProtonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//       (CCQESelect + "&&" + CCQENucEmitSelect +
//       "&&(abs(MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()-"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag())>10)").c_str(), DoLogx));

//     Selections1D["CCQEMuon_Pt"] =
//       (Selection1D("CCQEMuon_Pt","CCQEMuon_Pt",
//       "MuonProtonTransversity.MuonPt_MeV.Mag()",
//       ProtonMomBins, 0, ProtonMomMax,
//       CCQESelect.c_str(),
//       DoLogx));

//     Selections1D["CCQEHMProton_Pt"] =
//       (Selection1D("CCQEHMProton_Pt","CCQEHMProton_Pt",
//       "MuonProtonTransversity.HMProtonPt_MeV.Mag()",
//       ProtonMomBins, 0, ProtonMomMax,
//         CCQESelect.c_str(),
//         DoLogx));

//     Selections1D["CCQEHMProton_Phi"] =
//       (Selection1D("CCQEHMProton_Phi","CCQEHMProton_Phi",
//       "MuonProtonTransversity.HMProtonDirection.Phi()",
//       100, -1.0*M_PI, M_PI,
//       (CCQESelect + "&&" + FirstProtonNonZeroSelect).c_str(),
//       DoLogx));

//     Selections1D["CCQEHMProton_CTheta"] =
//       (Selection1D("CCQEHMProton_CTheta","CCQEHMProton_CTheta",
//       "MuonProtonTransversity.HMProtonDirection.CosTheta()",
//       100, -1, 1,
//       (CCQESelect + "&&" + FirstProtonNonZeroSelect).c_str(),
//       DoLogx));
//   }

//   void SetModeSelections(){

//     Selections1D["MEC"] =
//       (Selection1D("MEC","MEC",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + MECSelect).c_str(), DoLogx));

//     Selections1D["MECZoom"] =
//       (Selection1D("MECZoom","MECZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + MECSelect).c_str(), DoLogx));

//     Selections1D["CCQELike"] =
//       (Selection1D("CCQELike","CCQELike",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQE_MECSelect).c_str(), DoLogx));

//     Selections1D["CCQELikeZoom"] =
//       (Selection1D("CCQELikeZoom","CCQELikeZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQE_MECSelect).c_str(), DoLogx));

//     Selections1D["CC_not_QELike"] =
//       (Selection1D("CC_not_QELike","CC_not_QELike",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CC_not_QELikeSelect).c_str(), DoLogx));

//     Selections1D["CC_not_QELikeZoom"] =
//       (Selection1D("CC_not_QELikeZoom","CC_not_QELikeZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQE_MECSelect).c_str(), DoLogx));

//     Selections1D["CCRes"] =
//       (Selection1D("CCRes","CCRes",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCResSelect).c_str(), DoLogx));

//     Selections1D["CCResZoom"] =
//       (Selection1D("CCResZoom","CCResZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCResSelect).c_str(), DoLogx));

//     Selections1D["CCOther"] =
//       (Selection1D("CCOther","CCOther",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCOtherSelect).c_str(), DoLogx));

//     Selections1D["CCOtherZoom"] =
//       (Selection1D("CCOtherZoom","CCOtherZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCOtherSelect).c_str(), DoLogx));
//   }

//   void SetTopologyDefinedECal(){

//     std::string GotMuon = "(MuonProtonTransversity.MuonPDG==13)";
//     std::string HMTrackableAboveThresh = "(HMTrackableMomentum_MeV > 400.0)";
//     std::string Prong2 = "&&(NAboveThresholdTrackable == 1)";

//     Selections1D["ECalTop2Prong_CCQE"] =
//       (Selection1D("ECalTop2Prong_CCQE","ECalTop2Prong_CCQE",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==1)&&"
//         + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong_CCRes"] =
//       (Selection1D("ECalTop2Prong_CCRes","ECalTop2Prong_CCRes",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode>=11)&&"
//       "(MuonProtonTransversity.NeutConventionReactionCode<=13)&&"
//         + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong_CCMEC"] =
//       (Selection1D("ECalTop2Prong_CCMEC","ECalTop2Prong_CCMEC",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==2)&&"
//         + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong"] =
//       (Selection1D("ECalTop2Prong","ECalTop2Prong",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&" + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECal2ProngHMTrackableMom"] =
//       (Selection1D("ECal2ProngHMTrackableMom","ECal2ProngHMTrackableMom",
//       "MuonProtonTransversity.HMTrackableMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//       (GotMuon + "&&" + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong_CCQE_sm"] =
//       (Selection1D("ECalTop2Prong_CCQE_sm","ECalTop2Prong_CCQE_sm",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==1)&&"
//         + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong_CCRes_sm"] =
//       (Selection1D("ECalTop2Prong_CCRes_sm","ECalTop2Prong_CCRes_sm",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode>=11)&&"
//       "(MuonProtonTransversity.NeutConventionReactionCode<=13)&&"
//         + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong_CCMEC_sm"] =
//       (Selection1D("ECalTop2Prong_CCMEC_sm","ECalTop2Prong_CCMEC_sm",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&(MuonProtonTransversity.NeutConventionReactionCode==2)&&"
//         + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections1D["ECalTop2Prong_sm"] =
//       (Selection1D("ECalTop2Prong_sm","ECalTop2Prong_sm",
//       "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg",
//       30, DPhiTMin_deg, DPhiTMax_deg,
//       (GotMuon + "&&" + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));

//     Selections2D["ECalTop2Prong_sm_DeltaPhiTHMPvsMuonPt"] =
//       (Selection2D("ECalTop2Prong_sm_DeltaPhiTHMPvsMuonPt","ECalTop2Prong_sm_DeltaPhiTHMPvsMuonPt",
//         "MuonProtonTransversity.DeltaPhiT_HMTrackable_GSmear_deg:"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         8,200,600,
//         30, DPhiTMin_deg, DPhiTMax_deg,
//         (GotMuon + "&&" + HMTrackableAboveThresh + Prong2).c_str(), DoLogx));
//   }

//   void TargetMassSelection(){

//     static std::string const selccqe =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode==1)";
//     static std::string const selccqeNE =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode==1)"
//       "&&((MuonProtonTransversity.FirstProtonPDG!=2212)||"
//         "(MuonProtonTransversity.NFinalStateParticles!=2))";
//     static std::string const selmec =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode==2)";
//     static std::string const selres =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode>=11)"
//       "&&(MuonProtonTransversity.NeutConventionReactionCode<=13)";
//     static std::string const selcoh =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode==16)";
//     static std::string const selccother =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode>=17)"
//       "&&(MuonProtonTransversity.NeutConventionReactionCode<=23)";
//     static std::string const seldis =
//       "&&(MuonProtonTransversity.NeutConventionReactionCode==26)";

//     Selections1D["ReconTargetMassselccqe"] =
//       Selection1D("ReconTargetMass","ReconTargetMassselccqe",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + selccqe).c_str(),DoLogx);
//     Selections1D["ReconTargetMassselccqeNE"] =
//       Selection1D("ReconTargetMassNE","ReconTargetMassselccqeNE",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + selccqeNE).c_str(),DoLogx);
//     Selections1D["ReconTargetMassselmec"] =
//       Selection1D("ReconTargetMass","ReconTargetMassselmec",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + selmec).c_str(),DoLogx);

//     Selections1D["ReconTargetMassselres"] =
//       Selection1D("ReconTargetMass","ReconTargetMassselres",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + selres).c_str(),DoLogx);
//     Selections1D["ReconTargetMassselcoh"] =
//       Selection1D("ReconTargetMass","ReconTargetMassselcoh",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + selcoh).c_str(),DoLogx);
//     Selections1D["ReconTargetMassselccother"] =
//       Selection1D("ReconTargetMass","ReconTargetMassselccother",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + selccother).c_str(),DoLogx);
//     Selections1D["ReconTargetMassseldis"] =
//       Selection1D("ReconTargetMass","ReconTargetMassseldis",
//         "MuonProtonTransversity.ReconTargetMass",
//         100,0,2000,
//         (MuonValidSelect
//           + "&&" + FirstProtonValidSelect
//           + "&&" + FirstProtonNonZeroSelect
//           + seldis).c_str(),DoLogx);
//   }

//   void SetCCQESelections(){

//     Selections1D["CCQE"] =
//       (Selection1D("CCQE","CCQE",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQE_MEC"] =
//       (Selection1D("CCQE_MEC","CCQE_MEC",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQE_MECSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQECoplanarity"] =
//       (Selection1D("CCQECoplanarity","CCQECoplanarity",
//       "180.0 - MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENoGamma"] =
//       (Selection1D("CCQENoGamma","CCQENoGamma",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQESelect + "&&" + CCQEDPhiTSelect
//         + "&&(MuonProtonTransversity.NGamma==0)").c_str(), DoLogx));

//     Selections1D["CCQEZoom"] =
//       (Selection1D("CCQE All Zoom","CCQEZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQECoplanarityZoom"] =
//       (Selection1D("CCQECoplanarityZoom","CCQECoplanarityZoom",
//       "180.0 - MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, 174, 180,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENE"] =
//       (Selection1D("CCQE Nuclear Emission","CCQENE",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["NEZoom"] =
//       (Selection1D("Nuclear Emission Zoom","NEZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQENoNE"] =
//       (Selection1D("CCQE No Nuclear Emission","CCQENoNE",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQENoNEZoom"] =
//       (Selection1D("CCQE No Nuclear Emission Zoom","CCQENoNEZoom",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

// //DeltaPT
//     Selections1D["CCQEDPt"] =
//       (Selection1D("CCQE DPt","CCQEDPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQEDPt_MEC"] =
//       (Selection1D("CCQE_MEC DPt","CCQEDPt_MEC",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQE_MECSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENEDPt"] =
//       (Selection1D("CCQE Nuclear Emission DPt","CCQENEDPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQENoNEDPt"] =
//       (Selection1D("CCQE No Nuclear Emission DPt","CCQENoNEDPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

// //MUPt
//     Selections1D["CCQEMuPt"] =
//       (Selection1D("CCQE MuPt","CCQEMuPt",
//       "MuonProtonTransversity.MuonPt_MeV.Mag()",
//       500,0,700,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENEMuPt"] =
//       (Selection1D("CCQE Nuclear Emission MuPt","CCQENEMuPt",
//       "MuonProtonTransversity.MuonPt_MeV.Mag()",
//       500,0,700,
//       (CCQEDPhiTSelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections1D["CCQENoNEMuPt"] =
//       (Selection1D("CCQE No Nuclear Emission MuPt","CCQENoNEMuPt",
//       "MuonProtonTransversity.MuonPt_MeV.Mag()",
//       500,0,700,
//       (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

// //StrucNucleonPt
//     Selections1D["CCQEStruckNucleonZMomentum"] =
//       (Selection1D("CCQEStruckNucleonZMomentum","CCQEStruckNucleonZMomentum",
//       "MuonProtonTransversity.StruckNucleonDirection.fZ*"
//       "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//       101,-505,505,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQEStruckNucleonMomentum"] =
//       (Selection1D("CCQEStruckNucleonMomentum","CCQEStruckNucleonMomentum",
//       "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//       100,0,500,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQEStruckNucleonMomentum_MEC"] =
//       (Selection1D("CCQEStruckNucleonMomentum_MEC","CCQEStruckNucleonMomentum_MEC",
//       "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//       100,0,500,
//       (CCQE_MECSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQEStruckNucleon_BE"] =
//       (Selection1D("CCQEStruckNucleon_BE","CCQEStruckNucleon_BE",
//       "MuonProtonTransversity.StruckNucleon_BindingEnergy_MeV",
//       100,0,200,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx, true, 0.05));

//     Selections2D["CCQEStruckNucleon_BEvsMomentum"] =
//       (Selection2D("CCQEStruckNucleon_BEvsMomentum",
//         "CCQEStruckNucleon_BEvsMomentum",
//       "MuonProtonTransversity.StruckNucleon_BindingEnergy_MeV:"
//       "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//       50,0,100, 100,0,500,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

// //AlphaT
//     Selections1D["CCQEAlphaT"] =
//       (Selection1D("CCQE AlphaT","CCQEAlphaT",
//       "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQEAlphaT_MEC"] =
//       (Selection1D("CCQE_MEC AlphaT","CCQEAlphaT_MEC",
//       "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CCQE_MECSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

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

// //2D Distributions
//     Selections2D["DeltaPhiTHMPvsMuonPt"] =
//       (Selection2D("DeltaPhiTHMPvsMuonPt","DeltaPhiTHMPvsMuonPt",
//         "MuonProtonTransversity.DeltaPhiT_HMProton_deg:"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         500,0,1000,
//         DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["AlphaTHMPvsMuonPt"] =
//       (Selection2D("AlphaTHMPvsMuonPt","AlphaTHMPvsMuonPt",
//         "MuonProtonTransversity.DeltaAlphaT_HMProton_deg:"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         500,0,1000,
//         DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsMuonPt"] =
//       (Selection2D("DPTHMPvsMuonPt","DPTHMPvsMuonPt",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         500,0,1000,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DeltaPhiTHMPvsMuonPt_NoNE"] =
//       (Selection2D("DeltaPhiTHMPvsMuonPt_NoNE","DeltaPhiTHMPvsMuonPt_NoNE",
//         "MuonProtonTransversity.DeltaPhiT_HMProton_deg:"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         500,0,1000,
//         DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//         (CCQENoCCQENucEmitSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["AlphaTHMPvsMuonPt_NoNE"] =
//       (Selection2D("AlphaTHMPvsMuonPt_NoNE","AlphaTHMPvsMuonPt_NoNE",
//         "MuonProtonTransversity.DeltaAlphaT_HMProton_deg:"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         500,0,1000,
//         DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//         (CCQENoCCQENucEmitSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsMuonPt_NoNE"] =
//       (Selection2D("DPTHMPvsMuonPt_NoNE","DPTHMPvsMuonPt_NoNE",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.MuonPt_MeV.Mag()",
//         500,0,1000,
//         DPTBins, DPTMin, DPTMax,
//         (CCQENoCCQENucEmitSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsDAlphaT"] =
//       (Selection2D("DPTHMPvsDAlphaT","DPTHMPvsDAlphaT",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//         DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//       //en - emu
//     Selections2D["OmegaHMPvsDAlphaT"] =
//       (Selection2D("OmegaHMPvsDAlphaT","OmegaHMPvsDAlphaT",
//         "MuonProtonTransversity.IncNeutrinoMmtm.E() - MuonProtonTransversity.MuonMmtm.E():"
//         "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//         DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//         100, 0, 700,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPhiTHMPvsDAlphaT"] =
//       (Selection2D("DPhiTHMPvsDAlphaT","DPhiTHMPvsDAlphaT",
//         "MuonProtonTransversity.DeltaPhiT_HMProton_deg:"
//         "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//         DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//         DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DeltaPhiTHMPvsStruckNucleonZMomentum"] =
//       (Selection2D("DeltaPhiTHMPvsStruckNucleonZMomentum",
//         "DeltaPhiTHMPvsStruckNucleonZMomentum",
//         "MuonProtonTransversity.DeltaPhiT_HMProton_deg:"
//         "MuonProtonTransversity.StruckNucleonDirection.fZ*"
//         "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//         101,-505,505,
//         DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["AlphaTHMPvsStruckNucleonZMomentum"] =
//       (Selection2D("AlphaTHMPvsStruckNucleonZMomentum",
//         "AlphaTHMPvsStruckNucleonZMomentum",
//         "MuonProtonTransversity.DeltaAlphaT_HMProton_deg:"
//         "MuonProtonTransversity.StruckNucleonDirection.fZ*"
//         "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//         101,-505,505,
//         DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsStruckNucleonZMomentum"] =
//       (Selection2D("DPTHMPvsStruckNucleonZMomentum",
//         "DPTHMPvsStruckNucleonZMomentum",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.StruckNucleonDirection.fZ*"
//         "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//         101,-505,505,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["HMPMomVsStruckNucleonMom"] =
//       (Selection2D("HMPMomVsStruckNucleonMom",
//         "HMPMomVsStruckNucleonMom",
//         "MuonProtonTransversity.HMProtonMomentum_MeV:"
//         "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//         100, 0, 1000,
//         100, 0, 1000,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

// // MuonDirection
// // MuonMomentum_MeV
// // HMProtonDirection
// // HMProtonMomentum_MeV
// // MuonMmtm
// // StruckNucleonMomentum_MeV

//     Selections1D["LabMom"] =
//       (Selection1D("LabMom",
//         "LabMom",
//         // "MuonProtonTransversity.LabCCQEMomentum.Mag()",
//         "sqrt(MuonProtonTransversity.LabCCQEMomentum.fX*MuonProtonTransversity.LabCCQEMomentum.fX +"
//         "MuonProtonTransversity.LabCCQEMomentum.fY*MuonProtonTransversity.LabCCQEMomentum.fY +"
//         "MuonProtonTransversity.LabCCQEMomentum.fZ*MuonProtonTransversity.LabCCQEMomentum.fZ)",
//         100, 0, 1000,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["LabMomVsStruckNucleonMom"] =
//       (Selection2D("LabMomVsStruckNucleonMom",
//         "LabMomVsStruckNucleonMom",
//         "sqrt(MuonProtonTransversity.LabCCQEMomentum.fX*MuonProtonTransversity.LabCCQEMomentum.fX +"
//         "MuonProtonTransversity.LabCCQEMomentum.fY*MuonProtonTransversity.LabCCQEMomentum.fY +"
//         "MuonProtonTransversity.LabCCQEMomentum.fZ*MuonProtonTransversity.LabCCQEMomentum.fZ)"
//         // "MuonProtonTransversity.LabCCQEMomentum.Mag()"
//         ":"
//         "MuonProtonTransversity.StruckNucleonMomentum_MeV",
//         100, 0, 1000,
//         100, 0, 1000,
//         (CCQESelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsStruckNucleonPt"] =
//       (Selection2D("DPTHMPvsStruckNucleonPt",
//         "DPTHMPvsStruckNucleonPt",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag()",
//         100,0,700,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPTFirstPvsStruckNucleonPt"] =
//       (Selection2D("DPTFirstPvsStruckNucleonPt",
//         "DPTFirstPvsStruckNucleonPt",
//         "MuonProtonTransversity.DeltaPT_FirstProton_MeV.Mag():"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag()",
//         100,0,700,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsStruckNucleonPtNoNE"] =
//       (Selection2D("DPTHMPvsStruckNucleonPtNoNE",
//         "DPTHMPvsStruckNucleonPtNoNE",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag()",
//         100,0,700,
//         DPTBins, DPTMin, DPTMax,
//         (CCQEDPhiTSelect + "&&" + CCQENoCCQENucEmitSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsStruckNucleonPtNE"] =
//       (Selection2D("DPTHMPvsStruckNucleonPtNE",
//         "DPTHMPvsStruckNucleonPtNE",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag()",
//         100,0,700,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQENucEmitSelect).c_str(), DoLogx));

//     Selections2D["DPTHMPvsStruckNucleonPtNE_excludediagonal"] =
//       (Selection2D("DPTHMPvsStruckNucleonPtNE_excludediagonal",
//         "DPTHMPvsStruckNucleonPtNE_excludediagonal",
//         "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag():"
//         "MuonProtonTransversity.StruckNucleonPt_MeV.Mag()",
//         100,0,700,
//         DPTBins, DPTMin, DPTMax,
//         (CCQESelect + "&&" + CCQENucEmitSelect +
//           "&&(abs(MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()-"
//             "MuonProtonTransversity.StruckNucleonPt_MeV.Mag())>10)").c_str(),
//         DoLogx));
//   }

//   void SetERecSelections(){
// //ERec
//     Selections1D["CCQENuERec"] =
//       (Selection1D("CCQENuERec","CCQENuERec",
//       "MuonProtonTransversity.ERec_CCQE",
//       5000, 0, 5000,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENuMuPERec"] =
//       (Selection1D("CCQENuMuPERec","CCQENuMuPERec",
//       "MuonProtonTransversity.ERec_E_PMu",
//       5000, 0, 5000,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENuETrue"] =
//       (Selection1D("CCQENuETrue","CCQENuETrue",
//       "MuonProtonTransversity.IncNeutrinoMmtm.E()",
//       5000, 0, 5000,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENuERec_CCQENuETrue"] =
//       (Selection1D("CCQENuERec_CCQENuETrue","CCQENuERec_CCQENuETrue",
//       "MuonProtonTransversity.ERec_CCQE/MuonProtonTransversity.IncNeutrinoMmtm.E()",
//       1000, 0, 10,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQENuMuPERec_CCQENuETrue"] =
//       (Selection1D("CCQENuMuPERec_CCQENuETrue","CCQENuMuPERec_CCQENuETrue",
//       "MuonProtonTransversity.ERec_E_PMu/MuonProtonTransversity.IncNeutrinoMmtm.E()",
//       1000, 0, 10,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQE_mmtm_mu"] =
//       (Selection1D("CCQE_mmtm_mu","CCQE_mmtm_mu",
//       "MuonProtonTransversity.MuonMomentum_MeV",
//       500, 0, 5000,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CCQE_mmtm_p"] =
//       (Selection1D("CCQE_mmtm_p","CCQE_mmtm_p",
//       "MuonProtonTransversity.HMProtonMomentum_MeV",
//       500, 0, 5000,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections2D["CCQE_fsmmtm"] =
//       (Selection2D("CCQE_fsmmtm",
//         "CCQE_fsmmtm",
//         "MuonProtonTransversity.HMProtonMomentum_MeV:"
//         "MuonProtonTransversity.MuonMomentum_MeV",
//         100,0,1000,
//         100, 0, 1000,
//         (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//   }

//   void SetCC0PiSelections(){

//     Selections1D["CCQE_MEC_StruckNucleon_BE"] =
//       (Selection1D("CCQE_MEC_StruckNucleon_BE","CCQE_MEC_StruckNucleon_BE",
//       "MuonProtonTransversity.StruckNucleon_BindingEnergy_MeV",
//       100,0,200,
//       (CCQE_MECSelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx));

//     Selections1D["CC0Pi1P"] =
//       (Selection1D("CC0Pi1P","CC0Pi1P",
//       "MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CC0Pi1P_Select).c_str(), DoLogx));

//     Selections1D["CC0Pi1PDPt"] =
//       (Selection1D("CC0Pi1P DPt","CC0Pi1PDPt",
//       "MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CC0Pi1P_Select).c_str(), DoLogx));

//     Selections1D["CC0Pi1PAlphaT"] =
//       (Selection1D("CC0Pi1P AlphaT","CC0Pi1PAlphaT",
//       "MuonProtonTransversity.DeltaAlphaT_HMProton_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CC0Pi1P_Select).c_str(), DoLogx));
//   }

//   void FinalStateSpeciesSelections(std::string drawObject,
//     std::string selectString,
//     std::string prefix=""){

//     constexpr int NProtonsMax = 50;
//     constexpr int NGammaMax = 5;
//     constexpr int NNeutronsMax = 50;
//     constexpr int NPiPlusMax = 5;
//     constexpr int NPiZeroMax = 5;
//     constexpr int NPiMinusMax = 5;
//     constexpr int NFinalStateParticles = 100;

//     Selections1D[prefix + "NProtons"] =
//       (Selection1D("NProtons", prefix + "NProtons",
//       drawObject + ".NProtons",
//       NProtonsMax, 0, NProtonsMax, selectString.c_str(), DoLogx, true));

//     Selections1D[prefix + "NGamma"] =
//       (Selection1D("NGamma", prefix + "NGamma",
//       drawObject + ".NGamma",
//       NGammaMax, 0, NGammaMax, selectString.c_str(), DoLogx, true));

//     Selections1D[prefix + "NNeutron"] =
//       (Selection1D("NNeutrons", prefix + "NNeutron",
//       drawObject + ".NNeutrons",
//       NNeutronsMax, 0, NNeutronsMax, selectString.c_str(), DoLogx, true));

//     Selections1D[prefix + "NPiPlus"] =
//       (Selection1D("NPiPlus", prefix + "NPiPlus",
//       drawObject + ".NPiPlus",
//       NPiPlusMax, 0, NPiPlusMax, selectString.c_str(), DoLogx, true));

//     Selections1D[prefix + "NPiZero"] =
//       (Selection1D("NPiZero", prefix + "NPiZero",
//       drawObject + ".NPiZero",
//       NPiZeroMax, 0, NPiZeroMax, selectString.c_str(), DoLogx, true));

//     Selections1D[prefix + "NPiMinus"] =
//       (Selection1D("NPiMinus", prefix + "NPiMinus",
//       drawObject + ".NPiMinus",
//       NPiMinusMax, 0, NPiMinusMax, selectString.c_str(), DoLogx, true));

//     Selections1D[prefix + "NPrimaryParticles"] =
//       (Selection1D("NFinalStateParticles", prefix + "NPrimaryParticles",
//       drawObject + ".NFinalStateParticles",
//       NFinalStateParticles, 0, NFinalStateParticles, selectString.c_str(),
//       DoLogx, true));
//   }

//   void DetailsPiProdSelection(){
//     Selections1D["Res11"] =
//       (Selection1D("Res11","Res11",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["Res11Zoom"] =
//       (Selection1D("Res11Zoom","Res11Zoom",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["Res13"] =
//       (Selection1D("Res13","Res13",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["Res13Zoom"] =
//       (Selection1D("Res13Zoom","Res13Zoom",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (ResPPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["PiProdOther"] =
//       (Selection1D("PiProdOther","PiProdOther",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (OtherPiPlusSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["PiProdOtherZoom"] =
//       (Selection1D("PiProdOtherZoom","PiProdOtherZoom",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (OtherPiPlusSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

// //COH Momenta START
//     Selections1D["CohMuonMom"] =
//       (Selection1D("CohMuonMom","CohMuonMom",
//       "PionProductionTransversity.MuonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//         "(PionProductionTransversity.NeutConventionReactionCode==16)",
//         DoLogx));

//     Selections1D["CohHMPiPlusMom"] =
//       (Selection1D("CohHMPiPlusMom","CohHMPiPlusMom",
//       "PionProductionTransversity.HMPiPlusMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//         "(PionProductionTransversity.NeutConventionReactionCode==16)",
//         DoLogx));

//     Selections1D["CohMuon_Pt"] =
//       (Selection1D("CohMuon_Pt","CohMuon_Pt",
//       "PionProductionTransversity.MuonPt_MeV.Mag()",
//       ProtonMomBins, 0, ProtonMomMax,
//         "(PionProductionTransversity.NeutConventionReactionCode==16)",
//         DoLogx));

//     Selections1D["CohHMPiPlus_Pt"] =
//       (Selection1D("CohHMPiPlus_Pt","CohHMPiPlus_Pt",
//       "PionProductionTransversity.HMPiPlusPt_MeV.Mag()",
//       ProtonMomBins, 0, ProtonMomMax,
//         "(PionProductionTransversity.NeutConventionReactionCode==16)",
//         DoLogx));
//   }

//   void BasicPiProdSelections(){
//     Selections1D["CCPiProdInc"] =
//       (Selection1D("CCPiProdInc","CCPiProdInc",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CCSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["CCPiProdIncZoom"] =
//       (Selection1D("CCPiProdIncZoom","CCPiProdIncZoom",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CCSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["CohPiProd"] =
//       (Selection1D("CohPiProd","CohPiProd",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["CohPiProdZoom"] =
//       (Selection1D("Coherent Pion Production Zoom","CohPiProdZoom",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["NonCoh"] =
//       (Selection1D("NonCoh","NonCoh",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
//         + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str()
//       , DoLogx));

//     Selections1D["NonCohZoom"] =
//       (Selection1D("NonCohZoom","NonCohZoom",
//       "PionProductionTransversity.DeltaPhiT_HMPiPlus_deg",
//       DPhiTZoomxBins_deg, DPhiTZoomxMin_deg, DPhiTZoomxMax_deg,
//       ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
//         + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str(),
//       DoLogx));

//     Selections1D["CohPiProdPt"] =
//       (Selection1D("Coh Pt","CohPiProdPt",
//       "PionProductionTransversity.DeltaPT_HMPiPlus_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["CohPiProdAlphaT"] =
//       (Selection1D("Coh AlphaT","CohPiProdAlphaT",
//       "PionProductionTransversity.DeltaAlphaT_HMPiPlus_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       (CohPiPlusProdSelect + "&&" + FirstPiPlusValidSelect + "&&"
//         + FirstPiPlusNonZeroSelect).c_str(), DoLogx));

//     Selections1D["NonCohPt"] =
//       (Selection1D("NonCoh Pt","NonCohPt",
//       "PionProductionTransversity.DeltaPT_HMPiPlus_MeV.Mag()",
//       DPTBins, DPTMin, DPTMax,
//       ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
//         + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str(),
//       DoLogx));

//     Selections1D["NonCohAlphaT"] =
//       (Selection1D("NonCoh AlphaT","NonCohAlphaT",
//       "PionProductionTransversity.DeltaAlphaT_HMPiPlus_deg",
//       DAlphaTBins_deg, DAlphaTMin_deg, DAlphaTMax_deg,
//       ("(!"+CohPiPlusProdSelect + ")&&" + CCSelect + "&&"
//         + FirstPiPlusValidSelect + "&&" + FirstPiPlusNonZeroSelect).c_str(),
//       DoLogx));

//   }

//   void Multiplicity2DSels(){

//     Selections2D["NProtonsVsCCQEDphiT"] =
//       (Selection2D("NProtonsVsCCQEDphiT","NProtonsVsCCQEDphiT",
//       "MuonProtonTransversity.NProtons:MuonProtonTransversity.DeltaPhiT_HMProton_deg",
//       DPhiTBins_deg, DPhiTMin_deg, DPhiTMax_deg,
//       10, 0, 10,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx,
//       true));

//     Selections2D["NProtonsVsCCQEMuPt"] =
//       (Selection2D("NProtonsVsCCQEMuPt","NProtonsVsCCQEMuPt",
//       "MuonProtonTransversity.NProtons:MuonProtonTransversity.MuonPt_MeV.Mag()",
//       100, 0, 1000,
//       10, 0, 10,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx,
//       true));

//     Selections2D["NProtonsVsCCQEDPt"] =
//       (Selection2D("NProtonsVsCCQEDPt","NProtonsVsCCQEDPt",
//       "MuonProtonTransversity.NProtons:MuonProtonTransversity.DeltaPT_HMProton_MeV.Mag()",
//       100, 0, 1000,
//       10, 0, 10,
//       (CCQESelect + "&&" + CCQEDPhiTSelect).c_str(), DoLogx,
//       true));

//   }

  void InitSelections(){
//     (void)Verbosity;
//     InitCutStrings();
//     SetProtonMomentaSelections();
//     SetModeSelections();
//     SetTopologyDefinedECal();
//     TargetMassSelection();
//     SetCCQESelections();
//     FinalStateSpeciesSelections("MuonProtonTransversity",CCQESelect,"ccqe");
//     FinalStateSpeciesSelections("MuonProtonTransversity",CCQENoCCQENucEmitSelect,"ccqeNoNE");
//     FinalStateSpeciesSelections("MuonProtonTransversity",CCQENucEmitSelect,"ccqeNe");
//     DetailsPiProdSelection();
//     BasicPiProdSelections();
//     SetCC0PiSelections();
//     SetERecSelections();
//     Multiplicity2DSels();

//     ///Proton Momenta
//     Selections1D["FirstPMom"] =
//       (Selection1D("FirstPMom","FirstPMom",
//       "MuonProtonTransversity.FirstProtonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//       CCQEProtonMomentumSelect.c_str(), DoLogx));

//     Selections1D["HighestPMom"] =
//       (Selection1D("HighestPMom","HighestPMom",
//       "MuonProtonTransversity.HMProtonMomentum_MeV",
//       ProtonMomBins, 0, ProtonMomMax,
//       CCQEProtonMomentumSelect.c_str(), DoLogx));

  }


std::map<std::string,std::string> PreSelections;
void AddPreSelectionDescriptionXML(TXMLEngine &xmlengine,
  XMLDocPointer_t PreSelsNode){

  for(XMLNodePointer_t SelsDescriptNode = xmlengine.GetChild(PreSelsNode);
    (SelsDescriptNode != NULL);
    SelsDescriptNode = xmlengine.GetNext(SelsDescriptNode)){

    if(std::string(xmlengine.GetNodeName(SelsDescriptNode)) != "PreSelection"){
      std::cout << "[WARN]: Found unexpected node: "
        << xmlengine.GetNodeName(SelsDescriptNode)
        << " in the PreSelection element." << std::endl;
      continue;
    }
    std::string const &selName =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"Name");

    std::string const &selectionString =
      IOUtils::GetNodeContent(xmlengine, SelsDescriptNode);

    std::cout << "[INFO]: Found Preselection { Name: " << selName
    << ", Selection: " << selectionString << "}" << std::endl;

    if(!selName.length() || !selectionString.length()){
      return;
    }

    PreSelections[selName] = selectionString;
  }
}

namespace {
std::string CombineWeightStrings(const std::string& a,
  const std::string& b){

  if(!a.length()){
    if(!b.length()){

      return "";
    }
    return b;
  }
  if(!b.length()){
    return a;
  }
  return (a + "&&" + b);
}
}

void Add1DSelection(TXMLEngine &xmlengine,
  XMLDocPointer_t SelsDescriptNode){

    std::string const &selName =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"Name");

    if(selName.length() == 0){
      std::cout << "[ERROR]: Selection had no Name." << std::endl;
      return;
    }

    std::string const &selPreSelection =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"PreSelection");

    std::string const &selNBins =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"NBins");
    std::string const &selXBinLow =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinLow");
    std::string const &selXBinUpper =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinUpper");

    int NBins;
    double XBinLow;
    double XBinUpper;

    try{
      NBins = std::stoi(selNBins);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << selNBins << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinLow = std::stod(selXBinLow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << selXBinLow << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinUpper = std::stod(selXBinUpper);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse: " << selXBinUpper << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    bool valid, Invalidate;
    Invalidate = PGUtils::str2bool(
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,
       "InvalidateCache"), valid);
    if(!valid){
      Invalidate = true;
    }

    std::string selectionString = "";
    if(selPreSelection.length() && PreSelections.count(selPreSelection)){
      selectionString = CombineWeightStrings(PreSelections[selPreSelection],
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString"));
    } else {
      selectionString =
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString");
    }

    std::string const &drawString =
      IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "DrawString");

    if(!drawString.length() || !selectionString.length()){
      return;
    }

    Selections1D[selName] = Selection1D(
      selName,
      drawString,
      selectionString.c_str(),
      NBins,
      XBinLow,
      XBinUpper,
      Invalidate
      );

    std::cout << "[INFO]: Added Selection: " << Selections1D[selName]
      << std::endl;

}
void Add2DSelection(TXMLEngine &xmlengine,
  XMLDocPointer_t SelsDescriptNode){

    std::string const &selName =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"Name");

    if(selName.length() == 0){
      std::cout << "[ERROR]: Selection had no Name." << std::endl;
      return;
    }

    std::string const &selPreSelection =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"PreSelection");

    std::string const &selNXBins =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"NXBins");
    std::string const &selXBinLow =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinLow");
    std::string const &selXBinUpper =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"XBinUpper");

    std::string const &selNYBins =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"NYBins");
    std::string const &selYBinLow =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"YBinLow");
    std::string const &selYBinUpper =
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,"YBinUpper");

    int NXBins;
    double XBinLow;
    double XBinUpper;

    int NYBins;
    double YBinLow;
    double YBinUpper;

    try{
      NXBins = std::stoi(selNXBins);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse NXBins: " << selNXBins << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinLow = std::stod(selXBinLow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse XBinLow: " << selXBinLow << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      XBinUpper = std::stod(selXBinUpper);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse XBinUpper: " << selXBinUpper << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      NYBins = std::stoi(selNYBins);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse NYBins: " << selNYBins << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      YBinLow = std::stod(selYBinLow);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse YBinLow: " << selYBinLow << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    try{
      YBinUpper = std::stod(selYBinUpper);
    } catch (const std::invalid_argument& ia) {
      std::cout << "[ERROR]: Couldn't parse YBinUpper: " << selYBinUpper << "\n\t"
        << ia.what() << std::endl;
      return;
    }

    bool valid, Invalidate;
    Invalidate = PGUtils::str2bool(
      IOUtils::GetAttrValue(xmlengine, SelsDescriptNode,
       "InvalidateCache"), valid);
    if(!valid){
      Invalidate = true;
    }

    std::string selectionString = "";
    if(selPreSelection.length() && PreSelections.count(selPreSelection)){
      selectionString = CombineWeightStrings(PreSelections[selPreSelection],
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString"));
    } else {
      selectionString =
        IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "SelString");
    }

    std::string const &drawString =
      IOUtils::GetChildNodeContent(xmlengine, SelsDescriptNode, "DrawString");

    if(!drawString.length() || !selectionString.length()){
      return;
    }

    Selections2D[selName] = Selection2D(
      selName,
      drawString,
      selectionString.c_str(),
      NXBins,
      XBinLow,
      XBinUpper,
      NYBins,
      YBinLow,
      YBinUpper,
      Invalidate
      );

    std::cout << "[INFO]: Added Selection2D: " << Selections1D[selName]
      << std::endl;

}

void AddPlotSelectionDescriptionXML(TXMLEngine &xmlengine,
  XMLDocPointer_t PlotSelsNode){

  for(XMLNodePointer_t SelsDescriptNode = xmlengine.GetChild(PlotSelsNode);
    (SelsDescriptNode != NULL);
    SelsDescriptNode = xmlengine.GetNext(SelsDescriptNode)){

    if(std::string(xmlengine.GetNodeName(SelsDescriptNode)) == "Selection1D"){
      Add1DSelection(xmlengine, SelsDescriptNode);
    } else if(std::string(xmlengine.GetNodeName(SelsDescriptNode)) == "Selection2D"){
      Add2DSelection(xmlengine, SelsDescriptNode);
    } else {
      std::cout << "[WARN]: Found unexpected node: "
        << xmlengine.GetNodeName(SelsDescriptNode)
        << " in the Selection element." << std::endl;
      continue;
    }
  }
}

void InitSelectionsXML(std::string const &confFile){
  (void)Verbosity;

  TXMLEngine* xmlengine = new TXMLEngine;
  xmlengine->SetSkipComments();

  XMLNodePointer_t confSecNode =
    IOUtils::GetNamedChildElementOfDocumentRoot(*xmlengine,
      confFile, "Selections");
  if(!confSecNode){
    delete xmlengine; return;
  }

  for(XMLNodePointer_t PreSelsNode = xmlengine->GetChild(confSecNode);
    (PreSelsNode != NULL);
    PreSelsNode = xmlengine->GetNext(PreSelsNode)){
    if(std::string(xmlengine->GetNodeName(PreSelsNode)) != "PreSelections"){
      continue;
    }
    std::cout << "Found PreSelections element." << std::endl;
    AddPreSelectionDescriptionXML(*xmlengine,PreSelsNode);
    break;
  }

  for(XMLNodePointer_t PlotSelsNode = xmlengine->GetChild(confSecNode);
    (PlotSelsNode != NULL);
    PlotSelsNode = xmlengine->GetNext(PlotSelsNode)){
    if(std::string(xmlengine->GetNodeName(PlotSelsNode)) != "PlotSelections"){
      continue;
    }
    std::cout << "Found PlotSelections element." << std::endl;
    AddPlotSelectionDescriptionXML(*xmlengine,PlotSelsNode);
    break;

  }
  std::cout << "Done reading selections" << std::endl;
  delete xmlengine;
}
}
