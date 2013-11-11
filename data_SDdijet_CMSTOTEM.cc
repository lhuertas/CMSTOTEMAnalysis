
//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>

//OUR OWN CLASSES TO READ THE TREE
#include "MassParticles.h"
#include "MyBaseJet.h"
#include "MyBeamSpot.h"
#include "MyCaloJet.h"
#include "MyCastorDigi.h"
#include "MyCastorJet.h"
#include "MyCastorRecHit.h"
#include "MyDiJet.h"
#include "MyElectron.h"
#include "MyEvtId.h"
#include "MyFwdGap.h"
#include "MyGenJet.h"
#include "MyGenKin.h"
#include "MyGenMet.h"
#include "MyGenPart.h"
#include "MyHLTrig.h"
#include "MyJet.h"
#include "MyL1Trig.h"
#include "MyL1TrigOld.h"
//#include "MyMITEvtSel.h"
#include "MyMet.h"
#include "MyMuon.h"
#include "MyPFCand.h"
#include "MyPFJet.h"
#include "MyPUSumInfo.h"
#include "MyPart.h"
#include "MySimVertex.h"
#include "MyTracks.h"
#include "MyVertex.h"
#include "MyFSCHit.h"
#include "MyFSCDigi.h"

// TOTEM data formats
#include "T1Event.h"
#include "T2Event.h"
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpReconstructedProtonPair.h"
#include "RPRootDumpTrackInfo.h"
#include "RPRootDumpDigiInfo.h"
#include "RPRootDumpPatternInfo.h"

#include "analysis_tools.h"

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#define PI 3.141592653589793
using namespace std;

Double_t fFermiLike(Double_t *x, Double_t *par) {
  Double_t result = 0;
  result = 1.0/(TMath::Exp((par[0]-TMath::Sqrt(x[0]))/par[1]) + 1);
  return result;
}



void data_SDdijet_CMSTOTEM(string const& outputFileName = "data_SDdijet_CMSTOTEM.root", const Int_t nevt_max = -1){
  
  bool verbose = false;
  string treeName = "cms_totem";
  string jetCollName = "ak5PFJets";
  string jetCorrName = "ak5PFL2L3Residual";
  double ptJetMin = 15.0;
  double etaJetMax = 4.0;
  double etaMaxThreshold = 2.0;
  
  bool selectBunchCrossing = false;
  bool selectVertex = true;
  bool selectJets = false;
  bool selectEtaMax = false;
  bool selectEtaMin = false;
  bool selectZeroHitsT2Plus = false;
  bool selectZeroHitsT2Minus = false;
  bool selectSingleArmRecProton = true;
  bool selectDoubleArmRecProton = false;
  bool selectElastic = false;
  bool selectNonElastic = false;

  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;

  vector<string> hltPathNames;
  hltPathNames.push_back("HLT_L1DoubleEG3_FwdVeto_v1");
  hltPathNames.push_back("HLT_L1DoubleMu0_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20_RomanPotsOR_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20part1_v1");
  hltPathNames.push_back("HLT_L1DoubleJet24_v1");
  hltPathNames.push_back("HLT_L1DoubleJet20part2_v1");
  hltPathNames.push_back("HLT_L1Tech40_BPTXAND_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_1_v1");
  hltPathNames.push_back("HLT_L1Tech_HF9OR10_v1");
  hltPathNames.push_back("HLT_T1minbias_Tech55_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_2_v1");
  hltPathNames.push_back("HLT_L1Tech53_MB_3_v1");
  hltPathNames.push_back("HLT_RomanPots_Tech52_v1");
  hltPathNames.push_back("HLT_L1Tech54_ZeroBias_v1");
  hltPathNames.push_back("HLT_ZeroBias_v7");

  // Declaration of histograms
  map<string,TH1F*> histosTH1F;
  
  vector<string> selections;
  selections.push_back("All");
  selections.push_back("BunchCrossing");
  selections.push_back("HLT");
  selections.push_back("Vertex");
  selections.push_back("Jet");
  selections.push_back("EtaMax");
  selections.push_back("EtaMin");
  selections.push_back("ZeroHitsT2Plus");
  selections.push_back("ZeroHitsT2Minus");
  selections.push_back("SingleArmRP");
  selections.push_back("DoubleArmRP");
  selections.push_back("Elastic");
  selections.push_back("NonElastic");
  int nBinsEventSelection = selections.size();
  histosTH1F["EventSelection"] = new TH1F("EventSelection","EventSelection",nBinsEventSelection,0,nBinsEventSelection);
  for(size_t k = 0; k < selections.size(); ++k)
     histosTH1F["EventSelection"]->GetXaxis()->SetBinLabel( (k + 1), selections[k].c_str() );

  histosTH1F["bunchCrossingNumber"] = new TH1F("bunchCrossingNumber", "bunchCrossingNumber" , 3900 , 0 , 3900);

  histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
  histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

  int nBinsHLT = hltPathNames.size(); 
  histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);
  for(size_t k = 0; k < nBinsHLT; ++k) 
     histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1) , hltPathNames[k].c_str() );

  histosTH1F["vtx_zpos"] = new TH1F("vtx_zpos", "z(vtx)" , 150 , -30. , 30.);
  histosTH1F["vtx_xpos"] = new TH1F("vtx_xpos", "x(vtx)" , 150 , -1.5 , 1.5);
  histosTH1F["vtx_ypos"] = new TH1F("vtx_ypos", "y(vtx)" , 150 , -1.5 , 1.5);

  //histosTH1F["pt_gen"] = new TH1F("pt_gen" , "pt_gen;pt;nTracks" , 120 , 0 , 6);
  histosTH1F["track_pt"] = new TH1F("track_pt", "p_{T}(trk)" , 150 , 0. , 15.);
  histosTH1F["track_eta"] = new TH1F("track_eta", "#eta(trk)" , 200 , -5.2 , 5.2);
  histosTH1F["track_phi"] = new TH1F("track_phi", "#phi(trk)" , 200 , -M_PI , M_PI);
  histosTH1F["track_multiplicity"] = new TH1F("track_multiplicity", "n tracks" , 100 , 0 , 100);
  
  histosTH1F["jet_pt"] = new TH1F("jet_pt", "p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["jet_eta"] = new TH1F("jet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["jet_phi"] = new TH1F("jet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);

  histosTH1F["leadingJet_pt"] = new TH1F("leadingJet_pt", "p_{T}(jet)" , 20 , 0. , 200.);
  histosTH1F["leadingJet_eta"] = new TH1F("leadingJet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_phi"] = new TH1F("leadingJet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["leadingJet_pt_selected"] = new TH1F("leadingJet_pt_selected", "p_{T}(jet)" , 20 , 0. , 200.);
  histosTH1F["leadingJet_eta_selected"] = new TH1F("leadingJet_eta_selected", "#eta(jet)" , 20 , -5.2 , 5.2);

  histosTH1F["secondJet_pt"] = new TH1F("secondJet_pt", "p_{T}(jet)" , 20 , 0. , 200.);
  histosTH1F["secondJet_eta"] = new TH1F("secondJet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_phi"] = new TH1F("secondJet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_pt_selected"] = new TH1F("secondJet_pt_selected", "p_{T}(jet)" , 20 , 0. , 200.);
  histosTH1F["secondJet_eta_selected"] = new TH1F("secondJet_eta_selected", "#eta(jet)" , 20 , -5.2 , 5.2);

  histosTH1F["DeltaPtJet"] = new TH1F("Delta_pt_Jet", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet"] = new TH1F("Delta_eta_Jet", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet"] = new TH1F("Delta_phi_Jet", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["EtaJet_average"] = new TH1F("eta_Jet_average", "Average #eta(jet)" , 20 , -5.2 , 5.2);

  histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 20 , 0,11);
  histosTH1F["xi_plus_Reco"] = new TH1F("xi_plus", "#xi^{+}" , 20 , 0,1);
  histosTH1F["xi_minus_Reco"] = new TH1F("xi_minus", "#xi^{-}" , 20 , 0,1);
  histosTH1F["log_xi_plus_Reco"] = new TH1F("log_xi_plus", "Log #xi^{+}" , 20 , -3,0.5);
  histosTH1F["log_x_plus"] = new TH1F("log_x_plus", "Log x^{+}" , 20 , -4, 0);
  histosTH1F["log_x_minus"] = new TH1F("log_x_minus", "Log x^{-}" , 20 , -4, 0);
    

  histosTH1F["fscHit_energy"] = new TH1F("fscHit_energy", "FSC hit energy" , 150 , -100. , 200.);
  histosTH1F["fscHit_time"] = new TH1F("fscHit_time", "FSC hit time" , 150 , 0. , 300.);

  histosTH1F["t2_track_chi2Prob_zplus"] = new TH1F("t2_track_chi2Prob_zplus", "#chi^{2}" , 100 , 0. , 1.);
  histosTH1F["t2_track_entryX_zplus"] = new TH1F("t2_track_entryX_zplus", "x_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_entryY_zplus"] = new TH1F("t2_track_entryY_zplus", "y_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_multiplicity_zplus"] = new TH1F("t2_track_multiplicity_zplus", "n tracks" , 100 , 0 , 100);
  histosTH1F["t2_track_chi2Prob_zminus"] = new TH1F("t2_track_chi2Prob_zminus", "#chi^{2}" , 100 , 0. , 1.);
  histosTH1F["t2_track_entryX_zminus"] = new TH1F("t2_track_entryX_zminus", "x_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_entryY_zminus"] = new TH1F("t2_track_entryY_zminus", "y_{trk}" , 160 , -160. , 160.);
  histosTH1F["t2_track_multiplicity_zminus"] = new TH1F("t2_track_multiplicity_zminus", "n tracks" , 100 , 0 , 100);


  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};

  histosTH1F["proton_right_t"] = new TH1F("proton_right_t", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_cut"] = new TH1F("proton_right_t_cut", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_signal"] = new TH1F("proton_right_t_signal", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_true"] = new TH1F("proton_right_t_true", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_true_constbin"] = new TH1F("proton_right_t_true_constbin", "-t" , 20 , 0, 1);
  histosTH1F["proton_right_t_signal_constbin"] = new TH1F("proton_right_t_sigal_constbin", "-t" , 20 , 0, 1);
  histosTH1F["proton_right_t_signal_eff"] = new TH1F("proton_right_t_signal_eff", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_signal_effweight"] = new TH1F("proton_right_t_signal_effweight", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_signal_averagept_eff"] = new TH1F("proton_right_t_signal_averagept_eff", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_signal_averagept_effweight"] = new TH1F("proton_right_t_signal_averagept_effweight", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_halo"] = new TH1F("proton_right_t_halo", "-t" , 11 , tbins);
  histosTH1F["proton_right_t_halo_constbin"] = new TH1F("proton_right_t_halo_constbin", "-t" , 20 , 0, 1);
  histosTH1F["halo_right_pt30"] = new TH1F("halo_right_pt30", "-t halo" , 11 , tbins);
  histosTH1F["halo_right_pt30_constbin"] = new TH1F("halo_right_pt30_constbin", "-t halo" , 20 , 0, 1);
  histosTH1F["proton_right_chi2"] = new TH1F("proton_right_chi2", "#chi^{2}" , 20 , 0. , 100.);
  histosTH1F["proton_right_xi_signal"] = new TH1F("proton_right_xi_signal", "#xi Right RPs" , 20 , 0, 0.3);
  histosTH1F["proton_right_xi"] = new TH1F("proton_right_xi", "#xi Right RPs" , 20 , -0.1, 0.3);
  histosTH1F["proton_right_beta"] = new TH1F("proton_right_beta", "#beta Right RPs" , 20 , 0, 1);
  histosTH1F["proton_right_xi_cut"] = new TH1F("proton_right_xi_tcut", "#xi Right RPs" , 20 , 0, 0.3);
  histosTH1F["proton_right_logXi"] = new TH1F("proton_right_logXi","log(#xi)",20,-5.,0.);

  histosTH1F["proton_left_t"] = new TH1F("proton_left_t", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_cut"] = new TH1F("proton_left_t_cut", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal"] = new TH1F("proton_left_t_signal", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_true"] = new TH1F("proton_left_t_true", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_true_constbin"] = new TH1F("proton_left_t_true_constbin", "-t" , 20 , 0, 1);
  histosTH1F["proton_left_t_signal_constbin"] = new TH1F("proton_left_t_sigal_constbin", "-t" , 20 , 0, 1);
  histosTH1F["proton_left_t_signal_eff"] = new TH1F("proton_left_t_signal_eff", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal_effweight"] = new TH1F("proton_left_t_signal_effweight", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal_averagept_eff"] = new TH1F("proton_left_t_signal_averagept_eff", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_signal_averagept_effweight"] = new TH1F("proton_left_t_signal_averagept_effweight", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_halo"] = new TH1F("proton_left_t_halo", "-t" , 11 , tbins);
  histosTH1F["proton_left_t_halo_constbin"] = new TH1F("proton_left_t_halo_constbin", "-t" , 20 , 0, 1);
  histosTH1F["halo_left_pt30"] = new TH1F("halo_left_pt30", "-t halo" , 11 , tbins);
  histosTH1F["halo_left_pt30_constbin"] = new TH1F("halo_left_pt30_constbin", "-t halo" , 20 , 0, 1);
  histosTH1F["proton_left_chi2"] = new TH1F("proton_left_chi2", "#chi^{2}" , 20 , 0. , 100.);
  histosTH1F["proton_left_xi_signal"] = new TH1F("proton_left_xi_signal", "#xi Left RPs" , 20 , 0, 0.3);
  histosTH1F["proton_left_xi"] = new TH1F("proton_left_xi", "#xi Left RPs" , 20 , -0.1, 0.3);
  histosTH1F["proton_left_beta"] = new TH1F("proton_left_beta", "#beta Left RPs" , 20 , 0, 1);
  histosTH1F["proton_left_xi_cut"] = new TH1F("proton_left_xi_tcut", "#xi Left RPs" , 20 , 0, 0.3);
  histosTH1F["proton_left_logXi"] = new TH1F("proton_left_logXi","log(#xi)",20,-5.,0.);

  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);
  
  histosTH1F["xitotem_xicms_rightRPs"] = new TH1F("xitotem_xicms_rightRPs", "Right RPs" , 20 , -0.4 , 0.4);
  histosTH1F["xitotem_xicms_leftRPs"] = new TH1F("xitotem_xicms_leftRPs", "Left RPs" , 20 , -0.4 , 0.4);
  histosTH1F["xitotem_xicms_rightRPs_tcut"] = new TH1F("xitotem_xicms_rightRPs_tcut", "Right RPs" , 20 , -0.4 , 0.4);
  histosTH1F["xitotem_xicms_rightRPs_cut"] = new TH1F("xitotem_xicms_rightRPs_cut", "Right RPs" , 20 , -0.4 , 0.4);
  histosTH1F["xi_cms_totem_background_simulated"] = new TH1F("xitotem_xicms_rightRPs_simulated", "Right RPs" , 20 , -0.4 , 0.4);

  
  map<string,TH2F*> histosTH2F;
  histosTH2F["t2_track_multiplicity_vs_track_multiplicity"] = new TH2F("t2_track_multiplicity_vs_track_multiplicity","t2_track_multiplicity_vs_track_multiplicity", 100 , 0 , 100, 100 , 0 , 100);
  histosTH2F["t2_track_multiplicity_vs_leadingJet_pt"] = new TH2F("t2_track_multiplicity_vs_leadingJet_pt","t2_track_multiplicity_vs_leadingJet_pt", 150 , 0. , 150., 100 , 0 , 100);
  histosTH2F["t2_track_entryY_vs_entryX_zplus"] = new TH2F("t2_track_entryY_vs_entryX_zplus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);
  histosTH2F["t2_track_entryY_vs_entryX_zminus"] = new TH2F("t2_track_entryY_vs_entryX_zminus","t2_track_entryY_vs_entryX", 160 , -160. , 160., 160 , -160. , 160.);

  histosTH2F["proton_right_logXi_vs_pf_logXiPlus"] = new TH2F("proton_right_logXi_vs_pf_logXiPlus","proton_right_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_right_logXi_vs_pf_logXiMinus"] = new TH2F("proton_right_logXi_vs_pf_logXiMinus","proton_right_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_right_logXi_vs_t"] = new TH2F("proton_right_logXi_vs_t","proton_right_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
  histosTH2F["proton_right_xi_vs_pf_xiMinus"] = new TH2F("proton_right_xi_vs_pf_xiMinus","proton_right_xi_vs_pf_xiMinus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_right_xi_vs_pf_xiPlus"] = new TH2F("proton_right_xi_vs_pf_xiPlus","proton_right_xi_vs_pf_xiPlus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_right_t_vs_leadingJet_pt"] = new TH2F("proton_right_t_vs_leadingJet_pt","proton_right_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);

  histosTH2F["proton_left_logXi_vs_pf_logXiPlus"] = new TH2F("proton_left_logXi_vs_pf_logXiPlus","proton_left_logXi_vs_pf_logXiPlus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_pf_logXiMinus"] = new TH2F("proton_left_logXi_vs_pf_logXiMinus","proton_left_logXi_vs_pf_logXiMinus", 200, -5., 0., 200, -5., 0.);
  histosTH2F["proton_left_logXi_vs_t"] = new TH2F("proton_left_logXi_vs_t","proton_left_logXi_vs_t", 200, 0., 5., 200, -5., 0.);
  histosTH2F["proton_left_xi_vs_pf_xiMinus"] = new TH2F("proton_left_xi_vs_pf_xiMinus","proton_left_xi_vs_pf_xiMinus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_left_xi_vs_pf_xiPlus"] = new TH2F("proton_left_xi_vs_pf_xiPlus","proton_left_xi_vs_pf_xiPlus", 100, 0., 1., 100, 0., 1.);
  histosTH2F["proton_left_t_vs_leadingJet_pt"] = new TH2F("proton_left_t_vs_leadingJet_pt","proton_left_t_vs_leadingJet_pt", 150 , 0. , 150., 200 , 0. , 5.);
  
  double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;
  histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energy Vs Eta AllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energy Vs Eta Undefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energy Vs Eta Charged Hadron"] = new TH2F("energyVsEtaChargedHadron","energy Vs Eta Charged Hadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energy Vs Eta Electron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energy Vs Eta Muon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energy Vs Eta Gamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energy Vs Eta Neutral Hadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energy Vs Eta HadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energy Vs Eta GammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["xi+Vseta_max"] = new TH2F("xi+Vseta_max","#xi^{+] Vs #eta_{max}",200,0,5.2,nBinsEnergy,0,0.5);

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  //===================
  int i_tot = 0 , nevt_tot = 0;

  //vector<TString>* vfiles = new vector<TString>(1,"merged_reduced_8372_198903_LP_Jets1_1_test_v1.root"); 
//   vector<TString>* vfiles = new vector<TString>; 
//   vfiles->push_back( "CMS_Totem-LP_Jets1-198903_10_1_n7u-8369dc.root"); 

   const char *ext=".root";
 
   vector<TString>* vdirs = new vector<TString>; 
   vdirs->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/Jets1/");
   vdirs->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/Jets2/");
   vdirs->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/Jets1/");
   vdirs->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/Jets2/");
   
   vector<TString>* vfiles = new vector<TString>;
   for(vector<TString>::iterator itdirs = vdirs->begin(); itdirs != vdirs->end(); ++itdirs){
      TString& dirname = *itdirs;
       //vector<TString>* vfiles = new vector<TString>; 
      TSystemDirectory dir(dirname, dirname);
      TList *files = dir.GetListOfFiles();
      if (files) {
         TSystemFile *file;
         TString fname;
         TIter next(files);
         while ((file=(TSystemFile*)next())) {
             fname = file->GetName();
             if (!file->IsDirectory() && fname.EndsWith(ext)) {
                 TString root_file = dirname + string(fname.Data());
	         vfiles->push_back(root_file); cout<<root_file<<endl;      
             }
         }   
      } 
   }
 
  //Declaration of tree and its branches variables
//   TTree* tree = new TTree(treeName.c_str(),"");
  TTree* tree = NULL;
  MyEvtId*           evtId        = NULL;
  MyL1TrigOld*       l1Trig       = NULL;  
  MyHLTrig*          hltTrig      = NULL;
  //vector<MyGenPart>* genPart      = NULL;
  vector<MyTracks>*  track_coll   = NULL;
  vector<MyVertex>*  vertex_coll  = NULL;
  vector<MyPFJet>*   pfJet_coll   = NULL;
  vector<MyPFCand>*  pFlow_coll   = NULL;
  vector<MyFSCHit>*  fscHits_coll = NULL;
  vector<MyFSCDigi>* fscDigis_coll = NULL;
  //===================
  T2Event* t2_event = NULL;
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  RPRootDumpReconstructedProtonPair* rec_proton_pair  = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  map<unsigned int, RPRootDumpDigiInfo*> rp_digi_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_par_patterns_info;
  map<unsigned int, RPRootDumpPatternInfo*> rp_nonpar_patterns_info;
  map<unsigned int, std::vector<RPRootDumpTrackInfo>*> rp_multi_track_info;
  //===================  
  int n_evt_leadingJet = 0; 
  int n_evt_secondJet = 0; 
  int n_evt_leadingJet_selected = 0; 
  int n_evt_secondJet_selected = 0; 
  int n_evt_PF = 0; 

  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    //getting the tree form the current file
    tree = (TTree*) file->Get( treeName.c_str() );

    //ofstream ofs;
    //ofs.Open ("eff_ptJet2.txt");

    //Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

    //adding branches to the tree ----------------------------------------------------------------------
    tree->SetBranchAddress("cmsEvtUA",&evtId);
    tree->SetBranchAddress("cmsTrigUA",&l1Trig);
    tree->SetBranchAddress("cmsHLTTrigUA",&hltTrig);
    tree->SetBranchAddress("cmsTracksUA",&track_coll);
    tree->SetBranchAddress("cmsVerticesUA",&vertex_coll);
//    tree->SetBranchAddress("cmsPFJetsUA",&pfJet_coll);
     tree->SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll);
    tree->SetBranchAddress("cmsParticleFlowUA",&pFlow_coll);
    //tree->SetBranchAddress("cmsFSCHitsUA",&fscHits_coll);
    //tree->SetBranchAddress("cmsFSCDigisUA",&fscDigis_coll);
    tree->SetBranchAddress("branchT2EV.",&t2_event);
    tree->SetBranchAddress("rec_prot_left.",&rec_proton_left);
    tree->SetBranchAddress("rec_prot_right.",&rec_proton_right);
    tree->SetBranchAddress("rec_prot_pair.",&rec_proton_pair);
//     tree->SetBranchAddress("Evt",&evtId);
//     tree->SetBranchAddress("L1TrigOld",&l1Trig);
//     tree->SetBranchAddress("HLTrig",&hltTrig);
//     tree->SetBranchAddress("generalTracks",&track_coll); 
//     tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
//     tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
//     tree->SetBranchAddress("particleFlow",&pFlow_coll);
    //if(isMC) tree->SetBranchAddress("genPart",&genPart);
    std::vector<unsigned int> rp_list;
    rp_list.push_back(20); rp_list.push_back(21); rp_list.push_back(24); rp_list.push_back(25);
    rp_list.push_back(120); rp_list.push_back(121); rp_list.push_back(124); rp_list.push_back(125);
    char br_name[200];
    for (unsigned int a = 0; a < 2; ++a) {
       int s = 2;
       for (unsigned int r = 0; r < 6; r++) {
          unsigned int id = 100 * a + 10 * s + r;
          if( std::find(rp_list.begin(), rp_list.end(), id) == rp_list.end() ) continue;

          sprintf(br_name, "track_rp_%u.", id);
          std::cout << br_name << std::endl;
          tree->SetBranchAddress(br_name, &rp_track_info[id]);
       }
    } 
    
    /*//Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;*/
 
    int n_evt_Jet = 0;
    double eff_sum = 0; 
 
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){

      //printing the % of events done every 10k evts
      if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      bool passedHLT_jets = false;
      bool passedvtx = false;
      bool jet1_selected = false;
      bool jet2_selected = false;
      bool PF_eta_max = false;
      bool PF_eta_min = false;
      
      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = 1.;
      double event_weight_eff;
      double event_weight_averagept_eff;

      for (int itrig = 0 ; itrig < 128 ; ++itrig){
         if( l1Trig->PhysTrigWord[itrig] == 1) 
            histosTH1F["decisionPhysTrig"]->Fill( itrig, event_weight );
      }
        
      for (int itrig = 0 ; itrig < 64 ; ++itrig){
         if( l1Trig->TechTrigWord[itrig] == 1 )
            histosTH1F["decisionTechTrig"]->Fill( itrig, event_weight );
      }

      map<string,bool>::iterator it_hlt = (*hltTrig).HLTmap.begin();
      map<string,bool>::iterator it_hlt_end = (*hltTrig).HLTmap.end();
      for(; it_hlt != it_hlt_end; ++it_hlt){
         string const& hltName = it_hlt->first;
         vector<string>::const_iterator it_pos = find(hltPathNames.begin(),hltPathNames.end(),hltName);
         if(it_pos != hltPathNames.end()){
            size_t idx = it_pos - hltPathNames.begin();//cout <<hltName<<endl;
            if( hltName == "HLT_L1DoubleJet20part1_v1" || hltName == "HLT_L1DoubleJet20part2_v1"){
	      if( it_hlt->second == true ){
	         passedHLT_jets = true; 
		 histosTH1F["hltTrigFired"]->Fill( idx, event_weight );
	      }
	    }
	 }   
	         /*for(int ibin = 1; ibin <= histosTH1F["hltTrigFired"]->GetNbinsX(); ++ibin){
            if( hltName.c_str() != histosTH1F["hltTrigFired"]->GetXaxis()->GetBinLabel(ibin) ) continue;
            
            if( it_hlt->second ) 
               histosTH1F["hltTrigFired"]->Fill( histosTH1F["hltTrigFired"]->GetBinCenter( ibin ) );
         }*/ 
      }
      if(!passedHLT_jets) continue;
      
      //-------------------------------------------------------------------------------------------------
      //filling pt distribution for the generated particles
      //ie those from pythia generator, without reconstruction
      /*if(isMC){
        for(vector<MyGenPart>::iterator p=genPart->begin() ; p!=genPart->end() ; p++ )
          pt_gen->Fill(p->Pt());
      }*/
      
      //-------------------------------------------------------------------------------------------------
 	    
     // Vertices
      for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
        if (it_vtx!=vertex_coll->begin()) continue;
	if( it_vtx->ndof>4 ) passedvtx = true;   
      }
      if(!passedvtx) continue;
      
      // Tracks
      int n_tracks_selected = 0;
      for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
         histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
         histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
         histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );

         ++n_tracks_selected;
      }
      histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );

      
      
     //Jets with pt>30Gev and !eta!<2
      Double_t Jet1_pt; 
      Double_t Jet1_pz; 
      Double_t Jet1_E; 
      Double_t Jet2_pt; 
      Double_t Jet2_E; 
      Double_t Jet2_pz; 
      Double_t Jet1_eta; 
      Double_t Jet2_eta; 
      Double_t Jet1_phi; 
      Double_t Jet2_phi;
      Double_t eff; 
      Double_t averagept_eff; 

      ///Fit 
      TF1* func = new TF1("func", fFermiLike, 0., 20., 2);
      func->SetParameter(0,5.);
      func->SetParameter(1,0.5);      

       
      for(vector<MyPFJet>::iterator it_jet = pfJet_coll->begin() ; it_jet != pfJet_coll->end() ; ++it_jet){
         map<string,MyBaseJet>::iterator it_map = it_jet->mapjet.begin();
         for(; it_map != it_jet->mapjet.end(); ++it_map)
            if(verbose) cout << it_map->first << endl;

         MyBaseJet const& basejet = it_jet->mapjet[jetCorrName];
         histosTH1F["jet_pt"]->Fill( basejet.Pt(), event_weight  );
         if(basejet.Pt() > 0.) histosTH1F["jet_eta"]->Fill( basejet.Eta(), event_weight  );
         histosTH1F["jet_phi"]->Fill( basejet.Phi(), event_weight  );
      }
      if( pfJet_coll->size() > 0 ){
	 MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[jetCorrName];
         ++n_evt_leadingJet;
	 Jet1_E = leadingJet.E(); 
	 Jet1_pz = leadingJet.Pz(); 
	 Jet1_pt = leadingJet.Pt(); 
	 Jet1_eta = leadingJet.Eta(); 
	 Jet1_phi = leadingJet.Phi(); 
	 histosTH1F["leadingJet_pt"]->Fill(leadingJet.Pt(), event_weight  );
	 histosTH1F["leadingJet_eta"]->Fill( leadingJet.Eta(), event_weight  );
	 histosTH1F["leadingJet_phi"]->Fill( leadingJet.Phi(), event_weight  );
	 
	 if(Jet1_pt > 30. && Jet1_eta<4. && Jet1_eta>-4.){ jet1_selected = true;
	    ++n_evt_leadingJet_selected;
	 }
      }
      if(!jet1_selected) continue;
      
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet[jetCorrName];
         ++n_evt_secondJet;
         Jet2_E = secondJet.E(); 
         Jet2_pz = secondJet.Pz(); 
         Jet2_pt = secondJet.Pt(); 
	 Jet2_eta = secondJet.Eta(); 
	 Jet2_phi = secondJet.Phi(); 
         histosTH1F["secondJet_pt"]->Fill( secondJet.Pt(), event_weight  );  
	 histosTH1F["secondJet_eta"]->Fill( secondJet.Eta(), event_weight  );
	 histosTH1F["secondJet_phi"]->Fill( secondJet.Phi(), event_weight  );

	 if(Jet2_pt > 30. && Jet2_eta<4. && Jet2_eta>-4.){  jet2_selected = true;
	    ++n_evt_secondJet_selected;
	 }
      }
      if(!jet2_selected) continue;
      ++n_evt_Jet;
      double average_pt = (Jet1_pt+Jet2_pt)/2;
      eff = func->Eval(Jet2_pt);
      eff_sum += eff;
      averagept_eff = func->Eval(average_pt);
    
      double protons_correction = (1.067+1.066+1.060+1.063)/4; 

      event_weight_eff = protons_correction/eff;
      event_weight_averagept_eff = 1/averagept_eff;

      double x_plus = ((Jet1_E+Jet1_pz)+(Jet2_E+Jet2_pz))/8000;
      double x_minus = ((Jet1_E-Jet1_pz)+(Jet2_E-Jet2_pz))/8000;
       
      // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999;
      double cm = 8000;
      
      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 double pz = it_pfcand->Pz();
	 
         // HF eta rings 29, 30, 40, 41
         if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;
	 
	 bool Barrel_region = false;
         bool Endcap_region = false;
         bool Transition_region = false;
         bool Forward_region = false;
	 bool passed_threshold = true;
	 
         if( eta>=-1.4 && eta<=1.4 ) Barrel_region = true;
         else if( (eta>=-2.6 && eta<-1.4) || (eta>1.4 && eta<=2.6) ) Endcap_region = true;
         else if( (eta>=-3.2 && eta<-2.6) || (eta>2.6 && eta<=3.2) ) Transition_region = true;
         else if( eta<-3.2 || eta>3.2 ) Forward_region = true;
         else cout << "ERROR!!!!!!!!!" << endl;

	 //Applying threshold
	 if (Barrel_region == true){
	   //if(partType == MyPFCand::h0 && energy>=1.4) Treshold_Barrel1==true;
	   if(partType == MyPFCand::h0 && energy < 1.4) passed_threshold = false; 
	   if(partType == MyPFCand::gamma && energy<0.9) passed_threshold = false; 
	 }  
	 if (Endcap_region == true){
	   if(partType == MyPFCand::h0 && energy<2.7) passed_threshold = false; 
	   if(partType == MyPFCand::gamma && energy<2.5) passed_threshold = false; 
	 }  
	 if (Transition_region == true){
	   if(partType == MyPFCand::h0 && energy<3.8) passed_threshold = false; 
	   if(partType == MyPFCand::gamma && energy<2.5) passed_threshold = false; 
	   if(partType == MyPFCand::h_HF && energy<4) passed_threshold = false; 
	   if(partType == MyPFCand::egamma_HF && energy<3.5) passed_threshold = false; 
	 }  
	 if (Forward_region == true){
	   if(partType == MyPFCand::h_HF && energy<4) passed_threshold = false; 
	   if(partType == MyPFCand::egamma_HF && energy<3.5) passed_threshold = false; 
	 }  

	 if( !passed_threshold ) continue;
	 
	 soma1 += (energy + pz);
	 soma2 += (energy - pz);

         if (eta > eta_max) {eta_max = eta; PF_eta_max = true;} 
	 if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}

  
      }
      if(!PF_eta_max) continue;
      if(!PF_eta_min) continue;
      ++n_evt_PF;

      double xi_plus_Reco = soma1/cm;
      double xi_minus_Reco = soma2/cm;
      double delta_eta_maxmin = eta_max - eta_min;
      double correction = 1/0.8;

      histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
      histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
      histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
      histosTH1F["xi_plus_Reco"]->Fill( xi_plus_Reco, event_weight  );

      histosTH1F["xi_minus_Reco"]->Fill( xi_minus_Reco, event_weight  );
      histosTH1F["log_xi_plus_Reco"]->Fill( log10(xi_plus_Reco), event_weight  );
      histosTH2F["xi+Vseta_max"]->Fill( eta_max, xi_plus_Reco, event_weight );


      // TOTEM T2
      vector<double> const& t2_trk_entryX = t2_event->TrkEntryX;
      vector<double> const& t2_trk_entryY = t2_event->TrkEntryY;
      vector<double> const& t2_trk_entryZ =  t2_event->TrkEntryZ;
      vector<double> const& t2_trk_chiProb =  t2_event->TrkChiProb;

      int n_t2_tracks_selected = 0;
      int n_t2_tracks_selected_zplus = 0;
      int n_t2_tracks_selected_zminus = 0;
      size_t n_t2_tracks = t2_trk_chiProb.size();
      for(size_t i_t2_trk = 0; i_t2_trk < n_t2_tracks; ++i_t2_trk){
         double trk_entryZ = t2_trk_entryZ[i_t2_trk];
         int zside = ( trk_entryZ >= 0. ) ? 1 : -1;
         if( zside > 0 )
            histosTH1F["t2_track_chi2Prob_zplus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );
         else
            histosTH1F["t2_track_chi2Prob_zminus"]->Fill( t2_trk_chiProb[i_t2_trk], event_weight );

         // Select tracks
         if( t2_trk_chiProb[i_t2_trk] < 0.2 ) continue;

         ++n_t2_tracks_selected;
         if( zside > 0 ) ++n_t2_tracks_selected_zplus;
         else            ++n_t2_tracks_selected_zminus;

         if( zside > 0 ){
            histosTH1F["t2_track_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
            histosTH1F["t2_track_entryY_zplus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
            histosTH2F["t2_track_entryY_vs_entryX_zplus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
         } else{
            histosTH1F["t2_track_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], event_weight );
            histosTH1F["t2_track_entryY_zminus"]->Fill( t2_trk_entryY[i_t2_trk], event_weight );
            histosTH2F["t2_track_entryY_vs_entryX_zminus"]->Fill( t2_trk_entryX[i_t2_trk], t2_trk_entryY[i_t2_trk], event_weight );
         }
      }
      if( selectZeroHitsT2Plus && (n_t2_tracks_selected_zplus > 0) ) continue;
      histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Plus", event_weight );

      if( selectZeroHitsT2Minus && (n_t2_tracks_selected_zminus > 0) ) continue;
      histosTH1F["EventSelection"]->Fill( "ZeroHitsT2Minus", event_weight );

      bool proton_right_valid = rec_proton_right->valid;
      bool proton_left_valid = rec_proton_left->valid;
      if( selectSingleArmRecProton && (proton_right_valid && proton_left_valid) ) continue;
      histosTH1F["EventSelection"]->Fill( "SingleArmRP", event_weight );

      if( selectDoubleArmRecProton && !(proton_right_valid && proton_left_valid) ) continue;
      histosTH1F["EventSelection"]->Fill( "DoubleArmRP", event_weight );

      bool tag_elastic_top45_bot56 = elastic_top45_bot56(rp_track_info);      
      bool tag_elastic_bot45_top56 = elastic_bot45_top56(rp_track_info);      
      if( selectElastic && !(tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      histosTH1F["EventSelection"]->Fill( "Elastic", event_weight );

      if( selectNonElastic && (tag_elastic_top45_bot56 || tag_elastic_bot45_top56) ) continue;
      histosTH1F["EventSelection"]->Fill( "NonElastic", event_weight );

      // RP protons
      double xi_totem;
      double chi2_proton_right = rec_proton_right->chi2;
      double chindf_proton_right = rec_proton_right->chindf;
      double xi_proton_right = rec_proton_right->xi;//cout<<xi_proton_right<<endl;
      double t_proton_right = rec_proton_right->t;
      double beta_proton_right = x_minus/-xi_proton_right;
      bool xi_region_right = -xi_proton_right>=0.03 && -xi_proton_right<=0.1;
      bool t_region_right = -t_proton_right>=0.03 && -t_proton_right<=1;
      bool good_proton_right = proton_right_valid; // && (xi_proton_right < 0.);
      double chi2_proton_left = rec_proton_left->chi2;
      double chindf_proton_left = rec_proton_left->chindf;
      double xi_proton_left = rec_proton_left->xi;//cout<<xi_proton_left<<endl;
      double t_proton_left = rec_proton_left->t;
      double beta_proton_left = x_minus/-xi_proton_left;
      bool xi_region_left = -xi_proton_left>=0.03 && -xi_proton_left<=0.1;
      bool t_region_left = -t_proton_left>=0.03 && -t_proton_left<=1;
      bool good_proton_left = proton_left_valid; // && (xi_proton_left < 0.);
 
      
      if( good_proton_right &&  !good_proton_left){ 
         histosTH1F["proton_right_chi2"]->Fill( chi2_proton_right, event_weight );
         histosTH1F["proton_right_xi"]->Fill(-xi_proton_right, event_weight );
         histosTH1F["proton_right_t"]->Fill( -t_proton_right, event_weight );

         if (xi_region_right && t_region_right) {
 
             if (xi_minus_Reco>0.12){ 
    	        histosTH1F["proton_right_xi_cut"]->Fill(-xi_proton_right, event_weight );
	     } 

             histosTH1F["xitotem_xicms_rightRPs"]->Fill( xi_minus_Reco+xi_proton_right, event_weight );

             if(xi_minus_Reco+xi_proton_right>0){
               histosTH1F["proton_right_t_halo"]->Fill(-t_proton_right, event_weight );
               histosTH1F["proton_right_t_halo_constbin"]->Fill(-t_proton_right, event_weight );
             }

             histosTH1F["proton_right_logXi"]->Fill( log10(-xi_proton_right), event_weight );
             histosTH1F["pf_xiMinus_minus_proton_right_xi"]->Fill( (xi_minus_Reco + xi_proton_right), event_weight );
             histosTH2F["proton_right_logXi_vs_pf_logXiPlus"]->Fill( log10(xi_plus_Reco),log10(-xi_proton_right), event_weight );
             histosTH2F["proton_right_logXi_vs_pf_logXiMinus"]->Fill( log10(xi_minus_Reco),log10(-xi_proton_right), event_weight );
             histosTH2F["proton_right_xi_vs_pf_xiPlus"]->Fill( xi_plus_Reco, -xi_proton_right, event_weight );
             histosTH2F["proton_right_xi_vs_pf_xiMinus"]->Fill( xi_minus_Reco, -xi_proton_right, event_weight );
             histosTH2F["proton_right_logXi_vs_t"]->Fill( -t_proton_right, log10(-xi_proton_right), event_weight );
             histosTH1F["proton_right_t_cut"]->Fill( -t_proton_right, event_weight );
             histosTH1F["proton_right_logXi"]->Fill( log10(-xi_proton_right), event_weight );

             if (xi_minus_Reco+xi_proton_right>0) continue;
                histosTH1F["proton_right_t_signal"]->Fill( -t_proton_right, event_weight );
                histosTH1F["proton_right_t_signal_constbin"]->Fill( -t_proton_right, event_weight );
                histosTH1F["proton_right_t_signal_eff"]->Fill( -t_proton_right, eff );
                histosTH1F["proton_right_t_signal_effweight"]->Fill( -t_proton_right, event_weight_eff );
                histosTH1F["proton_right_t_signal_averagept_eff"]->Fill( -t_proton_right, averagept_eff );
                histosTH1F["proton_right_t_signal_averagept_effweight"]->Fill( -t_proton_right, event_weight_averagept_eff );
                histosTH1F["proton_right_xi_signal"]->Fill(-xi_proton_right, event_weight );
                histosTH1F["proton_right_beta"]->Fill(beta_proton_right, event_weight );
         
             if( pfJet_coll->size() > 0 ){
               MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
               histosTH2F["proton_right_t_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), -t_proton_right, event_weight );
             }
         }
      }
   
      if( good_proton_left &&  !good_proton_right){   
         histosTH1F["proton_left_chi2"]->Fill( chi2_proton_left, event_weight );
         histosTH1F["proton_left_xi"]->Fill(-xi_proton_left, event_weight );
         histosTH1F["proton_left_t"]->Fill( -t_proton_left, event_weight );

         if (xi_region_left && t_region_left) {

             if (xi_minus_Reco>0.12){
                histosTH1F["proton_left_xi_cut"]->Fill(-xi_proton_left, event_weight );
             }

             histosTH1F["xitotem_xicms_leftRPs"]->Fill( xi_minus_Reco+xi_proton_left, event_weight );

             if(xi_plus_Reco+xi_proton_left>0){
               histosTH1F["proton_left_t_halo"]->Fill(-t_proton_left, event_weight );
               histosTH1F["proton_left_t_halo_constbin"]->Fill(-t_proton_left, event_weight );
             }

/*             histosTH1F["proton_left_logXi"]->Fill( log10(-xi_proton_left), event_weight );
             histosTH1F["pf_xiMinus_minus_proton_left_xi"]->Fill( (xi_plus_Reco + xi_proton_left), event_weight );
             histosTH2F["proton_left_logXi_vs_pf_logXiPlus"]->Fill( log10(xi_plus_Reco),log10(-xi_proton_left), event_weight );
             histosTH2F["proton_left_logXi_vs_pf_logXiMinus"]->Fill( log10(xi_minus_Reco),log10(-xi_proton_left), event_weight );
             histosTH2F["proton_left_xi_vs_pf_xiPlus"]->Fill( xi_plus_Reco, -xi_proton_left, event_weight );
             histosTH2F["proton_left_xi_vs_pf_xiMinus"]->Fill( xi_minus_Reco, -xi_proton_left, event_weight );
             histosTH2F["proton_left_logXi_vs_t"]->Fill( -t_proton_left, log10(-xi_proton_left), event_weight );
             histosTH1F["proton_left_t_cut"]->Fill( -t_proton_left, event_weight );
             histosTH1F["proton_left_logXi"]->Fill( log10(-xi_proton_left), event_weight );
*/
             if (xi_plus_Reco+xi_proton_left>0) continue;
                histosTH1F["proton_left_t_signal"]->Fill( -t_proton_left, event_weight );
                histosTH1F["proton_left_t_signal_constbin"]->Fill( -t_proton_left, event_weight );
                histosTH1F["proton_left_t_signal_eff"]->Fill( -t_proton_left, eff );
                histosTH1F["proton_left_t_signal_effweight"]->Fill( -t_proton_left, event_weight_eff );
                histosTH1F["proton_left_t_signal_averagept_eff"]->Fill( -t_proton_left, averagept_eff );
                histosTH1F["proton_left_t_signal_averagept_effweight"]->Fill( -t_proton_left, event_weight_averagept_eff );
                histosTH1F["proton_left_xi_signal"]->Fill(-xi_proton_left, event_weight );
                histosTH1F["proton_left_beta"]->Fill(beta_proton_left, event_weight );

             if( pfJet_coll->size() > 0 ){
               MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[ jetCorrName ];
               histosTH2F["proton_left_t_vs_leadingJet_pt"]->Fill( leadingJet.Pt(), -t_proton_left, event_weight );
             }
         }
      }

      MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet[jetCorrName];
      histosTH1F["leadingJet_pt_selected"]->Fill(leadingJet.Pt(), event_weight  );
      histosTH1F["leadingJet_eta_selected"]->Fill( leadingJet.Eta(), event_weight  );


      MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet[jetCorrName];
      histosTH1F["secondJet_pt_selected"]->Fill( secondJet.Pt(), event_weight  );
      histosTH1F["secondJet_eta_selected"]->Fill( secondJet.Eta(), event_weight  );

      Double_t deltapt = abs(Jet1_pt - Jet2_pt);
      Double_t deltaeta = abs(Jet1_eta - Jet2_eta);
      Double_t deltaphi = abs(Jet1_phi - Jet2_phi);
      Double_t eta_average = 0.5*(Jet1_eta + Jet2_eta);
      histosTH1F["DeltaPtJet"]->Fill( deltapt, event_weight  );
      histosTH1F["DeltaEtaJet"]->Fill( deltaeta, event_weight  );
      histosTH1F["DeltaPhiJet"]->Fill( 1-(deltaphi)/PI, event_weight  );
      histosTH1F["EtaJet_average"]->Fill( eta_average, event_weight  );                    
      histosTH1F["log_x_minus"]->Fill( log10(x_minus), event_weight  );
      


    }//end loop for events
    file->Close();
    cout << "eff= "<<eff_sum/n_evt_Jet<<endl;
  
  }//end of loop over files

  //beam halo
  int halo_signal_pt30 = 0;
  for (int a = 1; a<=histosTH1F["proton_right_xi_cut"]->GetEntries(); a++){
       double xi_cms_right = histosTH1F["xi_minus_Reco"]->GetRandom();
       double xi_totem_right = histosTH1F["proton_right_xi_cut"]->GetRandom();
       histosTH1F["xi_cms_totem_background_simulated"]->Fill(xi_cms_right-xi_totem_right,1.0); 
       if (xi_cms_right-xi_totem_right<0.){++halo_signal_pt30;}
   }   
   cout<<"background_right:  "<<halo_signal_pt30<<endl;

  
   double nevents_t_halo_pt30[11];
   double halo_nosignal_pt30 = histosTH1F["proton_right_t_halo"]->GetEntries();
   double scale_pt30 = halo_signal_pt30/halo_nosignal_pt30;
   double nevents_data[11];
   double nevents_true[11];
   int nbins_const = histosTH1F["proton_right_t_halo_constbin"]->GetNbinsX();
   double nevents_t_halo_pt30_constbin[nbins_const];
   double nevents_data_constbin[nbins_const];
   double nevents_true_constbin[nbins_const];

   for (int b = 1; b<=11; b++){
       nevents_t_halo_pt30[b] = histosTH1F["proton_right_t_halo"]->GetBinContent(b)*scale_pt30;
       nevents_data[b] = histosTH1F["proton_right_t_signal"]->GetBinContent(b);
       nevents_true[b] = nevents_data[b]-nevents_t_halo_pt30[b]; 
       histosTH1F["halo_right_pt30"] -> SetBinContent(b, nevents_t_halo_pt30[b]);
       histosTH1F["proton_right_t_true"] -> SetBinContent(b, nevents_true[b]);
cout<<b<<"  data: "<<nevents_data[b]<<"  halo: "<<nevents_t_halo_pt30[b]<<"   true: "<<nevents_true[b]<<endl;
   }

   for (int c = 1; c<=nbins_const; c++){
       nevents_t_halo_pt30_constbin[c] = histosTH1F["proton_right_t_halo_constbin"]->GetBinContent(c)*scale_pt30;
       nevents_data_constbin[c] = histosTH1F["proton_right_t_signal_constbin"]->GetBinContent(c);
       nevents_true_constbin[c] = nevents_data_constbin[c]-nevents_t_halo_pt30_constbin[c]; 
       histosTH1F["halo_right_pt30_constbin"] -> SetBinContent(c, nevents_t_halo_pt30_constbin[c]);
       histosTH1F["proton_right_t_true_constbin"] -> SetBinContent(c, nevents_true_constbin[c]);
   }

  //output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();

  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();
                                  it_histo != histosTH1F.end(); ++it_histo)
     (*it_histo).second->Write();
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

 
  output->Close();

}
