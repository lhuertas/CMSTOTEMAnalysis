//STANDARD ROOT INCLUDES
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
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
#include <TRandom.h>

//STANDARD C++ INCLUDES
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

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
#include "RPRootDumpReconstructedProton.h"
#include "RPRootDumpTrackInfo.h"

#include "rp_aperture_config.h"
#include <exception>
//ROOUNFOLD CLASSES
#include "/storage1/lhuertas/Analysis/CMSTotem_RP/CMSTotem/Workspace/RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "/storage1/lhuertas/Analysis/CMSTotem_RP/CMSTotem/Workspace/RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "/storage1/lhuertas/Analysis/CMSTotem_RP/CMSTotem/Workspace/RooUnfold-1.1.1/src/RooUnfoldBinByBin.h"

#define PI 3.141592653589793
using namespace std;

Double_t beta_fit(Double_t *x, Double_t *par ){
  Double_t result = 0;
  result = par[0]+par[1]*x[0]+ par[2]*x[0]*x[0];
  return result;
}


double etaBinsHCALBoundaries[] = {-5.205, -4.903, -4.730,
                                  -4.552, -4.377, -4.204, -4.027, -3.853, -3.677, -3.503, -3.327, -3.152,
                                  -3.000, -2.868, -2.650, -2.500,
                                  -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479,
                                  -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.870, -0.783,
                                  -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
                                  0.000, 0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696,
                                  0.783, 0.870, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392,
                                  1.479, 1.566, 1.653, 1.740, 1.830, 1.930, 2.043, 2.172, 2.322,
                                  2.500, 2.650, 2.868, 3.000,
                                  3.152, 3.327, 3.503, 3.677, 3.853, 4.027, 4.204, 4.377, 4.552,
                                  4.730, 4.903, 5.205}; // 41 + 41 bins

class MyZeroBiasData {
   public:
      MyZeroBiasData() {}
      ~MyZeroBiasData() {}

      bool   proton_rec_left_valid;
      double proton_rec_left_t;
      double proton_rec_left_xi;
      bool   proton_rec_right_valid;
      double proton_rec_right_t;
      double proton_rec_right_xi;
      bool   vtx_valid;
      double vtx_ndof;
      double vtx_x;
      double vtx_y;
      double vtx_z;
      double leadingJet_pt;
      double leadingJet_eta;
      double secondJet_pt;
      double secondJet_eta;
      bool rp_track_valid_120; 
      bool rp_track_valid_121; 
      bool rp_track_valid_122; 
      bool rp_track_valid_123; 
      bool rp_track_valid_124;
      bool rp_track_valid_125;
      bool rp_track_valid_020;
      bool rp_track_valid_021;
      bool rp_track_valid_022;
      bool rp_track_valid_023;
      bool rp_track_valid_024;
      bool rp_track_valid_025;
      double rp_x_024;
      double rp_y_024;
      double rp_x_025;
      double rp_y_025;
      double rp_x_124;
      double rp_y_124;
      double rp_x_125;
      double rp_y_125;
      double xi_cms_plus; 
      double xi_cms_minus; 
};


void PYTHIA_reco_pu(string const& outputFileName = "/storage1/lhuertas/Analysis/CMSTotem/MC/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/pythia_reco_pu.root", const Int_t nevt_max = -1){
  
  bool verbose = false;
  string treeName = "evt";//"cms_totem";
  string jetCollName = "ak5PFJets";
  string jetCorrName = "ak5PFL2L3Residual";
  jetCorrName = "ak5PFL2L3"; 
  
  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


  // Declaration of histograms
  float xbins=50;
  float bin2=15;
//  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};
  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.16, 0.20, 0.25, 0.32, 0.42, 0.52, 0.65, 1.};
  map<string,TH1F*> histosTH1F;
//   histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
//   histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

/*  int nBinsHLT = hltPathNames.size(); 
  histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);*/
//   for(size_t k = 0; k < nBinsHLT; ++k) 
//      histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1) , hltPathNames[k].c_str() );

  double weight_st, xi_cms_st, xi_totem_st, xi_cms_minus_totem_st;
  TTree* small_tree = new TTree("small_tree","");
  small_tree->Branch("weight",&weight_st,"weight/D");
  small_tree->Branch("xi_cms",&xi_cms_st,"xi_cms/D");
  small_tree->Branch("xi_totem",&xi_totem_st,"xi_totem/D");
  small_tree->Branch("xi_cms_minus_totem",&xi_cms_minus_totem_st,"xi_cms_minus_totem/D");

  double weight_st_back, xi_cms_st_back, xi_totem_st_back, xi_cms_minus_totem_st_back;
  TTree* small_tree_back = new TTree("small_tree_back","");
  small_tree_back->Branch("weight",&weight_st_back,"weight/D");
  small_tree_back->Branch("xi_cms",&xi_cms_st_back,"xi_cms/D");
  small_tree_back->Branch("xi_totem",&xi_totem_st_back,"xi_totem/D");
  small_tree_back->Branch("xi_cms_minus_totem",&xi_cms_minus_totem_st_back,"xi_cms_minus_totem/D");

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

  histosTH1F["leadingJet_pt"] = new TH1F("leadingJet_pt", "p_{T}(jet)" , 20 , 0. , 230.);
  histosTH1F["leadingJet_eta"] = new TH1F("leadingJet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_phi"] = new TH1F("leadingJet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["leadingJet_pt_rec"] = new TH1F("leadingJet_pt_rec", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_gen"] = new TH1F("leadingJet_pt_gen", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_eta_rec"] = new TH1F("leadingJet_eta_rec", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_gen"] = new TH1F("leadingJet_eta_gen", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["leadingJet_pt_right_signal_kin_cut"] = new TH1F("leadingJet_pt_right_signal_kin_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_right_signal_kin_cut_bin2"] = new TH1F("leadingJet_pt_right_signal_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_left_signal_kin_cut_bin2"] = new TH1F("leadingJet_pt_left_signal_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_right_signal_kint_kinxi_cut"] = new TH1F("leadingJet_pt_right_signal_kint_kinxi_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_left_signal_kint_kinxi_cut"] = new TH1F("leadingJet_pt_left_signal_kint_kinxi_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_right_signal_kint_kinxi_cut_bin2"] = new TH1F("leadingJet_pt_right_signal_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_left_signal_kint_kinxi_cut_bin2"] = new TH1F("leadingJet_pt_left_signal_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_right_backg_kin_cut"] = new TH1F("leadingJet_pt_right_backg_kin_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_right_backg_kin_cut"] = new TH1F("leadingJet_pt_right_backg_kin_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_right_backg_kin_cut_bin2"] = new TH1F("leadingJet_pt_right_backg_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_left_backg_kin_cut_bin2"] = new TH1F("leadingJet_pt_left_backg_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_right_backg_kint_kinxi_cut"] = new TH1F("leadingJet_pt_right_backg_kint_kinxi_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_left_backg_kint_kinxi_cut"] = new TH1F("leadingJet_pt_left_backg_kint_kinxi_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["leadingJet_pt_right_backg_kint_kinxi_cut_bin2"] = new TH1F("leadingJet_pt_right_backg_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_pt_left_backg_kint_kinxi_cut_bin2"] = new TH1F("leadingJet_pt_left_backg_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["leadingJet_eta_right_signal_kin_cut"] = new TH1F("leadingJet_eta_right_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_left_signal_kin_cut"] = new TH1F("leadingJet_eta_left_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_right_signal_kint_kinxi_cut"] = new TH1F("leadingJet_eta_right_signal_kint_kinxi_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_left_signal_kint_kinxi_cut"] = new TH1F("leadingJet_eta_left_signal_kint_kinxi_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_right_backg_kin_cut"] = new TH1F("leadingJet_eta_right_backg_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_eta_right_backg_kint_kinxi_cut"] = new TH1F("leadingJet_eta_right_backg_kint_kinxi_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["leadingJet_pt_left_signal_kin_cut"] = new TH1F("leadingJet_pt_left_signal_kin_cut", "p_{T}(jet)" , 20 , 0. , 200.);
  histosTH1F["leadingJet_pt_left_backg_kin_cut"] = new TH1F("leadingJet_pt_left_backg_kin_cut", "p_{T}(jet)" , 20 , 0. , 200.);

  histosTH1F["secondJet_pt"] = new TH1F("secondJet_pt", "p_{T}(jet)" , 20 , 0. , 230.);
  histosTH1F["secondJet_eta"] = new TH1F("secondJet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_phi"] = new TH1F("secondJet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_pt_rec"] = new TH1F("secondJet_pt_rec", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_gen"] = new TH1F("secondJet_pt_gen", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_eta_rec"] = new TH1F("secondJet_eta_rec", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_gen"] = new TH1F("secondJet_eta_gen", "eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["secondJet_pt_right_signal_kin_cut"] = new TH1F("secondJet_pt_right_signal_kin_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["secondJet_pt_right_signal_kin_cut_bin2"] = new TH1F("secondJet_pt_right_signal_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_left_signal_kin_cut_bin2"] = new TH1F("secondJet_pt_left_signal_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_right_signal_kint_kinxi_cut"] = new TH1F("secondJet_pt_right_signal_kint_kinxi_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["secondJet_pt_left_signal_kint_kinxi_cut"] = new TH1F("secondJet_pt_left_signal_kint_kinxi_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["secondJet_pt_right_signal_kint_kinxi_cut_bin2"] = new TH1F("secondJet_pt_right_signal_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_left_signal_kint_kinxi_cut_bin2"] = new TH1F("secondJet_pt_left_signal_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_right_backg_kin_cut_bin2"] = new TH1F("secondJet_pt_right_backg_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_left_backg_kin_cut_bin2"] = new TH1F("secondJet_pt_left_backg_kin_cut_bin2", "p_{T}(jet)" , bin2 , 0. , 200.);
  histosTH1F["secondJet_pt_right_backg_kin_cut"] = new TH1F("secondJet_pt_right_backg_kin_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["secondJet_pt_left_backg_kin_cut"] = new TH1F("secondJet_pt_left_backg_kin_cut", "p_{T}(jet)" , xbins , 0. , 200.);
  histosTH1F["secondJet_pt_right_backg_kint_kinxi_cut"] = new TH1F("secondJet_pt_right_backg_kint_kinxi_cut", "p_{T}(jet)" , xbins, 0. , 200.);
  histosTH1F["secondJet_pt_left_backg_kint_kinxi_cut"] = new TH1F("secondJet_pt_left_backg_kint_kinxi_cut", "p_{T}(jet)" , xbins, 0. , 200.);
  histosTH1F["secondJet_pt_right_backg_kint_kinxi_cut_bin2"] = new TH1F("secondJet_pt_right_backg_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2, 0. , 200.);
  histosTH1F["secondJet_pt_left_backg_kint_kinxi_cut_bin2"] = new TH1F("secondJet_pt_left_backg_kint_kinxi_cut_bin2", "p_{T}(jet)" , bin2, 0. , 200.);
  histosTH1F["secondJet_eta_right_signal_kin_cut"] = new TH1F("secondJet_eta_right_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_left_signal_kin_cut"] = new TH1F("secondJet_eta_left_signal_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_right_signal_kint_kinxi_cut"] = new TH1F("secondJet_eta_right_signal_kint_kinxi_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_left_signal_kint_kinxi_cut"] = new TH1F("secondJet_eta_left_signal_kint_kinxi_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_right_backg_kin_cut"] = new TH1F("secondJet_eta_right_backg_kin_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_eta_right_backg_kint_kinxi_cut"] = new TH1F("secondJet_eta_right_backg_kint_kinxi_cut", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_pt_left_signal_kin_cut"] = new TH1F("secondJet_pt_left_signal_kin_cut", "p_{T}(jet)" , 20 , 0. , 200.);

  histosTH1F["mass_jets_rec"] = new TH1F("mass_jets_rec", "mass(jet)" , bin2 , 0 , 450);
  histosTH1F["mass_jets_gen"] = new TH1F("mass_jets_gen", "mass(jet)" , bin2 , 0 , 450);
  histosTH1F["DeltaPtJet_right_signal_kin_cut"] = new TH1F("Delta_pt_Jet_right_signal_kin_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaPtJet_left_signal_kin_cut"] = new TH1F("Delta_pt_Jet_left_signal_kin_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_right_signal_kin_cut"] = new TH1F("Delta_eta_Jet_right_signal_kin_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaEtaJet_left_signal_kin_cut"] = new TH1F("Delta_eta_Jet_left_signal_kin_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_right_signal_kin_cut"] = new TH1F("Delta_phi_Jet_right_signal_kin_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["DeltaPhiJet_left_signal_kin_cut"] = new TH1F("Delta_phi_Jet_left_signal_kin_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_right_signal_kin_cut"] = new TH1F("mass_Jet_right_signal_kin_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["Mass_Jet_left_signal_kin_cut"] = new TH1F("mass_Jet_left_signal_kin_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_right_backg_kin_cut"] = new TH1F("Delta_pt_Jet_right_backg_kin_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_right_backg_kin_cut"] = new TH1F("Delta_eta_Jet_right_backg_kin_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_right_backg_kin_cut"] = new TH1F("Delta_phi_Jet_right_backg_kin_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_right_backg_kin_cut"] = new TH1F("mass_Jet_right_backg_kin_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_left_backg_kin_cut"] = new TH1F("Delta_pt_Jet_left_backg_kin_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_left_backg_kin_cut"] = new TH1F("Delta_eta_Jet_left_backg_kin_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_left_backg_kin_cut"] = new TH1F("Delta_phi_Jet_left_backg_kin_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_left_backg_kin_cut"] = new TH1F("mass_Jet_left_backg_kin_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_right_signal_kint_kinxi_cut"] = new TH1F("Delta_pt_Jet_right_signal_kint_kinxi_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_right_signal_kint_kinxi_cut"] = new TH1F("Delta_eta_Jet_right_signal_kint_kinxi_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_right_signal_kint_kinxi_cut"] = new TH1F("Delta_phi_Jet_right_signal_kint_kinxi_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_right_signal_kint_kinxi_cut"] = new TH1F("mass_Jet_right_signal_kint_kinxi_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_left_signal_kint_kinxi_cut"] = new TH1F("Delta_pt_Jet_left_signal_kint_kinxi_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_left_signal_kint_kinxi_cut"] = new TH1F("Delta_eta_Jet_left_signal_kint_kinxi_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_left_signal_kint_kinxi_cut"] = new TH1F("Delta_phi_Jet_left_signal_kint_kinxi_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_left_signal_kint_kinxi_cut"] = new TH1F("mass_Jet_left_signal_kint_kinxi_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_right_backg_kint_kinxi_cut"] = new TH1F("Delta_pt_Jet_right_backg_kint_kinxi_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_right_backg_kint_kinxi_cut"] = new TH1F("Delta_eta_Jet_right_backg_kint_kinxi_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_right_backg_kint_kinxi_cut"] = new TH1F("Delta_phi_Jet_right_backg_kint_kinxi_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_right_backg_kint_kinxi_cut"] = new TH1F("mass_Jet_right_backg_kint_kinxi_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["DeltaPtJet_left_backg_kint_kinxi_cut"] = new TH1F("Delta_pt_Jet_left_backg_kint_kinxi_cut", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet_left_backg_kint_kinxi_cut"] = new TH1F("Delta_eta_Jet_left_backg_kint_kinxi_cut", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet_left_backg_kint_kinxi_cut"] = new TH1F("Delta_phi_Jet_left_backg_kint_kinxi_cut", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet_left_backg_kint_kinxi_cut"] = new TH1F("mass_Jet_left_backg_kint_kinxi_cut", "Mass(jet)" , 20 , 0 , 450);
  histosTH1F["EtaJet_average_right_signal_kin_cut"] = new TH1F("eta_Jet_average_right_signal_kin_cut", "Average #eta(jet)" , xbins , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_signal_kin_cut_bin2"] = new TH1F("eta_Jet_average_right_signal_kin_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_signal_kint_kinxi_cut"] = new TH1F("eta_Jet_average_right_signal_kint_kinxi_cut", "Average #eta(jet)" , xbins , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_signal_kint_kinxi_cut_bin2"] = new TH1F("eta_Jet_average_right_signal_kint_kinxi_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_backg_kin_cut"] = new TH1F("eta_Jet_average_right_backg_kin_cut", "Average #eta(jet)" , xbins , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_backg_kin_cut_bin2"] = new TH1F("eta_Jet_average_right_backg_kin_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_backg_kint_kinxi_cut"] = new TH1F("eta_Jet_average_right_backg_kint_kinxi_cut", "Average #eta(jet)" , xbins , -5.2 , 5.2);
  histosTH1F["EtaJet_average_right_backg_kint_kinxi_cut_bin2"] = new TH1F("eta_Jet_average_right_backg_kint_kinxi_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_signal_kin_cut"] = new TH1F("eta_Jet_average_left_signal_kin_cut", "Average #eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_signal_kin_cut_bin2"] = new TH1F("eta_Jet_average_left_signal_kin_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_signal_kint_kinxi_cut"] = new TH1F("eta_Jet_average_left_signal_kint_kinxi_cut", "Average #eta(jet)" , xbins , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_signal_kint_kinxi_cut_bin2"] = new TH1F("eta_Jet_average_left_signal_kint_kinxi_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_backg_kin_cut"] = new TH1F("eta_Jet_average_left_backg_kin_cut", "Average #eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_backg_kin_cut_bin2"] = new TH1F("eta_Jet_average_left_backg_kin_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_backg_kint_kinxi_cut"] = new TH1F("eta_Jet_average_left_backg_kint_kinxi_cut", "Average #eta(jet)" , xbins , -5.2 , 5.2);
  histosTH1F["EtaJet_average_left_backg_kint_kinxi_cut_bin2"] = new TH1F("eta_Jet_average_left_backg_kint_kinxi_cut_bin2", "Average #eta(jet)" , bin2 , -5.2 , 5.2);

  histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 50 , 0,11);
  histosTH1F["xi_plus_Reco"] = new TH1F("xi+", "#xi^{+}" , 82 , 0,4);
  histosTH1F["xi_minus_Reco"] = new TH1F("xi-", "#xi^{-}" , 82 , 0,4);
  histosTH1F["logxi_plus"] = new TH1F("logxi+", "Log #xi^{+}" , 82 , -3,0.5);
  histosTH1F["logxi_plus_gen"] = new TH1F("logxi+_gen", "Log #xi_{+}^{gen}" , 82 , -3,0.5);
  histosTH1F["logxi_minus_gen"] = new TH1F("logxi-_gen", "Log #xi_{-}^{gen}" , 82 , -3,0.5);
  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);

  histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_minus_rec"] = new TH1F("xi_proton_minus_rec", "xi_proton_minus_sel" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_plus_rec"] = new TH1F("xi_proton_plus_rec", "xi_proton_plus_sel" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_right_zb"] = new TH1F("xi_proton_right_zb", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_left_zb"] = new TH1F("xi_proton_left_zb", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right"] = new TH1F("xi_proton_signal_right", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kint"] = new TH1F("xi_proton_signal_right_kint", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kinxi_kint"] = new TH1F("xi_proton_signal_right_kinxi_kint", "xi_proton_right" , 25, -0.05, 0.2);
  histosTH1F["xi_minus_rec_proton_signal_right_kinxi_kint"] = new TH1F("xi_minus_rec_proton_signal_right_kinxi_kint", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_minus_rec_proton_signal_right_kinxi_kint_cut"] = new TH1F("xi_minus_rec_proton_signal_right_kinxi_kint_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_minus_rec_proton_backg_right_kinxi_kint"] = new TH1F("xi_minus_rec_proton_backg_right_kinxi_kint", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_minus_rec_proton_backg_right_kinxi_kint_cut"] = new TH1F("xi_minus_rec_proton_backg_right_kinxi_kint_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kint_bin25"] = new TH1F("xi_proton_signal_right_kint_bin25", "xi_proton_right" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kin_cut"] = new TH1F("xi_proton_signal_right_kin_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kin_cut_bin2"] = new TH1F("xi_proton_signal_right_kin_cut_bin2", "xi_proton_right" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kint_cut"] = new TH1F("xi_proton_signal_right_kint_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kint_cut_bin2"] = new TH1F("xi_proton_signal_right_kint_cut_bin2", "xi_proton_right" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kint_kinxi_cut"] = new TH1F("xi_proton_signal_right_kint_kinxi_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_right_kint_kinxi_cut_bin2"] = new TH1F("xi_proton_signal_right_kint_kinxi_cut_bin2", "xi_proton_right" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right"] = new TH1F("xi_proton_backg_right", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kin"] = new TH1F("xi_proton_backg_right_kin", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kin_cut"] = new TH1F("xi_proton_backg_right_kin_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kin_cut_bin2"] = new TH1F("xi_proton_backg_right_kin_cut_bin2", "xi_proton_right" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kint_cut"] = new TH1F("xi_proton_backg_right_kint_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kint_cut_bin2"] = new TH1F("xi_proton_backg_right_kint_cut_bin2", "xi_proton_right" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kint_kinxi_cut"] = new TH1F("xi_proton_backg_right_kint_kinxi_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kint_kinxi_cut_bin2"] = new TH1F("xi_proton_backg_right_kint_kinxi_cut_bin2", "xi_proton_right" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kint"] = new TH1F("xi_proton_backg_right_kint", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kinxi_kint"] = new TH1F("xi_proton_backg_right_kinxi_kint", "xi_proton_right" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kinxi_kint_sempeso"] = new TH1F("xi_proton_backg_right_kinxi_kint_sempeso", "xi_proton_right" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kinxi_kint_pesocuadrado"] = new TH1F("xi_proton_backg_right_kinxi_kint_pesocuadrado", "xi_proton_right" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_backg_right_kint_bin25"] = new TH1F("xi_proton_backg_right_kint_bin25", "xi_proton_right" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_right_both"] = new TH1F("xi_proton_right_both", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_both_right_kint"] = new TH1F("xi_proton_both_right_kint", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_both_right_kint_cut"] = new TH1F("xi_proton_both_right_kint_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_both_right_kin_cut"] = new TH1F("xi_proton_both_right_kin_cut", "xi_proton_right" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kin_cut_bin2"] = new TH1F("xi_proton_signal_left_kin_cut_bin2", "xi_proton_left" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kint_cut"] = new TH1F("xi_proton_signal_left_kint_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kint_cut_bin2"] = new TH1F("xi_proton_signal_left_kint_cut_bin2", "xi_proton_left" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kint_kinxi_cut"] = new TH1F("xi_proton_signal_left_kint_kinxi_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kint_kinxi_cut_bin2"] = new TH1F("xi_proton_signal_left_kint_kinxi_cut_bin2", "xi_proton_left" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kin_cut"] = new TH1F("xi_proton_backg_left_kin_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kin_cut_bin2"] = new TH1F("xi_proton_backg_left_kin_cut_bin2", "xi_proton_left" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kint_cut"] = new TH1F("xi_proton_backg_left_kint_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kint_cut_bin2"] = new TH1F("xi_proton_backg_left_kint_cut_bin2", "xi_proton_left" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kint_kinxi_cut"] = new TH1F("xi_proton_backg_left_kint_kinxi_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kint_kinxi_cut_bin2"] = new TH1F("xi_proton_backg_left_kint_kinxi_cut_bin2", "xi_proton_left" , bin2, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left"] = new TH1F("xi_proton_signal_left", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kint"] = new TH1F("xi_proton_signal_left_kint", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kinxi_kint"] = new TH1F("xi_proton_signal_left_kinxi_kint", "xi_proton_left" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kint_bin25"] = new TH1F("xi_proton_signal_left_kint_bin25", "xi_proton_left" , 25, -0.05, 0.2);
  histosTH1F["xi_proton_signal_left_kin_cut"] = new TH1F("xi_proton_signal_left_kin_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left"] = new TH1F("xi_proton_backg_left", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kin"] = new TH1F("xi_proton_backg_left_kin", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kin_cut"] = new TH1F("xi_proton_backg_left_kin_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_proton_backg_left_kint"] = new TH1F("xi_proton_backg_left_kint", "xi_proton_left" , xbins, -0.05, 0.2); 
  histosTH1F["xi_proton_backg_left_kinxi_kint"] = new TH1F("xi_proton_backg_left_kinxi_kint", "xi_proton_left" , 25, -0.05, 0.2); 
  histosTH1F["xi_proton_backg_left_kint_bin25"] = new TH1F("xi_proton_backg_left_kint_bin25", "xi_proton_left" , 25, -0.05, 0.2); 
  histosTH1F["xi_minus_cms_zb_kint"] = new TH1F("xi_minus_cms_zb_kint", "xi_cms_zb" , xbins, -0.05, 0.2);
  histosTH1F["xi_minus_cms_zb_kint_cut"] = new TH1F("xi_minus_cms_zb_kint_cut", "xi_cms_zb" , xbins, -0.05, 0.2);
  histosTH1F["xi_plus_rec_proton_signal_left_kinxi_kint"] = new TH1F("xi_plus_rec_proton_signal_left_kinxi_kint", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_plus_rec_proton_signal_left_kinxi_kint_cut"] = new TH1F("xi_plus_rec_proton_signal_left_kinxi_kint_cut", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_plus_rec_proton_backg_left_kinxi_kint"] = new TH1F("xi_plus_rec_proton_backg_left_kinxi_kint", "xi_proton_left" , xbins, -0.05, 0.2);
  histosTH1F["xi_plus_rec_proton_backg_left_kinxi_kint_cut"] = new TH1F("xi_plus_rec_proton_backg_left_kinxi_kint_cut", "xi_proton_left" , xbins, -0.05, 0.2);

  histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" , 11, tbins);
  histosTH1F["t_proton_plus_rec"] = new TH1F("t_proton_plus_rec", "t_proton_plus" , 11, tbins);
  histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" , 11, tbins);
  histosTH1F["t_proton_minus_rec"] = new TH1F("t_proton_minus_rec", "t_proton_minus" , 11, tbins);
  histosTH1F["t_proton_right_zb"] = new TH1F("t_proton_right_zb", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_right_zb_sel"] = new TH1F("t_proton_right_zb_sel", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_left_zb"] = new TH1F("t_proton_left_zb", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_left_zb_sel"] = new TH1F("t_proton_left_zb_sel", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_signal_right"] = new TH1F("t_proton_signal_right", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_signal_right_kin_cut"] = new TH1F("t_proton_signal_right_kin_cut", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_signal_right_kint_kinxi_cut"] = new TH1F("t_proton_signal_right_kint_kinxi_cut", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_signal_right_kint_kinxi_cut_bin"] = new TH1F("t_proton_signal_right_kint_kinxi_cut_bin", "t_proton_right" , xbins, 0,1);
  histosTH1F["t_proton_backg_right"] = new TH1F("t_proton_backg_right", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_backg_right_sel"] = new TH1F("t_proton_backg_right_sel", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_backg_right_kin_cut"] = new TH1F("t_proton_backg_right_kin_cut", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_backg_right_kint_kinxi_cut"] = new TH1F("t_proton_backg_right_kint_kinxi_cut", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_backg_right_kint_kinxi_cut_bin"] = new TH1F("t_proton_backg_right_kint_kinxi_cut_bin", "t_proton_right" , xbins, 0, 1);
  histosTH1F["t_proton_backg_right_kin"] = new TH1F("t_proton_backg_right_kin", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_signal_left"] = new TH1F("t_proton_signal_left", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_backg_left"] = new TH1F("t_proton_backg_left", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_right_both"] = new TH1F("t_proton_right_both", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_left_both"] = new TH1F("t_proton_left_both", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_both_right_kin_cut"] = new TH1F("t_proton_backg_both_kin_cut", "t_proton_minus" , 11, tbins);
  histosTH1F["t_proton_signal_left_kint_kinxi_cut"] = new TH1F("t_proton_signal_left_kint_kinxi_cut", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_signal_left_kint_kinxi_cut_bin"] = new TH1F("t_proton_signal_left_kint_kinxi_cut_bin", "t_proton_left" , xbins, 0,1);
  histosTH1F["t_proton_backg_left_kint_kinxi_cut"] = new TH1F("t_proton_backg_left_kint_kinxi_cut", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_backg_left_kint_kinxi_cut_bin"] = new TH1F("t_proton_backg_left_kint_kinxi_cut_bin", "t_proton_left" , xbins, 0, 1);
  histosTH1F["t_proton_signal_left_kin_cut"] = new TH1F("t_proton_signal_left_kin_cut", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_backg_left_sel"] = new TH1F("t_proton_backg_left_sel", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_backg_left_kin_cut"] = new TH1F("t_proton_backg_left_kin_cut", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_backg_left_kin"] = new TH1F("t_proton_backg_left_kin", "t_proton_left" , 11, tbins);

  Float_t bin[16] = {-0.4, -0.112, -0.096, -0.08, -0.064, -0.048, -0.032, -0.016, 0, 0.048, 0.112, 0.176, 0.24, 0.304, 0.368, 0.4};
  histosTH1F["xi_cms_totem_right_signal"] = new TH1F("xi_cms_totem_right_signal", "xi_cms_totem_right" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_signal_kin"] = new TH1F("xi_cms_totem_right_signal_kin", "xi_cms_totem_right" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_signal_kinxi_kint"] = new TH1F("xi_cms_totem_right_signal_kinxi_kint", "xi_cms_totem_right" , 15, bin);
  histosTH1F["xi_cms_totem_right_zb"] = new TH1F("xi_cms_totem_right_zb", "xi_cms_totem_right" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_zb_kin"] = new TH1F("xi_cms_totem_right_zb_kin", "xi_cms_totem_right" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_zb_kin_bin"] = new TH1F("xi_cms_totem_right_zb_kin_bin", "xi_cms_totem_right" , 50, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_zb_kinxi_kint"] = new TH1F("xi_cms_totem_right_zb_kinxi_kint", "xi_cms_totem_right" , 15, bin);
  histosTH1F["xi_cms_totem_right_zb_kinxi_kint_zb"] = new TH1F("xi_cms_totem_right_zb_kinxi_kint_zb", "xi_cms_totem_right" , 15, bin);
  histosTH1F["xi_cms_totem_right_zb_kinxi_kint_bin"] = new TH1F("xi_cms_totem_right_zb_kinxi_kint_bin", "xi_cms_totem_right" , 50, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_both"] = new TH1F("xi_cms_totem_right_both", "xi_cms_totem_right" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_right_both_kin"] = new TH1F("xi_cms_totem_right_both_kin", "xi_cms_totem_right" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_left_zb"] = new TH1F("xi_cms_totem_left_zb", "xi_cms_totem_left" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_left_signal"] = new TH1F("xi_cms_totem_left_signal", "xi_cms_totem_left" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_left_signal_kin"] = new TH1F("xi_cms_totem_left_signal_kin", "xi_cms_totem_left" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_left_zb_kin"] = new TH1F("xi_cms_totem_left_zb_kin", "xi_cms_totem_left" , xbins, -0.4, 0.4);
  histosTH1F["xi_cms_totem_left_signal_kinxi_kint"] = new TH1F("xi_cms_totem_left_signal_kinxi_kint", "xi_cms_totem_left" , 15, bin);
  histosTH1F["xi_cms_totem_left_zb_kinxi_kint"] = new TH1F("xi_cms_totem_left_zb_kinxi_kint", "xi_cms_totem_left" , 15, bin);

  histosTH1F["beta_proton_minus_signal"] = new TH1F("beta_proton_minus_signal", "beta_proton_minus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_minus_signal_kin_cut"] = new TH1F("beta_proton_minus_signal_kin_cut", "beta_proton_minus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_minus_signal_kin_cut_bin2"] = new TH1F("beta_proton_minus_signal_kin_cut_bin2", "beta_proton_minus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_minus_signal_kint_kinxi_cut"] = new TH1F("beta_proton_minus_signal_kint_kinxi_cut", "beta_proton_minus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_minus_signal_kint_kinxi_cut_bin2"] = new TH1F("beta_proton_minus_signal_kint_kinxi_cut_bin2", "beta_proton_minus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_minus_backg"] = new TH1F("beta_proton_minus_backg", "beta_proton_minus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_minus_backg_kin_cut"] = new TH1F("beta_proton_minus_backg_kin_cut", "beta_proton_minus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_minus_backg_kin_cut_bin2"] = new TH1F("beta_proton_minus_backg_kin_cut_bin2", "beta_proton_minus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_minus_backg_kint_kinxi_cut"] = new TH1F("beta_proton_minus_backg_kint_kinxi_cut", "beta_proton_minus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_minus_backg_kint_kinxi_cut_bin2"] = new TH1F("beta_proton_minus_backg_kint_kinxi_cut_bin2", "beta_proton_minus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_plus_signal"] = new TH1F("beta_proton_plus_signal", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_plus_signal_kin_cut"] = new TH1F("beta_proton_plus_signal_kin_cut", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_plus_backg"] = new TH1F("beta_proton_plus_backg", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_plus_backg_kin_cut"] = new TH1F("beta_proton_plus_backg_kin_cut", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_both_right"] = new TH1F("beta_proton_both_right", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_both_right_kin_cut"] = new TH1F("beta_proton_both_right_kin_cut", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_plus_signal_kin_cut_bin2"] = new TH1F("beta_proton_plus_signal_kin_cut_bin2", "beta_proton_plus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_plus_signal_kint_kinxi_cut"] = new TH1F("beta_proton_plus_signal_kint_kinxi_cut", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_plus_signal_kint_kinxi_cut_bin2"] = new TH1F("beta_proton_plus_signal_kint_kinxi_cut_bin2", "beta_proton_plus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_plus_backg_kin_cut_bin2"] = new TH1F("beta_proton_plus_backg_kin_cut_bin2", "beta_proton_plus" , bin2, 0, 0.9 );
  histosTH1F["beta_proton_plus_backg_kint_kinxi_cut"] = new TH1F("beta_proton_plus_backg_kint_kinxi_cut", "beta_proton_plus" , xbins, 0, 0.9 );
  histosTH1F["beta_proton_plus_backg_kint_kinxi_cut_bin2"] = new TH1F("beta_proton_plus_backg_kint_kinxi_cut_bin2", "beta_proton_plus" , bin2, 0, 0.9 );

  histosTH1F["log_x_parton_rec_selec_jj"] = new TH1F("log_x_parton_rec_selec_jj", "log_x_parton" , bin2 , -4, 0);
  histosTH1F["log_x_parton_rec_minus_jj"] = new TH1F("log_x_parton_rec_minus_jj", "log_x_parton" , bin2 , -4, 0);
  histosTH1F["log_x_parton_rec_plus_jj"] = new TH1F("log_x_parton_rec_plus_jj", "log_x_parton" , bin2 , -4, 0);
  histosTH1F["log_x_parton_gen_selec_jj"] = new TH1F("log_x_parton_gen_selec_jj", "log_x_parton" , bin2 , -4, 0);
  histosTH1F["log_x_parton_gen_minus_jj"] = new TH1F("log_x_parton_gen_minus_jj", "log_x_parton" , bin2 , -4, 0);
  histosTH1F["log_x_parton_gen_plus_jj"] = new TH1F("log_x_parton_gen_plus_jj", "log_x_parton" , bin2 , -4, 0);
  histosTH1F["log_x_parton_gen_minus_jj_recsel"] = new TH1F("log_x_parton_gen_minus_jj_recsel", "log_x_parton" , bin2, -4, 0);
  histosTH1F["log_x_parton_gen_plus_jj_recsel"] = new TH1F("log_x_parton_gen_plus_jj_recsel", "log_x_parton" , bin2, -4, 0);
  histosTH1F["log_x_parton_minus_signal_kin_cut"] = new TH1F("log_x_parton_minus_signal_kin_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_minus_signal_kin_cut_bin2"] = new TH1F("log_x_parton_minus_signal_kin_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_signal_kint_kinxi_cut"] = new TH1F("log_x_parton_minus_signal_kint_kinxi_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_minus_signal_kint_kinxi_cut_bin2"] = new TH1F("log_x_parton_minus_signal_kint_kinxi_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_signal_kin_cut"] = new TH1F("log_x_parton_plus_signal_kin_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_minus_backg_kin_cut"] = new TH1F("log_x_parton_minus_backg_kin_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_minus_backg_kin_cut_bin2"] = new TH1F("log_x_parton_minus_backg_kin_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_minus_backg_kint_kinxi_cut"] = new TH1F("log_x_parton_minus_backg_kint_kinxi_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_minus_backg_kint_kinxi_cut_bin2"] = new TH1F("log_x_parton_minus_backg_kint_kinxi_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_backg_kin_cut"] = new TH1F("log_x_parton_plus_backg_kin_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_plus_signal_kin_cut_bin2"] = new TH1F("log_x_parton_plus_signal_kin_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_signal_kint_kinxi_cut"] = new TH1F("log_x_parton_plus_signal_kint_kinxi_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_plus_signal_kint_kinxi_cut_bin2"] = new TH1F("log_x_parton_plus_signal_kint_kinxi_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_backg_kin_cut_bin2"] = new TH1F("log_x_parton_plus_backg_kin_cut_bin2", "log_x_parton" , bin2, -4, 0 );
  histosTH1F["log_x_parton_plus_backg_kint_kinxi_cut"] = new TH1F("log_x_parton_plus_backg_kint_kinxi_cut", "log_x_parton" , xbins, -4, 0 );
  histosTH1F["log_x_parton_plus_backg_kint_kinxi_cut_bin2"] = new TH1F("log_x_parton_plus_backg_kint_kinxi_cut_bin2", "log_x_parton" , bin2, -4, 0 );

  map<string,TH2F*> histosTH2F;
  double energyMin = -10.;
  double energyMax = 190.;
  int nBinsEnergy = 1000;
  histosTH2F["energyVsEtaAllTypes"] = new TH2F("energyVsEtaAllTypes","energyVsEtaAllTypes",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaUndefined"] = new TH2F("energyVsEtaUndefined","energyVsEtaUndefined",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaChargedHadron"] = new TH2F("energyVsEtaChargedHadron","energyVsEtaChargedHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaElectron"] = new TH2F("energyVsEtaElectron","energyVsEtaElectron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaMuon"] = new TH2F("energyVsEtaMuon","energyVsEtaMuon",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaGamma"] = new TH2F("energyVsEtaGamma","energyVsEtaGamma",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaNeutralHadron"] = new TH2F("energyVsEtaNeutralHadron","energyVsEtaNeutralHadron",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHF"] = new TH2F("energyVsEtaHadronHF","energyVsEtaHadronHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFEcalEnergy"] = new TH2F("energyVsEtaHadronHFEcalEnergy","energyVsEtaHadronHFEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaHadronHFNoEcalEnergy"] = new TH2F("energyVsEtaHadronHFNoEcalEnergy","energyVsEtaHadronHFNoEcalEnergy",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["energyVsEtaEGammaHF"] = new TH2F("energyVsEtaEGammaHF","energyVsEtaEGammaHF",82,etaBinsHCALBoundaries,nBinsEnergy,energyMin,energyMax);
  histosTH2F["xi_plus_reco_gen"] = new TH2F("xi+","xi+",82,0,0.5,82,0,0.5);
  histosTH2F["xi_minus_reco_gen"] = new TH2F("xi-","xi-",82,0,0.5,82,0,0.5);
  histosTH2F["logxi_plus_reco_gen"] = new TH2F("logxi+","xi+",82,-3,0.5,82,-3,0.5);
  histosTH2F["logxi_minus_reco_gen"] = new TH2F("logxi-","xi-",82,-3,0.5,82,-3,0.5);

  histosTH2F["rp_track_pos_y_vs_x_020"] = new TH2F("rp_track_pos_y_vs_x_020", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);
  histosTH2F["rp_track_pos_y_vs_x_120"] = new TH2F("rp_track_pos_y_vs_x_120", "y(RP track) vs x(RP track)" , 200, -10., 10., 500, -50., 50.);

  histosTH2F["proton_plus_xi_vs_t"] = new TH2F("proton_plus_xi_vs_t","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
  histosTH2F["proton_minus_xi_vs_t"] = new TH2F("proton_minus_xi_vs_t","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

  histosTH2F["proton_plus_xi_vs_t_accepted"] = new TH2F("proton_plus_xi_vs_t_accepted","proton_plus_xi_vs_t", 200, 0., 5., 200, 0., 1.);
  histosTH2F["proton_minus_xi_vs_t_accepted"] = new TH2F("proton_minus_xi_vs_t_accepted","proton_minus_xi_vs_t", 200, 0., 5., 200, 0., 1.);

  histosTH2F["pos_y_vs_x_proton_plus_accepted_020"] = new TH2F("pos_y_vs_x_proton_plus_accepted_020", "pos_y_vs_x_proton_plus" , 200, -0.05, 0.05, 200, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_accepted_120"] = new TH2F("pos_y_vs_x_proton_minus_accepted_120", "pos_y_vs_x_proton_minus" , 200, -0.05, 0.05, 200, -0.05, 0.05);

  histosTH2F["pos_y_vs_x_proton_plus_024_025"] = new TH2F("pos_y_vs_x_proton_plus_024_025", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_plus_024_025_accept"] = new TH2F("pos_y_vs_x_proton_plus_024_025_accept", "pos_y_vs_x_proton_plus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_124_125"] = new TH2F("pos_y_vs_x_proton_minus_124_125", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);
  histosTH2F["pos_y_vs_x_proton_minus_124_125_accept"] = new TH2F("pos_y_vs_x_proton_minus_124_125_accept", "pos_y_vs_x_proton_minus" , 400, -0.05, 0.05, 400, -0.05, 0.05);

  RooUnfoldResponse response_xminus (bin2, -4, 0,"unfolded_x_minus","unfolded_x_minus");
  histosTH2F["response_xminus_jj"] = (TH2F*) response_xminus.Hresponse();
  RooUnfoldResponse response_xplus (bin2, -4, 0,"unfolded_x_plus","unfolded_x_plus");
  histosTH2F["response_xplus_jj"] = (TH2F*) response_xplus.Hresponse();

  TH1F* event_selection = new TH1F("event_selection", "event_selection", 9, 0, 9);

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  gStyle->SetPalette(1);

  TH1F* rdm_evt_zb = new TH1F("rdm_evt_zb","rdm_evt_zb",50,0,10000000);
  //===================
  int i_tot = 0 , nevt_tot = 0;

  //MC Files
  const char *dirname="/storage1/lhuertas/Analysis/CMSTotem/MC/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/";
  const char *ext=".root";
  vector<TString>* vfiles = new vector<TString>; 
  
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

  //Declaration of tree and its branches variables
  TTree* tree = new TTree(treeName.c_str(),"");
  MyEvtId*           evtId        = NULL;
//   MyL1TrigOld*       l1Trig       = NULL;  
//   MyHLTrig*          hltTrig      = NULL;
  vector<MyGenPart>* genPart      = NULL;
  vector<MyTracks>*  track_coll   = NULL;
  vector<MyVertex>*  vertex_coll  = NULL;
  vector<MyPFJet>*   pfJet_coll   = NULL;
  vector<MyPFCand>*  pFlow_coll   = NULL;
  vector<MyGenJet>*   genJet_coll   = NULL;
  MyGenKin*  genKin   = NULL;
  //===============================

  //ZeroBias Files
  string treeNameZB = "cms_totem";
  TChain treeZB("cms_totem");
  vector<TString>* vdirs_zb = new vector<TString>;
  vdirs_zb->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/198902-8369_8371-V00-02-00/ZeroBias/");
  vdirs_zb->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/198903-8372-V00-02-00/ZeroBias/");
  vector<TString>* vfiles_zb = new vector<TString>;
  for(vector<TString>::iterator itdirs_zb = vdirs_zb->begin(); itdirs_zb != vdirs_zb->end(); ++itdirs_zb){
      TString& dirname_zb = *itdirs_zb;
      TSystemDirectory dir(dirname_zb, dirname_zb);
      TList *files = dir.GetListOfFiles();
      if (files) {
         TSystemFile *file;
         TString fname;
         TIter next(files);
         while ((file=(TSystemFile*)next())) {
             fname = file->GetName();
             if (!file->IsDirectory() && fname.EndsWith(ext)) {
                 TString root_file_zb = dirname_zb + string(fname.Data());
                 vfiles_zb->push_back(root_file_zb); cout<<root_file_zb<<endl;
                 treeZB.Add(root_file_zb);
             }
         }
      }
   }//cout<<"n_events_total ZB "<<treeZB.GetEntries()<<endl;

  //TTree* treeZB = new TTree(treeNameZB.c_str(),"");
  MyEvtId*           evtIdZB        = NULL;
  RPRootDumpReconstructedProton* rec_proton_left  = NULL;
  RPRootDumpReconstructedProton* rec_proton_right = NULL;
  map<unsigned int, RPRootDumpTrackInfo*> rp_track_info;
  vector<MyVertex>*  vertex_coll_ZB  = NULL;
  vector<MyPFCand>*  pFlow_coll_ZB   = NULL;
  vector<MyPFJet>*   pfJet_coll_ZB   = NULL;
//  TFile* fileZB = TFile::Open(fileNameZB.c_str(),"READ");
    
//  treeZB = (TTree*)fileZB->Get( treeNameZB.c_str() );
  int nevZB = int(treeZB.GetEntries());
  
  treeZB.SetBranchAddress("cmsEvtUA",&evtIdZB);
  treeZB.SetBranchAddress("rec_prot_left.",&rec_proton_left);
  treeZB.SetBranchAddress("rec_prot_right.",&rec_proton_right);
  treeZB.SetBranchAddress("cmsVerticesUA",&vertex_coll_ZB);
  treeZB.SetBranchAddress("cmsParticleFlowUA",&pFlow_coll_ZB);
  treeZB.SetBranchAddress("cmsak5PFJetsUA",&pfJet_coll_ZB);
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
          treeZB.SetBranchAddress(br_name, &rp_track_info[id]);
       }
  }

  std::vector<MyZeroBiasData> zeroBias; 
//  int nevt_max_ZB = 500000;
  for(int i_evt = 0; i_evt < nevZB && i_evt < nevt_max_corr; ++i_evt){
    
     if( ((i_evt+1) % 10000) == 0) cout <<int(double(i_evt+1)/1000)<<"k done"<<endl;

     //Filling the variables defined setting branches
     treeZB.GetEntry(i_evt);

     MyZeroBiasData mydata;
     mydata.proton_rec_left_valid = rec_proton_left->valid;
     mydata.proton_rec_left_t = rec_proton_left->t;
     mydata.proton_rec_left_xi = rec_proton_left->xi;
     mydata.proton_rec_right_valid = rec_proton_right->valid;
     mydata.proton_rec_right_t = rec_proton_right->t;
     mydata.proton_rec_right_xi = rec_proton_right->xi;
     
     MyVertex& primaryVertex = vertex_coll_ZB->at(0);
     /*for(vector<MyVertex>::iterator it_vtx_zb = vertex_coll_ZB->begin() ; it_vtx_zb != vertex_coll_ZB->end() ; ++it_vtx_zb){
	mydata.vtx_ndof = it_vtx_zb->ndof;    
        mydata.vtx_valid = it_vtx_zb->validity;
     }*/
     mydata.vtx_ndof = primaryVertex.ndof;
     mydata.vtx_valid = primaryVertex.validity;
     mydata.vtx_x = primaryVertex.x;
     mydata.vtx_y = primaryVertex.y;
     mydata.vtx_z = primaryVertex.z;

     if( pfJet_coll_ZB->size() > 0 ){
         MyBaseJet const& leadingJet = ( pfJet_coll_ZB->at(0) ).mapjet["ak5PFL2L3Residual"];
         mydata.leadingJet_pt = leadingJet.Pt();
         mydata.leadingJet_eta = leadingJet.Eta();
     }
     if( pfJet_coll_ZB->size() > 1 ){
         MyBaseJet const& secondJet = ( pfJet_coll_ZB->at(1) ).mapjet["ak5PFL2L3Residual"];
         mydata.secondJet_pt = secondJet.Pt();
         mydata.secondJet_eta = secondJet.Eta();
     }

     mydata.rp_track_valid_020 = rp_track_info[20]->valid;
     mydata.rp_track_valid_021 = rp_track_info[21]->valid;
     mydata.rp_track_valid_024 = rp_track_info[24]->valid;
     mydata.rp_track_valid_025 = rp_track_info[25]->valid;
     mydata.rp_track_valid_120 = rp_track_info[120]->valid;
     mydata.rp_track_valid_121 = rp_track_info[121]->valid;
     mydata.rp_track_valid_124 = rp_track_info[124]->valid;
     mydata.rp_track_valid_125 = rp_track_info[125]->valid;
     mydata.rp_x_024 = rp_track_info[24]->x;
     mydata.rp_y_024 = rp_track_info[24]->y;
     mydata.rp_x_025 = rp_track_info[25]->x;
     mydata.rp_y_025 = rp_track_info[25]->y;
     mydata.rp_x_124 = rp_track_info[124]->x;
     mydata.rp_y_124 = rp_track_info[124]->y;
     mydata.rp_x_125 = rp_track_info[125]->x;
     mydata.rp_y_125 = rp_track_info[125]->y;

     double sum_plus = 0;
     double sum_minus = 0;
     for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll_ZB->begin(); it_pfcand != pFlow_coll_ZB->end(); ++it_pfcand){
         int partType = it_pfcand->particleId;
         double eta = it_pfcand->Eta();
         double energy = it_pfcand->Energy();
         double pz = it_pfcand->Pz();         
         sum_plus += (energy + pz);
         sum_minus += (energy - pz);
     }
     mydata.xi_cms_plus = sum_plus/8000; 
     mydata.xi_cms_minus = sum_minus/8000; 

     zeroBias.push_back(mydata);
  }
  //fileZB->Close();

  cout << zeroBias.size() << " events analyzed" << endl;
  
 // for(vector<MyTOTEMData>::const_iterator it_zb = zeroBias.begin(); it_zb != zeroBias.end(); it_zb++){
 //    cout<<it_zb->proton_rec_right_t<<endl;  
 // }

  //================================================


  rp_aperture_config();

  gRandom->SetSeed(12345);

  double vtx = 0;
  double jet_gen = 0;
  double vtx_jetrec = 0;
  double vtx_jetrec_rp = 0;
  double vtx_jetrec_backg = 0;
  double vtx_jetrec_rp_proton_cut = 0;
  double vtx_jetrec_backg_proton_cut = 0;
  double nevents_total = 0; 
  double nevents_pf = 0;

  double nweight_total = 0; 
  double weight_total_leadingJet = 0; 
  double weight_total_secondJet = 0; 
  double weight_total_leadingJet_selected = 0; 
  double weight_total_secondJet_selected = 0; 
  double weight_total_Jet_selected = 0; 
  double weight_total_PF_selected = 0; 
  //starting Loop over files, stops at end of list of files or when reached nevt_max
  for(vector<TString>::iterator itfiles = vfiles->begin(); itfiles != vfiles->end() && i_tot < nevt_max_corr; ++itfiles){
  
    TFile* file = TFile::Open(*itfiles,"READ");
    
    //getting the tree form the current file
    tree = (TTree*) file->Get( treeName.c_str() );

    //Getting number of events
    int nev = int(tree->GetEntriesFast());
    nevt_tot += nev;
    cout <<"The current file has " << nev << " entries : "<< endl << *itfiles << endl;

    //adding branches to the tree ----------------------------------------------------------------------
    tree->SetBranchAddress("evtId",&evtId);
    tree->SetBranchAddress("generalTracks",&track_coll); 
    tree->SetBranchAddress("offlinePrimaryVertices",&vertex_coll);
    tree->SetBranchAddress("ak5PFJets",&pfJet_coll);
    tree->SetBranchAddress("particleFlow",&pFlow_coll);
    tree->SetBranchAddress("genKin",&genKin);
    tree->SetBranchAddress("genPart",&genPart);
    tree->SetBranchAddress("ak5GenJets",&genJet_coll);
  
    //starting loop over events, stops when reached end of file or nevt_max
    for(int i_evt = 0; i_evt < nev && i_tot < nevt_max_corr; ++i_evt , ++i_tot){
    
    //printing the % of events done every 10k evts
    if( ((i_tot+1) % 10000) == 0) cout <<int(double(i_tot+1)/1000)<<"k done"<<endl;
    
      //Filling the variables defined setting branches
      tree->GetEntry(i_evt);

      bool passedHLT = false;
      bool passedvtx = false;
      bool jet1_selected = false;
      bool jet2_selected = false;
      bool pz_proton_max = false;
      bool PF_eta_max = false;
      bool PF_eta_min = false;
      bool xi_negat_gen = false;
      bool xi_posit_gen = false;

      weight_st = -1.; weight_st_back = -1.;
      xi_cms_st = -999.; xi_totem_st = -999.; xi_cms_minus_totem_st = -999.;
      xi_cms_st_back = -999.; xi_totem_st_back= -999.; xi_cms_minus_totem_st_back = -999.;
      
      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = genKin->genWeight; 
  //    double event_weight = 1.0;
      nweight_total += event_weight; 
      ++nevents_total;
      
      weight_st = event_weight;
      weight_st_back = event_weight;
      //-------------------------------------------------------------------------------------------------
 	    
     // Vertices
      for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
      //  if (it_vtx!=vertex_coll->begin()) continue;
	if( it_vtx->ndof>4 ) passedvtx = true;   
        histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
        histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
        histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );
      }
      //if(!passedvtx) continue;
      MyVertex& primaryVertex = vertex_coll->at(0);
      double prim_vtx_r = sqrt( primaryVertex.x*primaryVertex.x + primaryVertex.y*primaryVertex.y );
      bool select_Vertex = ( !primaryVertex.fake && primaryVertex.validity && primaryVertex.ndof > 4);// && fabs( primaryVertex.z ) < 15.0 && prim_vtx_r < 2.0);
      if (select_Vertex) ++vtx;
      
      // Tracks
      int n_tracks_selected = 0;
      for(vector<MyTracks>::iterator it_trk = track_coll->begin() ; it_trk != track_coll->end() ; ++it_trk){
         histosTH1F["track_pt"]->Fill( it_trk->Pt(), event_weight );
         histosTH1F["track_eta"]->Fill( it_trk->Eta(), event_weight );
         histosTH1F["track_phi"]->Fill( it_trk->Phi(), event_weight );

         ++n_tracks_selected;
      }
      histosTH1F["track_multiplicity"]->Fill( n_tracks_selected, event_weight );


      //Jet generated level         
      double leadingJet_pt_gen = -999;
      double secondJet_pt_gen = -999;
      double Jet1_energy_gen;
      double Jet1_px_gen;
      double Jet1_py_gen;
      double Jet1_pz_gen;
      double Jet2_energy_gen;
      double Jet2_px_gen;
      double Jet2_py_gen;
      double Jet2_pz_gen;
      double Jet1_eta_gen;
      double Jet1_phi_gen;
      double Jet2_eta_gen;
      double Jet2_phi_gen;
      double jet_pt_gen;
      double jet_eta_gen;
      double jet_ene_gen;
      double jet_px_gen;
      double jet_py_gen;
      double jet_pz_gen;
      double jet_phi_gen;

         
      for(vector<MyGenJet>::iterator it_genjet = genJet_coll->begin(); it_genjet != genJet_coll->end(); ++it_genjet){
         jet_pt_gen = it_genjet->Pt();
         jet_eta_gen = it_genjet->Eta();
         jet_ene_gen = it_genjet->E();
         jet_px_gen = it_genjet->Px();
         jet_py_gen = it_genjet->Py();
         jet_pz_gen = it_genjet->Pz();
         jet_phi_gen = it_genjet->Phi();
         //cout<<i_evt<<" "<<pfJet_coll->size();
        // if (fabs(jet_eta_gen)>4.4) continue;
         
         if (jet_pt_gen>leadingJet_pt_gen){
             leadingJet_pt_gen = jet_pt_gen;
             Jet1_energy_gen = jet_ene_gen;
             Jet1_px_gen = jet_px_gen;
             Jet1_py_gen = jet_py_gen;
             Jet1_pz_gen = jet_pz_gen;
             Jet1_eta_gen = jet_eta_gen;
             Jet1_phi_gen = jet_phi_gen;
         }
         if (jet_pt_gen>secondJet_pt_gen && jet_pt_gen<leadingJet_pt_gen){
             secondJet_pt_gen = jet_pt_gen;
             Jet2_energy_gen = jet_ene_gen;
             Jet2_px_gen = jet_px_gen;
             Jet2_py_gen = jet_py_gen;
             Jet2_pz_gen = jet_pz_gen;
             Jet2_eta_gen = jet_eta_gen;
             Jet2_phi_gen = jet_phi_gen;
         }

      }
      //if(leadingJet_pt_gen < 30.) continue;
      //if(secondJet_pt_gen < 30.) continue;
      bool jet1_selected_gen = leadingJet_pt_gen > 30 && fabs(Jet1_eta_gen)<4.4;
      bool jet2_selected_gen = secondJet_pt_gen > 30 && fabs(Jet2_eta_gen)<4.4;
      bool selectJet_gen = jet1_selected_gen && jet2_selected_gen;

      double mass_jets_gen= sqrt(pow(Jet1_energy_gen+Jet2_energy_gen,2)-pow(Jet1_px_gen+Jet2_px_gen,2)-pow(Jet1_py_gen+Jet2_py_gen,2)-pow(Jet1_pz_gen+Jet2_pz_gen,2));
      double x_plus_gen = ((Jet1_energy_gen+Jet1_pz_gen)+(Jet2_energy_gen+Jet2_pz_gen))/8000;
      double x_minus_gen = ((Jet1_energy_gen-Jet1_pz_gen)+(Jet2_energy_gen-Jet2_pz_gen))/8000;
      double x_gen_selec = (x_minus_gen<x_plus_gen) ? x_minus_gen : x_plus_gen;
 

      //Jets with pt>30Gev and !eta!<2
      Double_t Jet1_pt; 
      Double_t Jet2_pt; 
      Double_t Jet1_eta; 
      Double_t Jet2_eta; 
      Double_t Jet1_phi, Jet1_px, Jet1_py, Jet1_pz, Jet1_energy; 
      Double_t Jet2_phi, Jet2_px, Jet2_py, Jet2_pz, Jet2_energy;
      
      for(vector<MyPFJet>::iterator it_jet = pfJet_coll->begin() ; it_jet != pfJet_coll->end() ; ++it_jet){
         map<string,MyBaseJet>::iterator it_map = it_jet->mapjet.begin();
         for(; it_map != it_jet->mapjet.end(); ++it_map)
            if(verbose) cout << it_map->first << endl;
         MyBaseJet const& basejet = it_jet->mapjet["ak5PFJets"];
         histosTH1F["jet_pt"]->Fill( basejet.Pt(), event_weight  );

         if(basejet.Pt() > 0.) histosTH1F["jet_eta"]->Fill( basejet.Eta(), event_weight  );
         histosTH1F["jet_phi"]->Fill( basejet.Phi(), event_weight  );
      }
      if( pfJet_coll->size() > 0 ){
	 MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet["ak5PFJets"];
         weight_total_leadingJet += event_weight;
	 Jet1_pt = leadingJet.Pt(); 
	 Jet1_eta = leadingJet.Eta(); 
	 Jet1_phi = leadingJet.Phi(); 
	 Jet1_px = leadingJet.Px(); 
	 Jet1_py = leadingJet.Py(); 
	 Jet1_pz = leadingJet.Pz(); 
	 Jet1_energy = leadingJet.E(); 
	 histosTH1F["leadingJet_pt"]->Fill( leadingJet.Pt(), event_weight  );
	 histosTH1F["leadingJet_eta"]->Fill( leadingJet.Eta(), event_weight  );
	 histosTH1F["leadingJet_phi"]->Fill( leadingJet.Phi(), event_weight  );
	 
 	 if(leadingJet.Pt() > 30. && fabs(leadingJet.Eta())<4.4 ){ 
           jet1_selected = true;
	    weight_total_leadingJet_selected += event_weight;}
      }
      
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet["ak5PFJets"];
         weight_total_secondJet += event_weight;
         Jet2_pt = secondJet.Pt(); 
	 Jet2_eta = secondJet.Eta(); 
	 Jet2_phi = secondJet.Phi(); 
         Jet2_px = secondJet.Px();
         Jet2_py = secondJet.Py();
         Jet2_pz = secondJet.Pz();
         Jet2_energy = secondJet.E();
         histosTH1F["secondJet_pt"]->Fill( secondJet.Pt(), event_weight  );  
	 histosTH1F["secondJet_eta"]->Fill( secondJet.Eta(), event_weight  );
	 histosTH1F["secondJet_phi"]->Fill( secondJet.Phi(), event_weight  );
	 
 	 if(secondJet.Pt() > 30. && fabs(secondJet.Eta())<4.4 ){  
            jet2_selected = true;
 	    weight_total_secondJet_selected += event_weight;}
      }
      bool selectJet_rec = jet1_selected && jet2_selected;
      //if( !selectJet) continue;
      double mass_jets= sqrt(pow(Jet1_energy+Jet2_energy,2)-pow(Jet1_px+Jet2_px,2)-pow(Jet1_py+Jet2_py,2)-pow(Jet1_pz+Jet2_pz,2));
      double x_plus = ((Jet1_energy+Jet1_pz)+(Jet2_energy+Jet2_pz))/8000;
      double x_minus = ((Jet1_energy-Jet1_pz)+(Jet2_energy-Jet2_pz))/8000;
      double x_selec = (x_minus<x_plus) ?  x_minus : x_plus;

      weight_total_Jet_selected += event_weight;
//       histosTH1F["DeltaPhiJet"]->Fill( deltaphi, event_weight  );	




      //GenPart
      double genEPlusPz = 0;
      double genEMinusPz = 0;
      double cm = 8000;
      double proton_pi = 4000;
      double proton_pz_plus=-999;
      double proton_px_plus = -999.;
      double proton_py_plus = -999.;
      double proton_energy_plus = 0.;
      double proton_pz_minus=999;
      double proton_px_minus = 999.;
      double proton_py_minus = 999.;
      double proton_energy_minus = 0.;
      double px_gen;
      double py_gen;
      double pz_gen;
      double energy_gen;
      double proton_pf;
      
      for(vector<MyGenPart>::iterator it_genpart = genPart->begin(); it_genpart != genPart->end(); ++it_genpart){
 
	 double eta_gen = it_genpart->Eta();
         int status = it_genpart->status;
         int id = it_genpart->pdgId;

         if (status == 1 && id == 2212) {
            energy_gen = it_genpart->Energy();
            px_gen = it_genpart->Px();
            py_gen = it_genpart->Py();
            pz_gen = it_genpart->Pz();
            proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
            double pz_cut = 0.7*proton_pi;
            if (fabs(pz_gen) > pz_cut){
                genEPlusPz += (energy_gen + pz_gen);
                genEMinusPz += (energy_gen - pz_gen);

                if (pz_gen > proton_pz_plus) {
                    proton_pz_plus = pz_gen; proton_energy_plus = energy_gen;
                    proton_px_plus = px_gen; proton_py_plus = py_gen;
                }
                if (pz_gen < proton_pz_minus) {
                   proton_pz_minus = pz_gen; proton_energy_minus = energy_gen;
                   proton_px_minus = px_gen; proton_py_minus = py_gen;
                }
            }
         }
      }
	 

      double xi_minus_gen = genEMinusPz/cm;
      double xi_proton_plus = -1.;
      double xi_proton_minus = -1.;
      double t_proton_plus = 0.;
      double t_proton_minus = 0.;
      double thx_proton_plus = 0.;
      double thy_proton_plus = 0.;
      double thx_proton_minus = 0.;
      double thy_proton_minus = 0.;


      bool proton_minus_rp_accept_120 = false;
      bool proton_minus_rp_accept_121 = false;
      bool proton_minus_rp_accept_122 = false;
      bool proton_minus_rp_accept_123 = false;
      bool proton_minus_rp_accept_124 = false;
      bool proton_minus_rp_accept_125 = false;
      bool proton_minus_rp_accept_020 = false;

      bool proton_plus_rp_accept_020 = false;
      bool proton_plus_rp_accept_021 = false;
      bool proton_plus_rp_accept_022 = false;
      bool proton_plus_rp_accept_023 = false;
      bool proton_plus_rp_accept_024 = false;
      bool proton_plus_rp_accept_025 = false;
      bool proton_plus_rp_accept_120 = false;

      std::map<int,std::vector<double> > proton_plus_pars;
      std::map<int,std::vector<double> > proton_minus_pars;

      bool fiducial_cut_rp_024= false; 
      bool fiducial_cut_rp_025 = false;
      if(proton_pz_plus > 0.){
         xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
         //t_proton_plus = -2*( (proton_pi*proton_energy_plus) - (proton_pi*proton_pz_plus) );
         TLorentzVector vec_pi(0.,0.,proton_pi,proton_pi);
         TLorentzVector vec_pf(proton_px_plus,proton_py_plus,proton_pz_plus,proton_energy_plus);
         TLorentzVector vec_t = (vec_pf - vec_pi);
         t_proton_plus = vec_t.Mag2();
         thx_proton_plus = atan(-proton_px_plus/proton_pi);
         thy_proton_plus = atan(proton_py_plus/proton_pi);

         //FIXME
         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_plus_rp_accept_020 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 20, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[20] = std::vector<double>(5,0.);
         proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
         proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
         proton_plus_pars[20][4] = out_xi;

         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[24] = std::vector<double>(5,0.);
         proton_plus_pars[24][0] = out_x; proton_plus_pars[24][1] = out_y;
         proton_plus_pars[24][2] = out_thx; proton_plus_pars[24][3] = out_thy;
         proton_plus_pars[24][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_plus_024_025"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );
         fiducial_cut_rp_024 = proton_plus_pars[24][0]>0 && proton_plus_pars[24][0]<0.006 && proton_plus_pars[24][1]>0.0084 && proton_plus_pars[24][1]<0.029;

         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[25] = std::vector<double>(5,0.);
         proton_plus_pars[25][0] = out_x; proton_plus_pars[25][1] = out_y;
         proton_plus_pars[25][2] = out_thx; proton_plus_pars[25][3] = out_thy;
         proton_plus_pars[25][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_plus_024_025"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );
         fiducial_cut_rp_025 = proton_plus_pars[25][0]>0 && proton_plus_pars[25][0]<0.006 && proton_plus_pars[25][1]<-0.0084 && proton_plus_pars[25][1]>-0.029;

         proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21);
         proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
         proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
         proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);
         if (fiducial_cut_rp_024 || fiducial_cut_rp_025){
             histosTH2F["pos_y_vs_x_proton_plus_024_025_accept"]->Fill( proton_plus_pars[24][0], proton_plus_pars[24][1] , event_weight );
             histosTH2F["pos_y_vs_x_proton_plus_024_025_accept"]->Fill( proton_plus_pars[25][0], proton_plus_pars[25][1] , event_weight );
         }
      }

      bool fiducial_cut_rp_124= false; 
      bool fiducial_cut_rp_125 = false;
      if(proton_pz_minus < 0.){
         xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
         // t_proton_minus = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) ); 
         TLorentzVector vec_pi(0.,0.,-proton_pi,proton_pi);
         TLorentzVector vec_pf(proton_px_minus,proton_py_minus,proton_pz_minus,proton_energy_minus);
         TLorentzVector vec_t = (vec_pf - vec_pi);
         t_proton_minus = vec_t.Mag2();
 
         thx_proton_minus = atan(-proton_px_minus/proton_pi);
         thy_proton_minus = atan(proton_py_minus/proton_pi);

         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_minus_rp_accept_120 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 120, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[120] = std::vector<double>(5,0.);
         proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
         proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
         proton_minus_pars[120][4] = out_xi;

         proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[124] = std::vector<double>(5,0.);
         proton_minus_pars[124][0] = out_x; proton_minus_pars[124][1] = out_y;
         proton_minus_pars[124][2] = out_thx; proton_minus_pars[124][3] = out_thy;
         proton_minus_pars[124][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_minus_124_125"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );
         fiducial_cut_rp_124 = proton_minus_pars[124][0]>0 && proton_minus_pars[124][0]<0.006 && proton_minus_pars[124][1]>0.0084 && proton_minus_pars[124][1]<0.027;

         proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[125] = std::vector<double>(5,0.);
         proton_minus_pars[125][0] = out_x; proton_minus_pars[125][1] = out_y;
         proton_minus_pars[125][2] = out_thx; proton_minus_pars[125][3] = out_thy;
         proton_minus_pars[125][4] = out_xi;
         histosTH2F["pos_y_vs_x_proton_minus_124_125"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );
         fiducial_cut_rp_125 = proton_minus_pars[125][0]>0 && proton_minus_pars[125][0]<0.006 && proton_minus_pars[125][1]<-0.0084 && proton_minus_pars[125][1]>-0.027;

         if (fiducial_cut_rp_124 || fiducial_cut_rp_125){
             histosTH2F["pos_y_vs_x_proton_minus_124_125_accept"]->Fill( proton_minus_pars[124][0], proton_minus_pars[124][1] , event_weight );
             histosTH2F["pos_y_vs_x_proton_minus_124_125_accept"]->Fill( proton_minus_pars[125][0], proton_minus_pars[125][1] , event_weight );
         }

         proton_minus_rp_accept_121 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 121);
         proton_minus_rp_accept_122 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 122);
         proton_minus_rp_accept_123 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 123);
         proton_minus_rp_accept_124 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 124);
         proton_minus_rp_accept_125 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 125);
         proton_minus_rp_accept_020 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 20);
      }


      // Particle-flow
      double soma1 = 0;
      double soma2 = 0;
      double eta_max=-999.;
      double eta_min=999.;
      
      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 double pz = it_pfcand->Pz();

	 // HF eta rings 29, 30, 40, 41
         if( ( (fabs(eta) <= 2.866) && (fabs(eta) > 3.152) ) || (fabs(eta) <= 4.730) ){

	    soma1 += (energy + pz);
	    soma2 += (energy - pz);

            if (eta > eta_max) { eta_max = eta; PF_eta_max = true;} 
	    if (eta < eta_min) { eta_min = eta; PF_eta_min = true;}

  
       	    histosTH2F["energyVsEtaAllTypes"]->Fill( eta, energy, event_weight );
	   
    	    if(partType == MyPFCand::X)
	     histosTH2F["energyVsEtaUndefined"]->Fill( eta, energy, event_weight );
	    else if(partType == MyPFCand::h)
	     histosTH2F["energyVsEtaChargedHadron"]->Fill( eta, energy, event_weight ); 
	    else if(partType == MyPFCand::e) 
	     histosTH2F["energyVsEtaElectron"]->Fill( eta, energy, event_weight );
	    else if(partType == MyPFCand::mu) 
	     histosTH2F["energyVsEtaMuon"]->Fill( eta, energy, event_weight );
	    else if(partType == MyPFCand::gamma) 
	     histosTH2F["energyVsEtaGamma"]->Fill( eta, energy, event_weight );
	    else if(partType == MyPFCand::h0) 
	     histosTH2F["energyVsEtaNeutralHadron"]->Fill( eta, energy, event_weight );
	    else if(partType == MyPFCand::h_HF){ 
	     histosTH2F["energyVsEtaHadronHF"]->Fill( eta, energy, event_weight );}
	    else if(partType == MyPFCand::egamma_HF) 
	     histosTH2F["energyVsEtaEGammaHF"]->Fill( eta, energy, event_weight );
         }
       }
       if(PF_eta_max && PF_eta_min) ++nevents_pf;
       weight_total_PF_selected += event_weight;


       double xi_plus_Reco = soma1/cm;
       double xi_minus_Reco = soma2/cm;
       double delta_eta_maxmin = eta_max - eta_min;  
       histosTH1F["Eta_max"]->Fill(eta_max, event_weight);
       histosTH1F["Eta_min"]->Fill(eta_min, event_weight);
       histosTH1F["Delta_eta_maxmin"]->Fill(delta_eta_maxmin, event_weight);
       
       //double correction = xi_plus_Reco/xi_plus_gen;
       //double resolution_before = (xi_plus_gen-xi_plus_Reco)/xi_plus_gen;
       //double xi_reconst = xi_plus_Reco/0.8;
       //double resolution_after = (xi_plus_gen-xi_reconst)/xi_plus_gen;

     //totem proton reconstructed
      float sigma_xi56 = 0.00720615 - 0.0418783*xi_proton_minus + 0.0999515*xi_proton_minus*xi_proton_minus; // sigma56 vs xi from Hubert
      float xi_proton_minus_rec = xi_proton_minus + gRandom->Gaus(0,sigma_xi56);
      double sigma_t56 = 0.233365*t_proton_minus - 0.0975751*t_proton_minus*t_proton_minus;  // sigma_t56 vs t from Hubert
      double t_proton_minus_rec = t_proton_minus + gRandom->Gaus(0,sigma_t56);
      float sigma_xi45=0.00714986 - 0.0408903*xi_proton_plus + 0.0965813*xi_proton_plus*xi_proton_plus;   // sigma45 vs xi from Hubert
      float xi_proton_plus_rec = xi_proton_plus + gRandom->Gaus(0,sigma_xi45);
      double sigma_t45=0.233365*t_proton_plus - 0.0975751*t_proton_plus*t_proton_plus;     // sigma_t45 vs t from Hubert
      double t_proton_plus_rec = t_proton_plus + gRandom->Gaus(0,sigma_t45);
      double proton_minus_beta = x_minus/xi_proton_minus_rec;
      double proton_plus_beta = x_plus/xi_proton_plus_rec;

  //    MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet["ak5PFJets"];
  //    MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet["ak5PFJets"];
      Double_t deltapt = abs(Jet1_pt - Jet2_pt);
      Double_t deltaeta = abs(Jet1_eta - Jet2_eta);
      Double_t deltaphi = abs(Jet1_phi - Jet2_phi);
      Double_t eta_average = 0.5*(Jet1_eta + Jet2_eta);


      //rp_accept
      bool proton_minus_rp_accept_mc = ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 && fiducial_cut_rp_124) || ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 && fiducial_cut_rp_125 );// || ( proton_minus_rp_accept_122 && proton_minus_rp_accept_123 );
      bool proton_plus_rp_accept_mc = ( proton_plus_rp_accept_020 && proton_plus_rp_accept_024 && fiducial_cut_rp_024 ) || ( proton_plus_rp_accept_021 && proton_plus_rp_accept_025 && fiducial_cut_rp_025 );// || ( proton_plus_rp_accept_022 && proton_plus_rp_accept_023 );
      bool proton_minus_kinec_accept_t_mc = fabs(t_proton_minus_rec)>=0.03 && fabs(t_proton_minus_rec)<=1;
      bool proton_minus_kinec_accept_xi_mc = xi_proton_minus_rec>=0.03 && xi_proton_minus_rec<=0.1;
      bool proton_plus_kinec_accept_t_mc = fabs(t_proton_plus_rec)>=0.03 && fabs(t_proton_plus_rec)<=1;
      bool proton_plus_kinec_accept_xi_mc = xi_proton_plus_rec>=0.03 && xi_proton_plus_rec<=0.1;


      //Access Zero-bias data
      //gRandom->SetSeed(12345);
      int i_evt_ZB = 0 + gRandom->Rndm()*(zeroBias.size());
      //rdm_evt_zb->Fill(i_evt_ZB, 1.);

      MyZeroBiasData const & zeroBiasData = zeroBias.at(i_evt_ZB);
      double xi_proton_right_zb = zeroBiasData.proton_rec_right_xi;      
      double t_proton_right_zb = zeroBiasData.proton_rec_right_t;
      bool valid_proton_right_zb = zeroBiasData.proton_rec_right_valid;      
      bool valid_proton_left_zb = zeroBiasData.proton_rec_left_valid;      
      double xi_proton_left_zb = zeroBiasData.proton_rec_left_xi;      
      double t_proton_left_zb = zeroBiasData.proton_rec_left_t;
      bool valid_vtx = zeroBiasData.vtx_valid;
      double ndof_vtx = zeroBiasData.vtx_ndof;
      bool vtx_zb = valid_vtx && ndof_vtx<4;
      //if (valid_vtx && ndof_vtx>4) continue;

      //fiducial cuts
      double rp_x_024_zb = zeroBiasData.rp_x_024;
      double rp_y_024_zb = zeroBiasData.rp_y_024;
      double rp_x_025_zb = zeroBiasData.rp_x_025;
      double rp_y_025_zb = zeroBiasData.rp_y_025;
      double rp_x_124_zb = zeroBiasData.rp_x_124;
      double rp_y_124_zb = zeroBiasData.rp_y_124;
      double rp_x_125_zb = zeroBiasData.rp_x_125;
      double rp_y_125_zb = zeroBiasData.rp_y_125;
      bool cut_rp_024_zb =  rp_x_024_zb>0 && rp_x_024_zb<6 && rp_y_024_zb>8.4 && rp_y_024_zb<29 ;
      bool cut_rp_025_zb =  rp_x_025_zb>0 && rp_x_025_zb<6 && rp_y_025_zb<-8.4 && rp_y_024_zb>-29 ;
      bool cut_rp_124_zb =  rp_x_124_zb>0 && rp_x_124_zb<6 && rp_y_124_zb>8.4 && rp_y_124_zb<27 ;
      bool cut_rp_125_zb =  rp_x_125_zb>0 && rp_x_125_zb<6 && rp_y_125_zb<-8.4 && rp_y_124_zb>-27 ;
      bool rp_track_valid_120_zb = zeroBiasData.rp_track_valid_120;
      bool rp_track_valid_020_zb = zeroBiasData.rp_track_valid_020;
      bool rp_track_valid_121_zb = zeroBiasData.rp_track_valid_121;
      bool rp_track_valid_021_zb = zeroBiasData.rp_track_valid_021;
      bool rp_track_valid_124_zb = zeroBiasData.rp_track_valid_124;
      bool rp_track_valid_024_zb = zeroBiasData.rp_track_valid_024;
      bool rp_track_valid_125_zb = zeroBiasData.rp_track_valid_125;
      bool rp_track_valid_025_zb = zeroBiasData.rp_track_valid_025;
      bool rp_track_accept_right_zb = ( rp_track_valid_120_zb && rp_track_valid_124_zb && cut_rp_124_zb )  || ( rp_track_valid_121_zb && rp_track_valid_125_zb && cut_rp_125_zb );
      bool rp_track_accept_left_zb = ( rp_track_valid_020_zb && rp_track_valid_024_zb && cut_rp_024_zb )  || ( rp_track_valid_021_zb && rp_track_valid_025_zb && cut_rp_025_zb );


      double leadingjet_pt_zb = zeroBiasData.leadingJet_pt;
      double leadingjet_eta_zb = zeroBiasData.leadingJet_eta;
      double secondjet_pt_zb = zeroBiasData.secondJet_pt;
      double secondjet_eta_zb = zeroBiasData.secondJet_eta;
      double proton_minus_beta_zb = x_minus/-xi_proton_right_zb;
      double proton_plus_beta_zb = x_plus/-xi_proton_left_zb;
      double xi_cms_plus_zb = zeroBiasData.xi_cms_plus;
      double xi_cms_minus_zb = zeroBiasData.xi_cms_minus;
      double xi_cms_minus = xi_cms_minus_zb+xi_minus_Reco;
      double xi_cms_plus = xi_cms_plus_zb+xi_plus_Reco;
      bool proton_right_kinec_xi_zb = -xi_proton_right_zb>=0.03 && -xi_proton_right_zb<=0.1;
      bool proton_right_kinec_t_zb = fabs(t_proton_right_zb)>=0.03 && fabs(t_proton_right_zb)<=1; 
      bool xi_cms_totem_cut_right_zb = xi_cms_minus+xi_proton_right_zb<=0;      
      bool proton_left_kinec_xi_zb = -xi_proton_left_zb>=0.03 && -xi_proton_left_zb<=0.1;
      bool proton_left_kinec_t_zb = -t_proton_left_zb>=0.03 && -t_proton_left_zb<=1; 
      bool xi_cms_totem_cut_left_zb = xi_cms_plus+xi_proton_left_zb<=0;      

      //decision where the proton in Totem acceptance came from
      bool xi_cms_totem_cut_minus =  xi_cms_minus-xi_proton_minus_rec<=0;
      bool xi_cms_totem_cut_plus =  xi_cms_plus-xi_proton_plus_rec<=0;
      bool signal_right = !valid_proton_right_zb && !valid_proton_left_zb && !proton_plus_rp_accept_mc && proton_minus_rp_accept_mc;
      bool signal_left = !valid_proton_right_zb && !valid_proton_left_zb && !proton_minus_rp_accept_mc && proton_plus_rp_accept_mc;
      bool backg_right = !proton_plus_rp_accept_mc && !proton_minus_rp_accept_mc && !valid_proton_left_zb && valid_proton_right_zb && rp_track_accept_right_zb && vtx_zb;
      bool backg_left = !proton_plus_rp_accept_mc && !proton_minus_rp_accept_mc && !valid_proton_right_zb && valid_proton_left_zb && rp_track_accept_left_zb && vtx_zb;
      bool rec_selection = false;
      bool gen_selection = false;

      if (selectJet_gen){
          ++jet_gen;
          gen_selection = true; 
          histosTH1F["leadingJet_pt_gen"]->Fill(leadingJet_pt_gen, event_weight);
          histosTH1F["secondJet_pt_gen"]->Fill(secondJet_pt_gen, event_weight);
          histosTH1F["leadingJet_eta_gen"]->Fill(Jet1_eta_gen, event_weight);
          histosTH1F["secondJet_eta_gen"]->Fill(Jet2_eta_gen, event_weight);
          histosTH1F["mass_jets_gen"]->Fill(mass_jets_gen, event_weight);
          histosTH1F["log_x_parton_gen_selec_jj"]->Fill(log10(x_gen_selec), event_weight);
          histosTH1F["log_x_parton_gen_minus_jj"]->Fill(log10(x_minus_gen), event_weight);
          histosTH1F["log_x_parton_gen_plus_jj"]->Fill(log10(x_plus_gen), event_weight);
      }

      if (select_Vertex && selectJet_rec){
         ++vtx_jetrec;
         rec_selection = true;
         histosTH1F["leadingJet_pt_rec"]->Fill(Jet1_pt, event_weight);
         histosTH1F["secondJet_pt_rec"]->Fill(Jet2_pt, event_weight);
         histosTH1F["leadingJet_eta_rec"]->Fill(Jet1_eta, event_weight);
         histosTH1F["secondJet_eta_rec"]->Fill(Jet2_eta, event_weight);
         histosTH1F["mass_jets_rec"]->Fill(mass_jets, event_weight);
         histosTH1F["log_x_parton_rec_minus_jj"]->Fill(log10(x_minus), event_weight);
         histosTH1F["log_x_parton_rec_plus_jj"]->Fill(log10(x_plus), event_weight);
         histosTH1F["log_x_parton_rec_selec_jj"]->Fill(log10(x_selec), event_weight);

         if (selectJet_gen){
            response_xminus.Fill (log10(x_minus), log10(x_minus_gen), 1.);
            response_xplus.Fill (log10(x_plus), log10(x_plus_gen), 1.);
            histosTH1F["log_x_parton_gen_minus_jj_recsel"]->Fill(log10(x_minus_gen), event_weight);
            histosTH1F["log_x_parton_gen_plus_jj_recsel"]->Fill(log10(x_plus_gen), event_weight);

         }

      }
     
      //Response Matrix 
      if ( !rec_selection && gen_selection ) response_xminus.Miss (log10(x_minus_gen), 1.);
      if ( rec_selection && !gen_selection ) response_xminus.Fake (log10(x_minus), 1.);

      if ( !rec_selection && gen_selection ) response_xplus.Miss (log10(x_plus_gen), 1.);
      if ( rec_selection && !gen_selection ) response_xplus.Fake (log10(x_plus), 1.);


      ///proton from MC in the RP acceptance
      if (select_Vertex && selectJet_rec && signal_right){
          ++vtx_jetrec_rp;
          histosTH1F["t_proton_signal_right"]->Fill( fabs(t_proton_minus_rec) , event_weight );
          histosTH1F["beta_proton_minus_signal"]->Fill( proton_minus_beta , event_weight );
          histosTH1F["xi_proton_signal_right"]->Fill( xi_proton_minus_rec , event_weight );
          histosTH1F["xi_cms_totem_right_signal"]->Fill( xi_cms_minus-xi_proton_minus_rec , event_weight );
          ///kinematic region
          if (proton_minus_kinec_accept_t_mc){
              histosTH1F["xi_proton_signal_right_kint"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["xi_proton_signal_right_kint_bin25"]->Fill( xi_proton_minus_rec , event_weight );}
          if (proton_minus_kinec_accept_t_mc && xi_proton_minus_rec<=0.1){
              histosTH1F["xi_proton_signal_right_kinxi_kint"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["xi_minus_rec_proton_signal_right_kinxi_kint"]->Fill( xi_cms_minus , event_weight );
              histosTH1F["xi_cms_totem_right_signal_kinxi_kint"]->Fill( xi_cms_minus-xi_proton_minus_rec , event_weight );
              xi_cms_st = xi_cms_minus;
              xi_totem_st = xi_proton_minus_rec;
              xi_cms_minus_totem_st = xi_cms_minus-xi_proton_minus_rec;
              small_tree->Fill();       
          }
          if (proton_minus_kinec_accept_t_mc && proton_minus_kinec_accept_xi_mc ) histosTH1F["xi_cms_totem_right_signal_kin"]->Fill( xi_cms_minus-xi_proton_minus_rec , event_weight );
          ///xi_cms -xi_totem cut 
          if (proton_minus_kinec_accept_t_mc && xi_cms_totem_cut_minus){
              histosTH1F["xi_proton_signal_right_kint_cut"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["xi_proton_signal_right_kint_cut_bin2"]->Fill( xi_proton_minus_rec , event_weight );}
          if (proton_minus_kinec_accept_t_mc && proton_minus_kinec_accept_xi_mc && xi_cms_totem_cut_minus){
              histosTH1F["xi_proton_signal_right_kin_cut"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["t_proton_signal_right_kin_cut"]->Fill( fabs(t_proton_minus_rec) , event_weight );
              histosTH1F["beta_proton_minus_signal_kin_cut"]->Fill( proton_minus_beta , event_weight );
              histosTH1F["log_x_parton_minus_signal_kin_cut"]->Fill( log10(x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_signal_kin_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_signal_kin_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_signal_kin_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["leadingJet_eta_right_signal_kin_cut"]->Fill( Jet1_eta, event_weight  );
              histosTH1F["secondJet_eta_right_signal_kin_cut"]->Fill( Jet2_eta, event_weight  );
              histosTH1F["DeltaPtJet_right_signal_kin_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_right_signal_kin_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_right_signal_kin_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_right_signal_kin_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_signal_right_kin_cut_bin2"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["beta_proton_minus_signal_kin_cut_bin2"]->Fill( proton_minus_beta , event_weight );
              histosTH1F["log_x_parton_minus_signal_kin_cut_bin2"]->Fill( log10(x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_signal_kin_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_signal_kin_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_signal_kin_cut_bin2"]->Fill( eta_average, event_weight  );
          }
          if (proton_minus_kinec_accept_t_mc && xi_proton_minus_rec<=0.1 && xi_cms_totem_cut_minus){
              ++vtx_jetrec_rp_proton_cut;
              histosTH1F["xi_minus_rec_proton_signal_right_kinxi_kint_cut"]->Fill( xi_cms_minus , event_weight );
              histosTH1F["xi_proton_signal_right_kint_kinxi_cut"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["t_proton_signal_right_kint_kinxi_cut"]->Fill( fabs(t_proton_minus_rec) , event_weight );
              histosTH1F["t_proton_signal_right_kint_kinxi_cut_bin"]->Fill( fabs(t_proton_minus_rec) , event_weight );
              histosTH1F["beta_proton_minus_signal_kint_kinxi_cut"]->Fill( proton_minus_beta , event_weight );
              histosTH1F["log_x_parton_minus_signal_kint_kinxi_cut"]->Fill( log10(x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_signal_kint_kinxi_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_signal_kint_kinxi_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_signal_kint_kinxi_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["leadingJet_eta_right_signal_kint_kinxi_cut"]->Fill( Jet1_eta, event_weight  );
              histosTH1F["secondJet_eta_right_signal_kint_kinxi_cut"]->Fill( Jet2_eta, event_weight  );
              histosTH1F["DeltaPtJet_right_signal_kint_kinxi_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_right_signal_kint_kinxi_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_right_signal_kint_kinxi_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_right_signal_kint_kinxi_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_signal_right_kint_kinxi_cut_bin2"]->Fill( xi_proton_minus_rec , event_weight );
              histosTH1F["beta_proton_minus_signal_kint_kinxi_cut_bin2"]->Fill( proton_minus_beta , event_weight );
              histosTH1F["log_x_parton_minus_signal_kint_kinxi_cut_bin2"]->Fill( log10(x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_signal_kint_kinxi_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_signal_kint_kinxi_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_signal_kint_kinxi_cut_bin2"]->Fill( eta_average, event_weight  );
         }

      }

      if (select_Vertex && selectJet_rec && signal_left){
          histosTH1F["t_proton_signal_left"]->Fill( fabs(t_proton_plus_rec) , event_weight );
          histosTH1F["beta_proton_plus_signal"]->Fill( proton_plus_beta , event_weight );
          histosTH1F["xi_proton_signal_left"]->Fill( xi_proton_plus_rec , event_weight );
          histosTH1F["xi_cms_totem_left_signal"]->Fill( xi_plus_Reco-xi_proton_plus_rec , event_weight );
          ///kinematic region
          if (proton_plus_kinec_accept_t_mc){
              histosTH1F["xi_proton_signal_left_kint"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["xi_proton_signal_left_kint_bin25"]->Fill( xi_proton_plus_rec , event_weight );}
          if (proton_plus_kinec_accept_t_mc && xi_proton_plus_rec<=0.1){
              histosTH1F["xi_proton_signal_left_kinxi_kint"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["xi_plus_rec_proton_signal_left_kinxi_kint"]->Fill( xi_cms_plus, event_weight );
              histosTH1F["xi_cms_totem_left_signal_kinxi_kint"]->Fill( xi_plus_Reco-xi_proton_plus_rec , event_weight );
          }
          if (proton_plus_kinec_accept_t_mc && proton_plus_kinec_accept_xi_mc) histosTH1F["xi_cms_totem_left_signal_kin"]->Fill( xi_plus_Reco-xi_proton_plus_rec , event_weight );
          ///xi_cms -xi_totem cut 
          if (proton_plus_kinec_accept_t_mc && xi_cms_totem_cut_plus){
              histosTH1F["xi_proton_signal_left_kint_cut"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["xi_proton_signal_left_kint_cut_bin2"]->Fill( xi_proton_plus_rec , event_weight );}
          if (proton_plus_kinec_accept_t_mc && proton_plus_kinec_accept_xi_mc  && xi_cms_totem_cut_plus){ 
              histosTH1F["xi_proton_signal_left_kin_cut"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["t_proton_signal_left_kin_cut"]->Fill( fabs(t_proton_plus_rec) , event_weight );
              histosTH1F["beta_proton_plus_signal_kin_cut"]->Fill( proton_plus_beta , event_weight );
              histosTH1F["log_x_parton_plus_signal_kin_cut"]->Fill( log10(x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_signal_kin_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_signal_kin_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_signal_kin_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["leadingJet_eta_left_signal_kin_cut"]->Fill( Jet1_eta, event_weight  );
              histosTH1F["secondJet_eta_left_signal_kin_cut"]->Fill( Jet2_eta, event_weight  );
              histosTH1F["DeltaPtJet_left_signal_kin_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_left_signal_kin_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_left_signal_kin_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_left_signal_kin_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_signal_left_kin_cut_bin2"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["beta_proton_plus_signal_kin_cut_bin2"]->Fill( proton_plus_beta , event_weight );
              histosTH1F["log_x_parton_plus_signal_kin_cut_bin2"]->Fill( log10(x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_signal_kin_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_signal_kin_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_signal_kin_cut_bin2"]->Fill( eta_average, event_weight  );
          }
          if (proton_plus_kinec_accept_t_mc && xi_proton_plus_rec<=0.1 && xi_cms_totem_cut_plus){
              histosTH1F["xi_plus_rec_proton_signal_left_kinxi_kint_cut"]->Fill( xi_cms_plus, event_weight );
              histosTH1F["xi_proton_signal_left_kint_kinxi_cut"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["t_proton_signal_left_kint_kinxi_cut"]->Fill( fabs(t_proton_plus_rec) , event_weight );
              histosTH1F["t_proton_signal_left_kint_kinxi_cut_bin"]->Fill( fabs(t_proton_plus_rec) , event_weight );
              histosTH1F["beta_proton_plus_signal_kint_kinxi_cut"]->Fill( proton_plus_beta , event_weight );
              histosTH1F["log_x_parton_plus_signal_kint_kinxi_cut"]->Fill( log10(x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_signal_kint_kinxi_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_signal_kint_kinxi_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_signal_kint_kinxi_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["leadingJet_eta_left_signal_kint_kinxi_cut"]->Fill( Jet1_eta, event_weight  );
              histosTH1F["secondJet_eta_left_signal_kint_kinxi_cut"]->Fill( Jet2_eta, event_weight  );
              histosTH1F["DeltaPtJet_left_signal_kint_kinxi_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_left_signal_kint_kinxi_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_left_signal_kint_kinxi_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_left_signal_kint_kinxi_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_signal_left_kint_kinxi_cut_bin2"]->Fill( xi_proton_plus_rec , event_weight );
              histosTH1F["beta_proton_plus_signal_kint_kinxi_cut_bin2"]->Fill( proton_plus_beta , event_weight );
              histosTH1F["log_x_parton_plus_signal_kint_kinxi_cut_bin2"]->Fill( log10(x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_signal_kint_kinxi_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_signal_kint_kinxi_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_signal_kint_kinxi_cut_bin2"]->Fill( eta_average, event_weight  );
         }
      }

      ///proton from ZB 
      if (select_Vertex && selectJet_rec && backg_right){
         ++vtx_jetrec_backg; 
         histosTH1F["t_proton_backg_right"]->Fill( -t_proton_right_zb , event_weight );
         histosTH1F["xi_proton_backg_right"]->Fill( -xi_proton_right_zb , event_weight );
         histosTH1F["beta_proton_minus_backg"]->Fill( proton_minus_beta_zb , event_weight );
         histosTH1F["xi_cms_totem_right_zb"]->Fill( xi_cms_minus+xi_proton_right_zb , event_weight );
         ///kinematic region
         if (proton_right_kinec_t_zb){
             histosTH1F["xi_proton_backg_right_kint"]->Fill( -xi_proton_right_zb , event_weight );
             histosTH1F["xi_proton_backg_right_kint_bin25"]->Fill( -xi_proton_right_zb , event_weight );}
         if (proton_right_kinec_t_zb && proton_right_kinec_xi_zb ) histosTH1F["xi_cms_totem_right_zb_kin"]->Fill( xi_cms_minus+xi_proton_right_zb , event_weight );
         if (proton_right_kinec_t_zb && -xi_proton_right_zb<0.1){
             histosTH1F["xi_proton_backg_right_kinxi_kint_sempeso"]->Fill( -xi_proton_right_zb , 1. );
             histosTH1F["xi_proton_backg_right_kinxi_kint_pesocuadrado"]->Fill( -xi_proton_right_zb , pow(event_weight,2) );
             histosTH1F["xi_proton_backg_right_kinxi_kint"]->Fill( -xi_proton_right_zb , event_weight );
             histosTH1F["xi_minus_rec_proton_backg_right_kinxi_kint"]->Fill( xi_cms_minus , event_weight );
             histosTH1F["xi_minus_cms_zb_kint"]->Fill( xi_cms_minus_zb , event_weight );
             histosTH1F["xi_cms_totem_right_zb_kinxi_kint"]->Fill( xi_cms_minus+xi_proton_right_zb , event_weight );
             histosTH1F["xi_cms_totem_right_zb_kinxi_kint_zb"]->Fill( xi_cms_minus_zb+xi_proton_right_zb , event_weight );
             xi_cms_st_back = xi_cms_minus;
             xi_totem_st_back = -xi_proton_right_zb;
             xi_cms_minus_totem_st_back = xi_cms_minus+xi_proton_right_zb;
             small_tree_back->Fill();
         }
         ///xi_cms-xi_totem cut
         if (proton_right_kinec_t_zb && xi_cms_totem_cut_right_zb ){
              histosTH1F["xi_proton_backg_right_kint_cut"]->Fill( -xi_proton_right_zb , event_weight );
              histosTH1F["xi_proton_backg_right_kint_cut_bin2"]->Fill( -xi_proton_right_zb , event_weight );}
         if (proton_right_kinec_t_zb && proton_right_kinec_xi_zb && xi_cms_totem_cut_right_zb){
              histosTH1F["xi_proton_backg_right_kin_cut"]->Fill( -xi_proton_right_zb , event_weight );
              histosTH1F["t_proton_backg_right_kin_cut"]->Fill( -t_proton_right_zb , event_weight );
              histosTH1F["beta_proton_minus_backg_kin_cut"]->Fill( proton_minus_beta_zb , event_weight );
              histosTH1F["log_x_parton_minus_backg_kin_cut"]->Fill(log10( x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_backg_kin_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_backg_kin_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["leadingJet_eta_right_backg_kin_cut"]->Fill( Jet1_eta, event_weight  );
              histosTH1F["secondJet_eta_right_backg_kin_cut"]->Fill( Jet2_eta, event_weight  );
              histosTH1F["EtaJet_average_right_backg_kin_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["DeltaPtJet_right_backg_kin_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_right_backg_kin_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_right_backg_kin_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_right_backg_kin_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_backg_right_kin_cut_bin2"]->Fill( -xi_proton_right_zb , event_weight );
              histosTH1F["beta_proton_minus_backg_kin_cut_bin2"]->Fill( proton_minus_beta_zb , event_weight );
              histosTH1F["log_x_parton_minus_backg_kin_cut_bin2"]->Fill(log10( x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_backg_kin_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_backg_kin_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_backg_kin_cut_bin2"]->Fill( eta_average, event_weight  );
         }
         if (proton_right_kinec_t_zb && -xi_proton_right_zb<0.1 && xi_cms_totem_cut_right_zb){
              ++vtx_jetrec_backg_proton_cut;
              histosTH1F["xi_minus_rec_proton_backg_right_kinxi_kint_cut"]->Fill( xi_cms_minus , event_weight );
              histosTH1F["xi_minus_cms_zb_kint_cut"]->Fill( xi_cms_minus_zb , event_weight );
              histosTH1F["xi_proton_backg_right_kint_kinxi_cut"]->Fill( -xi_proton_right_zb , event_weight );
              histosTH1F["t_proton_backg_right_kint_kinxi_cut"]->Fill( -t_proton_right_zb , event_weight );
              histosTH1F["t_proton_backg_right_kint_kinxi_cut_bin"]->Fill( -t_proton_right_zb , event_weight );
              histosTH1F["beta_proton_minus_backg_kint_kinxi_cut"]->Fill( proton_minus_beta_zb , event_weight );
              histosTH1F["log_x_parton_minus_backg_kint_kinxi_cut"]->Fill(log10( x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_backg_kint_kinxi_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_backg_kint_kinxi_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_backg_kint_kinxi_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["DeltaPtJet_right_backg_kint_kinxi_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_right_backg_kint_kinxi_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_right_backg_kint_kinxi_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_right_backg_kint_kinxi_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_backg_right_kint_kinxi_cut_bin2"]->Fill( -xi_proton_right_zb , event_weight );
              histosTH1F["beta_proton_minus_backg_kint_kinxi_cut_bin2"]->Fill( proton_minus_beta_zb , event_weight );
              histosTH1F["log_x_parton_minus_backg_kint_kinxi_cut_bin2"]->Fill(log10( x_minus) , event_weight );
              histosTH1F["leadingJet_pt_right_backg_kint_kinxi_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_right_backg_kint_kinxi_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_right_backg_kint_kinxi_cut_bin2"]->Fill( eta_average, event_weight  );
         }

      }

      if (select_Vertex && selectJet_rec && backg_left){
         histosTH1F["t_proton_backg_left"]->Fill( -t_proton_left_zb , event_weight );
         histosTH1F["xi_proton_backg_left"]->Fill( -xi_proton_left_zb , event_weight );
         histosTH1F["beta_proton_plus_backg"]->Fill( proton_plus_beta_zb , event_weight );
         histosTH1F["xi_cms_totem_left_zb"]->Fill( xi_plus_Reco+xi_proton_left_zb , event_weight );
         ///kinematic region
         if (proton_left_kinec_t_zb){ 
             histosTH1F["xi_proton_backg_left_kint"]->Fill( -xi_proton_left_zb , event_weight );
             histosTH1F["xi_proton_backg_left_kint_bin25"]->Fill( -xi_proton_left_zb , event_weight );}
         if (proton_left_kinec_t_zb && proton_left_kinec_xi_zb) histosTH1F["xi_cms_totem_left_zb_kin"]->Fill( xi_plus_Reco+xi_proton_left_zb , event_weight );
         if (proton_left_kinec_t_zb && -xi_proton_left_zb<0.1){ 
             histosTH1F["xi_proton_backg_left_kinxi_kint"]->Fill( -xi_proton_left_zb , event_weight );
             histosTH1F["xi_plus_rec_proton_backg_left_kinxi_kint"]->Fill( xi_cms_plus , event_weight );
             histosTH1F["xi_cms_totem_left_zb_kinxi_kint"]->Fill( xi_plus_Reco+xi_proton_left_zb , event_weight );}
         ///xi_cms-xi_totem cut
         if (proton_left_kinec_t_zb && xi_cms_totem_cut_left_zb ){
              histosTH1F["xi_proton_backg_left_kint_cut"]->Fill( -xi_proton_left_zb , event_weight );
              histosTH1F["xi_proton_backg_left_kint_cut_bin2"]->Fill( -xi_proton_left_zb , event_weight );}
         if (proton_left_kinec_t_zb && proton_left_kinec_xi_zb && xi_cms_totem_cut_left_zb){
              histosTH1F["xi_proton_backg_left_kin_cut"]->Fill( -xi_proton_left_zb , event_weight );
              histosTH1F["t_proton_backg_left_kin_cut"]->Fill( -t_proton_left_zb , event_weight );
              histosTH1F["beta_proton_plus_backg_kin_cut"]->Fill( proton_plus_beta_zb , event_weight );
              histosTH1F["log_x_parton_plus_backg_kin_cut"]->Fill( log10(x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_backg_kin_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_backg_kin_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_backg_kin_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["DeltaPtJet_left_backg_kin_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_left_backg_kin_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_left_backg_kin_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_left_backg_kin_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_backg_left_kin_cut_bin2"]->Fill( -xi_proton_left_zb , event_weight );
              histosTH1F["beta_proton_plus_backg_kin_cut_bin2"]->Fill( proton_plus_beta_zb , event_weight );
              histosTH1F["log_x_parton_plus_backg_kin_cut_bin2"]->Fill(log10( x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_backg_kin_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_backg_kin_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_backg_kin_cut_bin2"]->Fill( eta_average, event_weight  );
         }
         if (proton_left_kinec_t_zb && -xi_proton_left_zb<0.1 && xi_cms_totem_cut_left_zb){
              histosTH1F["xi_plus_rec_proton_backg_left_kinxi_kint_cut"]->Fill( xi_cms_plus , event_weight );
              histosTH1F["xi_proton_backg_left_kint_kinxi_cut"]->Fill( -xi_proton_left_zb , event_weight );
              histosTH1F["t_proton_backg_left_kint_kinxi_cut"]->Fill( -t_proton_left_zb , event_weight );
              histosTH1F["t_proton_backg_left_kint_kinxi_cut_bin"]->Fill( -t_proton_left_zb , event_weight );
              histosTH1F["beta_proton_plus_backg_kint_kinxi_cut"]->Fill( proton_plus_beta_zb , event_weight );
              histosTH1F["log_x_parton_plus_backg_kint_kinxi_cut"]->Fill(log10( x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_backg_kint_kinxi_cut"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_backg_kint_kinxi_cut"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_backg_kint_kinxi_cut"]->Fill( eta_average, event_weight  );
              histosTH1F["DeltaPtJet_left_backg_kint_kinxi_cut"]->Fill( deltapt, event_weight  );
              histosTH1F["DeltaEtaJet_left_backg_kint_kinxi_cut"]->Fill( deltaeta, event_weight  );
              histosTH1F["DeltaPhiJet_left_backg_kint_kinxi_cut"]->Fill( deltaphi, event_weight  );
              histosTH1F["Mass_Jet_left_backg_kint_kinxi_cut"]->Fill( mass_jets, event_weight  );

              histosTH1F["xi_proton_backg_left_kint_kinxi_cut_bin2"]->Fill( -xi_proton_left_zb , event_weight );
              histosTH1F["beta_proton_plus_backg_kint_kinxi_cut_bin2"]->Fill( proton_plus_beta_zb , event_weight );
              histosTH1F["log_x_parton_plus_backg_kint_kinxi_cut_bin2"]->Fill(log10( x_plus) , event_weight );
              histosTH1F["leadingJet_pt_left_backg_kint_kinxi_cut_bin2"]->Fill( Jet1_pt, event_weight  );
              histosTH1F["secondJet_pt_left_backg_kint_kinxi_cut_bin2"]->Fill( Jet2_pt, event_weight  );
              histosTH1F["EtaJet_average_left_backg_kint_kinxi_cut_bin2"]->Fill( eta_average, event_weight  );
         }
      }

      ///proton from MC and from ZB -> smallest xi
      if (proton_minus_rp_accept_mc && valid_proton_right_zb) {
         double t_proton_right_sel =  (xi_proton_minus_rec<-xi_proton_right_zb) ? fabs(t_proton_minus_rec) : -t_proton_right_zb;
         double xi_proton_right_sel =  (xi_proton_minus_rec<-xi_proton_right_zb) ? xi_proton_minus_rec : -xi_proton_right_zb;
         double beta_proton_right_sel = x_minus/xi_proton_right_sel;
         histosTH1F["t_proton_right_both"]->Fill( t_proton_right_sel , event_weight );
         histosTH1F["xi_proton_right_both"]->Fill( xi_proton_right_sel , event_weight );
         histosTH1F["xi_cms_totem_right_both"]->Fill( xi_cms_minus-xi_proton_right_sel , event_weight );
         histosTH1F["beta_proton_both_right"]->Fill( beta_proton_right_sel , event_weight );
         //kinematic region
          if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1) histosTH1F["xi_proton_both_right_kint"]->Fill( xi_proton_right_sel , event_weight );
          if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1 && xi_proton_right_sel>=0.03 && xi_proton_right_sel<=0.1) histosTH1F["xi_cms_totem_right_both_kin"]->Fill( xi_cms_minus-xi_proton_right_sel , event_weight );
          ///xi_cms -xi_totem cut 
          if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1 && xi_cms_minus-xi_proton_right_sel<=0) histosTH1F["xi_proton_both_right_kint_cut"]->Fill( xi_proton_right_sel , event_weight );
          if (t_proton_right_sel>=0.03 && t_proton_right_sel<=1 && xi_proton_right_sel>=0.03 && xi_proton_right_sel<=0.1 && xi_cms_minus-xi_proton_right_sel<=0){
              histosTH1F["xi_proton_both_right_kin_cut"]->Fill( xi_proton_right_sel , event_weight );
              histosTH1F["t_proton_both_right_kin_cut"]->Fill( t_proton_right_sel , event_weight );
              histosTH1F["beta_proton_both_right_kin_cut"]->Fill( beta_proton_right_sel , event_weight );
          }
      }

      
    }//end loop for events
   // cout <<"After the jet selection " << nevents_jets << " events  "<< endl;
   // cout <<"After GenPart selection " << nevents_gen << " events "<< endl;
   // cout <<"After PF selection " << nevents_pf << " events "<< endl;
   // cout <<"  "<< endl;

   file->Close();
 
  }//end of loop over files
  

  //output file
  TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  output->cd();
 
  event_selection->SetBinContent(1, nevents_total);
  event_selection->SetBinContent(2, vtx);
  event_selection->SetBinContent(3, nevents_pf); 
  event_selection->SetBinContent(4, jet_gen);
  event_selection->SetBinContent(5, vtx_jetrec); 
  event_selection->SetBinContent(6, vtx_jetrec_rp); 
  event_selection->SetBinContent(7, vtx_jetrec_rp_proton_cut); 
  event_selection->SetBinContent(8, vtx_jetrec_backg); 
  event_selection->SetBinContent(9, vtx_jetrec_backg_proton_cut);
  event_selection->GetXaxis()->SetBinLabel(1, "total");
  event_selection->GetXaxis()->SetBinLabel(2, "vtx");
  event_selection->GetXaxis()->SetBinLabel(3, "pf");
  event_selection->GetXaxis()->SetBinLabel(4, "jetgen");
  event_selection->GetXaxis()->SetBinLabel(5, "vtx_jetrec");
  event_selection->GetXaxis()->SetBinLabel(6, "vtx_jetrec_rp");
  event_selection->GetXaxis()->SetBinLabel(7, "vtx_jetrec_rp_kinproton_cut");
  event_selection->GetXaxis()->SetBinLabel(8, "vtx_jetrec_backg");
  event_selection->GetXaxis()->SetBinLabel(9, "vtx_jetrec_backg_proton_cut");
  event_selection->Write();

  // Jet Trigger
  float cross_section = 2.998e+10; //2.213e10; //pb cross section
  float luminity_HLT_L1Jet1_198902 = 0.015879;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet1_198903 = 0.008698;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198903
  float luminity_HLT_L1Jet2_198902 = 0.015879;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet2_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198903
  float luminity = luminity_HLT_L1Jet1_198902 + luminity_HLT_L1Jet1_198903 + luminity_HLT_L1Jet2_198902 + luminity_HLT_L1Jet2_198903;
 
  float n_events = luminity*cross_section;
  
  float f1 = (float) nevents_total;
  //float f1 = (float) nweight_total;
  Double_t scale = n_events/f1;
  cout << "Total events Pythia:   " << nevt_tot << "  total weights "<<nweight_total<< "   scale  "<< scale << endl;   
  
//   float f3 = (float) weight_total_PF_selected ; 
//   Double_t scale_PF = 1.0/f3;
  
  histosTH2F["xi_plus_reco_gen"]->SetOption("colz");
  histosTH2F["xi_minus_reco_gen"]->SetOption("colz");
  histosTH2F["logxi_plus_reco_gen"]->SetOption("colz");
  histosTH2F["logxi_minus_reco_gen"]->SetOption("colz");

  //rdm_evt_zb->Write();


  small_tree->Write();
  small_tree_back->Write();
  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();it_histo != histosTH1F.end(); ++it_histo){
     (*it_histo).second->Scale(scale);
     (*it_histo).second->Write();
  }
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

  output->Close();
}
