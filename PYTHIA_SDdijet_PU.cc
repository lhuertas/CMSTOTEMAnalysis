//STANDARD ROOT INCLUDES
#include <TROOT.h>
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

#include "rp_aperture_config.h"

#define PI 3.141592653589793
using namespace std;

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
};


void PYTHIA_SDdijet_PU(string const& outputFileName = "pythia_SDdijet_pu.root", const Int_t nevt_max = -1){
  
  bool verbose = false;
  string treeName = "evt";//"cms_totem";
  string jetCollName = "ak5PFJets";
  string jetCorrName = "ak5PFL2L3Residual";
  jetCorrName = "ak5PFL2L3"; 
  
  const Int_t nevt_max_corr = (nevt_max >= 0) ? nevt_max : 99999999;


  // Declaration of histograms
  map<string,TH1F*> histosTH1F;
//   histosTH1F["decisionPhysTrig"] = new TH1F("decisionPhysTrig", "decisionPhysTrig" , 128 , 0 , 128);
//   histosTH1F["decisionTechTrig"] = new TH1F("decisionTechTrig", "decisionTechTrig" , 64 , 0 , 64);

/*  int nBinsHLT = hltPathNames.size(); 
  histosTH1F["hltTrigFired"] = new TH1F("hltTrigFired", "hltTrigFired" , nBinsHLT , 0 , nBinsHLT);*/
//   for(size_t k = 0; k < nBinsHLT; ++k) 
//      histosTH1F["hltTrigFired"]->GetXaxis()->SetBinLabel( (k + 1) , hltPathNames[k].c_str() );

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
  histosTH1F["leadingJet_pt_selected"] = new TH1F("leadingJet_pt_selected", "p_{T}(jet)" , 20 , 0. , 230.);
  histosTH1F["leadingJet_eta_selected"] = new TH1F("leadingJet_eta_selected", "#eta(jet)" , 20 , -5.2 , 5.2);

  histosTH1F["secondJet_pt"] = new TH1F("secondJet_pt", "p_{T}(jet)" , 20 , 0. , 230.);
  histosTH1F["secondJet_eta"] = new TH1F("secondJet_eta", "#eta(jet)" , 20 , -5.2 , 5.2);
  histosTH1F["secondJet_phi"] = new TH1F("secondJet_phi", "#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["secondJet_pt_selected"] = new TH1F("secondJet_pt_selected", "p_{T}(jet)" , 20 , 0. , 230.);
  histosTH1F["secondJet_eta_selected"] = new TH1F("secondJet_eta_selected", "#eta(jet)" , 20 , -5.2 , 5.2);

  histosTH1F["DeltaPtJet"] = new TH1F("Delta_pt_Jet", "#Delta p_{T}(jet)" , 20 , 0. , 150.);
  histosTH1F["DeltaEtaJet"] = new TH1F("Delta_eta_Jet", "#Delta#eta(jet)" , 20 , 0 , 5.2);
  histosTH1F["DeltaPhiJet"] = new TH1F("Delta_phi_Jet", "#Delta#phi(jet)" , 20 , -M_PI , M_PI);
  histosTH1F["Mass_Jet"] = new TH1F("mass_Jet", "Mass(jet)" , 20 , 0 , 450);

  histosTH1F["Eta_max"] = new TH1F("Eta_max", "#eta^{max}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Eta_min"] = new TH1F("Eta_min", "#eta^{min}" , 82 , etaBinsHCALBoundaries);
  histosTH1F["Delta_eta_maxmin"] = new TH1F("Delta_eta_maxmin", "#eta^{max} - #eta^{min}" , 50 , 0,11);
  histosTH1F["xi_plus_Reco"] = new TH1F("xi+", "#xi^{+}" , 82 , 0,4);
  histosTH1F["xi_minus_Reco"] = new TH1F("xi-", "#xi^{-}" , 82 , 0,4);
  histosTH1F["logxi_plus"] = new TH1F("logxi+", "Log #xi^{+}" , 82 , -3,0.5);
  histosTH1F["logxi_plus_gen"] = new TH1F("logxi+_gen", "Log #xi_{+}^{gen}" , 82 , -3,0.5);
  histosTH1F["logxi_minus_gen"] = new TH1F("logxi-_gen", "Log #xi_{-}^{gen}" , 82 , -3,0.5);
  histosTH1F["correction"] = new TH1F("correction", "Correction factor" , 82 , 0,2);
  histosTH1F["resolution_after"] = new TH1F("resolution_after", "Resolution" , 82 , -2,2);
  histosTH1F["resolution_before"] = new TH1F("resolution_before", "Resolution" , 82 , -2,2);

  histosTH1F["pf_xiPlus_minus_proton_left_xi"] = new TH1F("pf_xiPlus_minus_proton_left_xi", "#xi diff." , 200 , -1. , 1.);
  histosTH1F["pf_xiMinus_minus_proton_right_xi"] = new TH1F("pf_xiMinus_minus_proton_right_xi", "#xi diff." , 200 , -1. , 1.);

  Float_t tbins[12] = { 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.22, 0.30, 0.40, 0.50, 0.65, 1.};
  histosTH1F["xi_proton_minus"] = new TH1F("xi_proton_minus", "xi_proton_minus" , 20, -0.1, 0.3);
  histosTH1F["xi_proton_plus"] = new TH1F("xi_proton_plus", "xi_proton_plus" , 20, -0.1, 0.3);
  histosTH1F["t_proton_plus"] = new TH1F("t_proton_plus", "t_proton_plus" , 11, tbins);
  histosTH1F["t_proton_minus"] = new TH1F("t_proton_minus", "t_proton_minus" , 11, tbins);
  histosTH1F["xi_proton_right_zb"] = new TH1F("xi_proton_right_zb", "xi_proton_right" , 20, -0.1, 0.3);
  histosTH1F["xi_proton_left_zb"] = new TH1F("xi_proton_left_zb", "xi_proton_left" , 20, -0.1, 0.3);
  histosTH1F["t_proton_right_zb"] = new TH1F("t_proton_right_zb", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_right_zb_sel"] = new TH1F("t_proton_right_zb_sel", "t_proton_right" , 11, tbins);
  histosTH1F["t_proton_left_zb"] = new TH1F("t_proton_left_zb", "t_proton_left" , 11, tbins);
  histosTH1F["t_proton_left_zb_sel"] = new TH1F("t_proton_left_zb_sel", "t_proton_left" , 11, tbins);


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

  for(map<string,TH1F*>::const_iterator it = histosTH1F.begin(); it != histosTH1F.end(); ++it)
      it->second->Sumw2();
  for(map<string,TH2F*>::const_iterator it = histosTH2F.begin(); it != histosTH2F.end(); ++it)
      it->second->Sumw2();
  gStyle->SetPalette(1);
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
  MyGenKin*  genKin   = NULL;
  //===============================

  //ZeroBias Files
  string treeNameZB = "cms_totem";
  TChain treeZB("cms_totem");
  vector<TString>* vdirs_zb = new vector<TString>;
  vdirs_zb->push_back("/storage1/lhuertas/Analysis/CMSTotem/data/MergedNtuples/HighBeta/ZeroBias/");
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
  vector<MyVertex>*  vertex_coll_ZB  = NULL;
//  TFile* fileZB = TFile::Open(fileNameZB.c_str(),"READ");
    
//  treeZB = (TTree*)fileZB->Get( treeNameZB.c_str() );
  int nevZB = int(treeZB.GetEntries());
  
  treeZB.SetBranchAddress("cmsEvtUA",&evtIdZB);
  treeZB.SetBranchAddress("rec_prot_left.",&rec_proton_left);
  treeZB.SetBranchAddress("rec_prot_right.",&rec_proton_right);
  treeZB.SetBranchAddress("cmsVerticesUA",&vertex_coll_ZB);

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

     zeroBias.push_back(mydata);
  }
  //fileZB->Close();

  cout << zeroBias.size() << " events analyzed" << endl;
  
 // for(vector<MyTOTEMData>::const_iterator it_zb = zeroBias.begin(); it_zb != zeroBias.end(); it_zb++){
 //    cout<<it_zb->proton_rec_right_t<<endl;  
 // }

  //================================================


  rp_aperture_config();

  
  double nevents_jets = 0; 
  double nevents_pf = 0; 
  double nevents_gen = 0; 
  double nevents_total = 0; 
  double events_jets = 0;
  double events_pf = 0;
  double events_gen = 0;

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
      
      //AT THIS POINT ON, CAN START USING THE VARIABLES LIKE TRACKS, VERTEX ETC !
      double event_weight = genKin->genWeight; 
     // double event_weight = 1.0;
      nweight_total += event_weight; 
      ++nevents_total;
      
      
      //-------------------------------------------------------------------------------------------------
 	    
     // Vertices
      for(vector<MyVertex>::iterator it_vtx = vertex_coll->begin() ; it_vtx != vertex_coll->end() ; ++it_vtx){
        if (it_vtx!=vertex_coll->begin()) continue;
	if( it_vtx->ndof>4 ) passedvtx = true;   
            //histosTH1F["vtx_zpos"]->Fill( it_vtx->z, event_weight );
            //histosTH1F["vtx_xpos"]->Fill( it_vtx->x, event_weight );
            //histosTH1F["vtx_ypos"]->Fill( it_vtx->y, event_weight );}
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
	 histosTH1F["leadingJet_pt"]->Fill( leadingJet.Pt(), event_weight  );
	 histosTH1F["leadingJet_eta"]->Fill( leadingJet.Eta(), event_weight  );
	 histosTH1F["leadingJet_phi"]->Fill( leadingJet.Phi(), event_weight  );
	 
 	 if(leadingJet.Pt() > 30. && leadingJet.Eta()<4. && leadingJet.Eta()>-4.){ 
            jet1_selected = true;
	    weight_total_leadingJet_selected += event_weight;
	    Jet1_pt = leadingJet.Pt(); 
	    Jet1_px = leadingJet.Px(); 
	    Jet1_py = leadingJet.Py(); 
	    Jet1_pz = leadingJet.Pz(); 
	    Jet1_energy = leadingJet.E(); 
	    Jet1_eta = leadingJet.Eta(); 
	    Jet1_phi = leadingJet.Phi(); }
      }
      if(!jet1_selected) continue;
      
      if( pfJet_coll->size() > 1 ){
	 MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet["ak5PFJets"];
         weight_total_secondJet += event_weight;
         Jet2_pt = secondJet.Pt(); 
	 Jet2_eta = secondJet.Eta(); 
	 Jet2_phi = secondJet.Phi(); 
         histosTH1F["secondJet_pt"]->Fill( secondJet.Pt(), event_weight  );  
	 histosTH1F["secondJet_eta"]->Fill( secondJet.Eta(), event_weight  );
	 histosTH1F["secondJet_phi"]->Fill( secondJet.Phi(), event_weight  );
	 
 	 if(secondJet.Pt() > 30. && secondJet.Eta()<4. && secondJet.Eta()>-4.){  
            jet2_selected = true;
	    weight_total_secondJet_selected += event_weight;
	    Jet2_pt = secondJet.Pt(); 
            Jet2_px = secondJet.Px();
            Jet2_py = secondJet.Py();
            Jet2_pz = secondJet.Pz();
            Jet2_energy = secondJet.E();
	    Jet2_eta = secondJet.Eta(); 
	    Jet2_phi = secondJet.Phi(); }
      }
      if(!jet2_selected) continue;
      weight_total_Jet_selected += event_weight;
//       histosTH1F["DeltaPhiJet"]->Fill( deltaphi, event_weight  );	
      ++nevents_jets;       
      double mass_jets= sqrt(pow(Jet1_energy+Jet2_energy,2)-pow(Jet1_px+Jet2_px,2)-pow(Jet1_py+Jet2_py,2)-pow(Jet1_pz+Jet2_pz,2));
      
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
	 
	 if (status != 1) continue; //final state for the particles
	 if (id != 2212) continue;
//	 if (eta_gen<4.9 && eta_gen>-4.9){ xi_posit_gen = true;
            energy_gen = it_genpart->Energy();
            px_gen = it_genpart->Px();
            py_gen = it_genpart->Py();
            pz_gen = it_genpart->Pz();
	    proton_pf = sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);  
            double pz_cut = 0.7*proton_pi;
            if (fabs(pz_gen) < pz_cut) continue;

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

      double xi_plus_gen = genEPlusPz/cm; //cout<<xi1_gen<<endl;
      double xi_minus_gen = genEMinusPz/cm;
      double xi_proton_plus = -1.;
      double xi_proton_minus = -1.;
      double t_proton_plus = 0.;
      double t_proton_minus = 0.;
      double thx_proton_plus = 0.;
      double thy_proton_plus = 0.;
      double thx_proton_minus = 0.;
      double thy_proton_minus = 0.;

      ++nevents_gen ;

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

      if(proton_pz_plus > 0.){
         xi_proton_plus =  ( 1 - (proton_pz_plus/proton_pi) );
         t_proton_plus = -2*( (proton_pi*proton_energy_plus) - (proton_pi*proton_pz_plus) );
         thx_proton_plus = atan(-proton_px_plus/proton_pi);
         thy_proton_plus = atan(proton_py_plus/proton_pi);

         //FIXME
         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_plus_rp_accept_020 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 20, out_x, out_thx, out_y, out_thy, out_xi);
         proton_plus_pars[20] = std::vector<double>(5,0.);
         proton_plus_pars[20][0] = out_x; proton_plus_pars[20][1] = out_y;
         proton_plus_pars[20][2] = out_thx; proton_plus_pars[20][3] = out_thy;
         proton_plus_pars[20][4] = out_xi;

         proton_plus_rp_accept_021 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 21);
         proton_plus_rp_accept_022 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 22);
         proton_plus_rp_accept_023 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 23);
         proton_plus_rp_accept_024 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 24);
         proton_plus_rp_accept_025 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 25);
         proton_plus_rp_accept_120 = protonRPDetected(0., thx_proton_plus, 0., thy_proton_plus, -xi_proton_plus, 120);
      }

      if(proton_pz_minus < 0.){
         xi_proton_minus = (proton_pz_minus < 0.) ? ( 1 + (proton_pz_minus/proton_pi) ) : -1.;
         t_proton_minus = -2*( (proton_pi*proton_energy_minus) + (proton_pi*proton_pz_minus) ); 

         thx_proton_minus = atan(-proton_px_minus/proton_pi);
         thy_proton_minus = atan(proton_py_minus/proton_pi);

         double out_x, out_thx, out_y, out_thy, out_xi;
         proton_minus_rp_accept_120 = protonRPDetected(0., thx_proton_minus, 0., thy_proton_minus, -xi_proton_minus, 120, out_x, out_thx, out_y, out_thy, out_xi);
         proton_minus_pars[120] = std::vector<double>(5,0.);
         proton_minus_pars[120][0] = out_x; proton_minus_pars[120][1] = out_y;
         proton_minus_pars[120][2] = out_thx; proton_minus_pars[120][3] = out_thy;
         proton_minus_pars[120][4] = out_xi;

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
/*      double cm = 8000;*/
      
      for(vector<MyPFCand>::iterator it_pfcand = pFlow_coll->begin(); it_pfcand != pFlow_coll->end(); ++it_pfcand){
	 int partType = it_pfcand->particleId;
	 double eta = it_pfcand->Eta();
	 double energy = it_pfcand->Energy();
	 double pz = it_pfcand->Pz();

	 // HF eta rings 29, 30, 40, 41
         if( ( (fabs(eta) >= 2.866) && (fabs(eta) < 3.152) ) || (fabs(eta) >= 4.730) ) continue;

	 soma1 += (energy + pz);
	 soma2 += (energy - pz);

         if (eta > eta_max) {eta_max = eta; PF_eta_max = true;} 
	 if (eta < eta_min) {eta_min = eta; PF_eta_min = true;}

  
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
       if(!PF_eta_max) continue;
       if(!PF_eta_min) continue;
       weight_total_PF_selected += event_weight;

       ++nevents_pf;

       double xi_plus_Reco = soma1/cm;
       double xi_minus_Reco = soma2/cm;
       double delta_eta_maxmin = eta_max - eta_min;  
      
       double correction = xi_plus_Reco/xi_plus_gen;
       double resolution_before = (xi_plus_gen-xi_plus_Reco)/xi_plus_gen;
       double xi_reconst = xi_plus_Reco/0.8;
       double resolution_after = (xi_plus_gen-xi_reconst)/xi_plus_gen;


      //rp_accept
      histosTH1F["xi_proton_minus"]->Fill( xi_proton_minus , event_weight );
      histosTH1F["t_proton_minus"]->Fill( fabs(t_proton_minus) , event_weight );

      bool proton_minus_rp_accept = ( proton_minus_rp_accept_120 && proton_minus_rp_accept_124 ) || ( proton_minus_rp_accept_121 && proton_minus_rp_accept_125 );

      //if(!proton_minus_rp_accept)continue;
      //   if (fabs(t_proton_minus)<0.03 || fabs(t_proton_minus)>1)continue;
      //   if (xi_proton_minus<0.03 || xi_proton_minus>0.1)continue;

      //Access Zero-bias data
      int i_evt_ZB = 0 + gRandom->Rndm()*(zeroBias.size());
      MyZeroBiasData const & zeroBiasData = zeroBias.at(i_evt_ZB);
      double xi_proton_right_zb = zeroBiasData.proton_rec_right_xi;      
      double t_proton_right_zb = zeroBiasData.proton_rec_right_t;
      bool valid_proton_right_zb = zeroBiasData.proton_rec_right_valid;      
      bool valid_proton_left_zb = zeroBiasData.proton_rec_left_valid;      
      double xi_proton_left_zb = zeroBiasData.proton_rec_left_xi;      
      double t_proton_left_zb = zeroBiasData.proton_rec_left_t;
      bool valid_vtx = zeroBiasData.vtx_valid;
      double ndof_vtx = zeroBiasData.vtx_ndof;
      if (valid_vtx && ndof_vtx>=4) continue;

  //    if (valid_proton_left_zb) continue;
      bool proton_right_zb = valid_proton_right_zb && (fabs(xi_proton_right_zb)>=0.03 && fabs(xi_proton_right_zb)<=0.1) && (fabs(t_proton_right_zb)>=0.03 && fabs(t_proton_right_zb)<=1);      
      bool proton_left_zb = valid_proton_left_zb && (fabs(xi_proton_left_zb)>=0.03 && fabs(xi_proton_left_zb)<=0.1) && (fabs(t_proton_left_zb)>=0.03 && fabs(t_proton_left_zb)<=1);      
      if (proton_right_zb) {
         histosTH1F["xi_proton_right_zb"]->Fill( xi_proton_right_zb , event_weight );
         histosTH1F["t_proton_right_zb"]->Fill( fabs(t_proton_right_zb) , event_weight );
         if (xi_minus_gen-xi_proton_right_zb<0) histosTH1F["t_proton_right_zb_sel"]->Fill( fabs(t_proton_right_zb) , event_weight );
      }

      if (proton_left_zb) {
         histosTH1F["xi_proton_left_zb"]->Fill( xi_proton_left_zb , event_weight );
         histosTH1F["t_proton_left_zb"]->Fill( fabs(t_proton_left_zb) , event_weight );
         if (xi_plus_gen-xi_proton_left_zb<0) histosTH1F["t_proton_left_zb_sel"]->Fill( fabs(t_proton_left_zb) , event_weight );
      }

      histosTH2F["xi_plus_reco_gen"]->Fill( xi_plus_gen, xi_plus_Reco, event_weight );
      histosTH2F["xi_minus_reco_gen"]->Fill( xi_minus_gen, xi_minus_Reco, event_weight );
      histosTH2F["logxi_plus_reco_gen"]->Fill( log10(xi_plus_gen), log10(xi_plus_Reco), event_weight );
      histosTH2F["logxi_minus_reco_gen"]->Fill( log10(xi_minus_gen), log10(xi_minus_Reco), event_weight );
      histosTH1F["Eta_max"]->Fill( eta_max, event_weight  );
      histosTH1F["Eta_min"]->Fill( eta_min, event_weight  );
      histosTH1F["Delta_eta_maxmin"]->Fill( delta_eta_maxmin, event_weight  );
      histosTH1F["xi_plus_Reco"]->Fill( xi_plus_Reco, event_weight  );
      histosTH1F["xi_minus_Reco"]->Fill( xi_minus_Reco, event_weight  );
      histosTH1F["logxi_plus"]->Fill( log10(xi_plus_Reco), event_weight  );

      MyBaseJet const& leadingJet = ( pfJet_coll->at(0) ).mapjet["ak5PFJets"];
      MyBaseJet const& secondJet = ( pfJet_coll->at(1) ).mapjet["ak5PFJets"];
      //...

      histosTH1F["secondJet_pt_selected"]->Fill( secondJet.Pt(), event_weight  );
      histosTH1F["secondJet_eta_selected"]->Fill( secondJet.Eta(), event_weight  );
      histosTH1F["leadingJet_pt_selected"]->Fill( leadingJet.Pt(), event_weight  );
      histosTH1F["leadingJet_eta_selected"]->Fill( leadingJet.Eta(), event_weight  );

      Double_t deltapt = abs(Jet1_pt - Jet2_pt);
      Double_t deltaeta = abs(Jet1_eta - Jet2_eta);
      Double_t deltaphi = abs(Jet1_phi - Jet2_phi);
      histosTH1F["DeltaPtJet"]->Fill( deltapt, event_weight  );
      histosTH1F["DeltaEtaJet"]->Fill( deltaeta, event_weight  );
      histosTH1F["DeltaPhiJet"]->Fill( 1-deltaphi/PI, event_weight  );

      histosTH1F["Mass_Jet"]->Fill( mass_jets, event_weight  );

      
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
 
  // Jet Trigger
  float cross_section = 2.213e10; //pb cross section for LP_Jets
  float luminity_HLT_L1Jet1_198902 = 0.015879;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet1_198903 = 0.008698;//pb ---- luminity for LP_Jets1_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198903
  float luminity_HLT_L1Jet2_198902 = 0.015879;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198902
  float luminity_HLT_L1Jet2_198903 = 0.008698;//pb ---- luminity for LP_Jets2_Run2012C-PromptReco-v1-HighBetaJuly2012-Run198903
  float luminity = luminity_HLT_L1Jet1_198902 + luminity_HLT_L1Jet1_198903 + luminity_HLT_L1Jet2_198902 + luminity_HLT_L1Jet2_198903;
 
  float n_events = luminity*cross_section;
  
  float f1 = (float) nevents_total;
//   float f1 = (float) nweight_total;
  Double_t scale = n_events/f1;
  
//   float f3 = (float) weight_total_PF_selected ; 
//   Double_t scale_PF = 1.0/f3;
  
  histosTH2F["xi_plus_reco_gen"]->SetOption("colz");
  histosTH2F["xi_minus_reco_gen"]->SetOption("colz");
  histosTH2F["logxi_plus_reco_gen"]->SetOption("colz");
  histosTH2F["logxi_minus_reco_gen"]->SetOption("colz");

  for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();it_histo != histosTH1F.end(); ++it_histo){
     (*it_histo).second->Scale(scale);
     (*it_histo).second->Write();
  }
  for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();
                                  it_histo != histosTH2F.end(); ++it_histo)
     (*it_histo).second->Write();

  output->Close();
}
