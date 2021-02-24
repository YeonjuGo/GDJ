//Author: Chris McGinn (2020.02.26)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
//#include "TLorentzVector.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"

//Local                                                                                   
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/keyHandler.h"
#include "include/ncollFunctions_5TeV.h"
#include "include/returnFileList.h"
#include "include/sampleHandler.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"
#include "include/etaPhiFunc.h"

const long MAXTREESIZE = 2000000000000; // set maximum tree size from 10 GB to 1862 GB, so that the code does not switch to a new file after 10 GB
int gdjNtuplePreProc_phoTaggedJetRaa(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  configParser config(inConfig_p);
  std::vector<std::string> necessaryParams = {"MCPREPROCDIRNAME",
    "OUTFILENAME",
    "CENTFILENAME",
    "ISPP",
    "ISMC"};  
  if(!checkEnvForParams(inConfig_p, necessaryParams)) return 1;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const bool isPP = inConfig_p->GetValue("ISPP", 0);
  const bool isMC = inConfig_p->GetValue("ISMC", 0);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const std::string inDirStr = inConfig_p->GetValue("MCPREPROCDIRNAME", "");
  const std::string inCentFileName = inConfig_p->GetValue("CENTFILENAME", "");
  if(!check.checkDir(inDirStr)){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }
  //if(!check.checkFileExt(inCentFileName, "txt")) return 1;

  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();

  //NColl weights
  std::vector<double> ncollWeights;
  std::vector<double> centCounts;
  for(Int_t cI = 0; cI < 100; ++cI){
    ncollWeights.push_back(findAvgNColl_Cent(cI, cI+1));
    centCounts.push_back(0.0);
  }
  for(Int_t cI = 0; cI < 100; ++cI){
    ncollWeights[99-cI] /= ncollWeights[0];
  }  

  std::vector<std::string> fileList = returnFileList(inDirStr, ".root");
  if(fileList.size() == 0){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' contains no root files. return 1" << std::endl;
    return 1;
  }
  TFile* inFile_p = new TFile(fileList[0].c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
  //TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");


  //  configParser inConfigGlobal(fileConfig_p);

  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  /////////////////////////////////////////
  // photon extra calibration
  const std::string inPhotonExtraCalibFileName = inConfig_p->GetValue("PHOEXTRACALIBFILE", "");
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  TFile* phoExtCaliFile_p = new TFile(inPhotonExtraCalibFileName.c_str(), "READ");
  if( phoExtCaliFile_p->IsZombie()){
    std::cout << "There is no photon extra scale factor file ( " << inPhotonExtraCalibFileName << " ). return 1" << std::endl;
    return 1;
  }
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const int centBins[] = {0,10,30,50,80};
  int nCentBins = sizeof(centBins)/sizeof(int) -1;
  if(isPP) nCentBins = 1;
  const int nCENTBINS = nCentBins; 
  TF1* f_phoExtCalib_barrel[nCENTBINS];
  TF1* f_phoExtCalib_endcap[nCENTBINS];
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  for(int icent = 0; icent < nCENTBINS; icent++){
    TString centLabel = Form("Cent%dto%d", centBins[icent], centBins[icent+1]);
    if(isPP) centLabel = "PP";
    f_phoExtCalib_barrel[icent] = (TF1*) phoExtCaliFile_p -> Get(Form("f_scale_vs_pt_%s_Eta0p00to1p37", centLabel.Data()));
    f_phoExtCalib_endcap[icent] = (TF1*) phoExtCaliFile_p -> Get(Form("f_scale_vs_pt_%s_Eta1p52to2p37", centLabel.Data()));
    if( f_phoExtCalib_barrel[icent]->IsZombie() || f_phoExtCalib_endcap[icent]->IsZombie() ){
      std::cout << "Can't retrieve the photon extra scale factor TF1. return 1" << std::endl;
      return 1;
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  /////////////////////////////////////////
  // JET extra calibration
  const std::string inJetExtraCalibFileName = inConfig_p->GetValue("JETEXTRACALIBFILE", "");
  TFile* jetExtCaliFile_p = new TFile(inJetExtraCalibFileName.c_str(), "READ");
  if( jetExtCaliFile_p->IsZombie()){
    std::cout << "There is no jet extra scale factor file ( " << inJetExtraCalibFileName << " ). return 1" << std::endl;
    return 1;
  }
  //const int centBins[] = {0,10,30,50,80};
  //int nCentBins = sizeof(centBins)/sizeof(int) -1;
  //if(isPP) nCentBins = 1;
  //const int nCENTBINS = nCentBins; 
  TF1* f_jetExtCalib[nCENTBINS];
  for(int icent = 0; icent < nCENTBINS; icent++){
    TString centLabel = Form("Cent%dto%d", centBins[icent], centBins[icent+1]);
    if(isPP) centLabel = "PP";
    f_jetExtCalib[icent] = (TF1*) jetExtCaliFile_p -> Get(Form("f_scale_vs_pt_%s_InclusiveJets", centLabel.Data()));
    if( f_jetExtCalib[icent]->IsZombie()){
      std::cout << "Can't retrieve the jet extra scale factor TF1. return 1" << std::endl;
      return 1;
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  /////////////////////////////////////////
  bool getTracks = false;
  if(checkEnvForParams(inConfig_p, {"GETTRACKS"})) getTracks = inConfig_p->GetValue("GETTRACKS", 0);

  bool getR2jets = false;
  if(checkEnvForParams(inConfig_p, {"GETR2JETS"})) getR2jets = inConfig_p->GetValue("GETR2JETS", 0);

  bool getR10jets = false;
  if(checkEnvForParams(inConfig_p, {"GETR10JETS"})) getR10jets = inConfig_p->GetValue("GETR10JETS", 0);

  bool getTruthParticle= false;
  if(checkEnvForParams(inConfig_p, {"GETTRUTHPARTICLE"})) getTruthParticle = inConfig_p->GetValue("GETTRUTHPARTICLE", 0);

  bool isTest = false;
  if(checkEnvForParams(inConfig_p, {"ISTEST"})) isTest = inConfig_p->GetValue("ISTEST", 0);

  std::string topOutDir = "output";
  if(checkEnvForParams(inConfig_p, {"OUTDIRNAME"})){
    topOutDir = inConfig_p->GetValue("OUTDIRNAME", "");  
    if(!check.checkDir(topOutDir)){
      std::cout << "Given output directory \'" << topOutDir << "\' does not exist. Please create it or remove param to use local default \'output\' directory. return 1" << std::endl;
      return 1;
    }
  }
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir(topOutDir);
  check.doCheckMakeDir(topOutDir + "/" + dateStr);

  std::string outFileName = inConfig_p->GetValue("OUTFILENAME", "");
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = topOutDir + "/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("gammaJetTree_p", "");

  std::vector<std::string> outBranchesToAdd = {"cent",
    "sampleTag",
    "sampleWeight",
    "ncollWeight",
    "fullWeight",
    "photon_pt_b4ExtraCalib",
    "akt4hi_em_xcalib_jet_pt_b4ExtraCalib",
    "truthPhotonPt",
    "truthPhotonEta",
    "truthPhotonPhi",
    "truthPhotonIso2",
    "truthPhotonIso3",
    "truthPhotonIso4"};

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  Int_t runNumber_, eventNumber_;
  UInt_t lumiBlock_;
  Bool_t passesToroid_;

  Bool_t is_pileup_;
  Bool_t is_oo_pileup_;

  Float_t pthat_;
  Int_t cent_;
  Int_t sampleTag_;
  Float_t sampleWeight_;
  Float_t ncollWeight_;
  Float_t fullWeight_;

  const Int_t nTreeParton_ = 2;
  Float_t treePartonPt_[nTreeParton_];
  Float_t treePartonEta_[nTreeParton_];
  Float_t treePartonPhi_[nTreeParton_];
  Int_t treePartonId_[nTreeParton_];

  Float_t actualInteractionsPerCrossing_, averageInteractionsPerCrossing_;
  Int_t nvert_;
  std::vector<float>* vert_x_p=nullptr;
  std::vector<float>* vert_y_p=nullptr;
  std::vector<float>* vert_z_p=nullptr;
  std::vector<int>* vert_type_p=nullptr;
  std::vector<int>* vert_ntrk_p=nullptr;

  Float_t fcalA_et_, fcalC_et_, fcalA_et_Cos2_, fcalC_et_Cos2_, fcalA_et_Sin2_, fcalC_et_Sin2_, fcalA_et_Cos3_, fcalC_et_Cos3_, fcalA_et_Sin3_, fcalC_et_Sin3_, fcalA_et_Cos4_, fcalC_et_Cos4_, fcalA_et_Sin4_, fcalC_et_Sin4_, evtPlane2Phi_, evtPlane3Phi_, evtPlane4Phi_;

  Int_t ntrk_;
  std::vector<float> *trk_pt_p=nullptr;
  std::vector<float> *trk_eta_p=nullptr;
  std::vector<float> *trk_phi_p=nullptr;
  std::vector<float> *trk_charge_p=nullptr;
  std::vector<bool> *trk_tight_primary_p=nullptr;
  std::vector<bool> *trk_minbias_p=nullptr;
  std::vector<float> *trk_d0_p=nullptr;
  std::vector<float> *trk_z0_p=nullptr;
  std::vector<float> *trk_vz_p=nullptr;
  std::vector<float> *trk_theta_p=nullptr;
  std::vector<int> *trk_nPixelHits_p=nullptr;
  std::vector<int> *trk_nSCTHits_p=nullptr;
  std::vector<int> *trk_nBlayerHits_p=nullptr;

  Int_t truth_n_;
  std::vector<float>* truth_charge_p=nullptr;
  std::vector<float>* truth_e_p=nullptr;
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_pdg_p=nullptr;
  std::vector<int>* truth_type_p=nullptr;
  std::vector<int>* truth_origin_p=nullptr;
  std::vector<int>* truth_status_p=nullptr;

  Int_t truthOut_n_;
  std::vector<float>* truthOut_charge_p=new std::vector<float>;
  std::vector<float>* truthOut_pt_p=new std::vector<float>;
  std::vector<float>* truthOut_eta_p=new std::vector<float>;
  std::vector<float>* truthOut_phi_p=new std::vector<float>;
  std::vector<float>* truthOut_e_p=new std::vector<float>;
  std::vector<float>* truthOut_pdg_p=new std::vector<float>;
  std::vector<int>* truthOut_type_p=new std::vector<int>;
  std::vector<int>* truthOut_origin_p=new std::vector<int>;
  std::vector<int>* truthOut_status_p=new std::vector<int>;

  Float_t truthPhotonPt_;
  Float_t truthPhotonEta_;
  Float_t truthPhotonPhi_;
  Float_t truthPhotonIso2_;
  Float_t truthPhotonIso3_;
  Float_t truthPhotonIso4_;

  Int_t akt2hi_jet_n_;
  std::vector<float>* akt2hi_em_xcalib_jet_m_p=nullptr;
  //std::vector<float>* akt2hi_em_xcalib_jet_clean_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_uncorrpt_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_uncorreta_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_uncorre_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt2hi_double_calib_jet_pt_p=nullptr;
  std::vector<int>* akt2hi_truthpos_p=nullptr;

  Int_t akt4hi_jet_n_;
  std::vector<float>* akt4hi_em_xcalib_jet_m_p=nullptr;
  //std::vector<float>* akt4hi_em_xcalib_jet_clean_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_b4ExtraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_uncorrpt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_uncorreta_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_uncorre_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt4hi_double_calib_jet_pt_p=nullptr;
  std::vector<int>* akt4hi_truthpos_p=nullptr;

  //JES, JER systematics
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_0_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_1_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_2_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_3_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_4_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_5_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_6_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_7_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_8_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_9_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_10_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_11_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_12_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_13_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_14_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_15_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_16_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_17_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_0_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_1_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_2_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_3_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_4_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_5_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_6_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_7_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_8_p=nullptr;

  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_0_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_1_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_2_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_3_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_4_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_5_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_6_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_7_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_8_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_9_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_10_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_11_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_12_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_13_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_14_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_15_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_16_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JES_17_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_0_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_1_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_2_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_3_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_4_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_5_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_6_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_7_extraCalib_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_sys_JER_8_extraCalib_p=nullptr;

  Int_t akt10hi_jet_n_;
  std::vector<float>* akt10hi_em_xcalib_jet_m_p=nullptr;
  //std::vector<float>* akt10hi_em_xcalib_jet_clean_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_uncorrpt_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_uncorreta_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt10hi_em_xcalib_jet_uncorre_p=nullptr;
  std::vector<float>* akt10hi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt10hi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt10hi_constit_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt10hi_constit_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt10hi_double_calib_jet_pt_p=nullptr;
  std::vector<int>* akt10hi_truthpos_p=nullptr;

  Int_t akt2to10hi_jet_n_;
  std::vector<float>* akt2to10hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt2to10hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt2to10hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt2to10hi_em_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt2to10hi_em_xcalib_jet_m_p=nullptr;
  std::vector<int>* akt2to10hi_truthpos_p=nullptr;

  Int_t photon_n_;
  std::vector<float>* photon_pt_extraCalib_p=nullptr;
  std::vector<float>* photon_pt_sys1_extraCalib_p=nullptr;
  std::vector<float>* photon_pt_sys2_extraCalib_p=nullptr;
  std::vector<float>* photon_pt_sys3_extraCalib_p=nullptr;
  std::vector<float>* photon_pt_sys4_extraCalib_p=nullptr;
  //std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_pt_sys1_p=nullptr;
  std::vector<float>* photon_pt_sys2_p=nullptr;
  std::vector<float>* photon_pt_sys3_p=nullptr;
  std::vector<float>* photon_pt_sys4_p=nullptr;
  std::vector<float>* photon_pt_b4ExtraCalib_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;
  std::vector<bool>* photon_loose_p=nullptr;
  std::vector<unsigned int>* photon_isem_p=nullptr;
  std::vector<int>* photon_convFlag_p=nullptr;
  std::vector<float>* photon_Rconv_p=nullptr;
  std::vector<float>* photon_etcone20_p=nullptr;
  std::vector<float>* photon_etcone30_p=nullptr;
  std::vector<float>* photon_etcone40_p=nullptr;
  std::vector<float>* photon_topoetcone20_p=nullptr;
  std::vector<float>* photon_topoetcone30_p=nullptr;
  std::vector<float>* photon_topoetcone40_p=nullptr;
  //std::vector<float>* photon_etcone20_ptcorrected_p=nullptr;
  //std::vector<float>* photon_etcone30_ptcorrected_p=nullptr;
  //std::vector<float>* photon_etcone40_ptcorrected_p=nullptr;
  std::vector<float>* photon_Rhad1_p=nullptr;
  std::vector<float>* photon_Rhad_p=nullptr;
  std::vector<float>* photon_e277_p=nullptr;
  std::vector<float>* photon_Reta_p=nullptr;
  std::vector<float>* photon_Rphi_p=nullptr;
  std::vector<float>* photon_weta1_p=nullptr;
  std::vector<float>* photon_weta2_p=nullptr;
  std::vector<float>* photon_wtots1_p=nullptr;
  std::vector<float>* photon_f1_p=nullptr;
  std::vector<float>* photon_f3_p=nullptr;
  std::vector<float>* photon_fracs1_p=nullptr;
  std::vector<float>* photon_DeltaE_p=nullptr;
  std::vector<float>* photon_Eratio_p=nullptr;

  Int_t akt2_truth_jet_n_;
  std::vector<float>* akt2_truth_jet_pt_p=nullptr;
  std::vector<float>* akt2_truth_jet_eta_p=nullptr;
  std::vector<float>* akt2_truth_jet_phi_p=nullptr;
  std::vector<float>* akt2_truth_jet_e_p=nullptr;
  std::vector<float>* akt2_truth_jet_m_p=nullptr;
  std::vector<int>* akt2_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt2_truth_jet_recopos_p=nullptr;
  Int_t akt4_truth_jet_n_;
  std::vector<float>* akt4_truth_jet_pt_p=nullptr;
  std::vector<float>* akt4_truth_jet_eta_p=nullptr;
  std::vector<float>* akt4_truth_jet_phi_p=nullptr;
  std::vector<float>* akt4_truth_jet_e_p=nullptr;
  std::vector<float>* akt4_truth_jet_m_p=nullptr;
  std::vector<int>* akt4_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt4_truth_jet_recopos_p=nullptr;
  Int_t akt10_truth_jet_n_;
  std::vector<float>* akt10_truth_jet_pt_p=nullptr;
  std::vector<float>* akt10_truth_jet_eta_p=nullptr;
  std::vector<float>* akt10_truth_jet_phi_p=nullptr;
  std::vector<float>* akt10_truth_jet_e_p=nullptr;
  std::vector<float>* akt10_truth_jet_m_p=nullptr;
  std::vector<int>* akt10_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt10_truth_jet_recopos_p=nullptr;
  Int_t akt2to10_truth_jet_n_;
  std::vector<float>* akt2to10_truth_jet_pt_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_eta_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_phi_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_e_p=nullptr;
  std::vector<float>* akt2to10_truth_jet_m_p=nullptr;
  std::vector<int>* akt2to10_truth_jet_partonid_p=nullptr;
  std::vector<int>* akt2to10_truth_jet_recopos_p=nullptr;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  outTree_p->SetMaxTreeSize(MAXTREESIZE);
  outTree_p->Branch("runNumber", &runNumber_, "runNumber/I");
  outTree_p->Branch("eventNumber", &eventNumber_, "eventNumber/I");
  outTree_p->Branch("lumiBlock", &lumiBlock_, "lumiBlock/i");
  outTree_p->Branch("passesToroid", &passesToroid_, "passesToroid/i");
  outTree_p->Branch("is_pileup", &is_pileup_, "is_pileup/i");
  outTree_p->Branch("is_oo_pileup", &is_oo_pileup_, "is_oo_pileup/i");

  if(isMC) outTree_p->Branch("pthat", &pthat_, "pthat/F");
  if(!isPP) outTree_p->Branch("cent", &cent_, "cent/I");
  outTree_p->Branch("sampleTag", &sampleTag_, "sampleTag/I");
  if(isMC) outTree_p->Branch("sampleWeight", &sampleWeight_, "sampleWeight/F");
  if(!isPP && isMC) outTree_p->Branch("ncollWeight", &ncollWeight_, "ncollWeight/F");
  if(isMC) outTree_p->Branch("fullWeight", &fullWeight_, "fullWeight/F");

  if(isMC){
    outTree_p->Branch("treePartonPt", treePartonPt_, ("treePartonPt[" + std::to_string(nTreeParton_) + "]/F").c_str());
    outTree_p->Branch("treePartonEta", treePartonEta_, ("treePartonEta[" + std::to_string(nTreeParton_) + "]/F").c_str());
    outTree_p->Branch("treePartonPhi", treePartonPhi_, ("treePartonPhi[" + std::to_string(nTreeParton_) + "]/F").c_str());
    outTree_p->Branch("treePartonId", treePartonId_, ("treePartonId[" + std::to_string(nTreeParton_) + "]/I").c_str());
  }

  outTree_p->Branch("actualInteractionsPerCrossing", &actualInteractionsPerCrossing_, "actualInteractionsPerCrossing/F");
  outTree_p->Branch("averageInteractionsPerCrossing", &averageInteractionsPerCrossing_, "averageInteractionsPerCrossing/F");
  outTree_p->Branch("nvert", &nvert_, "nvert/I");

  outTree_p->Branch("vert_x", &vert_x_p);
  outTree_p->Branch("vert_y", &vert_y_p);
  outTree_p->Branch("vert_z", &vert_z_p);
  outTree_p->Branch("vert_type", &vert_type_p);
  outTree_p->Branch("vert_ntrk", &vert_ntrk_p);

  outTree_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  outTree_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F"); 
  outTree_p->Branch("fcalA_et_Cos2", &fcalA_et_Cos2_, "fcalA_et_Cos2/F");
  outTree_p->Branch("fcalC_et_Cos2", &fcalC_et_Cos2_, "fcalC_et_Cos2/F");
  outTree_p->Branch("fcalA_et_Sin2", &fcalA_et_Sin2_, "fcalA_et_Sin2/F");
  outTree_p->Branch("fcalC_et_Sin2", &fcalC_et_Sin2_, "fcalC_et_Sin2/F");

  outTree_p->Branch("fcalA_et_Cos3", &fcalA_et_Cos3_, "fcalA_et_Cos3/F");
  outTree_p->Branch("fcalC_et_Cos3", &fcalC_et_Cos3_, "fcalC_et_Cos3/F");
  outTree_p->Branch("fcalA_et_Sin3", &fcalA_et_Sin3_, "fcalA_et_Sin3/F");
  outTree_p->Branch("fcalC_et_Sin3", &fcalC_et_Sin3_, "fcalC_et_Sin3/F");

  outTree_p->Branch("fcalA_et_Cos4", &fcalA_et_Cos4_, "fcalA_et_Cos4/F");
  outTree_p->Branch("fcalC_et_Cos4", &fcalC_et_Cos4_, "fcalC_et_Cos4/F");
  outTree_p->Branch("fcalA_et_Sin4", &fcalA_et_Sin4_, "fcalA_et_Sin4/F");
  outTree_p->Branch("fcalC_et_Sin4", &fcalC_et_Sin4_, "fcalC_et_Sin4/F");

  outTree_p->Branch("evtPlane2Phi", &evtPlane2Phi_, "evtPlane2Phi/F");
  outTree_p->Branch("evtPlane3Phi", &evtPlane3Phi_, "evtPlane3Phi/F");
  outTree_p->Branch("evtPlane4Phi", &evtPlane4Phi_, "evtPlane4Phi/F");

  if(getTracks){
    outTree_p->Branch("ntrk", &ntrk_, "ntrk/I");
    outTree_p->Branch("trk_pt", &trk_pt_p);
    outTree_p->Branch("trk_eta", &trk_eta_p);
    outTree_p->Branch("trk_phi", &trk_phi_p);

    outTree_p->Branch("trk_charge", &trk_charge_p);
    outTree_p->Branch("trk_tight_primary", &trk_tight_primary_p);
    outTree_p->Branch("trk_minbias", &trk_minbias_p);
    outTree_p->Branch("trk_d0", &trk_d0_p);
    outTree_p->Branch("trk_z0", &trk_z0_p);
    outTree_p->Branch("trk_vz", &trk_vz_p);
    outTree_p->Branch("trk_theta", &trk_theta_p);
    outTree_p->Branch("trk_nPixelHits", &trk_nPixelHits_p);
    outTree_p->Branch("trk_nSCTHits", &trk_nSCTHits_p);
    outTree_p->Branch("trk_nBlayerHits", &trk_nBlayerHits_p);
  }

  if(isMC){
    if(getTruthParticle){
      //outTree_p->Branch("truth_n", &truth_n_, "truth_n/I");
      //outTree_p->Branch("truth_charge", &truth_charge_p);
      //outTree_p->Branch("truth_e", &truth_e_p);
      //outTree_p->Branch("truth_pt", &truth_pt_p);
      //outTree_p->Branch("truth_eta", &truth_eta_p);
      //outTree_p->Branch("truth_phi", &truth_phi_p);
      //outTree_p->Branch("truth_pdg", &truth_pdg_p);
      //outTree_p->Branch("truth_type", &truth_type_p);
      //outTree_p->Branch("truth_origin", &truth_origin_p);
      //outTree_p->Branch("truth_status", &truth_status_p);

      outTree_p->Branch("truth_n", &truthOut_n_, "truth_n/I");
      outTree_p->Branch("truth_charge", &truthOut_charge_p);
      outTree_p->Branch("truth_pt", &truthOut_pt_p);
      outTree_p->Branch("truth_eta", &truthOut_eta_p);
      outTree_p->Branch("truth_phi", &truthOut_phi_p);
      outTree_p->Branch("truth_pdg", &truthOut_pdg_p);
      outTree_p->Branch("truth_e", &truthOut_e_p);
      outTree_p->Branch("truth_type", &truthOut_type_p);
      outTree_p->Branch("truth_origin", &truthOut_origin_p);
      outTree_p->Branch("truth_status", &truthOut_status_p);
    }

    outTree_p->Branch("truthPhotonPt", &truthPhotonPt_, "truthPhotonPt/F");
    outTree_p->Branch("truthPhotonPhi", &truthPhotonPhi_, "truthPhotonPhi/F");
    outTree_p->Branch("truthPhotonEta", &truthPhotonEta_, "truthPhotonEta/F");
    outTree_p->Branch("truthPhotonIso2", &truthPhotonIso2_, "truthPhotonIso2/F");
    outTree_p->Branch("truthPhotonIso3", &truthPhotonIso3_, "truthPhotonIso3/F");
    outTree_p->Branch("truthPhotonIso4", &truthPhotonIso4_, "truthPhotonIso4/F");
  }

  if(getR2jets){
    outTree_p->Branch("akt2hi_jet_n", &akt2hi_jet_n_, "akt2hi_jet_n/I");
    outTree_p->Branch("akt2hi_em_xcalib_jet_m", &akt2hi_em_xcalib_jet_m_p);
    //outTree_p->Branch("akt2hi_em_xcalib_jet_clean", &akt2hi_em_xcalib_jet_clean_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_pt", &akt2hi_em_xcalib_jet_pt_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_uncorrpt", &akt2hi_em_xcalib_jet_uncorrpt_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_eta", &akt2hi_em_xcalib_jet_eta_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_uncorreta", &akt2hi_em_xcalib_jet_uncorreta_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_phi", &akt2hi_em_xcalib_jet_phi_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_e", &akt2hi_em_xcalib_jet_e_p);
    outTree_p->Branch("akt2hi_em_xcalib_jet_uncorre", &akt2hi_em_xcalib_jet_uncorre_p);
    outTree_p->Branch("akt2hi_constit_xcalib_jet_pt", &akt2hi_constit_xcalib_jet_pt_p);
    outTree_p->Branch("akt2hi_constit_xcalib_jet_eta", &akt2hi_constit_xcalib_jet_eta_p);
    outTree_p->Branch("akt2hi_constit_xcalib_jet_phi", &akt2hi_constit_xcalib_jet_phi_p);
    outTree_p->Branch("akt2hi_constit_xcalib_jet_e", &akt2hi_constit_xcalib_jet_e_p);
    outTree_p->Branch("akt2hi_double_calib_jet_pt", &akt2hi_double_calib_jet_pt_p);
    if(isMC) outTree_p->Branch("akt2hi_truthpos", &akt2hi_truthpos_p);
  }

  outTree_p->Branch("akt4hi_jet_n", &akt4hi_jet_n_, "akt4hi_jet_n/I");
  outTree_p->Branch("akt4hi_em_xcalib_jet_m", &akt4hi_em_xcalib_jet_m_p);
  //outTree_p->Branch("akt4hi_em_xcalib_jet_clean", &akt4hi_em_xcalib_jet_clean_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_extraCalib_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_pt_b4ExtraCalib", &akt4hi_em_xcalib_jet_pt_b4ExtraCalib_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_uncorrpt", &akt4hi_em_xcalib_jet_uncorrpt_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_uncorreta", &akt4hi_em_xcalib_jet_uncorreta_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_e", &akt4hi_em_xcalib_jet_e_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_uncorre", &akt4hi_em_xcalib_jet_uncorre_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_pt", &akt4hi_constit_xcalib_jet_pt_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_eta", &akt4hi_constit_xcalib_jet_eta_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_phi", &akt4hi_constit_xcalib_jet_phi_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_e", &akt4hi_constit_xcalib_jet_e_p);
  outTree_p->Branch("akt4hi_double_calib_jet_pt", &akt4hi_double_calib_jet_pt_p);
  if(isMC) outTree_p->Branch("akt4hi_truthpos", &akt4hi_truthpos_p);

  if(isMC){
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_0", &akt4hi_em_xcalib_jet_pt_sys_JES_0_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_1", &akt4hi_em_xcalib_jet_pt_sys_JES_1_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_2", &akt4hi_em_xcalib_jet_pt_sys_JES_2_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_3", &akt4hi_em_xcalib_jet_pt_sys_JES_3_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_4", &akt4hi_em_xcalib_jet_pt_sys_JES_4_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_5", &akt4hi_em_xcalib_jet_pt_sys_JES_5_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_6", &akt4hi_em_xcalib_jet_pt_sys_JES_6_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_7", &akt4hi_em_xcalib_jet_pt_sys_JES_7_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_8", &akt4hi_em_xcalib_jet_pt_sys_JES_8_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_9", &akt4hi_em_xcalib_jet_pt_sys_JES_9_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_10", &akt4hi_em_xcalib_jet_pt_sys_JES_10_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_11", &akt4hi_em_xcalib_jet_pt_sys_JES_11_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_12", &akt4hi_em_xcalib_jet_pt_sys_JES_12_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_13", &akt4hi_em_xcalib_jet_pt_sys_JES_13_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_14", &akt4hi_em_xcalib_jet_pt_sys_JES_14_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_15", &akt4hi_em_xcalib_jet_pt_sys_JES_15_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_16", &akt4hi_em_xcalib_jet_pt_sys_JES_16_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JES_17", &akt4hi_em_xcalib_jet_pt_sys_JES_17_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_0", &akt4hi_em_xcalib_jet_pt_sys_JER_0_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_1", &akt4hi_em_xcalib_jet_pt_sys_JER_1_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_2", &akt4hi_em_xcalib_jet_pt_sys_JER_2_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_3", &akt4hi_em_xcalib_jet_pt_sys_JER_3_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_4", &akt4hi_em_xcalib_jet_pt_sys_JER_4_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_5", &akt4hi_em_xcalib_jet_pt_sys_JER_5_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_6", &akt4hi_em_xcalib_jet_pt_sys_JER_6_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_7", &akt4hi_em_xcalib_jet_pt_sys_JER_7_extraCalib_p);
    outTree_p->Branch("akt4hi_em_xcalib_jet_pt_sys_JER_8", &akt4hi_em_xcalib_jet_pt_sys_JER_8_extraCalib_p);
  }

  if(getR10jets){  
    if(!isPP){  
      outTree_p->Branch("akt10hi_jet_n", &akt10hi_jet_n_, "akt10hi_jet_n/I");
      outTree_p->Branch("akt10hi_em_xcalib_jet_m", &akt10hi_em_xcalib_jet_m_p);
      //outTree_p->Branch("akt10hi_em_xcalib_jet_clean", &akt10hi_em_xcalib_jet_clean_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_pt", &akt10hi_em_xcalib_jet_pt_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_uncorrpt", &akt10hi_em_xcalib_jet_uncorrpt_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_eta", &akt10hi_em_xcalib_jet_eta_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_uncorreta", &akt10hi_em_xcalib_jet_uncorreta_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_phi", &akt10hi_em_xcalib_jet_phi_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_e", &akt10hi_em_xcalib_jet_e_p);
      outTree_p->Branch("akt10hi_em_xcalib_jet_uncorre", &akt10hi_em_xcalib_jet_uncorre_p);
      outTree_p->Branch("akt10hi_constit_xcalib_jet_pt", &akt10hi_constit_xcalib_jet_pt_p);
      outTree_p->Branch("akt10hi_constit_xcalib_jet_eta", &akt10hi_constit_xcalib_jet_eta_p);
      outTree_p->Branch("akt10hi_constit_xcalib_jet_phi", &akt10hi_constit_xcalib_jet_phi_p);
      outTree_p->Branch("akt10hi_constit_xcalib_jet_e", &akt10hi_constit_xcalib_jet_e_p);
      outTree_p->Branch("akt10hi_double_calib_jet_pt", &akt10hi_double_calib_jet_pt_p);
      if(isMC) outTree_p->Branch("akt10hi_truthpos", &akt10hi_truthpos_p);
    }
    outTree_p->Branch("akt2to10hi_jet_n", &akt2to10hi_jet_n_, "akt2to10hi_jet_n/I");
    outTree_p->Branch("akt2to10hi_em_xcalib_jet_pt", &akt2to10hi_em_xcalib_jet_pt_p);
    outTree_p->Branch("akt2to10hi_em_xcalib_jet_eta", &akt2to10hi_em_xcalib_jet_eta_p);
    outTree_p->Branch("akt2to10hi_em_xcalib_jet_phi", &akt2to10hi_em_xcalib_jet_phi_p);
    outTree_p->Branch("akt2to10hi_em_xcalib_jet_e", &akt2to10hi_em_xcalib_jet_e_p);
    outTree_p->Branch("akt2to10hi_em_xcalib_jet_m", &akt2to10hi_em_xcalib_jet_m_p);
    if(isMC) outTree_p->Branch("akt2to10hi_truthpos", &akt2to10hi_truthpos_p);
  }

  outTree_p->Branch("photon_n", &photon_n_, "photon_n/I");
  //outTree_p->Branch("photon_pt", &photon_pt_p);
  //outTree_p->Branch("photon_pt_sys1", &photon_pt_sys1_p);
  //outTree_p->Branch("photon_pt_sys2", &photon_pt_sys2_p);
  //outTree_p->Branch("photon_pt_sys3", &photon_pt_sys3_p);
  //outTree_p->Branch("photon_pt_sys4", &photon_pt_sys4_p);
  outTree_p->Branch("photon_pt", &photon_pt_extraCalib_p);
  if(isMC){
    outTree_p->Branch("photon_pt_sys1", &photon_pt_sys1_extraCalib_p);
    outTree_p->Branch("photon_pt_sys2", &photon_pt_sys2_extraCalib_p);
    outTree_p->Branch("photon_pt_sys3", &photon_pt_sys3_extraCalib_p);
    outTree_p->Branch("photon_pt_sys4", &photon_pt_sys4_extraCalib_p);
  }
  outTree_p->Branch("photon_pt_b4ExtraCalib", &photon_pt_b4ExtraCalib_p);
  outTree_p->Branch("photon_eta", &photon_eta_p);
  outTree_p->Branch("photon_phi", &photon_phi_p);
  outTree_p->Branch("photon_tight", &photon_tight_p);
  outTree_p->Branch("photon_loose", &photon_loose_p);
  outTree_p->Branch("photon_isem", &photon_isem_p);
  outTree_p->Branch("photon_convFlag", &photon_convFlag_p);
  outTree_p->Branch("photon_Rconv", &photon_Rconv_p);
  outTree_p->Branch("photon_etcone20", &photon_etcone20_p);
  outTree_p->Branch("photon_etcone30", &photon_etcone30_p);
  outTree_p->Branch("photon_etcone40", &photon_etcone40_p);
  outTree_p->Branch("photon_topoetcone20", &photon_topoetcone20_p);
  outTree_p->Branch("photon_topoetcone30", &photon_topoetcone30_p);
  outTree_p->Branch("photon_topoetcone40", &photon_topoetcone40_p);
  //outTree_p->Branch("photon_etcone20_ptcorrected", &photon_etcone20_ptcorrected_p);
  //outTree_p->Branch("photon_etcone30_ptcorrected", &photon_etcone30_ptcorrected_p);
  //outTree_p->Branch("photon_etcone40_ptcorrected", &photon_etcone40_ptcorrected_p);
  outTree_p->Branch("photon_Rhad1", &photon_Rhad1_p);
  outTree_p->Branch("photon_Rhad", &photon_Rhad_p);
  outTree_p->Branch("photon_e277", &photon_e277_p);
  outTree_p->Branch("photon_Reta", &photon_Reta_p);
  outTree_p->Branch("photon_Rphi", &photon_Rphi_p);
  outTree_p->Branch("photon_weta1", &photon_weta1_p);
  outTree_p->Branch("photon_weta2", &photon_weta2_p);
  outTree_p->Branch("photon_weta2", &photon_weta2_p);
  outTree_p->Branch("photon_wtots1", &photon_wtots1_p);
  outTree_p->Branch("photon_f1", &photon_f1_p);
  outTree_p->Branch("photon_f3", &photon_f3_p);
  outTree_p->Branch("photon_fracs1", &photon_fracs1_p);
  outTree_p->Branch("photon_DeltaE", &photon_DeltaE_p);
  outTree_p->Branch("photon_Eratio", &photon_Eratio_p);

  if(isMC){
    if(getR2jets){
      outTree_p->Branch("akt2_truth_jet_n", &akt2_truth_jet_n_, "akt2_truth_jet_n/I");
      outTree_p->Branch("akt2_truth_jet_pt", &akt2_truth_jet_pt_p);
      outTree_p->Branch("akt2_truth_jet_eta", &akt2_truth_jet_eta_p);
      outTree_p->Branch("akt2_truth_jet_phi", &akt2_truth_jet_phi_p);
      outTree_p->Branch("akt2_truth_jet_e", &akt2_truth_jet_e_p);
      outTree_p->Branch("akt2_truth_jet_m", &akt2_truth_jet_m_p);
      outTree_p->Branch("akt2_truth_jet_partonid", &akt2_truth_jet_partonid_p);
      outTree_p->Branch("akt2_truth_jet_recopos", &akt2_truth_jet_recopos_p);
    }

    outTree_p->Branch("akt4_truth_jet_n", &akt4_truth_jet_n_, "akt4_truth_jet_n/I");
    outTree_p->Branch("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
    outTree_p->Branch("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
    outTree_p->Branch("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
    outTree_p->Branch("akt4_truth_jet_e", &akt4_truth_jet_e_p);
    outTree_p->Branch("akt4_truth_jet_m", &akt4_truth_jet_m_p);
    outTree_p->Branch("akt4_truth_jet_partonid", &akt4_truth_jet_partonid_p);
    outTree_p->Branch("akt4_truth_jet_recopos", &akt4_truth_jet_recopos_p);

    if(getR10jets){
      if(!isPP){
        outTree_p->Branch("akt10_truth_jet_n", &akt10_truth_jet_n_, "akt10_truth_jet_n/I");
        outTree_p->Branch("akt10_truth_jet_pt", &akt10_truth_jet_pt_p);
        outTree_p->Branch("akt10_truth_jet_eta", &akt10_truth_jet_eta_p);
        outTree_p->Branch("akt10_truth_jet_phi", &akt10_truth_jet_phi_p);
        outTree_p->Branch("akt10_truth_jet_e", &akt10_truth_jet_e_p);
        outTree_p->Branch("akt10_truth_jet_m", &akt10_truth_jet_m_p);
        outTree_p->Branch("akt10_truth_jet_partonid", &akt10_truth_jet_partonid_p);
        outTree_p->Branch("akt10_truth_jet_recopos", &akt10_truth_jet_recopos_p);
      }

      outTree_p->Branch("akt2to10_truth_jet_n", &akt2to10_truth_jet_n_, "akt2to10_truth_jet_n/I");
      outTree_p->Branch("akt2to10_truth_jet_pt", &akt2to10_truth_jet_pt_p);
      outTree_p->Branch("akt2to10_truth_jet_eta", &akt2to10_truth_jet_eta_p);
      outTree_p->Branch("akt2to10_truth_jet_phi", &akt2to10_truth_jet_phi_p);
      outTree_p->Branch("akt2to10_truth_jet_e", &akt2to10_truth_jet_e_p);
      outTree_p->Branch("akt2to10_truth_jet_m", &akt2to10_truth_jet_m_p);
      outTree_p->Branch("akt2to10_truth_jet_partonid", &akt2to10_truth_jet_partonid_p);
      outTree_p->Branch("akt2to10_truth_jet_recopos", &akt2to10_truth_jet_recopos_p);
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::map<int, unsigned long long> tagToCounts;
  std::map<int, double> tagToXSec;
  std::map<int, double> tagToFilterEff;
  std::vector<int> minPthats;
  std::map<int, int> minPthatToTag;
  std::map<int, double> tagToWeight;

  std::map<std::string, std::vector<std::string> > configMap;
  configMap["MCPREPROCDIRNAME"] = {inDirStr};


  //Basic pre-processing for output config
  std::vector<std::string> listOfBranchesIn, listOfBranchesHLT, listOfBranchesHLTPre;
  std::vector<std::string> listOfBranchesOut = getVectBranchList(outTree_p);
  ULong64_t totalNEntries = 0;
  for(auto const & file : fileList){
    inFile_p = new TFile(file.c_str(), "READ");
    inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    totalNEntries += inTree_p->GetEntries();
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");

    if(!isPP){
      inTree_p->SetBranchStatus("*", 0);
      inTree_p->SetBranchStatus("fcalA_et", 1);
      inTree_p->SetBranchStatus("fcalC_et", 1);

      inTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
      inTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);

      for(Long64_t entry = 0; entry < inTree_p->GetEntries(); ++entry){
        inTree_p->GetEntry(entry);

        cent_ = centTable.GetCent(fcalA_et_ + fcalC_et_);
        ++(centCounts[cent_]);
      }
    }

    std::map<std::string, std::string> tempConfigMap = GetMapFromEnv(inConfig_p);

    if(isMC){
      sampleHandler sHandler;

      std::string inDataSetName = inConfig_p->GetValue("INDATASET", "");
      if(inDataSetName.size() == 0 || !sHandler.Init(inDataSetName)){
        std::cout << "GDJMCNTUPLEPREPROC ERROR - Given input \'" << file << "\' contains INDATASET \'" << inDataSetName << "\' that is not valid. return 1" << std::endl;

        inFile_p->Close();
        delete inFile_p;

        outFile_p->Close();
        delete outFile_p;

        return 1;
      }

      sampleTag_ = sHandler.GetTag();

      if(tagToCounts.count(sampleTag_) == 0){
        tagToCounts[sampleTag_] = inTree_p->GetEntries();
        tagToXSec[sampleTag_] = sHandler.GetXSection();
        tagToFilterEff[sampleTag_] = sHandler.GetFilterEff();
        int minPthat = sHandler.GetMinPthat();
        minPthats.push_back(minPthat);
        minPthatToTag[sHandler.GetMinPthat()] = sampleTag_;
      }
      else tagToCounts[sampleTag_] += inTree_p->GetEntries();
    }

    for(auto const & val : tempConfigMap){
      bool isFound = false;
      std::vector<std::string> tempVect = configMap[val.first];
      for(unsigned int vI = 0; vI < tempVect.size(); ++vI){
        if(isStrSame(tempVect[vI], val.second)){
          isFound = true;
          break;
        }
      }

      if(!isFound){
        tempVect.push_back(val.second);
        configMap[val.first] = tempVect;
      }      
    }

    std::vector<std::string> tempBranches = getVectBranchList(inTree_p);
    for(auto const & branch : tempBranches){
      if(branch.size() >= 4){
        if(branch.substr(0,4).find("HLT_") != std::string::npos){
          if(branch.find("presc") != std::string::npos){
            if(!vectContainsStr(branch, &listOfBranchesHLTPre)) listOfBranchesHLTPre.push_back(branch);
          }
          else{
            if(!vectContainsStr(branch, &listOfBranchesHLT)) listOfBranchesHLT.push_back(branch);
          }	    

          continue;
        }
      }

      if(!vectContainsStr(branch, &listOfBranchesIn)){
        listOfBranchesIn.push_back(branch);
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  if(isMC){
    std::sort(std::begin(minPthats), std::end(minPthats));
    std::cout << "Min pthats: " << std::endl;
    for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
      int tag = minPthatToTag[minPthats[pI]];
      double weight = tagToXSec[tag]*tagToFilterEff[tag]/(double)tagToCounts[tag];
      tagToWeight[tag] = weight;
    }

    for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
      if(pI == 0) continue;

      int tag = minPthatToTag[minPthats[pI]];
      int tag0 = minPthatToTag[minPthats[0]];

      tagToWeight[tag] /= tagToWeight[tag0];
    }

    int tag0 = minPthatToTag[minPthats[0]];
    tagToWeight[tag0] = 1.0;

    for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
      int tag = minPthatToTag[minPthats[pI]];

      std::cout << " " << pI << "/" << minPthats.size() << ": " << minPthats[pI] << ", " << tagToWeight[tag] << std::endl;
    }
  }

  std::vector<Bool_t*> hltVect;
  std::vector<Float_t*> hltPreVect;
  hltVect.reserve(listOfBranchesHLT.size());
  hltPreVect.reserve(listOfBranchesHLTPre.size());

  for(auto const & branchHLT : listOfBranchesHLT){
    hltVect.push_back(new Bool_t(false));

    outTree_p->Branch(branchHLT.c_str(), hltVect[hltVect.size()-1], (branchHLT + "/O").c_str());
  }
  for(auto const & branchHLTPre : listOfBranchesHLTPre){
    hltPreVect.push_back(new Float_t(0.0));

    outTree_p->Branch(branchHLTPre.c_str(), hltPreVect[hltPreVect.size()-1], (branchHLTPre + "/F").c_str());
  }

  bool allBranchesGood = true;
  //for(auto const & branchIn : listOfBranchesIn){
  //  bool containsBranch = vectContainsStr(branchIn, &listOfBranchesOut);
  //  allBranchesGood = allBranchesGood && containsBranch;

  //  if(!containsBranch) std::cout << "GDJMCNTUPLEPREPROC ERROR - \'" << branchIn << "\' branch is not part of output tree. please fix macro, return 1" << std::endl;
  //}

  bool allBranchesGood2 = true;
  for(auto const & branchOut : listOfBranchesOut){
    if(vectContainsStr(branchOut, &outBranchesToAdd)) continue;//Dont want to throw error on branches we added

    bool containsBranch = vectContainsStr(branchOut, &listOfBranchesIn);
    allBranchesGood2 = allBranchesGood2 && containsBranch;

    if(!containsBranch) std::cout << "GDJMCNTUPLEPREPROC ERROR - \'" << branchOut << "\' branch is not part of input tree. please fix macro, return 1" << std::endl;
  }

  if(!allBranchesGood || !allBranchesGood2){
    outFile_p->Close();
    delete outFile_p;
    return 1;
  }

  if(!isPP && isMC){
    for(unsigned int cI = 0; cI < ncollWeights.size(); ++cI){      
      if(centCounts[cI] < 0.5) continue;
      ncollWeights[cI] /= centCounts[cI];
    }
  }

  sampleHandler sHandler;

  ULong64_t nDiv = TMath::Max((ULong64_t)1, totalNEntries/20);
  ULong64_t currTotalEntries = 0;
  UInt_t nFile = 0;
  for(auto const & file : fileList){
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    TFile* inFile_p = new TFile(file.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");
    std::string inDataSetName = inConfig_p->GetValue("INDATASET", "");

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(isMC){
      if(inDataSetName.size() == 0 || !sHandler.Init(inDataSetName)){
        std::cout << "GDJMCNTUPLEPREPROC ERROR - Given input \'" << file << "\' contains INDATASET \'" << inDataSetName << "\' that is not valid. return 1" << std::endl;

        inFile_p->Close();
        delete inFile_p;

        outFile_p->Close();
        delete outFile_p;

        return 1;
      }

      sampleTag_ = sHandler.GetTag();
      sampleWeight_ = tagToWeight[sampleTag_];
    }

    for(unsigned int bI = 0; bI < listOfBranchesHLT.size(); ++bI){
      inTree_p->SetBranchAddress(listOfBranchesHLT[bI].c_str(), hltVect[bI]);
    }

    for(unsigned int bI = 0; bI < listOfBranchesHLTPre.size(); ++bI){
      inTree_p->SetBranchAddress(listOfBranchesHLTPre[bI].c_str(), hltPreVect[bI]);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    inTree_p->SetBranchAddress("runNumber", &runNumber_);
    inTree_p->SetBranchAddress("eventNumber", &eventNumber_);
    inTree_p->SetBranchAddress("lumiBlock", &lumiBlock_);

    if(isMC){
      inTree_p->SetBranchAddress("pthat", &pthat_);

      inTree_p->SetBranchAddress("treePartonPt", treePartonPt_);
      inTree_p->SetBranchAddress("treePartonEta", treePartonEta_);
      inTree_p->SetBranchAddress("treePartonPhi", treePartonPhi_);
      inTree_p->SetBranchAddress("treePartonId", treePartonId_);
    }

    inTree_p->SetBranchAddress("actualInteractionsPerCrossing", &actualInteractionsPerCrossing_);
    inTree_p->SetBranchAddress("averageInteractionsPerCrossing", &averageInteractionsPerCrossing_);

    inTree_p->SetBranchAddress("nvert", &nvert_);
    inTree_p->SetBranchAddress("vert_x", &vert_x_p);
    inTree_p->SetBranchAddress("vert_y", &vert_y_p);
    inTree_p->SetBranchAddress("vert_z", &vert_z_p);
    inTree_p->SetBranchAddress("vert_type", &vert_type_p);
    inTree_p->SetBranchAddress("vert_ntrk", &vert_ntrk_p);

    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);
    inTree_p->SetBranchAddress("fcalA_et_Cos2", &fcalA_et_Cos2_);
    inTree_p->SetBranchAddress("fcalC_et_Cos2", &fcalC_et_Cos2_);
    inTree_p->SetBranchAddress("fcalA_et_Sin2", &fcalA_et_Sin2_);
    inTree_p->SetBranchAddress("fcalC_et_Sin2", &fcalC_et_Sin2_);

    inTree_p->SetBranchAddress("fcalA_et_Cos3", &fcalA_et_Cos3_);
    inTree_p->SetBranchAddress("fcalC_et_Cos3", &fcalC_et_Cos3_);
    inTree_p->SetBranchAddress("fcalA_et_Sin3", &fcalA_et_Sin3_);
    inTree_p->SetBranchAddress("fcalC_et_Sin3", &fcalC_et_Sin3_);

    inTree_p->SetBranchAddress("fcalA_et_Cos4", &fcalA_et_Cos4_);
    inTree_p->SetBranchAddress("fcalC_et_Cos4", &fcalC_et_Cos4_);
    inTree_p->SetBranchAddress("fcalA_et_Sin4", &fcalA_et_Sin4_);
    inTree_p->SetBranchAddress("fcalC_et_Sin4", &fcalC_et_Sin4_);

    inTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi_);
    inTree_p->SetBranchAddress("evtPlane3Phi", &evtPlane3Phi_);
    inTree_p->SetBranchAddress("evtPlane4Phi", &evtPlane4Phi_);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(getTracks){
      inTree_p->SetBranchAddress("ntrk", &ntrk_);
      inTree_p->SetBranchAddress("trk_pt", &trk_pt_p);
      inTree_p->SetBranchAddress("trk_eta", &trk_eta_p);
      inTree_p->SetBranchAddress("trk_phi", &trk_phi_p);

      inTree_p->SetBranchAddress("trk_charge", &trk_charge_p);
      inTree_p->SetBranchAddress("trk_tight_primary", &trk_tight_primary_p);
      inTree_p->SetBranchAddress("trk_minbias", &trk_minbias_p);
      inTree_p->SetBranchAddress("trk_d0", &trk_d0_p);
      inTree_p->SetBranchAddress("trk_z0", &trk_z0_p);
      inTree_p->SetBranchAddress("trk_vz", &trk_vz_p);
      inTree_p->SetBranchAddress("trk_theta", &trk_theta_p);
      inTree_p->SetBranchAddress("trk_nPixelHits", &trk_nPixelHits_p);
      inTree_p->SetBranchAddress("trk_nSCTHits", &trk_nSCTHits_p);
      inTree_p->SetBranchAddress("trk_nBlayerHits", &trk_nBlayerHits_p);
    }

    if(isMC && getTruthParticle){
      inTree_p->SetBranchAddress("truth_n", &truth_n_);
      inTree_p->SetBranchAddress("truth_charge", &truth_charge_p);
      inTree_p->SetBranchAddress("truth_e", &truth_e_p);
      inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
      inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
      inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
      inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);
      inTree_p->SetBranchAddress("truth_type", &truth_type_p);
      inTree_p->SetBranchAddress("truth_origin", &truth_origin_p);
      inTree_p->SetBranchAddress("truth_status", &truth_status_p);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(getR2jets){
      inTree_p->SetBranchAddress("akt2hi_jet_n", &akt2hi_jet_n_);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_m", &akt2hi_em_xcalib_jet_m_p);
      //inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_clean", &akt2hi_em_xcalib_jet_clean_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_pt", &akt2hi_em_xcalib_jet_pt_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_uncorrpt", &akt2hi_em_xcalib_jet_uncorrpt_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_eta", &akt2hi_em_xcalib_jet_eta_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_uncorreta", &akt2hi_em_xcalib_jet_uncorreta_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_phi", &akt2hi_em_xcalib_jet_phi_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_e", &akt2hi_em_xcalib_jet_e_p);
      inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_uncorre", &akt2hi_em_xcalib_jet_uncorre_p);
      inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_pt", &akt2hi_constit_xcalib_jet_pt_p);
      inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_eta", &akt2hi_constit_xcalib_jet_eta_p);
      inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_phi", &akt2hi_constit_xcalib_jet_phi_p);
      inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_e", &akt2hi_constit_xcalib_jet_e_p);
      inTree_p->SetBranchAddress("akt2hi_double_calib_jet_pt", &akt2hi_double_calib_jet_pt_p);
      if(isMC) inTree_p->SetBranchAddress("akt2hi_truthpos", &akt2hi_truthpos_p);
    }

    inTree_p->SetBranchAddress("akt4hi_jet_n", &akt4hi_jet_n_);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_m", &akt4hi_em_xcalib_jet_m_p);
    //inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_clean", &akt4hi_em_xcalib_jet_clean_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_b4ExtraCalib_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_uncorrpt", &akt4hi_em_xcalib_jet_uncorrpt_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_uncorreta", &akt4hi_em_xcalib_jet_uncorreta_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_e", &akt4hi_em_xcalib_jet_e_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_uncorre", &akt4hi_em_xcalib_jet_uncorre_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_pt", &akt4hi_constit_xcalib_jet_pt_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_eta", &akt4hi_constit_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_phi", &akt4hi_constit_xcalib_jet_phi_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_e", &akt4hi_constit_xcalib_jet_e_p);
    inTree_p->SetBranchAddress("akt4hi_double_calib_jet_pt", &akt4hi_double_calib_jet_pt_p);
    if(isMC) inTree_p->SetBranchAddress("akt4hi_truthpos", &akt4hi_truthpos_p);

    if(isMC){
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_0", &akt4hi_em_xcalib_jet_pt_sys_JES_0_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_1", &akt4hi_em_xcalib_jet_pt_sys_JES_1_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_2", &akt4hi_em_xcalib_jet_pt_sys_JES_2_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_3", &akt4hi_em_xcalib_jet_pt_sys_JES_3_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_4", &akt4hi_em_xcalib_jet_pt_sys_JES_4_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_5", &akt4hi_em_xcalib_jet_pt_sys_JES_5_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_6", &akt4hi_em_xcalib_jet_pt_sys_JES_6_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_7", &akt4hi_em_xcalib_jet_pt_sys_JES_7_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_8", &akt4hi_em_xcalib_jet_pt_sys_JES_8_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_9", &akt4hi_em_xcalib_jet_pt_sys_JES_9_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_10", &akt4hi_em_xcalib_jet_pt_sys_JES_10_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_11", &akt4hi_em_xcalib_jet_pt_sys_JES_11_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_12", &akt4hi_em_xcalib_jet_pt_sys_JES_12_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_13", &akt4hi_em_xcalib_jet_pt_sys_JES_13_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_14", &akt4hi_em_xcalib_jet_pt_sys_JES_14_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_15", &akt4hi_em_xcalib_jet_pt_sys_JES_15_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_16", &akt4hi_em_xcalib_jet_pt_sys_JES_16_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JES_17", &akt4hi_em_xcalib_jet_pt_sys_JES_17_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_0", &akt4hi_em_xcalib_jet_pt_sys_JER_0_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_1", &akt4hi_em_xcalib_jet_pt_sys_JER_1_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_2", &akt4hi_em_xcalib_jet_pt_sys_JER_2_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_3", &akt4hi_em_xcalib_jet_pt_sys_JER_3_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_4", &akt4hi_em_xcalib_jet_pt_sys_JER_4_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_5", &akt4hi_em_xcalib_jet_pt_sys_JER_5_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_6", &akt4hi_em_xcalib_jet_pt_sys_JER_6_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_7", &akt4hi_em_xcalib_jet_pt_sys_JER_7_p);
      inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt_sys_JER_8", &akt4hi_em_xcalib_jet_pt_sys_JER_8_p);
    } 

    if(getR10jets){
      if(!isPP){
        inTree_p->SetBranchAddress("akt10hi_jet_n", &akt10hi_jet_n_);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_m", &akt10hi_em_xcalib_jet_m_p);
        //inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_clean", &akt10hi_em_xcalib_jet_clean_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_pt", &akt10hi_em_xcalib_jet_pt_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_uncorrpt", &akt10hi_em_xcalib_jet_uncorrpt_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_eta", &akt10hi_em_xcalib_jet_eta_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_uncorreta", &akt10hi_em_xcalib_jet_uncorreta_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_phi", &akt10hi_em_xcalib_jet_phi_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_e", &akt10hi_em_xcalib_jet_e_p);
        inTree_p->SetBranchAddress("akt10hi_em_xcalib_jet_uncorre", &akt10hi_em_xcalib_jet_uncorre_p);
        inTree_p->SetBranchAddress("akt10hi_constit_xcalib_jet_pt", &akt10hi_constit_xcalib_jet_pt_p);
        inTree_p->SetBranchAddress("akt10hi_constit_xcalib_jet_eta", &akt10hi_constit_xcalib_jet_eta_p);
        inTree_p->SetBranchAddress("akt10hi_constit_xcalib_jet_phi", &akt10hi_constit_xcalib_jet_phi_p);
        inTree_p->SetBranchAddress("akt10hi_constit_xcalib_jet_e", &akt10hi_constit_xcalib_jet_e_p);
        inTree_p->SetBranchAddress("akt10hi_double_calib_jet_pt", &akt10hi_double_calib_jet_pt_p);
        if(isMC) inTree_p->SetBranchAddress("akt10hi_truthpos", &akt10hi_truthpos_p);
      }

      inTree_p->SetBranchAddress("akt2to10hi_jet_n", &akt2to10hi_jet_n_);
      inTree_p->SetBranchAddress("akt2to10hi_em_xcalib_jet_pt", &akt2to10hi_em_xcalib_jet_pt_p);
      inTree_p->SetBranchAddress("akt2to10hi_em_xcalib_jet_eta", &akt2to10hi_em_xcalib_jet_eta_p);
      inTree_p->SetBranchAddress("akt2to10hi_em_xcalib_jet_phi", &akt2to10hi_em_xcalib_jet_phi_p);
      inTree_p->SetBranchAddress("akt2to10hi_em_xcalib_jet_e", &akt2to10hi_em_xcalib_jet_e_p);
      inTree_p->SetBranchAddress("akt2to10hi_em_xcalib_jet_m", &akt2to10hi_em_xcalib_jet_m_p);
      if(isMC) inTree_p->SetBranchAddress("akt2to10hi_truthpos", &akt2to10hi_truthpos_p);
    }

    inTree_p->SetBranchAddress("photon_n", &photon_n_);
    inTree_p->SetBranchAddress("photon_pt", &photon_pt_b4ExtraCalib_p);
    if(isMC){
      inTree_p->SetBranchAddress("photon_pt_sys1", &photon_pt_sys1_p);
      inTree_p->SetBranchAddress("photon_pt_sys2", &photon_pt_sys2_p);
      inTree_p->SetBranchAddress("photon_pt_sys3", &photon_pt_sys3_p);
      inTree_p->SetBranchAddress("photon_pt_sys4", &photon_pt_sys4_p);
    }
    inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
    inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
    inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
    inTree_p->SetBranchAddress("photon_loose", &photon_loose_p);
    inTree_p->SetBranchAddress("photon_isem", &photon_isem_p);
    inTree_p->SetBranchAddress("photon_convFlag", &photon_convFlag_p);
    inTree_p->SetBranchAddress("photon_Rconv", &photon_Rconv_p);
    inTree_p->SetBranchAddress("photon_etcone20", &photon_etcone20_p);
    inTree_p->SetBranchAddress("photon_etcone30", &photon_etcone30_p);
    inTree_p->SetBranchAddress("photon_etcone40", &photon_etcone40_p);
    inTree_p->SetBranchAddress("photon_topoetcone20", &photon_topoetcone20_p);
    inTree_p->SetBranchAddress("photon_topoetcone30", &photon_topoetcone30_p);
    inTree_p->SetBranchAddress("photon_topoetcone40", &photon_topoetcone40_p);
    //inTree_p->SetBranchAddress("photon_etcone20_ptcorrected", &photon_etcone20_ptcorrected_p);
    //inTree_p->SetBranchAddress("photon_etcone30_ptcorrected", &photon_etcone30_ptcorrected_p);
    //inTree_p->SetBranchAddress("photon_etcone40_ptcorrected", &photon_etcone40_ptcorrected_p);
    inTree_p->SetBranchAddress("photon_Rhad1", &photon_Rhad1_p);
    inTree_p->SetBranchAddress("photon_Rhad", &photon_Rhad_p);
    inTree_p->SetBranchAddress("photon_e277", &photon_e277_p);
    inTree_p->SetBranchAddress("photon_Reta", &photon_Reta_p);
    inTree_p->SetBranchAddress("photon_Rphi", &photon_Rphi_p);
    inTree_p->SetBranchAddress("photon_weta1", &photon_weta1_p);
    inTree_p->SetBranchAddress("photon_weta2", &photon_weta2_p);
    inTree_p->SetBranchAddress("photon_wtots1", &photon_wtots1_p);
    inTree_p->SetBranchAddress("photon_f1", &photon_f1_p);
    inTree_p->SetBranchAddress("photon_f3", &photon_f3_p);
    inTree_p->SetBranchAddress("photon_fracs1", &photon_fracs1_p);
    inTree_p->SetBranchAddress("photon_DeltaE", &photon_DeltaE_p);
    inTree_p->SetBranchAddress("photon_Eratio", &photon_Eratio_p);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(isMC){
      if(getR2jets){
        inTree_p->SetBranchAddress("akt2_truth_jet_n", &akt2_truth_jet_n_);
        inTree_p->SetBranchAddress("akt2_truth_jet_pt", &akt2_truth_jet_pt_p);
        inTree_p->SetBranchAddress("akt2_truth_jet_eta", &akt2_truth_jet_eta_p);
        inTree_p->SetBranchAddress("akt2_truth_jet_phi", &akt2_truth_jet_phi_p);
        inTree_p->SetBranchAddress("akt2_truth_jet_e", &akt2_truth_jet_e_p);
        inTree_p->SetBranchAddress("akt2_truth_jet_m", &akt2_truth_jet_m_p);
        inTree_p->SetBranchAddress("akt2_truth_jet_partonid", &akt2_truth_jet_partonid_p);
        inTree_p->SetBranchAddress("akt2_truth_jet_recopos", &akt2_truth_jet_recopos_p);
      }

      inTree_p->SetBranchAddress("akt4_truth_jet_n", &akt4_truth_jet_n_);
      inTree_p->SetBranchAddress("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_e", &akt4_truth_jet_e_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_m", &akt4_truth_jet_m_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_partonid", &akt4_truth_jet_partonid_p);
      inTree_p->SetBranchAddress("akt4_truth_jet_recopos", &akt4_truth_jet_recopos_p);

      if(getR10jets){
        if(!isPP){
          inTree_p->SetBranchAddress("akt10_truth_jet_n", &akt10_truth_jet_n_);
          inTree_p->SetBranchAddress("akt10_truth_jet_pt", &akt10_truth_jet_pt_p);
          inTree_p->SetBranchAddress("akt10_truth_jet_eta", &akt10_truth_jet_eta_p);
          inTree_p->SetBranchAddress("akt10_truth_jet_phi", &akt10_truth_jet_phi_p);
          inTree_p->SetBranchAddress("akt10_truth_jet_e", &akt10_truth_jet_e_p);
          inTree_p->SetBranchAddress("akt10_truth_jet_m", &akt10_truth_jet_m_p);
          inTree_p->SetBranchAddress("akt10_truth_jet_partonid", &akt10_truth_jet_partonid_p);
          inTree_p->SetBranchAddress("akt10_truth_jet_recopos", &akt10_truth_jet_recopos_p);
        }

        inTree_p->SetBranchAddress("akt2to10_truth_jet_n", &akt2to10_truth_jet_n_);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_pt", &akt2to10_truth_jet_pt_p);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_eta", &akt2to10_truth_jet_eta_p);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_phi", &akt2to10_truth_jet_phi_p);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_e", &akt2to10_truth_jet_e_p);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_m", &akt2to10_truth_jet_m_p);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_partonid", &akt2to10_truth_jet_partonid_p);
        inTree_p->SetBranchAddress("akt2to10_truth_jet_recopos", &akt2to10_truth_jet_recopos_p);
      }
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    // START EVENT LOOP!

    ULong64_t temp_nEntries = inTree_p->GetEntries();
    if(isTest) temp_nEntries = 10;
    const ULong64_t nEntries = temp_nEntries;
    std::cout << "total entry = " << nEntries << std::endl;
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
      if(currTotalEntries%nDiv == 0) std::cout << " Entry " << currTotalEntries << "/" << totalNEntries << "... (File " << nFile << "/" << fileList.size() << ")"  << std::endl;
      inTree_p->GetEntry(entry);

      int icentBin = 0;
      if(!isPP){
        cent_ = centTable.GetCent(fcalA_et_ + fcalC_et_);
        ncollWeight_ = ncollWeights[cent_];
        for (; cent_>=centBins[icentBin+1] && icentBin<nCENTBINS; ++icentBin);
      }
      else ncollWeight_ = 1.0;

      fullWeight_ = sampleWeight_*ncollWeight_;
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      if(isMC){
        truthPhotonPt_ = -999.;     
        truthPhotonPhi_ = -999.;      
        truthPhotonEta_ = -999.;
        truthPhotonIso2_ = -999.;
        truthPhotonIso3_ = -999.;
        truthPhotonIso4_ = -999.;

        truthOut_charge_p->clear();
        truthOut_pt_p->clear();
        truthOut_eta_p->clear();
        truthOut_phi_p->clear();
        truthOut_e_p->clear();
        truthOut_pdg_p->clear();
        truthOut_type_p->clear();
        truthOut_origin_p->clear();
        truthOut_status_p->clear();
        truthOut_n_ = 0;

        for(unsigned int tI = 0; tI < truth_pt_p->size(); ++tI){
          if(truth_status_p->at(tI) != 1) continue; 

          truthOut_charge_p->push_back(truth_charge_p->at(tI));
          truthOut_pt_p->push_back(truth_pt_p->at(tI));
          truthOut_eta_p->push_back(truth_eta_p->at(tI));
          truthOut_phi_p->push_back(truth_phi_p->at(tI));
          truthOut_e_p->push_back(truth_e_p->at(tI));
          truthOut_pdg_p->push_back(truth_pdg_p->at(tI));
          truthOut_type_p->push_back(truth_type_p->at(tI));
          truthOut_origin_p->push_back(truth_origin_p->at(tI));
          truthOut_status_p->push_back(truth_status_p->at(tI));   
          ++truthOut_n_; 

          if(truth_type_p->at(tI) != 14) continue;
          if(truth_origin_p->at(tI) != 37) continue;
          if(truth_pdg_p->at(tI) != 22) continue;
          if(truth_status_p->at(tI) != 1) continue;

          float genEtSum2 = 0;
          float genEtSum3 = 0;
          float genEtSum4 = 0;
          for(unsigned int tI2 = 0; tI2 <truth_pt_p->size() ; ++tI2){
            if( tI2==tI) continue;
            if( truth_status_p->at(tI2) != 1) continue;
            if( (truth_type_p->at(tI2) == 14 && truth_origin_p->at(tI2)==37 && truth_pdg_p->at(tI2) == 22)) continue;

            //TLorentzVector temp;
            //temp.SetPtEtaPhiE(truth_pt_p->at(tI2), truth_eta_p->at(tI2), truth_phi_p->at(tI2), truth_e_p->at(tI2));
            ROOT::Math::PtEtaPhiEVector temp(truth_pt_p->at(tI2), truth_eta_p->at(tI2), truth_phi_p->at(tI2), truth_e_p->at(tI2));
            Float_t dR = getDR(truth_eta_p->at(tI), truth_phi_p->at(tI), truth_eta_p->at(tI2), truth_phi_p->at(tI2));
            if(dR < 0.4) genEtSum4 += temp.Et();
            if(dR < 0.3) genEtSum3 += temp.Et();
            if(dR < 0.2) genEtSum2 += temp.Et();
          }

          if(truthPhotonPt_ > 0){
            if(doGlobalDebug) std::cout << "WARNING: MORE THAN ONE GEN LEVEL PHOTON FOUND IN FILE, ENTRY: " << file << ", " << entry << std::endl;
            if(truth_pt_p->at(tI) > truthPhotonPt_){
              truthPhotonPt_ = truth_pt_p->at(tI);
              truthPhotonPhi_ = truth_phi_p->at(tI);
              truthPhotonEta_ = truth_eta_p->at(tI);
              truthPhotonIso2_ = genEtSum2; 
              truthPhotonIso3_ = genEtSum3; 
              truthPhotonIso4_ = genEtSum4; 
            }
          }
          else{
            truthPhotonPt_ = truth_pt_p->at(tI);
            truthPhotonPhi_ = truth_phi_p->at(tI);
            truthPhotonEta_ = truth_eta_p->at(tI);
            truthPhotonIso2_ = genEtSum2; 
            truthPhotonIso3_ = genEtSum3; 
            truthPhotonIso4_ = genEtSum4; 
          }
        }//truth particle loop
        //if(truthPhotonPt_ > 0 && truthPhotonPhi_ < -100) 
        if(isTest)
          std::cout << "truthPt = " << truthPhotonPt_ << ", truthEta = " << truthPhotonEta_ << ", truthPhi = " << truthPhotonPhi_ << ", truthIso = " << truthPhotonIso3_ << std::endl;

      }//if(isMC)?

      //////////////////////////////////////////////////
      // photon extra scale
      photon_pt_extraCalib_p->clear();
      //photon energy systematic is only varied in MC! not data!
      if(isMC){
        photon_pt_sys1_extraCalib_p->clear(); 
        photon_pt_sys2_extraCalib_p->clear(); 
        photon_pt_sys3_extraCalib_p->clear(); 
        photon_pt_sys4_extraCalib_p->clear(); 
      }

      for(unsigned int pI = 0; pI < photon_pt_b4ExtraCalib_p->size(); ++pI){
        float pt_b4Calib = photon_pt_b4ExtraCalib_p->at(pI);
        //photon_pt_b4ExtraCalib_p->push_back(pt_b4Calib);
        float calibFactor = 1.;
        float phoAbsEta = abs(photon_eta_p->at(pI));
        if(phoAbsEta <= 1.37) 
          calibFactor = 1./f_phoExtCalib_barrel[icentBin]->Eval(pt_b4Calib);
        else if(phoAbsEta >= 1.52 && phoAbsEta <= 2.37) 
          calibFactor = 1./f_phoExtCalib_endcap[icentBin]->Eval(pt_b4Calib);
        photon_pt_extraCalib_p->push_back(calibFactor*pt_b4Calib); 
        if(isTest)
          std::cout << "photon pt before and after : " << pt_b4Calib << ", " << photon_pt_extraCalib_p->at(pI) << ", calibration factor = " << calibFactor << std::endl;
        
        //std::cout << "pt before : "  << photon_pt_p->at(pI) << std::endl;
        if(isMC){
          photon_pt_sys1_extraCalib_p->push_back(photon_pt_sys1_p->at(pI)*calibFactor); 
          photon_pt_sys2_extraCalib_p->push_back(photon_pt_sys2_p->at(pI)*calibFactor); 
          photon_pt_sys3_extraCalib_p->push_back(photon_pt_sys3_p->at(pI)*calibFactor); 
          photon_pt_sys4_extraCalib_p->push_back(photon_pt_sys4_p->at(pI)*calibFactor); 
          //photon_pt_sys1_extraCalib_p->at(pI) = photon_pt_sys1_p->at(pI)*calibFactor; 
          //photon_pt_sys2_extraCalib_p->at(pI) = photon_pt_sys2_p->at(pI)*calibFactor; 
          //photon_pt_sys3_extraCalib_p->at(pI) = photon_pt_sys3_p->at(pI)*calibFactor; 
          //photon_pt_sys4_extraCalib_p->at(pI) = photon_pt_sys4_p->at(pI)*calibFactor; 

          if(isTest){
            float temp_sys1 = photon_pt_sys1_p->at(pI);
            float temp_sys2 = photon_pt_sys2_p->at(pI);
            float temp_sys3 = photon_pt_sys3_p->at(pI);
            float temp_sys4 = photon_pt_sys4_p->at(pI);
            std::cout << "photon:: before calibration pt, sys1,2,3,4 = " << photon_pt_b4ExtraCalib_p->at(pI) << ", " << temp_sys1  << ", " << temp_sys2  << ", " << temp_sys3  << ", " << temp_sys4  << std::endl; 
            std::cout << "photon:: after  calibration pt, sys1,2,3,4 = " << photon_pt_extraCalib_p->at(pI) << ", " << photon_pt_sys1_extraCalib_p->at(pI)  << ", " << photon_pt_sys2_extraCalib_p->at(pI)  << ", " << photon_pt_sys3_extraCalib_p->at(pI)  << ", " << photon_pt_sys4_extraCalib_p->at(pI) << std::endl; 
          }
        }

      }

      //////////////////////////////////////////////////
      // jet extra scale 
      akt4hi_em_xcalib_jet_pt_extraCalib_p->clear();
      //jet energy systematic is only varied in MC! not in data!
      if(isMC){
        akt4hi_em_xcalib_jet_pt_sys_JES_0_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_1_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_2_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_3_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_4_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_5_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_6_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_7_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_8_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_9_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_10_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_11_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_12_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_13_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_14_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_15_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_16_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JES_17_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_0_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_1_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_2_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_3_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_4_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_5_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_6_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_7_extraCalib_p->clear();
        akt4hi_em_xcalib_jet_pt_sys_JER_8_extraCalib_p->clear();
      }

      for(unsigned int pI = 0; pI < akt4hi_em_xcalib_jet_pt_b4ExtraCalib_p->size(); ++pI){
        float pt_b4Calib = akt4hi_em_xcalib_jet_pt_b4ExtraCalib_p->at(pI);
        float calibFactor = 1./f_jetExtCalib[icentBin]->Eval(pt_b4Calib);
        akt4hi_em_xcalib_jet_pt_extraCalib_p->push_back(calibFactor*pt_b4Calib); 
        if(isTest)
          std::cout << "jet pt before and after : " << pt_b4Calib << ", " << akt4hi_em_xcalib_jet_pt_extraCalib_p->at(pI) << ", calibration factor = " << calibFactor << std::endl;
        //std::cout << "pt before : "  << photon_pt_p->at(pI) << std::endl;
        if(isMC){
          akt4hi_em_xcalib_jet_pt_sys_JES_0_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_0_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_1_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_1_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_2_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_2_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_3_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_3_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_4_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_4_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_5_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_5_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_6_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_6_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_7_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_7_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_8_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_8_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_9_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_9_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_10_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_10_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_11_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_11_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_12_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_12_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_13_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_13_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_14_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_14_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_15_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_15_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_16_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_16_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JES_17_extraCalib_p->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_17_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_0_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_0_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_1_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_1_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_2_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_2_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_3_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_3_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_4_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_4_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_5_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_5_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_6_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_6_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_7_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_7_p->at(pI)*calibFactor);
          akt4hi_em_xcalib_jet_pt_sys_JER_8_extraCalib_p ->push_back(akt4hi_em_xcalib_jet_pt_sys_JES_8_p->at(pI)*calibFactor);

          if(isTest){
            float temp_sys1 = akt4hi_em_xcalib_jet_pt_sys_JES_0_p->at(pI);
            float temp_sys2 = akt4hi_em_xcalib_jet_pt_sys_JES_1_p->at(pI);
            float temp_sys3 = akt4hi_em_xcalib_jet_pt_sys_JER_0_p->at(pI);
            float temp_sys4 = akt4hi_em_xcalib_jet_pt_sys_JER_1_p->at(pI);
            std::cout << "jet:: before calibration pt, sys1,2,3,4 = " << pt_b4Calib << ", " << temp_sys1  << ", " << temp_sys2  << ", " << temp_sys3  << ", " << temp_sys4  << std::endl; 
            std::cout << "jet:: after  calibration pt, sys1,2,3,4 = " << akt4hi_em_xcalib_jet_pt_extraCalib_p->at(pI) << ", " << akt4hi_em_xcalib_jet_pt_sys_JES_0_extraCalib_p->at(pI)  << ", " << akt4hi_em_xcalib_jet_pt_sys_JES_1_extraCalib_p->at(pI)  << ", " << akt4hi_em_xcalib_jet_pt_sys_JER_0_extraCalib_p->at(pI)  << ", " << akt4hi_em_xcalib_jet_pt_sys_JER_1_extraCalib_p->at(pI) << std::endl; 
          }
        }

      }

      outTree_p->Fill();
      ++currTotalEntries;
    }

    inFile_p->Close();
    delete inFile_p;

    ++nFile;
  }

  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;

  TEnv outConfig;
  for(auto const & val : configMap){  
    for(unsigned int vI = 0; vI < val.second.size(); ++vI){
      std::string configVal = val.first;
      if(val.second.size() != 1) configVal = configVal + "_" + std::to_string(vI);
      outConfig.SetValue(configVal.c_str(), val.second[vI].c_str());
    }    
  }

  outConfig.Write("config", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete inConfig_p;

  std::cout << "GDJMCNTUPLEPREPROC COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjNtuplePreProc_phoTaggedJetRaa.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjNtuplePreProc_phoTaggedJetRaa(argv[1]);
  return retVal;
}
