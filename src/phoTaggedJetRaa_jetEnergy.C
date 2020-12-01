//Author: Chris McGinn (2020.02.19)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TDirectoryFile.h"
#include "TEnv.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/photonUtil.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"
#include "include/returnFileList.h"
#include "/direct/usatlas+u/goyeonju/phoTaggedJetRaa/include/yjUtility.h"

void fillTH1(TH1F* inHist_p, Float_t fillVal, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal, weight);
  }
  return;
}

void fillTH2(TH2F* inHist_p, Float_t fillVal1, Float_t fillVal2, Float_t weight = -1.0)
{
  if(weight < 0) inHist_p->Fill(fillVal1, fillVal2);
  else{
    if(inHist_p->GetSumw2()->fN == 0) inHist_p->Sumw2();
    inHist_p->Fill(fillVal1, fillVal2, weight);
  }
  return;
}

int phoTaggedJetRaa_jetEnergy(std::string inConfigFileName)
{
  const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101
  TRandom3* randGen_p = new TRandom3(randSeed);

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INDIRNAME",
                                              "VERSION",
                                              "CENTFILENAME",
                          "ISPP",
					      "ISMC",
					      "JETR",
					      "GAMMAEXCLUSIONDR",
		"PHOTONSELECTION",
        "PHOGENMATCHINGDR",
        "PHOISOCONESIZE",
        "GENISOCUT",
        "ISOCUT",
        "DOPTCORRECTEDISO",
        "DOCENTCORRECTEDISO",
        "BKGISOGAP",
					      "CENTBINS",
                          "DOPTBINSSUBINCONFIG",
                          "GAMMAPTBINSSUB",
					      "NGAMMAPTBINSSUB",
					      "GAMMAPTBINSSUBLOW",
					      "GAMMAPTBINSSUBHIGH",
					      "GAMMAPTBINSSUBDOLOG",
                          "NPHOETABINS",
                          "ETABINS_I",
                          "ETABINS_F",
					      "NJTETABINS",
					      "JTETABINSLOW",
					      "JTETABINSHIGH",
					      "JTETABINSDOABS",
                          "DOJTPTBINSCONFIG",
					      "NJTPTBINS",
					      "JTPTBINS",
					      "JTPTBINSLOW",
					      "JTPTBINSHIGH",
					      "JTPTBINSDOLOG",
					      "NDPHIBINS",
					      "DPHIBINSLOW",
					      "DPHIBINSHIGH",
					      "GAMMAJTDPHI"  
					      };

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
 
  std::string inCentFileName = config_p->GetValue("CENTFILENAME", "");
  std::string version = config_p->GetValue("VERSION", "temp");
  std::string inGRLFileName = config_p->GetValue("GRLFILENAME", "");

  ////photon configuration
  const bool doPtCorrectedIso = config_p->GetValue("DOPTCORRECTEDISO", 1);
  const bool doCentCorrectedIso = config_p->GetValue("DOCENTCORRECTEDISO", 1);
  const Float_t isoCut = config_p->GetValue("ISOCUT", 3);
  //const Float_t genIsoCut = config_p->GetValue("GENISOCUT", 5);
  //const Float_t phoGenMatchingDR = config_p->GetValue("PHOGENMATCHINGDR", 0.2);
  const Float_t phoIsoConeSize = config_p->GetValue("PHOISOCONESIZE", 3);
  std::string label_phoIsoConeSize = Form("%d",(int)(phoIsoConeSize));

  const std::string nMaxEvtStr = config_p->GetValue("NEVT", "");
  ULong64_t nMaxEvt = 0;
  if(nMaxEvtStr.size() != 0){
    nMaxEvt = std::stol(nMaxEvtStr);
  }

  const int jetR = config_p->GetValue("JETR", 4);
  if(jetR != 2 && jetR != 4){
    std::cout << "Given parameter jetR, \'" << jetR << "\' is not \'2\' or \'4\'. return 1" << std::endl;
    return 1;
  }
  const std::string jetRStr = prettyString(((double)jetR)/10., 1, false);
  const double gammaExclusionDR = config_p->GetValue("GAMMAEXCLUSIONDR", 0.7);
  
  const bool isPP = config_p->GetValue("ISPP", 1);
  const bool isMC = config_p->GetValue("ISMC", 0);
  if(!isMC) return 1;

  ////////////////////////////////////////
  // output file name
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + version); // check dated output subdir exists; if not create

  std::string systStr = "PP";
  if(!isPP) systStr = "PbPb";

  std::string outFileName = "output/" + version + "/phoTagJetRaa_jetEnergy_" + systStr + "Data_" + version + ".root ";
  if(isMC)
    outFileName = "output/" + version + "/phoTagJetRaa_jetEnergy_" + systStr + "MC_" + version + ".root ";
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  ///////////////////////////////////////////////////////////////
  // Centrality binning
  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();

  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
  if(!isMC && !check.checkFileExt(inGRLFileName, "xml")) return 1; // GRL File xml

  const Int_t nMaxSubBins = 60;
  const Int_t nMaxCentBins = 10;
  const Int_t nMaxPtBins = 200;
  const Int_t nMaxEtaPhiBins = 100;
  Int_t nCentBins = 1;
  
  std::vector<int> centBins;
  std::vector<std::string> centBinsStr = {systStr};
  std::map<std::string, std::string> binsToLabelStr;
  if(!isPP){
    centBins = strToVectI(config_p->GetValue("CENTBINS", "0,10,30,80"));
    nCentBins = centBins.size()-1;

    if(!goodBinning(inConfigFileName, nMaxCentBins, nCentBins, "CENTBINS")) return 1;

    centBinsStr.clear();
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      centBinsStr.push_back("Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]));

      binsToLabelStr[centBinsStr[cI]] = std::to_string(centBins[cI]) + "-" + std::to_string(centBins[cI+1]) + "%";
    }
  }
  else binsToLabelStr[centBinsStr[0]] = "pp";
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  ///////////////////////////////////////////////////////////////
  //photon pT SUB bins handling
  const bool doPtBinSubInConfig = config_p->GetValue("DOPTBINSSUBINCONFIG", 0);
  const Int_t nGammaPtBinsSub = config_p->GetValue("NGAMMAPTBINSSUB", 4);
  if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaPtBinsSub, "NGAMMAPTBINSSUB")) return 1;
  const Float_t gammaPtBinsSubLow = config_p->GetValue("GAMMAPTBINSSUBLOW", 80.0);
  const Float_t gammaPtBinsSubHigh = config_p->GetValue("GAMMAPTBINSSUBHIGH", 280.0);
  const Bool_t gammaPtBinsSubDoLog = config_p->GetValue("GAMMAPTBINSSUBDOLOG", 1);

  Double_t gammaPtBinsSub[nMaxSubBins+1];
  if(!doPtBinSubInConfig){
      if(gammaPtBinsSubDoLog) getLogBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
      else getLinBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
  } else{
      std::vector<int> ptBins_vec = strToVectI(config_p->GetValue("GAMMAPTBINSSUB", "50,60,70,90,130,1000"));
      Int_t tempPtbinSize =  ptBins_vec.size();
      for(Int_t ipt=0;ipt<tempPtbinSize;++ipt){
          gammaPtBinsSub[ipt] = ptBins_vec.at(ipt);
      }
  }

  std::vector<std::string> gammaPtBinsSubStr;
  for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
    gammaPtBinsSubStr.push_back(Form("GammaPt%dto%d",(int)gammaPtBinsSub[pI],(int)gammaPtBinsSub[pI+1]));
    ReplaceStringInPlace(gammaPtBinsSubStr[pI], ".", "p");
    binsToLabelStr[gammaPtBinsSubStr[pI]] = prettyString(gammaPtBinsSub[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[pI+1], 1, false);
  }
  gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(nGammaPtBinsSub));

  binsToLabelStr[gammaPtBinsSubStr[gammaPtBinsSubStr.size()-1]] = prettyString(gammaPtBinsSub[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false);
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////////
    //photon eta SUB bins handling
  const int nPhoEtaBins = config_p->GetValue("NPHOETABINS", 2);
  std::vector<float> etaBins_i = strToVectF(config_p->GetValue("ETABINS_I", ""));
  std::vector<float> etaBins_f = strToVectF(config_p->GetValue("ETABINS_F", ""));
  std::vector<std::string> etaBinsStr;
  for(Int_t ieta=0;ieta<nPhoEtaBins;++ieta){
    etaBinsStr.push_back(Form("Eta%.2fto%.2f",etaBins_i[ieta],etaBins_f[ieta])); 
    ReplaceStringInPlace(etaBinsStr[ieta], ".", "p");
  }
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////////
    //jet pT bins handling
  const bool doJtPtBinInConfig = config_p->GetValue("DOJTPTBINSCONFIG", 1);
  Int_t nJtPtBins_temp = config_p->GetValue("NJTPTBINS", 20);
  if(!goodBinning(inConfigFileName, nMaxPtBins, nJtPtBins_temp, "NJTPTBINS")) return 1;
  const Float_t jtPtBinsLow = config_p->GetValue("JTPTBINSLOW", 30.0);
  const Float_t jtPtBinsHigh = config_p->GetValue("JTPTBINSHIGH", 230.0);
  const Bool_t jtPtBinsDoLog = config_p->GetValue("JTPTBINSDOLOG", 1);
  Double_t jtPtBins[nMaxPtBins+1];
  if(!doJtPtBinInConfig){
      if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins_temp, jtPtBins);
      else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins_temp, jtPtBins);
  } else{
      std::vector<int> ptBins_vec = strToVectI(config_p->GetValue("JTPTBINS", "40,50,60,70,90,130,200"));
      nJtPtBins_temp =  ptBins_vec.size()-1;
        cout << " === JET PT BINS DEFINED IN CONFIG FILES === "  << endl;
      for(Int_t ipt=0;ipt<nJtPtBins_temp+1;++ipt){
          jtPtBins[ipt] = ptBins_vec.at(ipt);
          cout << "JET PT BIN " << ipt << " = " <<  jtPtBins[ipt] << endl;
      }
  }
  const Int_t nJtPtBins = nJtPtBins_temp; 

  std::string jtPtBinsGlobalStr = "GlobalJtPt0";
  std::string jtPtBinsGlobalLabel = prettyString(jtPtBinsLow,1,false) + " < p_{T,jet} < " + prettyString(jtPtBinsHigh,1,false);
  binsToLabelStr[jtPtBinsGlobalStr] = jtPtBinsGlobalLabel;   

  std::string multiJtCutGlobalStr = "MultiJt0";
  std::string multiJtCutGlobalLabel = "N_{Jet,Reco.} >= 2";
  binsToLabelStr[multiJtCutGlobalStr] = multiJtCutGlobalLabel;   

  std::vector<std::string> jtPtBinsStr;
  for(Int_t pI = 0; pI < nJtPtBins; ++pI){
    if(pI < 10) jtPtBinsStr.push_back("JtPt0" + std::to_string(pI));
    else jtPtBinsStr.push_back("JtPt" + std::to_string(pI));

    binsToLabelStr[jtPtBinsStr[pI]] = prettyString(jtPtBins[pI], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[pI+1], 1, false);
  }
  jtPtBinsStr.push_back("JtPt" + std::to_string(nJtPtBins));
  binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] = prettyString(jtPtBins[0], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[nJtPtBins], 1, false);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////////
    //jet eta bins handling
  const Int_t nJtEtaBins = config_p->GetValue("NJTETABINS", 14);
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nJtEtaBins, "NJTETABINS")) return 1;
  const Float_t jtEtaBinsLow = config_p->GetValue("JTETABINSLOW", 0.0);
  const Float_t jtEtaBinsHigh = config_p->GetValue("JTETABINSHIGH", 2.8);
  const Bool_t jtEtaBinsDoAbs = config_p->GetValue("JTETABINSDOABS", 1);
  if(jtEtaBinsDoAbs && jtEtaBinsLow < 0){
    std::cout << "ERROR - config \'" << inConfigFileName << "\' contains jtEtaBinsLow \'" << jtEtaBinsLow << "\' less than 0 despite requested jtEtaBinsDoAbs \'" << jtEtaBinsDoAbs << "\'. return 1" << std::endl;
    return 1;
  }
  Double_t jtEtaBins[nMaxEtaPhiBins+1];
  getLinBins(jtEtaBinsLow, jtEtaBinsHigh, nJtEtaBins, jtEtaBins);

    ///////////////////////////////////////////////////////////////
    //photon-jet dphi bins handling
  const Double_t gammaJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);  

  const int nRatio = 100; 
  const float ratioMin = 0; 
  const float ratioMax = 3; 

  const int nJet = 100; 
  const float jetMin = 40; 
  const float jetMax = 200; 
  ///////////////////////////////////////////////
  // output file and histogram definition
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* pthat_p = nullptr;
  TH1F* pthat_Unweighted_p = nullptr;
  TH1F* centrality_p = nullptr;
  TH1F* centrality_Unweighted_p = nullptr;

  TH1F* h1F_genMatchedRecoPt[nMaxCentBins]; // reco gen-matched
  TH1F* h1F_recoMatchedGenPt[nMaxCentBins]; // gen reco-matched
  TH1F* h1F_recoMatchedGenPt_finerBin[nMaxCentBins]; // gen reco-matched
  TH1F* h1F_recoPt_noMatching[nMaxCentBins]; // reco not gen-matched
  TH1F* h1F_genPt_noMatching[nMaxCentBins]; // gen not reco-matched
  TH2F* h2F_reco_over_gen_ratio_vs_genPt[nMaxCentBins]; 
  TH2F* h2F_reco_over_gen_ratio_vs_genPt_quarkJet[nMaxCentBins]; 
  TH2F* h2F_reco_over_gen_ratio_vs_genPt_gluonJet[nMaxCentBins]; 
  TH1F* h1F_genPt_vs_quarkFraction[nMaxCentBins]; // gen not reco-matched
  TH2F* h2F_genPt_recoPt[nMaxCentBins]; //response matrix
  TH2F* h2F_genPt_recoPt_split[nMaxCentBins]; //response matrix
  TH1F* h1F_genMatchedRecoPt_split[nMaxCentBins]; // reco gen-matched

  if(isMC){
    pthat_p = new TH1F(("pthat_" + systStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    pthat_Unweighted_p = new TH1F(("pthat_Unweighted_" + systStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    centerTitles({pthat_p, pthat_Unweighted_p});
  }
  
  if(!isPP){
    centrality_p = new TH1F(("centrality_" + systStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
    centerTitles(centrality_p);

    if(isMC){
      centrality_Unweighted_p = new TH1F(("centrality_Unweighted_" + systStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
      centerTitles(centrality_Unweighted_p);
    }
  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


  ////////////////////////////////////////////////
  // define histogram 
  for(Int_t cI = 0; cI < nCentBins; ++cI){
      h1F_genMatchedRecoPt[cI] = new TH1F(("h1F_genMatchedRecoPt_" + centBinsStr[cI]).c_str(), Form(";Truth-matched Reco p_{T}^{jet};Entries%s",""), nJtPtBins, jtPtBins); 
      h1F_genMatchedRecoPt_split[cI] = new TH1F(("h1F_genMatchedRecoPt_split_" + centBinsStr[cI]).c_str(), Form(";Truth-matched Reco p_{T}^{jet};Entries%s",""), nJtPtBins, jtPtBins); 
      h1F_recoMatchedGenPt[cI] = new TH1F(("h1F_recoMatchedGenPt_" + centBinsStr[cI]).c_str(), Form(";Reco-matched Gen p_{T}^{jet};Entries%s",""), nJtPtBins, jtPtBins); 
      h1F_recoMatchedGenPt_finerBin[cI] = new TH1F(("h1F_recoMatchedGenPt_finerBin_" + centBinsStr[cI]).c_str(), Form(";Reco-matched Gen p_{T}^{jet};Entries%s",""), nJet, jetMin, jetMax); 
      h1F_recoPt_noMatching[cI] = new TH1F(("h1F_recoPt_noMatching_" + centBinsStr[cI]).c_str(), Form(";Reco p_{T}^{jet};Entries%s",""), nJtPtBins, jtPtBins);
      h1F_genPt_noMatching[cI] = new TH1F(("h1F_genPt_noMatching_" + centBinsStr[cI]).c_str(), Form(";Gen p_{T}^{jet};Entries%s",""), nJtPtBins, jtPtBins); 
      h2F_genPt_recoPt[cI] = new TH2F(("h2F_genPt_recoPt_" + centBinsStr[cI]).c_str(), Form(";Reco p_{T}^{jet}%s;Gen p_{T}^{jet}",""), nJtPtBins, jtPtBins, nJtPtBins, jtPtBins); // x-axis: reco, y-axis: gen
      h2F_genPt_recoPt_split[cI] = new TH2F(("h2F_genPt_recoPt_split_" + centBinsStr[cI]).c_str(), Form(";Reco p_{T}^{jet}%s;Gen p_{T}^{jet}",""), nJtPtBins, jtPtBins, nJtPtBins, jtPtBins); // x-axis: reco, y-axis: gen
      h2F_reco_over_gen_ratio_vs_genPt[cI] = new TH2F(("h2F_reco_over_gen_ratio_vs_genPt_" + centBinsStr[cI]).c_str(), Form(";Truth p_{T}^{jet};p_{T}^{reco jet}/p_{T}^{truth jet}%s",""), nJet, jetMin, jetMax, nRatio, ratioMin, ratioMax);

      h2F_reco_over_gen_ratio_vs_genPt_quarkJet[cI] = new TH2F(("h2F_reco_over_gen_ratio_vs_genPt_quarkJet_" + centBinsStr[cI]).c_str(), Form(";Truth p_{T}^{jet};p_{T}^{reco jet}/p_{T}^{truth jet}%s",""), nJet, jetMin, jetMax, nRatio, ratioMin, ratioMax);
      h2F_reco_over_gen_ratio_vs_genPt_gluonJet[cI] = new TH2F(("h2F_reco_over_gen_ratio_vs_genPt_gluonJet_" + centBinsStr[cI]).c_str(), Form(";Truth p_{T}^{jet};p_{T}^{reco jet}/p_{T}^{truth jet}%s",""), nJet, jetMin, jetMax, nRatio, ratioMin, ratioMax);
      h1F_genPt_vs_quarkFraction[cI] = new TH1F(("h1F_genPt_quarkFraction_" + centBinsStr[cI]).c_str(), Form(";Reco p_{T}^{jet}%s;Fraction",""), nJet, jetMin, jetMax); 
  }


  ////////////////////////////////////////////////
  // retrieve input file 
  const std::string inDirStr = config_p->GetValue("INDIRNAME", "");
  std::vector<std::string> fileList = returnFileList(inDirStr, ".root");
  if(fileList.size() == 0){
      std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' contains no root files. return 1" << std::endl;
      return 1;
  }

  TChain* inTree_p = new TChain("gammaJetTree_p");
  for(auto const & file : fileList){
      std::cout << file.c_str() << std::endl;
      inTree_p->Add(file.c_str());
  }

  inTree_p->SetBranchStatus("*", 0);
  outFile_p->cd();
  
  //Grab the hltbranches for some basic prescale checks
  std::vector<std::string> listOfBranches = getVectBranchList(inTree_p);
  std::vector<std::string> hltList;
  std::vector<std::string> hltListPres;
  std::string hltStr = "HLT_";
  std::string prescaleStr = "_prescale";

  //HLTLists are built
  for(auto const & branchStr : listOfBranches){
    if(branchStr.size() < hltStr.size()) continue;
    if(!isStrSame(branchStr.substr(0, hltStr.size()), hltStr)) continue;

    if(branchStr.find(prescaleStr) != std::string::npos) hltListPres.push_back(branchStr);
    else hltList.push_back(branchStr);
  }

  bool allHLTPrescalesFound = true;
  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    std::string modStr = hltList[hI] + prescaleStr;
    if(!vectContainsStr(modStr, &hltListPres)){
      allHLTPrescalesFound = false;
      std::cout << "HLT " << hltList[hI] << " has no prescale. return 1" << std::endl;
    }
  }
  if(!allHLTPrescalesFound) return 1;

  float hltPrescaleDelta = 0.01;
  std::vector<bool*> hltVect;
  std::vector<float*> hltPrescaleVect;
  Int_t runNumber;
  UInt_t lumiBlock;
  Float_t pthat;
  Float_t sampleWeight;
  Float_t ncollWeight;
  Float_t fullWeight;
  Float_t fcalA_et, fcalC_et;
  std::vector<float>* vert_z_p=nullptr;
  
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;

  Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;
  
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;  
  std::vector<unsigned int>* photon_isem_p=nullptr;
  std::vector<float>* photon_etcone_p=nullptr;
  
  std::vector<float>* aktRhi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_uncorrpt_p=nullptr;
  std::vector<float>* aktRhi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_uncorreta_p=nullptr;
  std::vector<float>* aktRhi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* aktRhi_em_xcalib_jet_phi_p=nullptr;
  std::vector<int>* aktRhi_truthpos_p=nullptr;

  std::vector<float>* aktR_truth_jet_pt_p=nullptr;
  std::vector<float>* aktR_truth_jet_eta_p=nullptr;
  std::vector<float>* aktR_truth_jet_phi_p=nullptr;
  std::vector<int>* aktR_truth_jet_partonid_p=nullptr;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inTree_p->SetBranchStatus("*", 0);

  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    hltVect.push_back(new bool(false));
    hltPrescaleVect.push_back(new float(0.0));

    inTree_p->SetBranchStatus(hltList[hI].c_str(), 1);
    inTree_p->SetBranchStatus(hltListPres[hI].c_str(), 1);
  }  

  inTree_p->SetBranchStatus("runNumber", 1);
  inTree_p->SetBranchStatus("lumiBlock", 1);

  if(isMC){
    inTree_p->SetBranchStatus("pthat", 1);
    inTree_p->SetBranchStatus("sampleWeight", 1);
    if(!isPP) inTree_p->SetBranchStatus("ncollWeight", 1);
    inTree_p->SetBranchStatus("fullWeight", 1);

    inTree_p->SetBranchStatus("truth_pt", 1);
    inTree_p->SetBranchStatus("truth_eta", 1);
    inTree_p->SetBranchStatus("truth_phi", 1);
    inTree_p->SetBranchStatus("truth_pdg", 1);

    inTree_p->SetBranchStatus("truthPhotonPt", 1);
    inTree_p->SetBranchStatus("truthPhotonEta", 1);
    inTree_p->SetBranchStatus("truthPhotonPhi", 1);
    //inTree_p->SetBranchStatus(("truthPhotonIso"+label_phoIsoConeSize).c_str(), 1);
  }
  
  if(!isPP){
    inTree_p->SetBranchStatus("fcalA_et", 1);
    inTree_p->SetBranchStatus("fcalC_et", 1);
  }

  inTree_p->SetBranchStatus("vert_z", 1);
  
  inTree_p->SetBranchStatus("photon_pt", 1);
  inTree_p->SetBranchStatus("photon_eta", 1);
  inTree_p->SetBranchStatus("photon_phi", 1);
  inTree_p->SetBranchStatus("photon_tight", 1);
  inTree_p->SetBranchStatus("photon_isem", 1);
  inTree_p->SetBranchStatus(("photon_etcone"+label_phoIsoConeSize+"0").c_str(), 1);
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorrpt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_pt").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorreta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_eta").c_str(), 1);
  inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), 1);
  
  if(isMC){
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_truthpos").c_str(), 1);
    
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_pt").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_eta").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_phi").c_str(), 1);
    inTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "_truth_jet_partonid").c_str(), 1);
  }
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    inTree_p->SetBranchAddress(hltList[hI].c_str(), hltVect[hI]);
    inTree_p->SetBranchAddress(hltListPres[hI].c_str(), hltPrescaleVect[hI]);
  }

  inTree_p->SetBranchAddress("runNumber", &runNumber);
  inTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
  if(isMC){
    inTree_p->SetBranchAddress("pthat", &pthat);
    inTree_p->SetBranchAddress("sampleWeight", &sampleWeight);
    if(!isPP) inTree_p->SetBranchAddress("ncollWeight", &ncollWeight);
    inTree_p->SetBranchAddress("fullWeight", &fullWeight);

    inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
    inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
    inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);
    inTree_p->SetBranchAddress("truth_pdg", &truth_pdg_p);

    inTree_p->SetBranchAddress("truthPhotonPt", &truthPhotonPt);
    inTree_p->SetBranchAddress("truthPhotonEta", &truthPhotonEta);
    inTree_p->SetBranchAddress("truthPhotonPhi", &truthPhotonPhi);
    //inTree_p->SetBranchAddress(("truthPhotonIso"+label_phoIsoConeSize).c_str(), &truthPhotonIso);
  }

  if(!isPP){
    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
  }

  inTree_p->SetBranchAddress("vert_z", &vert_z_p);
  
  inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
  inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
  inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
  inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
  inTree_p->SetBranchAddress("photon_isem", &photon_isem_p);
  inTree_p->SetBranchAddress(("photon_etcone"+label_phoIsoConeSize+"0").c_str(), &photon_etcone_p);

  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), &aktRhi_em_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorrpt").c_str(), &aktRhi_em_xcalib_jet_uncorrpt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_pt").c_str(), &aktRhi_constit_xcalib_jet_pt_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), &aktRhi_em_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_uncorreta").c_str(), &aktRhi_em_xcalib_jet_uncorreta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_constit_xcalib_jet_eta").c_str(), &aktRhi_constit_xcalib_jet_eta_p);
  inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), &aktRhi_em_xcalib_jet_phi_p);

  if(isMC){
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_truthpos").c_str(), &aktRhi_truthpos_p);    

    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_pt").c_str(), &aktR_truth_jet_pt_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_eta").c_str(), &aktR_truth_jet_eta_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_phi").c_str(), &aktR_truth_jet_phi_p);
    inTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "_truth_jet_partonid").c_str(), &aktR_truth_jet_partonid_p);
  }


  Double_t recoJtPtMin = 100000.;
  
  ULong64_t nEntriesTemp = inTree_p->GetEntries();
  if(doGlobalDebug) nEntriesTemp = 2000;
  if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, nMaxEvt);
  const ULong64_t nEntries = nEntriesTemp;
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);
  std::cout << "Processing " << nEntries << " events..." << std::endl;

  bool didOneFireMiss = false;
  std::vector<int> skippedCent;
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
 
  /////////////////////////////////////////////////////////////////////
  // EVENT LOOP 
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

    bool isEvenEvt = false;
    if(entry%2 == 0) isEvenEvt = true;

    double vert_z = vert_z_p->at(0);
    vert_z /= 1000.;
    if(vert_z <= -15. || vert_z >= 15.) continue;      

    if(!didOneFireMiss && !isMC){//only check this once per input
      //check at least one of the purported selection triggers fired
      bool oneFire = false;
      for(unsigned int hI = 0; hI < hltVect.size(); ++hI){
	if(*(hltVect[hI])){
	  oneFire = true;
	  break;
	}
      }

      if(!oneFire){
	std::cout << "WARNING - YOU HAVE EVENTS w/ NO TRIGGERS!!!" << std::endl;
	didOneFireMiss = true;
      }
    }
    
    //Check prescale is 1 in case i made a mistake on first unprescaled
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    for(unsigned int hI = 0; hI < hltPrescaleVect.size(); ++hI){
      if(TMath::Abs((*(hltPrescaleVect[hI])) - 1.0) > hltPrescaleDelta){
	std::cout << "WARNING - prescale for \'" << hltList[hI] << "\' has non-unity value, \'" << (*(hltPrescaleVect[hI])) << "\'." << std::endl;
      }
    }
      if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    Int_t centPos = -1;
    Double_t cent = -1;
    if(!isPP){
      cent = centTable.GetCent(fcalA_et + fcalC_et);
      centPos = ghostPos(centBins, cent, true, doGlobalDebug);
    }
    else centPos = 0;

    if(centPos < 0){
      bool vectContainsCent = vectContainsInt((Int_t)cent, &skippedCent);

      if(!vectContainsCent){
	std::cout << "phoTaggedJetRaa_jetEnergy Warning - Skipping centrality \'" << (Int_t)cent << "\' as given centrality binning is \'" << centBins[0] << "-" << centBins[centBins.size()-1] << "\'. if this is incorrect please fix." << std::endl;
	skippedCent.push_back((Int_t)cent);
      }
      continue;
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
      
    if(!isPP){
      fillTH1(centrality_p, cent, fullWeight);
      if(isMC) centrality_Unweighted_p->Fill(cent);
    }

    if(isMC){
      fillTH1(pthat_p, pthat, fullWeight);
      pthat_Unweighted_p->Fill(pthat);
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    double leadingPhoPt = 0;
    int leadingPhoIndex = -1;
    /////////////////////////////////////////////////////////////////////
    // PHOTON LOOP
    for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
      if(photon_pt_p->at(pI) < gammaPtBinsSub[0]) continue;
      if(photon_pt_p->at(pI) >= gammaPtBinsSub[nGammaPtBinsSub]) continue;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      Float_t etaValMain = TMath::Abs(photon_eta_p->at(pI));
      if(etaValMain <= etaBins_i[0]) continue;
      if(etaValMain >= etaBins_f[nPhoEtaBins-1]) continue;
      if(etaValMain>= 1.37 &&  etaValMain < 1.52) continue;

      if(leadingPhoPt < photon_pt_p->at(pI)){
          leadingPhoPt = photon_pt_p->at(pI);
          leadingPhoIndex = pI;
      }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

      //if(isMC){
      //    if(truthPhotonPt<=0) continue; //prompt photons
      //    if(truthPhotonIso>genIsoCut) continue; // truth isolation condition
      //    if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) > phoGenMatchingDR) continue;

      //    if(leadingPhoPt_genMatchedReco < photon_pt_p->at(pI)){
      //        leadingPhoPt_genMatchedReco = photon_pt_p->at(pI);
      //        leadingPhoIndex_genMatchedReco = pI;
      //    }
      //}
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    if(leadingPhoIndex == -1) continue;
    //if(leadingPhoIndex_genMatchedReco == -1) isGoodGenMatchedRecoPhoton = false;

    //Int_t ptPos = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photon_pt_p->at(leadingPhoIndex), true, doGlobalDebug);
    //Int_t etaPos = ghostPos(nGammaEtaBinsSub, gammaEtaBinsSub, etaValSub, true, doGlobalDebug);
    int tempEtaPos = -1;
    Float_t etaValMain = TMath::Abs(photon_eta_p->at(leadingPhoIndex));
    for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
        if(etaValMain>=etaBins_i[eI] && etaValMain<etaBins_f[eI]) tempEtaPos=eI;
    }
    if(tempEtaPos==-1) continue; //eta cut

    //Int_t ptPos_genMatchedReco = -1;
    //int tempEtaPos_genMatchedReco = -1;
    //if(isGoodGenMatchedRecoPhoton){
    //    ptPos_genMatchedReco = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photon_pt_p->at(leadingPhoIndex_genMatchedReco), true, doGlobalDebug);
    //    Float_t etaValMain_genMatchedReco = TMath::Abs(photon_eta_p->at(leadingPhoIndex_genMatchedReco));
    //    for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
    //        if(etaValMain_genMatchedReco>=etaBins_i[eI] && etaValMain_genMatchedReco<etaBins_f[eI]) tempEtaPos_genMatchedReco=eI;
    //    }
    //}
    //if(tempEtaPos_genMatchedReco == -1) isGoodGenMatchedRecoPhoton = false;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    /////////////////// photon isolation correction
    float correctedIso = photon_etcone_p->at(leadingPhoIndex);

    if(doPtCorrectedIso && doCentCorrectedIso)
        correctedIso = getCorrectedPhotonIsolation(isPP, photon_etcone_p->at(leadingPhoIndex), photon_pt_p->at(leadingPhoIndex), photon_eta_p->at(leadingPhoIndex), cent);

    if(!(photon_tight_p->at(leadingPhoIndex)==1 && correctedIso < isoCut)) continue;


    /////////////////////////////////////////////////////////////////////
    // Truth jet loop
    for(unsigned int tjI = 0; tjI < aktR_truth_jet_pt_p->size(); ++tjI){
        if(aktR_truth_jet_eta_p->at(tjI) <= jtEtaBinsLow) continue;
        if(aktR_truth_jet_eta_p->at(tjI) >= jtEtaBinsHigh) continue;

        Float_t dR = getDR(aktR_truth_jet_eta_p->at(tjI), aktR_truth_jet_phi_p->at(tjI), photon_eta_p->at(leadingPhoIndex), photon_phi_p->at(leadingPhoIndex));
        if(dR < gammaExclusionDR) continue;
        Float_t dPhi = TMath::Abs(getDPHI(aktR_truth_jet_phi_p->at(tjI), photon_phi_p->at(leadingPhoIndex)));
        if(dPhi < gammaJtDPhiCut) continue;

        fillTH1(h1F_genPt_noMatching[centPos], aktR_truth_jet_pt_p->at(tjI), fullWeight);

    }// truth jet loop

    /////////////////////////////////////////////////////////////////////
    // Reco JET LOOP 
    for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
        if(aktRhi_em_xcalib_jet_eta_p->at(jI) <= jtEtaBinsLow) continue;
        if(aktRhi_em_xcalib_jet_eta_p->at(jI) >= jtEtaBinsHigh) continue;

        Float_t dR = getDR(aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), photon_eta_p->at(leadingPhoIndex), photon_phi_p->at(leadingPhoIndex));
        if(dR < gammaExclusionDR) continue;
        Float_t dPhi = TMath::Abs(getDPHI(aktRhi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(leadingPhoIndex)));
        if(dPhi < gammaJtDPhiCut) continue;

        fillTH1(h1F_recoPt_noMatching[centPos], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);

        int truthPos = aktRhi_truthpos_p->at(jI);
        if(truthPos < 0) continue; 
        if(aktR_truth_jet_partonid_p->at(truthPos) == -1) continue; 

        double truthJetPt = aktR_truth_jet_pt_p->at(truthPos); 
        float truthJetEta = aktR_truth_jet_eta_p->at(truthPos); 
        if(truthJetEta <= jtEtaBinsLow) continue;
        if(truthJetEta >= jtEtaBinsHigh) continue;

        double recoJetPt = aktRhi_em_xcalib_jet_pt_p->at(jI); 

        fillTH1(h1F_genMatchedRecoPt[centPos], recoJetPt, fullWeight);
        fillTH1(h1F_recoMatchedGenPt[centPos], truthJetPt, fullWeight);
        fillTH1(h1F_recoMatchedGenPt_finerBin[centPos], truthJetPt, fullWeight);
        fillTH2(h2F_reco_over_gen_ratio_vs_genPt[centPos], truthJetPt, recoJetPt/truthJetPt, fullWeight);
        fillTH2(h2F_genPt_recoPt[centPos], recoJetPt, truthJetPt, fullWeight); // x-axis: reco, y-axis: gen
        if(isEvenEvt) 
            fillTH2(h2F_genPt_recoPt_split[centPos], recoJetPt, truthJetPt, fullWeight); // x-axis: reco, y-axis: gen
        else 
            fillTH1(h1F_genMatchedRecoPt_split[centPos], recoJetPt, fullWeight);

        if(aktR_truth_jet_partonid_p->at(truthPos) == 21){     
            fillTH2(h2F_reco_over_gen_ratio_vs_genPt_gluonJet[centPos], truthJetPt, recoJetPt/truthJetPt, fullWeight);
        } else {
            fillTH2(h2F_reco_over_gen_ratio_vs_genPt_quarkJet[centPos], truthJetPt, recoJetPt/truthJetPt, fullWeight);
            fillTH1(h1F_genPt_vs_quarkFraction[centPos], truthJetPt, fullWeight);
        }

    }// reco jet loop
      
  } // END OF EVENT LOOP
  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
 // inFile_p->Close();
 // delete inFile_p;

  outFile_p->cd();


  ///////////////////////////////////////////////////////////
  // Write histograms in the output file
  for(Int_t cI = 0; cI < nCentBins; ++cI){
      h1F_genMatchedRecoPt[cI]->Write("", TObject::kOverwrite);
      h1F_recoMatchedGenPt[cI]->Write("", TObject::kOverwrite);
      h1F_recoMatchedGenPt_finerBin[cI]->Write("", TObject::kOverwrite);
      h1F_recoPt_noMatching[cI]->Write("", TObject::kOverwrite);
      h1F_genPt_noMatching[cI]->Write("", TObject::kOverwrite);
      h2F_reco_over_gen_ratio_vs_genPt[cI]->Write("", TObject::kOverwrite);
      h2F_reco_over_gen_ratio_vs_genPt_quarkJet[cI]->Write("", TObject::kOverwrite);
      h2F_reco_over_gen_ratio_vs_genPt_gluonJet[cI]->Write("", TObject::kOverwrite);
      h1F_genPt_vs_quarkFraction[cI]->Write("", TObject::kOverwrite);
      h2F_genPt_recoPt[cI]->Write("", TObject::kOverwrite);
      h2F_genPt_recoPt_split[cI]->Write("", TObject::kOverwrite);
      h1F_genMatchedRecoPt_split[cI]->Write("", TObject::kOverwrite);
  }

  //runNumber_p->Write("", TObject::kOverwrite);

  if(!isPP){
    centrality_p->Write("", TObject::kOverwrite);
    if(isMC) centrality_Unweighted_p->Write("", TObject::kOverwrite);
  }

  if(isMC){
    pthat_p->Write("", TObject::kOverwrite);
    pthat_Unweighted_p->Write("", TObject::kOverwrite);
  }


    ///////////////////////////////////////////////////////////
    // delete histograms

  for(Int_t cI = 0; cI < nCentBins; ++cI){
     delete h1F_genMatchedRecoPt[cI];
     delete h1F_recoMatchedGenPt[cI];
     delete h1F_recoMatchedGenPt_finerBin[cI];
     delete h1F_recoPt_noMatching[cI];
     delete h1F_genPt_noMatching[cI];
     delete h2F_reco_over_gen_ratio_vs_genPt[cI];
     delete h2F_reco_over_gen_ratio_vs_genPt_quarkJet[cI];
     delete h2F_reco_over_gen_ratio_vs_genPt_gluonJet[cI];
     delete h1F_genPt_vs_quarkFraction[cI];
     delete h2F_genPt_recoPt[cI];
     delete h2F_genPt_recoPt_split[cI];
     delete h1F_genMatchedRecoPt_split[cI];
  }
  
  //delete runNumber_p;

  if(isMC){
    delete pthat_p;
    delete pthat_Unweighted_p;
  }

  if(!isPP){
    delete centrality_p;
    if(isMC) delete centrality_Unweighted_p;
  }

  config_p->SetValue("RECOJTPTMIN", prettyString(recoJtPtMin, 1, false).c_str());
  config_p->Write("config", TObject::kOverwrite);

  TEnv labelEnv;
  for(auto const & lab : binsToLabelStr){
    labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
  }
  labelEnv.Write("label", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  std::cout << "phoTaggedJetRaa_jetEnergy COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/phoTaggedJetRaa_jetEnergy.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += phoTaggedJetRaa_jetEnergy(argv[1]);
  return retVal;
}
