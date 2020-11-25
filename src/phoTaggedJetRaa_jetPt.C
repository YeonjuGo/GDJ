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
//#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TTree.h"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"

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

int phoTaggedJetRaa_jetPt(std::string inConfigFileName, bool doBkgPhoton = false)
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
					                          "MIXFILENAME",
                                              "ISPP",
					                          "ISMC",
					      "DOMIX",
                          "NMIX",
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

  std::vector<std::string> mixParams = {"DOMIXCENT",
					"NMIXCENTBINS",
					"MIXCENTBINSLOW",
					"MIXCENTBINSHIGH",
					"DOMIXPSI2",
					"NMIXPSI2BINS",
					"MIXPSI2BINSLOW",
					"MIXPSI2BINSHIGH",
					"DOMIXVZ",
					"NMIXVZBINS",
					"MIXVZBINSLOW",
					"MIXVZBINSHIGH"};

  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
 
  std::string inCentFileName = config_p->GetValue("CENTFILENAME", "");
  std::string version = config_p->GetValue("VERSION","temp");
  std::string inGRLFileName = config_p->GetValue("GRLFILENAME", "");
  std::string inMixFileName = config_p->GetValue("MIXFILENAME", "");

  const int photonSelection = config_p->GetValue("PHOTONSELECTION", 0); // 0: tight, isolated, 1:tight, non-isolated, 2:non-tight, isolated, 3: non-tight, non-isolated, 4: non-isolated, 5: non-tight
  const bool doPtCorrectedIso = config_p->GetValue("DOPTCORRECTEDISO", 1);
  const bool doCentCorrectedIso = config_p->GetValue("DOCENTCORRECTEDISO", 1);
  const Float_t bkgIsoGap = config_p->GetValue("BKGISOGAP", 2);
  const Float_t isoCut = config_p->GetValue("ISOCUT", 3);
  const Float_t genIsoCut = config_p->GetValue("GENISOCUT", 5);
  const Float_t phoIsoConeSize = config_p->GetValue("PHOISOCONESIZE", 3);
  std::string label_phoIsoConeSize = Form("%d",(int)(phoIsoConeSize));
  std::cout << "label_phoIsoConeSize = " << label_phoIsoConeSize << std::endl;
  const Float_t phoGenMatchingDR = config_p->GetValue("PHOGENMATCHINGDR", 0.2);

  const unsigned long long nMix = config_p->GetValue("NMIX", 10);

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

  const bool doMix = config_p->GetValue("DOMIX", 0);

  if(doMix){
    if(!checkEnvForParams(config_p, mixParams)) return 1;
  } 
  
  //Mixing categories, temp hardcoding
  //Centrality, percent level
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const Int_t nMaxMixBins = 200;
  bool doMixCent = false;
  Int_t nMixCentBins = -1;
  Float_t mixCentBinsLow = 0;
  Float_t mixCentBinsHigh = 100;
  Double_t mixCentBins[nMaxMixBins+1];

  bool doMixPsi2 = false;
  Int_t nMixPsi2Bins = -1;
  Float_t mixPsi2BinsLow = -TMath::Pi()/2.;
  Float_t mixPsi2BinsHigh = TMath::Pi()/2.;
  Double_t mixPsi2Bins[nMaxMixBins+1];

  bool doMixVz = false;
  Int_t nMixVzBins = -1;
  Float_t mixVzBinsLow = -15.;
  Float_t mixVzBinsHigh = 15.;
  Double_t mixVzBins[nMaxMixBins+1];

  std::vector<std::vector<unsigned long long> > mixVect;
  std::vector<std::vector<unsigned long long> > keyVect;
  
  if(doMix){
    if(!check.checkFileExt(inMixFileName, ".root")) return 1;

    doMixCent = (bool)config_p->GetValue("DOMIXCENT", 0);

    if(doMixCent){      
      nMixCentBins = config_p->GetValue("NMIXCENTBINS", 10);
      if(nMixCentBins > nMaxMixBins){
	std::cout << "phoTaggedJetRaa_jetPt ERROR - centrality mixing bins \'" << nMixCentBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixCentBinsLow = config_p->GetValue("MIXCENTBINSLOW", 0.0);
      mixCentBinsHigh = config_p->GetValue("MIXCENTBINSHIGH", 80.0);
      getLinBins(mixCentBinsLow, mixCentBinsHigh, nMixCentBins, mixCentBins);

      mixVect.push_back({});
      for(Int_t mI = 0; mI < nMixCentBins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);
	keyVect.push_back({(unsigned long long)mI});
      }
    }

    doMixPsi2 = (bool)config_p->GetValue("DOMIXPSI2", 0);

    if(doMixPsi2){
      nMixPsi2Bins = config_p->GetValue("NMIXPSI2BINS", 16);
      if(nMixPsi2Bins > nMaxMixBins){
	std::cout << "phoTaggedJetRaa_jetPt ERROR - centrality mixing bins \'" << nMixPsi2Bins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixPsi2BinsLow = config_p->GetValue("MIXPSI2BINSLOW", 0.0);
      mixPsi2BinsHigh = config_p->GetValue("MIXPSI2BINSHIGH", TMath::Pi()/2.0);
      getLinBins(mixPsi2BinsLow, mixPsi2BinsHigh, nMixPsi2Bins, mixPsi2Bins);

      std::vector<std::vector<unsigned long long> > tempKeyVect;
      mixVect.push_back({});
      for(Int_t mI = 0; mI < nMixPsi2Bins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);

	if(keyVect.size() != 0){
	  for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
	    std::vector<unsigned long long> tempVect = keyVect[vI];
	    tempVect.push_back((unsigned long long)mI);
	    tempKeyVect.push_back(tempVect);
	  }
	}
	else tempKeyVect.push_back({(unsigned long long)mI});
      }

      keyVect = tempKeyVect;
    }

    doMixVz = (bool)config_p->GetValue("DOMIXVZ", 0);

    if(doMixVz){
      nMixVzBins = config_p->GetValue("NMIXVZBINS", 10);
      if(nMixVzBins > nMaxMixBins){
	std::cout << "phoTaggedJetRaa_jetPt ERROR - centrality mixing bins \'" << nMixVzBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
	return 1;
      }
      mixVzBinsLow = config_p->GetValue("MIXVZBINSLOW", -15.0);
      mixVzBinsHigh = config_p->GetValue("MIXVZBINSHIGH", 15.0);    
      getLinBins(mixVzBinsLow, mixVzBinsHigh, nMixVzBins, mixVzBins);

      std::vector<std::vector<unsigned long long> > tempKeyVect;
      mixVect.push_back({});
      for(Int_t mI = 0; mI < nMixVzBins; ++mI){
	mixVect[mixVect.size()-1].push_back(mI);

	if(keyVect.size() != 0){
	  for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
	    std::vector<unsigned long long> tempVect = keyVect[vI];
	    tempVect.push_back((unsigned long long)mI);
	    tempKeyVect.push_back(tempVect);
	  }
	}
	else tempKeyVect.push_back({(unsigned long long)mI});
      }

      keyVect = tempKeyVect;
    }

    if(doMixCent){
      std::cout << "MIXING IN CENTRALItY: " << std::endl;
      for(Int_t cI = 0; cI < nMixCentBins; ++cI){
	std::cout << " BIN " << cI << ": " << mixCentBins[cI] << "-" << mixCentBins[cI+1] << std::endl;
      }      
    }

    if(doMixPsi2){
      std::cout << "MIXING IN PSI2: " << std::endl;
      for(Int_t cI = 0; cI < nMixPsi2Bins; ++cI){
	std::cout << " BIN " << cI << ": " << mixPsi2Bins[cI] << "-" << mixPsi2Bins[cI+1] << std::endl;
      }      
    }

    if(doMixVz){
      std::cout << "MIXING IN VZ: " << std::endl;
      for(Int_t cI = 0; cI < nMixVzBins; ++cI){
	std::cout << " BIN " << cI << ": " << mixVzBins[cI] << "-" << mixVzBins[cI+1] << std::endl;
      }      
    }

  }

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  keyHandler runLumiKey("runLumiHandler");//for lumi estimation
  runLumiKey.Init({1000000, 10000});//runnumbers, then lumi
  
  keyHandler keyBoy("mixingHandler");//For Mixing
  std::map<unsigned long long, std::vector<std::vector<ROOT::Math::PtEtaPhiMVector> > > mixingMap;
  std::map<unsigned long long, unsigned long long> mixingMapCounter, signalMapCounter;
  if(doMix){
    std::vector<unsigned long long> sizes;
    for(unsigned int vI = 0; vI < mixVect.size(); ++vI){
      sizes.push_back(mixVect[vI].size()+1);
    }

    //    return 0;
    
    keyBoy.Init(sizes);    

    for(unsigned int vI = 0; vI < keyVect.size(); ++vI){
      unsigned long long key = keyBoy.GetKey(keyVect[vI]);
      mixingMap[key] = {};
      mixingMap[key].reserve(40);

      mixingMapCounter[key] = 0;
      signalMapCounter[key] = 0;
    }
  }  
  
  const bool isPP = config_p->GetValue("ISPP", 1);
  const bool isMC = config_p->GetValue("ISMC", 0);

  ////////////////////////////////////////
  // output file name
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + version); // check dated output subdir exists; if not create

  std::string systStr = "PP";
  if(!isPP) systStr = "PbPb";
  std::string capStr = "";
  cout << "doBkgPhoton = " << doBkgPhoton << endl;
  if(doBkgPhoton) capStr = version + "_bkgPhoton";
  else capStr = version + "sigPhoton";

  std::string outFileName = "output/" + version + "/phoTagJetRaa_jetPt_" + systStr + "Data_" + capStr + ".root ";
  if(isMC)
      outFileName = "output/" + version + "/phoTagJetRaa_jetPt_" + systStr + "MC_" + capStr + ".root ";

  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////////
    // Centrality binning
  if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
  if(!isMC && !check.checkFileExt(inGRLFileName, "xml")) return 1; // GRL File xml
  if(doMix && !check.checkFileExt(inMixFileName, "root")) return 1; // Check mixed file exists if mixing is requested
  
  const Int_t nMaxSubBins = 60;
  const Int_t nMaxCentBins = 10;
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


  const Int_t nMaxPtBins = 200;
  const Int_t nMaxEtaPhiBins = 100;
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

  const Int_t nPhiBins = config_p->GetValue("NPHIBINS", 32);
  if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nPhiBins, "NPHIBINS")) return 1;
  Double_t phiBins[nMaxEtaPhiBins+1];
  getLinBins(-TMath::Pi()+0.01, TMath::Pi()+0.01, nPhiBins, phiBins);

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

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
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
    //photon-jet dphi bins handling
  const Int_t nDPhiBins = config_p->GetValue("NDPHIBINS", 16);
  const Float_t dPhiBinsLow = config_p->GetValue("DPHIBINSLOW", 0.0);
  const Float_t dPhiBinsHigh = config_p->GetValue("DPHIBINSHIGH", TMath::Pi());
  const Double_t gammaJtDPhiCut = mathStringToNum(config_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const std::string gammaJtDPhiStr = "DPhi0";
  std::string gammaJtDPhiLabel = returnAllCapsString(config_p->GetValue("GAMMAJTDPHI", ""));
  if(gammaJtDPhiLabel.find("PI") != std::string::npos) gammaJtDPhiLabel.replace(gammaJtDPhiLabel.find("PI"), 2, "#pi");
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  gammaJtDPhiLabel = "|#Delta#phi_{#gamma,jet}| > " + gammaJtDPhiLabel;
  binsToLabelStr[gammaJtDPhiStr] = gammaJtDPhiLabel;

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  cout << "dphi" << endl; 

  ///////////////////////////////////////////////
  // output file
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1F* runNumber_p = nullptr;
  TH1F* pthat_p = nullptr;
  TH1F* pthat_Unweighted_p = nullptr;
  TH1F* centrality_p = nullptr;
  TH1F* centrality_Unweighted_p = nullptr;

  TH1F* h1F_nPhoton[nMaxCentBins][nPhoEtaBins];
  TH1F* h1F_nMix[nMaxCentBins][nPhoEtaBins];
  TH1F* h1F_nJetPerPhoton[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub];
  TH1F* h1F_jetPt_raw[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub];
  TH1F* h1F_jetPt_raw_mix[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub];
  TH1F* h1F_dphi_raw[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub];
  TH1F* h1F_dphi_raw_mix[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub];

  // MC
  TH1F* h1F_nPhoton_gen[nMaxCentBins][nPhoEtaBins];
  TH1F* h1F_nPhoton_genMatchedReco[nMaxCentBins][nPhoEtaBins];
  TH1F* h1F_genPho_genJet_jetPt[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub+1];
  TH1F* h1F_genPho_genJet_dphi[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub+1];
  TH1F* h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub+1];
  TH1F* h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[nMaxCentBins][nPhoEtaBins][nGammaPtBinsSub+1];
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

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
      for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
          h1F_nPhoton[cI][eI] = new TH1F(("h1F_nPhoton_" + centBinsStr[cI] + "_" + etaBinsStr[eI]).c_str(), Form(";%s;N_{#gamma} per each bin","#gamma E_{T} [GeV]"), nGammaPtBinsSub, gammaPtBinsSub);
          h1F_nMix[cI][eI] = new TH1F(("h1F_nMix_" + centBinsStr[cI] + "_" + etaBinsStr[eI]).c_str(), Form(";%s;number of mixed events","#gamma E_{T} [GeV]"), nGammaPtBinsSub, gammaPtBinsSub);

          if(isMC){
              h1F_nPhoton_gen[cI][eI] = new TH1F(("h1F_nPhoton_gen_" + centBinsStr[cI] + "_" + etaBinsStr[eI]).c_str(), Form(";%s;N_{#gamma} per each bin","#gamma E_{T} [GeV]"), nGammaPtBinsSub, gammaPtBinsSub);
              h1F_nPhoton_genMatchedReco[cI][eI] = new TH1F(("h1F_nPhoton_genMatchedReco_" + centBinsStr[cI] + "_" + etaBinsStr[eI]).c_str(), Form(";%s;N_{#gamma} per each bin","#gamma E_{T} [GeV]"), nGammaPtBinsSub, gammaPtBinsSub);
          }

          for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
              TString titleStr_jetPt = "#gamma-tagged Jet p_{T} [GeV/c]";
                TString titleStr_Y = "1/N_{#gamma} dN_{#gamma,jet}/d";
                h1F_jetPt_raw[cI][eI][pI] = new TH1F(("h1F_jetPt_raw_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";%s;%sp_{T}",titleStr_jetPt.Data(),titleStr_Y.Data()), nJtPtBins, jtPtBins);
                h1F_dphi_raw[cI][eI][pI] = new TH1F(("h1F_dphi_raw_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";#Delta#phi_{#gamma,jet};%s#Delta#phi_{#gamma,jet}",titleStr_Y.Data()), nDPhiBins, dPhiBinsLow, dPhiBinsHigh);
                
                if(isMC){
                    h1F_genPho_genJet_jetPt[cI][eI][pI] = new TH1F(("h1F_genPho_genJet_jetPt_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";%s;%sp_{T}",titleStr_jetPt.Data(),titleStr_Y.Data()), nJtPtBins, jtPtBins);
                    h1F_genPho_genJet_dphi[cI][eI][pI] = new TH1F(("h1F_genPho_genJet_dphi_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";#Delta#phi_{#gamma,jet};%s#Delta#phi_{#gamma,jet}",titleStr_Y.Data()), nDPhiBins, dPhiBinsLow, dPhiBinsHigh);
                    h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[cI][eI][pI] = new TH1F(("h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";%s;%sp_{T}",titleStr_jetPt.Data(),titleStr_Y.Data()), nJtPtBins, jtPtBins);
                    h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[cI][eI][pI] = new TH1F(("h1F_genMatchedRecoPho_genMatchedRecoJet_dphi_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";#Delta#phi_{#gamma,jet};%s#Delta#phi_{#gamma,jet}",titleStr_Y.Data()), nDPhiBins, dPhiBinsLow, dPhiBinsHigh);

                    if(pI==nGammaPtBinsSub-1){
                        h1F_genPho_genJet_jetPt[cI][eI][pI+1] = new TH1F(("h1F_genPho_genJet_jetPt_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI+1]).c_str(), Form(";%s;%sp_{T}",titleStr_jetPt.Data(),titleStr_Y.Data()), nJtPtBins, jtPtBins);
                        h1F_genPho_genJet_dphi[cI][eI][pI+1] = new TH1F(("h1F_genPho_genJet_dphi_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI+1]).c_str(), Form(";#Delta#phi_{#gamma,jet};%s#Delta#phi_{#gamma,jet}",titleStr_Y.Data()), nDPhiBins, dPhiBinsLow, dPhiBinsHigh);
                        h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[cI][eI][pI+1] = new TH1F(("h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI+1]).c_str(), Form(";%s;%sp_{T}",titleStr_jetPt.Data(),titleStr_Y.Data()), nJtPtBins, jtPtBins);
                        h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[cI][eI][pI+1] = new TH1F(("h1F_genMatchedRecoPho_genMatchedRecoJet_dphi_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI+1]).c_str(), Form(";#Delta#phi_{#gamma,jet};%s#Delta#phi_{#gamma,jet}",titleStr_Y.Data()), nDPhiBins, dPhiBinsLow, dPhiBinsHigh);
                    } 

                }

                h1F_nJetPerPhoton[cI][eI][pI] = new TH1F(("h1F_nJetPerPhoton_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), ";;N_{jet} per photon", 30, 0, 30);
                if(doMix){
                    h1F_jetPt_raw_mix[cI][eI][pI] = new TH1F(("h1F_jetPt_raw_mix_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";%s;1/N_{#gamma} dN_{#gamma,jet}/dp_{T}",titleStr_jetPt.Data()), nJtPtBins, jtPtBins);
                    h1F_dphi_raw_mix[cI][eI][pI] = new TH1F(("h1F_dphi_raw_mix_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI]).c_str(), Form(";#Delta#phi_{#gamma,jet};%s#Delta#phi_{#gamma,jet}",titleStr_Y.Data()), nDPhiBins, dPhiBinsLow, dPhiBinsHigh);
                
              }
          }
      }
  }

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
  inTree_p->SetBranchStatus("runNumber", 1);
  Int_t runMin = inTree_p->GetMinimum("runNumber");
  Int_t runMax = inTree_p->GetMaximum("runNumber");
  Int_t nRunBins = runMax - runMin;
  Float_t runMinF = ((Float_t)runMin) - 0.5;
  Float_t runMaxF = ((Float_t)runMax) + 0.5;

  outFile_p->cd();
  runNumber_p = new TH1F(("runNumber_" + systStr + "_h").c_str(), ";Run;Counts", nRunBins+1, runMinF, runMaxF);

  std::map<unsigned long long, bool> runLumiIsFired;
  std::map<int, int> runLumiCounter;
  std::map<int, int> runLumiTotal;
  if(!isMC){
    std::ifstream inFile(inGRLFileName.c_str());
    std::string tempStr;

    std::string currRunStr = "";
    
    while(std::getline(inFile, tempStr)){
      if(tempStr.find("<Run") != std::string::npos){
	tempStr.replace(0, tempStr.find(">")+1, "");
	tempStr.replace(tempStr.rfind("<"), tempStr.size(), "");
	
	currRunStr = tempStr;

	if(runLumiCounter.count(std::stoi(currRunStr)) != 0) std::cout << "Warning - counts found already for run \'" << currRunStr << "\'" << std::endl;
	else{
	  runLumiCounter[std::stoi(currRunStr)] = 0;
	  runLumiTotal[std::stoi(currRunStr)] = 0;
	}
      }
      else if(currRunStr.size() != 0 && tempStr.find("<LB") != std::string::npos){
	tempStr.replace(0, tempStr.find("\"")+1, "");
	tempStr.replace(tempStr.rfind("\""), tempStr.size(), "");
	std::string firstNumStr = tempStr.substr(0, tempStr.find("\""));
	std::string secondNumStr = tempStr;
	while(secondNumStr.find("\"") != std::string::npos){
	  secondNumStr.replace(0, secondNumStr.find("\"")+1, "");
	}
	while(firstNumStr.size() < 3){firstNumStr = "0" + firstNumStr;}
	while(secondNumStr.size() < 3){secondNumStr = "0" + secondNumStr;}

	Int_t firstNum = std::stoi(firstNumStr);
	Int_t secondNum = std::stoi(secondNumStr);
	
	for(Int_t nI = firstNum; nI <= secondNum; ++nI){
	  ++(runLumiTotal[std::stoi(currRunStr)]);
	  runLumiIsFired[runLumiKey.GetKey({(unsigned long long)std::stoi(currRunStr), (unsigned long long)nI})] = false;
	}
      }
    }
    
    inFile.close();
  }
  
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
  Float_t evtPlane2Phi;
  std::vector<float>* vert_z_p=nullptr;
  
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<int>* truth_pdg_p=nullptr;

  Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;
  Float_t truthPhotonIso;
  
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
  
  TFile* mixFile_p = nullptr;
  TTree* mixTree_p = nullptr;
  if(doMix){
    mixFile_p = new TFile(inMixFileName.c_str(), "READ");
    mixTree_p = (TTree*)mixFile_p->Get("gammaJetTree_p");

    mixTree_p->SetBranchStatus("*", 0);
    mixTree_p->SetBranchStatus("vert_z", 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), 1);
    mixTree_p->SetBranchStatus(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), 1);

    mixTree_p->SetBranchAddress("vert_z", &vert_z_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_pt").c_str(), &aktRhi_em_xcalib_jet_pt_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_eta").c_str(), &aktRhi_em_xcalib_jet_eta_p);
    mixTree_p->SetBranchAddress(("akt" + std::to_string(jetR) + "hi_em_xcalib_jet_phi").c_str(), &aktRhi_em_xcalib_jet_phi_p);

    if(!isPP){
      mixTree_p->SetBranchStatus("fcalA_et", 1);
      mixTree_p->SetBranchStatus("fcalC_et", 1);
      
      mixTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
      mixTree_p->SetBranchAddress("fcalC_et", &fcalC_et);

      if(doMixPsi2){
	mixTree_p->SetBranchStatus("evtPlane2Phi", 1);	
	mixTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi);
      }
    }

    ///////////////////////////////////////////////////////////////////
    // make vector sets for jets in mixing events 
    ULong64_t nEntriesTemp = mixTree_p->GetEntries();
    if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, (ULong64_t)nMaxEvt*10);
    const ULong64_t nMixEntries = nEntriesTemp;

    for(ULong64_t entry = 0; entry < nMixEntries; ++entry){
      mixTree_p->GetEntry(entry);

      double vert_z = vert_z_p->at(0);
      vert_z /= 1000.;
      if(vert_z <= -15. || vert_z >= 15.) continue;      
      
      Double_t cent = -1;
      unsigned long long centPos = 0;
      unsigned long long psi2Pos = 0;
      if(!isPP){
	cent = centTable.GetCent(fcalA_et + fcalC_et);
	if(cent < mixCentBinsLow || cent >= mixCentBinsHigh) continue;
	if(doMixCent) centPos = ghostPos(nMixCentBins, mixCentBins, cent);

	if(doMixPsi2){
	  if(evtPlane2Phi > TMath::Pi()/2) evtPlane2Phi -= TMath::Pi();
	  else if(evtPlane2Phi < -TMath::Pi()/2) evtPlane2Phi += TMath::Pi();

	  psi2Pos = ghostPos(nMixPsi2Bins, mixPsi2Bins, evtPlane2Phi);
	}	
      }      

      unsigned long long vzPos = 0;
      if(doMixVz) vzPos = ghostPos(nMixVzBins, mixVzBins, vert_z);
      
      std::vector<unsigned long long> eventKeyVect;
      if(doMixCent) eventKeyVect.push_back(centPos);
      if(doMixPsi2) eventKeyVect.push_back(psi2Pos);
      if(doMixVz) eventKeyVect.push_back(vzPos);
      
      unsigned long long key = keyBoy.GetKey(eventKeyVect);//, vzPos, evtPlanePos});
      
      std::vector<ROOT::Math::PtEtaPhiMVector> jets;
      for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
	if(aktRhi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
	if(aktRhi_em_xcalib_jet_pt_p->at(jI) >= jtPtBinsHigh) continue;
	if(aktRhi_em_xcalib_jet_eta_p->at(jI) <= jtEtaBinsLow) continue;
	if(aktRhi_em_xcalib_jet_eta_p->at(jI) >= jtEtaBinsHigh) continue;

	//TLorentzVector temp;
    ROOT::Math::PtEtaPhiMVector temp(aktRhi_em_xcalib_jet_pt_p->at(jI), aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), 0.0);

	jets.push_back(temp);
      }

      mixingMap[key].push_back(jets);
      ++(mixingMapCounter[key]);
      ++(signalMapCounter[key]);   
    }
    
    mixFile_p->Close();
    delete mixFile_p;

    unsigned long long minKey = 0;
    unsigned long long minimumVal = 9999999;
    for(auto const & mixes : mixingMapCounter){
      if(mixes.second < minimumVal){
	minKey = mixes.first;
	minimumVal = mixes.second;
      }
    }

    std::cout << "MINIMUM NUMBER TO MIX, CORRESPONDING KEY: " << minimumVal << ", " << minKey << std::endl;
  }

    ///////////////////////////////////////////////////////////////////
    // import branches from input file 
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
    inTree_p->SetBranchStatus(("truthPhotonIso"+label_phoIsoConeSize).c_str(), 1);
  }
  
  if(!isPP){
    inTree_p->SetBranchStatus("fcalA_et", 1);
    inTree_p->SetBranchStatus("fcalC_et", 1);
    if(doMixPsi2) inTree_p->SetBranchStatus("evtPlane2Phi", 1);
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
    inTree_p->SetBranchAddress(("truthPhotonIso"+label_phoIsoConeSize).c_str(), &truthPhotonIso);
  }

  if(!isPP){
    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et);
    if(doMixPsi2) inTree_p->SetBranchAddress("evtPlane2Phi", &evtPlane2Phi);
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
  }


  Double_t recoJtPtMin = 100000.;
  
  ULong64_t nEntriesTemp = inTree_p->GetEntries();
  if(doGlobalDebug) nEntriesTemp = 2000;
  if(nMaxEvtStr.size() != 0) nEntriesTemp = TMath::Min(nEntriesTemp, nMaxEvt);
  const ULong64_t nEntries = nEntriesTemp;
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::vector<std::vector<Double_t> > gammaCountsPerPtCent;
  for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
    gammaCountsPerPtCent.push_back({});

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      gammaCountsPerPtCent[pI].push_back(0.0);
    }
  }
  std::cout << "Processing " << nEntries << " events..." << std::endl;

  bool didOneFireMiss = false;
  std::vector<int> skippedCent;
 
  /////////////////////////////////////////////////////////////////////
  // EVENT LOOP 
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    inTree_p->GetEntry(entry);

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
	std::cout << "phoTaggedJetRaa_jetPt Warning - Skipping centrality \'" << (Int_t)cent << "\' as given centrality binning is \'" << centBins[0] << "-" << centBins[centBins.size()-1] << "\'. if this is incorrect please fix." << std::endl;
	skippedCent.push_back((Int_t)cent);
      }
      continue;
    }

    if(!isMC) fullWeight = -1.0;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    unsigned long long tempKey = runLumiKey.GetKey({(unsigned long long)runNumber, (unsigned long long)lumiBlock});    
    runLumiIsFired[tempKey] = true;
     
     
    runNumber_p->Fill(runNumber);	
    //fillTH1(runNumber_p, runNumber, fullWeight);	
    if(!isPP){
      if(isMC){
          fillTH1(centrality_p, cent, fullWeight);
          centrality_Unweighted_p->Fill(cent);
      } else{ centrality_p->Fill(cent); }
    }

    if(isMC){
      fillTH1(pthat_p, pthat, fullWeight);
      pthat_Unweighted_p->Fill(pthat);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    /////////////////////////////////////////////////////////////////////
    // Truth photon and Truth jets
    if(isMC){
        bool isGoodTruthPhoton = true;
        if(truthPhotonPt < gammaPtBinsSub[0]) isGoodTruthPhoton = false;
        if(truthPhotonPt >= gammaPtBinsSub[nGammaPtBinsSub]) isGoodTruthPhoton = false;
        int tempEtaPos_gen = -1;
        Float_t eta_gen_abs= TMath::Abs(truthPhotonEta);
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            if(eta_gen_abs>=etaBins_i[eI] && eta_gen_abs<etaBins_f[eI]) tempEtaPos_gen=eI;
        }
        if(tempEtaPos_gen<0) isGoodTruthPhoton = false;
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

        if(isGoodTruthPhoton){ 
            Int_t ptPos_genPho = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, truthPhotonPt, true, doGlobalDebug);
            fillTH1(h1F_nPhoton_gen[centPos][tempEtaPos_gen],truthPhotonPt,fullWeight);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
            for(unsigned int tjI = 0; tjI < aktR_truth_jet_pt_p->size(); ++tjI){
                if(aktR_truth_jet_eta_p->at(tjI) <= jtEtaBinsLow) continue;
                if(aktR_truth_jet_eta_p->at(tjI) >= jtEtaBinsHigh) continue;

                Float_t dPhi = TMath::Abs(getDPHI(aktR_truth_jet_phi_p->at(tjI), truthPhotonPhi));

                fillTH1(h1F_genPho_genJet_dphi[centPos][tempEtaPos_gen][ptPos_genPho], dPhi, fullWeight);
                fillTH1(h1F_genPho_genJet_dphi[centPos][tempEtaPos_gen][nGammaPtBinsSub], dPhi, fullWeight);

                if(dPhi >= gammaJtDPhiCut){
                    fillTH1(h1F_genPho_genJet_jetPt[centPos][tempEtaPos_gen][ptPos_genPho],aktR_truth_jet_pt_p->at(tjI), fullWeight);
                    fillTH1(h1F_genPho_genJet_jetPt[centPos][tempEtaPos_gen][nGammaPtBinsSub],aktR_truth_jet_pt_p->at(tjI), fullWeight);
                }

            }//truth jet loop
        }// isGoodTruthPhoton
    }//isMC 

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    double leadingPhoPt = 0;
    int leadingPhoIndex = -1; 
    double leadingPhoPt_genMatchedReco = 0;
    int leadingPhoIndex_genMatchedReco = -1; 
    bool isGoodGenMatchedRecoPhoton = true;
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

      if(isMC){
          if(truthPhotonPt<=0) continue; //prompt photons
          if(truthPhotonIso>genIsoCut) continue; // truth isolation condition
          if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) > phoGenMatchingDR) continue;

          if(leadingPhoPt_genMatchedReco < photon_pt_p->at(pI)){ 
              leadingPhoPt_genMatchedReco = photon_pt_p->at(pI);
              leadingPhoIndex_genMatchedReco = pI; 
          }
      }
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    if(leadingPhoIndex == -1) continue;
    if(isMC && leadingPhoIndex_genMatchedReco == -1) isGoodGenMatchedRecoPhoton = false;

    Int_t ptPos = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photon_pt_p->at(leadingPhoIndex), true, doGlobalDebug);
    int tempEtaPos = -1;
    Float_t etaValMain = TMath::Abs(photon_eta_p->at(leadingPhoIndex));
    for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
        if(etaValMain>=etaBins_i[eI] && etaValMain<etaBins_f[eI]) tempEtaPos=eI;
    }
    if(tempEtaPos==-1) continue; //eta cut

    Int_t ptPos_genMatchedReco = -1;
    int tempEtaPos_genMatchedReco = -1;
    if(isMC && isGoodGenMatchedRecoPhoton){
        ptPos_genMatchedReco = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photon_pt_p->at(leadingPhoIndex_genMatchedReco), true, doGlobalDebug);
        Float_t etaValMain_genMatchedReco = TMath::Abs(photon_eta_p->at(leadingPhoIndex_genMatchedReco));
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            if(etaValMain_genMatchedReco>=etaBins_i[eI] && etaValMain_genMatchedReco<etaBins_f[eI]) tempEtaPos_genMatchedReco=eI;
        }
    }
    if(isMC && tempEtaPos_genMatchedReco == -1) isGoodGenMatchedRecoPhoton = false;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    /////////////////// photon isolation correction
    float correctedIso = photon_etcone_p->at(leadingPhoIndex);
    if(doPtCorrectedIso && doCentCorrectedIso)
        correctedIso = getCorrectedPhotonIsolation(isPP, photon_etcone_p->at(leadingPhoIndex), photon_pt_p->at(leadingPhoIndex), photon_eta_p->at(leadingPhoIndex), cent);
    else if(doPtCorrectedIso && !doCentCorrectedIso)
        correctedIso = getPtCorrectedPhotonIsolation(photon_etcone_p->at(leadingPhoIndex), photon_pt_p->at(leadingPhoIndex), photon_eta_p->at(leadingPhoIndex));

    if(!doBkgPhoton){
        if(!(photon_tight_p->at(leadingPhoIndex)==1 && correctedIso < isoCut)) continue;
    } else{
        if(photonSelection==1){
            if(!(photon_tight_p->at(leadingPhoIndex)==1 && correctedIso > isoCut+bkgIsoGap)) continue;
        } else if(photonSelection==2){
            if(!(photon_tight_p->at(leadingPhoIndex)==0 && correctedIso < isoCut)) continue;
        } else if(photonSelection==3){
            if(!(photon_tight_p->at(leadingPhoIndex)==0 && correctedIso > isoCut+bkgIsoGap)) continue;
        } else if(photonSelection==4){
            if(!(correctedIso > isoCut+bkgIsoGap)) continue;
        } else if(photonSelection==5){
            if(!(photon_tight_p->at(leadingPhoIndex)==0)) continue;
        }
    }

    if(isMC){ 
        fillTH1(h1F_nPhoton[centPos][tempEtaPos],photon_pt_p->at(leadingPhoIndex),fullWeight);
        if(isGoodGenMatchedRecoPhoton) fillTH1(h1F_nPhoton_genMatchedReco[centPos][tempEtaPos_genMatchedReco],photon_pt_p->at(leadingPhoIndex_genMatchedReco),fullWeight);
    } else { h1F_nPhoton[centPos][tempEtaPos]->Fill(photon_pt_p->at(leadingPhoIndex)); }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    /////////////////////////////////////////////////////////////////////
    // JET LOOP 
    int multCounter = 0;
    for(unsigned int jI = 0; jI < aktRhi_em_xcalib_jet_pt_p->size(); ++jI){
        if(aktRhi_em_xcalib_jet_eta_p->at(jI) <= jtEtaBinsLow) continue;
        if(aktRhi_em_xcalib_jet_eta_p->at(jI) >= jtEtaBinsHigh) continue;

        Float_t dR = getDR(aktRhi_em_xcalib_jet_eta_p->at(jI), aktRhi_em_xcalib_jet_phi_p->at(jI), photon_eta_p->at(leadingPhoIndex), photon_phi_p->at(leadingPhoIndex));
        if(dR < gammaExclusionDR) continue;

        if(recoJtPtMin > aktRhi_em_xcalib_jet_pt_p->at(jI)) recoJtPtMin = aktRhi_em_xcalib_jet_pt_p->at(jI);

        if(aktRhi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
        if(aktRhi_em_xcalib_jet_pt_p->at(jI) >= jtPtBinsHigh) continue;

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
        Float_t dPhi = TMath::Abs(getDPHI(aktRhi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(leadingPhoIndex)));

        bool isGoodGenMatchedRecoJet = true;

        if(isMC){ 
            fillTH1(h1F_dphi_raw[centPos][tempEtaPos][ptPos], dPhi, fullWeight);
            if(dPhi >= gammaJtDPhiCut){
                fillTH1(h1F_jetPt_raw[centPos][tempEtaPos][ptPos], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
                ++multCounter;
            }

            if(isGoodGenMatchedRecoPhoton){
                if(aktRhi_truthpos_p->at(jI) < 0) isGoodGenMatchedRecoJet = false;
                if(isGoodGenMatchedRecoJet){
                    if(aktR_truth_jet_pt_p->at(aktRhi_truthpos_p->at(jI)) < 10) continue; 
                    Float_t dPhi_genMatchedReco = TMath::Abs(getDPHI(aktRhi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(leadingPhoIndex_genMatchedReco)));
                    fillTH1(h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[centPos][tempEtaPos_genMatchedReco][ptPos_genMatchedReco], dPhi_genMatchedReco, fullWeight);
                    fillTH1(h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[centPos][tempEtaPos_genMatchedReco][nGammaPtBinsSub], dPhi_genMatchedReco, fullWeight);

                    if(dPhi_genMatchedReco >= gammaJtDPhiCut){
                        fillTH1(h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[centPos][tempEtaPos_genMatchedReco][ptPos_genMatchedReco], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
                        fillTH1(h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[centPos][tempEtaPos_genMatchedReco][nGammaPtBinsSub], aktRhi_em_xcalib_jet_pt_p->at(jI), fullWeight);
                    }
                }
            } 
        } else { //data
            h1F_dphi_raw[centPos][tempEtaPos][ptPos]->Fill(dPhi);
            if(dPhi >= gammaJtDPhiCut){
                h1F_jetPt_raw[centPos][tempEtaPos][ptPos]->Fill(aktRhi_em_xcalib_jet_pt_p->at(jI));
                ++multCounter;
            }
        }
    } // END OF JET LOOP

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    if(isMC)
        fillTH1(h1F_nJetPerPhoton[centPos][tempEtaPos][ptPos], multCounter, fullWeight);
    else 
        h1F_nJetPerPhoton[centPos][tempEtaPos][ptPos]->Fill(multCounter);

    ///////////////////////////////////////////////////////////////////
    // Minbias JET LOOP for event mixing / in each photon loop!
	if(doMix){
	  unsigned long long mixCentPos = 0;
	  unsigned long long mixPsi2Pos = 0;
	  if(!isPP){
	    if(doMixCent) mixCentPos = ghostPos(nMixCentBins, mixCentBins, cent);
	    if(doMixPsi2){
	      if(evtPlane2Phi > TMath::Pi()/2) evtPlane2Phi -= TMath::Pi();
	      else if(evtPlane2Phi < -TMath::Pi()/2) evtPlane2Phi += TMath::Pi();
	      mixPsi2Pos = ghostPos(nMixPsi2Bins, mixPsi2Bins, evtPlane2Phi);
	    }
	  }

	  unsigned long long mixVzPos = 0;
	  if(doMixVz) mixVzPos = ghostPos(nMixVzBins, mixVzBins, vert_z);
	  
	  std::vector<unsigned long long> eventKeyVect;
	  if(doMixCent) eventKeyVect.push_back(mixCentPos);
	  if(doMixPsi2) eventKeyVect.push_back(mixPsi2Pos);
	  if(doMixVz) eventKeyVect.push_back(mixVzPos);
	  
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

	  unsigned long long key = keyBoy.GetKey(eventKeyVect);
	  unsigned long long maxPos = mixingMap[key].size();
	  if(maxPos == 0){
	    std::cout << "WHOOPS NO AVAILABLE MIXED EVENT. bailing" << std::endl;
	    std::cout << key << ", " << mixCentPos << ", " << cent << std::endl;
	    return 1;
	  }

	  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 


      //std::cout << "key, Number of mixing events MAX: " << key << ", " << maxPos << std::endl; 
	  //std::cout << "CENT: " << mixCentPos << ", " << cent << std::endl;
	  //std::cout << "PSI2: " << mixPsi2Pos << ", " << evtPlane2Phi << std::endl;
	  //std::cout << "VZ: " << mixVzPos << ", " << vert_z << std::endl;
      for(unsigned long long jetPos=0; jetPos < nMix; ++jetPos){
          std::vector<ROOT::Math::PtEtaPhiMVector> const& jets = mixingMap[key][jetPos];
          if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE, JETS SIZE, MAX, CHOSEN: " << __FILE__ << ", " << __LINE__ << ", " << jets.size() << ", " << maxPos << ", " << jetPos << std::endl; 

          for(unsigned int jI = 0; jI < jets.size(); ++jI){

              if(jets[jI].Pt()< jtPtBinsLow) continue;
              if(jets[jI].Eta() <= jtEtaBinsLow) continue;
              if(jets[jI].Eta() >= jtEtaBinsHigh) continue;

              Float_t dR = getDR(jets[jI].Eta(), jets[jI].Phi(), photon_eta_p->at(leadingPhoIndex), photon_phi_p->at(leadingPhoIndex));
              if(dR < gammaExclusionDR) continue;

              Float_t dPhi = TMath::Abs(getDPHI(jets[jI].Phi(), photon_phi_p->at(leadingPhoIndex)));
              if(isMC){
                  fillTH1(h1F_dphi_raw_mix[centPos][tempEtaPos][ptPos], dPhi, fullWeight);
                  if(dPhi >= gammaJtDPhiCut)
                      fillTH1(h1F_jetPt_raw_mix[centPos][tempEtaPos][ptPos],  jets[jI].Pt(), fullWeight);
              } else{
                  h1F_dphi_raw_mix[centPos][tempEtaPos][ptPos]->Fill(dPhi);
                  if(dPhi >= gammaJtDPhiCut)
                      h1F_jetPt_raw_mix[centPos][tempEtaPos][ptPos]->Fill(jets[jI].Pt());
              }
          }// END OF MIXING JET LOOP EACH EVENTS (MINBIAS)
          fillTH1(h1F_nMix[centPos][tempEtaPos],photon_pt_p->at(leadingPhoIndex),fullWeight);
	  }// the number of mixing eventsLOOP 

    } // doMix
  } // END OF EVENT LOOP
  

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  outFile_p->cd();


  ///////////////////////////////////////////////////////////
  // Write histograms in the output file
  for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
          for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
              h1F_jetPt_raw[cI][eI][pI]->Write("", TObject::kOverwrite);
              h1F_dphi_raw[cI][eI][pI]->Write("", TObject::kOverwrite);
              if(doMix){
                  h1F_jetPt_raw_mix[cI][eI][pI]->Write("", TObject::kOverwrite);
                  h1F_dphi_raw_mix[cI][eI][pI]->Write("", TObject::kOverwrite);
              }
              h1F_nJetPerPhoton[cI][eI][pI]->Write("", TObject::kOverwrite);
                
              if(isMC){
                  h1F_genPho_genJet_jetPt[cI][eI][pI]->Write("", TObject::kOverwrite);
                  h1F_genPho_genJet_dphi[cI][eI][pI]->Write("", TObject::kOverwrite);
                  h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[cI][eI][pI]->Write("", TObject::kOverwrite);
                  h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[cI][eI][pI]->Write("", TObject::kOverwrite);
              }
          }//pt
          if(doMix) h1F_nMix[cI][eI]->Write("", TObject::kOverwrite);
          h1F_nPhoton[cI][eI]->Write("", TObject::kOverwrite);
          if(isMC){
              h1F_nPhoton_gen[cI][eI]->Write("", TObject::kOverwrite);
              h1F_nPhoton_genMatchedReco[cI][eI]->Write("", TObject::kOverwrite);
          }
      }
  }

  runNumber_p->Write("", TObject::kOverwrite);

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
      for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
          for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
              delete h1F_jetPt_raw[cI][eI][pI];
              delete h1F_dphi_raw[cI][eI][pI];
              if(doMix){
                  delete h1F_jetPt_raw_mix[cI][eI][pI];
                  delete h1F_dphi_raw_mix[cI][eI][pI];
              }
              delete h1F_nJetPerPhoton[cI][eI][pI];
              if(isMC){
                 delete h1F_genPho_genJet_jetPt[cI][eI][pI];
                 delete h1F_genPho_genJet_dphi[cI][eI][pI];
                 delete h1F_genMatchedRecoPho_genMatchedRecoJet_jetPt[cI][eI][pI];
                 delete h1F_genMatchedRecoPho_genMatchedRecoJet_dphi[cI][eI][pI];
              }
          }//pt
         if(doMix) delete h1F_nMix[cI][eI];
         delete h1F_nPhoton[cI][eI];
          if(isMC){
              delete h1F_nPhoton_gen[cI][eI];
              delete h1F_nPhoton_genMatchedReco[cI][eI];
          }
      }
  }
  
  delete runNumber_p;

  if(isMC){
    delete pthat_p;
    delete pthat_Unweighted_p;
  }

  if(!isPP){
    delete centrality_p;
    if(isMC) delete centrality_Unweighted_p;
  }

  config_p->Write("config", TObject::kOverwrite);

  TEnv labelEnv;
  for(auto const & lab : binsToLabelStr){
    labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
  }
  labelEnv.Write("label", TObject::kOverwrite);

  outFile_p->Close();
  delete outFile_p;

  delete randGen_p;
  
  std::cout << "phoTaggedJetRaa_jetPt COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/phoTaggedJetRaa_jetPt.exe <inConfigFileName> <doBKGphoton>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += phoTaggedJetRaa_jetPt(argv[1], atoi(argv[2]));
  return retVal;
}
