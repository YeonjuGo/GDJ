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
#include "include/configParser.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/keyHandler.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"

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

int phoTaggedJetRaa(std::string inConfigFileName)
{
    const Int_t randSeed = 5573; // from coin flips -> binary number 1010111000101
    TRandom3* randGen_p = new TRandom3(randSeed);

    checkMakeDir check;
    if(!check.checkFileExt(inConfigFileName, ".txt")) return 1;

    globalDebugHandler gDebug;
    const bool doGlobalDebug = gDebug.GetDoGlobalDebug();

    configParser config(inConfigFileName);

    std::vector<std::string> necessaryParams = {"INFILENAME",
        "OUTFILENAME",
        "CENTFILENAME",
        "MIXFILENAME",
        "DOMIX",
        "ISPP",
        "ISMC",
        "CENTBINS",
        "NGAMMAPTBINS",
        "GAMMAPTBINSLOW",
        "GAMMAPTBINSHIGH",
        "GAMMAPTBINSDOLOG",
        "NETABINS",
        "ETABINSLOW",
        "ETABINSHIGH",
        "ETABINSDOABS",
        "NPHIBINS",
        "NGAMMAPTBINSSUB",
        "GAMMAPTBINSSUBLOW",
        "GAMMAPTBINSSUBHIGH",
        "GAMMAPTBINSDOLOG",
        "NETABINSSUB",
        "ETABINSSUBLOW",
        "ETABINSSUBHIGH",
        "ETABINSSUBDOABS",
        "NJTPTBINS",
        "JTPTBINSLOW",
        "JTPTBINSHIGH",
        "JTPTBINSDOLOG",
        "NDPHIBINS",
        "GAMMAJTDPHI",  
        "DOGAMMAJTDPHICUT",  
        "NXJBINS",
        "XJBINSLOW",
        "XJBINSHIGH",
        "NMASSBINS",
        "MASSBINSLOW",
        "MASSBINSHIGH",
        "DOJETMAXPTCUT",
        "DRPHOTONJET",
        "DRJETCONE",
        "DOQUARKJET",
        "DOGLUONJET",
        "DOBKGPHOTONS",
        "DOTHISPHOPTBINS",
        "DOTHISJETPTBINS",
        "PHOPTBINS",
        "JETPTBINS",
        "ISOCUT"};

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


    if(!config.ContainsParamSet(necessaryParams)) return 1;  

    std::string inROOTFileName = config.GetConfigVal("INFILENAME");
    std::string inCentFileName = config.GetConfigVal("CENTFILENAME");
    std::string outFileName = config.GetConfigVal("OUTFILENAME");
    std::string inGRLFileName = config.GetConfigVal("GRLFILENAME");
    std::string inMixFileName = config.GetConfigVal("MIXFILENAME");
    std::string jetDR = config.GetConfigVal("DRJETCONE");
    const bool doMix = std::stoi(config.GetConfigVal("DOMIX"));
    const bool doJetMaxPtCut = std::stoi(config.GetConfigVal("DOJETMAXPTCUT"));
    const bool doBackgroundPhotons = std::stoi(config.GetConfigVal("DOBKGPHOTONS"));
    const Float_t DRPhotonJetCut = std::stof(config.GetConfigVal("DRPHOTONJET"));
    const Float_t doDPhiCut = std::stof(config.GetConfigVal("DOGAMMAJTDPHICUT"));
    const Float_t isoCut = std::stof(config.GetConfigVal("ISOCUT"));

    //const bool doThisPhoPtBins= std::stoi(config.GetConfigVal("DOTHISPHOPTBINS"));
    const bool doThisJetPtBins= std::stoi(config.GetConfigVal("DOTHISPHOPTBINS"));
    if(doMix && !config.ContainsParamSet(mixParams)) return 1;

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

        doMixCent = (bool)std::stoi(config.GetConfigVal("DOMIXCENT"));

        if(doMixCent){      
            nMixCentBins = std::stoi(config.GetConfigVal("NMIXCENTBINS"));
            if(nMixCentBins > nMaxMixBins){
                std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixCentBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
                return 1;
            }
            mixCentBinsLow = std::stof(config.GetConfigVal("MIXCENTBINSLOW"));
            mixCentBinsHigh = std::stof(config.GetConfigVal("MIXCENTBINSHIGH"));    
            getLinBins(mixCentBinsLow, mixCentBinsHigh, nMixCentBins, mixCentBins);

            mixVect.push_back({});
            for(Int_t mI = 0; mI < nMixCentBins; ++mI){
                mixVect[mixVect.size()-1].push_back(mI);
                keyVect.push_back({(unsigned long long)mI});
            }
        }

        doMixPsi2 = (bool)std::stoi(config.GetConfigVal("DOMIXPSI2"));

        if(doMixPsi2){
            nMixPsi2Bins = std::stoi(config.GetConfigVal("NMIXPSI2BINS"));
            if(nMixPsi2Bins > nMaxMixBins){
                std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixPsi2Bins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
                return 1;
            }
            mixPsi2BinsLow = std::stof(config.GetConfigVal("MIXPSI2BINSLOW"));
            mixPsi2BinsHigh = std::stof(config.GetConfigVal("MIXPSI2BINSHIGH"));    
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

        doMixVz = (bool)std::stoi(config.GetConfigVal("DOMIXVZ"));

        if(doMixVz){
            nMixVzBins = std::stoi(config.GetConfigVal("NMIXVZBINS"));
            if(nMixVzBins > nMaxMixBins){
                std::cout << "GDJNTUPLETOHIST ERROR - centrality mixing bins \'" << nMixVzBins << "\' exceeds max \'" << nMaxMixBins << "\'. return 1" << std::endl;
                return 1;
            }
            mixVzBinsLow = std::stof(config.GetConfigVal("MIXVZBINSLOW"));
            mixVzBinsHigh = std::stof(config.GetConfigVal("MIXVZBINSHIGH"));    
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

    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    keyHandler runLumiKey("runLumiHandler");//for lumi estimation
    runLumiKey.Init({1000000, 10000});//runnumbers, then lumi

    keyHandler keyBoy("mixingHandler");//For Mixing
    std::map<unsigned long long, std::vector<std::vector<TLorentzVector> > > mixingMap;
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
        }
    }  

    const bool isMC = std::stoi(config.GetConfigVal("ISMC"));
   // bool doQuarkJet = false; 
   // bool doGluonJet = false; 
   // if(isMC){
   //     doQuarkJet = std::stoi(config.GetConfigVal("DOQUARKJET"));
   //     doGluonJet = std::stoi(config.GetConfigVal("DOGLUONJET"));
   // }

    if(!check.checkFileExt(inROOTFileName, "root")) return 1; // Check input is valid ROOT file
    if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
    if(!isMC && !check.checkFileExt(inGRLFileName, "xml")) return 1; // GRL File xml
    if(doMix && !check.checkFileExt(inMixFileName, "root")) return 1; // Check mixed file exists if mixing is requested

    const std::string dateStr = getDateStr();
    check.doCheckMakeDir("output"); // check output dir exists; if not create
    check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

    centralityFromInput centTable(inCentFileName);
    if(doGlobalDebug) centTable.PrintTableTex();

    if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
    outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    const bool isPP = std::stoi(config.GetConfigVal("ISPP"));
    const Int_t nMaxSubBins = 10;
    const Int_t nMaxCentBins = 10;
    Int_t nCentBins = 1;

    std::string systStr = "PP";
    if(!isPP) systStr = "PbPb";

    std::vector<int> centBins;
    std::vector<std::string> centBinsStr = {systStr};
    std::map<std::string, std::string> binsToLabelStr;
    if(!isPP){
        centBins = strToVectI(config.GetConfigVal("CENTBINS"));
        nCentBins = centBins.size()-1;

        if(!goodBinning(inConfigFileName, nMaxCentBins, nCentBins, "CENTBINS")) return 1;

        centBinsStr.clear();
        for(Int_t cI = 0; cI < nCentBins; ++cI){
            centBinsStr.push_back("Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]));

            binsToLabelStr[centBinsStr[cI]] = std::to_string(centBins[cI]) + "-" + std::to_string(centBins[cI+1]) + "%";
        }
    }
    else binsToLabelStr[centBinsStr[0]] = "pp";

    const Int_t nMaxPtBins = 300;
    const Int_t nGammaPtBins = std::stoi(config.GetConfigVal("NGAMMAPTBINS"));
    if(!goodBinning(inConfigFileName, nMaxPtBins, nGammaPtBins, "NGAMMAPTBINS")) return 1;
    const Float_t gammaPtBinsLow = std::stof(config.GetConfigVal("GAMMAPTBINSLOW"));
    const Float_t gammaPtBinsHigh = std::stof(config.GetConfigVal("GAMMAPTBINSHIGH"));
    const Bool_t gammaPtBinsDoLog = std::stof(config.GetConfigVal("GAMMAPTBINSDOLOG"));
    Double_t gammaPtBins[nMaxPtBins+1];
    if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
    else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
    std::vector<std::string> genGammaPtBinsStr, recoGammaPtBinsStr;
    for(Int_t pI = 0; pI < nGammaPtBins; ++pI){
        genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(pI));
        recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(pI));

        binsToLabelStr[genGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
        binsToLabelStr[recoGammaPtBinsStr[pI]] = prettyString(gammaPtBins[pI], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[pI+1], 1, false);
    }
    genGammaPtBinsStr.push_back("GenGammaPt" + std::to_string(nGammaPtBins));
    recoGammaPtBinsStr.push_back("RecoGammaPt" + std::to_string(nGammaPtBins));

    binsToLabelStr[genGammaPtBinsStr[genGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Gen. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);
    binsToLabelStr[recoGammaPtBinsStr[recoGammaPtBinsStr.size()-1]] = prettyString(gammaPtBins[0], 1, false) + " < Reco. p_{T,#gamma} < " + prettyString(gammaPtBins[nGammaPtBins], 1, false);

    const Int_t nMaxEtaPhiBins = 100;
    const Int_t nEtaBins = std::stoi(config.GetConfigVal("NETABINS"));
    if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nEtaBins, "NETABINS")) return 1;
    const Float_t etaBinsLow = std::stof(config.GetConfigVal("ETABINSLOW"));
    const Float_t etaBinsHigh = std::stof(config.GetConfigVal("ETABINSHIGH"));
    const Bool_t etaBinsDoAbs = std::stoi(config.GetConfigVal("ETABINSDOABS"));
    if(etaBinsDoAbs && etaBinsLow < 0){
        std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsLow \'" << etaBinsLow << "\' less than 0 despite requested etaBinsDoAbs \'" << etaBinsDoAbs << "\'. return 1" << std::endl;
        return 1;
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    Double_t etaBins[nMaxEtaPhiBins+1];
    getLinBins(etaBinsLow, etaBinsHigh, nEtaBins, etaBins);

    const Int_t nPhiBins = std::stoi(config.GetConfigVal("NPHIBINS"));
    if(!goodBinning(inConfigFileName, nMaxEtaPhiBins, nPhiBins, "NPHIBINS")) return 1;
    Double_t phiBins[nMaxEtaPhiBins+1];
    getLinBins(-TMath::Pi()+0.01, TMath::Pi()+0.01, nPhiBins, phiBins);

    //Pt sub bins handling
    //std::vector<int> phoPtBins;
    Int_t nGammaPtBinsSub = std::stoi(config.GetConfigVal("NGAMMAPTBINSSUB"));
    //if(doThisPhoPtBins){
    //    gammaPtBinsSub = strToVectI(config.GetConfigVal("PHOPTBINS"));
    //    nGammaPtBinsSub = gammaPtBinsSub.size()-1;
    //} else{
    nGammaPtBinsSub = std::stoi(config.GetConfigVal("NGAMMAPTBINSSUB"));
    if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaPtBinsSub, "NGAMMAPTBINSSUB")) return 1;
    const Float_t gammaPtBinsSubLow = std::stof(config.GetConfigVal("GAMMAPTBINSSUBLOW"));
    const Float_t gammaPtBinsSubHigh = std::stof(config.GetConfigVal("GAMMAPTBINSSUBHIGH"));
    if(gammaPtBinsSubLow < gammaPtBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubLow \'" << gammaPtBinsSubLow << "\' less than gammaPtBinsLow \'" << gammaPtBinsLow << "\'. return 1" << std::endl;
    if(gammaPtBinsSubHigh > gammaPtBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubHigh \'" << gammaPtBinsSubHigh << "\' greater than gammaPtBinsHigh \'" << gammaPtBinsHigh << "\'. return 1" << std::endl;
    if(gammaPtBinsSubLow < gammaPtBinsLow || gammaPtBinsSubHigh > gammaPtBinsHigh) return 1;
    const Bool_t gammaPtBinsSubDoLog = std::stof(config.GetConfigVal("GAMMAPTBINSDOLOG"));
    Double_t gammaPtBinsSub[nMaxSubBins+1];
    //Double_t gammaPtBinsSub[nMaxSubBins+1];
    if(gammaPtBinsSubDoLog) getLogBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
    else getLinBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
    //}

    std::vector<std::string> gammaPtBinsSubStr;
    for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
        gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(pI));

        binsToLabelStr[gammaPtBinsSubStr[pI]] = prettyString(gammaPtBinsSub[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[pI+1], 1, false);
    }
    gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(nGammaPtBinsSub));

    binsToLabelStr[gammaPtBinsSubStr[gammaPtBinsSubStr.size()-1]] = prettyString(gammaPtBinsSub[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false);

    //Eta sub bins handling
    const Int_t nEtaBinsSub = std::stoi(config.GetConfigVal("NETABINSSUB"));
    if(!goodBinning(inConfigFileName, nMaxSubBins, nEtaBinsSub, "NETABINSSUB")) return 1;
    const Float_t etaBinsSubLow = std::stof(config.GetConfigVal("ETABINSSUBLOW"));
    const Float_t etaBinsSubHigh = std::stof(config.GetConfigVal("ETABINSSUBHIGH"));
    if(etaBinsSubLow < etaBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than etaBinsLow \'" << etaBinsLow << "\'. return 1" << std::endl;
    if(etaBinsSubHigh > etaBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubHigh \'" << etaBinsSubHigh << "\' greater than etaBinsHigh \'" << etaBinsHigh << "\'. return 1" << std::endl;
    if(etaBinsSubLow < etaBinsLow || etaBinsSubHigh > etaBinsHigh) return 1;
    const	Bool_t etaBinsSubDoAbs = std::stoi(config.GetConfigVal("ETABINSSUBDOABS"));
    if(etaBinsSubDoAbs && etaBinsSubLow < 0){
        std::cout << "ERROR - config \'" << inConfigFileName << "\' contains etaBinsSubLow \'" << etaBinsSubLow << "\' less than 0 despite requested etaBinsSubDoAbs \'" << etaBinsSubDoAbs << "\'. return 1" << std::endl;
        return 1;
    }
    Double_t etaBinsSub[nMaxSubBins+1];
    getLinBins(etaBinsSubLow, etaBinsSubHigh, nEtaBinsSub, etaBinsSub);
    std::vector<std::string> etaBinsSubStr;
    std::string preStr = "";
    if(etaBinsSubDoAbs) preStr = "Abs";
    for(Int_t eI = 0; eI < nEtaBinsSub; ++eI){
        etaBinsSubStr.push_back(preStr + "Eta" + std::to_string(eI));

        if(etaBinsSubDoAbs) binsToLabelStr[etaBinsSubStr[eI]] = prettyString(etaBinsSub[eI], 2, false) + "<|#eta_{#gamma}|<" + prettyString(etaBinsSub[eI+1], 2, false);
        else binsToLabelStr[etaBinsSubStr[eI]] = prettyString(etaBinsSub[eI], 2, false) + " < #eta_{#gamma} < " + prettyString(etaBinsSub[eI+1], 2, false);
    }
    etaBinsSubStr.push_back(preStr + "Eta" + std::to_string(nEtaBinsSub));
    if(etaBinsSubDoAbs) binsToLabelStr[etaBinsSubStr[etaBinsSubStr.size()-1]] = prettyString(etaBinsSub[0], 2, false) + "<|#eta_{#gamma}|<" + prettyString(etaBinsSub[nEtaBinsSub], 2, false);
    else binsToLabelStr[etaBinsSubStr[etaBinsSubStr.size()-1]] = prettyString(etaBinsSub[0], 2, false) + " < #eta_{#gamma} < " + prettyString(etaBinsSub[nEtaBinsSub], 2, false);


    //jet pt binning
    //const Int_t nMaxCentBins = 10;
    //std::vector<std::string> centBinsStr = {systStr};

        //if(!goodBinning(inConfigFileName, nMaxCentBins, nCentBins, "CENTBINS")) return 1;

  //      centBinsStr.clear();
  //      for(Int_t cI = 0; cI < nCentBins; ++cI){
  //          centBinsStr.push_back("Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]));

  //          binsToLabelStr[centBinsStr[cI]] = std::to_string(centBins[cI]) + "-" + std::to_string(centBins[cI+1]) + "%";
  //      }
  //  else binsToLabelStr[centBinsStr[0]] = "pp";

    Int_t nJtPtBins;
    std::vector<int> jetBins_vec;
    Double_t jtPtBins[nMaxPtBins+1];
    const Float_t jtPtBinsLow = std::stof(config.GetConfigVal("JTPTBINSLOW"));
    const Float_t jtPtBinsHigh = std::stof(config.GetConfigVal("JTPTBINSHIGH"));
    const Bool_t jtPtBinsDoLog = std::stof(config.GetConfigVal("JTPTBINSDOLOG"));

    if(doThisJetPtBins){
        jetBins_vec = strToVectI(config.GetConfigVal("JETPTBINS"));
        nJtPtBins = jetBins_vec.size()-1;
        for(int pI=0; pI < nJtPtBins; ++pI){
            jtPtBins[pI] = jetBins_vec[pI]; 
        }
    } else {
        nJtPtBins = std::stoi(config.GetConfigVal("NJTPTBINS"));
        if(!goodBinning(inConfigFileName, nMaxPtBins, nJtPtBins, "NJTPTBINS")) return 1;
        //Double_t jtPtBins[nMaxPtBins+1];
        if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
        else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
    }



    std::string jtPtBinsGlobalStr = "GlobalJtPt0";
    std::string jtPtBinsGlobalLabel = "";
    if(doThisJetPtBins) jtPtBinsGlobalLabel = "p_{T,jet} > " + std::to_string(jetBins_vec[0]); 
    else jtPtBinsGlobalLabel = "p_{T,jet} > " + prettyString(jtPtBinsLow,1,false);
    binsToLabelStr[jtPtBinsGlobalStr] = jtPtBinsGlobalLabel;   

    std::string multiJtCutGlobalStr = "MultiJt0";
    std::string multiJtCutGlobalLabel = "N_{Jet,Reco.} >= 2";
    binsToLabelStr[multiJtCutGlobalStr] = multiJtCutGlobalLabel;   

    std::vector<std::string> jtPtBinsStr;
    for(Int_t pI = 0; pI < nJtPtBins; ++pI){
        if(pI < 10) jtPtBinsStr.push_back("JtPt0" + std::to_string(pI));
        else jtPtBinsStr.push_back("JtPt" + std::to_string(pI));

        if(doThisJetPtBins) binsToLabelStr[jtPtBinsStr[pI]] = " p_{T,Jet} > " + std::to_string(jetBins_vec[pI]); 
        else binsToLabelStr[jtPtBinsStr[pI]] = " p_{T,Jet} > " + prettyString(jtPtBins[pI], 1, false);
        //binsToLabelStr[jtPtBinsStr[pI]] = " p_{T,Jet} > " + prettyString(jtPtBins[pI], 1, false);
        //binsToLabelStr[jtPtBinsStr[pI]] = prettyString(jtPtBins[pI], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[pI+1], 1, false);
    }
    jtPtBinsStr.push_back("JtPt" + std::to_string(nJtPtBins));
    binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] =" p_{T,Jet} > " + prettyString(jtPtBins[0], 1, false);
    //binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] = prettyString(jtPtBins[0], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[nJtPtBins], 1, false);


    //const Int_t nDPhiBins = std::stoi(config.GetConfigVal("NDPHIBINS"));
    const Double_t gammaJtDPhiCut = mathStringToNum(config.GetConfigVal("GAMMAJTDPHI"));  
    const std::string gammaJtDPhiStr = "DPhi0";
    std::string gammaJtDPhiLabel = returnAllCapsString(config.GetConfigVal("GAMMAJTDPHI"));
    if(gammaJtDPhiLabel.find("PI") != std::string::npos) gammaJtDPhiLabel.replace(gammaJtDPhiLabel.find("PI"), 2, "#pi");
    gammaJtDPhiLabel = "|#Delta#phi_{#gamma,jet}| > " + gammaJtDPhiLabel;
    //binsToLabelStr[gammaJtDPhiStr] = gammaJtDPhiLabel;

    const Int_t nMassBins = std::stoi(config.GetConfigVal("NMASSBINS"));
    const Float_t massBinsLow = std::stof(config.GetConfigVal("MASSBINSLOW"));
    const Float_t massBinsHigh = std::stof(config.GetConfigVal("MASSBINSHIGH"));
    Double_t massBins[nMaxPtBins+1];
    getLinBins(massBinsLow, massBinsHigh, nMassBins, massBins);

    const Int_t nXJBins = std::stoi(config.GetConfigVal("NXJBINS"));
    const Float_t xjBinsLow = std::stof(config.GetConfigVal("XJBINSLOW"));
    const Float_t xjBinsHigh = std::stof(config.GetConfigVal("XJBINSHIGH"));
    Double_t xjBins[nMaxPtBins+1];
    getLinBins(xjBinsLow, xjBinsHigh, nXJBins, xjBins);

    ////////////////////////////////////////////////////////////
    // define histograms  

    TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
    TH1F* runNumber_p = nullptr;
    TH1F* lumiFractionPerRun_p = nullptr;
    TH1F* pthat_p = nullptr;
    TH1F* pthat_Unweighted_p = nullptr;
    TH1F* centrality_p = nullptr;
    TH1F* centrality_Unweighted_p = nullptr;

    //DATA
    TH2F* PhotonTaggedJetPtVsPhotonPt_CentDep[nMaxCentBins];
    TH1F* PhotonTaggedJetPt_CentPhoPtDep[nMaxCentBins][nMaxSubBins+1];
    TH1F* PhotonTaggedJetPtPerPho_CentPhoPtDep[nMaxCentBins][nMaxSubBins+1];

    //MC
    TH2F* GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[nMaxCentBins]; // gen matching & use gen values
    TH2F* GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[nMaxCentBins];// gen matching
    TH1F* GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[nMaxCentBins][nMaxSubBins+1];// gen matching & use gen values
    TH1F* GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[nMaxCentBins][nMaxSubBins+1];// gen matching
    TH1F* PhotonTaggedGenJetPt_CentPhoPtDep[nMaxCentBins][nMaxSubBins+1];// gen matching & use gen values
    TH1F* PhotonTaggedGenMatchedJetPt_CentPhoPtDep[nMaxCentBins][nMaxSubBins+1];// gen matching
    //TH1F* photonEff_CentDep_Den[nMaxCentBins];// gen matching
    //TH1F* photonEff_CentDep_Num[nMaxCentBins];// gen matching
    //TH1F* photonEff_CentDep_Eff[nMaxCentBins];// gen matching


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


    ///////////////////////////////////////////////////////////
    // set histograms 
    float nPtBins_2D = 220;
    float nPtmin_2D = 30;
    float nPtmax_2D = 250;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

        PhotonTaggedJetPtVsPhotonPt_CentDep[cI] = new TH2F(("PhotonTaggedJetPtVsPhotonPt_CentDep" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma p_{T} [GeV];Jet p_{T} [GeV]", nPtBins_2D, nPtmin_2D, nPtmax_2D, nPtBins_2D, nPtmin_2D, nPtmax_2D);

        centerTitles({PhotonTaggedJetPtVsPhotonPt_CentDep[cI]});
        setSumW2({PhotonTaggedJetPtVsPhotonPt_CentDep[cI]});

        if(isMC){
            GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[cI] = new TH2F(("GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma p_{T} [GeV];Jet p_{T} [GeV]", nPtBins_2D, nPtmin_2D, nPtmax_2D, nPtBins_2D, nPtmin_2D, nPtmax_2D);
            GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[cI] = new TH2F(("GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep" + centBinsStr[cI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma p_{T} [GeV];Jet p_{T} [GeV]", nPtBins_2D, nPtmin_2D, nPtmax_2D, nPtBins_2D, nPtmin_2D, nPtmax_2D);

            centerTitles({GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[cI],GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[cI]});
            setSumW2({GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[cI],GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[cI]});
            
        }

        for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
            PhotonTaggedJetPt_CentPhoPtDep[cI][pI] = new TH1F(("PhotonTaggedJetPt_CentPhoPtDep" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];N_{#gamma,jet}", nJtPtBins, jtPtBins);
            PhotonTaggedJetPtPerPho_CentPhoPtDep[cI][pI] = new TH1F(("PhotonTaggedJetPtPerPho_CentPhoPtDep" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{1}{N_{#gamma}}#frac{N_{#gamma,jet}}{dp_{T}^{jet}}", nJtPtBins, jtPtBins);

            centerTitles({PhotonTaggedJetPt_CentPhoPtDep[cI][pI],PhotonTaggedJetPtPerPho_CentPhoPtDep[cI][pI]});
            setSumW2({PhotonTaggedJetPt_CentPhoPtDep[cI][pI],PhotonTaggedJetPtPerPho_CentPhoPtDep[cI][pI]});

            if(isMC){
                GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[cI][pI] = new TH1F(("GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
                GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI] = new TH1F(("GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
                PhotonTaggedGenJetPt_CentPhoPtDep[cI][pI] = new TH1F(("PhotonTaggedGenJetPt_CentPhoPtDep" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);
                PhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI] = new TH1F(("PhotonTaggedGenMatchedJetPt_CentPhoPtDep" + centBinsStr[cI] + "_" + gammaPtBinsSubStr[pI] + "_" + gammaJtDPhiStr + "_h").c_str(), ";#gamma-tagged Jet p_{T} [GeV];#frac{N_{#gamma,jet}}{N_{#gamma}}", nJtPtBins, jtPtBins);

                centerTitles({GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[cI][pI],GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI]});
                setSumW2({GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[cI][pI],GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI]});

                centerTitles({PhotonTaggedGenJetPt_CentPhoPtDep[cI][pI],PhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI]});
                setSumW2({PhotonTaggedGenJetPt_CentPhoPtDep[cI][pI],PhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI]});
            }
        }//photon pt loop
    }//centrailty loop

    ///////////////////////////////////////////////////////////
    // import input trees 

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    TFile* inFile_p = new TFile(inROOTFileName.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("runNumber", 1);
    Int_t runMin = inTree_p->GetMinimum("runNumber");
    Int_t runMax = inTree_p->GetMaximum("runNumber");
    Int_t nRunBins = runMax - runMin;
    Float_t runMinF = ((Float_t)runMin) - 0.5;
    Float_t runMaxF = ((Float_t)runMax) + 0.5;

    outFile_p->cd();
    runNumber_p = new TH1F(("runNumber_" + systStr + "_h").c_str(), ";Run;Counts", nRunBins+1, runMinF, runMaxF);
    if(!isMC) lumiFractionPerRun_p = new TH1F(("lumiFractionPerRun_" + systStr + "_h").c_str(), ";Run;Fraction of Lumiblocks", nRunBins+1, runMinF, runMaxF);
    centerTitles(runNumber_p);
    if(!isMC) centerTitles(lumiFractionPerRun_p);

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

    inFile_p->cd();

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

    float treePartonPt[2];
    float treePartonEta[2];
    float treePartonPhi[2];
    int treePartonId[2];

    //treePartonPt:treePartonEta:treePartonPhi:treePartonId:truth_pt:truth_eta:truth_phi:truth_pdg:akt4_truth_jet_pt:akt4_truth_jet_eta:akt4_truth_jet_phi:truthPhotonPt

    Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;

    std::vector<float>* photon_pt_p=nullptr;
    std::vector<float>* photon_eta_p=nullptr;
    std::vector<float>* photon_phi_p=nullptr;
    std::vector<bool>* photon_tight_p=nullptr;  
    std::vector<bool>* photon_loose_p=nullptr;  
    std::vector<float>* photon_etcone30_p=nullptr;

    std::vector<float>* akthi_em_xcalib_jet_pt_p=nullptr;
    std::vector<float>* akthi_em_xcalib_jet_eta_p=nullptr;
    std::vector<float>* akthi_em_xcalib_jet_phi_p=nullptr;
    std::vector<int>* akthi_truthpos_p=nullptr;

    std::vector<float>* akt_truth_jet_pt_p=nullptr;
    std::vector<float>* akt_truth_jet_eta_p=nullptr;
    std::vector<float>* akt_truth_jet_phi_p=nullptr;

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

        inTree_p->SetBranchStatus("treePartonPt", 1);
        inTree_p->SetBranchStatus("treePartonEta", 1);
        inTree_p->SetBranchStatus("treePartonPhi", 1);
        inTree_p->SetBranchStatus("treePartonId", 1);

        inTree_p->SetBranchStatus("truthPhotonPt", 1);
        inTree_p->SetBranchStatus("truthPhotonEta", 1);
        inTree_p->SetBranchStatus("truthPhotonPhi", 1);
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
    inTree_p->SetBranchStatus("photon_loose", 1);
    inTree_p->SetBranchStatus("photon_etcone30", 1);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    inTree_p->SetBranchStatus(("akt"+jetDR+"hi_em_xcalib_jet_pt").c_str(), 1);
    inTree_p->SetBranchStatus(("akt"+jetDR+"hi_em_xcalib_jet_eta").c_str(), 1);
    inTree_p->SetBranchStatus(("akt"+jetDR+"hi_em_xcalib_jet_phi").c_str(), 1);

    if(isMC){
        inTree_p->SetBranchStatus(("akt"+jetDR+"hi_truthpos").c_str(), 1);

        inTree_p->SetBranchStatus(("akt"+jetDR+"_truth_jet_pt").c_str(), 1);
        inTree_p->SetBranchStatus(("akt"+jetDR+"_truth_jet_eta").c_str(), 1);
        inTree_p->SetBranchStatus(("akt"+jetDR+"_truth_jet_phi").c_str(), 1);
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

        inTree_p->SetBranchAddress("treePartonPt", &treePartonPt);
        inTree_p->SetBranchAddress("treePartonEta", &treePartonEta);
        inTree_p->SetBranchAddress("treePartonPhi", &treePartonPhi);
        inTree_p->SetBranchAddress("treePartonId", &treePartonId);

        inTree_p->SetBranchAddress("truthPhotonPt", &truthPhotonPt);
        inTree_p->SetBranchAddress("truthPhotonEta", &truthPhotonEta);
        inTree_p->SetBranchAddress("truthPhotonPhi", &truthPhotonPhi);
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
    inTree_p->SetBranchAddress("photon_loose", &photon_loose_p);
    inTree_p->SetBranchAddress("photon_etcone30", &photon_etcone30_p);

    inTree_p->SetBranchAddress(("akt"+jetDR+"hi_em_xcalib_jet_pt").c_str(), &akthi_em_xcalib_jet_pt_p);
    inTree_p->SetBranchAddress(("akt"+jetDR+"hi_em_xcalib_jet_eta").c_str(), &akthi_em_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress(("akt"+jetDR+"hi_em_xcalib_jet_phi").c_str(), &akthi_em_xcalib_jet_phi_p);

    if(isMC){
        inTree_p->SetBranchAddress(("akt"+jetDR+"hi_truthpos").c_str(), &akthi_truthpos_p);    

        inTree_p->SetBranchAddress(("akt"+jetDR+"_truth_jet_pt").c_str(), &akt_truth_jet_pt_p);
        inTree_p->SetBranchAddress(("akt"+jetDR+"_truth_jet_eta").c_str(), &akt_truth_jet_eta_p);
        inTree_p->SetBranchAddress(("akt"+jetDR+"_truth_jet_phi").c_str(), &akt_truth_jet_phi_p);
    }


    Double_t recoJtPtMin = 100000.;
    const ULong64_t nEntries = inTree_p->GetEntries();
    const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

    //variable to count the number of photons in a given centrality and photon pt bin 
    std::vector<std::vector<Double_t> > gammaCountsPerPtCent;
    std::vector<std::vector<Double_t> > gammaJetCountsPerPtCent;
    std::vector<std::vector<Double_t> > genMatchedGammaCountsPerPtCent;
    for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
        gammaCountsPerPtCent.push_back({});
        gammaJetCountsPerPtCent.push_back({});
        genMatchedGammaCountsPerPtCent.push_back({});

        for(Int_t cI = 0; cI < nCentBins; ++cI){
            gammaCountsPerPtCent[pI].push_back(0.0);
            gammaJetCountsPerPtCent[pI].push_back(0.0);
            genMatchedGammaCountsPerPtCent[pI].push_back(0.0);
        }
    }

    //variable to count the number of events in a given centrality bin
    std::vector<Double_t> eventCountsPerCent;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        eventCountsPerCent.push_back(0.0);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    std::cout << "Processing " << nEntries << " events..." << std::endl;

    bool didOneFireMiss = false;
    std::vector<int> skippedCent;
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

    //bool gluonJetFlag = 0; // 0 is for quark jet
    //unsigned int treeJetIndex = -1;

    ///////////////////////////////////////////////////////////
    // Event loop! 
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
        if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
        inTree_p->GetEntry(entry);

        double vert_z = vert_z_p->at(0);
        vert_z /= 1000.;
        if(vert_z <= -15. || vert_z >= 15.) continue;      
        //    if(vert_z <= vzMixBinsLow || vert_z >= vzMixBinsHigh) continue;

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
                std::cout << "gdjNTupleToHist Warning - Skipping centrality \'" << (Int_t)cent << "\' as given centrality binning is \'" << centBins[0] << "-" << centBins[centBins.size()-1] << "\'. if this is incorrect please fix." << std::endl;
                skippedCent.push_back((Int_t)cent);
            }

            continue;
        }

        if(!isMC) fullWeight = -1.0;

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

        unsigned long long tempKey = runLumiKey.GetKey({(unsigned long long)runNumber, (unsigned long long)lumiBlock});    
        runLumiIsFired[tempKey] = true;

        fillTH1(runNumber_p, runNumber, fullWeight);	
        if(!isPP){
            fillTH1(centrality_p, cent, fullWeight);
            if(isMC) centrality_Unweighted_p->Fill(cent);
        }

        if(isMC){
            fillTH1(pthat_p, pthat, fullWeight);
            pthat_Unweighted_p->Fill(pthat);
        }

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 

        bool isPhotonJetEvent = false;
        ///////////////////////////////////////////////////////////
        // photon loop 
        for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
            if(doBackgroundPhotons){
                if(photon_loose_p->at(pI)) continue;
            } else{
                if(!photon_tight_p->at(pI)) continue;
            }
            if(photon_pt_p->at(pI) < gammaPtBins[0]) continue;
            //if(photon_pt_p->at(pI) >= gammaPtBins[nGammaPtBins]) continue;

            //Isolation as taken from internal note of 2015 data analysis      
            if(photon_etcone30_p->at(pI) > isoCut) continue;

            Float_t etaValMain = photon_eta_p->at(pI);
            Float_t etaValSub = etaValMain;
            if(etaBinsDoAbs) etaValMain = TMath::Abs(etaValMain);
            if(etaBinsSubDoAbs) etaValSub = TMath::Abs(etaValSub);

            if(etaValMain <= etaBins[0]) continue;
            if(etaValMain >= etaBins[nEtaBins]) continue;

            Int_t ptPos = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photon_pt_p->at(pI), true, doGlobalDebug);
            //Int_t etaPos = ghostPos(nEtaBinsSub, etaBinsSub, etaValSub, true, doGlobalDebug);

            if(ptPos < 0) continue;

            float phoPt = photon_pt_p->at(pI);

            //count the number of photons in a given centrality and pt bins
            if(!isMC){
                ++(gammaCountsPerPtCent[ptPos][centPos]);
                ++(gammaCountsPerPtCent[nGammaPtBinsSub][centPos]);
            }
            else{
                gammaCountsPerPtCent[ptPos][centPos] += fullWeight;
                gammaCountsPerPtCent[nGammaPtBinsSub][centPos] += fullWeight;
                if(truthPhotonPt>0){
                    genMatchedGammaCountsPerPtCent[ptPos][centPos] += fullWeight;
                    genMatchedGammaCountsPerPtCent[nGammaPtBinsSub][centPos] += fullWeight;
                }
            }

            ///////////////////////////////////////////////////////////
            // jet loop 
            for(unsigned int jI = 0; jI < akthi_em_xcalib_jet_pt_p->size(); ++jI){

                //jet eta cut
                if(akthi_em_xcalib_jet_eta_p->at(jI) <= etaBinsLow) continue;
                if(akthi_em_xcalib_jet_eta_p->at(jI) >= etaBinsHigh) continue;

                // dR(photon,jet) cut
                Float_t dR = getDR(akthi_em_xcalib_jet_eta_p->at(jI), akthi_em_xcalib_jet_phi_p->at(jI), photon_eta_p->at(pI), photon_phi_p->at(pI));
                if(dR < DRPhotonJetCut) continue;

                // reco jet min pt check //?? revisit here!
                if(recoJtPtMin > akthi_em_xcalib_jet_pt_p->at(jI)) recoJtPtMin = akthi_em_xcalib_jet_pt_p->at(jI);

                //jet min and max pt cut
                if(akthi_em_xcalib_jet_pt_p->at(jI) < jtPtBinsLow) continue;
                if(doJetMaxPtCut && akthi_em_xcalib_jet_pt_p->at(jI) >= jtPtBinsHigh) continue;

                // dphi(photon, jet) cut    
                Float_t dPhi = TMath::Abs(getDPHI(akthi_em_xcalib_jet_phi_p->at(jI), photon_phi_p->at(pI)));
                if(doDPhiCut && dPhi < gammaJtDPhiCut) continue;

                isPhotonJetEvent = true;
                float jetPt =  akthi_em_xcalib_jet_pt_p->at(jI);
                fillTH1(PhotonTaggedJetPt_CentPhoPtDep[centPos][ptPos], jetPt, fullWeight);
                fillTH1(PhotonTaggedJetPt_CentPhoPtDep[centPos][nGammaPtBinsSub], jetPt, fullWeight);
                fillTH1(PhotonTaggedJetPtPerPho_CentPhoPtDep[centPos][ptPos], jetPt, fullWeight);
                fillTH1(PhotonTaggedJetPtPerPho_CentPhoPtDep[centPos][nGammaPtBinsSub], jetPt, fullWeight);

                fillTH2(PhotonTaggedJetPtVsPhotonPt_CentDep[centPos], phoPt, jetPt, fullWeight);

                if(isMC){
                    if(truthPhotonPt>0){// gen prompt photon pt!
                        int genJetIndex = akthi_truthpos_p->at(jI);
                        if(genJetIndex >=0){
                            float truthJetPt = akt_truth_jet_pt_p->at(genJetIndex);
                            fillTH1(GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[centPos][ptPos], jetPt, fullWeight);
                            fillTH1(GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[centPos][nGammaPtBinsSub], jetPt, fullWeight);
                            fillTH1(GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[centPos][ptPos], truthJetPt, fullWeight);
                            fillTH1(GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[centPos][nGammaPtBinsSub], truthJetPt, fullWeight);

                            fillTH2(GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[centPos], phoPt, jetPt, fullWeight);
                            fillTH2(GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[centPos], truthPhotonPt, truthJetPt, fullWeight);
                        }//if(gen jet matched)
                    }//if(gen photon matched)

                    //for total reco photons, no gen matching for photons
                    int genJetIndex = akthi_truthpos_p->at(jI);
                    if(genJetIndex >=0){
                        float truthJetPt = akt_truth_jet_pt_p->at(genJetIndex);
                        fillTH1(PhotonTaggedGenMatchedJetPt_CentPhoPtDep[centPos][ptPos], jetPt, fullWeight);
                        fillTH1(PhotonTaggedGenMatchedJetPt_CentPhoPtDep[centPos][nGammaPtBinsSub], jetPt, fullWeight);
                        fillTH1(PhotonTaggedGenJetPt_CentPhoPtDep[centPos][ptPos], truthJetPt, fullWeight);
                        fillTH1(PhotonTaggedGenJetPt_CentPhoPtDep[centPos][nGammaPtBinsSub], truthJetPt, fullWeight);
                    }
                    
                }//MC
            }//jet loop
            if(isPhotonJetEvent){
                ++(gammaJetCountsPerPtCent[ptPos][centPos]);
                ++(gammaJetCountsPerPtCent[nGammaPtBinsSub][centPos]);
            }
        }//photon loop
        if(isPhotonJetEvent){
            ++(eventCountsPerCent[centPos]);
        }
    }//event loop

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    inFile_p->Close();
    delete inFile_p;

    outFile_p->cd();

    //Pre-write and delete some of these require some mods
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
            //print number of photons! (normalization factor)
            std::cout << "gammaCountsPerPtCent (cent" << cI << ", photonPt" << pI << ") : " << gammaCountsPerPtCent[pI][cI] << std::endl;
            if(isMC) std::cout << "genMatchedGammaCountsPerPtCent (cent" << cI << ", photonPt" << pI << ") : " << genMatchedGammaCountsPerPtCent[pI][cI] << std::endl;
            
            PhotonTaggedJetPtPerPho_CentPhoPtDep[cI][pI]->Scale(1./gammaCountsPerPtCent[pI][cI]);
        }
    }
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
            //print number of photons! (normalization factor)
            std::cout << "gammaJetCountsPerPtCent (cent" << cI << ", photonPt" << pI << ") : " << gammaJetCountsPerPtCent[pI][cI] << std::endl;
        }
    }
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        std::cout << "eventCountsPerPtCent (cent" << cI << ") : " << eventCountsPerCent[cI] << std::endl;
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
    else{
        for(auto const & iter : runLumiIsFired){
            if(iter.second){
                unsigned long long run = (runLumiKey.InvertKey(iter.first))[0];

                ++(runLumiCounter[run]);
            }
        }

        double num = 0.0;
        double denom = 0.0;
        for(auto const & iter : runLumiTotal){
            num += ((double)runLumiCounter[iter.first]);
            denom += ((double)iter.second);
            double val = ((double)runLumiCounter[iter.first])/((double)iter.second);
            int binVal = lumiFractionPerRun_p->FindBin(iter.first);
            lumiFractionPerRun_p->SetBinContent(binVal, val);
            lumiFractionPerRun_p->SetBinError(binVal, 0.0);
        }

        std::cout << "Total lumiblocks: " << num << "/" << denom << "=" << num/denom << std::endl;

        lumiFractionPerRun_p->Write("", TObject::kOverwrite);
    }



    ///////////////////////////////////////////////////////////
    // Write histograms in the output file 
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        TDirectoryFile* centDir_p = (TDirectoryFile*)outFile_p->mkdir(centBinsStr[cI].c_str());
        centDir_p->cd();
        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
        PhotonTaggedJetPtVsPhotonPt_CentDep[cI]->Write("", TObject::kOverwrite);
        if(isMC){
            GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[cI]->Write("", TObject::kOverwrite);
            GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[cI]->Write("", TObject::kOverwrite);
        }

        for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
            PhotonTaggedJetPt_CentPhoPtDep[cI][pI]->Write("", TObject::kOverwrite);
            PhotonTaggedJetPtPerPho_CentPhoPtDep[cI][pI]->Write("", TObject::kOverwrite);
            if(isMC){
                GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[cI][pI]->Write("", TObject::kOverwrite);
                GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI]->Write("", TObject::kOverwrite);
                PhotonTaggedGenJetPt_CentPhoPtDep[cI][pI]->Write("", TObject::kOverwrite);
                PhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI]->Write("", TObject::kOverwrite);
            }
        }
        centDir_p->Close();
        delete centDir_p;
    }

    //Before we do a bunch of deletions lets auto-generate some stats tables
    if(!isMC){
        std::string lStr = "l | ";
        for(int gI = 0; gI < nGammaPtBinsSub; ++gI){
            lStr = lStr + "l | ";
        }
        lStr = lStr + "l";
        std::string multiColumn = std::to_string(nGammaPtBinsSub+2);

        std::string tempDPhiStr = gammaJtDPhiLabel.substr(gammaJtDPhiLabel.find(">")+2, gammaJtDPhiLabel.size());
        tempDPhiStr.replace(tempDPhiStr.find("#"),1,"\\");
        std::string jtPtStr = prettyString(jtPtBins[0], 1, false) + " < \\ptj < " + prettyString(jtPtBins[nJtPtBins], 1, false) + "; |\\dphigj| > " + tempDPhiStr;

        std::cout << "\\begin{table}" << std::endl;
        std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
        std::cout << "\\vspace{-0.6cm}" << std::endl;
        std::cout << "\\begin{center}" << std::endl;
        std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
        std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With Exactly 2 Jets (" << jtPtStr << ")} \\\\" << std::endl;
        std::string columnStr = " Centrality &";
        for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
            columnStr = columnStr + " " + prettyString(gammaPtBinsSub[pI], 1, false) + " < \\ptg < " + prettyString(gammaPtBinsSub[pI+1], 1, false) + " &";
        }
        columnStr = columnStr + " " + prettyString(gammaPtBinsSub[0], 1, false) + " < \\ptg < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false) + " \\\\ \\hline" ;

        std::cout << columnStr << std::endl;

        for(Int_t cI = 0; cI < nCentBins; ++cI){ 
            std::string outStr = binsToLabelStr[centBinsStr[cI]];
            if(outStr.find("%") != std::string::npos){
                outStr.replace(outStr.find("%"), 1, "\\%");
            }

            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                outStr = outStr + " & " + std::to_string(((int)PhotonTaggedJetPt_CentPhoPtDep[cI][pI]->GetBinContent(3)));
            }

            std::cout << outStr << " \\\\" << std::endl;
        }  
        std::cout << "\\end{tabular}" << std::endl;
        std::cout << "\\end{center}" << std::endl;
        std::cout << "\\end{table}" << std::endl;

        std::cout << std::endl;
        std::cout << "\\begin{table}" << std::endl;
        std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
        std::cout << "\\vspace{-0.6cm}" << std::endl;
        std::cout << "\\begin{center}" << std::endl;
        std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
        std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With 2 or More Jets (" << jtPtStr << ")} \\\\" << std::endl;
        std::cout << columnStr << std::endl;

        for(Int_t cI = 0; cI < nCentBins; ++cI){ 
            std::string outStr = binsToLabelStr[centBinsStr[cI]];
            if(outStr.find("%") != std::string::npos){
                outStr.replace(outStr.find("%"), 1, "\\%");
            }

            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                int countNJet = 0;

                for(Int_t bIX = 2; bIX < PhotonTaggedJetPt_CentPhoPtDep[cI][pI]->GetXaxis()->GetNbins(); ++bIX){
                    countNJet += PhotonTaggedJetPt_CentPhoPtDep[cI][pI]->GetBinContent(bIX+1);
                }

                outStr = outStr + " & " + std::to_string(countNJet);
            }

            std::cout << outStr << " \\\\" << std::endl;
        }  
        std::cout << "\\end{tabular}" << std::endl;
        std::cout << "\\end{center}" << std::endl;
        std::cout << "\\end{table}" << std::endl;

        std::cout << std::endl;
        std::cout << "\\begin{table}" << std::endl;
        std::cout << "\\fontsize{9}{9}\\selectfont" << std::endl;
        std::cout << "\\vspace{-0.6cm}" << std::endl;
        std::cout << "\\begin{center}" << std::endl;
        std::cout << "\\hspace*{-2cm}\\begin{tabular}{ " << lStr << " }" << std::endl;
        std::cout << "\\multicolumn{" << multiColumn << "}{c}{Event Counts With 1 or More Jets (" << jtPtStr << ")} \\\\" << std::endl;
        std::cout << columnStr << std::endl;

        for(Int_t cI = 0; cI < nCentBins; ++cI){ 
            std::string outStr = binsToLabelStr[centBinsStr[cI]];
            if(outStr.find("%") != std::string::npos){
                outStr.replace(outStr.find("%"), 1, "\\%");
            }

            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                int countNJet = 0;

                for(Int_t bIX = 1; bIX < PhotonTaggedJetPt_CentPhoPtDep[cI][pI]->GetXaxis()->GetNbins(); ++bIX){
                    countNJet += PhotonTaggedJetPt_CentPhoPtDep[cI][pI]->GetBinContent(bIX+1);
                }

                outStr = outStr + " & " + std::to_string(countNJet);
            }

            std::cout << outStr << " \\\\" << std::endl;
        }  
        std::cout << "\\end{tabular}" << std::endl;
        std::cout << "\\end{center}" << std::endl;
        std::cout << "\\end{table}" << std::endl;
    }

    delete runNumber_p;

    if(isMC){
        delete pthat_p;
        delete pthat_Unweighted_p;
    }
    else delete lumiFractionPerRun_p;

    if(!isPP){
        delete centrality_p;
        if(isMC) delete centrality_Unweighted_p;
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////
    // delete histograms
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        delete PhotonTaggedJetPtVsPhotonPt_CentDep[cI];
        if(isMC){
            delete GenMatchedPhotonTaggedGenJetPtVsGenPhotonPt_CentDep[cI];
            delete GenMatchedPhotonTaggedGenMatchedJetPtVsGenMatchedPhotonPt_CentDep[cI];
        }

        for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
            delete PhotonTaggedJetPt_CentPhoPtDep[cI][pI];
            delete PhotonTaggedJetPtPerPho_CentPhoPtDep[cI][pI];
            if(isMC){
                delete GenMatchedPhotonTaggedGenJetPt_CentPhoPtDep[cI][pI];
                delete GenMatchedPhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI];
                delete PhotonTaggedGenJetPt_CentPhoPtDep[cI][pI];
                delete PhotonTaggedGenMatchedJetPt_CentPhoPtDep[cI][pI];
            }
        }
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::map<std::string, std::string> configMap = config.GetConfigMap(); //grab the config in map form          
    configMap["RECOJTPTMIN"] = prettyString(recoJtPtMin, 1, false);

    TEnv configEnv; // We will convert to tenv and store in file
    for(auto const & val : configMap){
        configEnv.SetValue(val.first.c_str(), val.second.c_str()); //Fill out the map
    }
    configEnv.Write("config", TObject::kOverwrite);  

    TEnv labelEnv;
    for(auto const & lab : binsToLabelStr){
        labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
    }
    labelEnv.Write("label", TObject::kOverwrite);

    outFile_p->Close();
    delete outFile_p;

    delete randGen_p;

    std::cout << "phoTaggedJetRaa COMPLETE. return 0." << std::endl;
    return 0;
}

int main(int argc, char* argv[])
{
    if(argc != 2){
        std::cout << "Usage: ./bin/phoTaggedJetRaa.exe <inConfigFileName>" << std::endl;
        std::cout << "TO DEBUG:" << std::endl;
        std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
        std::cout << "TO TURN OFF DEBUG:" << std::endl;
        std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
        std::cout << "return 1." << std::endl;
        return 1;
    }

    int retVal = 0;
    retVal += phoTaggedJetRaa(argv[1]);
}
