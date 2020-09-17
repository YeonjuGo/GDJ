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

int phoTaggedJetRaa_photonEff(std::string inConfigFileName)
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
        "PHOGENMATCHINGDR",
        "PHOISOCONESIZE",
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
        "GAMMAJTDPHI",  
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
        "DOLOOSEID",
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
    std::string phoGenMatchingDR = config.GetConfigVal("DRJETCONE");
    const bool doMix = std::stoi(config.GetConfigVal("DOMIX"));
    const bool doLooseID = std::stoi(config.GetConfigVal("DOLOOSEID"));
    const Float_t isoCut = std::stof(config.GetConfigVal("ISOCUT"));
    const Float_t phoIsoConeSize = std::stof(config.GetConfigVal("PHOISOCONESIZE"));
    std::string label_phoIsoConeSize = Form("%d",(int)(phoIsoConeSize*100));
    std::cout << "label_phoIsoConeSize = " << label_phoIsoConeSize << std::endl;

    if(doMix && !config.ContainsParamSet(mixParams)) return 1;

    //Mixing categories, temp hardcoding
    //Centrality, percent level

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::vector<std::vector<unsigned long long> > mixVect;
    std::vector<std::vector<unsigned long long> > keyVect;

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
    outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + "_photonEff.root";

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
    const Int_t nGammaPtBinsSub = std::stoi(config.GetConfigVal("NGAMMAPTBINSSUB"));
    if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaPtBinsSub, "NGAMMAPTBINSSUB")) return 1;
    const Float_t gammaPtBinsSubLow = std::stof(config.GetConfigVal("GAMMAPTBINSSUBLOW"));
    const Float_t gammaPtBinsSubHigh = std::stof(config.GetConfigVal("GAMMAPTBINSSUBHIGH"));
    if(gammaPtBinsSubLow < gammaPtBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubLow \'" << gammaPtBinsSubLow << "\' less than gammaPtBinsLow \'" << gammaPtBinsLow << "\'. return 1" << std::endl;
    if(gammaPtBinsSubHigh > gammaPtBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubHigh \'" << gammaPtBinsSubHigh << "\' greater than gammaPtBinsHigh \'" << gammaPtBinsHigh << "\'. return 1" << std::endl;
    if(gammaPtBinsSubLow < gammaPtBinsLow || gammaPtBinsSubHigh > gammaPtBinsHigh) return 1;
    const Bool_t gammaPtBinsSubDoLog = std::stof(config.GetConfigVal("GAMMAPTBINSDOLOG"));
    Double_t gammaPtBinsSub[nMaxSubBins+1];
    if(gammaPtBinsSubDoLog) getLogBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
    else getLinBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
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

    const Int_t nJtPtBins = std::stoi(config.GetConfigVal("NJTPTBINS"));
    if(!goodBinning(inConfigFileName, nMaxPtBins, nJtPtBins, "NJTPTBINS")) return 1;
    const Float_t jtPtBinsLow = std::stof(config.GetConfigVal("JTPTBINSLOW"));
    const Float_t jtPtBinsHigh = std::stof(config.GetConfigVal("JTPTBINSHIGH"));
    const Bool_t jtPtBinsDoLog = std::stof(config.GetConfigVal("JTPTBINSDOLOG"));
    Double_t jtPtBins[nMaxPtBins+1];
    if(jtPtBinsDoLog) getLogBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);
    else getLinBins(jtPtBinsLow, jtPtBinsHigh, nJtPtBins, jtPtBins);

    std::string jtPtBinsGlobalStr = "GlobalJtPt0";
    std::string jtPtBinsGlobalLabel = "p_{T,jet} > " + prettyString(jtPtBinsLow,1,false);
    //std::string jtPtBinsGlobalLabel = prettyString(jtPtBinsLow,1,false) + " < p_{T,jet} < " + prettyString(jtPtBinsHigh,1,false);
    binsToLabelStr[jtPtBinsGlobalStr] = jtPtBinsGlobalLabel;   

    std::string multiJtCutGlobalStr = "MultiJt0";
    std::string multiJtCutGlobalLabel = "N_{Jet,Reco.} >= 2";
    binsToLabelStr[multiJtCutGlobalStr] = multiJtCutGlobalLabel;   

    std::vector<std::string> jtPtBinsStr;
    for(Int_t pI = 0; pI < nJtPtBins; ++pI){
        if(pI < 10) jtPtBinsStr.push_back("JtPt0" + std::to_string(pI));
        else jtPtBinsStr.push_back("JtPt" + std::to_string(pI));

        binsToLabelStr[jtPtBinsStr[pI]] = " p_{T,Jet} > " + prettyString(jtPtBins[pI], 1, false);
        //binsToLabelStr[jtPtBinsStr[pI]] = prettyString(jtPtBins[pI], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[pI+1], 1, false);
    }
    jtPtBinsStr.push_back("JtPt" + std::to_string(nJtPtBins));
    binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] =" p_{T,Jet} > " + prettyString(jtPtBins[0], 1, false);
    //binsToLabelStr[jtPtBinsStr[jtPtBinsStr.size()-1]] = prettyString(jtPtBins[0], 1, false) + " < p_{T,Jet} < " + prettyString(jtPtBins[nJtPtBins], 1, false);


    ////////////////////////////////////////////////////////////
    // define histograms  

    TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
    TH1F* runNumber_p = nullptr;
    TH1F* pthat_p = nullptr;
    TH1F* pthat_Unweighted_p = nullptr;
    TH1F* centrality_p = nullptr;
    TH1F* centrality_Unweighted_p = nullptr;

    const int nPhoEtaBins = 2;
    const double etaBins_i[] = {0.0, 1.52};
    const double etaBins_f[] = {1.37, 2.37};
    std::vector<std::string> etaBinsStr;
    for(Int_t eI = 0; eI < nPhoEtaBins; ++eI)
        etaBinsStr.push_back("Eta" + std::to_string(etaBins_i[eI]) + "to" + std::to_string(etaBins_f[eI]));

    const int nPhoPtBins = 7;
    const double ptBins_i[] = {50, 60, 80, 100, 120, 150, 200};
    const double ptBins_f[] = {60, 80, 100, 120, 150, 200, 250};
    std::vector<std::string> ptBinsStr;
    for(Int_t pI = 0; pI < nPhoPtBins; ++pI)
        ptBinsStr.push_back("PhotonPt" + std::to_string(ptBins_i[pI]) + "to" + std::to_string(ptBins_f[pI]));
    ptBinsStr.push_back("PhotonPt50to250");
    //ptBinsStr.push_back("PhotonPt" + "50" + "to" + "250");
    const float minIso = -10;
    const float maxIso = 50;
    const int nIso = maxIso - minIso;

    TH1F* photonEff_TOT_CentDep_Den[nMaxCentBins][nPhoEtaBins];// gen matching
    TH1F* photonEff_TOT_CentDep_Num[nMaxCentBins][nPhoEtaBins];// gen matching
    TH1F* photonEff_TOT_CentDep_Eff[nMaxCentBins][nPhoEtaBins];// gen matching
    TH1F* photonEff_ID_CentDep_Num[nMaxCentBins][nPhoEtaBins];// gen matching
    TH1F* photonEff_ID_CentDep_Eff[nMaxCentBins][nPhoEtaBins];// gen matching
    TH1F* photonEff_ISO_CentDep_Num[nMaxCentBins][nPhoEtaBins];// gen matching
    TH1F* photonEff_ISO_CentDep_Eff[nMaxCentBins][nPhoEtaBins];// gen matching
    TH2F* photon_truthIso_vs_recoIso[nMaxCentBins][nPhoEtaBins][nPhoPtBins+1];// gen matching
    //TH1F* photonEff_ISO_CentDep_Den[nMaxCentBins];// gen matching
    //TH1F* photonEff_ID_CentDep_Den[nMaxCentBins];// gen matching


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
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

            if(isMC){
                photonEff_TOT_CentDep_Den[cI][eI] = new TH1F(("photonEff_TOT_CentDep_Den" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBins, gammaPtBins);
                photonEff_TOT_CentDep_Num[cI][eI] = new TH1F(("photonEff_TOT_CentDep_Num" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBins, gammaPtBins);
                photonEff_TOT_CentDep_Eff[cI][eI] = new TH1F(("photonEff_TOT_CentDep_Eff" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Total Efficiency", nGammaPtBins, gammaPtBins);

                photonEff_ID_CentDep_Num[cI][eI] = new TH1F(("photonEff_ID_CentDep_Num" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBins, gammaPtBins);
                photonEff_ID_CentDep_Eff[cI][eI] = new TH1F(("photonEff_ID_CentDep_Eff" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Identification Efficiency", nGammaPtBins, gammaPtBins);

                photonEff_ISO_CentDep_Num[cI][eI] = new TH1F(("photonEff_ISO_CentDep_Num" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBins, gammaPtBins);
                photonEff_ISO_CentDep_Eff[cI][eI] = new TH1F(("photonEff_ISO_CentDep_Eff" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Isolation Efficiency", nGammaPtBins, gammaPtBins);

                for(Int_t pI = 0; pI < nPhoPtBins+1; ++pI){
                    photon_truthIso_vs_recoIso[cI][eI][pI] = new TH2F(("photon_truthIso_vs_recoIso" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + ptBinsStr[pI] + "_h").c_str(), ";#gamma gen-level isolation [GeV];#gamma reco-level isolation [GeV]", nIso, minIso, maxIso, nIso, minIso, maxIso);
                    centerTitles({photon_truthIso_vs_recoIso[cI][eI][pI]});
                    setSumW2({photon_truthIso_vs_recoIso[cI][eI][pI]});

                }

                centerTitles({photonEff_ID_CentDep_Num[cI][eI],photonEff_ID_CentDep_Eff[cI][eI]});
                setSumW2({photonEff_ID_CentDep_Num[cI][eI],photonEff_ID_CentDep_Eff[cI][eI]});
                centerTitles({photonEff_ISO_CentDep_Num[cI][eI],photonEff_ISO_CentDep_Eff[cI][eI]});
                setSumW2({photonEff_ISO_CentDep_Num[cI][eI],photonEff_ISO_CentDep_Eff[cI][eI]});
                centerTitles({photonEff_TOT_CentDep_Den[cI][eI],photonEff_TOT_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Eff[cI][eI]});
                setSumW2({photonEff_TOT_CentDep_Den[cI][eI],photonEff_TOT_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Eff[cI][eI]});
            }//MC
        }//eta loop

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
    centerTitles(runNumber_p);

    inFile_p->cd();

    //Grab the hltbranches for some basic prescale checks
    std::vector<std::string> listOfBranches = getVectBranchList(inTree_p);
    //std::vector<std::string> hltList;
    //std::vector<std::string> hltListPres;
    //std::string hltStr = "HLT_";
    //std::string prescaleStr = "_prescale";

    ////HLTLists are built
    //for(auto const & branchStr : listOfBranches){
    //    if(branchStr.size() < hltStr.size()) continue;
    //    if(!isStrSame(branchStr.substr(0, hltStr.size()), hltStr)) continue;

    //    if(branchStr.find(prescaleStr) != std::string::npos) hltListPres.push_back(branchStr);
    //    else hltList.push_back(branchStr);
    //}

    //bool allHLTPrescalesFound = true;
    //for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    //    std::string modStr = hltList[hI] + prescaleStr;
    //    if(!vectContainsStr(modStr, &hltListPres)){
    //        allHLTPrescalesFound = false;
    //        std::cout << "HLT " << hltList[hI] << " has no prescale. return 1" << std::endl;
    //    }
    //}
    //if(!allHLTPrescalesFound) return 1;

    //float hltPrescaleDelta = 0.01;
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
    std::vector<int>* truth_type_p=nullptr;
    std::vector<int>* truth_origin_p=nullptr;

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
    std::vector<float>* photon_etcone_p=nullptr;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    inTree_p->SetBranchStatus("*", 0);

    //for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    //    hltVect.push_back(new bool(false));
    //    hltPrescaleVect.push_back(new float(0.0));

    //    inTree_p->SetBranchStatus(hltList[hI].c_str(), 1);
    //    inTree_p->SetBranchStatus(hltListPres[hI].c_str(), 1);
    //}  

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
        inTree_p->SetBranchStatus("truth_type", 1);
        inTree_p->SetBranchStatus("truth_origin", 1);

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
    }

    inTree_p->SetBranchStatus("vert_z", 1);

    inTree_p->SetBranchStatus("photon_pt", 1);
    inTree_p->SetBranchStatus("photon_eta", 1);
    inTree_p->SetBranchStatus("photon_phi", 1);
    inTree_p->SetBranchStatus("photon_tight", 1);
    inTree_p->SetBranchStatus("photon_loose", 1);
    inTree_p->SetBranchStatus(("photon_etcone"+label_phoIsoConeSize).c_str(), 1);

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    //inTree_p->SetBranchStatus(("akt"+jetDR+"hi_em_xcalib_jet_pt").c_str(), 1);
    //inTree_p->SetBranchStatus(("akt"+jetDR+"hi_em_xcalib_jet_eta").c_str(), 1);
    //inTree_p->SetBranchStatus(("akt"+jetDR+"hi_em_xcalib_jet_phi").c_str(), 1);

    //if(isMC){
    //    inTree_p->SetBranchStatus(("akt"+jetDR+"hi_truthpos").c_str(), 1);

    //    inTree_p->SetBranchStatus(("akt"+jetDR+"_truth_jet_pt").c_str(), 1);
    //    inTree_p->SetBranchStatus(("akt"+jetDR+"_truth_jet_eta").c_str(), 1);
    //    inTree_p->SetBranchStatus(("akt"+jetDR+"_truth_jet_phi").c_str(), 1);
    //}

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    //for(unsigned int hI = 0; hI < hltList.size(); ++hI){
    //    inTree_p->SetBranchAddress(hltList[hI].c_str(), hltVect[hI]);
    //    inTree_p->SetBranchAddress(hltListPres[hI].c_str(), hltPrescaleVect[hI]);
    //}

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
        inTree_p->SetBranchAddress("truth_type", &truth_type_p);
        inTree_p->SetBranchAddress("truth_origin", &truth_origin_p);

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
    }

    inTree_p->SetBranchAddress("vert_z", &vert_z_p);

    inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
    inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
    inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
    inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
    inTree_p->SetBranchAddress("photon_loose", &photon_loose_p);

    inTree_p->SetBranchAddress(("photon_etcone"+label_phoIsoConeSize).c_str(), &photon_etcone_p);

    //inTree_p->SetBranchAddress(("akt"+jetDR+"hi_em_xcalib_jet_pt").c_str(), &akthi_em_xcalib_jet_pt_p);
    //inTree_p->SetBranchAddress(("akt"+jetDR+"hi_em_xcalib_jet_eta").c_str(), &akthi_em_xcalib_jet_eta_p);
    //inTree_p->SetBranchAddress(("akt"+jetDR+"hi_em_xcalib_jet_phi").c_str(), &akthi_em_xcalib_jet_phi_p);

    //if(isMC){
    //    inTree_p->SetBranchAddress(("akt"+jetDR+"hi_truthpos").c_str(), &akthi_truthpos_p);    

    //    inTree_p->SetBranchAddress(("akt"+jetDR+"_truth_jet_pt").c_str(), &akt_truth_jet_pt_p);
    //    inTree_p->SetBranchAddress(("akt"+jetDR+"_truth_jet_eta").c_str(), &akt_truth_jet_eta_p);
    //    inTree_p->SetBranchAddress(("akt"+jetDR+"_truth_jet_phi").c_str(), &akt_truth_jet_phi_p);
    //}


    //Double_t recoJtPtMin = 100000.;
    const ULong64_t nEntries = inTree_p->GetEntries();
    const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

    //variable to count the number of events in a given centrality bin
    std::vector<Double_t> eventCountsPerCent;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        eventCountsPerCent.push_back(0.0);
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
    std::cout << "Processing " << nEntries << " events..." << std::endl;

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

        //unsigned long long tempKey = runLumiKey.GetKey({(unsigned long long)runNumber, (unsigned long long)lumiBlock});    
        //runLumiIsFired[tempKey] = true;

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

        ///////////////////////////////////////////////////////////
        // photon loop 

        for(unsigned int pI = 0; pI < photon_pt_p->size(); ++pI){
            if(truthPhotonPt<=0) continue; //prompt isolated photons?!
            if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) > phoGenMatchingDR) continue; 
            // if(abs(truthPhotonEta)>2.4) continue;
            //for (; abs(truthPhotonEta)>=etaBins_i[tempEtaPos] && tempEtaPos<etaBins_f[tempEtaPos] ; ++tempEtaPos);
            int tempEtaPos = -1;
            for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
                if(abs(truthPhotonEta)>=etaBins_i[eI] && abs(truthPhotonEta)<etaBins_f[eI]) tempEtaPos=eI; 
            }
            if(tempEtaPos==-1) continue; //eta cut

            if(photon_pt_p->at(pI) < gammaPtBins[0]) continue; // min gamma pT cut
            //if(photon_pt_p->at(pI) >= gammaPtBins[nGammaPtBins]) continue;

            fillTH1(photonEff_TOT_CentDep_Den[centPos][tempEtaPos],truthPhotonPt,fullWeight);

            if(doLooseID){
                if(photon_loose_p->at(pI))
                    fillTH1(photonEff_ID_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
            } else{
                if(photon_tight_p->at(pI))
                    fillTH1(photonEff_ID_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
            } 

            if(photon_etcone_p->at(pI) <= isoCut){
                fillTH1(photonEff_ISO_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
                if(doLooseID && photon_loose_p->at(pI))
                    fillTH1(photonEff_TOT_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
                else if(!doLooseID && photon_tight_p->at(pI))
                    fillTH1(photonEff_TOT_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
            }
        }//photon loop
    }//event loop

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    inFile_p->Close();
    delete inFile_p;

    outFile_p->cd();

    //Pre-write and delete some of these require some mods
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            photonEff_ID_CentDep_Eff[cI][eI]->Divide(photonEff_ID_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Den[cI][eI],1.,1.,"B");
            photonEff_ISO_CentDep_Eff[cI][eI]->Divide(photonEff_ISO_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Den[cI][eI],1.,1.,"B");
            photonEff_TOT_CentDep_Eff[cI][eI]->Divide(photonEff_TOT_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Den[cI][eI],1.,1.,"B");
        }
    }

    if(!isPP){
        centrality_p->Write("", TObject::kOverwrite);
        if(isMC) centrality_Unweighted_p->Write("", TObject::kOverwrite);
    }

    if(isMC){
        pthat_p->Write("", TObject::kOverwrite);
        pthat_Unweighted_p->Write("", TObject::kOverwrite);
    }


    ///////////////////////////////////////////////////////////
    // Write histograms in the output file 
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
            photonEff_ID_CentDep_Num[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_ID_CentDep_Eff[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_ISO_CentDep_Num[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_ISO_CentDep_Eff[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_TOT_CentDep_Num[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_TOT_CentDep_Den[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_TOT_CentDep_Eff[cI][eI]->Write("", TObject::kOverwrite);
        }
    }

    if(isMC){
        delete pthat_p;
        delete pthat_Unweighted_p;
    }

    if(!isPP){
        delete centrality_p;
        if(isMC) delete centrality_Unweighted_p;
    }

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////
    // delete histograms
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            delete photonEff_ID_CentDep_Num[cI][eI];
            delete photonEff_ID_CentDep_Eff[cI][eI];
            delete photonEff_ISO_CentDep_Num[cI][eI];
            delete photonEff_ISO_CentDep_Eff[cI][eI];
            delete photonEff_TOT_CentDep_Num[cI][eI];
            delete photonEff_TOT_CentDep_Den[cI][eI];
            delete photonEff_TOT_CentDep_Eff[cI][eI];
        }
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::map<std::string, std::string> configMap = config.GetConfigMap(); //grab the config in map form          
    //configMap["RECOJTPTMIN"] = prettyString(recoJtPtMin, 1, false);

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

    std::cout << "phoTaggedJetRaa_photonEff COMPLETE. return 0." << std::endl;
    return 0;
}

int main(int argc, char* argv[])
{
    if(argc != 2){
        std::cout << "Usage: ./bin/phoTaggedJetRaa_photonEff.exe <inConfigFileName>" << std::endl;
        std::cout << "TO DEBUG:" << std::endl;
        std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
        std::cout << "TO TURN OFF DEBUG:" << std::endl;
        std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
        std::cout << "return 1." << std::endl;
        return 1;
    }

    int retVal = 0;
    retVal += phoTaggedJetRaa_photonEff(argv[1]);
}
