//Author: Chris McGinn (2020.02.19)
//Modified by: Yeonju Go (2020.09.18)
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
#include "TF1.h"

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
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/photonUtil.h"
#include "include/treeUtil.h"
#include "include/toStringWithPrecision.h"
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

int phoTaggedJetRaa_photonEff(std::string inConfigFileName)
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
                                                "SYSTEMATIC",
                                                "CENTFILENAME",
                                                "ISPP",
                                                "ISMC",
        "PHOGENMATCHINGDR",
        "PHOISOCONESIZE",
        "GENISOCUT",
        "ISOCUT",
        "DOPTCORRECTEDISO",
        "DOCENTCORRECTEDISO",
        "DOLOOSEID_PHOTONEFF",
        "CENTBINS",
        "DOPTBINSINCONFIG",
        "PTBINS",
        "NGAMMAPTBINS",
        "GAMMAPTBINSLOW",
        "GAMMAPTBINSHIGH",
        "GAMMAPTBINSDOLOG",
        "DOPTBINSSUBINCONFIG",
        "GAMMAPTBINSSUB",
        "NGAMMAPTBINSSUB",
        "GAMMAPTBINSSUBLOW",
        "GAMMAPTBINSSUBHIGH",
        "GAMMAPTBINSSUBDOLOG"
        };


    if(!checkEnvForParams(config_p, necessaryParams)) return 1;

    std::string inCentFileName = config_p->GetValue("CENTFILENAME","");
    std::string version = config_p->GetValue("VERSION","temp");
    std::string systematic = config_p->GetValue("SYSTEMATIC","nominal");
    const bool doLooseID = config_p->GetValue("DOLOOSEID_PHOTONEFF",0);
    const bool doPtCorrectedIso = config_p->GetValue("DOPTCORRECTEDISO", 1);
    const bool doCentCorrectedIso = config_p->GetValue("DOCENTCORRECTEDISO", 1);
    const bool doPtBinInConfig = config_p->GetValue("DOPTBINSINCONFIG", 1);
    const bool doPtBinSubInConfig = config_p->GetValue("DOPTBINSSUBINCONFIG", 1);
    const Float_t isoCut = config_p->GetValue("ISOCUT", 3);
    const Float_t genIsoCut = config_p->GetValue("GENISOCUT", 5);
    const Float_t phoGenMatchingDR = config_p->GetValue("PHOGENMATCHINGDR", 0.2);
    const Float_t phoIsoConeSize = config_p->GetValue("PHOISOCONESIZE",3);
    std::string label_phoIsoConeSize = Form("%d",(int)(phoIsoConeSize));
    std::cout << "label_phoIsoConeSize = " << label_phoIsoConeSize << std::endl;

    const bool isPP = config_p->GetValue("ISPP", 1);
    const bool isMC = config_p->GetValue("ISMC", 1);
    if(!isMC) return 1;

    ////////////////////////////////////////
    // output file name
    check.doCheckMakeDir("output"); // check output dir exists; if not create
    check.doCheckMakeDir("output/" + version); // check dated output subdir exists; if not create

    std::string systStr = "PP";
    if(!isPP) systStr = "PbPb";
    std::string capStr = version + "_" + systematic;
    std::string outFileName = "output/" + version + "/phoTagJetRaa_photonEfficiency_" + systStr + "MC_" + capStr + ".root ";

    centralityFromInput centTable(inCentFileName);
    if(doGlobalDebug) centTable.PrintTableTex();
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////////
    // Centrality binning 
    if(!check.checkFileExt(inCentFileName, "txt")) return 1; // Check centrality table is valid TXT file   
    const Int_t nMaxSubBins = 20;
    const Int_t nMaxCentBins = 15;
    Int_t nCentBins = 1;

    std::vector<int> centBins;
    std::vector<std::string> centBinsStr = {systStr};
    std::map<std::string, std::string> binsToLabelStr;
    Double_t centBinsArr[nMaxSubBins+1];
    if(!isPP){
        centBins = strToVectI(config_p->GetValue("CENTBINS", "0,10,30,80"));
        nCentBins = centBins.size()-1;

        for(Int_t icent=0;icent<nCentBins+1;++icent){
            centBinsArr[icent] = centBins.at(icent);
        }

        if(!goodBinning(inConfigFileName, nMaxCentBins, nCentBins, "CENTBINS")) return 1;

        centBinsStr.clear();
        for(Int_t cI = 0; cI < nCentBins+1; ++cI){
            if(cI==nCentBins){ 
                centBinsStr.push_back("Cent" + std::to_string(centBins[0]) +"to" + std::to_string(centBins[nCentBins]));
                binsToLabelStr[centBinsStr[cI]] = Form("%d-%d%s", centBins[0], centBins[nCentBins], "%");
            } else{ 
                centBinsStr.push_back("Cent" + std::to_string(centBins[cI]) + "to" + std::to_string(centBins[cI+1]));
                binsToLabelStr[centBinsStr[cI]] = std::to_string(centBins[cI]) + "-" + std::to_string(centBins[cI+1]) + "%";
            }
        }
    } 
    else binsToLabelStr[centBinsStr[0]] = "pp";
    Int_t nCentBins_withIncBin = nCentBins+1;
    if(isPP) nCentBins_withIncBin = 1;

    ///////////////////////////////////////////////////////////////
    // photon pT main binning
    const Int_t nMaxPtBins = 300;
    const Int_t nGammaPtBins = config_p->GetValue("NGAMMAPTBINS",10);
    if(!goodBinning(inConfigFileName, nMaxPtBins, nGammaPtBins, "NGAMMAPTBINS")) return 1;
    const Float_t gammaPtBinsLow = config_p->GetValue("GAMMAPTBINSLOW",50.0);
    const Float_t gammaPtBinsHigh = config_p->GetValue("GAMMAPTBINSHIGH",1000.0);
    const Bool_t gammaPtBinsDoLog = config_p->GetValue("GAMMAPTBINSDOLOG",1);
    Double_t gammaPtBins[nMaxPtBins+1];
    if(!doPtBinInConfig){
        if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
        else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
    } else{
        std::vector<int> ptBins_vec = strToVectI(config_p->GetValue("PTBINS", "50,55,60,70,90,130,1000"));
        Int_t tempPtbinSize =  ptBins_vec.size();
        for(Int_t ipt=0;ipt<tempPtbinSize;++ipt){
            gammaPtBins[ipt] = ptBins_vec.at(ipt);
        }
    }
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

    ///////////////////////////////////////////////////////////////
    //photon pT SUB bins handling
    const Int_t nGammaPtBinsSub = config_p->GetValue("NGAMMAPTBINSSUB", 9);
    if(!goodBinning(inConfigFileName, nMaxSubBins, nGammaPtBinsSub, "NGAMMAPTBINSSUB")) return 1;
    const Float_t gammaPtBinsSubLow = config_p->GetValue("GAMMAPTBINSSUBLOW", 50);
    const Float_t gammaPtBinsSubHigh = config_p->GetValue("GAMMAPTBINSSUBHIGH", 1000);
    if(gammaPtBinsSubLow < gammaPtBinsLow) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubLow \'" << gammaPtBinsSubLow << "\' less than gammaPtBinsLow \'" << gammaPtBinsLow << "\'. return 1" << std::endl;
    if(gammaPtBinsSubHigh > gammaPtBinsHigh) std::cout << "ERROR - config \'" << inConfigFileName << "\' contains gammaPtBinsSubHigh \'" << gammaPtBinsSubHigh << "\' greater than gammaPtBinsHigh \'" << gammaPtBinsHigh << "\'. return 1" << std::endl;
    if(gammaPtBinsSubLow < gammaPtBinsLow || gammaPtBinsSubHigh > gammaPtBinsHigh) return 1;
    const Bool_t gammaPtBinsSubDoLog = config_p->GetValue("GAMMAPTBINSDOLOG", 1);
    Double_t gammaPtBinsSub[nGammaPtBinsSub+1];
    if(!doPtBinSubInConfig){
    if(gammaPtBinsSubDoLog) getLogBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
    else getLinBins(gammaPtBinsSubLow, gammaPtBinsSubHigh, nGammaPtBinsSub, gammaPtBinsSub);
    } else{
        std::vector<int> ptBins_vec = strToVectI(config_p->GetValue("GAMMAPTBINSSUB", "50,55,60,70,90,130,1000" ));
        Int_t tempPtbinSize =  ptBins_vec.size();
        for(Int_t ipt=0;ipt<tempPtbinSize;++ipt){
            gammaPtBinsSub[ipt] = ptBins_vec.at(ipt);
        }
    }
    std::vector<std::string> gammaPtBinsSubStr;
    for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
        gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(pI));
        binsToLabelStr[gammaPtBinsSubStr[pI]] = prettyString(gammaPtBinsSub[pI], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[pI+1], 1, false);
    }
    gammaPtBinsSubStr.push_back("GammaPt" + std::to_string(nGammaPtBinsSub));
    binsToLabelStr[gammaPtBinsSubStr[gammaPtBinsSubStr.size()-1]] = prettyString(gammaPtBinsSub[0], 1, false) + " < p_{T,#gamma} < " + prettyString(gammaPtBinsSub[nGammaPtBinsSub], 1, false);


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

    /////////////////////////////////////////
    // isolation binning
    const float minIso = -30;
    const float maxIso = 30;
    const int nIso = (maxIso - minIso)*4;

    ////////////////////////////////////////////////////////////
    // define histograms  

    TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
    TH1F* runNumber_p = nullptr;
    TH1F* pthat_p = nullptr;
    TH1F* pthat_Unweighted_p = nullptr;
    TH1F* centrality_p = nullptr;
    TH1F* centrality_Unweighted_p = nullptr;

    TH1F* photonEff_TOT_CentDep_Den[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photonEff_TOT_CentDep_Num[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photonEff_TOT_CentDep_Eff[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photonEff_ID_CentDep_Num[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photonEff_ID_CentDep_Eff[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photonEff_ISO_CentDep_Num[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photonEff_ISO_CentDep_Eff[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photon_arith_meanRecoIso_vs_pt[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photon_gaus_meanRecoIso_vs_pt[nMaxCentBins+1][nPhoEtaBins];// gen matching
    TH1F* photon_arith_meanRecoIso_vs_cent[nGammaPtBinsSub+1][nPhoEtaBins];// gen matching
    TH1F* photon_gaus_meanRecoIso_vs_cent[nGammaPtBinsSub+1][nPhoEtaBins];// gen matching

    TH2F* h2F_photonEff_ISO_Den[nPhoEtaBins];// gen matching
    TH2F* h2F_photonEff_ISO_Num[nPhoEtaBins];// gen matching
    TH2F* h2F_photonEff_ISO_Eff[nPhoEtaBins];// gen matching
    TH2F* photon_truthIso_vs_recoIso[nMaxCentBins+1][nPhoEtaBins][nGammaPtBinsSub+1];// gen matching
    TH1F* h1F_photon_recoIso[nMaxCentBins+1][nPhoEtaBins][nGammaPtBinsSub+1];// gen matching
    TF1* fit_photon_recoIso[nMaxCentBins+1][nPhoEtaBins][nGammaPtBinsSub+1];// gen matching
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    //TString ptLabel[nGAMMAPTBINS+1];
    //for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
    //    ptLabel[pI] = Form("%d < E_{T} " );
    //}
    ///////////////////////////////////////////////////////////
    // set histograms 
    for(Int_t cI = 0; cI < nCentBins_withIncBin; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
            std::cout << "cI = " << cI << ", eI = " << eI << std::endl;

            photonEff_TOT_CentDep_Den[cI][eI] = new TH1F(("photonEff_TOT_CentDep_Den_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBinsSub, gammaPtBinsSub);
            photonEff_TOT_CentDep_Num[cI][eI] = new TH1F(("photonEff_TOT_CentDep_Num_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBinsSub, gammaPtBinsSub);
            photonEff_TOT_CentDep_Eff[cI][eI] = new TH1F(("photonEff_TOT_CentDep_Eff_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Total Efficiency", nGammaPtBinsSub, gammaPtBinsSub);
            
            if(cI==0 && !isPP){
                h2F_photonEff_ISO_Den[eI]  = new TH2F(("h2F_photonEff_ISO_Den_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Centrality [%]", 10, 50, 250, nCentBins, centBinsArr);
                h2F_photonEff_ISO_Num[eI]  = new TH2F(("h2F_photonEff_ISO_Num_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Centrality [%]", 10, 50, 250, nCentBins, centBinsArr);
                h2F_photonEff_ISO_Eff[eI]  = new TH2F(("h2F_photonEff_ISO_Eff_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Centrality [%]", 10, 50, 250, nCentBins, centBinsArr);
                centerTitles({h2F_photonEff_ISO_Eff[eI], h2F_photonEff_ISO_Den[eI], h2F_photonEff_ISO_Num[eI]});
                setSumW2({h2F_photonEff_ISO_Eff[eI], h2F_photonEff_ISO_Den[eI], h2F_photonEff_ISO_Num[eI]});
            }

            photonEff_ID_CentDep_Num[cI][eI] = new TH1F(("photonEff_ID_CentDep_Num_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBinsSub, gammaPtBinsSub);
            photonEff_ID_CentDep_Eff[cI][eI] = new TH1F(("photonEff_ID_CentDep_Eff_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Identification Efficiency", nGammaPtBinsSub, gammaPtBinsSub);

            photonEff_ISO_CentDep_Num[cI][eI] = new TH1F(("photonEff_ISO_CentDep_Num_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Arbitrary normalized", nGammaPtBinsSub, gammaPtBinsSub);
            photonEff_ISO_CentDep_Eff[cI][eI] = new TH1F(("photonEff_ISO_CentDep_Eff_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";Gen #gamma p_{T} [GeV];Isolation Efficiency", nGammaPtBinsSub, gammaPtBinsSub);
            photon_arith_meanRecoIso_vs_pt[cI][eI] = new TH1F(("photon_arith_meanRecoIso_vs_pt_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";p_{T}^{#gamma,reco} [GeV];<Iso E_{T}^{reco}>", nGammaPtBinsSub, gammaPtBinsSub);

            photon_gaus_meanRecoIso_vs_pt[cI][eI] = new TH1F(("photon_gaus_meanRecoIso_vs_pt_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_h").c_str(), ";p_{T}^{#gamma,reco} [GeV];<Iso E_{T}^{reco}>", nGammaPtBinsSub, gammaPtBinsSub);

            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                photon_truthIso_vs_recoIso[cI][eI][pI] = new TH2F(("photon_truthIso_vs_recoIso_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";#gamma gen-level isolation [GeV];#gamma reco-level isolation [GeV]", nIso, minIso, maxIso, nIso, minIso, maxIso);
                h1F_photon_recoIso[cI][eI][pI] = new TH1F(("h1F_photon_recoIso_" + centBinsStr[cI] + "_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";#gamma reco-level isolation [GeV];Entries", nIso, minIso, maxIso);
                centerTitles({photon_truthIso_vs_recoIso[cI][eI][pI],h1F_photon_recoIso[cI][eI][pI]});
                setSumW2({photon_truthIso_vs_recoIso[cI][eI][pI],h1F_photon_recoIso[cI][eI][pI]});

                if(cI==0){ 
                    photon_arith_meanRecoIso_vs_cent[pI][eI] = new TH1F(("photon_arith_meanRecoIso_vs_cent_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";Centrality [%];<Iso E_{T}^{reco}>", nCentBins, centBinsArr);
                    photon_gaus_meanRecoIso_vs_cent[pI][eI] = new TH1F(("photon_gaus_meanRecoIso_vs_cent_" + etaBinsStr[eI] + "_" + gammaPtBinsSubStr[pI] + "_h").c_str(), ";Centrality [%];<Iso E_{T}^{reco}>", nCentBins, centBinsArr);
                    centerTitles({photon_arith_meanRecoIso_vs_cent[pI][eI], photon_gaus_meanRecoIso_vs_cent[pI][eI]});
                    setSumW2({photon_arith_meanRecoIso_vs_cent[pI][eI], photon_gaus_meanRecoIso_vs_cent[pI][eI]});

                }
            }
            centerTitles({photonEff_ID_CentDep_Num[cI][eI],photonEff_ID_CentDep_Eff[cI][eI]});
            setSumW2({photonEff_ID_CentDep_Num[cI][eI],photonEff_ID_CentDep_Eff[cI][eI]});
            centerTitles({photonEff_ISO_CentDep_Num[cI][eI],photonEff_ISO_CentDep_Eff[cI][eI]});
            setSumW2({photonEff_ISO_CentDep_Num[cI][eI],photonEff_ISO_CentDep_Eff[cI][eI]});
            centerTitles({photonEff_TOT_CentDep_Den[cI][eI],photonEff_TOT_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Eff[cI][eI]});
            setSumW2({photonEff_TOT_CentDep_Den[cI][eI],photonEff_TOT_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Eff[cI][eI]});
            centerTitles({photon_arith_meanRecoIso_vs_pt[cI][eI],photon_gaus_meanRecoIso_vs_pt[cI][eI]});
            setSumW2({photon_arith_meanRecoIso_vs_pt[cI][eI],photon_gaus_meanRecoIso_vs_pt[cI][eI]});
            
        }//eta loop
    }//centrailty loop
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////
    // pthat and centrality reweighting histograms
    pthat_p = new TH1F(("pthat_" + systStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    pthat_Unweighted_p = new TH1F(("pthat_Unweighted_" + systStr + "_h").c_str(), ";p_{T} Hat;Counts", 250, 35, 535);
    centerTitles({pthat_p, pthat_Unweighted_p});

    if(!isPP){
        centrality_p = new TH1F(("centrality_" + systStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
        centerTitles(centrality_p);

        centrality_Unweighted_p = new TH1F(("centrality_Unweighted_" + systStr + "_h").c_str(), ";Centrality (%);Counts", 100, -0.5, 99.5);
        centerTitles(centrality_Unweighted_p);
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;


    ///////////////////////////////////////////////////////////
    // import input trees 
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

    //outFile_p->cd();
    runNumber_p = new TH1F(("runNumber_" + systStr + "_h").c_str(), ";Run;Counts", nRunBins+1, runMinF, runMaxF);
    centerTitles(runNumber_p);
    //inFile_p->cd();

    //Grab the hltbranches for some basic prescale checks
    std::vector<std::string> listOfBranches = getVectBranchList(inTree_p);
    Int_t runNumber;
    UInt_t lumiBlock;
    Float_t pthat;
    Float_t sampleWeight;
    Float_t ncollWeight;
    Float_t fullWeight;
    Float_t fcalA_et, fcalC_et;
    std::vector<float>* vert_z_p=nullptr;

    float treePartonPt[2];
    float treePartonEta[2];
    float treePartonPhi[2];
    int treePartonId[2];

    Float_t truthPhotonPt, truthPhotonPhi, truthPhotonEta;
    Float_t truthPhotonIso;

    std::vector<float>* photon_pt_p=nullptr;
    std::vector<float>* photon_eta_p=nullptr;
    std::vector<float>* photon_phi_p=nullptr;
    std::vector<bool>* photon_tight_p=nullptr;  
    std::vector<bool>* photon_loose_p=nullptr;  
    std::vector<float>* photon_etcone_p=nullptr;

    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    inTree_p->SetBranchStatus("*", 0);
    inTree_p->SetBranchStatus("runNumber", 1);
    inTree_p->SetBranchStatus("lumiBlock", 1);

    if(isMC){
        inTree_p->SetBranchStatus("pthat", 1);
        inTree_p->SetBranchStatus("sampleWeight", 1);
        if(!isPP) inTree_p->SetBranchStatus("ncollWeight", 1);
        inTree_p->SetBranchStatus("fullWeight", 1);

        inTree_p->SetBranchStatus("treePartonPt", 1);
        inTree_p->SetBranchStatus("treePartonEta", 1);
        inTree_p->SetBranchStatus("treePartonPhi", 1);
        inTree_p->SetBranchStatus("treePartonId", 1);

        inTree_p->SetBranchStatus("truthPhotonPt", 1);
        inTree_p->SetBranchStatus("truthPhotonEta", 1);
        inTree_p->SetBranchStatus("truthPhotonPhi", 1);
        inTree_p->SetBranchStatus(("truthPhotonIso"+label_phoIsoConeSize).c_str(), 1);
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
    inTree_p->SetBranchStatus(("photon_etcone"+label_phoIsoConeSize+"0").c_str(), 1);
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    inTree_p->SetBranchAddress("runNumber", &runNumber);
    inTree_p->SetBranchAddress("lumiBlock", &lumiBlock);
    if(isMC){
        inTree_p->SetBranchAddress("pthat", &pthat);
        inTree_p->SetBranchAddress("sampleWeight", &sampleWeight);
        if(!isPP) inTree_p->SetBranchAddress("ncollWeight", &ncollWeight);
        inTree_p->SetBranchAddress("fullWeight", &fullWeight);

        inTree_p->SetBranchAddress("treePartonPt", &treePartonPt);
        inTree_p->SetBranchAddress("treePartonEta", &treePartonEta);
        inTree_p->SetBranchAddress("treePartonPhi", &treePartonPhi);
        inTree_p->SetBranchAddress("treePartonId", &treePartonId);

        inTree_p->SetBranchAddress("truthPhotonPt", &truthPhotonPt);
        inTree_p->SetBranchAddress("truthPhotonEta", &truthPhotonEta);
        inTree_p->SetBranchAddress("truthPhotonPhi", &truthPhotonPhi);
        inTree_p->SetBranchAddress(("truthPhotonIso"+label_phoIsoConeSize).c_str(), &truthPhotonIso);
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
    inTree_p->SetBranchAddress(("photon_etcone"+label_phoIsoConeSize+"0").c_str(), &photon_etcone_p);


    //variable to count the number of events in a given centrality bin
    std::vector<Double_t> eventCountsPerCent;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
        eventCountsPerCent.push_back(0.0);
    }
    std::vector<int> skippedCent;

    ULong64_t nEntries_ = inTree_p->GetEntries();
    if(doGlobalDebug) nEntries_ = 2000;
    const ULong64_t nEntries = nEntries_;
    const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////
    // Event loop! 
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
        if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
        inTree_p->GetEntry(entry);

        double vert_z = vert_z_p->at(0);
        vert_z /= 10.;
        if(vert_z <= -15. || vert_z >= 15.) continue;      

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
            double photonPt = photon_pt_p->at(pI);
            if(photonPt < gammaPtBinsSub[0]) continue; // min gamma pT cut
            if(photonPt >= gammaPtBinsSub[nGammaPtBinsSub]) continue;
            Int_t ptPos = ghostPos(nGammaPtBinsSub, gammaPtBinsSub, photonPt, true, doGlobalDebug);
            if(doGlobalDebug) std::cout << "ptPos = " <<  ptPos << std::endl;

            if(truthPhotonPt<=0) continue; //prompt isolated photons?!
            if(getDR(photon_eta_p->at(pI), photon_phi_p->at(pI), truthPhotonEta, truthPhotonPhi) > phoGenMatchingDR) continue; 

            // kinematics
            int tempEtaPos = -1;
            for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
                if(abs(truthPhotonEta)>=etaBins_i[eI] && abs(truthPhotonEta)<etaBins_f[eI]) tempEtaPos=eI; 
            }
            if(tempEtaPos==-1) continue; //eta cut

            /////////////////////////////////////
            // isolation correction
            float correctedIso = photon_etcone_p->at(pI);
            if(doPtCorrectedIso && doCentCorrectedIso)
                correctedIso = getCorrectedPhotonIsolation(isPP, photon_etcone_p->at(pI), photonPt, photon_eta_p->at(pI), cent);
            else if(doPtCorrectedIso && !doCentCorrectedIso) 
                correctedIso = getPtCorrectedPhotonIsolation(photon_etcone_p->at(pI), photonPt, photon_eta_p->at(pI));

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
            /////////////////////////////////
            // fill 2D truth iso vs. reco iso (before reco isolation cut 
            if(doLooseID && photon_loose_p->at(pI)==1){
                fillTH2(photon_truthIso_vs_recoIso[centPos][tempEtaPos][ptPos],truthPhotonIso,correctedIso,fullWeight);
                fillTH2(photon_truthIso_vs_recoIso[centPos][tempEtaPos][nGammaPtBinsSub],truthPhotonIso,correctedIso,fullWeight);
                if(!isPP){
                fillTH2(photon_truthIso_vs_recoIso[nCentBins][tempEtaPos][ptPos],truthPhotonIso,correctedIso,fullWeight);
                fillTH2(photon_truthIso_vs_recoIso[nCentBins][tempEtaPos][nGammaPtBinsSub],truthPhotonIso,correctedIso,fullWeight);
                }
            } else if(!doLooseID && photon_tight_p->at(pI)==1){
                fillTH2(photon_truthIso_vs_recoIso[centPos][tempEtaPos][ptPos],truthPhotonIso,correctedIso,fullWeight);
                fillTH2(photon_truthIso_vs_recoIso[centPos][tempEtaPos][nGammaPtBinsSub],truthPhotonIso,correctedIso,fullWeight);
                if(!isPP){
                fillTH2(photon_truthIso_vs_recoIso[nCentBins][tempEtaPos][ptPos],truthPhotonIso,correctedIso,fullWeight);
                fillTH2(photon_truthIso_vs_recoIso[nCentBins][tempEtaPos][nGammaPtBinsSub],truthPhotonIso,correctedIso,fullWeight);
                }
            }
            if(truthPhotonIso>genIsoCut) continue;

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
            /////////////////////////////////
            // fill total 1D truth photon pT
            fillTH1(photonEff_TOT_CentDep_Den[centPos][tempEtaPos],truthPhotonPt,fullWeight);
            if(!isPP){
                fillTH1(photonEff_TOT_CentDep_Den[nCentBins][tempEtaPos],truthPhotonPt,fullWeight);
                fillTH2(h2F_photonEff_ISO_Den[tempEtaPos],truthPhotonPt,cent,fullWeight);
            }

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
            /////////////////////////////////
            // fill identified 1D truth photon pT
            if(doLooseID){
                if(photon_loose_p->at(pI)){
                    fillTH1(photonEff_ID_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
                    if(!isPP){
                        fillTH1(photonEff_ID_CentDep_Num[nCentBins][tempEtaPos],truthPhotonPt,fullWeight);
                    }
                }
            } else{
                if(photon_tight_p->at(pI)){
                    fillTH1(photonEff_ID_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
                    if(!isPP){
                        fillTH1(photonEff_ID_CentDep_Num[nCentBins][tempEtaPos],truthPhotonPt,fullWeight);
                    }
                }
            } 

        if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl; 
            /////////////////////////////////
            // fill isolated (or identified && isolated) 1D truth photon pT
            if(correctedIso > isoCut) continue;

            fillTH1(photonEff_ISO_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
            if(!isPP){
                fillTH1(photonEff_ISO_CentDep_Num[nCentBins][tempEtaPos],truthPhotonPt,fullWeight);
                fillTH2(h2F_photonEff_ISO_Num[tempEtaPos],truthPhotonPt,cent,fullWeight);
            }
            if(doLooseID && photon_loose_p->at(pI)){
                fillTH1(photonEff_TOT_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
                if(!isPP)
                fillTH1(photonEff_TOT_CentDep_Num[nCentBins][tempEtaPos],truthPhotonPt,fullWeight);
            } else if(!doLooseID && photon_tight_p->at(pI)){
                fillTH1(photonEff_TOT_CentDep_Num[centPos][tempEtaPos],truthPhotonPt,fullWeight);
                if(!isPP)
                fillTH1(photonEff_TOT_CentDep_Num[nCentBins][tempEtaPos],truthPhotonPt,fullWeight);
            }
        }//photon loop
    }//event loop
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    outFile_p->cd();
    //Pre-write and delete some of these require some mods
    for(Int_t cI = 0; cI < nCentBins_withIncBin; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            photonEff_ID_CentDep_Eff[cI][eI]->Divide(photonEff_ID_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Den[cI][eI],1.,1.,"B");
            photonEff_ISO_CentDep_Eff[cI][eI]->Divide(photonEff_ISO_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Den[cI][eI],1.,1.,"B");
            photonEff_TOT_CentDep_Eff[cI][eI]->Divide(photonEff_TOT_CentDep_Num[cI][eI],photonEff_TOT_CentDep_Den[cI][eI],1.,1.,"B");
            if(cI==1){
                h2F_photonEff_ISO_Eff[eI]->Divide(h2F_photonEff_ISO_Num[eI],h2F_photonEff_ISO_Den[eI],1.,1.,"B");
            }
            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                h1F_photon_recoIso[cI][eI][pI] = (TH1F*) photon_truthIso_vs_recoIso[cI][eI][pI]->ProjectionY();
                fit_photon_recoIso[cI][eI][pI] = cleverGaus(h1F_photon_recoIso[cI][eI][pI], Form("fit_%s", h1F_photon_recoIso[cI][eI][pI]->GetName()),1.0); 
            }
        }
    }

    // energy scale vs. pT
    for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
        for(Int_t cI = 0; cI < nCentBins_withIncBin; ++cI){
            for(Int_t pI = 0; pI < nGammaPtBinsSub; ++pI){
                photon_arith_meanRecoIso_vs_pt[cI][eI]->SetBinContent(pI+1,h1F_photon_recoIso[cI][eI][pI]->GetMean());
                photon_arith_meanRecoIso_vs_pt[cI][eI]->SetBinError(pI+1,h1F_photon_recoIso[cI][eI][pI]->GetMeanError());
                photon_gaus_meanRecoIso_vs_pt[cI][eI]->SetBinContent(pI+1,fit_photon_recoIso[cI][eI][pI]->GetParameter(1));
                photon_gaus_meanRecoIso_vs_pt[cI][eI]->SetBinError(pI+1,fit_photon_recoIso[cI][eI][pI]->GetParError(1));
            }
        }
    }
    // energy scale vs. cent
    for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
        for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
            for(Int_t cI = 0; cI < nCentBins_withIncBin; ++cI){
                photon_arith_meanRecoIso_vs_cent[pI][eI]->SetBinContent(cI+1,h1F_photon_recoIso[cI][eI][pI]->GetMean());
                photon_arith_meanRecoIso_vs_cent[pI][eI]->SetBinError(cI+1,h1F_photon_recoIso[cI][eI][pI]->GetMeanError());
                photon_gaus_meanRecoIso_vs_cent[pI][eI]->SetBinContent(cI+1,fit_photon_recoIso[cI][eI][pI]->GetParameter(1));
                photon_gaus_meanRecoIso_vs_cent[pI][eI]->SetBinError(cI+1,fit_photon_recoIso[cI][eI][pI]->GetParError(1));
            }
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
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    ///////////////////////////////////////////////////////////
    // Write histograms in the output file 
    for(Int_t cI = 0; cI < nCentBins_withIncBin; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            photonEff_ID_CentDep_Num[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_ID_CentDep_Eff[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_ISO_CentDep_Num[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_ISO_CentDep_Eff[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_TOT_CentDep_Num[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_TOT_CentDep_Den[cI][eI]->Write("", TObject::kOverwrite);
            photonEff_TOT_CentDep_Eff[cI][eI]->Write("", TObject::kOverwrite);
            if(cI==1){
                h2F_photonEff_ISO_Num[eI]->Write("", TObject::kOverwrite);
                h2F_photonEff_ISO_Den[eI]->Write("", TObject::kOverwrite);
                h2F_photonEff_ISO_Eff[eI]->Write("", TObject::kOverwrite);
            }
            photon_arith_meanRecoIso_vs_pt[cI][eI]->Write("", TObject::kOverwrite);
            photon_gaus_meanRecoIso_vs_pt[cI][eI]->Write("", TObject::kOverwrite);
            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                photon_truthIso_vs_recoIso[cI][eI][pI]->Write("", TObject::kOverwrite);
                h1F_photon_recoIso[cI][eI][pI]->Write("", TObject::kOverwrite);
                if(cI==0) photon_arith_meanRecoIso_vs_cent[pI][eI]->Write("", TObject::kOverwrite);
                if(cI==0) photon_gaus_meanRecoIso_vs_cent[pI][eI]->Write("", TObject::kOverwrite);
            }
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
    for(Int_t cI = 0; cI < nCentBins_withIncBin; ++cI){
        for(Int_t eI = 0; eI < nPhoEtaBins; ++eI){
            delete photonEff_ID_CentDep_Num[cI][eI];
            delete photonEff_ID_CentDep_Eff[cI][eI];
            delete photonEff_ISO_CentDep_Num[cI][eI];
            delete photonEff_ISO_CentDep_Eff[cI][eI];
            delete photonEff_TOT_CentDep_Num[cI][eI];
            delete photonEff_TOT_CentDep_Den[cI][eI];
            delete photonEff_TOT_CentDep_Eff[cI][eI];
            if(cI==1){
                delete h2F_photonEff_ISO_Num[eI];
                delete h2F_photonEff_ISO_Den[eI];
                delete h2F_photonEff_ISO_Eff[eI];
            }
            delete photon_arith_meanRecoIso_vs_pt[cI][eI];
            delete photon_gaus_meanRecoIso_vs_pt[cI][eI];
            for(Int_t pI = 0; pI < nGammaPtBinsSub+1; ++pI){
                delete photon_truthIso_vs_recoIso[cI][eI][pI];
                delete h1F_photon_recoIso[cI][eI][pI];
                if(cI==0) delete photon_arith_meanRecoIso_vs_cent[pI][eI];
                if(cI==0) delete photon_gaus_meanRecoIso_vs_cent[pI][eI];
            }
        }
    }
    if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    config_p->Write("config", TObject::kOverwrite);
    TEnv labelEnv;
    for(auto const & lab : binsToLabelStr){
        labelEnv.SetValue(lab.first.c_str(), lab.second.c_str());
    }
    labelEnv.Write("label", TObject::kOverwrite);

    outFile_p->Close();
    delete outFile_p;

    delete randGen_p;

    std::cout << "phoTaggedJetRaa_photonEfficiency COMPLETE. return 0." << std::endl;
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
