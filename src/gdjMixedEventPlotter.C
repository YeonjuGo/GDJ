//Author: Chris McGinn (2020.04.20)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

//Local
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/kirchnerPalette.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"

void plotMixClosure(const bool doGlobalDebug, std::map<std::string, std::string>* labelMap, std::map<std::string, std::string>* configMap, std::string dateStr,  std::vector<TH1F*> hists_p)
{
  std::vector<std::string> legStrs = {"Raw", "Mixed", "Raw - Mixed"};
  if(legStrs.size() != hists_p.size()) return;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  const bool isMC = std::stoi((*configMap)["ISMC"]);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string preStr = hists_p[0]->GetName();
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  while(preStr.find("/") != std::string::npos){preStr.replace(0, preStr.find("/")+1, "");}
  preStr.replace(preStr.find("_"), preStr.size(), "");

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  std::string labelStr = hists_p[0]->GetName();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  labelStr.replace(0, labelStr.find("_")+1, "");
  labelStr.replace(labelStr.rfind("_h"), 2, "");
  std::vector<std::string> preLabels;
  if(isMC) preLabels.push_back("#bf{ATLAS Internal} Monte Carlo");
  else preLabels.push_back("#bf{ATLAS Internal} Data");
  while(labelStr.find("_") != std::string::npos){
    std::string centStr = labelStr.substr(0, labelStr.find("_"));
    preLabels.push_back(centStr);
    labelStr.replace(0, labelStr.find("_")+1, "");
  }
  if(labelStr.size() != 0) preLabels.push_back(labelStr);
    
  Double_t padSplit = 0.35;
  Double_t leftMargin = 0.11;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TLine* line_p = new TLine();
  line_p->SetLineStyle(2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
  TPad* pads_p[2];
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);
  canv_p->cd();
  
  pads_p[0] = new TPad("pad0", "", 0.0, padSplit, 1.0, 1.0);
  pads_p[0]->SetTopMargin(0.01);
  pads_p[0]->SetRightMargin(0.01);
  pads_p[0]->SetLeftMargin(leftMargin);
  pads_p[0]->SetBottomMargin(0.001);

  canv_p->cd();
  pads_p[0]->Draw("SAME");
  pads_p[0]->cd();
  canv_p->cd();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, padSplit);
  pads_p[1]->SetTopMargin(0.001);
  pads_p[1]->SetRightMargin(0.01);
  pads_p[1]->SetLeftMargin(leftMargin);
  pads_p[1]->SetBottomMargin(leftMargin/padSplit);

  canv_p->cd();
  pads_p[1]->Draw("SAME");
  pads_p[1]->cd();
  canv_p->cd();

  
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  kirchnerPalette kPal;
  const int nMarkers = 4;
  int markers[nMarkers] = {24,25,28,46};
  const int nColors = 4;
  int colors[nColors] = {0,1,3,4};

  TLegend* leg_p = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(unsigned int cI = 0; cI < hists_p.size(); ++cI){
    leg_p->AddEntry(hists_p[cI], legStrs[cI].c_str(), "P L");
    
    canv_p->cd();
    pads_p[0]->cd();

    hists_p[cI]->SetMarkerSize(1);
    hists_p[cI]->SetMarkerStyle(markers[cI%nMarkers]);
    hists_p[cI]->SetMarkerColor(kPal.getColor(colors[cI%nColors]));
    hists_p[cI]->SetLineColor(kPal.getColor(colors[cI%nColors]));

    
    if(preStr.find("DPhi") != std::string::npos) hists_p[cI]->SetMaximum(0.45);
    else if(preStr.find("XJ") != std::string::npos) hists_p[cI]->SetMaximum(0.30);
    hists_p[cI]->SetMinimum(-0.04);

    hists_p[cI]->GetXaxis()->SetTitleFont(43);
    hists_p[cI]->GetYaxis()->SetTitleFont(43);
    hists_p[cI]->GetXaxis()->SetTitleSize(14);
    hists_p[cI]->GetYaxis()->SetTitleSize(14);

    hists_p[cI]->GetXaxis()->SetLabelFont(43);
    hists_p[cI]->GetYaxis()->SetLabelFont(43);
    hists_p[cI]->GetXaxis()->SetLabelSize(12);
    hists_p[cI]->GetYaxis()->SetLabelSize(12);
    
    hists_p[cI]->GetYaxis()->SetNdivisions(505);
    hists_p[cI]->GetXaxis()->SetNdivisions(505);

    hists_p[cI]->GetYaxis()->SetTitle("N_{#gamma,jet}/N_{#gamma}");

    hists_p[cI]->GetYaxis()->SetTitleOffset(1.6);
    
    if(cI == 0) hists_p[cI]->DrawCopy("HIST E1 P");
    else hists_p[cI]->DrawCopy("HIST E1 P SAME");
    gStyle->SetOptStat(0);
    
    canv_p->cd();
    pads_p[1]->cd();

    hists_p[cI]->SetMaximum(0.04);
    hists_p[cI]->SetMinimum(-0.04);

    hists_p[cI]->GetXaxis()->SetTitleOffset(2.8);
    
    if(cI == 0) hists_p[cI]->DrawCopy("HIST E1 P");
    else hists_p[cI]->DrawCopy("HIST E1 P SAME");    
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  canv_p->cd();
  pads_p[0]->cd();
  leg_p->Draw("SAME");

  labelStr = "";

  const double xPos = 0.18;
  double yPos = 0.9;

  std::string preLabelSaveStr = "";
  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    if(preLabels[pI].find("ATLAS") != std::string::npos){
      if(isMC) preLabelSaveStr = preLabelSaveStr + "MC_";
      else preLabelSaveStr = preLabelSaveStr + "DATA_";
    }
    else preLabelSaveStr = preLabelSaveStr + preLabels[pI] + "_";
  }

  for(unsigned int pI = 0; pI < preLabels.size(); ++pI){
    std::string preStr = "";
    if(preLabels[pI].find("Cent") != std::string::npos) preStr = "Pb+Pb";
    if(labelMap->count(preLabels[pI]) != 0){
      if(preStr.size() != 0) preLabels[pI] = preStr + " " + (*labelMap)[preLabels[pI]];
      else preLabels[pI] = (*labelMap)[preLabels[pI]];
	
    }

    if(labelStr.size() + preLabels[pI].size() > 40){
      if(labelStr.find(";") != std::string::npos) labelStr.replace(labelStr.rfind(";"), labelStr.size(), "");
      label_p->DrawLatex(xPos, yPos, labelStr.c_str());
      yPos -= 0.085;
      labelStr = "";
    }
    labelStr = labelStr + preLabels[pI] + "; ";
  }
  if(labelStr.find(";") != std::string::npos) labelStr.replace(labelStr.rfind(";"), labelStr.size(), "");
  
  label_p->DrawLatex(xPos, yPos, labelStr.c_str());
  
  
  canv_p->cd();
  pads_p[0]->cd();

  line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 0.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 0.0);

  canv_p->cd();
  pads_p[1]->cd();

  line_p->DrawLine(hists_p[0]->GetBinLowEdge(1), 0.0, hists_p[0]->GetBinLowEdge(hists_p[0]->GetXaxis()->GetNbins()+1), 0.0);

  canv_p->cd();
  std::string saveName = "pdfDir/" + dateStr + "/" + preStr + "_" + preLabelSaveStr + dateStr + ".pdf";
  quietSaveAs(canv_p, saveName);
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  delete leg_p;
  delete pads_p[0];
  delete pads_p[1];
  delete canv_p;

  delete line_p;
  delete label_p;
  
  return;
}

int gdjMixedEventPlotter(std::string inFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inFileName, "root")) return 1;

  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir("pdfDir/" + dateStr);
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* config_p = (TEnv*)inFile_p->Get("config");
  TEnv* label_p = (TEnv*)inFile_p->Get("label");
  configParser labels(label_p);
  configParser configs(config_p);
  std::map<std::string, std::string> labelMap = labels.GetConfigMap();
  std::map<std::string, std::string> configMap = configs.GetConfigMap();
  
  std::vector<std::string> necessaryParams = {"ISMC",
					      "NGAMMAPTBINSSUB",
					      "ISPP"};

  std::vector<std::string> pbpbParams = {"CENTBINS"};
  
  if(!configs.ContainsParamSet(necessaryParams)) return 1;

  const int nGammaPtBinsSub = std::stoi(configs.GetConfigVal("NGAMMAPTBINSSUB"));
  const bool isPP = std::stoi(configs.GetConfigVal("ISPP"));
  
  int nCentBins = 1;
  std::vector<std::string> centBins = {"PP"};
  if(!isPP){
    if(!configs.ContainsParamSet(pbpbParams)) return 1;    
    
    centBins = strToVect(configs.GetConfigVal("CENTBINS"));
    nCentBins = centBins.size()-1;
  }

  std::vector<std::string> observables1 = {"JtDPhi", "JtXJ", "JtPt", "JtMult", "MultiJtPt", "MultiJtXJ"};
  std::vector<std::string> backStr1 = {"_GlobalJtPt0", "_GlobalJtPt0_DPhi0", "_DPhi0", "_GlobalJtPt0_DPhi0", "_DPhi0_MultiJt0", "_GlobalJtPt0_DPhi0_MultiJt0"};
  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "PP";
    if(!isPP) centStr = "Cent" + centBins[cI] + "to" + centBins[cI+1];
    
    for(Int_t gI = 0; gI < nGammaPtBinsSub+1; ++gI){
      for(unsigned int oI = 0; oI < observables1.size(); ++oI){     
	std::string rawName = centStr + "/photon" + observables1[oI] + "VCentPt_" + centStr + "_GammaPt" + std::to_string(gI) + backStr1[oI] + "_h";
	
	TH1F* raw_p = (TH1F*)inFile_p->Get(rawName.c_str());
	rawName.replace(rawName.find("photon"), std::string("photon").size(), "photonMix");
	TH1F* mix_p = (TH1F*)inFile_p->Get(rawName.c_str());
	rawName.replace(rawName.find("photonMix"), std::string("photonMix").size(), "photonSub");
	TH1F* sub_p = (TH1F*)inFile_p->Get(rawName.c_str());

	if(doGlobalDebug) std::cout << "OBSERVABLE: " << observables1[oI] << std::endl;
	plotMixClosure(doGlobalDebug, &labelMap, &configMap, dateStr, {raw_p, mix_p, sub_p});
      }
    }
  }

  inFile_p->Close();
  delete inFile_p;

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "GDJMIXEDEVENTPLOTTER COMPLETE. return 0" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjMixedEventPlotter.exe <inFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjMixedEventPlotter(argv[1]);
  return retVal;
}