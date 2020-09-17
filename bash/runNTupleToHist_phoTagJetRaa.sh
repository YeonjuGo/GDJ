#!/bin/bash

export DOGLOBALDEBUGROOT=0
#2020Sep04
./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi34_doJetPtBins.txt >& ./log/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi34_doJetPtBins.log & sleep 1
./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi34_doJetPtBins.txt>& ./log/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi34_doJetPtBins.log & sleep 1
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi34.txt>& ./log/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi34.log & sleep 1
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi34.txt>& ./log/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi34.log & sleep 1
#
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi00.txt>& ./log/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi00.log & sleep 1
##./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi00.txt>& ./log/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi00.log & sleep 1
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi00.txt>& ./log/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi00.log & sleep 1
##./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi00.txt>& ./log/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi00.log & sleep 1

##2020Sep04
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi34.txt>& ./log/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi34.log & sleep 1
##./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi34.txt>& ./log/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi34.log & sleep 1
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi34.txt>& ./log/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi34.log & sleep 1
##./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi34.txt>& ./log/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi34.log & sleep 1
#
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi00.txt>& ./log/phoTagJetRaa_PPData_noMix_tightID_Iso3_dphi00.log & sleep 1
##./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi00.txt>& ./log/phoTagJetRaa_PbPbData_noMix_tightID_Iso8_dphi00.log & sleep 1
#./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi00.txt>& ./log/phoTagJetRaa_PPMC_noMix_tightID_Iso3_dphi00.log & sleep 1
##./bin/phoTaggedJetRaa.exe input/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi00.txt>& ./log/phoTagJetRaa_PbPbMC_noMix_tightID_Iso8_dphi00.log & sleep 1
wait $(jobs -p)
#./bin/gdjQuarkGluonJetsPlotter.exe input/comp_QuarkGluonJets_dijet_PbPbMC_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjQuarkGluonJetsPlotter.exe input/comp_QuarkGluonJets_dijet_PPMC_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_noMix.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_MC_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_MC_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_MC_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataMCPlotter_dijet.exe input/comp_DataMC_dijet_PP_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataMCPlotter_dijet.exe input/comp_DataMC_dijet_PbPb_noMix_jetMinPt30_jetCone4.txt
