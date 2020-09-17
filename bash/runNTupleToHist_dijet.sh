#!/bin/bash

export DOGLOBALDEBUGROOT=0
#2020Aug21
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPData_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbData_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPData_dijet_noMix_jetMinPt30_jetCone4_bkgPhotons.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbData_dijet_noMix_jetMinPt30_jetCone4_bkgPhotons.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt30_jetCone4_bkgPhotons.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt30_jetCone4_bkgPhotons.txt & sleep 1
#wait $(jobs -p)
./bin/gdjSignalBackgroundPhotonPlotter_dijet.exe input/comp_SignalBackgroundPhotons_PbPb_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
./bin/gdjSignalBackgroundPhotonPlotter_dijet.exe input/comp_SignalBackgroundPhotons_pp_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
./bin/gdjSignalBackgroundPhotonPlotter_dijet.exe input/comp_SignalBackgroundPhotons_PbPbMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
./bin/gdjSignalBackgroundPhotonPlotter_dijet.exe input/comp_SignalBackgroundPhotons_ppMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_SignalBackgroundPhotons_PbPbMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_SignalBackgroundPhotons_ppMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1

#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPData_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPData_dijet_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPData_dijet_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbData_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbData_dijet_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbData_dijet_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt30_jetCone4_quarkJet.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PPMC_dijet_noMix_jetMinPt30_jetCone4_gluonJet.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt30_jetCone4_quarkJet.txt & sleep 1
#./bin/gdjNTupleToHist_dijet.exe input/ntupleToHist_PbPbMC_dijet_noMix_jetMinPt30_jetCone4_gluonJet.txt & sleep 1
#wait $(jobs -p)
#./bin/gdjQuarkGluonJetsPlotter.exe input/comp_QuarkGluonJets_dijet_PbPbMC_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjQuarkGluonJetsPlotter.exe input/comp_QuarkGluonJets_dijet_PPMC_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_MC_noMix_jetMinPt40_jetCone4.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_MC_noMix_jetMinPt30_jetCone2.txt & sleep 1
#./bin/gdjDataPbPbPPPlotter.exe input/comp_PbPbpp_dijet_MC_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataMCPlotter_dijet.exe input/comp_DataMC_dijet_PP_noMix_jetMinPt30_jetCone4.txt & sleep 1
#./bin/gdjDataMCPlotter_dijet.exe input/comp_DataMC_dijet_PbPb_noMix_jetMinPt30_jetCone4.txt
