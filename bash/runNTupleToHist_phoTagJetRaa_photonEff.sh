#!/bin/bash

export DOGLOBALDEBUGROOT=0
#2020Sep04
#./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PbPbMC_tightID_Iso8.txt >& ./log/phoTagJetRaa_photonEff_PbPbMC_tightID_Iso8.log & sleep 1
#./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PbPbMC_tightID_Iso5.txt >& ./log/phoTagJetRaa_photonEff_PbPbMC_tightID_Iso5.log & sleep 1
#./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PbPbMC_tightID_Iso3.txt >& ./log/phoTagJetRaa_photonEff_PbPbMC_tightID_Iso3.log & sleep 1
#./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PbPbMC_looseID_Iso8.txt >& ./log/phoTagJetRaa_photonEff_PbPbMC_looseID_Iso8.log & sleep 1
#./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PPMC_tightID_Iso3.txt >& ./log/phoTagJetRaa_photonEff_PPMC_tightID_Iso3.log & sleep 1
#./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PPMC_looseID_Iso3.txt >& ./log/phoTagJetRaa_photonEff_PPMC_looseID_Iso3.log & sleep 1
./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PPMC_tightID_Iso5.txt >& ./log/phoTagJetRaa_photonEff_PPMC_tightID_Iso5.log & sleep 1
./bin/phoTaggedJetRaa_photonEff.exe input/phoTagJetRaa_photonEff_PPMC_tightID_Iso8.txt >& ./log/phoTagJetRaa_photonEff_PPMC_tightID_Iso8.log & sleep 1
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
