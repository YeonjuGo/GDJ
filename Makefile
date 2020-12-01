CXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

define GDJDIRERR
 GDJDIR is not set at all. Please set this environment variable to point to your build - this should be either
export GDJDIR=$(PWD)
or
source setEnv.sh
if you have made appropriate changes.
For more, see README for full setup recommendations
endef

ifndef GDJDIR
$(error "$(GDJDIRERR)")	
endif

INCLUDE=-I$(GDJDIR)
LIB=-L$(GDJDIR)/lib
ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(GDJDIR)/bin
MKDIR_LIB=mkdir -p $(GDJDIR)/lib
MKDIR_OBJ=mkdir -p $(GDJDIR)/obj
MKDIR_OUTPUT=mkdir -p $(GDJDIR)/output
MKDIR_PDF=mkdir -p $(GDJDIR)/pdfDir

all: mkdirBin mkdirLib mkdirObj mkdirOutput mkdirPdf obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o lib/libATLASGDJ.so bin/gdjNtuplePreProc.exe bin/gdjNTupleToHist.exe bin/gdjNTupleToHist_dijet.exe bin/gdjHistDumper.exe bin/gdjGammaJetResponsePlot.exe bin/gdjMixedEventPlotter.exe bin/gdjQuarkGluonJetsPlotter.exe bin/gdjResponsePlotter.exe bin/gdjDataMCPlotter_dijet.exe bin/gdjSignalBackgroundPhotonPlotter_dijet.exe bin/gdjDataPbPbPPPlotter.exe bin/gdjDataMCRawPlotter.exe bin/grlToTex.exe bin/testKeyHandler.exe bin/testSampleHandler.exe bin/phoTaggedJetRaa_photonEff.exe bin/phoTaggedJetRaa_photonEnergy.exe bin/phoTaggedJetRaa_photonPurity.exe bin/phoTaggedJetRaa.exe bin/phoTaggedJetRaa_jetPt_photonEffPurCorrected.exe bin/phoTaggedJetRaa_jetPt.exe bin/phoTaggedJetRaa_jetEnergy.exe bin/gdjNtuplePreProc_onlyPhoton.exe bin/gdjNtuplePreProc_phoTaggedJetRaa.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirLib:
	$(MKDIR_LIB)

mkdirObj:
	$(MKDIR_OBJ)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

obj/centralityFromInput.o: src/centralityFromInput.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/centralityFromInput.C -o obj/centralityFromInput.o $(INCLUDE) $(ROOT)

obj/checkMakeDir.o: src/checkMakeDir.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/checkMakeDir.C -o obj/checkMakeDir.o $(INCLUDE)

obj/configParser.o: src/configParser.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/configParser.C -o obj/configParser.o $(INCLUDE) $(ROOT)

obj/globalDebugHandler.o: src/globalDebugHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/globalDebugHandler.C -o obj/globalDebugHandler.o $(ROOT) $(INCLUDE)

obj/keyHandler.o: src/keyHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/keyHandler.C -o obj/keyHandler.o $(INCLUDE)

obj/sampleHandler.o: src/sampleHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/sampleHandler.C -o obj/sampleHandler.o $(ROOT) $(INCLUDE)

lib/libATLASGDJ.so:
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libATLASGDJ.so obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o $(ROOT) $(INCLUDE)

bin/gdjNtuplePreProc.exe: src/gdjNtuplePreProc.C
	$(CXX) $(CXXFLAGS) src/gdjNtuplePreProc.C -o bin/gdjNtuplePreProc.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNtuplePreProc_onlyPhoton.exe: src/gdjNtuplePreProc_onlyPhoton.C
	$(CXX) $(CXXFLAGS) src/gdjNtuplePreProc_onlyPhoton.C -o bin/gdjNtuplePreProc_onlyPhoton.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNtuplePreProc_phoTaggedJetRaa.exe: src/gdjNtuplePreProc_phoTaggedJetRaa.C
	$(CXX) $(CXXFLAGS) src/gdjNtuplePreProc_phoTaggedJetRaa.C -o bin/gdjNtuplePreProc_phoTaggedJetRaa.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNTupleToHist.exe: src/gdjNTupleToHist.C
	$(CXX) $(CXXFLAGS) src/gdjNTupleToHist.C -o bin/gdjNTupleToHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNTupleToHist_dijet.exe: src/gdjNTupleToHist_dijet.C
	$(CXX) $(CXXFLAGS) src/gdjNTupleToHist_dijet.C -o bin/gdjNTupleToHist_dijet.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHistDumper.exe: src/gdjHistDumper.C
	$(CXX) $(CXXFLAGS) src/gdjHistDumper.C -o bin/gdjHistDumper.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjGammaJetResponsePlot.exe: src/gdjGammaJetResponsePlot.C
	$(CXX) $(CXXFLAGS) src/gdjGammaJetResponsePlot.C -o bin/gdjGammaJetResponsePlot.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjMixedEventPlotter.exe: src/gdjMixedEventPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjMixedEventPlotter.C -o bin/gdjMixedEventPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjResponsePlotter.exe: src/gdjResponsePlotter.C
	$(CXX) $(CXXFLAGS) src/gdjResponsePlotter.C -o bin/gdjResponsePlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjSignalBackgroundPhotonPlotter_dijet.exe: src/gdjSignalBackgroundPhotonPlotter_dijet.C
	$(CXX) $(CXXFLAGS) src/gdjSignalBackgroundPhotonPlotter_dijet.C -o bin/gdjSignalBackgroundPhotonPlotter_dijet.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjQuarkGluonJetsPlotter.exe: src/gdjQuarkGluonJetsPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjQuarkGluonJetsPlotter.C -o bin/gdjQuarkGluonJetsPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjDataMCPlotter_dijet.exe: src/gdjDataMCPlotter_dijet.C
	$(CXX) $(CXXFLAGS) src/gdjDataMCPlotter_dijet.C -o bin/gdjDataMCPlotter_dijet.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjDataPbPbPPPlotter.exe: src/gdjDataPbPbPPPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjDataPbPbPPPlotter.C -o bin/gdjDataPbPbPPPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjDataMCRawPlotter.exe: src/gdjDataMCRawPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjDataMCRawPlotter.C -o bin/gdjDataMCRawPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/grlToTex.exe: src/grlToTex.C
	$(CXX) $(CXXFLAGS) src/grlToTex.C -o bin/grlToTex.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/testKeyHandler.exe: src/testKeyHandler.C
	$(CXX) $(CXXFLAGS) src/testKeyHandler.C -o bin/testKeyHandler.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/testSampleHandler.exe: src/testSampleHandler.C
	$(CXX) $(CXXFLAGS) src/testSampleHandler.C -o bin/testSampleHandler.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa.exe: src/phoTaggedJetRaa.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa.C -o bin/phoTaggedJetRaa.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa_photonEff.exe: src/phoTaggedJetRaa_photonEff.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa_photonEff.C -o bin/phoTaggedJetRaa_photonEff.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa_photonPurity.exe: src/phoTaggedJetRaa_photonPurity.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa_photonPurity.C -o bin/phoTaggedJetRaa_photonPurity.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa_photonEnergy.exe: src/phoTaggedJetRaa_photonEnergy.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa_photonEnergy.C -o bin/phoTaggedJetRaa_photonEnergy.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa_jetEnergy.exe: src/phoTaggedJetRaa_jetEnergy.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa_jetEnergy.C -o bin/phoTaggedJetRaa_jetEnergy.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa_jetPt.exe: src/phoTaggedJetRaa_jetPt.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa_jetPt.C -o bin/phoTaggedJetRaa_jetPt.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/phoTaggedJetRaa_jetPt_photonEffPurCorrected.exe: src/phoTaggedJetRaa_jetPt_photonEffPurCorrected.C
	$(CXX) $(CXXFLAGS) src/phoTaggedJetRaa_jetPt_photonEffPurCorrected.C -o bin/phoTaggedJetRaa_jetPt_photonEffPurCorrected.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ


clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bash/*~
	rm -f bash/#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f include/*~
	rm -f include/#*#
	rm -f input/*~
	rm -f input/#*#
	rm -f lib/*.so
	rm -rf lib
	rm -f obj/*.o
	rm -rf obj
	rm -f src/*~
	rm -f src/#*#
