AUTOMAKE_OPTIONS = foreign
AM_CXXFLAGS = -pthread
AM_CPPFLAGS = \
  -I$(top_srcdir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include
AM_LDFLAGS = \
  -pthread \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64
pkginclude_HEADERS = \
  InModuleTracks.h \
  InModuleTrack.h \
  InModuleTrackv1.h \
  InModuleTrackContainer.h \
  InModuleTrackContainerv1.h \
  FullTrack.h \
  FullTrackv1.h \
  FullTrackContainer.h \
  FullTrackContainerv1.h \
  FullTrackVertex.h \
  FullTrackVertexv1.h \
  FullTrackVertexContainer.h \
  FullTrackVertexContainerv1.h \
  FullTrackConnector.h \
  FullTrackVertexer.h \
  FinalTrack.h \
  FinalTrackv1.h \
  FinalTrackContainer.h \
  FinalTrackContainerv1.h \
  FinalTrackVertex.h \
  FinalTrackVertexv1.h \
  FinalTrackVertexContainer.h \
  FinalTrackVertexContainerv1.h \
  FinalTrackVertexer.h \
  TpcPolyTrack.h \
  TpcPolyTrackv1.h \
  TpcPolyTrackContainer.h \
  TpcPolyTrackContainerv1.h \
  TpcPolyTrackReco.h \
  TpcPolyCluster.h \
  TpcPolyClusterv1.h \
  TpcPolyClusterContainer.h \
  TpcPolyClusterContainerv1.h \
  TpcPolyClusterizer.h \
  TpcPolyClusterTrack.h \
  TpcPolyClusterTrackv1.h \
  TpcPolyClusterTrackContainer.h \
  TpcPolyClusterTrackContainerv1.h \
  TpcPolyHelixFitter.h \
  Fitter.h \
  IdealPadMap.h \
  TpcPadMap.h \
  TpcPadMapv1.h \
  TpcPadMapBuilder.h
lib_LTLIBRARIES = \
  libInModuleTracks.la
libInModuleTracks_la_SOURCES = \
  InModuleTracks.cc \
  InModuleTrack.cc \
  InModuleTrackv1.cc \
  InModuleTrackContainer.cc \
  InModuleTrackContainerv1.cc \
  FullTrack.cc \
  FullTrackv1.cc \
  FullTrackContainer.cc \
  FullTrackContainerv1.cc \
  FullTrackVertex.cc \
  FullTrackVertexv1.cc \
  FullTrackVertexContainer.cc \
  FullTrackVertexContainerv1.cc \
  FullTrackConnector.cc \
  FullTrackVertexer.cc \
  FinalTrack.cc \
  FinalTrackv1.cc \
  FinalTrackContainer.cc \
  FinalTrackContainerv1.cc \
  FinalTrackVertex.cc \
  FinalTrackVertexv1.cc \
  FinalTrackVertexContainer.cc \
  FinalTrackVertexContainerv1.cc \
  FinalTrackVertexer.cc \
  TpcPolyTrack.cc \
  TpcPolyTrackv1.cc \
  TpcPolyTrackContainer.cc \
  TpcPolyTrackContainerv1.cc \
  TpcPolyTrackReco.cc \
  TpcPolyHelixFitter.cc \
  TpcPolyCluster.cc \
  TpcPolyClusterv1.cc \
  TpcPolyClusterContainer.cc \
  TpcPolyClusterContainerv1.cc \
  TpcPolyClusterizer.cc \
  TpcPolyClusterTrack.cc \
  TpcPolyClusterTrackv1.cc \
  TpcPolyClusterTrackContainer.cc \
  TpcPolyClusterTrackContainerv1.cc \
  Fitter.cc \
  IdealPadMap.cc \
  TpcPadMap.cc \
  TpcPadMapv1.cc \
  TpcPadMapBuilder.cc \
  InModuleTrackDict.cc
libInModuleTracks_la_LIBADD = \
  -lphool \
  -lfun4all \
  -lffamodules \
  -lcdbobjects \
  -lg4eval \
  -lg4eval_io \
  -lphg4hit \
  -lmbd_io \
  -lg4detectors_io \
  -lmvtx \
  -lmvtx_io \
  -lintt \
  -lintt_io \
  -lmicromegas_io \
  -lSubsysReco \
  -ltpc_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lgenfit2 \
  -lgenfit2exp \
  -lPHGenFit \
  -lphfield \
  -lPHGarfield
BUILT_SOURCES = \
  testexternals.cc \
  InModuleTrackDict.cc
noinst_PROGRAMS = \
  testexternals
testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libInModuleTracks.la
InModuleTrackDict.cc: InModuleTrack.h InModuleTrackv1.h InModuleTrackContainer.h InModuleTrackContainerv1.h FullTrack.h FullTrackv1.h FullTrackContainer.h FullTrackContainerv1.h FullTrackVertex.h FullTrackVertexv1.h FullTrackVertexContainer.h FullTrackVertexContainerv1.h FinalTrack.h FinalTrackv1.h FinalTrackContainer.h FinalTrackContainerv1.h FinalTrackVertex.h FinalTrackVertexv1.h FinalTrackVertexContainer.h FinalTrackVertexContainerv1.h TpcPolyTrack.h TpcPolyTrackv1.h TpcPolyTrackContainer.h TpcPolyTrackContainerv1.h TpcPolyCluster.h TpcPolyClusterv1.h TpcPolyClusterContainer.h TpcPolyClusterContainerv1.h TpcPolyClusterTrack.h TpcPolyClusterTrackv1.h TpcPolyClusterTrackContainer.h TpcPolyClusterTrackContainerv1.h Fitter.h IdealPadMap.h TpcPadMap.h TpcPadMapv1.h LinkDef.h
	rootcling -f $@ -c $(AM_CPPFLAGS) $^
testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@
clean-local:
	rm -f $(BUILT_SOURCES)
	rm -f InModuleTrackDict_rdict.pcm
