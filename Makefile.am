AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  InModuleTracks.h \
  InModuleTrack.h \
  InModuleTrackHit.h \
  InModuleTrackContainer.h

lib_LTLIBRARIES = \
  libInModuleTracks.la

libInModuleTracks_la_SOURCES = \
  InModuleTracks.cc \
  InModuleTrack.cc \
  InModuleTrackHit.cc \
  InModuleTrackContainer.cc \
  InModuleTrackDict.cc

libInModuleTracks_la_LIBADD = \
  -lphool \
  -lfun4all \
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
  -ltrackeralign \
  -ltrackbase_historic_io

BUILT_SOURCES = \
  testexternals.cc \
  InModuleTrackDict.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libInModuleTracks.la

InModuleTrackDict.cc: InModuleTrack.h InModuleTrackHit.h InModuleTrackContainer.h LinkDef.h
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