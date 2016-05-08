bin_PROGRAMS = biobloommaker

biobloommaker_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloommaker_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)

biobloommaker_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz \
	-lsdsl
	
biobloommaker_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloommaker_SOURCES = BioBloomMaker.cpp \
	BloomFilterGenerator.h BloomFilterGenerator.cpp \
	WindowedFileParser.h WindowedFileParser.cpp \
	BloomMapGenerator.h BloomMapGenerator.cpp \
	SimMat.h \
	SpacedSeedIndex.h
