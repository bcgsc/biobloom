bin_PROGRAMS = biobloommaker biobloommimaker

biobloommaker_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloommaker_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)

biobloommaker_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl
	
biobloommaker_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloommaker_SOURCES = BioBloomMaker.cpp \
	BloomFilterGenerator.h BloomFilterGenerator.cpp



biobloommimaker_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloommimaker_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)

biobloommimaker_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl
	
biobloommimaker_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloommimaker_SOURCES = BioBloomMIMaker.cpp \
	MIBFGen.hpp \
	Options.cpp Options.h
