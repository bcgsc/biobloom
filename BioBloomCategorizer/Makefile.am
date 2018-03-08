bin_PROGRAMS = biobloomcategorizer

biobloomcategorizer_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloomcategorizer_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)
	
biobloomcategorizer_LDADD = $(top_builddir)/Common/libcommon.a -lz \
	-lsdsl
	
biobloomcategorizer_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloomcategorizer_SOURCES = BioBloomCategorizer.cpp \
	ResultsManager.hpp \
	BioBloomClassifier.h BioBloomClassifier.cpp \
	MIBFClassifier.hpp \
	Options.h Options.cpp
