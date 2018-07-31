bin_PROGRAMS = biobloomcategorizer biobloommicategorizer

biobloomcategorizer_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloomcategorizer_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)
	
biobloomcategorizer_LDADD = $(top_builddir)/Common/libcommon.a \
	-lsdsl
	
biobloomcategorizer_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloomcategorizer_SOURCES = BioBloomCategorizer.cpp \
	ResultsManager.hpp \
	BioBloomClassifier.h BioBloomClassifier.cpp \
	MIBFClassifier.hpp \
	Options.h Options.cpp


biobloommicategorizer_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloommicategorizer_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)
	
biobloommicategorizer_LDADD = $(top_builddir)/Common/libcommon.a \
	-lsdsl
	
biobloommicategorizer_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloommicategorizer_SOURCES = BioBloomMICategorizer.cpp \
	ResultsManager.hpp \
	MIBFClassifier.hpp \
	Options.h Options.cpp
