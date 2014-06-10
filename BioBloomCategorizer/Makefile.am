bin_PROGRAMS = biobloomcategorizer

biobloomcategorizer_CXXFLAGS = $(OPENMP_CXXFLAGS) $(AM_CXXFLAGS)

biobloomcategorizer_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)
	
biobloomcategorizer_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
	
biobloomcategorizer_LDFLAGS = $(OPENMP_CXXFLAGS)

biobloomcategorizer_SOURCES = BioBloomCategorizer.cpp \
	MultiFilter.h MultiFilter.cpp \
	ResultsManager.h ResultsManager.cpp \
	BioBloomClassifier.h BioBloomClassifier.cpp
