bin_PROGRAMS = biobloomcategorizer

biobloomcategorizer_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer
	
biobloomcategorizer_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

biobloomcategorizer_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz

biobloomcategorizer_SOURCES = BioBloomCategorizer.cpp \
	MultiFilter.h MultiFilter.cpp \
	ResultsManager.h ResultsManager.cpp \
	BioBloomClassifier.h BioBloomClassifier.cpp