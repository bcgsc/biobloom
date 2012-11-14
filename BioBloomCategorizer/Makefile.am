bin_PROGRAMS = BioBloomCategorizer

BioBloomCategorizer_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BioBloomCategorizer_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

BioBloomCategorizer_SOURCES = BioBloomCategorizer.cpp \
	MultiFilter.h MultiFilter.cpp \
	BioBloomClassifier.h BioBloomClassifier.cpp