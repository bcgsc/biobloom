bin_PROGRAMS = BioBloomMaker

BioBloomMaker_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BioBloomMaker_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a

BioBloomMaker_SOURCES = BioBloomMaker.cpp \
	bloom_filter.hpp \
	BloomFilterGenerator.h BloomFilterGenerator.cpp \
	WindowedFileParser.h WindowedFileParser.cpp