check_PROGRAMS = BloomFilterTests \
	BloomFilterCategorizerTests \
	BloomFilterMakerTests \
	BloomFilterInfoTests \
	SeqEvalTests \
	ntHashTests

BloomFilterTests_LDADD = $(top_builddir)/Common/libcommon.a -lz
BloomFilterTests_SOURCES = BloomFilterTests.cpp

BloomFilterInfoTests_LDADD = $(top_builddir)/Common/libcommon.a -lz
BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp

BloomFilterMakerTests_LDADD = $(top_builddir)/Common/libcommon.a -lz
BloomFilterMakerTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
BloomFilterMakerTests_SOURCES = BloomFilterMakerTests.cpp
BloomFilterMakerTests_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BloomFilterCategorizerTests_LDADD = $(top_builddir)/Common/libcommon.a -lz
BloomFilterCategorizerTests_SOURCES = BloomFilterCategorizerTests.cpp
BloomFilterCategorizerTests_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

SeqEvalTests_LDADD = $(top_builddir)/Common/libcommon.a
SeqEvalTests_SOURCES = SeqEvalTests.cpp
SeqEvalTests_CPPFLAGS = -I$(top_srcdir)/Common

ntHashTests_LDADD = $(top_builddir)/Common/libcommon.a -lz -lsdsl
ntHashTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
ntHashTests_SOURCES = ntHashTests.cpp
ntHashTests_CPPFLAGS = -I$(top_srcdir)/Common
