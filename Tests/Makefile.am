bin_PROGRAMS = ReadProcessorTests \
	BloomFilterTests \
	BloomFilterCategorizerTests \
	BloomFilterMakerTests \
	WindowedParser \
	BloomFilterInfoTests \
	SeqEvalTests

ReadProcessorTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
ReadProcessorTests_SOURCES = ReadProcessorTests.cpp

BloomFilterTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
BloomFilterTests_SOURCES = BloomFilterTests.cpp

BloomFilterInfoTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp

WindowedParser_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
WindowedParser_SOURCES = WindowedParserTests.cpp
WindowedParser_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BloomFilterMakerTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
BloomFilterMakerTests_SOURCES = BloomFilterMakerTests.cpp
BloomFilterMakerTests_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BloomFilterCategorizerTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz
BloomFilterCategorizerTests_SOURCES = BloomFilterCategorizerTests.cpp
BloomFilterCategorizerTests_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

SeqEvalTests_LDADD = $(top_builddir)/Common/libcommon.a
SeqEvalTests_SOURCES = SeqEvalTests.cpp
SeqEvalTests_CPPFLAGS = -I$(top_srcdir)/Common
