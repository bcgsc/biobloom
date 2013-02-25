bin_PROGRAMS = BloomFilterTests \
	BloomFilterCategorizerTests \
	BloomFilterMakerTests \
	WindowedParser \
	BloomFilterInfoTests \
	ReadProcessorTests

ReadProcessorTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
ReadProcessorTests_SOURCES = ReadProcessorTests.cpp

BloomFilterTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
BloomFilterTests_SOURCES = BloomFilterTests.cpp

BloomFilterInfoTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp

WindowedParser_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
WindowedParser_SOURCES = WindowedParserTests.cpp
WindowedParser_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BloomFilterMakerTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
BloomFilterMakerTests_SOURCES = BloomFilterMakerTests.cpp
BloomFilterMakerTests_CPPFLAGS = -I$(top_srcdir)/BioBloomMaker \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer

BloomFilterCategorizerTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a
BloomFilterCategorizerTests_SOURCES = BloomFilterCategorizerTests.cpp
BloomFilterCategorizerTests_CPPFLAGS = -I$(top_srcdir)/BioBloomCategorizer \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer