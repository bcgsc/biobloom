bin_PROGRAMS = bbr

bbr_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

bbr_CPPFLAGS = -I$(top_srcdir)/BioBloomRecruit \
	-I$(top_srcdir)/Common \
	-I$(top_srcdir)/DataLayer \
	-I$(top_srcdir)

bbr_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Common/libcommon.a -lz \
	$(top_builddir)/bloomfilter/libbloomfilter.a \
	-lsdsl
	
bbr_LDFLAGS = $(OPENMP_CXXFLAGS)

bbr_SOURCES = BioBloomRecruit.cpp \
	MIBloomTag.cpp MIBloomTag.h \
	MIRecurit.cpp MIRecurit.h
	
