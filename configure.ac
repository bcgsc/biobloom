#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(BIOBLOOMTOOLS, 3.0.0b, cjustin@bcgsc.ca)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([BioBloomMaker/BioBloomMaker.cpp])
AC_CONFIG_HEADER([config.h])

AC_PROG_RANLIB

# Checks for programs.
# make sure -O3 by default for speed
if test -z $CXXFLAGS; then
    CXXFLAGS='-O3'
fi
AC_PROG_CXX

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Options to configure.
# Boost
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
			[specify directory for the boost header files]))
if test "$with_boost" -a -d "$with_boost"; then
	boost_cppflags="-isystem$with_boost"
fi

# Find the absolute paths.
srcdir=$(cd $srcdir; pwd)

# Set compiler flags.
boost_ver=1.58.0
boost_ver_dir=1_58_0
AC_SUBST(CPPFLAGS,"-I$srcdir $boost_cppflags $CPPFLAGS -isystem$srcdir/$boost_ver_dir")

# make sure -std=c++0x is used for sdsl
AC_SUBST(CPPFLAGS,"-std=c++0x $CPPFLAGS")

# Check for Boost.
if test $ac_cv_header_boost_unordered_unordered_map_hpp != yes; then
	AC_MSG_ERROR([Requires the Boost C++ libraries, which may
	be downloaded from here: http://www.boost.org/users/download/
	The following commands will download and install Boost:
	cd $srcdir
	wget http://downloads.sourceforge.net/project/boost/boost/$boost_ver/$boost_ver_dir.tar.bz2
	tar jxf $boost_ver_dir.tar.bz2
	cd -])
fi

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers
AC_CHECK_HEADERS([boost/unordered/unordered_map.hpp])
AC_CHECK_HEADERS([boost/property_tree/ptree.hpp])
AC_CHECK_HEADERS([boost/property_tree/ini_parser.hpp])
AC_CHECK_HEADERS([boost/tuple/tuple.hpp])
AC_CHECK_HEADERS([zlib.h])

# TODO: CHECK FOR SDSL Headers

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, '-Wall -Wextra -Werror')

# Checks for typedefs, structures, and compiler characteristics.

AC_CONFIG_FILES([
		Makefile
		Common/Makefile
		DataLayer/Makefile
		BioBloomMaker/Makefile
		BioBloomCategorizer/Makefile
		bloomfilter/Makefile
        Tests/Makefile
])
AC_OUTPUT

