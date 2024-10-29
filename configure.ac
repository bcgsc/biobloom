#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(BIOBLOOMTOOLS, 2.3.4, cjustin@bcgsc.ca)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([BioBloomMaker/BioBloomMaker.cpp])
AC_CONFIG_HEADER([config.h])

#see http://stackoverflow.com/questions/8215785/how-can-i-use-ac-revision-with-git
AC_REVISION([m4_esyscmd([./tools/configure.commit])])
AC_DEFINE([GIT_REVISION], "[m4_esyscmd([./tools/configure.commit])]", [git revision information])

AC_PROG_RANLIB

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AC_CHECK_LIB([dl], [dlopen])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Set compiler flags. (TODO Probably doing some bad stuff here)
AC_SUBST(CXXFLAGS,"-std=c++17 $CXXFLAGS")
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Werror $AM_CXXFLAGS")

# Options to configure.
# Boost
AC_ARG_WITH(boost, AS_HELP_STRING([--with-boost=PATH],
			[specify directory for the boost header files]))
if test "$with_boost" -a -d "$with_boost"; then
	boost_cppflags="-isystem$with_boost"
fi

AC_ARG_WITH(sparsehash, AS_HELP_STRING([--with-sparsehash=PATH],
			[specify directory for the google sparsehash header files]))
if test "$with_sparsehash" -a -d "$with_sparsehash"; then
	sparsehash_cppflags="-isystem$with_sparsehash"
fi

AC_ARG_WITH(sdsl, AS_HELP_STRING([--with-sdsl=PATH],
			[specify directory for the sdsl header files]))
if test "$with_sdsl" -a -d "$with_sdsl"; then
	sdsl_cppflags="-isystem$with_sdsl/include"
	sdsl_ldflags="-L$with_sdsl/lib"
fi

# Find the absolute paths.
srcdir=$(cd $srcdir; pwd)

# Set compiler flags.
AC_SUBST(CPPFLAGS,"-I$srcdir $boost_cppflags $sparsehash_cppflags $sdsl_cppflags $CPPFLAGS")
AC_SUBST(LDFLAGS,"$sdsl_ldflags $LDFLAGS")

# Checks for header files.
AC_LANG([C++])	#sets the auto conf to look for c++ headers
AC_CHECK_HEADERS([boost/property_tree/ptree.hpp])
AC_CHECK_HEADERS([boost/property_tree/ini_parser.hpp])
AC_CHECK_HEADERS([btl_bloomfilter/BloomFilter.hpp])
AC_CHECK_HEADERS([zlib.h])

AC_CHECK_HEADERS([google/dense_hash_map])
if test $ac_cv_header_google_dense_hash_map != yes; then
	AC_MSG_ERROR([Requires google sparsehash])
fi

AC_CHECK_HEADERS([sdsl/bit_vector_il.hpp])
if test $ac_cv_header_sdsl_bit_vector_il_hpp != yes; then
	AC_MSG_ERROR([Requires sdsl])
fi

# Check for Boost.
if test $ac_cv_header_boost_property_tree_ptree_hpp != yes; then
	AC_MSG_ERROR([Requires Boost])
fi

# Check for OpenMP.
AC_OPENMP
if test -z $OPENMP_CXXFLAGS; then
	OPENMP_CXXFLAGS=-Wno-unknown-pragmas
fi

# Checks for typedefs, structures, and compiler characteristics.

AC_CONFIG_FILES([
		Makefile
		Common/Makefile
		BioBloomMaker/Makefile
		BioBloomCategorizer/Makefile
		Tests/Makefile
])
AC_OUTPUT

