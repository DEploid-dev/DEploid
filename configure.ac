AC_INIT([dEploid], esyscmd([sh -c "cat version | tr -d '\n'"]),[https://github.com/mcveanlab/DEploid/issues])

#AM_INIT_AUTOMAKE([subdir-objects -Wall foreign])
AM_INIT_AUTOMAKE([1.9 tar-ustar subdir-objects -Wall foreign])

# Use -O3 as default optimization level
: ${CXXFLAGS="-O3"}

# Suppress output at compilation
m4_ifndef([AM_SILENT_RULES], [m4_define([AM_SILENT_RULES],[])])
AM_SILENT_RULES([yes])

# Load macros in 'm4'-dir
#AC_CONFIG_MACRO_DIR([m4])

AC_PROG_INSTALL
AC_PREREQ
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_RANLIB

# Check for C++11
#AX_CXX_COMPILE_STDCXX_11(,mandatory)

# Checks for libraries
AC_CHECK_LIB(cppunit,TestCase,[])

# Checks for header files.
AC_HEADER_STDC
AC_LANG(C++)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

## Enable Link-time optimization if supported (gcc only)
#if test x$CXX = xg++; then
  #AX_CHECK_COMPILE_FLAG([-flto], [OPT_CXXFLAGS="$OPT_CXXFLAGS -flto"], [], [-Werror])
  #AC_SUBST(OPT_CXXFLAGS)
#fi

# Enable Warnings
AX_CXXFLAGS_WARN_ALL
#AX_CHECK_COMPILE_FLAG([-Wextra], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Wextra"], [], [-Werror])
## AX_CHECK_COMPILE_FLAG([-Wshadow], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Wshadow"], [], [-Werror])
#AX_CHECK_COMPILE_FLAG([-Wnon-virtual-dtor], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Wnon-virtual-dtor"], [], [-Werror])
##AX_CHECK_COMPILE_FLAG([-Wold-style-cast], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Wold-style-cast"], [], [-Werror])
#AX_CHECK_COMPILE_FLAG([-Wcast-align], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Wcast-align"], [], [-Werror])
#AX_CHECK_COMPILE_FLAG([-Wunused], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Wunused"], [], [-Werror])
#AX_CHECK_COMPILE_FLAG([-Woverloaded-virtual], [OPT_CXXFLAGS="$OPT_CXXFLAGS -Woverloaded-virtual"], [], [-Werror])
#AX_CHECK_COMPILE_FLAG([-pedantic], [OPT_CXXFLAGS="$OPT_CXXFLAGS -pedantic"], [], [-Werror])

m4_include([m4/ac_doxygen.m4])
DX_INIT_DOXYGEN(PfDeconv, Doxyfile, html/)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
