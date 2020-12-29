AC_INIT([dEploid], esyscmd([sh -c "cat VERSION | tr -d '\n'"]),[https://github.com/mcveanlab/DEploid/issues])

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
AC_CHECK_HEADERS([zlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

## Enable Link-time optimization if supported (gcc only)
#if test x$CXX = xg++; then
  #AX_CHECK_COMPILE_FLAG([-flto], [OPT_CXXFLAGS="$OPT_CXXFLAGS -flto"], [], [-Werror])
  #AC_SUBST(OPT_CXXFLAGS)
#fi

# bash complete
AC_ARG_WITH([bash-completion-dir],
    AS_HELP_STRING([--with-bash-completion-dir[=PATH]],
        [Install the bash auto-completion script in this directory. @<:@default=yes@:>@]),
    [],
    [with_bash_completion_dir=yes])

if test "x$with_bash_completion_dir" = "xyes"; then
    PKG_CHECK_MODULES([BASH_COMPLETION], [bash-completion >= 2.0],
        [BASH_COMPLETION_DIR="`pkg-config --variable=completionsdir bash-completion`"],
        [BASH_COMPLETION_DIR="$datadir/bash-completion/completions"])
else
    BASH_COMPLETION_DIR="$with_bash_completion_dir"
fi

AC_SUBST([BASH_COMPLETION_DIR])
AM_CONDITIONAL([ENABLE_BASH_COMPLETION],[test "x$with_bash_completion_dir" != "xno"])


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


AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
   then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AM_CONDITIONAL([HAVE_DOXYGEN],
[test -n "$DOXYGEN"])AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([docs/doxygen/Doxyfile])])


AC_CONFIG_FILES([Makefile docs/doxygen/Makefile])
AC_OUTPUT
