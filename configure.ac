#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([MulWick], 1.0, [shiozaki.toru@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cc])
AM_CONFIG_HEADER([config.h])
AC_F77_LIBRARY_LDFLAGS

# for the time being...
CXXFLAGS="$CXXFLAGS -std=c++0x -O0 -g"
# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

AC_CONFIG_MACRO_DIR([m4])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile prep/Makefile])
AC_OUTPUT
