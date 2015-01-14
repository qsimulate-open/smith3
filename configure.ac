#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([SMITH], 3.0, [shiozaki@northwestern.edu], [], [http://www.nubakery.org])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cc])
AM_CONFIG_HEADER([config.h])

# Does not take any time anyways
CXXFLAGS="$CXXFLAGS -std=c++11 -O0 -g"
# Checks for programs.
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])

AC_CONFIG_FILES([Makefile prep/Makefile])
AC_OUTPUT
