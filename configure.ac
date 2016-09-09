# switch off caching, since that leads to insurmountable trouble if
# you call ./configure subsequently on different systems without
# clearing the cache
define([AC_CACHE_LOAD], )
define([AC_CACHE_SAVE], )

AC_INIT
AC_CONFIG_HEADER([include/base/FerretConfig.h:include/base/FerretConfig.h.in])
AC_CONFIG_FILES(Makefile)

# Compute the canonical target system type variable "target", and its
# three individual parts: "target_cpu", "target_vendor", and "target_os"
AC_CANONICAL_TARGET

# The AC_SUBST(foo) command creates an output variable from a shell variable.
# This command also causes AC_OUTPUT to replace all instances of @VAR@ with
# the value of the shell variable VAR in the output files (e.g. one or more Makefiles).
AC_SUBST(host)

# Let the user specify the location of boost using
# --with-boost=/path/to/boost.  Then /path/to/boost/include will be
# appended to the CXXFLAGS variable when performing the compile test
# below.
AC_ARG_WITH(boost,
  AS_HELP_STRING([--with-boost=PATH],[Boost headers should be in PATH/include]),
                 [BOOST_INCLUDE_PATH="-I${withval}/include"],
                 [BOOST_INCLUDE_PATH=""])

# Let the user's $BOOST_DIR environment variable, if set, dictate the
# location of Boost.
if test "x$BOOST_DIR" != x; then
  BOOST_INCLUDE_PATH="-I${BOOST_DIR}/include"
fi

# We are going to test using the C++ compiler.
AC_LANG_PUSH([C++])

# Test including required boost header files.
old_CXXFLAGS="$CXXFLAGS"

# Append the Boost header location to CXXFLAGS for testing.
if test "x$BOOST_INCLUDE_PATH" != x; then
  CXXFLAGS="$CXXFLAGS $BOOST_INCLUDE_PATH"
fi

# AC_TRY_COMPILE([test program's headers],
#                [test program's body],
#                [action if success],
#                [action if fail])
AC_TRY_COMPILE(
  [@%:@include <boost/math/special_functions/legendre.hpp>],
  [],
  [AC_DEFINE(FERRET_HAVE_BOOST_MATH_SPECIAL_FUNCTIONS, 1, [Flag indicating that Boost special functions are available])
   AC_MSG_RESULT([Boost special math function header found])
   AC_SUBST(BOOST_INCLUDE_PATH)],
  [AC_MSG_RESULT([Boost special math function header not found, some Kernels may be disabled])]
)

# Reset CXXFLAGS to their original state
CXXFLAGS="$old_CXXFLAGS"
AC_LANG_POP([C++])

# Generate files and clean up after ourselves.
AC_OUTPUT
rm -f conftest* config.cache

# Local Variables:
# mode: m4
# End:
