# Shell script

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

# This file should be sourced by configure, and sets the following
# environment variables corresponding to the recommended settings for a
# given OS/CPU/compiler combination:
#
# cppflags_default       # Base CPPFLAGS                     (default: "")

# cflags_default         # Base CFLAGS                       (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging    (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization (default: "-O")
# cflags_default_hot     # Optimization for specific files   (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling    (default: "-g")
# cflags_default_omp     # Added to $CFLAGS for OpenMP       (default: "")

# cxxflags_default       # Base CXXFLAGS                       (default: "")
# cxxflags_default_dbg   # Added to $CXXFLAGS for debugging    (default: "-g")
# cxxflags_default_opt   # Added to $CXXFLAGS for optimization (default: "-O")
# cxxflags_default_hot   # Optimization for specific files     (default: "-O")
# cxxflags_default_prf   # Added to $CXXFLAGS for profiling    (default: "-g")
# cxxflags_default_omp   # Added to $CXXFLAGS for OpenMP       (default: "")
# cxxflags_default_std   # C++ standard variant                (default: "")

# fcflags_default        # Base FCFLAGS                       (default: "")
# fcflags_default_dbg    # Added to $FCFLAGS for debugging    (default: "-g")
# fcflags_default_opt    # Added to $FCFLAGS for optimization (default: "-O")
# fcflags_default_hot    # Optimization for specific files    (default: "-O")
# fcflags_default_prf    # Added to $FCFLAGS for profiling    (default: "-g")
# fcflags_default_omp    # Added to $FCFLAGS for OpenMP       (default: "")
#
# ldflags_default        # Base LDFLAGS                       (default: "")
# ldflags_default_dbg    # Added to $LDFLAGS for debugging    (default: "-g")
# ldflags_default_opt    # Added to $LDFLAGS for optimization (default: "-O")
# ldflags_default_prf    # Added to $LDFLAGS for profiling    (default: "-g")
# ldflags_rpath          # Added to $LDFLAGS for shared libs  (default: "")

# libs_default           # Base LIBS                          (default: "")
# libs_default_dbg       # Added to $LIBS for debugging       (default: "")
# libs_default_opt       # Added to $LIBS for optimization    (default: "")
# libs_default_prf       # Added to $LIBS for profiling       (default: "")

# Two other environment variable strings are defined, containing possibly
# more detailed compiler information:
#
# cs_ac_cc_version      # Compiler version string, 1 line max.
# cs_ac_cc_version_full # Compiler version string, 10 lines max.

# The sourcing approach and some tests are borrowed from the HDF5 configure
# environment.
#
# We choose to source this script rather than use a more classical m4 macro
# for this functionality, so that a user may more easily modify
# default compiler options or port to a new machine without requiring
# any advanced knowledge of autoconf or m4 macros, or installation of
# an autoconf environment on the target machine.

# Initialize local variables
#---------------------------

outfile=cs_ac_env-tmp

# Some compilers require a file to compile even for version info.

cat > conftest.c <<\_______EOF
int main()
{
  return 0;
}
_______EOF

# Compiler version info may be localization dependent (especially for gcc)

save_LANG=$LANG
unset LANG;

# Libraries only added on special cases (such as IBM Blue Gene),
# so initialize empty variables here

libs_default=""
libs_default_dbg=""
libs_default_opt=""
libs_default_prf=""

##################
#                #
#  Preprocessor  #
#                #
##################


# Default pre-processor flags (not too dependent on compiler)
#----------------------------

case "$host_os" in
  *)
    cppflags_default=""
    ;;
esac


################
#              #
#  C compiler  #
#              #
################

cs_ac_cc_version=unknown
cs_cc_compiler_known=no

cflags_default_prf="-g"

# Are we using gcc ?
#-------------------

cs_gcc=no

if test "x$GCC" = "xyes"; then

  # Intel and Pathscale compilers may pass as GCC but
  # may be recognized by version string

  if test -n "`$CC --version | grep icc`" ; then
    cs_gcc=icc
  elif test -n "`$CC --version | grep clang`" ; then
    cs_gcc=clang
  elif test -n "`$CC --version 2>&1 | grep PathScale`" ; then
    cs_gcc=pathcc
  else
    cs_gcc=gcc
  fi

fi

if test "x$cs_gcc" = "xgcc"; then

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  cs_ac_cc_version=`$CC --version 2>&1 | head -1`
  cs_cc_compiler_known=yes

  # Practical version info for option setting
  cs_cc_version="`$CC -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"
  cs_cc_vendor=`echo $cs_cc_version |sed 's/\([a-z]*\).*/\1/'`
  cs_cc_version=`echo $cs_cc_version |sed 's/[-a-z]//g'`

  if test "x" = "x$cs_cc_vendor" -a "x" != "x$cs_cc_version"; then
    cs_cc_vendor=gcc
  fi
  if test "-" != "$cs_cc_vendor-$cs_cc_version"; then
    echo "compiler '$CC' is GNU $cs_cc_vendor-$cs_cc_version"
  fi

  # Some version numbers
  cs_cc_vers_major=`echo $cc_version | cut -f1 -d.`
  cs_cc_vers_minor=`echo $cc_version | cut -f2 -d.`
  cs_cc_vers_patch=`echo $cc_version | cut -f3 -d.`
  test -n "$cs_cc_vers_major" || cs_cc_vers_major=0
  test -n "$cs_cc_vers_minor" || cs_cc_vers_minor=0
  test -n "$cs_cc_vers_patch" || cs_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-std=c99 -fms-extensions -funsigned-char -pedantic -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wfloat-equal"
  cflags_default_dbg="-g"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_omp="-fopenmp"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in

    *i?86|*x86_64)
      cflags_default_opt="-funroll-loops -O2 -Wuninitialized"
      case "$host_cpu" in
        i686)
          cflags_default_opt="$cflags_default_opt -march=i686"
          ;;
      esac
      ;;

  esac

  # Modify default flags depending on gcc version (as older versions
  # may not handle all flags)

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-[45]*)
      ;;
    *)
      cflags_default="$cflags_default -Wmisleading-indentation -Wduplicated-cond"
      ;;
  esac

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-4.[012345678]*)
      ;;
    *)
      cflags_default="$cflags_default -fdiagnostics-color=auto"
      ;;
  esac

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-4.[012]*)
      cflags_default_omp=""
      ;;
  esac

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-4.[01234]*)
      ;;
    *)
      cflags_default_opt="$cflags_default_opt -fexcess-precision=fast"
      cflags_default_hot="$cflags_default_hot -fexcess-precision=fast"
      ;;
  esac

  case "$host_os" in
    mingw32)
      cflags_default="`echo $cflags_default | sed -e 's/-std=c99/-std=gnu99/g'`"
      cflags_default="$cflags_default -Wno-format -Wno-pedantic-ms-format"
      ;;
  esac

# Otherwise, are we using icc ?
#------------------------------

elif test "x$cs_gcc" = "xicc"; then

  cs_cc_version=`echo $CC --version | grep icc |sed 's/[a-zA-Z()]//g'`

  echo "compiler '$CC' is Intel ICC"

  # Version strings for logging purposes and known compiler flag
  $CC -V conftest.c > $outfile 2>&1
  cs_ac_cc_version=`$CC --version 2>&1 | head -1`
  cs_cc_compiler_known=yes

  # Some version numbers
  cs_cc_vers_major=`echo $cs_ac_cc_version | cut -f 3 -d" " | cut -f1 -d.`
  cs_cc_vers_minor=`echo $cs_ac_cc_version | cut -f 3 -d" " | cut -f2 -d.`
  cs_cc_vers_patch=`echo $cs_ac_cc_version | cut -f 3 -d" " | cut -f3 -d.`
  test -n "$cs_cc_vers_major" || cs_cc_vers_major=0
  test -n "$cs_cc_vers_minor" || cs_cc_vers_minor=0
  test -n "$cs_cc_vers_patch" || cs_cc_vers_patch=0

  # Default compiler flags
  # (temporarily disable "operands evaluated in unspecified order" remark -- 981)
  cflags_default="-std=c99 -restrict -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -wd981"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_omp="-qopenmp"
  case "$cs_cc_vers_major" in
    1[0123456])
      cflags_default_omp="-openmp"
      ;;
  esac

# Otherwise, are we using clang ?
#--------------------------------

elif test "x$cs_gcc" = "xclang"; then

  cs_cc_version=`echo $CC --version | grep clang | cut -f 3 -d ' '`

  echo "compiler '$CC' is clang"

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  cs_ac_cc_version=`$CC --version 2>&1 | head -1`
  cs_cc_compiler_known=yes

  # Default compiler flags
  # (temporarily disable "operands evaluated in unspecified order" remark -- 981)
  cflags_default="-std=c99 -funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_omp="-fopenmp=libomp"

# Otherwise, are we using pathcc ?
#---------------------------------

elif test "x$cs_gcc" = "xpathcc"; then

  $CC --version 2>&1 | grep 'PathScale' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is PathScale C compiler"

    # Version strings for logging purposes and known compiler flag
    $CC --version > $outfile 2>&1
    cs_ac_cc_version=`grep -i Compiler $outfile`
    cs_cc_compiler_known=yes

    # Default compiler flags
    cflags_default="-c99 -noswitcherror"
    cflags_default="-std=c99 -funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wunused-value"
    cflags_default_dbg="-g"
    cflags_default_opt="-O2"
    cflags_default_hot="-Ofast"
    cflags_default_omp="-openmp"

  fi

fi

# Otherwise, are we using pgcc ?
#-------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Portland Group pgcc"

    # Version strings for logging purposes and known compiler flag
    $CC -V > $outfile 2>&1
    cs_ac_cc_version=`grep -i pgcc $outfile`
    cs_cc_compiler_known=yes

    # Default compiler flags
    cflags_default="-c99 -noswitcherror"
    cflags_default_dbg="-g -Mbounds"
    cflags_default_opt="-O2"
    cflags_default_hot="-fast"
    cflags_default_omp="-mp"

  fi

fi

# Otherwise, are we using xlc ?
#------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -qversion 2>&1 | grep 'XL C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is IBM XL C compiler"

    # Version strings for logging purposes and known compiler flag
    $CC -qversion > $outfile 2>&1
    cs_ac_cc_version=`grep 'XL C' $outfile`
    cs_cc_compiler_known=yes
    cs_linker_set=yes

    # Default compiler flags
    cflags_default="-qlanglvl=stdc99 -q64"
    cflags_default_opt="-O3"
    cflags_default_hot="-O3"
    cflags_default_dbg="-g -qfullpath"
    cflags_default_omp="-qsmp=omp -qthreaded"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O3"
    ldflags_default_dbg="-g -qfullpath"
    ldflags_rpath="-R"

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null
    if test "$?" = "0" ; then
      # Default compiler flags (we assume that MPI wrappers are used)
      if test -f /bgsys/drivers/ppcfloor/cnk/bin/bgq_kernel.elf ; then
        cs_ibm_bg_type="Q"
        cppflags_default=""
        cflags_default=""                    # "-qlanglvl=extc99" by default
        cflags_default_opt="-g -O3"
        cflags_default_hot="-g -O3 -qhot"
        cflags_default_dbg="-g"
      fi
    fi

  fi
fi

# Otherwise, are we using the Cray compiler ?
#------------------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'Cray C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Cray C compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cc_version=`$CC -V 2>&1 | grep "Cray C" | head -1`
    cs_cc_compiler_known=yes
    cs_linker_set=yes

    # Default compiler flags
    cflags_default=""                        # "-h c99" by default
    cflags_default_opt="-O2"
    cflags_default_hot="-O3"
    cflags_default_dbg="-g"
    cflags_default_omp="-h omp"              # default: use "-h noomp" to disable

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O2"
    ldflags_default_dbg="-g"

  fi
fi

# Compiler still not identified
#------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  case "$host_os" in

    *)

      # Generic
      #--------

      cflags_default=""
      cflags_default_opt="-O"
      cflags_default_dbg="-g"
      cflags_default_omp=""
      ;;

  esac

fi

if test "x$cflags_default_hot" = "x" ; then
  cflags_default_hot=$cflags_default_opt
fi

if test -f $outfile ; then
  cs_ac_cc_version_full=`sed -e '11,$d' $outfile`
  rm -f $outfile
fi


##################
#                #
#  C++ compiler  #
#                #
##################

cs_ac_cxx_version=unknown
cs_cxx_compiler_known=no

cxxflags_default_prf="-g"

# Are we using g++ ?
#-------------------

cs_gxx=no

if test "x$GXX" = "xyes"; then

  # Intel and Pathscale compilers may pass as GXX but
  # may be recognized by version string

  if test -n "`$CXX --version | grep icpc`" ; then
    cs_gxx=icpc
  elif test -n "`$CXX --version | grep icc`" ; then
    cs_gxx=icc
  elif test -n "`$CXX --version | grep clang`" ; then
    cs_gxx=clang
  elif test -n "`$CXX --version 2>&1 | grep PathScale`" ; then
    cs_gxx=pathCC
  else
    cs_gxx=g++
  fi

fi

if test "x$cs_gxx" = "xg++"; then

  # Version strings for logging purposes and known compiler flag
  $CXX -v > $outfile 2>&1
  cs_ac_cxx_version=`$CXX --version 2>&1 | head -1`
  cs_cxx_compiler_known=yes

  # Practical version info for option setting
  cs_cxx_version="`$CXX -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"
  cs_cxx_vendor=`echo $cxx_version |sed 's/\([a-z]*\).*/\1/'`
  cs_cxx_version=`echo $cs_cxx_version |sed 's/[-a-z]//g'`

  if test "x" = "x$cs_cxx_vendor" -a "x" != "x$cs_cxx_version"; then
    cs_cxx_vendor=g++
  fi
  if test "-" != "$cs_cxx_vendor-$cs_cxx_version"; then
    echo "compiler '$CXX' is GNU $cs_cxx_vendor-$cs_cxx_version"
  fi

  # Some version numbers
  cs_cxx_vers_major=`echo $cs_cxx_version | cut -f1 -d.`
  cs_cxx_vers_minor=`echo $cs_cxx_version | cut -f2 -d.`
  cs_cxx_vers_patch=`echo $cs_cxx_version | cut -f3 -d.`
  test -n "$cs_cxx_vers_major" || cs_cxx_vers_major=0
  test -n "$cs_cxx_vers_minor" || cs_cxx_vers_minor=0
  test -n "$cs_cxx_vers_patch" || cs_cxx_vers_patch=0

  # Default compiler flags
  cxxflags_default="-W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wunused -Wfloat-equal"
  cxxflags_default_dbg="-g"
  cxxflags_default_opt="-O2"
  cxxflags_default_hot="-O3"
  cxxflags_default_omp="-fopenmp"
  cxxflags_default_std="-funsigned-char"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in

    *i?86|*x86_64)
      cxxflags_default_opt="-funroll-loops -O2 -Wuninitialized"
      case "$host_cpu" in
        i686)
          cxxflags_default_opt="$cxxflags_default_opt -march=i686"
          ;;
      esac
      ;;

  esac

  # Modify default flags depending on g++ version (as older versions
  # may not handle all flags)

  case "$cs_cxx_vendor-$cs_cxx_version" in
    g++-[45]*)
      ;;
    *)
      cxxflags_default="$cxxflags_default -Wmisleading-indentation -Wduplicated-cond"
      ;;
  esac

  case "$cs_cxx_vendor-$cs_cxx_version" in
    g++-4.[012345678]*)
      ;;
    *)
      cxxflags_default="$cxxflags_default -fdiagnostics-color=auto"
      ;;
  esac

  case "$cs_cxx_vendor-$cs_cxx_version" in
    g++-2.*|g++-3*|g++-4.[012]*)
      cxxflags_default_omp=""
      ;;
  esac

# Otherwise, are we using icc ?
#------------------------------

elif test "x$cs_gxx" = "xicpc" -o "x$cs_gxx" = "xicc"; then

  if test "x$cs_gxx" = "xicpc"; then
    cs_cxx_version=`echo $CXX --version | grep icpc |sed 's/[a-zA-Z()]//g'`
  else
    cs_cxx_version=`echo $CXX --version | grep icc |sed 's/[a-zA-Z()]//g'`
  fi

  echo "compiler '$CXX' is Intel ICC"

  # Version strings for logging purposes and known compiler flag
  $CXX -V conftest.c > $outfile 2>&1
  cs_ac_cxx_version=`$CXX --version 2>&1 | head -1`
  cs_cxx_compiler_known=yes

  cs_cxx_vers_major=`echo $cs_ac_cxx_version | cut -f 3 -d" " | cut -f1 -d.`
  cs_cxx_vers_minor=`echo $cs_ac_cxx_version | cut -f 3 -d" " | cut -f2 -d.`
  cs_cxx_vers_patch=`echo $cs_ac_cxx_version | cut -f 3 -d" " | cut -f3 -d.`
  test -n "$cs_cxx_vers_major" || cs_cxx_vers_major=0
  test -n "$cs_cxx_vers_minor" || cs_cxx_vers_minor=0
  test -n "$cs_cxx_vers_patch" || cs_cxx_vers_patch=0

  # Default compiler flags
  cxxflags_default="-std=c++11 -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cxxflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cxxflags_default_opt="-O2"
  cxxflags_default_hot="-O3"
  cxxflags_default_omp="-qopenmp"
  cxxflags_default_std="-funsigned-char"

  case "$cs_cxx_vers_major" in
    1[0123456])
      cxxflags_default_omp="-openmp"
      ;;
  esac

# Otherwise, are we using clang ?
#--------------------------------

elif test "x$cs_gxx" = "xclang"; then

  cs_cc_version=`echo $CXX --version | grep clang | cut -f 3 -d ' '`

  echo "compiler '$CXX' is clang"

  # Version strings for logging purposes and known compiler flag
  $CXX -v > $outfile 2>&1
  cs_ac_cxx_version=`$CXX --version 2>&1 | head -1`
  cs_cxx_compiler_known=yes

  # Default compiler flags
  # (temporarily disable "operands evaluated in unspecified order" remark -- 981)
  cxxlags_default="-Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cxxflags_default_dbg="-g -O0"
  cxxflags_default_opt="-O2"
  cxxflags_default_hot="-O3"
  cxxflags_default_omp="-fopenmp=libomp"
  cxxflags_default_std="-funsigned-char"

# Otherwise, are we using pgcc ?
#-------------------------------

else

  $CXX -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CXX' is Portland Group pgCC"

    # Version strings for logging purposes and known compiler flag
    $CXX -V conftest.c > $outfile 2>&1
    cs_ac_cxx_version=`grep -i pgcc $outfile`
    cs_cxx_compiler_known=yes

    # Default compiler flags
    cxxflags_default="-noswitcherror"
    cxxflags_default_dbg="-g -Mbounds"
    cxxflags_default_opt="-fast -fastsse"
    cxxflags_default_omp="-mp"
    cxxflags_default_std="-Xa"

  fi

fi

# Otherwise, are we using xlc++ ?
#--------------------------------

if test "x$cs_cxx_compiler_known" != "xyes" ; then

  $CXX -qversion 2>&1 | grep 'XL C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CXX' is IBM XL C/C++ compiler"

    # Version strings for logging purposes and known compiler flag
    $CXX -qversion > $outfile 2>&1
    cs_ac_cxx_version=`grep 'XL C' $outfile`
    cs_cxx_compiler_known=yes

    # Default compiler flags
    cxxflags_default="-q64 -qlanglvl=redefmac"
    cxxflags_default_opt="-O3"
    cxxflags_default_hot="-O3"
    cxxflags_default_dbg="-g"
    cxxflags_default_omp="-qsmp=omp -qthreaded"
    cxxflags_default_std=""

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null
    if test "$?" = "0" ; then
      # Default compiler flags (we assume that MPI wrappers are used)
      if test -f /bgsys/drivers/ppcfloor/cnk/bin/bgq_kernel.elf ; then
        cxxflags_default="-qlanglvl=redefmac"
        cxxflags_default_opt="-g -O3"
        cxxflags_default_hot="-g -O3 -qhot"
        cxxflags_default_dbg="-g"
      fi
    fi

  fi
fi

# Otherwise, are we using the Cray compiler ?
#------------------------------------------

if test "x$cs_cxx_compiler_known" != "xyes" ; then

  $CXX -V 2>&1 | grep 'Cray C++' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CXX' is Cray C++"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cxx_version=`$CXX -V 2>&1 | grep "Cray C++" | head -1`
    cs_cxx_compiler_known=yes

    # Default compiler flags
    cxxflags_default=""                        # "-h c99" by default
    cxxflags_default_opt="-O2"
    cxxflags_default_hot="-O3"
    cxxflags_default_dbg="-g"
    cfxxlags_default_omp="-h omp"              # default: use "-h noomp" to disable
    cxxflags_default_std=""

  fi

fi

# Otherwise, are we using pathcc ?
#---------------------------------

if test "x$cs_cxx_compiler_known" != "xyes" ; then

  $CXX --version 2>&1 | grep 'PathScale' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CXX' is PathScale C++ compiler"

    # Version strings for logging purposes and known compiler flag
    $CXX --version > $outfile 2>&1
    cs_ac_cxx_version=`grep -i Compiler $outfile`
    cs_cxx_compiler_known=yes

    # Default compiler flags
    cxxflags_default="-W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wunused-value"
    cxxflags_default_dbg="-g"
    cxxflags_default_opt="-O2"
    cxxflags_default_hot="-Ofast"
    cxxflags_default_omp="-openmp"
    cxxflags_default_std=""

  fi

fi

# Compiler still not identified
#------------------------------

if test "x$cs_cxx_compiler_known" != "xyes" ; then

  case "$host_os" in

    *)

      # Generic
      #--------

      cxxflags_default=""
      cxxflags_default_opt="-O"
      cxxflags_default_dbg="-g"
      cxxflags_default_omp=""
      cxxflags_default_std=""
      ;;

  esac

fi

if test "x$cxxflags_default_hot" = "x" ; then
  cxxflags_default_hot=$cxxflags_default_opt
fi

if test -f $outfile ; then
  cs_ac_cxx_version_full=`sed -e '11,$d' $outfile`
  rm -f $outfile
fi


######################
#                    #
#  Fortran compiler  #
#                    #
######################

cs_ac_fc_version=unknown
cs_fc_compiler_known=no

fcflags_default_prf="-g"

# Are we using gfortran ?
#------------------------

cs_gfortran=no

# Are we using gfortran ?
#------------------------

$FC --version 2>&1 | grep 'GNU Fortran' > /dev/null

if test "$?" = "0" ; then

  cs_fc_version=`echo $FC --version  |sed 's/[a-zA-Z()]//g'`
  cs_fc_version="`$FC -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"

  cs_fc_compiler_known=yes
  cs_gfortran=gfortran

  # Version strings for logging purposes and known compiler flag
  $FC -v > $outfile 2>&1
  cs_ac_fc_version=`$FC --version 2>&1 | head -1`
  cs_fc_compiler_known=yes

  # Practical version info for option setting
  cs_fc_version="`$FC -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"
  cs_fc_vendor=`echo $fc_version |sed 's/\([a-z]*\).*/\1/'`
  cs_fc_version=`echo $cs_fc_version |sed 's/[-a-z]//g'`

  if test "x" = "x$cs_fc_vendor" -a "x" != "x$cs_fc_version"; then
    cs_fc_vendor=gfortran
  fi
  if test "-" != "$cs_fc_vendor-$cs_fc_version"; then
    echo "compiler '$FC' is GNU $cs_fc_vendor-$cs_fc_version"
  else
    echo "compiler '$FC' is gfortran"
  fi

  # Some version numbers
  cs_fc_vers_major=`echo $cc_version | cut -f1 -d.`
  cs_fc_vers_minor=`echo $cc_version | cut -f2 -d.`
  cs_fc_vers_patch=`echo $cc_version | cut -f3 -d.`
  test -n "$cs_fc_vers_major" || cs_fc_vers_major=0
  test -n "$cs_fc_vers_minor" || cs_fc_vers_minor=0
  test -n "$cs_fc_vers_patch" || cs_fc_vers_patch=0

  # Default compiler flags
  fcflags_default="-x f95-cpp-input -Wall -pedantic-errors -std=f2003"
  fcflags_default_dbg="-g -fcheck=bounds"
  fcflags_default_opt="-O"
  fcflags_default_hot="-O2"
  fcflags_default_omp="-fopenmp"

  if test "xgfortran" = "x$cs_fc_vendor"; then
    case "$cs_fc_version" in
      4.[234]*)
        fcflags_default_dbg="`echo $fcflags_default_dbg | sed -e 's/-fcheck=bounds/-fbounds-check/g'`"
        ;;
      4.[56789]*)
        ;;
      *)
        fcflags_default="$fcflags_default -fdiagnostics-color=auto"
       ;;
    esac
  fi

fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  # Are we using ifort ?
  #---------------------

  $FC --version 2>&1 | grep 'IFORT' > /dev/null
  if test "$?" = "0" ; then

    cs_fc_version=`echo $FC --version | grep ifort |sed 's/[a-zA-Z()]//g'`

    echo "compiler '$FC' is Intel Fortran"

    # Version strings for logging purposes and known compiler flag
    $FC -V > $outfile 2>&1
    cs_ac_fc_version=`$FC --version 2>&1 | head -1`
    cs_fc_compiler_known=yes

    cs_fc_vers_major=`echo $cs_ac_fc_version | cut -f 3 -d" " | cut -f1 -d.`
    cs_fc_vers_minor=`echo $cs_ac_fc_version | cut -f 3 -d" " | cut -f2 -d.`
    cs_fc_vers_patch=`echo $cs_ac_fc_version | cut -f 3 -d" " | cut -f3 -d.`
    test -n "$cs_fc_vers_major" || cs_fc_vers_major=0
    test -n "$cs_fc_vers_minor" || cs_fc_vers_minor=0
    test -n "$cs_fc_vers_patch" || cs_fc_vers_patch=0

    # Default compiler flags
    # (temporarily disable "unused variable" remark -- 7712)
    fcflags_default="-cpp -fpic -warn -diag-disable 7712"
    fcflags_default_dbg="-g -O0 -traceback -check all -check nopointer -fpe0 -ftrapuv"
    fcflags_default_opt="-O2"
    fcflags_default_hot="-O3"
    fcflags_default_omp="-qopenmp"

    case "$cs_cxx_vers_major" in
      1[0123456])
        fcflags_default_omp="-openmp"
        ;;
    esac

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  # Are we using pgf95 ?
  #---------------------

  $FC -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$FC' is Portland Group Fortran compiler"

    # Version strings for logging purposes and known compiler flag
    $FC -V > $outfile 2>&1
    cs_ac_cc_version=`grep -i pgf $outfile`
    cs_fc_compiler_known=yes

    # Default compiler flags
    fcflags_default="-Mpreprocess -noswitcherror"
    fcflags_default_dbg="-g -Mbounds"
    fcflags_default_opt="-O2"
    fcflags_default_hot="-fast"
    fcflags_default_omp="-mp"

  fi

fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  # Are we using xlf ?
  #-------------------

  $FC -qversion 2>&1 | grep 'XL Fortran' > /dev/null

  if test "$?" = "0" ; then

    echo "compiler '$FC' is IBM XL Fortran compiler"

    # Version strings for logging purposes and known compiler flag
    $FC -qversion > $outfile 2>&1
    cs_ac_fc_version=`grep 'XL Fortran' $outfile`
    cs_fc_compiler_known=yes

    fcflags_default="-q64 -qextname -qsuffix=cpp=f90"
    fcflags_default_dbg="-g"
    fcflags_default_opt="-O3"
    fcflags_default_omp="-qsmp=omp -qthreaded"

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null

    if test "$?" = "0" ; then

      # Default compiler flags (we assume that MPI wrappers are used)
      if test -f /bgsys/drivers/ppcfloor/cnk/bin/bgq_kernel.elf ; then
        fcflags_default="-qextname -qsuffix=cpp=f90"
        fcflags_default_dbg="-g -qcheck"
        fcflags_default_opt="-g -O3"
        fcflags_default_hot="-g -O3 -qhot"
      fi
    fi

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  # Are we using the Cray compiler ?
  #-------------------------------

  $FC -V 2>&1 | grep 'Cray Fortran' > /dev/null

  if test "$?" = "0" ; then

    echo "compiler '$FC' is Cray Fortran compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_fc_version=`$FC -V 2>&1 | grep "Cray Fortran" | head -1`
    cs_fc_compiler_known=yes

    fcflags_default="-eF -em -J."
    fcflags_default_dbg="-g"
    fcflags_default_opt="-O2"
    fcflags_default_omp="-h omp"              # default: use "-h noomp" to disable

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  # Are we using pathf95 ?
  #---------------------

  $FC --version 2>&1 | grep 'PathScale' > /dev/null

  if test "$?" = "0" ; then

    echo "compiler '$FC' is PathScale Fortran compiler"

    # Version strings for logging purposes and known compiler flag
    $FC --version > $outfile 2>&1
    cs_ac_fc_version=`grep 'PathScale' $outfile`
    cs_fc_compiler_known=yes

    fcflags_default="-Wall -Wno-unused -cpp"
    fcflags_default_dbg="-g -ffortran-bounds-check"
    fcflags_default_opt="-O"
    fcflags_default_hot="-fast"
    fcflags_default_omp="-openmp"

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  case "$host_os" in

    *)

      # Generic
      #--------

      fcflags_default=""
      fcflags_default_opt="-O"
      fcflags_default_hot="-O"
      fcflags_default_dbg="-g"
      ;;

  esac

fi

if test "x$fcflags_default_hot" = "x" ; then
  fcflags_default_hot=$fcflags_default_opt
fi

if test -f $outfile ; then
  cs_ac_fc_version_full=`sed -e '11,$d' $outfile`
fi


############
#          #
#  Linker  #
#          #
############


# Default linker flags
#---------------------

ldflags_default_prf="-g"

if test "x$cs_linker_set" != "xyes" ; then

  case "$host_os" in

    linux*)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_rpath="-Wl,-rpath -Wl,"
      if test "x$cs_gfortran" = "xgfortran"; then
        libgfortran_path=`$FC --print-file-name=libgfortran.so`
        libgfortran_dir=`dirname $libgfortran_path`
        unset libgfortran_path
        if test "`echo $libgfortran_dir | cut -c1-4`" != "/usr" ; then
          ldflags_rpath="${ldflags_rpath}${libgfortran_dir}"
        fi
        libgfortran_path=`$FC --print-file-name=libgfortran.so`
        libgfortran_dir=`dirname $libgfortran_path`
        unset libgfortran_path
        if test "`echo $libgfortran_dir | cut -c1-4`" != "/usr" ; then
          ldflags_rpath="${ldflags_rpath}${libgfortran_dir}"
        fi
        unset libgfortran_dir
      fi
      ;;

    *)
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ;;

  esac

fi

# Finish

export LANG=$save_LANG

# Clean temporary files

rm -f conftest* a.out $outfile

