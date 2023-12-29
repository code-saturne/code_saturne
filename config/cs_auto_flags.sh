# Shell script

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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
# cppflags_default       # Base CPPFLAGS                       (default: "")

# cflags_default         # Base CFLAGS                         (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging      (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization   (default: "-O")
# cflags_default_hot     # Optimization for specific files     (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling      (default: "-g")
# cflags_default_omp     # Added to $CFLAGS for OpenMP         (default: "")
# cflags_default_omp_ad  # Added to $CFLAGS for OpenMP offload (default: "")
# cflags_default_shared  # Added to $CFLAGS for shared libs    (default: "-fPIC -DPIC")

# cxxflags_default       # Base CXXFLAGS                       (default: "")
# cxxflags_default_dbg   # Added to $CXXFLAGS for debugging    (default: "-g")
# cxxflags_default_opt   # Added to $CXXFLAGS for optimization (default: "-O")
# cxxflags_default_hot   # Optimization for specific files     (default: "-O")
# cxxflags_default_prf   # Added to $CXXFLAGS for profiling    (default: "-g")
# cxxflags_default_omp   # Added to $CXXFLAGS for OpenMP       (default: "")
# cxxflags_default_omp_ad  # Added to $CXXFLAGS for OpenMP offload (default: "")
# cxxflags_default_std   # C++ standard variant                (default: "")
# cxxflags_default_shared  # Added to $CXXFLAGS for shared libs (default: "-fPIC -DPIC")

# fcflags_default        # Base FCFLAGS                       (default: "")
# fcflags_default_dbg    # Added to $FCFLAGS for debugging    (default: "-g")
# fcflags_default_opt    # Added to $FCFLAGS for optimization (default: "-O")
# fcflags_default_hot    # Optimization for specific files    (default: "-O")
# fcflags_default_prf    # Added to $FCFLAGS for profiling    (default: "-g")
# fcflags_default_omp    # Added to $FCFLAGS for OpenMP       (default: "")
# fclags_default_shared  # Added to $FCFLAGS for shared libs  (default: "-fPIC -DPIC")
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

# cs_ac_cc_version       # C compiler version string, 1 line max.
# cs_ac_cc_version_full  # C compiler version string, 10 lines max.
# cs_ac_cxx_version      # C++ compiler version string, 1 line max.
# cs_ac_cxx_version_full # C++ compiler version string, 10 lines max.
# cs_ac_fc_version       # Fortran compiler version string, 1 line max.
# cs_ac_fc_version_full  # Fortran compiler version string, 10 lines max.

# The sourcing approach and some tests are borrowed from the HDF5 configure
# environment.

# Remarks:
#--------

# We choose to source this script rather than use a more classical m4 macro
# for this functionality, so that a user may more easily modify
# default compiler options or port to a new machine without requiring
# any advanced knowledge of autoconf or m4 macros, or installation of
# an autoconf environment on the target machine.

# Some associated settings relative to library creation are also done in
# `build-aux/cs_link_library.py`, which may also need to be adapted for
# ports to some machines, using different linker options.

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

# Libraries only added on special cases),
# so initialize empty variables here

libs_default=""
libs_default_dbg=""
libs_default_opt=""
libs_default_prf=""

############################
#                          #
#  Shared library options  #
#                          #
############################

# Options to generate position-independent code for shared libraries.
# can be modified later if necessary, but usually common to most compilers.

case "$host_os" in
  darwin*)
    cflags_default_shared="-fPIC -DPIC"
    fcflags_default_shared="-fPIC -DPIC"
    cxxflags_default_shared="-fPIC -DPIC"
    ldflags_default_shared="-dynamiclib -undefined dynamic_lookup"
    ldflags_default_soname="-install_name @rpath/"
    ;;
  *)
    cflags_default_shared="-fPIC -DPIC"
    fcflags_default_shared="-fPIC -DPIC"
    cxxflags_default_shared="-fPIC -DPIC"
    ldflags_default_shared="-shared"
    ldflags_default_soname="-Wl,-soname -Wl,"
    ;;
esac

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

  # Intel and LLVM compilers may pass as GCC but
  # may be recognized by version string
  # NVIDIA compiler defines __GNUC__ which confuses configure

  cs_ac_cc_version=`$CC $user_CFLAGS --version 2>&1 | head -1`

  if test -n "`echo $cs_ac_cc_version | grep ICC`" ; then
    cs_gcc=icc
  elif test -n "`echo $cs_ac_cc_version | grep -e DPC++ -e oneAPI`" ; then
    cs_gcc=oneAPI
  elif test -n "`echo $cs_ac_cc_version | grep clang`" ; then
    cs_gcc=clang
  elif test -n "`echo $cs_ac_cc_version | grep Cray | grep -v GCC`" ; then
    cs_gcc=cray
  elif test -n "`echo $cs_ac_cc_version | grep FCC`" ; then
    cs_gcc=fujitsu
  elif test -n "`echo $cs_ac_cc_version | grep Arm`" ; then
    cs_gcc=arm
  elif test -n "`$CC $user_CFLAGS --version 2>&1 | grep NVIDIA`" ; then
    cs_gcc=no
  else
    cs_gcc=gcc
  fi

fi

if test "x$cs_gcc" = "xgcc"; then

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
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
  cflags_default="-funsigned-char -W -Wall -Werror=shadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Werror=missing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wfloat-equal -Werror=implicit-function-declaration"
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
    gcc-[4]*)
      cflags_default="$cflags_default -std=c11"
      ;;
    gcc-[5]*)
      ;;
    *)
      cflags_default="$cflags_default -Wmisleading-indentation -Wduplicated-cond"
      ;;
  esac

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-4.[012345678]*)
      ;;
    *)
      cflags_default="$cflags_default -fdiagnostics-color=auto -Werror=format-security"
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

# Otherwise, are we using ICC Classic ?
#--------------------------------------

elif test "x$cs_gcc" = "xicc" ; then

  cs_cc_version=`echo $cs_ac_cc_version | grep ICC |sed 's/[a-zA-Z()]//g'`
  echo "compiler '$CC' is Intel ICC Classic"

  # Version strings for logging purposes and known compiler flag
  $CC $user_CFLAGS -V conftest.c > $outfile 2>&1
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
  cflags_default="-std=c11 -restrict -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -wd981"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_omp="-qopenmp"
  case "$cs_cc_vers_major" in
    1[0123456])
      cflags_default_omp="-openmp"
      ;;
  esac

# Otherwise, are we using Intel LLVM DPC++/C++ Compiler ?
#--------------------------------------------------------

elif test "x$cs_gcc" = "xoneAPI" ; then

  cs_cc_version=`echo $cs_ac_cc_version | grep -e DPC++ -e oneAPI |sed 's/[a-zA-Z()+/]//g'`
  echo "compiler '$CC' is oneAPI DPC++/C++ Compiler"

  # Version strings for logging purposes and known compiler flag
  $CC $user_CFLAGS -V conftest.c > $outfile 2>&1
  cs_cc_compiler_known=yes

  # Some version numbers
  cs_cc_vers_major=`echo $cs_ac_cc_version | cut -f 5 -d" " | cut -f1 -d.`
  cs_cc_vers_minor=`echo $cs_ac_cc_version | cut -f 5 -d" " | cut -f2 -d.`
  cs_cc_vers_patch=`echo $cs_ac_cc_version | cut -f 5 -d" " | cut -f3 -d.`
  test -n "$cs_cc_vers_major" || cs_cc_vers_major=0
  test -n "$cs_cc_vers_minor" || cs_cc_vers_minor=0
  test -n "$cs_cc_vers_patch" || cs_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -Wno-unused-command-line-argument"
  cflags_default_dbg="-g -O0"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_omp="-fiopenmp"
  cflags_default_omp_ad="-fopenmp-targets=spir64"

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
  cflags_default="-funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_omp="-fopenmp=libomp"

fi

# Otherwise, are we using nvc/pgcc ?
#---------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC --version 2>&1 | grep 'NVIDIA' > /dev/null
  if test "$?" = "0" ; then
    $CC -V 2>&1 | grep 'Compilers and Tools' > /dev/null
    if test "$?" = "0" ; then

      echo "compiler '$CC' is NVIDIA compiler"

      # Version strings for logging purposes and known compiler flag
      $CC -V > $outfile 2>&1
      cs_ac_cc_version=`grep pgcc $outfile | head -1`
      if test "$cs_ac_cc_version" = "" ; then
        cs_ac_cc_version=`grep nvc $outfile | head -1`
      fi
      cs_cc_compiler_known=yes

      # Default compiler flags
      cflags_default=""
      cflags_default_dbg="-g -Mbounds"
      cflags_default_opt="-O2"
      cflags_default_hot="-fast"
      cflags_default_omp="-mp=gpu"

    fi
  fi

fi

# Otherwise, are we using xlc ?
#------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -qversion 2>&1 | grep 'XL C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is IBM XL C/C++ compiler"

    # Version strings for logging purposes and known compiler flag
    $CC -qversion > $outfile 2>&1
    cs_ac_cc_version=`grep 'XL C' $outfile`
    cs_cc_compiler_known=yes
    cs_linker_set=yes

    # Default compiler flags
    cflags_default=""
    cflags_default_opt="-O3"
    cflags_default_hot="-O3"
    cflags_default_dbg="-g -qfullpath"
    cflags_default_omp="-qsmp=omp -qthreaded"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O3"
    ldflags_default_dbg="-g -qfullpath"
    ldflags_rpath="-R"

  fi
fi

# Otherwise, are we using the Cray compiler ?
#--------------------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'Cray C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Cray C compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cc_version=`$CC -V 2>&1 | grep "Cray C" | head -1`
    cs_cc_compiler_known=yes
    cs_linker_set=yes

    # Default compiler flags
    cflags_default="-h std=c11"
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

# Otherwise, are we using the Fujitsu compiler ?
#-----------------------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'Fujitsu C/C++' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Fujitsu C compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cc_version=`$CC -V 2>&1 | grep "Fujitsu C/C++" | head -1`
    cs_cc_compiler_known=yes
    cs_linker_set=yes

    # Default compiler flags
    cflags_default="-x c11 -fPIC"
    cflags_default_opt="-O2"
    cflags_default_hot="-O3"
    cflags_default_dbg="-g"
    cflags_default_omp="-Kopenmp"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O2"
    ldflags_default_dbg="-g"

  fi
fi

# Otherwise, are we using the Arm compiler ?
#-------------------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  $CC -v 2>&1 | grep 'Arm C/C++/Fortran' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Arm C compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cc_version=`$CC -v 2>&1 | grep "Arm C/C++/Fortran" | head -1`
    cs_cc_compiler_known=yes
    cs_linker_set=yes

    # Default compiler flags
    cflags_default="-std=c11 -fPIC"
    cflags_default_opt="-O2"
    cflags_default_hot="-O3"
    cflags_default_dbg="-g"
    cflags_default_omp="-fopenmp"

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

  # Intel and LLVM compilers may pass as GXX but
  # may be recognized by version string

  cs_ac_cxx_version=`$CXX $user_CXXFLAGS --version 2>&1 | head -1`

  if test -n "`echo $cs_ac_cxx_version | grep ICC`" ; then
    cs_gxx=icpc
  elif test -n "`echo $cs_ac_cxx_version | grep -e DPC++ -e oneAPI`" ; then
    cs_gxx=oneAPI
  elif test -n "`echo $cs_ac_cxx_version | grep clang`" ; then
    cs_gxx=clang
  elif test -n "`echo $cs_ac_cxx_version | grep Cray | grep -v GCC`" ; then
    cs_gxx=cray
  elif test -n "`echo $cs_ac_cxx_version | grep FCC`" ; then
    cs_gxx=fujitsu
  elif test -n "`echo $cs_ac_cxx_version | grep Arm`" ; then
    cs_gxx=arm
  elif test -n "`$CXX $user_CFLAGS --version 2>&1 | grep NVIDIA`" ; then
    cs_gxx=no
  else
    cs_gxx=g++
  fi

fi

if test "x$cs_gxx" = "xg++"; then

  # Version strings for logging purposes and known compiler flag
  $CXX -v > $outfile 2>&1
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
    g++-4.[012345678]*)
      cxxflags_default="$cxxflags_default -std=c++1y"
      ;;
    *)
      cxxflags_default="$cxxflags_default -std=c++14"
      ;;
  esac

  # Warning flags as are available by compiler version

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

  # No OpenMP in very old gcc versions

  case "$cs_cxx_vendor-$cs_cxx_version" in
    g++-2.*|g++-3*|g++-4.[012]*)
      cxxflags_default_omp=""
      ;;
  esac

# Otherwise, are we using ICPC Classic ?
#---------------------------------------

elif test "x$cs_gxx" = "xicpc"; then

  cs_cxx_version=`echo $cs_ac_cxx_version | grep ICC |sed 's/[a-zA-Z()]//g'`
  echo "compiler '$CXX' is Intel ICPC Classic"

  # Version strings for logging purposes and known compiler flag
  $CXX $user_CXXFLAGS -V conftest.c > $outfile 2>&1
  cs_cxx_compiler_known=yes

  cs_cxx_vers_major=`echo $cs_ac_cxx_version | cut -f 3 -d" " | cut -f1 -d.`
  cs_cxx_vers_minor=`echo $cs_ac_cxx_version | cut -f 3 -d" " | cut -f2 -d.`
  cs_cxx_vers_patch=`echo $cs_ac_cxx_version | cut -f 3 -d" " | cut -f3 -d.`
  test -n "$cs_cxx_vers_major" || cs_cxx_vers_major=0
  test -n "$cs_cxx_vers_minor" || cs_cxx_vers_minor=0
  test -n "$cs_cxx_vers_patch" || cs_cxx_vers_patch=0

  # Default compiler flags
  cxxflags_default="-std=c++14 -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
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

# Otherwise, are we using ICC NextGen ?  This is a deprecated beta version.
#--------------------------------------

elif test "x$cs_gxx" = "xoneAPI"; then

  cs_cxx_version=`echo $cs_ac_cxx_version | grep -e DPC++ -e oneAPI |sed 's/[a-zA-Z()+/]//g'`
  echo "compiler '$CXX' is oneAPI DPC++/C++ Compiler"

  # Version strings for logging purposes and known compiler flag
  $CXX $user_CXXFLAGS -V conftest.c > $outfile 2>&1
  cs_cxx_compiler_known=yes

  cs_cxx_vers_major=`echo $cs_ac_cxx_version | cut -f 5 -d" " | cut -f1 -d.`
  cs_cxx_vers_minor=`echo $cs_ac_cxx_version | cut -f 5 -d" " | cut -f2 -d.`
  cs_cxx_vers_patch=`echo $cs_ac_cxx_version | cut -f 5 -d" " | cut -f3 -d.`
  test -n "$cs_cxx_vers_major" || cs_cxx_vers_major=0
  test -n "$cs_cxx_vers_minor" || cs_cxx_vers_minor=0
  test -n "$cs_cxx_vers_patch" || cs_cxx_vers_patch=0

  # Default compiler flags
  cxxflags_default="-Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -Wno-unused-command-line-argument"
  cxxflags_default_dbg="-g -O0"
  cxxflags_default_opt="-O2"
  cxxflags_default_hot="-O3"
  cxxflags_default_omp="-fiopenmp"
  cxxflags_default_omp_ad="-fopenmp-targets=spir64"
  cxxflags_default_std="-funsigned-char"

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
  cxxflags_default="-Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cxxflags_default_dbg="-g -O0"
  cxxflags_default_opt="-O2"
  cxxflags_default_hot="-O3"
  cxxflags_default_omp="-fopenmp=libomp"
  cxxflags_default_std="-funsigned-char"

# Otherwise, are we using pgc++/nvc++ ?
#--------------------------------------

else

  $CXX -V 2>&1 | grep 'NVIDIA' > /dev/null
  if test "$?" = "0" ; then
    $CXX -V 2>&1 | grep 'Compilers and Tools' > /dev/null
    if test "$?" = "0" ; then

      echo "compiler '$CXX' is NVIDIA compiler"

      # Version strings for logging purposes and known compiler flag
      $CXX -V conftest.c > $outfile 2>&1
      cs_ac_cxx_version=`grep pgc++ $outfile | head -1`
      if test "$cs_ac_cxx_version" = "" ; then
        cs_ac_cxx_version=`grep nvc++ $outfile | head -1`
      fi
      cs_cxx_compiler_known=yes

      # Default compiler flags
      cxxflags_default=""
      cxxflags_default_dbg="-g -Mbounds"
      cxxflags_default_opt="-O2"
      cxxflags_default_hot="-fast"
      cxxflags_default_omp="-mp=gpu"
      cxxflags_default_std=""

    fi
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
    cxxflags_default="-h std=c++14"
    cxxflags_default_opt="-O2"
    cxxflags_default_hot="-O3"
    cxxflags_default_dbg="-g"
    cfxxlags_default_omp="-h omp"              # default: use "-h noomp" to disable
    cxxflags_default_std=""

  fi

fi

# Otherwise, are we using the Fujitsu compiler ?
#------------------------------------------

if test "x$cs_cxx_compiler_known" != "xyes" ; then

  $CXX -V 2>&1 | grep 'Fujitsu C/C++' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CXX' is Fujitsu C++"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cxx_version=`$CXX -V 2>&1 | grep "Fujitsu C/C++" | head -1`
    cs_cxx_compiler_known=yes

    # Default compiler flags
    cxxflags_default="-x c++14 -fPIC"
    cxxflags_default_opt="-O2"
    cxxflags_default_hot="-O3"  # Bug observed when -O3 is used
    cxxflags_default_dbg="-g"
    cfxxlags_default_omp="-Kopenmp"
    cxxflags_default_std=""
  fi
fi

# Otherwise, are we using the Arm compiler ?
#------------------------------------------

if test "x$cs_cxx_compiler_known" != "xyes" ; then

  $CXX -v 2>&1 | grep 'Arm C/C++/Fortran' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CXX' is Arm C++"

    # Version strings for logging purposes and known compiler flag
    cs_ac_cxx_version=`$CXX -v 2>&1 | grep "Arm C/C++/Fortran" | head -1`
    cs_cxx_compiler_known=yes

    # Default compiler flags
    cxxflags_default="-std=c++14 -fPIC"
    cxxflags_default_opt="-O2"
    cxxflags_default_hot="-O3"
    cxxflags_default_dbg="-g"
    cfxxlags_default_omp="-fopenmp"
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

cs_ac_fc_version=`$FC $user_FCFLAGS --version 2>&1 | head -1`

# Are we using gfortran ?
#------------------------

echo $cs_ac_fc_version | grep 'GNU Fortran' > /dev/null

if test "$?" = "0" ; then

  cs_fc_version=`echo $FC --version  |sed 's/[a-zA-Z()]//g'`
  cs_fc_version="`$FC -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"

  cs_fc_compiler_known=yes
  cs_gfortran=gfortran

  # Version strings for logging purposes and known compiler flag
  $FC -v > $outfile 2>&1
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
  fcflags_default="-x f95-cpp-input -Wall -pedantic-errors -std=f2008"
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

  # Are we using IFORT Classic ?
  #-----------------------------

  echo $cs_ac_fc_version | grep 'ifort' > /dev/null
  if test "$?" = "0" ; then

    cs_fc_version=`echo $cs_ac_fc_version | sed 's/[a-zA-Z()]//g'`
    cs_fc_vendor='Intel IFORT Classic'

    echo "compiler '$FC' is Intel Fortran Classic"

    # Version strings for logging purposes and known compiler flag
    $FC $user_FCLAGS -V > $outfile 2>&1
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

  # Are we using Intel LLVM FORTRAN Compiler ?
  #-------------------------------------------

  echo $cs_ac_fc_version | grep 'ifx' > /dev/null
  if test "$?" = "0" ; then

    cs_fc_version=`echo $cs_ac_fc_version | sed 's/[a-zA-Z()+/]//g'`

    echo "compiler '$FC' is Intel oneAPI Fortran"

    # Version strings for logging purposes and known compiler flag
    $FC $user_FCLAGS -V > $outfile 2>&1
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
    fcflags_default_omp="-fiopenmp"

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  cs_ac_fc_version=""

  # Are we using pgfortran/nvfortran ?
  #-----------------------------------

  $FC -V 2>&1 | grep 'NVIDIA' > /dev/null
  if test "$?" = "0" ; then
    $FC -V 2>&1 | grep 'Compilers and Tools' > /dev/null
    if test "$?" = "0" ; then

      echo "compiler '$FC' is NVIDIA compiler"

      # Version strings for logging purposes and known compiler flag
      $FC -V > $outfile 2>&1
      cs_ac_fc_version=`grep pgf $outfile | head -1`
      if test "$cs_ac_fc_version" = "" ; then
        cs_ac_fc_version=`grep nvfortran $outfile | head -1`
      fi
      cs_fc_compiler_known=yes

      # Default compiler flags
      fcflags_default="-Mpreprocess -noswitcherror"
      fcflags_default_dbg="-g -Mbounds"
      fcflags_default_opt="-O2"
      fcflags_default_hot="-fast"
      fcflags_default_omp="-mp=gpu"

    fi
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

  # Are we using the Fujitsu compiler ?
  #-------------------------------

  $FC -V 2>&1 | grep 'Fujitsu Fortran' > /dev/null

  if test "$?" = "0" ; then

    echo "compiler '$FC' is Fujitsu Fortran compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_fc_version=`$FC -V 2>&1 | grep "Fujitsu Fortran" | head -1`
    cs_fc_compiler_known=yes

    fcflags_default="-cpp -fPIC"
    fcflags_default_dbg="-g"
    fcflags_default_opt="-O2"
    fcflags_default_omp="-Kopenmp"

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  # Are we using the Arm compiler ?
  #--------------------------------

  $FC -v 2>&1 | grep 'Arm C/C++/Fortran' > /dev/null

  if test "$?" = "0" ; then

    echo "compiler '$FC' is Arm Fortran compiler"

    # Version strings for logging purposes and known compiler flag
    cs_ac_fc_version=`$FC -v 2>&1 | grep "Arm C/C++/Fortran" | head -1`
    cs_fc_compiler_known=yes

    fcflags_default="-cpp -fPIC"
    fcflags_default_dbg="-g"
    fcflags_default_opt="-O2"
    fcflags_default_omp="-fopenmp"

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


###################
#                 #
#  CUDA compiler  #
#                 #
###################

if test "x$NVCC" != "x" ; then

  # Practical version info
  cs_nvcc_version="`$NVCC --version 2>&1 |grep 'release' |\
                    sed 's/.*release \([-a-z0-9\.]*\).*/\1/'`"

  echo "compiler '${NVCC}' version ${cs_nvcc_version}"

  cs_ac_nvcc_version="`$NVCC --version 2>&1 |grep 'release' | head -1`"

  # Default compiler flags
  nvccflags_default=""
  nvccflags_default_dbg="-g -G"
  nvccflags_default_opt="-O2"
  nvccflags_default_prf="-O2 -g -lineinfo"

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
          ldflags_rpath_add="${ldflags_rpath}${libgfortran_dir}"
        fi
        libgfortran_path=`$FC --print-file-name=libgfortran.so`
        libgfortran_dir=`dirname $libgfortran_path`
        unset libgfortran_path
        if test "`echo $libgfortran_dir | cut -c1-4`" != "/usr" ; then
          ldflags_rpath_add="${ldflags_rpath}${libgfortran_dir}"
        fi
        unset libgfortran_dir
      fi
      ;;

    darwin*)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_rpath="-Wl,-rpath -Wl,"
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

