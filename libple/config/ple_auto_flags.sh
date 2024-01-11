# Shell script

# Copyright (C) 2005-2023 EDF

# This file is part of the PLE software package.  For license
# information, see the COPYING file in the top level directory of the
# PLE source distribution.

# This file should be sourced by configure, and sets the following
# environment variables corresponding to the recommended settings for a
# given OS/CPU/compiler combination:
#
# cppflags_default       # Base CPPFLAGS                      (default: "")

# cflags_default         # Base CFLAGS                        (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging     (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization  (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling     (default: "-g")
# cflags_default_omp     # Added to $CFLAGS for OpenMP        (default: "")
# cflags_default_shared  # Added to $CFLAGS for shared libs   (default: "-fPIC -DPIC")
#
# ldflags_default        # Base LDFLAGS                       (default: "")
# ldflags_default_dbg    # Added to $LDFLAGS for debugging    (default: "-g")
# ldflags_default_opt    # Added to $LDFLAGS for optimization (default: "-O")
# ldflags_default_prf    # Added to $LDFLAGS for profiling    (default: "-g")
# ldflags_rpath          # Added to $LDFLAGS for shared libs  (default: "")

# Two other environment variable strings are defined, containing possibly
# more detailed compiler information:
#
# ple_ac_cc_version      # Compiler version string, 1 line max.
# ple_ac_cc_version_full # Compiler version string, 10 lines max.

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

outfile=ple_ac_cc_env-tmp

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

# Default pre-processor flags (not too dependent on compiler)
#----------------------------

case "$host_os" in
  *)
    cppflags_default=""
    ;;
esac

# Default compiler flags
#-----------------------

cflags_default_prf="-g"
ldflags_default_prf="-g"

# Options to generate position-independent code for shared libraries.
# can be modified later if necessary, but usually common to most compilers.

case "$host_os" in
  darwin*)
    cflags_default_shared="-fPIC -DPIC"
    ldflags_default_shared="-dynamiclib -undefined dynamic_lookup"
    ldflags_default_soname="-install_name @rpath/"
    ;;
  *)
    cflags_default_shared="-fPIC -DPIC"
    ldflags_default_shared="-shared"
    ldflags_default_soname="-Wl,-soname -Wl,"
    ;;
esac

# Compiler info (will be determined in following tests)

ple_ac_cc_version=unknown
ple_cc_compiler_known=no

# Are we using gcc ?
#-------------------

ple_gcc=no

if test "x$GCC" = "xyes"; then

  # Intel and LLVM compilers may pass as GCC but
  # may be recognized by version string

  ple_ac_cc_version=`$CC $user_CFLAGS --version 2>&1 | head -1`

  if test -n "`echo $ple_ac_cc_version | grep ICC`" ; then
    ple_gcc=icc
  elif test -n "`echo $ple_ac_cc_version | grep ICX`" ; then
    ple_gcc=icx
  elif test -n "`echo $ple_ac_cc_version | grep -e DPC++ -e oneAPI`" ; then
    ple_gcc=oneAPI
  elif test -n "`echo $ple_ac_cc_version | grep clang`" ; then
    ple_gcc=clang
  elif test -n "`echo $ple_ac_cc_version | grep Cray | grep -v GCC`" ; then
    ple_gcc=cray
  elif test -n "`echo $ple_ac_cc_version | grep FCC`" ; then
    ple_gcc=fujitsu
  elif test -n "`echo $ple_ac_cc_version | grep Arm`" ; then
    ple_gcc=arm
  elif test -n "`$CC $user_CFLAGS --version 2>&1 | grep NVIDIA`" ; then
    ple_gcc=no
  else
    ple_gcc=gcc
  fi

fi

if test "x$ple_gcc" = "xgcc"; then

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  ple_cc_compiler_known=yes

  # Practical version info for option setting
  ple_cc_version="`$CC -v 2>&1 |grep 'gcc version' |\
                  sed 's/.*gcc version \([-a-z0-9\.]*\).*/\1/'`"
  ple_cc_vendor=`echo $ple_cc_version |sed 's/\([a-z]*\).*/\1/'`
  ple_cc_version=`echo $ple_cc_version |sed 's/[-a-z]//g'`

  if test "x" = "x$ple_cc_vendor" -a "x" != "x$ple_cc_version"; then
    ple_cc_vendor=gcc
  fi
  if test "-" != "$ple_cc_vendor-$ple_cc_version"; then
    echo "compiler '$CC' is GNU $ple_cc_vendor-$ple_cc_version"
  fi

  # Some version numbers
  ple_cc_vers_major=`echo $ple_cc_version | cut -f1 -d.`
  ple_cc_vers_minor=`echo $ple_cc_version | cut -f2 -d.`
  ple_cc_vers_patch=`echo $ple_cc_version | cut -f3 -d.`
  test -n "$ple_cc_vers_major" || ple_cc_vers_major=0
  test -n "$ple_cc_vers_minor" || ple_cc_vers_minor=0
  test -n "$ple_cc_vers_patch" || ple_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wfloat-equal -Werror=implicit-function-declaration"
  cflags_default_dbg="-g"
  cflags_default_opt="-O2"
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

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-[4]*)
      cflags_default="$cflags_default -std=c11"
      ;;
    gcc-[5]*)
      ;;
    *)
      cflags_default="$cflags_default -Wmisleading-indentation -Wduplicated-cond"
      ;;
  esac

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-4.[012345678]*)
      ;;
    *)
      cflags_default="$cflags_default -fdiagnostics-color=auto -Werror=format-security"
      ;;
  esac

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-4.[012]*)
      cflags_default_omp=""
      ;;
  esac

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-4.[01234]*)
      ;;
    *)
      cflags_default_opt="$cflags_default_opt -fexcess-precision=fast"
      ;;
  esac

# Otherwise, are we using ICC Classic ?
#--------------------------------------

elif test "x$ple_gcc" = "xicc" ; then

  ple_cc_version=`echo $ple_ac_cc_version | grep ICC |sed 's/[a-zA-Z()]//g'`
  echo "compiler '$CC' is Intel ICC Classic"

  # Version strings for logging purposes and known compiler flag
  $CC $user_CFLAGS -V conftest.c > $outfile 2>&1
  ple_cc_compiler_known=yes

  # Some version numbers
  ple_cc_vers_major=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f1 -d.`
  ple_cc_vers_minor=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f2 -d.`
  ple_cc_vers_patch=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f3 -d.`
  test -n "$ple_cc_vers_major" || ple_cc_vers_major=0
  test -n "$ple_cc_vers_minor" || ple_cc_vers_minor=0
  test -n "$ple_cc_vers_patch" || ple_cc_vers_patch=0

  # Default compiler flags
  # (temporarily disable "operands evaluated in unspecified order" remark -- 981)
  cflags_default="-std=c99 -restrict -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -wd981"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_omp="-qopenmp"

  case "$ple_cc_vers_major" in
    1[0123456])
      cflags_default_omp="-openmp"
      ;;
  esac

# Otherwise, are we using ICC NextGen ?  This is a deprecated beta version.
#--------------------------------------

elif test "x$ple_gcc" = "xicx" ; then

  ple_cc_version=`echo $ple_ac_cc_version | grep ICX |sed 's/[a-zA-Z()]//g'`
  echo "compiler '$CC' is Intel ICC NextGen"

  # Version strings for logging purposes and known compiler flag
  $CC $user_CFLAGS -V conftest.c > $outfile 2>&1
  ple_cc_compiler_known=yes

  # Some version numbers
  ple_cc_vers_major=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f1 -d.`
  ple_cc_vers_minor=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f2 -d.`
  ple_cc_vers_patch=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f3 -d.`
  test -n "$ple_cc_vers_major" || ple_cc_vers_major=0
  test -n "$ple_cc_vers_minor" || ple_cc_vers_minor=0
  test -n "$ple_cc_vers_patch" || ple_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_omp="-qopenmp"

# Otherwise, are we using Intel LLVM DPC++/C++ Compiler (OneAPI) ?
#-----------------------------------------------------------------

elif test "x$ple_gcc" = "xoneAPI" ; then

  ple_cc_version=`echo $ple_ac_cc_version | grep -e DPC++ -e oneAPI |sed 's/[a-zA-Z()+/]//g'`
  echo "compiler '$CC' is oneAPI DPC++/C++ Compiler"

  # Version strings for logging purposes and known compiler flag
  $CC $user_CFLAGS -V conftest.c > $outfile 2>&1
  ple_cc_compiler_known=yes

  # Some version numbers
  ple_cc_vers_major=`echo $ple_ac_cc_version | cut -f 5 -d" " | cut -f1 -d.`
  ple_cc_vers_minor=`echo $ple_ac_cc_version | cut -f 5 -d" " | cut -f2 -d.`
  ple_cc_vers_patch=`echo $ple_ac_cc_version | cut -f 5 -d" " | cut -f3 -d.`
  test -n "$ple_cc_vers_major" || ple_cc_vers_major=0
  test -n "$ple_cc_vers_minor" || ple_cc_vers_minor=0
  test -n "$ple_cc_vers_patch" || ple_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -Wno-unused-command-line-argument"
  cflags_default_dbg="-g -O0"
  cflags_default_opt="-O2"
  cflags_default_omp="-fiopenmp"

# Otherwise, are we using clang ?
#--------------------------------

elif test "x$ple_gcc" = "xclang"; then

  ple_cc_version=`echo $CC --version | grep clang | cut -f 3 -d ' '`

  echo "compiler '$CC' is clang"

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  ple_ac_cc_version=`$CC --version 2>&1 | head -1`
  ple_cc_compiler_known=yes

  # Default compiler flags
  # (temporarily disable "operands evaluated in unspecified order" remark -- 981)
  cflags_default="-std=c99 -funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0"
  cflags_default_opt="-O2"
  cflags_default_omp="-fopenmp=libomp"

fi

# Otherwise, are we using nvc/pgcc ?
#---------------------------------

if test "x$ple_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'NVIDIA' > /dev/null
  if test "$?" = "0" ; then
    $CC -V 2>&1 | grep 'Compilers and Tools' > /dev/null
    if test "$?" = "0" ; then

      echo "compiler '$CC' is NVIDIA compiler"

      # Version strings for logging purposes and known compiler flag
      $CC -V > $outfile 2>&1
      ple_ac_cc_version=`grep pgcc $outfile | head -1`
      if test "$ple_ac_cc_version" = "" ; then
        ple_ac_cc_version=`grep nvc $outfile | head -1`
      fi
      ple_cc_compiler_known=yes

      # Default compiler flags
      cflags_default=""
      cflags_default_dbg="-g -Mbounds"
      cflags_default_opt="-O2"
      cflags_default_omp="-mp"

    fi
  fi

fi

# Otherwise, are we using xlc ?
#------------------------------

if test "x$ple_cc_compiler_known" != "xyes" ; then

  $CC -qversion 2>&1 | grep 'XL C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is IBM XL C compiler"

    # Version strings for logging purposes and known compiler flag
    $CC -qversion > $outfile 2>&1
    ple_ac_cc_version=`grep 'XL C' $outfile`
    ple_cc_compiler_known=yes
    ple_linker_set=yes

    # Default compiler flags
    cflags_default=""
    cflags_default_opt="-O3"
    cflags_default_dbg="-g -qfullpath"
    cflags_default_omp="-qsmp=omp -qthreaded"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O3"
    ldflags_default_dbg="-g -qfullpath"

  fi
fi

# Otherwise, are we using the Cray compiler ?
#--------------------------------------------

if test "x$ple_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'Cray C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Cray C compiler"

    # Version strings for logging purposes and known compiler flag
    ple_ac_cc_version=`$CC -V 2>&1 | grep "Cray C" | head -1`
    ple_cc_compiler_known=yes
    ple_linker_set=yes

    # Default compiler flags
    cflags_default=""                        # "-h c99" by default
    cflags_default_opt="-O2"
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

if test "x$ple_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'Fujitsu C/C++' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Fujitsu C compiler"

    # Version strings for logging purposes and known compiler flag
    ple_ac_cc_version=`$CC -V 2>&1 | grep "Fujitsu C" | head -1`
    ple_cc_compiler_known=yes
    ple_linker_set=yes

    # Default compiler flags
    cflags_default="-x c11 -fPIC"
    cflags_default_opt="-O2"
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

if test "x$ple_cc_compiler_known" != "xyes" ; then

  $CC -v 2>&1 | grep 'Arm C/C++/Fortran' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Arm C compiler"

    # Version strings for logging purposes and known compiler flag
    ple_ac_cc_version=`$CC -v 2>&1 | grep "Arm C/C++/Fortran" | head -1`
    ple_cc_compiler_known=yes
    ple_linker_set=yes

    # Default compiler flags
    cflags_default="-std=c11 -fPIC"
    cflags_default_opt="-O2"
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

if test "x$ple_cc_compiler_known" != "xyes" ; then

  case "$host_os" in

    SUPER-UX* | superux*)

     # Native NEC SX vectorizing C compiler (sxmpicc)
     #-------------------------------------

     $CC -V conftest.c 2>&1 | grep 'NEC' | grep 'SX' > /dev/null
     if test "$?" = "0" ; then

       echo "compiler '$CC' is NEC SX compiler"

       # Version strings for logging purposes and known compiler flag
       $CC -V conftest.c > $outfile 2>&1
       ple_ac_cc_version=`grep ccom $outfile`
       ple_cc_compiler_known=yes
       ple_linker_set=yes

       # Default compiler flags
       cflags_default="-Kc99 -pvctl,loopcnt=2147483647"
       cflags_default_opt=""
       cflags_default_dbg=""
       cflags_default_omp=""

       # Default linker flags
       ldflags_default=""
       ldflags_default_opt="-O"
       ldflags_default_dbg="-g"

     fi
     ;;

    *)

      # Unknown
      #--------

      cflags_default=""
      cflags_default_opt="-O"
      cflags_default_dbg="-g"
      cflags_default_omp=""
      ;;

  esac

fi


if test -f $outfile ; then
  ple_ac_cc_version_full=`sed -e '11,$d' $outfile`
  rm -f $outfile
fi

# Default linker flags
#---------------------

if test "x$ple_linker_set" != "xyes" ; then

  case "$host_os" in

    linux*)
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

