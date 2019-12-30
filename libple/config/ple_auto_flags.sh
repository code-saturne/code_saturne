# Shell script

# Copyright (C) 2005-2017 EDF

# This file is part of the PLE software package.  For license
# information, see the COPYING file in the top level directory of the
# PLE source distribution.

# This file should be sourced by configure, and sets the following
# environment variables corresponding to the recommended settings for a
# given OS/CPU/compiler combination:
#
# cppflags_default       # Base CPPFLAGS                     (default: "")
#
# cflags_default         # Base CFLAGS                       (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging    (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling    (default: "-g")
# cflags_default_omp     # Added to $CFLAGS for OpenMP       (default: "")
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

cflags_default_prf="-g"
ldflags_default_prf="-g"

ple_ac_cc_version=unknown
ple_cc_compiler_known=no

# Are we using gcc ?
#-------------------

ple_gcc=no

if test "x$GCC" = "xyes"; then

  # Intel and Pathscale compilers may pass as GCC but
  # may be recognized by version string

  if test -n "`$CC --version | grep icc`" ; then
    ple_gcc=icc
  elif test -n "`$CC --version | grep clang`" ; then
    ple_gcc=clang
  elif test -n "`$CC --version 2>&1 | grep PathScale`" ; then
    ple_gcc=pathcc
  else
    ple_gcc=gcc
  fi

fi

if test "x$ple_gcc" = "xgcc"; then

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  ple_ac_cc_version=`$CC --version 2>&1 | head -1`
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
  cflags_default="-std=c99 -funsigned-char -pedantic -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused"
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
    gcc-[45]*)
      ;;
    *)
      cflags_default="$cflags_default -Wmisleading-indentation -Wduplicated-cond"
      ;;
  esac

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-4.[012345678]*)
      ;;
    *)
      cflags_default="$cflags_default -fdiagnostics-color=auto"
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

  case "$host_os" in
    mingw64)
      cflags_default="`echo $cflags_default | sed -e 's/-std=c99/-std=gnu99/g'`"
      ;;
  esac

# Otherwise, are we using icc ?
#------------------------------

elif test "x$ple_gcc" = "xicc"; then

  ple_cc_version=`echo $CC --version | grep icc |sed 's/[a-zA-Z()]//g'`

  echo "compiler '$CC' is Intel ICC"

  # Version strings for logging purposes and known compiler flag
  $CC -V conftest.c > $outfile 2>&1
  ple_ac_cc_version=`$CC --version 2>&1 | head -1`
  ple_cc_compiler_known=yes

  # Some version numbers
  ple_cc_vers_major=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f1 -d.`
  ple_cc_vers_minor=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f2 -d.`
  ple_cc_vers_patch=`echo $ple_ac_cc_version | cut -f 3 -d" " | cut -f3 -d.`
  test -n "$ple_cc_vers_major" || ple_cc_vers_major=0
  test -n "$ple_cc_vers_minor" || ple_cc_vers_minor=0
  test -n "$ple_cc_vers_patch" || ple_cc_vers_patch=0

  # Default compiler flags
  cflags_default="-strict-ansi -std=c99 -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_omp="-qopenmp"
  case "$ple_cc_vers_major" in
    1[0123456])
      cflags_default_omp="-openmp"
      ;;
  esac

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
  cflags_default="-std=c99 -funsigned-char -Wall -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0"
  cflags_default_opt="-O2"
  cflags_default_omp="-openmp"

# Otherwise, are we using pathcc ?
#---------------------------------

elif test "x$ple_gcc" = "xpathcc"; then

  $CC --version 2>&1 | grep 'PathScale' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is PathScale C compiler"

    # Version strings for logging purposes and known compiler flag
    $CC --version > $outfile 2>&1
    ple_ac_cc_version=`grep -i Compiler $outfile`
    ple_cc_compiler_known=yes

    # Default compiler flags
    cflags_default="-c99 -noswitcherror"
    cflags_default="-std=c99 -funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wunused-value"
    cflags_default_dbg="-g"
    cflags_default_opt="-O2"
    cflags_default_omp="-openmp"

  fi

fi

# Otherwise, are we using pgcc ?
#-------------------------------

if test "x$ple_cc_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Portland Group pgcc"

    # Version strings for logging purposes and known compiler flag
    $CC -V > $outfile 2>&1
    ple_ac_cc_version=`grep -i pgcc $outfile`
    ple_cc_compiler_known=yes

    # Default compiler flags
    cflags_default="-c99 -noswitcherror"
    cflags_default_dbg="-g -Mbounds"
    cflags_default_opt="-O2"
    cflags_default_omp="-mp"

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
    cflags_default="-qlanglvl=stdc99 -q64"
    cflags_default_opt="-O3"
    cflags_default_dbg="-g -qfullpath"
    cflags_default_omp="-qsmp=omp -qthreaded"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O3"
    ldflags_default_dbg="-g -qfullpath"

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null
    if test "$?" = "0" ; then
      # Default compiler flags (we assume that MPI wrappers are used)
      ple_ibm_bg_type=`grep 'Blue Gene' $outfile | sed -e 's/.*Blue Gene\/\([A-Z]\).*/\1/'`
      if test "x$ple_ibm_bg_type" = "xP" ; then
        cppflags_default="-I/bgsys/drivers/ppcfloor/arch/include"
        cflags_default="-qlanglvl=extc99"
        cflags_default_opt="-O3"
        cflags_default_dbg="-g"
      elif test -f /bgsys/drivers/ppcfloor/cnk/bin/bgq_kernel.elf ; then
        ple_ibm_bg_type="Q"
        cppflags_default=""
        cflags_default=""                    # "-qlanglvl=extc99" by default
        cflags_default_opt="-g -O3"
        cflags_default_dbg="-g"
      fi
    fi

  fi
fi

# Otherwise, are we using the Cray compiler ?
#------------------------------------------

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

