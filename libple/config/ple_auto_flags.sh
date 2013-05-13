# Shell script

# Copyright (C) 2005-2013 EDF

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
# cflags_default_prf     # Added to $CFLAGS for profiling    (default: "")
#
# ldflags_default        # Base LDFLAGS                       (default: "")
# ldflags_default_dbg    # Added to $LDFLAGS for debugging    (default: "-g")
# ldflags_default_opt    # Added to $LDFLAGS for optimization (default: "-O")
# ldflags_default_prf    # Added to $LDFLAGS for profiling    (default: "")
#
# ple_disable_shared     # Disable shared librairies          (default: "")

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

ple_ac_cc_version=unknown

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

# Are we using gcc ?
#-------------------

ple_gcc=no

if test "x$GCC" = "xyes"; then

  # Intel compiler passes as GCC but may be recognized by version string
  if test -n "`$CC --version | grep icc`" ; then
    ple_gcc=icc
  elif test -n "`$CC --version 2>&1 | grep PathScale`" ; then
    ple_gcc=pathcc
  elif test -n "`$CC --version 2>&1 | grep Open64`" ; then
    ple_gcc=open64
  else
    ple_gcc=gcc
  fi

fi

if test "x$ple_gcc" = "xgcc"; then

  # Version strings for logging purposes and known compiler flag
  $CC -v > $outfile 2>&1
  ple_ac_cc_version=`$CC --version 2>&1 | head -1`
  ple_compiler_known=yes

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
  cflags_default_prf="-pg"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in

    *i?86|*x86_64)
      cflags_default_opt="-funroll-loops -O2 -Wuninitialized"
      case "$host_cpu" in
        i686)
          case "$ple_cc_vendor-$ple_cc_version" in
            gcc-3*|gcc-4*)
              cflags_default_opt="$cflags_default_opt -march=i686"
              ;;
          esac
          ;;
      esac
      ;;

  esac

  # Modify default flags depending on gcc version (as older versions
  # may not handle all flags)

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-3.*|gcc-4.*)
      cflags_default="`echo $cflags_default | sed -e 's/-ansi/-std=c99/g'`"
      cflags_default="$cflags_default -Wfloat-equal"
      ;;

  esac

  case "$ple_cc_vendor-$ple_cc_version" in
    gcc-3*|gcc-4.[01234]*)
      ;;
    *)
      cflags_default_opt="$cflags_default_opt -fexcess-precision=fast"
      ;;
  esac

  case "$host_os" in
    mingw32)
    cflags_default="$cflags_default -Wno-format -Wno-pedantic-ms-format"
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
  ple_compiler_known=yes

  # Default compiler flags
  cflags_default="-strict-ansi -std=c99 -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_prf="-p"

# Otherwise, are we using pathcc ?
#---------------------------------

elif test "x$ple_cc_compiler_known" != "xyes" ; then

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
    cflags_default_prf=""

  fi

# Otherwise, are we using opencc ?
#---------------------------------

elif test "x$ple_gcc" = "xopen64"; then

  $CC --version 2>&1 | grep 'Open64' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Open64 C compiler"

    # Version strings for logging purposes and known compiler flag
    $CC --version > $outfile 2>&1
    ple_ac_cc_version=`grep -i Compiler $outfile`
    ple_cc_compiler_known=yes

    # Default compiler flags
    cflags_default="-std=c99"
    cflags_default="-std=c99 -funsigned-char -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused -Wunused-value"
    cflags_default_dbg="-g"
    cflags_default_opt="-O2"
    cflags_default_prf=""

  fi

fi

# Otherwise, are we using pgcc ?
#-------------------------------

if test "x$ple_compiler_known" != "xyes" ; then

  $CC -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Portland Group pgcc"

    # Version strings for logging purposes and known compiler flag
    $CC -V > $outfile 2>&1
    ple_ac_cc_version=`grep -i pgcc $outfile`
    ple_compiler_known=yes

    # Default compiler flags
    cflags_default="-c99"
    cflags_default_dbg="-g -Mbounds"
    cflags_default_opt="-O2"
    cflags_default_prf="-Mprof=func,lines"

  fi

fi

# Otherwise, are we using xlc ?
#------------------------------

if test "x$ple_compiler_known" != "xyes" ; then

  $CC -qversion 2>&1 | grep 'XL C' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is IBM XL C compiler"

    # Version strings for logging purposes and known compiler flag
    $CC -qversion > $outfile 2>&1
    ple_ac_cc_version=`grep 'XL C' $outfile`
    ple_compiler_known=yes
    ple_linker_set=yes

    # Default compiler flags
    cflags_default="-qlanglvl=stdc99 -q64"
    cflags_default_opt="-O2"
    cflags_default_dbg="-g"
    cflags_default_prf="-pg"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O2"
    ldflags_default_dbg="-g"
    ldflags_default_prf="-pg"

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null
    if test "$?" = "0" ; then
      # Default compiler flags
      ple_ibm_bg_type=`grep 'Blue Gene' $outfile | sed -e 's/.*Blue Gene\/\([A-Z]\).*/\1/'`
      if test "$ple_ibm_bg_type" = "L" ; then
        cppflags_default="-I/bgl/BlueLight/ppcfloor/bglsys/include"
        cflags_default="-g -qmaxmem=-1 -qarch=440d -qtune=440"
        cflags_default_opt="-O2"
        cflags_default_dbg=""
        ldflags_default="-L/bgl/BlueLight/ppcfloor/bglsys/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts -lnss_files -lnss_dns -lresolv"
        ple_disable_shared=yes # Disable shared libraries
      elif test "$ple_ibm_bg_type" = "P" ; then
        cppflags_default="-I/bgsys/drivers/ppcfloor/comm/include"
        cflags_default="-g -qmaxmem=-1 -qarch=450d -qtune=450"
        cflags_default_opt="-O1"
        cflags_default_dbg=""
        ldflags_default=""
        ple_disable_shared=yes # Disable shared libraries
      else
        ple_ibm_bg_type="Q"
        cppflags_default=""
        cflags_default=""
        cflags_default_opt="-O3"
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
    cflags_default_prf="-h profile_generate" # resulting code must be run under CrayPat

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O2"
    ldflags_default_dbg="-g"
    ldflags_default_prf="-h profile_generate"

  fi
fi

# Compiler still not identified
#------------------------------

if test "x$ple_compiler_known" != "xyes" ; then

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
        ple_compiler_known=yes
        ple_linker_set=yes

        # Default compiler flags
        cflags_default="-Kc99 -pvctl,loopcnt=2147483647"
        cflags_default_opt=""
        cflags_default_dbg=""
        cflags_default_prf=""

        # Default linker flags
        ldflags_default=""
        ldflags_default_opt="-O"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-pg"

      fi
      ;;

    hpux*)

      # Native HP-UX C compiler
      #------------------------

      $CC -V conftest.c 2>&1 | grep 'HP' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is HP compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        ple_ac_cc_version=`grep ccom $outfile`
        ple_compiler_known=yes
        ple_linker_set=yes

        # Default compiler flags
        cflags_default="-AC99 +e"
        cflags_default_opt="+O2"
        cflags_default_dbg="-g"
        cflags_default_prf="-G" # -G for gprof, -p for prof

        # Default linker flags
        ldflags_default="+FPVZOUD +U77"
        ldflags_default_opt="+O1"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-fbexe"

        if test "$host_cpu" = "ia64" ; then
          cflags_default="$cflags_default +DD64"
          ldflags_default="$ldflags_default +DD64"
        else
          cflags_default="$cflags_default +DA2.0w"
          ldflags_default="$ldflags_default +DA2.0w"
        fi

      fi
      ;;

    solaris2.*)

      # Sun Workshop compiler
      #----------------------

      $CC -V 2>&1 | grep 'WorkShop Compilers' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is Sun WorkShop Compilers"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        ple_ac_cc_version=`grep cc $outfile`
        ple_compiler_known=yes

        # Default compiler flags
        cflags_default="-Xa -Xc99"
        cflags_default_opt="-xO2"
        cflags_default_dbg="-g"
        cflags_default_prf="-pg"

     fi
     ;;

    *)

      # Unknown
      #--------

      cflags_default=""
      cflags_default_opt="-O"
      cflags_default_dbg="-g"
      cflags_default_prf=""

      ;;

  esac

fi

# Default linker flags
#---------------------

if test "x$ple_linker_set" != "xyes" ; then

  case "$host_os" in

    linux*)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_default_prf="-pg"
      ;;

    solaris2.*)
      ldflags_default_opt=""
      ldflags_default_dbg="-g"
      ldflags_default_prf=""
      ;;

    *)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_default_prf="-pg"
      ;;

  esac

fi

# Finish

export LANG=$save_LANG

if test -f $outfile ; then 
  ple_ac_cc_version_full=`sed -e '11,$d' $outfile`
fi

# Clean temporary files

rm -f conftest* a.out $outfile

