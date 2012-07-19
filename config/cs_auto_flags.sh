# Shell script

# This file is part of the Code_Saturne Kernel, element of the
# Code_Saturne CFD tool.
#
# Copyright (C) 2009-2010 EDF S.A., France
#
# The Code_Saturne Kernel is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#
# The Code_Saturne Kernel is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public Licence
# along with the Code_Saturne Preprocessor; if not, write to the
# Free Software Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# This file should be sourced by configure, and sets the following
# environment variables corresponding to the recommended settings for a
# given OS/CPU/compiler combination:
#
# cppflags_default       # Base CPPFLAGS                     (default: "")

# cflags_default         # Base CFLAGS                       (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging    (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization (default: "-O")
# cflags_default_hot     # Optimization for specific files   (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling    (default: "")
# cflags_default_ext     # Added to $CFLAGS for extended     (default: "")
#                        # precision if available

# fcflags_default        # Base FCFLAGS                       (default: "")
# fcflags_default_dbg    # Added to $FCFLAGS for debugging    (default: "-g")
# fcflags_default_opt    # Added to $FCFLAGS for optimization (default: "-O")
# fcflags_default_hot    # Optimization for specific files    (default: "-O")
# fcflags_default_prf    # Added to $FCFLAGS for profiling    (default: "")
#
# ldflags_default        # Base LDFLAGS                       (default: "")
# ldflags_default_dbg    # Added to $LDFLAGS for debugging    (default: "-g")
# ldflags_default_opt    # Added to $LDFLAGS for optimization (default: "-O")
# ldflags_default_prf    # Added to $LDFLAGS for profiling    (default: "")
# ldflags_rpath          # Added to $LDFLAGS for shared libs  (default: "")

# libs_default           # Base LIBS                          (default: "")
# libs_default_dbg       # Added to $LDFLAGS for debugging    (default: "")
# libs_default_opt       # Added to $LDFLAGS for optimization (default: "")
# libs_default_prf       # Added to $LDFLAGS for profiling    (default: "")

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

# Are we using gcc ?
#-------------------

cs_gcc=no

if test "x$GCC" = "xyes"; then

  # Intel compiler passes as GCC but may be recognized by version string
  if test -n "`$CC --version | grep icc`" ; then
    cs_gcc=icc
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
  cflags_default="-ansi -funsigned-char -pedantic -W -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs -Wunused"
  cflags_default_dbg="-g"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_prf="-pg"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in

    *i?86|*x86_64)
      cflags_default_opt="-funroll-loops -O2 -Wuninitialized"
      case "$host_cpu" in
        i686)
          case "$cs_cc_vendor-$cs_cc_version" in
            gcc-2.9[56]*|gcc-3*|gcc-4*)
              cflags_default_opt="$cflags_default_opt -march=i686"
              ;;
          esac
          ;;

      esac
      ;;

    *alphaev6|*alphaev67|*alphaev68|*alphaev7)
      cflags_default_opt="-mcpu=ev6 -O"
      ;;

  esac

  # Modify default flags depending on gcc version (as older versions
  # may not handle all flags)

  case "$cs_cc_vendor-$cs_cc_version" in

    gcc-2.9[56]*)
      cflags_default="$cflags_default -Wno-long-long"
      ;;

    gcc-3.*|gcc-4.*)
      cflags_default="`echo $cflags_default | sed -e 's/-ansi/-std=c99/g'`"
      cflags_default="$cflags_default -Wfloat-equal"
      ;;

  esac

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-3*|gcc-4.[01234]*)
      ;;
    *)
      cflags_default_opt="$cflags_default_opt -fexcess-precision=fast"
      cflags_default_hot="$cflags_default_hot -fexcess-precision=fast"
      ;;
  esac

  case "$host_os" in
    *cygwin)
    cflags_default="`echo $cflags_default | sed -e 's/c99/gnu99/g'`"
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

  # Default compiler flags
  # (temporarily disable "operands evaluated in unspecified order" remark -- 981)
  cflags_default="-strict-ansi -std=c99 -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused -wd981"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp64 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_prf="-p"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in
    *ia64)
      cflags_default_opt="-O2 -mcpu=itanium2-p9000"
      cflags_default_hot="-O3 -mcpu=itanium2-p9000"
      cflags_default_ext="-fp-model extended"
      ;;
  esac

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
    cflags_default_prf="-Mprof=func,lines"

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
    cflags_default_dbg="-g"
    cflags_default_prf="-pg"

    # Default  linker flags
    ldflags_default=""
    ldflags_default_opt="-O3"
    ldflags_default_dbg="-g"
    ldflags_default_prf="-pg"

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null
    if test "$?" = "0" ; then
      # Default compiler flags (we assume that MPI wrappers are used)
      cs_ibm_bg_type=`grep 'Blue Gene' $outfile | sed -e 's/.*Blue Gene\/\([A-Z]\).*/\1/'`
      if test "x$cs_ibm_bg_type" = "xL" ; then
        cppflags_default="-I/bgl/BlueLight/ppcfloor/bglsys/include"
        cflags_default="-qlanglvl=stdc99"
        cflags_default_opt="-O3"
        cflags_default_hot="-O3 -qhot"
        cflags_default_dbg="-g"
        ldflags_default="-Wl,-allow-multiple-definition"
      elif test "x$cs_ibm_bg_type" = "xP" ; then
        cppflags_default="-I/bgsys/drivers/ppcfloor/arch/include"
        cflags_default="-qlanglvl=extc99"
        cflags_default_opt="-O3"
        cflags_default_hot="-O3 -qhot"
        cflags_default_dbg="-g"
        ldflags_default="-Wl,-allow-multiple-definition"
      else
        cs_ibm_bg_type="Q"
        cppflags_default=""
        cflags_default=""                    # "-qlanglvl=extc99" by default
        cflags_default_opt="-O3"
        cflags_default_hot="-O3 -qhot"
        cflags_default_dbg="-g"
        ldflags_default="-Wl,-allow-multiple-definition"
      fi
    fi

  fi
fi

# Compiler still not identified
#------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  case "$host_os" in

    osf*)

      # Native Compaq Tru64 Unix C compiler
      #------------------------------------

      $CC -V 2>&1 | grep 'Compaq Tru64' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is Compaq Tru64 compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        cs_ac_cc_version=`grep 'Compaq C' $outfile`
        cs_cc_compiler_known=yes
        cs_linker_set=yes

        # Default compiler flags
        case "$host_cpu" in
          alphaev6|alphaev67|alphaev68|alphaev7)
            cflags_default="-arch host -tune host -ansi_alias -std -check_bounds -trapuv -check -msg_enable alignment -msg_enable noansi -msg_enable performance -portable -msg_enable c_to_cxx"
            cflags_default_opt="-O"
            cflags_default_hot="-O"
            cflags_default_dbg="-g"
            cflags_default_prf="-pg"
          ;;
        esac

        # Default  linker flags
        ldflags_default="-Wl,-call_shared"
        ldflags_default_opt="-O0"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-pg"
        ldflags_rpath="-Wl,-rpath -Wl,"
      fi
      ;;

    SUPER-UX* | superux*)

      # Native NEC SX vectorizing C compiler (sxmpicc)
      #-------------------------------------

      $CC -V conftest.c 2>&1 | grep 'NEC' | grep 'SX' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is NEC SX compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        cs_ac_cc_version=`grep ccom $outfile`
        cs_cc_compiler_known=yes
        cs_linker_set=yes

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

    irix5.*|irix6.*)

      # Native SGI IRIX C compiler
      #---------------------------

      $CC -version 2>&1 | grep 'MIPSpro' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is MIPSpro compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -version > $outfile 2>&1
        cs_ac_cc_version=`grep MIPSpro $outfile`
        cs_cc_compiler_known=yes

        # Default compiler flags
        cflags_default="-c99 -64"
        cflags_default_opt="-O2 -woff 1521,1552,1096"
        cflags_default_dbg="-g -woff 1429,1521,1209 -fullwarn"
        cflags_default_prf="$cflags_default_opt"

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
        cs_ac_cc_version=`grep ccom $outfile`
        cs_cc_compiler_known=yes
        cs_linker_set=yes

        # Default compiler flags
        cflags_default="-AC99 +e"
        cflags_default_opt="+O2"
        cflags_default_hot="+O3"
        cflags_default_dbg="-g"
        cflags_default_prf="-G"

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
        cs_ac_cc_version=`grep cc $outfile`
        cs_cc_compiler_known=yes

        # Default compiler flags
        cflags_default="-Xa -Xc99"
        cflags_default_opt="-xO2"
        cflags_default_hot="-xO3"
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

if test "x$cflags_default_hot" = "x" ; then
  cflags_default_hot=$cflags_default_opt
fi

if test -f $outfile ; then
  cs_ac_cc_version_full=`sed -e '11,$d' $outfile`
  rm -f $outfile
fi


######################
#                    #
#  Fortran compiler  #
#                    #
######################

cs_ac_fc_version=unknown
cs_fc_compiler_known=no

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

  echo "compiler '$FC' is gfortran"

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
  fi

  # Some version numbers
  cs_fc_vers_major=`echo $cc_version | cut -f1 -d.`
  cs_fc_vers_minor=`echo $cc_version | cut -f2 -d.`
  cs_fc_vers_patch=`echo $cc_version | cut -f3 -d.`
  test -n "$cs_fc_vers_major" || cs_fc_vers_major=0
  test -n "$cs_fc_vers_minor" || cs_fc_vers_minor=0
  test -n "$cs_fc_vers_patch" || cs_fc_vers_patch=0

  # Default compiler flags
  fcflags_default="-x f95-cpp-input -Wall -Wno-unused"
  fcflags_default_dbg="-g -fbounds-check"
  fcflags_default_opt="-O"
  fcflags_default_hot="-O2"
  fcflags_default_prf="-pg"

  # Deactivate bounds checking with older gfortran 4.1
  # to avoid issue in turrij.f90
  case "$cs_fc_vers_major-$cs_fc_vers_minor" in
    gcc-4.[01]*)
     fcflags_default_dbg="-g"
      ;;
  esac

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

    # Default compiler flags
    # (temporarily disable "unused variable" remark -- 7712)
    fcflags_default="-cpp -fpic -warn -diag-disable 7712"
    fcflags_default_dbg="-g -O0 -traceback -check all -fpe0 -ftrapuv"
    fcflags_default_opt="-O2"
    fcflags_default_hot="-O3"
    fcflags_default_prf="-p"

    # Modify default flags on certain systems

    case "$host-os-$host_cpu" in
      *ia64)
        fcflags_default_opt="-O2 -mcpu=itanium2-p9000"
        fcflags_default_hot="-O3 -mcpu=itanium2-p9000"
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
    fcflags_default_prf="-Mprof=func,lines"

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
    fcflags_default_prf="-pg"

    # Adjust options for IBM Blue Gene cross-compiler

    grep 'Blue Gene' $outfile > /dev/null

    if test "$?" = "0" ; then

      # Default compiler flags (we assume that MPI wrappers are used)
      cs_ibm_bg_type=`grep 'Blue Gene' $outfile | sed -e 's/.*Blue Gene\/\([A-Z]\).*/\1/'`
      if test "$cs_ibm_bg_type" = "L" ; then
        fcflags_default="-qextname -qsuffix=cpp=f90"
        fcflags_default_dbg="-g"
        fcflags_default_opt="-O3"
        fcflags_default_hot="-O3 -qhot"
      elif test "$cs_ibm_bg_type" = "P" ; then
        fcflags_default="-qextname -qsuffix=cpp=f90"
        fcflags_default_dbg="-g"
        fcflags_default_opt="-O3"
        fcflags_default_hot="-O3 -qhot"
      elif test "x$cs_ibm_bg_type" = "xQ" ; then
        fcflags_default="-qextname -qsuffix=cpp=f90"
        fcflags_default_dbg="-g -qcheck"
        fcflags_default_opt="-O3"
        fcflags_default_hot="-O3 -qhot"
      fi
    fi

  fi
fi

if test "x$cs_fc_compiler_known" != "xyes" ; then

  case "$host_os" in

    SUPER-UX* | superux*)

      # Native NEC SX vectorizing Fortran compiler (sxmpif90)
      #-------------------------------------------

      $FC -V conftest.f 2>&1 | grep 'NEC' | grep 'SX' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$FC' is NEC SX compiler"

        # Version strings for logging purposes and known compiler flag
        $FC -V conftest.f > $outfile 2>&1
        cs_ac_fc_version=`grep ccom $outfile`
        cs_fc_compiler_known=yes

        # Default compiler flags
        fcflags_default="-Ep -C hopt -I. -Wf,-pvctl,loopcnt=2147483647"
        fcflags_default_opt=""
        fcflags_default_hot=""
        fcflags_default_dbg=""
        fcflags_default_prf=""

      fi
      ;;

    hpux*)

      # Native HP-UX Fortran compiler
      #------------------------------

      $FC +version 2>&1 | grep 'HP' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$FC' is HP compiler"

        # Version strings for logging purposes and known compiler flag
        $FC -V conftest.f > $outfile 2>&1
        cs_ac_fc_version=`$FC +version`
        cs_fc_compiler_known=yes

        # Default compiler flags
        fcflags_default="+cpp=yes +fp_exception +FPVZOUD +U77 +DA2.0W +noppu"
        fcflags_default_opt="+O1"
        fcflags_default_hot="+O2"
        fcflags_default_dbg="-g"
        fcflags_default_prf="-G"

      fi
      ;;

    *)

      # Unknown
      #--------

      fcflags_default=""
      fcflags_default_opt="-O"
      fcflags_default_hot="-O"
      fcflags_default_dbg="-g"
      fcflags_default_prf=""
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

if test "x$cs_linker_set" != "xyes" ; then

  case "$host_os" in

    linux*)
      ldflags_default="-Wl,-export-dynamic"
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g"
      ldflags_default_prf="-pg"
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

    osf*)
      ldflags_default=""
      ldflags_default_opt="-O"
      ldflags_default_dbg="-g3"
      ldflags_default_prf="-pg"
      ;;

    irix5.*|irix6.*)
      ldflags_default="-64 -Wl,-woff,85"
      ldflags_default_opt=""
      ldflags_default_dbg="-g"
      ldflags_default_prf="-p"
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

# Clean temporary files

rm -f conftest* a.out $outfile

