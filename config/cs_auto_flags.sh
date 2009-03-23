# Shell script

# This file is part of the Code_Saturne Kernel, element of the
# Code_Saturne CFD tool.
#
# Copyright (C) 2009 EDF S.A., France
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
# cflags_default_omp     # Added to $CFLAGS for OpenMP       (default: "")

# fcflags_default        # Base FCFLAGS                       (default: "")
# fcflags_default_dbg    # Added to $FCFLAGS for debugging    (default: "-g")
# fcflags_default_opt    # Added to $FCFLAGS for optimization (default: "-O")
# fcflags_default_hot    # Optimization for specific files    (default: "-O")
# fcflags_default_prf    # Added to $FCFLAGS for profiling    (default: "")
# fcflags_default_omp    # Added to $FCFLAGS for OpenMP       (default: "")
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

# cs_disable_shared      # Disable shared librairies          (default: "")

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
  cs_cc_vendor=`echo $cc_version |sed 's/\([a-z]*\).*/\1/'`
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
  cflags_default_omp="-fopenmp"

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
      cflags_default="$cflags_default -Wfloat-equal -Wpadded"
      ;;

  esac

  case "$cs_cc_vendor-$cs_cc_version" in
    gcc-2.*|gcc-3*|gcc-4.[012]*)
      cflags_default_omp=""
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
  cflags_default="-strict-ansi -std=c99 -funsigned-char -Wall -Wcheck -Wshadow -Wpointer-arith -Wmissing-prototypes -Wuninitialized -Wunused"
  cflags_default_dbg="-g -O0 -traceback -w2 -Wp66 -ftrapuv"
  cflags_default_opt="-O2"
  cflags_default_hot="-O3"
  cflags_default_prf="-p"
  cflags_default_omp="-openmp"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in
    *ia64)
      cflags_default_opt="-O2 -mcpu=itanium2-p9000"
      ;;
  esac

# Otherwise, are we using pgcc ?
#-------------------------------

else

  $CC -V 2>&1 | grep 'The Portland Group' > /dev/null
  if test "$?" = "0" ; then

    echo "compiler '$CC' is Portland Group pgcc"

    # Version strings for logging purposes and known compiler flag
    $CC -V conftest.c > $outfile 2>&1
    cs_ac_cc_version=`grep -i pgcc $outfile`
    cs_cc_compiler_known=yes

    # Default compiler flags
    cflags_default="-Xa -fPIC"
    cflags_default_dbg="-g"
    cflags_default_opt="-fast -fastsse"
    cflags_default_prf="-Mprof=func,lines"
    cflags_default_omp="-mp"

  fi

fi

# Compiler still not identified
#------------------------------

if test "x$cs_cc_compiler_known" != "xyes" ; then

  case "$host_os" in

    linux* | none)

      # IBM Blue Gene
      #--------------

      $CC -qversion 2>&1 | grep 'XL C' | grep 'Blue Gene' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is IBM XL C compiler for Blue Gene"

        # Version strings for logging purposes and known compiler flag
        $CC -qversion > $outfile 2>&1
        cs_ac_cc_version=`grep 'XL C' $outfile`
        cs_compiler_known=yes
        cs_linker_set=yes

        # Default compiler flags
        if test -d /bgl/BlueLight/ppcfloor/bglsys/include ; then
          cppflags_default="-I/opt/ibmmath/essl/4.2/include -I/bgl/BlueLight/ppcfloor/bglsys/include"
          cflags_default="-g -qmaxmem=-1 -qarch=440d -qtune=440"
          cflags_default_opt="-O3"
          cflags_default_hot="-O3 -qhot"
          cflags_default_dbg=""
        elif test -d /bgsys/drivers/ppcfloor/comm/include ; then
          cppflags_default="-I/opt/ibmmath/essl/4.4/include -I/bgsys/drivers/ppcfloor/comm/include"
          cflags_default="-g -qmaxmem=-1 -qarch=450d -qtune=450"
          cflags_default_opt="-O3"
          cflags_default_hot="-O3 -qhot"
          cflags_default_dbg=""
        else
          cppflags_default=""
          cflags_default=""
          cflags_default_opt="-O3"
          cflags_default_hot="-O3"
          cflags_default_dbg="-g"
        fi
        cflags_default_prf="-pg"
        cflags_default_omp="-qsmp=omp"

        # Default  linker flags
        ldflags_default=""
        ldflags_default_opt="-O3"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-pg"

      fi
      ;;

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
            cflags_default_omp="-omp"
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

    SUPER-UX*)

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
        cflags_default="-Kc99 -ftrace --pvctl,loopcnt=2147483647"
        cflags_default_opt=""
        cflags_default_dbg=""
        cflags_default_prf=""
        cflags_default_omp=""

        # Default linker flags
        ldflags_default="-ftrace -Wf,-pvctl,loopcnt=2147483647 -f90lib"
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
        cflags_default_opt="-O2 -woff 1429,1521"
        cflags_default_dbg="-g -woff 1429,1521,1209 -fullwarn"
        cflags_default_prf="-O0"

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
        cflags_default="-Aa +e +DA2.0W"
        cflags_default_opt="+O2"
        cflags_default_hot="+O3 +Oinline=Orient3D_split,Orient3D_normalize,Orient3D_set_maxvalue"
        cflags_default_dbg="-g"
        cflags_default_prf="-G"
        cflags_default_omp="+Oopenmp" # most pragmas require +O3

        # Default linker flags
        ldflags_default="+DA2.0W +FPVZOUD +U77"
        ldflags_default_opt="+O1"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-fbexe"
        cflags_default_omp="+Oopenmp"

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
        cflags_default="-Xa"
        cflags_default_opt="-xO2"
        cflags_default_dbg="-g"
        cflags_default_prf="-pg"
        cflags_default_omp="-xopenmp"

     fi
     ;;

    *)

      # Unknown
      #--------

      cflags_default=""
      cflags_default_opt="-O"
      cflags_default_dbg="-g"
      cflags_default_prf=""
      cflags_default_omp=""
      ;;

  esac

fi

if test "x$cflags_default_hot" = "x" ; then
  cflags_default_hot=$cflags_default
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
  fcflags_default="-x f95-cpp-input"
  fcflags_default_dbg="-g"
  fcflags_default_opt="-O"
  fcflags_default_prf="-pg"
  fcflags_default_omp="-fopenmp"

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
    fcflags_default="-cpp -fpic -warn"
    fcflags_default_dbg="-g -O0 -traceback -check all -fpe0 -ftrapuv"
    fcflags_default_opt="-O2"
    fcflags_default_hot="-O3"
    fcflags_default_prf="-p"
    fcflags_default_omp="-openmp"

    # Modify default flags on certain systems

    case "$host-os-$host_cpu" in
      *ia64)
        fcflags_default_opt="-O2 -mcpu=itanium2-p9000"
        fcflags_default_hot="-O3 -mcpu=itanium2-p9000"
        ;;
    esac

  fi

  case "$host_os" in

    linux* | none)

      # IBM Blue Gene
      #--------------

      $FC -qversion 2>&1 | grep 'XL Fortran' | grep 'Blue Gene' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$FC' is IBM XL Fortran compiler for Blue Gene"

        # Version strings for logging purposes and known compiler flag
        $FC -qversion > $outfile 2>&1
        cs_ac_fc_version=`grep 'XL Fortran' $outfile`
        cs_fc_compiler_known=yes

        # Default compiler flags
        if test -d /bgl/BlueLight/ppcfloor/bglsys/include ; then
          fcflags_default="-g -qmaxmem=-1 -qarch=440d -qtune=440 -qextname -qsuffix=cpp=f90"
          fcflags_default_opt="-O3"
          fcflags_default_hot="-O3 -qhot"
        elif test -d /bgsys/drivers/ppcfloor/comm/include ; then
          fcflags_default="-g -qmaxmem=-1 -qarch=450d -qtune=450 -qextname -qsuffix=cpp=f90"
          fcflags_default_opt="-O3"
          fcflags_default_hot="-O3 -qhot"
        else
          fcflags_default="-qsuffix=cpp=f90"
          fcflags_default_opt="-O3"
        fi
        fcflags_default_prf="-pg"
        fcflags_default_omp="-qsmp=omp"

      fi
      ;;

    SUPER-UX*)

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
        fcflags_default="-Ep -C hopt -ftrace -I. -Wf,-pvctl,loopcnt=2147483647"
        fcflags_default_opt=""
        fcflags_default_dbg=""
        fcflags_default_prf=""
        fcflags_default_omp=""

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
        fcflags_default_omp="+Oopenmp" # most pragmas require +O3

      fi
      ;;

    *)

      # Unknown
      #--------

      fcflags_default=""
      fcflags_default_opt="-O"
      fcflags_default_dbg="-g"
      fcflags_default_prf=""
      ;;

  esac

fi

if test "x$fcflags_default_hot" = "x" ; then
  fcflags_default_hot=$fcflags_default
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
      ldflags_default="-rdynamic"
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
      ldflags_default_prf=""
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

# Additional libraries and options for specific systems
#------------------------------------------------------

if test -d /bgl/BlueLight/ppcfloor/bglsys ; then #  For Blue Gene/L

  bg_sys_ldflags="-L/bgl/BlueLight/ppcfloor/bglsys/lib"
  bg_sys_libs="-lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts -lnss_files -lnss_dns -lresolv"
  bg_trace="/bgl/local/lib/libmpitrace.a"

  ldflags_default="${ldflags_default} -Wl,-allow-multiple-definition -L/opt/ibmcmp/xlmass/bg/4.3/blrts_lib -L/opt/ibmmath/essl/4.2/lib ${bg_sys_ldflags}"
  libs_default="-lmass -lmassv -lesslbg ${bg_trace} ${bg_sys_libs}"
  cs_disable_shared=yes

elif test -d /bgsys/drivers/ppcfloor/comm/include ; then #  For Blue Gene/P

  bg_sys_ldflags="-L/bgsys/drivers/ppcfloor/comm/lib -L/bgsys/drivers/ppcfloor/runtime/SPI"
  bg_sys_libs="-lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk -lSPI.cna -lrt -lpthread"
  bg_trace="/bgsys/local/tools_ibm/lib/libmpitrace.a"

  ldflags_default="${ldflags_default} -Wl,-allow-multiple-definition -L/opt/ibmcmp/xlmass/bg/4.4/bglib -L/opt/ibmmath/essl/4.4/lib ${bg_sys_ldflags}"
  libs_default="-lmass -lmassv -lesslsmpbg ${bg_trace} ${bg_sys_libs}"
  cs_disable_shared=yes

fi

# Finish

export LANG=$save_LANG

# Clean temporary files

rm -f conftest* a.out $outfile

