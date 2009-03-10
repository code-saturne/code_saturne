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
#
# cflags_default         # Base CFLAGS                       (default: "")
# cflags_default_dbg     # Added to $CFLAGS for debugging    (default: "-g")
# cflags_default_opt     # Added to $CFLAGS for optimization (default: "-O")
# cflags_default_prf     # Added to $CFLAGS for profiling    (default: "")
#
# fcflags_default        # Base FCFLAGS                       (default: "")
# fcflags_default_dbg    # Added to $FCFLAGS for debugging    (default: "-g")
# fcflags_default_opt    # Added to $FCFLAGS for optimization (default: "-O")
# fcflags_default_prf    # Added to $FCFLAGS for profiling    (default: "")
#
# ldflags_default        # Base LDFLAGS                       (default: "")
# ldflags_default_dbg    # Added to $LDFLAGS for debugging    (default: "-g")
# ldflags_default_opt    # Added to $LDFLAGS for optimization (default: "-O")
# ldflags_default_prf    # Added to $LDFLAGS for profiling    (default: "")

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
  cflags_default_dbg="-g -O0 -traceback"
  cflags_default_opt="-O3"
  cflags_default_prf="-p"

  # Modify default flags on certain systems

  case "$host-os-$host_cpu" in
    *ia64)
      cflags_default_opt="-O3 -mcpu=itanium2-p9000"
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

        # Default compiler flags
        case "$host_cpu" in
          alphaev6|alphaev67|alphaev68|alphaev7)
            cflags_default="-arch host -tune host -ansi_alias -std -check_bounds -trapuv -check -msg_enable alignment -msg_enable noansi -msg_enable performance -portable -msg_enable c_to_cxx"
            cflags_default_opt="-O"
            cflags_default_dbg="-g"
            cflags_default_prf="-pg"
          ;;
        esac

      fi
      ;;

    uxpv*)

      # Native Fujitsu vectorizing C compiler (tested on VPP5000)
      #---------------------------------------

      $CC -V 2>&1 | grep 'Fujitsu' > /dev/null
      if test "$?" = "0" ; then

        echo "compiler '$CC' is Fujitsu compiler"

        # Version strings for logging purposes and known compiler flag
        $CC -V conftest.c > $outfile 2>&1
        cs_ac_cc_version=`grep ccom $outfile`
        cs_cc_compiler_known=yes
        cs_linker_set=yes

        # Default compiler flags
        cflags_default="-KA64 -Kvp"
        cflags_default_opt="-O"
        cflags_default_dbg="-g -Kargchk -w4"
        cflags_default_prf="-pg"

        # Default linker flags
        ldflags_default="-Kvp -Kargchk -Wl,-S1:d"
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
        cflags_default_dbg="-g"
        cflags_default_prf="-G"

        # Default linker flags
        ldflags_default="+DA2.0W"
        ldflags_default_opt="+O2"
        ldflags_default_dbg="-g"
        ldflags_default_prf="-fbexe"

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

cs_gfortran=gfortran

if test "x$cs_gfortran" = "xgfortran"; then

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

fi

# Compiler still not identified
#------------------------------

if test "x$cs_fc_compiler_known" != "xyes" ; then

  case "$host_os" in

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

# Finish

export LANG=$save_LANG

# Clean temporary files

rm -f conftest* a.out $outfile

