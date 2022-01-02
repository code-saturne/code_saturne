dnl Copyright (C) 2004-2022 EDF
dnl
dnl This file is part of the PLE software package.  For license
dnl information, see the COPYING file in the top level directory of the
dnl PLE source distribution.

# PLE_AC_CHECK_SIZEOF(TYPE, [PREFIX])
#------------------------------------
# get type sizes
# Optionnaly, the corresponding SIZEOF definition may be prefixed
# by the PREFIX variable

AC_DEFUN([PLE_AC_CHECK_SIZEOF],[

if test "$1" = "" ; then
  AC_MSG_ERROR([configure test cannot be run])
fi

AC_REQUIRE([PLE_AC_CONFIG_PUBL_INIT])dnl

ple_ac_lcname=`echo "$1" | sed y/' *'/'_p'/`
ple_ac_lower='abcdefghijklmnopqrstuvwxyz'
ple_ac_upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
if test "$2" != "" ; then
  ple_ac_szname_prefix=`echo $2 | sed y/$ple_ac_lower/$ple_ac_upper/`_
else
  ple_ac_szname_prefix=""
fi
ple_ac_szname_postfix=`echo $ple_ac_lcname | sed y/$ple_ac_lower/$ple_ac_upper/`
ple_ac_szname="${ple_ac_szname_prefix}SIZEOF_${ple_ac_szname_postfix}"
unset ple_ac_lower
unset ple_ac_upper
unset ple_ac_szname_prefix
unset ple_ac_szname_postfix

AC_CHECK_SIZEOF($1)
eval ple_ac_sizeof=\$ac_cv_sizeof_$ple_ac_lcname
if test "$ple_ac_sizeof" != "" -a "$ple_ac_sizeof" != "0"; then
  PLE_AC_CONFIG_PUBL_DEFINE([$ple_ac_szname], [$ple_ac_sizeof],
                            [The size of a '$1', as computed by sizeof.])
else
  PLE_AC_CONFIG_PUBL_SET([$ple_ac_szname], [no],
                         [The size of a '$1', as computed by sizeof.])
fi

unset ple_ac_lcname
unset ple_ac_szname
unset ple_ac_sizeof

/bin/rm -f conftest*])dnl
