dnl----------------------------------------------------------------------------
dnl   This file is part of the Code_Saturne Kernel, element of the
dnl   Code_Saturne CFD tool.
dnl
dnl   Copyright (C) 2009 EDF S.A., France
dnl
dnl   The Code_Saturne Kernel is free software; you can redistribute it
dnl   and/or modify it under the terms of the GNU General Public License
dnl   as published by the Free Software Foundation; either version 2 of
dnl   the License, or (at your option) any later version.
dnl
dnl   The Code_Saturne Kernel is distributed in the hope that it will be
dnl   useful, but WITHOUT ANY WARRANTY; without even the implied warranty
dnl   of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl   GNU General Public License for more details.
dnl
dnl   You should have received a copy of the GNU General Public Licence
dnl   along with the Code_Saturne Preprocessor; if not, write to the
dnl   Free Software Foundation, Inc.,
dnl   51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
dnl-----------------------------------------------------------------------------

# CS_AC_CONFIG_INFO_INIT([OUTPUT FILE NAME])
#-----------------------
# Write main config info file header.

AC_DEFUN([CS_AC_CONFIG_INFO_INIT], [

# First arg is output file name
if test "$1" = "" ; then
  cs_ac_config_info="cs-config"
else
  cs_ac_config_info=$1
fi

outfile="$cs_ac_config_info"-tmp

rm -f $outfile

AC_MSG_NOTICE([initializing $cs_ac_config_info])
])

# CS_AC_CONFIG_INFO_EXTRA([extra config info strings])
#------------------------
# Write extra info to config info file header.

AC_DEFUN([CS_AC_CONFIG_INFO_EXTRA], [

AC_REQUIRE([CS_AC_CONFIG_INFO_INIT])dnl

outfile="$cs_ac_config_info"-tmp

echo "# Configuration:" >> $outfile
echo "# --------------"    >> $outfile
echo "$1" >> $outfile

])

# CS_AC_CONFIG_INFO_CC([version], [version_full])
#---------------------
# Write available compiler info to config info file header.

AC_DEFUN([CS_AC_CONFIG_INFO_CC], [

AC_REQUIRE([CS_AC_CONFIG_INFO_INIT])dnl

outfile="$cs_ac_config_info"-tmp

if test "$1" != "" -o "$2" != "" ; then
  echo "" >> $outfile
  echo "# C compiler used for build: $1" >> $outfile
  echo "# --------------------------"    >> $outfile
  if test "$2" != "" ; then
    echo "$2" | sed 's/^/# /' >> $outfile
  fi
else
  AC_MSG_NOTICE([C compiler version info unavailable for configuration file])
fi
])

# CS_AC_CONFIG_INFO_FC([version], [version_full])
#---------------------
# Write available compiler info to config info file header.

AC_DEFUN([CS_AC_CONFIG_INFO_FC], [

AC_REQUIRE([CS_AC_CONFIG_INFO_INIT])dnl

outfile="$cs_ac_config_info"-tmp

if test "$1" != "" -o "$2" != "" ; then
  echo "" >> $outfile
  echo "# Fortran compiler used for build: $1" >> $outfile
  echo "# --------------------------------"    >> $outfile
  if test "$2" != "" ; then
    echo "$2" | sed 's/^/# /' >> $outfile
  fi
else
  AC_MSG_NOTICE([Fortran compiler version info unavailable for configuration file])
fi
])

# CS_AC_CONFIG_INFO_FINALIZE
#---------------------------
# Write main config info file body (script part).

AC_DEFUN([CS_AC_CONFIG_INFO_FINALIZE], [

AC_REQUIRE([CS_AC_CONFIG_INFO_INIT])dnl

outfile="$cs_ac_config_info"-tmp

AC_MSG_NOTICE([closing $cs_ac_config_info])

diff $outfile $cs_ac_config_info > /dev/null 2>&1
if test $? -eq 0 ; then
  AC_MSG_NOTICE([$cs_ac_config_info is unchanged])
  rm -f $outfile
else
  mv $outfile $cs_ac_config_info
  chmod +x $cs_ac_config_info
fi

])dnl

