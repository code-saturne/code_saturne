#!/bin/sh
#============================================================================
#
#     This file is part of the Code_Saturne Kernel, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2008 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne Kernel is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne Kernel is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#============================================================================
#
# A priori detection of the rank of a MPI process
#================================================
#
# In order to get the MPI rank from a script launched
# by mpirun (or prun, mpijob, or equivalent):
#
# MPI_RANK=`$CS_HOME/bin/mpi_rank.sh $@`
#
# Mainly useful to launch MPMD applications
# like coupling within MPI 1.2 environment
# which does not have a command like mpiexec
# (or some appschemes as LAM)

# On HP AlphaServer cluster (CCRT: Chrome)
if [ "$RMS_RANK" != "" ] ; then
  MPI_RANK="$RMS_RANK"

# On Opteron or Itanium cluster under Linux with Slurm (CCRT: Platine, Tantale)
elif [ "$SLURM_PROCID" != "" ] ; then
  MPI_RANK="$SLURM_PROCID"

# On Opteron cluster under Linux MPICH-GM (EDF Chatou cluster)
elif [ "$GMPI_ID" != "" ] ; then
  MPI_RANK="$GMPI_ID"

# On EDF MFTT cluster under Linux with MPI Scali (SCI network)
elif [ "$SCAMPI_PROCESS_PARAM" != "" ] ; then
  MPI_RANK=`echo $SCAMPI_PROCESS_PARAM | cut -f4 -d' '`

# For MPICH2 (it also provides the mpiexec command)
elif [ "$PMI_RANK" != "" ] ; then
  MPI_RANK="$PMI_RANK"

# For Open MPI (it also provides the mpiexec command)
elif [ "$OMPI_MCA_ns_nds_vpid" != "" ] ; then
  MPI_RANK="$OMPI_MCA_ns_nds_vpid"

# For LAM 7.1 (an appscheme can also be used)
elif [ "$LAMRANK" != "" ] ; then
  MPI_RANK="$LAMRANK"

# On IBM machine with MVAPICH (following a test "Computing on Demand")
elif [ "$SMPIRUN_RANK" != "" ] ; then
  MPI_RANK="$MPIRUN_RANK"

# On cluster with HP-MPI (CCRT: Tantale, old system)
elif [ "$MPI_PROC" != "" ] ; then
  MPI_RANK=`echo $MPI_PROC | cut -f 2 -d,`

# For "standard" MPICH 1.2 (with "usual" chp4 communication)
else
  MPI_RANK=0
  next_arg_is_rank=""
  for arg in "$@" ; do
    if [ "$arg" = "-p4rmrank" ] ; then
      next_arg_is_rank="1"
    elif [ "$next_arg_is_rank" = "1" ] ; then
      MPI_RANK="$arg"
      next_arg_is_rank="0"
    fi
  done

# End of known cases
fi

# Output of the obtained rank

echo "$MPI_RANK"

