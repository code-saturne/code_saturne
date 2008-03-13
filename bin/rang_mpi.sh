#!/bin/sh
#============================================================================
#
#                    Code_Saturne version 1.3
#                    ------------------------
#
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
# Détection a priori du rang MPI d'un processus
#==============================================
#
# Pour récupérer le rang MPI depuis un script lancé
# par mpirun (ou prun, mpijob, ou autre équivalent) :
#
# MPI_RANK=`$CS_HOME/bin/rang_mpi.sh $@`
#
# Utile surtout pour lancer des applications
# MPMD de type couplage avec des environnements
# MPI 1.2 ne proposant pas l'équivalent de la
# commande mpiexec (ou des appschemas LAM)

# Pour Cluster HP AlphaServer (Chrome au CCRT)
if [ "$RMS_RANK" != "" ] ; then
  MPI_RANK="$RMS_RANK"

# Pour Cluster Opteron ou Itanium sous Linux avec Slurm (Argent et Tantale au CCRT)
elif [ "$SLURM_PROCID" != "" ] ; then
  MPI_RANK="$SLURM_PROCID"

# Pour Cluster Opteron sous Linux MPICH-GM (Cluster Chatou)
elif [ "$GMPI_ID" != "" ] ; then
  MPI_RANK="$GMPI_ID"

# Pour Cluster MFTT sous Linux avec MPI Scali (réseau SCI)
elif [ "$SCAMPI_PROCESS_PARAM" != "" ] ; then
  MPI_RANK=`echo $SCAMPI_PROCESS_PARAM | cut -f4 -d' '`

# Pour MPICH2 (qui fournit aussi la commande mpiexec)
elif [ "$PMI_RANK" != "" ] ; then
  MPI_RANK="$PMI_RANK" 

# Pour Open MPI 1.0 (qui fournit aussi la commande mpiexec)
elif [ "$OMPI_MCA_ns_nds_vpid" != "" ] ; then
  MPI_RANK="$OMPI_MCA_ns_nds_vpid"

# Pour LAM 7.1 (où l'on peut aussi utiliser un appschema)
elif [ "$LAMRANK" != "" ] ; then
  MPI_RANK="$LAMRANK"

# Pour machine IBM sous MVAPICH (d'après test Computing on Demand)
elif [ "$SMPIRUN_RANK" != "" ] ; then
  MPI_RANK="$MPIRUN_RANK"

# Pour Cluster utilisant HP-MPI (ancien système Tantale au CCRT)
elif [ "$MPI_PROC" != "" ] ; then
  MPI_RANK=`echo $MPI_PROC | cut -f 2 -d,`

# Pour MPICH 1.2 "standard" (avec communication chp4 "usuelle")
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

# Fin des cas connus
fi

# Affichage du rang obtenu

echo "$MPI_RANK"

