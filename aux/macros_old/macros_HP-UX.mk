#============================================================================
#
#                    Code_Saturne version 1.3
#                    ------------------------
#
#
#     This file is part of the Code_Saturne Kernel, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2007 EDF S.A., France
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
# Macros du Makefile Code_Saturne pour HP-UX
############################################
#
# Macros pour BFT
#----------------

BFT_HOME       = /home/saturne/opt/bft-1.0.5/arch/HP-UX
BFT_INC        = -I$(BFT_HOME)/include
BFT_LDFLAGS    = -L$(BFT_HOME)/lib -lbft -L/home/saturne/opt/zlib-1.2.1/arch/HP-UX/lib -lz

# Macro pour FVM
#---------------

FVM_HOME        =/home/saturne/opt/fvm-0.8.0/arch/HP-UX

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macro pour MPI
#---------------

# Option MPI
MPI             =0
MPE             =0
MPE_COMM        =1
MPI_INC         =
MPI_LIB         =


# Macro pour Sockets
#-------------------

# Option Socket
SOCKET          =1
SOCKET_INC      =
SOCKET_LIB      =

# Macro pour XML
#---------------

# Option XML
XML             =1

XML_HOME = /home/saturne/opt/libxml2-2.6.19

XML_INC  =-I$(XML_HOME)/include/libxml2
XML_LIB  =-L$(XML_HOME)/arch/HP-UX/lib -lxml2

# Macro pour BLAS
#----------------

# Option BLAS
BLAS            =0
BLAS_INC        =
BLAS_CFLAGS     =
BLAS_LDFLAGS    =

# Preprocesseur
#--------------

PREPROC         =
PREPROCFLAGS    =


# Compilateur C natif
#---------------------
# Option `-Ae' pour extended ANSI (tremai.c)
# Option `+e' pour les `long long' de HDF
# O2 poutot que O3 (option par defaut)

CCOMP                  = /opt/ansic/bin/cc
CCOMPFLAGSDEF          = -Ae +e +DA2.0W
CCOMPFLAGS             = $(CCOMPFLAGSDEF) +O2
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) +O3 +Oinline=Orient3D_split,Orient3D_normalize,Orient3D_set_maxvalue
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) +O2
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) +O2
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) +O0
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g
CCOMPFLAGSPROF         = -G
CCOMPFLAGSVERS         = -V


# Compilateur FORTRAN
#--------------------
# Option `+FPVZOUD' pour les arrets en cas d'exception
# O2 (defaut) est plus rapide que O3 (et a priori avec moins de risques)

FTNCOMP                = f90
FTNCOMPFLAGSDEF        = +fp_exception +FPVZOUD +U77 +DA2.0W +noppu
FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) +O1
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) +O2
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) +O2
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) +O2
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) +O0
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g
FTNCOMPFLAGSPROF       = -G
FTNCOMPFLAGSVERS       = +version

FTNPREPROCOPT          =

# Linker
#-------

# Linker

LDEDL           = f90
LDEDLFLAGS      = +O1 +FPVZOUD +U77 +DA2.0W
LDEDLFLAGSLO    = +O0 +FPVZOUD +U77 +DA2.0W
LDEDLFLAGSDBG   = -g  +FPVZOUD +U77 +DA2.0W
LDEDLFLAGSPROF  = -G
LDEDLFLAGSVERS  =
LDEDLRPATH      =


# Positionnement des variables pour le pre-processeur
#----------------------------------------------------
#
# _POSIX_SOURCE          : utilisation des fonctions standard POSIX

VARDEF          = -D_POSIX_SOURCE


# Librairies a "linker"
#----------------------

# Librairies de base toujours prises en compte

LIBBASIC = $(FVM_LDFLAGS) $(BFT_LDFLAGS) -lm -lF90 -L/opt/fortran90/lib -lU77

# Librairies en mode sans option

LIBOPT   =

# Librairies en mode optimisation reduite

LIBLO    =

# Librairies en mode DEBUG

LIBDBG   =

# Librairie en mode ElectricFence (malloc debugger)

LIBEF    =-L/home/saturne/opt/ElectricFence-2.0.5/lib/HP-UX -lefence

# Liste eventuelle des fichiers a compiler avec des options particulieres
#------------------------------------------------------------------------

# Sous la forme :
# LISTE_OPT_PART = fic_1.c fic_2.c \
#                fic_3.F
# Obligatoirement O1 ou inferieur : clptur.F condli.F lagent.F
#
# paquet 70% cpu promav gradrc gradco prodsc
# paquet 10% cpu jacobi prcpol bilsc2 ;
#
#  Pour le fichier c, il s'agit de permettre l'inline (non ansi)

LISTE_OPT_PART1 = bilsc2.F gradco.F gradrc.F jacobi.F prcpol.F prodsc.F prods2.F prods3.F promav.F cs_lagrang.c cs_matrix.c cs_sles.c cs_blas.c cs_benchmark.c
LISTE_OPT_PART2 =
LISTE_OPT_PART3 =

