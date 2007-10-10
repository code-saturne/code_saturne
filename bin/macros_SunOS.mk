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
# Macros du Makefile Code_Saturne pour SunOS
############################################
#
# Macros pour BFT
#----------------

BFT_HOME       = /home/saturne/opt/bft-1.0.5/arch/SunOS
BFT_INC        = -I$(BFT_HOME)/include
BFT_LDFLAGS    = -L$(BFT_HOME)/lib -lbft

# Macro pour FVM
#---------------

FVM_HOME        =/home/saturne/opt/fvm-0.8.0/arch/SunOS

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
SOCKET          =0
SOCKET_INC      =
SOCKET_LIB      =

# Macro pour XML
#---------------

# Option XML
XML             =1

XML_HOME = /home/saturne/opt/libxml2-2.6.19

XML_INC  =-I$(XML_HOME)/include/libxml2
XML_LIB  =-L$(XML_HOME)/arch/SunOS/lib -lxml2

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
#--------------------
# WorkShop Compilers 5.0 98/12/15 C 5.0

CCOMP                  = $(CC)
CCOMPFLAGSDEF          = -Xa
CCOMPFLAGS             = $(CCOMPFLAGSDEF) -xO2
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -xO3 -xinline=Orient3D_split
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -xO2
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -xO2
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -xO0
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g
CCOMPFLAGSPROF         = -p
CCOMPFLAGSVERS         = -V



# Compilateur FORTRAN
#--------------------

FTNCOMP                = f77
FTNCOMPFLAGSDEF        = 
FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) -O2
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) -O0
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g
FTNCOMPFLAGSPROF       = -G
FTNCOMPFLAGSVERS       = -V

FTNPREPROCOPT          =

# Linker
#-------

# Linker


LDEDL           = f77
LDEDLFLAGS      = -O
LDEDLFLAGSLO    = -O
LDEDLFLAGSEF    = -g
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -G
LDEDLFLAGSVERS  = -V
LDEDLRPATH      = -R


# Positionnement des variables pour le pre-processeur
#----------------------------------------------------
#
# _POSIX_SOURCE          : utilisation des fonctions standard POSIX

VARDEF          = -D_POSIX_SOURCE


# Librairies a "linker"
#----------------------

# Librairies de base toujours prises en compte

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm -lmalloc 

# Librairies en mode sans option

LIBOPT   =

# Librairies en mode optimisation reduite

LIBLO    =

# Librairies en mode DEBUG

LIBDBG   =

# Librairie en mode ElectricFence (malloc debugger)

LIBEF    =-lefence

# Liste eventuelle des fichiers a compiler avec des options particulieres
#------------------------------------------------------------------------

# Sous la forme :
# LISTE_OPT_PART = fic_1.c fic_2.c \
#                fic_3.F
#  Pour le fichier c, il s'agit de permettre l'inline (non ansi)

LISTE_OPT_PART1 = cs_lagrang.c cs_matrix.c cs_sles.c cs_blas.c cs_benchmark.c
LISTE_OPT_PART2 =
LISTE_OPT_PART3 =
