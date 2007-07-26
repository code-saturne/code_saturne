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
# Macros du Makefile Code_Saturne pour Bull Novascale
#####################################################
#
# Macro pour BFT
#---------------

BFT_HOME       = /home/saturne/Saturne/opt/bft-1.0.4/arch/Linux_IA64
BFT_INC        = -I$(BFT_HOME)/include
BFT_LDFLAGS    = -L$(BFT_HOME)/lib -lbft

# Macro pour FVM
#---------------

FVM_HOME       = /home/saturne/Saturne/opt/fvm-0.8.0/arch/Linux_IA64
FVM_INC        = -I$(FVM_HOME)/include
FVM_LDFLAGS    = -L$(FVM_HOME)/lib -lfvm

# Macro pour MPI
#---------------

# Option MPI
MPI             =1
MPE             =0
MPE_COMM        =0
#
#Les bibliothèques MPI sont directement appelées par mpicc et mpif77
MPI_INC         =
MPI_LIB         =-lmpi -lrt -ldl -lpmi -lcpuset

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

XML_HOME = /home/saturne/Saturne/opt/libxml2-2.6.19

XML_INC  =-I$(XML_HOME)/include/libxml2
XML_LIB  =-L$(XML_HOME)/arch/Linux/lib -lxml2

# Macro pour BLAS
#----------------

# Option BLAS
BLAS            =1
BLAS_CFLAGS     =-D_CS_HAVE_MKL
BLAS_INC        =-I/applications/intel/cmkl/9.0.018/include
BLAS_LDFLAGS    =-L/applications/intel/cmkl/9.0.018/lib/64 -lmkl -lmkl_blacs_intelmpi20 -lmkl_ipf -lguide -lpthread 

# Preprocesseur
#--------------

PREPROC         =
PREPROCFLAGS    =

# Compilateur C
#--------------

CCOMP                  = icc

CCOMPFLAGSDEF          = -fpic -std=c99 -strict-ansi -Wall -Wcheck -Wmissing-prototypes \
			 -Wuninitialized -Wshadow -funsigned-char -Wpointer-arith \
			 -mtune=itanium2-p9000

CCOMPFLAGS             = $(CCOMPFLAGSDEF) -O2
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -O3
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -O3
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -O3
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -O1            
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g -O0 -traceback
CCOMPFLAGSPROF         = -p
CCOMPFLAGSVERS         = -V

# Compilateur FORTRAN
#--------------------
#  Profiling gprof : -p
#  Attention, par defaut on a -Mbounds avec mpif77 (l'inverse de mpicc)

FTNCOMP                = ifort

FTNCOMPFLAGSDEF        = -fpic -mtune=itanium2-p9000 -warn

FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) -O2
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) -O1
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g -O0 -traceback
FTNCOMPFLAGSPROF       = -p
FTNCOMPFLAGSVERS       = -V

FTNPREPROCOPT          =

# Linker
#-------

# Linker

LDEDL           = ifort
LDEDLFLAGS      = -O -nofor_main
LDEDLFLAGSLO    = -O0 -nofor_main
LDEDLFLAGSDBG   = -g -nofor_main
LDEDLFLAGSPROF  = -p
LDEDLFLAGSVERS  = -V
LDEDLRPATH      = -Wl,-rpath -Wl,


# Positionnement des variables pour le pre-processeur
#----------------------------------------------------
#
# _POSIX_SOURCE          : utilisation des fonctions standard POSIX
#
VARDEF          = -D_POSIX_SOURCE


# Librairies a "linker"
#----------------------

# Librairies de base toujours prises en compte

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm

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
#
# paquet 70% cpu promav gradrc gradco prodsc
# paquet 10% cpu jacobi prcpol bilsc2 ;
#
#  Pour les fortrans, les listes ci-dessous servent a differencier
#	les options d'optimisation
#

LISTE_OPT_PART1 = gradco.F gradrc.F jacobi.F prcpol.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

