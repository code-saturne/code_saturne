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
# Macros du Makefile Code_Saturne pour HP Proliant DL585
########################################################
#
# Macros pour BFT
#----------------

BFT_HOME       = /home/cont002/saturne/opt/bft-1.0.6/arch/Linux
BFT_INC        = -I$(BFT_HOME)/include
BFT_LDFLAGS    = -L$(BFT_HOME)/lib -lbft

# Macros pour FVM
#----------------

FVM_HOME       = /home/cont002/saturne/opt/fvm-0.10.0/arch/Linux
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

XML_HOME = /usr

XML_INC  =-I$(XML_HOME)/include/libxml2
XML_LIB  =-L$(XML_HOME)/lib -lxml2

# Macro pour BLAS
#----------------

# Option BLAS
BLAS            =1

BLAS_INC        =-I/applications/atlas/include
BLAS_CFLAGS     =-D_CS_HAVE_CBLAS
BLAS_LDFLAGS    =-L/applications/atlas/lib -lcblas -latlas -lg2c

# Preprocesseur
#--------------

PREPROC         =
PREPROCFLAGS    =


# Compilateur
#------------

CCOMP                  = mpicc
#CCOMPFLAGSDEF          = -Xa -Ktrap=fp   Bug compilateur PGI 6.2 si -Ktrap en meme temps que fastsse
CCOMPFLAGSDEF          = -Xa
CCOMPFLAGS             = $(CCOMPFLAGSDEF) -O1 
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -O2 -fast -fastsse  
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -O2 -fast -fastsse
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -O1
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -O0            
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g -Mbounds
CCOMPFLAGSPROF         = -pg
CCOMPFLAGSVERS         = -V



# Compilateur FORTRAN
#--------------------
#  Profiling gprof : -pg
#  Attention, par defaut on a -Mbounds avec mpif77 (l'inverse de mpicc)

FTNCOMP                = mpif77
#FTNCOMPFLAGSDEF        = -Ktrap=fp -fastsse     Bug compilateur PGI 6.2 si -Ktrap en meme temps que fastsse
FTNCOMPFLAGSDEF        = -fastsse
FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) -O -Mnobounds
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) -O2 -Mnobounds
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) -O2 -Mnobounds
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) -O -Mnobounds
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) -O -Mnobounds
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g -O0
FTNCOMPFLAGSPROF       = -pg
FTNCOMPFLAGSVERS       = -V

FTNPREPROCOPT          =

# Linker
#-------

# Linker

LDEDL           = mpif77
LDEDLFLAGS      = -O -Mnomain
LDEDLFLAGSLO    = -O0 -Mnomain
LDEDLFLAGSDBG   = -g -Mnomain
LDEDLFLAGSPROF  = -pg
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
#    prodsc est 4 fois plus rapide en O6 qu'en O2
#    bilsc2 plus rapide en O1
#    pour les autres, on privilegie l'O2, qui est suppose plus fiable
#      mais fait perdre un  peu de temps (2% de perte par rapport a 
#      gradco O3, gradrc jacobi prcpol promav O5) 
#
#  Pour les fortrans, les listes ci-dessous servent a differencier
#	les options d'optimisation
#

LISTE_OPT_PART1 = gradco.F gradrc.F jacobi.F prcpol.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

