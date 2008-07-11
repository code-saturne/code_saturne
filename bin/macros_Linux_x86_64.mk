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
# Macros du Makefile Code_Saturne pour Linux_x86_64
###################################################
#
# Macro pour BFT
#---------------

BFT_HOME        =/home/saturne/opt/bft-1.0.6/arch/Linux_x86_64

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macro pour FVM
#---------------

FVM_HOME        =/home/saturne/opt/fvm-0.10.0/arch/Linux_x86_64

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macro pour MPI
#---------------

# Option MPI
MPI             =1
MPE             =0
MPE_COMM        =0

# Pour Open MPI sur saturne
MPI_HOME        =/home/saturne/opt/openmpi-1.2.6/arch/Linux_x86_64
MPI_INC         =-I$(MPI_HOME)/include
MPI_LIB         =-pthread -L$(MPI_HOME)/lib -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

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

XML_HOME =

XML_INC  =-I/usr/include/libxml2
XML_LIB  =-lxml2

# Macro pour BLAS
#----------------

# Option BLAS
BLAS            =1
BLAS_HOME       =
BLAS_INC        =
BLAS_CFLAGS     =-D_CS_HAVE_CBLAS
BLAS_LDFLAGS    =-lcblas -latlas

# Macro pour gettext
#-------------------

# Option gettext
NLS				=0


# Preprocesseur
#--------------

PREPROC         =
PREPROCFLAGS    =


# Compilateur C
#--------------

CCOMP                  = /home/saturne/opt/gcc-4.3.1/arch/Linux_x86_64/bin/gcc

CCOMPFLAGSDEF          = -std=c99 -funsigned-char -pedantic -W -Wall -Wshadow \
                         -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings \
                         -Wstrict-prototypes -Wmissing-prototypes \
                         -Wmissing-declarations -Wnested-externs -Wno-uninitialized 

CCOMPFLAGS             = $(CCOMPFLAGSDEF) -O -Wno-unused
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -O2
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -O2
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -O0
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -O0            
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g3            
CCOMPFLAGSPROF         = -pg
CCOMPFLAGSVERS         = -v            


# Compilateur FORTRAN 
#--------------------
#  Profiling gprof : -pg -a

FTNCOMP                = /home/saturne/opt/gcc-4.3.1/arch/Linux_x86_64/bin/gfortran

FTNCOMPFLAGSDEF        = -I.

FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) -O1
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) -O2
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) -O0
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) -O0
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g
FTNCOMPFLAGSPROF       = -pg
FTNCOMPFLAGSVERS       = -v

FTNPREPROCOPT          =

# Linker
#-------

# Linker

LDEDL           = $(FTNCOMP)
LDEDLFLAGS      = -O
LDEDLFLAGSLO    = -O0
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -pg
LDEDLFLAGSVERS  = -v
LDEDLRPATH      = -rdynamic -Wl,-rpath -Wl,/home/saturne/opt/gcc-4.3.1/arch/Linux_x86_64/lib64:

# If libxml2 must be installed under CS_ROOT, it may be necessary to add
# $(CS_ROOT)/opt/libxml2-2.6.19/arch/$NOM_ARCH/lib
# to LDEDLRPATH

# Options pour librairie partagee
#--------------------------------

BUILD_SO        =1
CCFLAGSSO       =-fPIC
FTNFLAGSSO      =-fPIC
LDEDLFLAGSSO    =-Wl,-soname -Wl,libcs14.so -fPIC -shared


# Positionnement des variables pour le pre-processeur
#----------------------------------------------------
#
# _POSIX_SOURCE          : utilisation des fonctions standard POSIX

VARDEF          = -D_POSIX_SOURCE


# Librairies a "linker"
#----------------------

# Librairies de base toujours prises en compte

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm -lpthread

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
#
# paquet 70% cpu promav gradrc prodsc
# paquet 10% cpu bilsc2 ;
#    prodsc est 4 fois plus rapide en O3 qu'en O2
#    bilsc2 plus rapide en O1
#    pour les autres, on privilegie l'O2, qui est suppose plus fiable
#      mais fait perdre un  peu de temps
#
#  Pour les fortrans, les listes ci-dessous servent a differencier
#	les options d'optimisation
#
#
#  Temporairement, gradmc en O1 pour eviter un bug d'optim potentiel
#       avec gcc 3.3.2 (resolu en 3.3.3)
#

LISTE_OPT_PART1 = gradrc.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 = gradmc.F

