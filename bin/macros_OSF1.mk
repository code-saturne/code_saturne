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
# Macros du Makefile Code_Saturne pour HP SC
############################################
#
# Macros pour BFT
#----------------

BFT_HOME       = /home/saturne/Saturne/opt/bft-1.0.4/arch/OSF1
BFT_INC        = -I$(BFT_HOME)/include
BFT_LDFLAGS    = -L$(BFT_HOME)/lib -lbft

# Macro pour FVM
#---------------

FVM_HOME        =/home/saturne/Saturne/opt/fvm-0.8.0/arch/OSF1

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macro pour MPI
#---------------

# Option MPI
MPI             =1
MPE             =0
MPE_COMM        =1
#
# Pour scali
MPI_INC         =-I/usr/include
MPI_LIB         =-L/usr/lib -lmpi


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
XML_LIB  =-L$(XML_HOME)/arch/OSF1/lib -lxml2

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


# Compilateur C natif on evite de trapper les underflows ? -fptm u
#--------------------

CCOMP                  = $(CC)
CCOMPFLAGSDEF          = -ansi_alias -std \
                         -trapuv \
                         -check -msg_enable alignment -msg_enable noansi \
                         -msg_enable performance -portable -msg_enable c_to_cxx
CCOMPFLAGS             = $(CCOMPFLAGSDEF) -O -arch host -tune host   
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -O            
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -O  
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -O  
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -O0            
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g -check_bounds
CCOMPFLAGSPROF         = -pg
CCOMPFLAGSVERS         = -V



# Compilateur FORTRAN on evite de trapper les underflows ? -check underflow
#--------------------
#  Profiling gprof : -pg -a

FTNCOMP                = f90
FTNCOMPFLAGSDEF        = -warn declarations -std 
FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) -O -arch host -tune host 
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) -O
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) -O
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) -O
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) -O
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g -O0 -check bounds
FTNCOMPFLAGSPROF       = -pg
FTNCOMPFLAGSVERS       = -version

FTNPREPROCOPT          =

# Linker
#-------

# Linker

LDEDL           = $(CC)
LDEDLFLAGS      = -O
LDEDLFLAGSLO    = -O0
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -pg
LDEDLFLAGSVERS  = -V
LDEDLRPATH      = -Wl,-call_shared -Wl,-rpath -Wl,


# Positionnement des variables pour le pre-processeur
#----------------------------------------------------
#
# _POSIX_SOURCE          : utilisation des fonctions standard POSIX
#
# Sur OSF1, on utilise de préférence _OSF_SOURCE ; ceci entraîne
# la définition de _POSIX_SOURCE et de _XOPEN_SOURCE (mais pas de
# _X_OPEN_SOURCE_EXTENDED), mais ces macros ne doivent pas être définies
# directement dans ce cas : les fichiers include s'en chargent mais
# doivent positionner certaines variables dans un ordre bien défini.

VARDEF          = -D_OSF_SOURCE


# Librairies a "linker"
#----------------------

# Librairies de base toujours prises en compte

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm -lUfor -lfor

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

