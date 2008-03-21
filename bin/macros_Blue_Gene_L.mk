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
# Macros du Makefile Code_Saturne pour Blue Gene
################################################
#
# Chemins système
#----------------

BGL_SYS  = /bgl/BlueLight/ppcfloor/bglsys

# Macro pour BFT
#---------------

BFT_HOME        =/gpfs2/home/saturne/opt/bft-1.0.6/arch/bgl

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macro pour FVM
#---------------

FVM_HOME        =/gpfs2/home/saturne/opt/fvm-0.10.0/arch/bgl

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macro pour MPI
#---------------

# Option MPI
MPI             =1
MPE             =0
MPE_COMM        =0

# Pour MPI BlueGene
MPI_HOME        =
MPI_INC         = -I$(BGL_SYS)/include
MPI_LIB  = 

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

XML_HOME = /gpfs2/home/saturne/opt/libxml2-2.6.19

XML_INC  =-I$(XML_HOME)/include/libxml2
XML_LIB  =-L$(XML_HOME)/arch/bgl/lib -lxml2

# Macro pour BLAS
#----------------

# Option BLAS
BLAS            =1
ESSL            =1 # librairie ESSL IBM avec extension BLAS
BLAS_INC        =-I/opt/ibmmath/essl/4.2/include
BLAS_CFLAGS     =-D_CS_HAVE_ESSL
BLAS_LDFLAGS    =


# Preprocesseur
#--------------

PREPROC         =
PREPROCFLAGS    =


# Compilateur C natif
#--------------------

CCOMP                  = blrts_xlc

CCOMPFLAGSDEF          = -g -qmaxmem=-1 -qarch=440d -qtune=440
#CCOMPFLAGSDEF          = -g -qmaxmem=-1 -qarch=440d -qtune=440 -qflttrap=enable:overflow:zerodivide -qsigtrap=xl_trcedump
#CCOMPFLAGSDEF          = -g -qmaxmem=-1 -qarch=440d -qtune=440 -qsource -qlist

CCOMPFLAGS             = $(CCOMPFLAGSDEF) -O3
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -O3 -qhot
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -O3 -qhot
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -O3 -qhot
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -O0
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g
CCOMPFLAGSPROF         = -pg
CCOMPFLAGSVERS         = -v


# Compilateur FORTRAN
#--------------------
#  Profiling gprof : -pg -a

FTNCOMP                = blrts_xlf

FTNCOMPFLAGSDEF        = -g -qmaxmem=-1 -qarch=440d -qtune=440 -qextname
#FTNCOMPFLAGSDEF        = -g -qmaxmem=-1 -qarch=440d -qtune=440 -qextname -qflttrap=enable:overflow:zerodivide -qsigtrap=xl_trcedump
#FTNCOMPFLAGSDEF        = -g -qmaxmem=-1 -qarch=440d -qtune=440 -qextname -qsource -qlist

FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF) -O3
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF) -O3 -qhot
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF) -O3 -qhot
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF) -O3 -qhot
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF) -O0
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF) -g
FTNCOMPFLAGSPROF       = -pg
FTNCOMPFLAGSVERS       = -v

FTNPREPROCOPT          = -WF,

# Linker
#-------

# Linker

LDEDL           = blrts_xlf -qflttrap=enable:overflow:zerodivide
#LDEDL          = blrts_xlf -qflttrap=enable:overflow:zerodivide -qsigtrap=xl_trcedump
LDEDLFLAGS      = -O3
LDEDLFLAGSLO    = -O0
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -pg
LDEDLFLAGSVERS  = -v
LDEDLRPATH      =


# Positionnement des variables pour le pre-processeur
#----------------------------------------------------
#
# _POSIX_SOURCE          : utilisation des fonctions standard POSIX

VARDEF          = -D_POSIX_SOURCE


# Librairies a "linker"
#----------------------

# Zlib utilisee par HDF5
ZLIB     = -L/gpfs2/home/saturne/opt/zlib-1.2.1/arch/bgl/lib -lz

# Librairies IBM
MASS     = -L/opt/opt/ibmcmp/xlmass/bg/4.3/blrts_lib -lmass -lmassv
LIBMAT   = /bgl/local/lib/libmpitrace.a
ESSL     = /opt/ibmmath/essl/4.2/lib/libesslbg.a
EXIT     = /bgl/local/lib/libexit.a

# Librairies de base toujours prises en compte

LIBBASIC = $(ZLIB)\
-Wl,-allow-multiple-definition $(MASS) $(ESSL) $(LIBMAT) -L$(BGL_SYS)/lib -lmpich.rts -lmsglayer.rts -lrts.rts -ldevices.rts -lnss_files -lnss_dns -lresolv

# Librairies en mode sans option

LIBOPT   =

# Librairies en mode optimisation reduite

LIBLO    =

# Librairies en mode DEBUG

LIBDBG   =

# Librairie en mode ElectricFence (malloc debugger)

LIBEF    =

# Liste eventuelle des fichiers a compiler avec des options particulieres
#------------------------------------------------------------------------

# Sous la forme :
# LISTE_OPT_PART = fic_1.c fic_2.c \
#                fic_3.F
#
# paquet 70% cpu promav gradrc gradco prodsc
# paquet 10% cpu jacobi prcpol bilsc2 ;
#    option -qhot recommande pour ces sous-programmes
#    pour les autres, on privilegie l'O3, qui est suppose plus fiable
#      mais fait perdre un  peu de temps
#
#  Pour les fortrans, les listes ci-dessous servent a differencier
#	les options d'optimisation
#

LISTE_OPT_PART1 = gradco.F gradrc.F jacobi.F prcpol.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

