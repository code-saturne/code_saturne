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
# Macros for Makefile under HP Proliant DL585
#############################################
#
# Macros for BFT
#---------------

BFT_HOME        =/home/cont002/saturne/opt/bft-1.0.7/arch/Linux

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macros for FVM
#---------------

FVM_HOME        =/home/cont002/saturne/opt/fvm-0.11.0/arch/Linux

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macros for MPI
#---------------

# MPI support
MPI             =1
MPE             =0
MPE_COMM        =0

#Les bibliothèques MPI sont directement appelées par mpicc et mpif77
MPI_BIN         =
MPI_INC         =
MPI_LIB         =


# Macros for Sockets
#-------------------

# Sockets support
SOCKET          =1
SOCKET_INC      =
SOCKET_LIB      =

# Macros for XML
#---------------

# XML support
XML             =1

XML_HOME =

XML_INC  =-I/usr/include/libxml2
XML_LIB  =-lxml2

# Macros for BLAS
#----------------

# BLAS support
BLAS            =1
BLAS_HOME       =/applications/atlas
BLAS_INC        =-I$(BLAS_HOME)/include
BLAS_CFLAGS     =-D_CS_HAVE_CBLAS
BLAS_LDFLAGS    =-L$(BLAS_HOME)/lib -lcblas -latlas -lg2c

# Macros for gettext
#-------------------

# gettext support
NLS             =0

# Set CS_LANG to FR to have French translation
CS_LANG         =


# Preprocessor
#-------------

PREPROC         =
PREPROCFLAGS    =


# C compiler
#-----------

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


# Fortran compiler
#-----------------
#  Profiling gprof : -pg
#  Attention, par defaut on a -Mbounds avec mpif77 (l'inverse de mpicc)

FTNCOMP                = mpif77

FTNCOMPFLAGSDEF        = -fastsse
#FTNCOMPFLAGSDEF        = -Ktrap=fp -fastsse     Bug compilateur PGI 6.2 si -Ktrap en meme temps que fastsse

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


# Set preprocessor variables
#---------------------------
#
# _POSIX_SOURCE          : POSIX standard functions
#
VARDEF          = -D_POSIX_SOURCE

# Libraries to link against
#--------------------------

# Base libraries (always used)

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm

# Libraries in production mode

LIBOPT   =

# Libraries in reduced optimization mode

LIBLO    =

# Libraries in DEBUG mode

LIBDBG   =

# Library in ElectricFence (malloc debugger) mode

LIBEF    =-lefence

# Optional lists of files to compile with specific options
#---------------------------------------------------------

# In the form:
# LISTE_OPT_PART = fic_1.c fic_2.c \
#                fic_3.F
#
# 70% cpu promav gradrc prodsc
# 10% cpu bilsc2 ;
#    option O3 recommended for these subroutines
#    for others, we prefer O2, less risky, but slightly slower
#
# The file lists below correspond to different optimization levels
#

LISTE_OPT_PART1 = gradrc.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

