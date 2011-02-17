#============================================================================
#
#                    Code_Saturne version 1.3
#                    ------------------------
#
#
#     This file is part of the Code_Saturne Kernel, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
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
# Macros for Makefile under Darwin
##################################
#
# Macros for BFT
#---------------

BFT_HOME        =/Users/saturne/Code_Saturne/1.3/opt/bft-1.1/arch/Darwin

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macros for FVM
#---------------

FVM_HOME        =/Users/saturne/Code_Saturne/1.3/opt/fvm-0.15/arch/Darwin

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macros for MPI
#---------------

# MPI support
MPI             =1
MPE             =0
MPE_COMM        =0

# We use MPI wrappers, so nothing to add here
MPI_HOME        =
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
XML             =0

XML_HOME =

XML_INC  =-I/usr/include/libxml2
XML_LIB  =-lxml2

# Macros for BLAS
#----------------

# BLAS support
BLAS            =0
BLAS_HOME       =
BLAS_INC        =-I/usr/include
BLAS_CFLAGS     =-D_CS_HAVE_CBLAS
BLAS_LDFLAGS    =-lcblas -latlas

# Macros for gettext
#-------------------

# gettext support
NLS             =0

# Set CS_LANG to FR to have French translation
CS_LANG         =


# C compiler
#-----------

CCOMP                  = mpicc

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


# Fortran compiler
#-----------------

FTNCOMP                = mpif77

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

LDEDL           = $(FTNCOMP)
LDEDLFLAGS      = -O
LDEDLFLAGSLO    = -O0
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -pg
LDEDLFLAGSVERS  = -v
LDEDLRPATH      = -rdynamic -Wl,-rpath -Wl,:


# Set preprocessor variables
#---------------------------

VARDEF          = -D_POSIX_SOURCE

# Libraries to link against
#--------------------------

# Base libraries (always used)

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm -lpthread -lz

# Libraries in production mode

LIBOPT   =

# Libraries in reduced optimization mode

LIBLO    =

# Libraries in DEBUG mode

LIBDBG   =

# Library in ElectricFence (malloc debugger) mode

LIBEF    =

# Optional lists of files to compile with specific options
#---------------------------------------------------------

# In the form:
# LISTE_OPT_PART = fic_1.c fic_2.c \
#                fic_3.F
#
# 70% cpu promav gradrc gradco prodsc
# 10% cpu jacobi prcpol bilsc2 ;
#    option O3 recommended for these subroutines
#    for others, we prefer O2, less risky, but slightly slower
#
# The file lists below correspond to different optimization levels
#

LISTE_OPT_PART1 = gradco.F gradrc.F jacobi.F prcpol.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

