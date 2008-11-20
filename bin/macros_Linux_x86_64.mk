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
# Macros for Makefile under Linux x86_64
########################################
#
# Macros for BFT
#---------------

BFT_HOME        =/home/saturne/opt/bft-1.0.8/arch/Linux_x86_64

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macros for FVM
#---------------

FVM_HOME        =/home/saturne/opt/fvm-0.12.0/arch/Linux_x86_64

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macros for MPI
#---------------

# MPI support
MPI             =1
MPE             =0
MPE_COMM        =0

# For Open MPI on saturne
MPI_HOME        =/home/saturne/opt/openmpi-1.2.6/arch/Linux_x86_64
MPI_BIN         =$(MPI_HOME)/bin
MPI_INC         =-I$(MPI_HOME)/include
MPI_LIB         =-pthread -L$(MPI_HOME)/lib -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl

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


# Preprocessor
#-------------

PREPROC         =
PREPROCFLAGS    =


# C compiler
#-----------

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


# Fortran compiler
#-----------------
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

# Shared library options
#-----------------------

BUILD_SO        =1
CCFLAGSSO       =-fPIC
FTNFLAGSSO      =-fPIC
LDEDLFLAGSSO    =-Wl,-soname -Wl,libcs14.so -fPIC -shared


# Set preprocessor variables
#---------------------------
#
# _POSIX_SOURCE          : POSIX standard functions

VARDEF          = -D_POSIX_SOURCE

# Libraries to link against
#--------------------------

# Base libraries (always used)

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm -lpthread

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
#  Temporarily, gradmc is compiled with O1 to bypass a potential optimization bug
#       with gcc 3.3.2 (resolved with 3.3.3)
#

LISTE_OPT_PART1 = gradmc.F gradrc.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

