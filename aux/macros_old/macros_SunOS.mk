#============================================================================
#
#                    Code_Saturne version 1.3
#                    ------------------------
#
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
# Macros for Makefile under Solaris
###################################
#
# Macros for BFT
#---------------

BFT_HOME        =/home/saturne/opt/bft-1.0.7/arch/SunOS

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macros for FVM
#---------------

FVM_HOME        =/home/saturne/opt/fvm-0.11.0/arch/SunOS

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macros for MPI
#---------------

# MPI support
MPI             =1
MPE             =0
MPE_COMM        =0

# For Open MPI on saturne
MPI_HOME        =/home/saturne/opt/openmpi-1.2.7/arch/SunOS
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

# Option XML
XML             =1

XML_HOME =

XML_INC  =-I/usr/include/libxml2
XML_LIB  =-lxml2

# Macro pour BLAS
#----------------

# BLAS support
BLAS            =0
BLAS_HOME       =/home/saturne/opt/atlas-3.8.2/arch/SunOS
BLAS_INC        =-I$(BLAS_HOME)/include
BLAS_CFLAGS     =-D_CS_HAVE_CBLAS
BLAS_LDFLAGS    =-L$(BLAS_HOME)/lib -lcblas -latlas

# Macros for gettext
#-------------------

# gettext support
NLS				=0

# Set CS_LANG to FR to have French translation
CS_LANG         =FR


# Preprocessor
#-------------

PREPROC         =
PREPROCFLAGS    =


# C compiler
#-----------

CCOMP                  = $(CC)

CCOMPFLAGSDEF          = -Xa # -Xc99 may be use on up-to-date compiler versions

CCOMPFLAGS             = $(CCOMPFLAGSDEF) -xO2
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -xO3 -xinline=Orient3D_split
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -xO2
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -xO2
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -xO0
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g
CCOMPFLAGSPROF         = -p
CCOMPFLAGSVERS         = -V


# Fortran compiler
#-----------------
#  Profiling gprof : -pg -a

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
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -G
LDEDLFLAGSVERS  = -V
LDEDLRPATH      = -R


# Set preprocessor variables
#---------------------------
#
# _POSIX_SOURCE          : POSIX standard functions

VARDEF          = -D_POSIX_SOURCE

# Libraries to link against
#--------------------------

# Base libraries (always used)

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) -lm -lmalloc 

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
# 70% cpu promav gradrc gradco prodsc
# 10% cpu jacobi prcpol bilsc2 ;
#    option O3 recommended for these subroutines
#    for others, we prefer O2, less risky, but slightly slower
#
# The file lists below correspond to different optimization levels
#

LISTE_OPT_PART1 = gradco.F gradrc.F jacobi.F prcpol.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 = gradmc.F

