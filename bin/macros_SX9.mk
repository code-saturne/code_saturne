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
# Macros for Makefile under NEC SX9
###################################
#
# Macros for BFT
#---------------

BFT_HOME        =/home/saturne/opt/bft-1.0.8/arch/SX9

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macros for FVM
#---------------

FVM_HOME        =/home/saturne/opt/fvm-0.13.0/arch/SX9

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
SOCKET          =0
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
BLAS            =1
BLAS_HOME       =
BLAS_INC        =-I/usr/include
BLAS_CFLAGS     =-DHAVE_CBLAS
BLAS_LDFLAGS    =-lcblas -latlas

# Macros for dlopen
#------------------

# Dynamic library loader support

DLOPEN          =0


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

CCOMP                  = sxmpicc

CCOMPFLAGSDEF          = -Kc99 -ftrace -D__uxpvp__ -pvctl,loopcnt=2147483647

CCOMPFLAGS             = $(CCOMPFLAGSDEF)
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF)
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF)
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF)
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF)
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF)
CCOMPFLAGSPROF         = 
CCOMPFLAGSVERS         = 

# Fortran compiler
#-----------------


FTNCOMP                = sxmpif90

FTNCOMPFLAGSDEF        = -Ep -C hopt -ftrace -D__uxpvp__ -I. -Wf,-pvctl,loopcnt=2147483647

FTNCOMPFLAGS           = $(FTNCOMPFLAGSDEF)
FTNCOMPFLAGSOPTPART1   = $(FTNCOMPFLAGSDEF)
FTNCOMPFLAGSOPTPART2   = $(FTNCOMPFLAGSDEF)
FTNCOMPFLAGSOPTPART3   = $(FTNCOMPFLAGSDEF)
FTNCOMPFLAGSLO         = $(FTNCOMPFLAGSDEF)
FTNCOMPFLAGSDBG        = $(FTNCOMPFLAGSDEF)
FTNCOMPFLAGSPROF       = 
FTNCOMPFLAGSVERS       =

FTNPREPROCOPT          =

# Linker
#-------

# Linker

LDEDL           = $(FTNCOMP)
LDEDLFLAGS      = -ftrace
LDEDLFLAGSLO    =
LDEDLFLAGSDBG   = 
LDEDLFLAGSPROF  =
LDEDLFLAGSVERS  =
LDEDLRPATH      = -Wf,-pvctl,loopcnt=2147483647


# Archiver
#---------

AR              = sxar
ARFLAGS         = cr


# Set preprocessor variables
#---------------------------
#
# _POSIX_SOURCE          : POSIX standard functions

VARDEF          = 

# Libraries to link against
#--------------------------

# Base libraries (always used)

LIBBASIC = $(BFT_LDFLAGS) $(FVM_LDFLAGS) 

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
#                fic_3.f90
#
# 70% cpu promav gradrc prodsc
# 10% cpu bilsc2 ;
#    option O3 recommended for these subroutines
#    for others, we prefer O2, less risky, but slightly slower
#
# The file lists below correspond to different optimization levels
#

LISTE_OPT_PART1 = gradmc.f90 gradrc.f90 promav.f90 cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.f90 prods2.f90 prods3.f90 cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

