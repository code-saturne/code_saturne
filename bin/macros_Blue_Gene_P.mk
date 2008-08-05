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
# Macros for Makefile under Blue Gene/P
#######################################
#
# System paths
#-------------

BGL_SYS  = /bgsys/drivers/ppcfloor/comm

# Macros for BFT
#---------------

BFT_HOME        =/gpfs/home/saturne/opt/bft-1.0.6/arch/bgp

BFT_INC         =-I$(BFT_HOME)/include
BFT_LDFLAGS     =-L$(BFT_HOME)/lib -lbft

# Macros for FVM
#---------------

FVM_HOME        =/gpfs/home/saturne/opt/fvm-0.10.0/arch/bgp

FVM_INC         =-I$(FVM_HOME)/include
FVM_LDFLAGS     =-L$(FVM_HOME)/lib -lfvm

# Macros for MPI
#---------------

# MPI support
MPI             =1
MPE             =0
MPE_COMM        =0

# For Blue Gene MPI
MPI_HOME        =
MPI_INC         = -I$(BGL_SYS)/include
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
XML             =1

XML_HOME = /gpfs/home/saturne/opt/libxml2-2.3.32/arch/bgp

XML_INC  =-I$(XML_HOME)/include/libxml2
XML_LIB  =-L$(XML_HOME)/lib -lxml2

# Macros for BLAS
#----------------

# BLAS support
BLAS            =1
ESSL            =1 # IBM ESSL library with BLAS extension
BLAS_INC        =-I/opt/ibmmath/essl/4.3/include
BLAS_CFLAGS     =-D_CS_HAVE_ESSL
BLAS_LDFLAGS    =

# Macros for gettext
#-------------------

# gettext support
NLS                                =0


# Preprocessor
#-------------

PREPROC         =
PREPROCFLAGS    =


# C compiler
#-----------

CCOMP                  = bgxlc

CCOMPFLAGSDEF          =

#CCOMPFLAGSDEF          = -g -qmaxmem=-1 -qarch=450d -qtune=450
#CCOMPFLAGSDEF          = -g -qmaxmem=-1 -qarch=450d -qtune=450 -qflttrap=enable:overflow:zerodivide -qsigtrap=xl_trcedump
#CCOMPFLAGSDEF          = -g -qmaxmem=-1 -qarch=450d -qtune=450 -qsource -qlist

CCOMPFLAGS             = $(CCOMPFLAGSDEF) -O3
CCOMPFLAGSOPTPART1     = $(CCOMPFLAGSDEF) -O3 -qhot
CCOMPFLAGSOPTPART2     = $(CCOMPFLAGSDEF) -O3 -qhot
CCOMPFLAGSOPTPART3     = $(CCOMPFLAGSDEF) -O3 -qhot
CCOMPFLAGSLO           = $(CCOMPFLAGSDEF) -O0
CCOMPFLAGSDBG          = $(CCOMPFLAGSDEF) -g
CCOMPFLAGSPROF         = -pg
CCOMPFLAGSVERS         = -v


# Fortran compiler
#-----------------
#  Profiling gprof : -pg -a

FTNCOMP                = bgxlf

FTNCOMPFLAGSDEF        = -qextname

#FTNCOMPFLAGSDEF        = -g -qmaxmem=-1 -qarch=450d -qtune=450 -qextname
#FTNCOMPFLAGSDEF        = -g -qmaxmem=-1 -qarch=450d -qtune=450 -qextname -qflttrap=enable:overflow:zerodivide -qsigtrap=xl_trcedump
#FTNCOMPFLAGSDEF        = -g -qmaxmem=-1 -qarch=450d -qtune=450 -qextname -qsource -qlist

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

LDEDL           = bgxlf_r
#LDEDL           = bgxlf_r -qflttrap=enable:overflow:zerodivide -qsigtrap=xl_trcedump
LDEDLFLAGS      = -O3
LDEDLFLAGSLO    = -O0
LDEDLFLAGSDBG   = -g
LDEDLFLAGSPROF  = -pg
LDEDLFLAGSVERS  = -v
LDEDLRPATH      =


# Set preprocessor variables
#---------------------------
#
# _POSIX_SOURCE          : POSIX standard functions

VARDEF          = -D_POSIX_SOURCE

# Libraries to link against
#--------------------------

# Zlib used by HDF5
ZLIB     = -L/bgsys/local/tools_ibm/lib -lz

# IBM libraries
SYSLIBS  = -L/bgsys/drivers/ppcfloor/comm/lib -lmpich.cnk -ldcmfcoll.cnk -ldcmf.cnk -L/bgsys/drivers/ppcfloor/runtime/SPI -lSPI.cna -lrt -lpthread
MASS     = -L/opt/ibmcmp/xlmass/bg/4.4/bglib -lmass -lmassv
ESSL     = -L/opt/ibmmath/essl/4.3/lib -lesslbg -lesslsmpbg
TRACE    = /bgsys/local/tools_ibm/lib/libmpitrace.a

# Base libraries (always used)

LIBBASIC = $(ZLIB) -Wl,--allow-multiple-definition $(MASS) $(ESSL) $(TRACE) $(SYSLIBS)

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
#    option -qhot recommended for these subroutines
#    for others, we prefer O3, less risky, but slightly slower
#
# The file lists below correspond to different optimization levels
#

LISTE_OPT_PART1 = gradco.F gradrc.F jacobi.F prcpol.F promav.F cs_matrix.c cs_sles.c
LISTE_OPT_PART2 = prodsc.F prods2.F prods3.F cs_blas.c cs_benchmark.c
LISTE_OPT_PART3 =

