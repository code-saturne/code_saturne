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
# Makefile for Code_Saturne
###########################
#
# Environment Macros
# ==================

SHELL           =/bin/sh

VERSION         ="1.4.0"

TERM            =xterm



# Project configuration macros
# ============================

# Base project directories

DIR_CS         = $(CS_HOME)
BIN_CS         = $(DIR_CS)/bin
INC_CS         = $(DIR_CS)/include
SRC_CS         = $(DIR_CS)/src
USR_CS         = $(DIR_CS)/users
ARCH_CS        = $(DIR_CS)/arch/$(NOM_ARCH)
ARCHBIN_CS     = $(ARCH_CS)/bin
ARCHINC_CS     = $(ARCH_CS)/include
ARCHLIB_CS     = $(ARCH_CS)/lib
ARCHLOCALE_CS  = $(ARCH_CS)/share/locale
ARCHUSR_CS     = $(ARCH_CS)/share/users
BLD_CS         = $(DIR_CS)/build

# Names of the base library and Code_Saturne executable file

NOM_LIB_CS_PRE = "saturne"
NOM_EXE        = "cs14.exe"
NOM_PACKAGE    = "ncs"


# OS-specific macros
# ==================

include $(BIN_CS)/macros_$(NOM_ARCH).mk



# User macros
#============

# Definition of the chosen mode (EF, DBG ...)

MYOPT           = $(OPT)

# Definition of variables for compilation

MYCCOMPFLAGS    = $(OPTC)
MYFTNCOMPFLAGS  = $(OPTF)

# Definition of variables for linking

MYLDEDLFLAGS    = $(OPTLD)
MYLIBS          = $(OPTLIBS)


# Macros defining libraries to link
#===================================

LIBRARYOPT    = $(LIBOPT)

LIBRARYLO     = $(LIBLO)

LIBRARYDBG    = $(LIBDBG)

LIBRARYEF     = $(LIBEF)

LIBRARY       = $(LIBBASIC)



# Presentation macros
#--------------------

BO = `tput bold`
SO = `tput smso`
NO = `tput rmso; tput sgr0`



# Definition of objects to create
#================================

.SUFFIXES:
.SUFFIXES: .c .f90 .F90 .o

.c.o:
	@echo $<
	$(CCOMP)   $(INC_CS_I) $(ALLCCOMPFLAGS) -c $<

.f90.o:
	@echo $<
	$(FTNCOMP) $(ALLFCOMPFLAGS) -c $<

.F90.o:
	@echo $<
	$(FTNCOMP) $(ALLFCOMPFLAGS) -c $<


# Target definitions
#===================


default:
	@dirtrav=`pwd` ;\
	 if   [ `echo "$${dirtrav}" | grep "$(BLD_CS)"` ] ;  \
	 then make lib ;\
	 else make link ;\
	 fi


saturne: link


lib: privtitre privconfig privcompil privlib privfin


lib_install: privtitre privconfig privcompil privlib privlibinst privincinst privusrinst privfin


link: privtitre privconfig privcompil privlink privfin


compversion: privtitre privcompil privlink privfin


comp: privtitre privcomp privfin


privcomp: privconfig privcompil


#
# echo_mpi_bin must be called using "make -f", so it must not include recursive
#              calls to make, as they would not find the Makefile.

echo_mpi_bin:
	@echo $(MPI_BIN)

#
# privcompil recognizes the working directory (module and OPT),
#            makes a symbolic link to the source files (in case of build),
#            generates the dependency lists,
#            and sets OBJECTS and compiler options.
#

privcompil:
	@dirtrav=`pwd` ;\
	 dirprec=`dirname $${dirtrav}` ;\
	 if   [ "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	 then opttrav="" ;\
	      modtrav=`basename $${dirtrav}` ;\
	 else dirprec1=`dirname $${dirprec}` ;\
	      if   [ "$${dirprec1}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	      then opttrav=`basename $${dirtrav}` ;\
	           modtrav=`basename $${dirprec}` ;\
	      else opttrav="$(MYOPT)" ;\
	           modtrav="" ;\
	      fi ;\
	 fi ;\
	 if   [ "$${modtrav}" = "BASE" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/atmo $(INC_CS)/cfbl $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/ctwr $(INC_CS)/elec $(INC_CS)/fuel $(INC_CS)/lagr $(INC_CS)/pprt $(INC_CS)/rayt $(INC_CS)/mati" ;\
	      src_cs_all="$(USR_CS)/base $(SRC_CS)/base $(USR_CS)/pprt $(SRC_CS)/pprt" ;\
	 elif [ "$${modtrav}" = "ATMO" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/pprt $(INC_CS)/atmo" ;\
	      src_cs_all="$(USR_CS)/atmo $(SRC_CS)/atmo" ;\
	 elif [ "$${modtrav}" = "CFBL" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/pprt $(INC_CS)/cfbl" ;\
	      src_cs_all="$(USR_CS)/cfbl $(SRC_CS)/cfbl" ;\
	 elif [ "$${modtrav}" = "COGZ" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/pprt $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/rayt" ;\
	      src_cs_all="$(USR_CS)/cogz $(SRC_CS)/cogz" ;\
	 elif [ "$${modtrav}" = "CPLV" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/lagr $(INC_CS)/pprt $(INC_CS)/rayt";\
	      src_cs_all="$(USR_CS)/cplv $(SRC_CS)/cplv" ;\
	 elif [ "$${modtrav}" = "CTWR" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/pprt $(INC_CS)/ctwr";\
	      src_cs_all="$(USR_CS)/ctwr $(SRC_CS)/ctwr" ;\
	 elif [ "$${modtrav}" = "FUEL" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/fuel $(INC_CS)/lagr $(INC_CS)/pprt $(INC_CS)/rayt";\
	      src_cs_all="$(USR_CS)/fuel $(SRC_CS)/fuel" ;\
	 elif [ "$${modtrav}" = "ELEC" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/elec $(INC_CS)/pprt $(INC_CS)/rayt" ;\
	      src_cs_all="$(USR_CS)/elec $(SRC_CS)/elec" ;\
	 elif [ "$${modtrav}" = "LAGR" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/lagr $(INC_CS)/pprt $(INC_CS)/rayt" ;\
	      src_cs_all="$(USR_CS)/lagr $(SRC_CS)/lagr" ;\
	 elif [ "$${modtrav}" = "RAYT" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/cfbl $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/fuel $(INC_CS)/lagr $(INC_CS)/pprt $(INC_CS)/rayt" ;\
	      src_cs_all="$(USR_CS)/rayt $(SRC_CS)/rayt" ;\
	 elif [ "$${modtrav}" = "MATI" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/cfbl $(INC_CS)/cogz $(INC_CS)/cplv $(INC_CS)/lagr $(INC_CS)/mati $(INC_CS)/pprt $(INC_CS)/rayt" ;\
	      src_cs_all="$(SRC_CS)/mati" ;\
	 elif [ "$${modtrav}" = "" ] ;\
	 then inc_cs_all="$(ARCHINC_CS)" ;\
	      src_cs_all="" ;\
	 fi ;\
	 if   [ "$${modtrav}" != "LINK" -a "$${modtrav}" != "" ] ;\
	 then dirtrav=`pwd` ;\
	      dirprec=`dirname $${dirtrav}` ;\
	      dirprec1=`dirname $${dirprec}` ;\
	      if   [   "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)" \
	           -o "$${dirprec1}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	      then for dir in `echo $${inc_cs_all} $${src_cs_all}` ;\
	           do \
		        for fich in $${dir}/*.* ;\
		        do \
		             if [ ! -d "$${fich}" -a -f $${fich} ] ;\
			     then fichloc=`basename $${fich}` ;\
			          if [ ! -h "$${fichloc}" ] ;\
			          then ln -s $${fich} . ;\
			          fi ;\
			     fi ;\
		        done ;\
	           done ;\
	      fi ;\
	 fi ;\
	 make privdepend ;\
	 echo "***  Compilation" ;\
	 cflags="$(VARDEF)" ;\
	 dir_inc_cs="$(BFT_INC) $(FVM_INC)" ;\
	 if [ "$${modtrav}" = "" ]  ;\
     then cflags="$${cflags} -I$(ARCHINC_CS)" ;\
          fflags="$${fflags} -I$(ARCHINC_CS)" ;\
	 fi ; \
	 if   [ "${MPI}" -eq "1" ] ;\
	 then dir_inc_cs="$${dir_inc_cs} $(MPI_INC)" ;\
	      cflags="$${cflags} -DHAVE_MPI " ;\
	 fi ; \
       if [ "${MPE}" -eq "1" ] ;\
       then dir_inc_cs="$${dir_inc_cs} $(MPE_INC)" ;\
            cflags="$${cflags} -DHAVE_MPE ";\
	   if [ "${MPE_COMM}" -eq "1" ] ; \
               then cflags="$${cflags} -D_CS_COMM_PROFILING ";\
         fi ; \
       fi ; \
	 if   [ "${SOCKET}" -eq "1" ] ;\
	 then dir_inc_cs="$${dir_inc_cs} $(SOCKET_INC)" ;\
	      cflags="$${cflags} -DHAVE_SOCKET " ;\
	 fi ; \
	 if   [ "${BLAS}" -eq "1" ] ;\
	 then dir_inc_cs="$${dir_inc_cs} $(BLAS_INC)" ;\
	      cflags="$${cflags} $(BLAS_CFLAGS) " ;\
	 fi ; \
	 if   [ "${XML}" -eq "1" ] ;\
	 then dir_inc_cs="$${dir_inc_cs} $(XML_INC)" ;\
	      cflags="$${cflags} -DHAVE_XML " ;\
	 fi ; \
	 if   [ "${NLS}" -eq "1" ] ;\
	 then cflags="$${cflags} -DENABLE_NLS -DPACKAGE=\\\"$(NOM_PACKAGE)\\\" -DLOCALEDIR=\\\"$(ARCHLOCALE_CS)\\\"" ;\
	 fi ; \
	 if   [ "${DLOPEN}" -eq "1" ] ;\
	 then cflags="$${cflags} -DHAVE_DLOPEN " ;\
	 fi ; \
	 if   [ "${CS_LANG}" = "FR" ] ;\
	 then fflags="$${fflags} $(FTNPREPROCOPT)-D_CS_LANG_FR" ;\
	 fi ; \
	 if   [ "$(MYOPT)" = "LO" ] ;\
	 then make privdependobj OBJECTS="`cat .liste_mkfiles.o`" \
	                         INC_CS_I="$${dir_inc_cs}" \
	           ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DDEBUG   $(CCOMPFLAGSLO)"  \
	           ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGSLO)" ;\
	 elif [ "$(MYOPT)" = "DBG" -o "$(MYOPT)" = "EF" ] ;\
	 then make privdependobj OBJECTS="`cat .liste_mkfiles.o`" \
	                         INC_CS_I="$${dir_inc_cs}" \
	           ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DDEBUG   $(CCOMPFLAGSDBG)"  \
	           ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGSDBG)"  ;\
	 elif [ "$(MYOPT)" = "PROF" ] ;\
	 then \
	      if   [ -n "$(LISTE_OPT_PART1)" ] ;\
	      then make privdependoptpart LISTE_OPT_PART="$(LISTE_OPT_PART1)" ;\
	           make privdependobj OBJECTS="`cat .liste_mkoptpart.o`" \
				 INC_CS_I="$${dir_inc_cs}" \
		   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DNDEBUG  $(CCOMPFLAGSOPTPART1) $(CCOMPFLAGSPROF)" \
		   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGSOPTPART1) $(FTNCOMPFLAGSPROF)"  ;\
	      fi ;\
	      if   [ -n "$(LISTE_OPT_PART2)" ] ;\
	      then make privdependoptpart LISTE_OPT_PART="$(LISTE_OPT_PART2)" ;\
	           make privdependobj OBJECTS="`cat .liste_mkoptpart.o`" \
				 INC_CS_I="$${dir_inc_cs}" \
		   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DNDEBUG  $(CCOMPFLAGSOPTPART2) $(CCOMPFLAGSPROF)" \
		   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGSOPTPART2) $(FTNCOMPFLAGSPROF)"  ;\
	      fi ;\
	      if   [ -n "$(LISTE_OPT_PART3)" ] ;\
	      then make privdependoptpart LISTE_OPT_PART="$(LISTE_OPT_PART3)" ;\
	           make privdependobj OBJECTS="`cat .liste_mkoptpart.o`" \
				 INC_CS_I="$${dir_inc_cs}" \
		   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DNDEBUG  $(CCOMPFLAGSOPTPART3) $(CCOMPFLAGSPROF)" \
		   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGSOPTPART3) $(FTNCOMPFLAGSPROF)"  ;\
	      fi ;\
              make privdependobj OBJECTS="`cat .liste_mkfiles.o`" \
                                 INC_CS_I="$${dir_inc_cs}" \
                   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DNDEBUG  $(CCOMPFLAGS) $(CCOMPFLAGSPROF)" \
                   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGS) $(FTNCOMPFLAGSPROF)"  ;\
	 elif [ "$(MYOPT)" = "VERS" ] ;\
	 then make privdependobj OBJECTS="`cat .liste_mkfiles.o`" \
	                         INC_CS_I="$${dir_inc_cs}" \
	           ALLCCOMPFLAGS=" $${cflags}  \
				    -DNDEBUG  $(CCOMPFLAGSVERS)"  \
	           ALLFCOMPFLAGS=" $${fflags}  \
				              $(FTNCOMPFLAGSVERS)" ;\
	 elif [ "$(MYOPT)" = "" ] ;\
	 then if   [ -n "$(LISTE_OPT_PART1)" ] ;\
	      then make privdependoptpart LISTE_OPT_PART="$(LISTE_OPT_PART1)" ;\
	           make privdependobj OBJECTS="`cat .liste_mkoptpart.o`" \
                         INC_CS_I="$${dir_inc_cs}" \
		   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
			            -DNDEBUG  $(CCOMPFLAGSOPTPART1)" \
		   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGSOPTPART1)"  ;\
	      fi ;\
	      if   [ -n "$(LISTE_OPT_PART2)" ] ;\
	      then make privdependoptpart LISTE_OPT_PART="$(LISTE_OPT_PART2)" ;\
	           make privdependobj OBJECTS="`cat .liste_mkoptpart.o`" \
                         INC_CS_I="$${dir_inc_cs}" \
		   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				   -DNDEBUG  $(CCOMPFLAGSOPTPART2)" \
		   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				         $(FTNCOMPFLAGSOPTPART2)"  ;\
	      fi ;\
	      if   [ -n "$(LISTE_OPT_PART3)" ] ;\
	      then make privdependoptpart LISTE_OPT_PART="$(LISTE_OPT_PART3)" ;\
	           make privdependobj OBJECTS="`cat .liste_mkoptpart.o`" \
				 INC_CS_I="$${dir_inc_cs}" \
		   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
			     	    -DNDEBUG  $(CCOMPFLAGSOPTPART3)" \
		   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				         $(FTNCOMPFLAGSOPTPART3)"  ;\
	      fi ;\
	      make privdependobj OBJECTS="`cat .liste_mkfiles.o`" \
                                 INC_CS_I="$${dir_inc_cs}" \
                   ALLCCOMPFLAGS=" $${cflags} $(MYCCOMPFLAGS) \
				    -DNDEBUG  $(CCOMPFLAGS)" \
                   ALLFCOMPFLAGS=" $${fflags} $(MYFTNCOMPFLAGS) \
				              $(FTNCOMPFLAGS)"  ;\
	 fi
	@echo "***  Compilation finished"   ; echo

#
# privdepend genere les listes de dependance
#

privdepend:
	@mkfiles="" ;\
	 liste="" ;\
	for file in *.c *.f90 *.F90 ;\
	   do if [ "$${file}" != "*.c" -a "$${file}" != "*.f90" -a "$${file}" != "*.F90" ] ;\
	      then mkfiles="$${mkfiles} `basename $${file}`" ;\
	      fi ;\
	   done ;\
	   echo "$${mkfiles}" | sed -e "s/\.c/.o/g" \
	                            -e "s/\.f90/.o/g" \
	                            -e "s/\.F90/.o/g" \
                                  -e "s/ \*\.o//g"> .liste_mkfiles.o
#
# privdependoptpart genere les listes de dependance

privdependoptpart:
	@mkfiles="" ;\
	 if   [ -n "$(LISTE_OPT_PART)" ] ;\
	 then \
	      for file in $(LISTE_OPT_PART) ;\
	      do \
	           if [ -n "$${file}" ] ;\
	           then if [ -r "$${file}" ] ;\
		          then mkfiles="$${mkfiles} `basename $${file}`" ;\
	                fi ;\
		     fi ;\
	      done ;\
	 fi ;\
	 if [ -f .liste_mkoptpart.o ] ;\
	 then \rm -f .liste_mkoptpart.o ;\
	 fi ;\
	 echo "$${mkfiles}" | sed -e "s/\.c/.o/g" \
	                          -e "s/\.f90/.o/g" \
	                          -e "s/\.F90/.o/g" > .liste_mkoptpart.o


privdependobj: $(OBJECTS)


privtitre:
	@echo
	@echo "***  Code_Saturne $(VERSION) compilation/link"
	@echo "***  System: $(NOM_ARCH)"
	@echo


#
# privlib generates the library 
#           its name is based on the working directory
#

privlib:
	@echo "***  Creating library"
	@dirtrav=`pwd` ;\
	 dirprec=`dirname $${dirtrav}` ;\
	 if   [ "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	 then opttrav="" ;\
	      modtrav=`basename $${dirtrav}` ;\
	 else dirprec1=`dirname $${dirprec}` ;\
	      if   [ "$${dirprec1}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	      then opttrav=`basename $${dirtrav}` ;\
	           modtrav=`basename $${dirprec}` ;\
	      else opttrav="$(MYOPT)" ;\
	           modtrav="" ;\
	      fi ;\
	 fi ;\
	 \rm -f "lib$(NOM_LIB_CS_PRE)$${modtrav}$${opttrav}.a" ;\
	 $(AR) $(ARFLAGS) "lib$(NOM_LIB_CS_PRE)$${modtrav}$${opttrav}.a" *.o 
	@echo "***  Library creation finished" ; echo

#
# privlibinst installs the library
#           its name is based on the working directory
#

privlibinst:
	@echo "***  Installing library"
	@if [ ! -d "$(DIR_CS)/arch" ] ; \
	 then mkdir "$(DIR_CS)/arch" ; \
	 fi ;\
         if [ ! -d "$(ARCH_CS)" ] ; \
	 then mkdir "$(ARCH_CS)" ; \
	 fi ;\
         if [ ! -d "$(ARCHBIN_CS)" ] ; \
	 then mkdir "$(ARCHBIN_CS)" ; \
	 fi ;\
         if [ ! -d "$(ARCHLIB_CS)" ] ; \
	 then mkdir "$(ARCHLIB_CS)" ; \
	 fi 
	@dirtrav=`pwd` ;\
	 dirprec=`dirname $${dirtrav}` ;\
	 if   [ "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	 then opttrav="" ;\
	      modtrav=`basename $${dirtrav}` ;\
	 else dirprec1=`dirname $${dirprec}` ;\
	      if   [ "$${dirprec1}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	      then opttrav=`basename $${dirtrav}` ;\
	           modtrav=`basename $${dirprec}` ;\
	      else opttrav="$(MYOPT)" ;\
	           modtrav="" ;\
	      fi ;\
	 fi ;\
	 mv "lib$(NOM_LIB_CS_PRE)$${modtrav}$${opttrav}.a" "$(ARCHLIB_CS)" 
	@echo "***  Library installation finished" ; echo


privincinst:
	@echo "***  Installing headers"
	@if [ ! -d "$(DIR_CS)/arch" ] ; \
	 then mkdir "$(DIR_CS)/arch" ; \
	 fi ;\
         if [ ! -d "$(ARCH_CS)" ] ; \
	 then mkdir "$(ARCH_CS)" ; \
	 fi ;\
         if [ ! -d "$(ARCHINC_CS)" ] ; \
	 then mkdir "$(ARCHINC_CS)" ; \
	 fi
	@dirtrav=`pwd` ;\
	 dirprec=`dirname $${dirtrav}` ;\
	 if   [ "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	 then modtrav=`basename $${dirtrav}` ;\
	 else dirprec1=`dirname $${dirprec}` ;\
	      if   [ "$${dirprec1}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	      then modtrav=`basename $${dirprec}` ;\
	      else modtrav="" ;\
	      fi ;\
	 fi ;\
	 if   [ "$${modtrav}" = "BASE" ] ;\
	 then inc_cs_all="$(INC_CS)/base $(INC_CS)/pprt" ;\
         elif [ "$${modtrav}" = "ATMO" ] ;\
	 then inc_cs_all="$(INC_CS)/atmo" ;\
         elif [ "$${modtrav}" = "CFBL" ] ;\
	 then inc_cs_all="$(INC_CS)/cfbl" ;\
         elif [ "$${modtrav}" = "COGZ" ] ;\
	 then inc_cs_all="$(INC_CS)/cogz" ;\
         elif [ "$${modtrav}" = "CPLV" ] ;\
	 then inc_cs_all="$(INC_CS)/cplv" ;\
         elif [ "$${modtrav}" = "CTWR" ] ;\
	 then inc_cs_all="$(INC_CS)/ctwr" ;\
         elif [ "$${modtrav}" = "ELEC" ] ;\
	 then inc_cs_all="$(INC_CS)/elec" ;\
         elif [ "$${modtrav}" = "FUEL" ] ;\
	 then inc_cs_all="$(INC_CS)/fuel" ;\
         elif [ "$${modtrav}" = "LAGR" ] ;\
	 then inc_cs_all="$(INC_CS)/lagr" ;\
         elif [ "$${modtrav}" = "RAYT" ] ;\
	 then inc_cs_all="$(INC_CS)/rayt" ;\
         elif [ "$${modtrav}" = "MATI" ] ;\
	 then inc_cs_all="$(INC_CS)/mati" ;\
	 fi ;\
	 for dir in $${inc_cs_all} ;\
	 do \
	     for h in $${dir}/*.h ; \
	     do \
		     cp -f "$${h}" "$(ARCHINC_CS)" ; \
		 done ; \
	 done
	@echo "***  Headers installation finished" ; echo


privusrinst:
	@echo "***  Installing user files"
	@if [ ! -d "$(DIR_CS)/arch" ] ; \
	 then mkdir "$(DIR_CS)/arch" ; \
	 fi ;\
     if [ ! -d "$(ARCH_CS)" ] ; \
	 then mkdir "$(ARCH_CS)" ; \
	 fi ;\
	 if [ ! -d "$(ARCH_CS)/share" ] ; \
	 then mkdir "$(ARCH_CS)/share" ; \
	 fi ;\
     if [ ! -d "$(ARCHUSR_CS)" ] ; \
	 then mkdir "$(ARCHUSR_CS)" ; \
	 fi
	@dirtrav=`pwd` ;\
	 dirprec=`dirname $${dirtrav}` ;\
	 if   [ "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	 then modtrav=`basename $${dirtrav}` ;\
	 else dirprec1=`dirname $${dirprec}` ;\
	      if   [ "$${dirprec1}" = "$(BLD_CS)/$(NOM_ARCH)" ] ;\
	      then modtrav=`basename $${dirprec}` ;\
	      else modtrav="" ;\
	      fi ;\
	 fi ;\
	 if   [ "$${modtrav}" = "BASE" ] ;\
	 then user_cs_all="$(USR_CS)/base $(USR_CS)/pprt" ;\
         elif [ "$${modtrav}" = "ATMO" ] ;\
	 then user_cs_all="$(USR_CS)/atmo" ;\
         elif [ "$${modtrav}" = "CFBL" ] ;\
	 then user_cs_all="$(USR_CS)/cfbl" ;\
         elif [ "$${modtrav}" = "COGZ" ] ;\
	 then user_cs_all="$(USR_CS)/cogz" ;\
         elif [ "$${modtrav}" = "CPLV" ] ;\
	 then user_cs_all="$(USR_CS)/cplv" ;\
         elif [ "$${modtrav}" = "CTWR" ] ;\
	 then user_cs_all="$(USR_CS)/ctwr" ;\
         elif [ "$${modtrav}" = "ELEC" ] ;\
	 then user_cs_all="$(USR_CS)/elec" ;\
         elif [ "$${modtrav}" = "FUEL" ] ;\
	 then user_cs_all="$(USR_CS)/fuel" ;\
         elif [ "$${modtrav}" = "LAGR" ] ;\
	 then user_cs_all="$(USR_CS)/lagr" ;\
         elif [ "$${modtrav}" = "RAYT" ] ;\
	 then user_cs_all="$(USR_CS)/rayt" ;\
	 fi ;\
	 for dir in $${user_cs_all} ;\
	 do \
	     moddir=`basename $${dir}` ; \
	     if [ ! -d "$(ARCHUSR_CS)/$${moddir}" ] ; \
		 then mkdir "$(ARCHUSR_CS)/$${moddir}" ; \
		 fi ; \
	     for f in $${dir}/*.f90 ; \
	     do \
		     cp -f "$${f}" "$(ARCHUSR_CS)/$${moddir}/" ; \
			 finst=$(ARCHUSR_CS)/$${moddir}/`basename $${f}` ; \
             oldvers="^\!VERS\n" ; \
             newvers="\!VERS                  Code_Saturne version $(VERSION)\n\!                      --------------------------" ; \
             sed "N;s/$${oldvers}/$${newvers}/" "$${finst}" > "$${finst}.tmp" && mv "$${finst}.tmp" "$${finst}" ; \
		 done ; \
	 done
	@echo "***  User files installation finished" ; echo


privmodlink:
	@liste_obj="" ;\
	 for file in *.o ;\
	 do	if [ "$${file}" != "*.o" ] ;\
		then liste_obj="$${liste_obj} $${file}" ;\
		fi ;\
	 done ;\
	 librariesloc=""; \
	 if   [ "$(MYOPT)" = "EF" ] ;\
	 then \
	   for lib in BASE ATMO CFBL COGZ CPLV CTWR ELEC FUEL LAGR MATI RAYT ;\
	   do \
	       librariesloc="$${librariesloc} -l$(NOM_LIB_CS_PRE)$${lib}DBG" ;\
	   done ;\
	    librariesloc="$${librariesloc} $(LIBRARYEF)" ;\
	 else \
	   for lib in BASE ATMO CFBL COGZ CPLV CTWR ELEC FUEL LAGR MATI RAYT ;\
	   do \
	       librariesloc="$${librariesloc} -l$(NOM_LIB_CS_PRE)$${lib}$(MYOPT)" ;\
	   done ;\
	 fi ;\
	 if [ "$(MYOPT)" = "VERS" ] ;\
	 then librariesloc="" ;\
	 fi ;\
	 fvm_config="$(FVM_HOME)/bin/fvm-config" ;\
	 librariessup="$(FVM_LDFLAGS)" ;\
	 if [ "${MPI}" -eq "1" ] ;\
	   then librariessup="$${librariessup} -lfvm_coupl" ;\
	 fi ;\
	 librariessup="$${librariessup} $(BFT_LDFLAGS)" ;\
	 librariessup="$${librariessup} `$${fvm_config} --ldflags cgns`" ;\
	 librariessup="$${librariessup} `$${fvm_config} --libs cgns`" ;\
	 librariessup="$${librariessup} `$${fvm_config} --ldflags med`" ;\
	 librariessup="$${librariessup} `$${fvm_config} --libs med`" ;\
	 librariessup="$${librariessup} `$${fvm_config} --ldflags hdf5`" ;\
	 librariessup="$${librariessup} `$${fvm_config} --libs hdf5`" ;\
	 if [ "$(LDEDLRPATH)" != "" ] ;\
	 then libpath="" ;\
	     for str in $${librariessup} $(LIBRARY);\
	     do  libdir="" ;\
	         if [ "`echo $${str} | cut -c1-2`" = "-L" ] ;\
	         then libdir=`echo $${str} | cut -c3-` ;\
	         elif [ "`echo $${str} | cut -c1-6`" = "-Wl,-L" ] ;\
	         then libdir=`echo $${str} | cut -c7-` ;\
	         fi ;\
	         if [ "$${libdir}" != "" ] ;\
               then islib=`find $${libdir} -name "*.so" |wc -w` ;\
		     if [ $$islib -gt 0 ] ;\
		     then \
		         if [ "$${libpath}" = "" ] ;\
	    	         then libpath="$(LDEDLRPATH)$${libdir}" ;\
	    	         else libpath="$${libpath}:$${libdir}" ;\
		         fi ;\
		     fi ;\
	         fi ;\
	     done ;\
	     librariessup="$${librariessup} $${libpath}" ;\
         fi ;\
         libraries="-L$(ARCHLIB_CS) $${librariesloc} $${librariesloc} $${librariessup} $(LIBRARY)" ;\
	     link_cmd="$(LDEDL) -o $(NOM_EXE) $(ALLLDEDLFLAGS)" ;\
	     link_cmd="$${link_cmd} $${liste_obj} $${libraries} $(MYLIBS)" ;\
	     echo "$${link_cmd}" ;\
	     $${link_cmd}
	@echo ; echo "***  Linking finished; library:" ; echo
	@ls -l $(NOM_EXE)


privlink:
	@echo "***  Linking"
	@library=" $(LIBRARY) " ;\
	 if   [ "$(MYOPT)" = "LO" ] ;\
	 then library="$${library} $(LIBRARYLO) " ;\
	 elif [ "$(MYOPT)" = "DBG" -o "$(MYOPT)" = "EF" ]  ;\
	 then library="$${library} $(LIBRARYDBG) " ;\
	 elif [ "$(MYOPT)" = "PROF" ]  ;\
	 then library="$${library} $(LIBRARYOPT) " ;\
	 elif [ "$(MYOPT)" = "VERS" ]  ;\
	 then library="$${library}                 " ;\
	 elif [ "$(MYOPT)" = "" ]  ;\
	 then library="$${library} $(LIBRARYOPT) " ;\
	 fi ; \
	 if [ "${MPI}" -eq "1" ] ;\
	 then \
	   if [ "${MPE}" -eq "1" ] ;\
	     then library="${MPE_LIB} $${library} " ;\
	   fi ; \
           library="${MPI_LIB} $${library} " ;\
	 fi ; \
	 if [ "${SOCKET}" -eq "1" ] ;\
	 then library="$${library} ${SOCKET_LIB} " ;\
	 fi ; \
	 if [ "${BLAS}" -eq "1" ] ;\
	 then library="${BLAS_LDFLAGS} $${library} " ;\
	 fi ; \
	 if [ "${XML}" -eq "1" ] ;\
	 then library="${XML_LIB} $${library} " ;\
	 fi ; \
	 if   [ "$(MYOPT)" = "LO" ] ;\
	 then make privmodlink \
	           ALLLDEDLFLAGS="$(MYLDEDLFLAGS) $(LDEDLFLAGSLO)" \
	                 LIBRARY="$${library}" ;\
	 elif [ "$(MYOPT)" = "DBG" ]  ;\
	 then make privmodlink \
	           ALLLDEDLFLAGS="$(MYLDEDLFLAGS) $(LDEDLFLAGSDBG)" \
	                 LIBRARY="$${library}" ;\
	 elif [ "$(MYOPT)" = "EF" ]  ;\
	 then make privmodlink \
	           ALLLDEDLFLAGS="$(MYLDEDLFLAGS) $(LDEDLFLAGSDBG)" \
	                 LIBRARY="$${library}" ;\
	 elif [ "$(MYOPT)" = "PROF" ] ;\
	 then make privmodlink \
	           ALLLDEDLFLAGS="$(MYLDEDLFLAGS) $(LDEDLFLAGS) $(LDEDLFLAGSPROF)" \
	                 LIBRARY="$${library}" ;\
	 elif [ "$(MYOPT)" = "VERS" ] ;\
	 then make privmodlink \
	           ALLLDEDLFLAGS="$(MYLDEDLFLAGS) $(LDEDLFLAGSVERS)" \
	                 LIBRARY="$${library}" ;\
	 elif [ "$(MYOPT)" = "" ] ;\
	 then make privmodlink \
	           ALLLDEDLFLAGS="$(MYLDEDLFLAGS)  $(LDEDLFLAGS)" \
	                 LIBRARY="$${library}" ;\
	 fi


all $(NOM_EXE): privtitre privcomp privlink privfin

#
# privconfig used cd commands so that grep does not suffer from
#   diverse possibilities  ~foo, $HOME, ../..
#   as such, $(BLD_CS) must be accessible ... 
#   if it is not and we are a user, it means we have made an error
# awk is used to check if the dirbuild string is part of the dirtrav string

privconfig:
	@dirtrav=`pwd` ;\
	 if   [   -d "$(BLD_CS)" -a   -x "$(BLD_CS)" ] ;\
	 then dirbuild=$(BLD_CS) ;\
	      findirtrav=`echo "$${dirtrav}" | awk -F"$${dirbuild}" '{print $$2}'` ;\
	      if   [ "$${findirtrav}" != "" ] ;\
	      then echo "***  \"make\" in a library generation directory" ;\
	           make privconflib  ;\
	      else echo "***  \"make\" in a user directory" ;\
	      fi ;\
	 elif [   -d "$(BLD_CS)" -a ! -x "$(BLD_CS)" ] ;\
	 then echo "***  \"make\" in a user directory" ;\
	 elif [ ! -d "$(BLD_CS)"                   ] ;\
	 then echo "***  \"make\" in a user directory" ;\
	      echo "              ($(BO)$(BLD_CS)$(NO) does not exist)";\
	      tput init ;\
	 fi

#
# privconflib checks that the OPT option is consistent with
#           the working directory
#

privconflib:
	@dirtrav=`pwd` ;\
	 module=`basename $${dirtrav}` ;\
	 if   [ "$${dirtrav}" = "$(BLD_CS)/$(NOM_ARCH)/$${module}" ] ;\
	 then opttrav="" ;\
	      modtrav=`basename $${dirtrav}` ;\
	      dirprec=`dirname $${dirtrav}` ;\
	 else opttrav=`basename $${dirtrav}` ;\
	      dirprec=`dirname $${dirtrav}` ;\
	      module=`basename $${dirprec}` ;\
	      if   [ "$${dirprec}" = "$(BLD_CS)/$(NOM_ARCH)/$${module}" ] ;\
	      then modtrav=`basename $${dirprec}` ;\
	           dirprec=`dirname $${dirprec}` ;\
	      else echo ;\
	           echo "$(SO)CONFIGURATION ERROR$(NO)" ;\
	           echo ;\
	           echo "$(BO)You must be in directory:$(NO)" ;\
	           echo "$(BO)$(BLD_CS)/$(NOM_ARCH)/$${module}$(NO)" ;\
	           echo "$(BO)or in one of its sub-directories$(NO)" ;\
	           echo ;\
	           tput init ;\
	           exit 1 ;\
	      fi ;\
	 fi ;\
	 dirprec=`dirname $${dirprec}` ;\
	 dirlib=`dirname $${dirprec}` ;\
	 if   [ ! "$$dirlib" = "$(DIR_CS)" ] ;\
	 then echo ;\
	      echo "$(SO)CONFIGURATION ERROR$(NO)" ;\
	      echo ;\
	      echo "$(BO) You must be in the general directory: $(DIR_CS)$(NO)" ;\
	      echo ;\
	      tput init ;\
	      exit 1 ;\
	 fi ;\
	 if   [ ! "$${modtrav}" = "$${module}" ] ;\
	 then echo ;\
	      echo "$(SO)CONFIGURATION ERROR$(NO)" ;\
	      echo ;\
	      echo "$(BO)You must be in directory: $${module}$(NO)" ;\
	      echo ;\
	      tput init ;\
	      exit 1 ;\
	 fi ;\
	 if   [ ! "$${opttrav}" = "$(MYOPT)" ] ;\
	 then echo ;\
	      echo "$(SO)CONFIGURATION ERROR$(NO)" ;\
	      echo ;\
	      if   [ "$(MYOPT)" = "" ] ;\
	      then echo "$(BO)In directory $${opttrav}$(NO)" ;\
	           echo "$(BO) you must use make option: OPT=$${opttrav}$(NO)" ;\
	      else echo "$(BO)You must be in directory $(MYOPT)$(NO)" ;\
	           echo "$(BO) to use option OPT=$(MYOPT)$(NO)" ;\
	      fi ;\
	      echo ;\
	      tput init ;\
	      exit 1 ;\
	 fi ;\
	 if [ "$${modtrav}" = "BASE" ] ;\
	 then INC_CS="toto" ;\
	 fi 


privfin:
	@echo


clean purge rmobj: privconfig
	@\rm -f *.o *.inl core \
	      .liste_mkfiles.* .liste_mkoptpart.o so_locations
	@for file in *.c *.f90 *.F90 *.h ;\
	 do \
	      if [ -h "$${file}" ] ;\
	      then \rm -f $${file} ;\
	      fi ;\
	 done


rmexe: privconfig
	@\rm -f $(NOM_EXE) *.out

rmlib: privconfig
	@\rm -f "$(NOM_LIB_CS_PRE)$(MYOPT).a"

new: privconfig clean rmexe privfin


help info: privtitre
	@echo
	@echo "Usage : $(BO)make$(NO)                     (in $(BO)build/$(NOM_ARCH)/$(NO)$(SO)<MOD>$(NO)$(BO)/$(NO))"
	@echo "   or : $(BO)make$(NO)          $(BO)OPT=$(NO)$(SO)<OPT>$(NO)  (in $(BO)build/$(NOM_ARCH)/$(NO)$(SO)<MOD>$(NO)$(BO)/$(NO)$(SO)<OPT>$(NO)$(BO)/$(NO))"
	@echo
	@echo "   or : $(BO)make$(NO) $(SO)<target>$(NO) $(BO)OPT=$(NO)$(SO)<OPT>$(NO)  (in a user directory)"
	@echo "                                 for Code_Saturne link"
	@echo
	@echo "     with:"
	@echo "     $(SO)<target>$(NO) : $(BO)saturne$(NO)     : Code_Saturne compilation and link"
	@echo "                $(BO)all$(NO)         : same"
	@echo "                $(BO)lib$(NO)         : library creation"
	@echo "                $(BO)lib_install$(NO) : library creation and installation"
	@echo "                $(BO)comp$(NO)        : compilation only"
	@echo "                $(BO)link$(NO)        : link of Code_Saturne"
	@echo "                $(BO)saturne$(NO)     : same"
	@echo "                $(BO)clean$(NO)       : delete object files"
	@echo "                $(BO)purge$(NO)       : same"
	@echo "                $(BO)rmlib$(NO)       : delete library"
	@echo "                $(BO)rmobj$(NO)       : same"
	@echo "                $(BO)rmexe$(NO)       : delete executable"
	@echo "                $(BO)new$(NO)         : clean + rmlib + rmexe"
	@echo "                $(BO)help$(NO)        : this message"
	@echo "                $(BO)info$(NO)        : same"
	@echo
	@echo "     $(SO)<OPT>$(NO)    : $(BO)LO$(NO)       : option for reduced optimization"
	@echo "                $(BO)DBG$(NO)      : option for debug"
	@echo "                $(BO)EF$(NO)       : option for debug with ElectricFence"
	@echo "                $(BO)VERS$(NO)     : option for compiler versions"
	@echo "                $(BO)PROF$(NO)     : option for profiling"
	@echo
	@echo "     and the additional user options:"
	@echo "     $(SO)<OPTC>$(NO)    : C compiler options"
	@echo "     $(SO)<OPTF>$(NO)    : Fortran compiler options"
	@echo "     $(SO)<OPTLD>$(NO)   : linker options"
	@echo "     $(SO)<OPTLIBS>$(NO) : additional library options"
	@echo
	@echo "   Environment configuration queries (for run script):"
	@echo "     $(BO)make echo_mpi_bin$(NO) : echo MPI environment executables path"
	@echo

