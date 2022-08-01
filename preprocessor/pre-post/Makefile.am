## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

# Main part

AM_CPPFLAGS = -I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/util \
$(HDF5_CPPFLAGS) $(HDF5_CPPFLAGS_MPI) \
$(MED_CPPFLAGS) $(CGNS_CPPFLAGS)
AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT)

# Headers

noinst_HEADERS = \
ecs_table_comm.h \
ecs_table_post_cgns.h \
ecs_table_post_ens.h \
ecs_table_post.h \
ecs_table_post_med.h \
ecs_maillage_ncs.h \
ecs_maillage_post.h \
ecs_maillage_pre.h \
ecs_med.h \
ecs_med_priv.h \
ecs_post_cgns.h \
ecs_post_cgns_priv.h \
ecs_post_ens.h \
ecs_post_ens_priv.h \
ecs_post.h \
ecs_post_med.h \
ecs_pre_ccm.h \
ecs_pre_cgns.h \
ecs_pre_ens.h \
ecs_pre_gambit.h \
ecs_pre_gmsh.h \
ecs_pre.h \
ecs_pre_ideas.h \
ecs_pre_med.h \
ecs_pre_nopo.h

# Library source files

noinst_LIBRARIES = libprepostecs.a
libprepostecs_a_LIBADD =

libprepostecs_a_CPPFLAGS = \
-DPKGLIBDIR=\"'$(pkglibdir)'\" \
$(AM_CPPFLAGS)

libprepostecs_a_SOURCES = \
ecs_table_comm.c \
ecs_table_post.c \
ecs_table_post_cgns.c \
ecs_table_post_ens.c \
ecs_table_post_med.c \
ecs_maillage_ncs.c \
ecs_maillage_post.c \
ecs_maillage_pre.c \
ecs_med.c \
ecs_post.c \
ecs_post_cgns.c \
ecs_post_ens.c \
ecs_post_med.c \
ecs_pre.c \
ecs_pre_cgns.c \
ecs_pre_ens.c \
ecs_pre_gambit.c \
ecs_pre_gmsh.c \
ecs_pre_ideas.c \
ecs_pre_med.c \
ecs_pre_nopo.c

local_libraries_add =

# CCMIO (plugin or linked)

if HAVE_CCM

ecs_pre_ccm_CPPFLAGS = \
-I$(top_srcdir)/preprocessor/base -I$(top_srcdir)/preprocessor/util \
$(CCM_CPPFLAGS) $(CCM_CXXFLAGS) $(MPI_CPPFLAGS)
ecs_pre_ccm_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)

ecs_pre_ccm = $(top_srcdir)/preprocessor/pre-post/ecs_pre_ccm.c
ecs_pre_ccm.o: $(ecs_pre_ccm)
	$(CC) -c $(DEFS) $(DEFAULT_INCLUDES) $(ecs_pre_ccm_CPPFLAGS) $(ecs_pre_ccm_CFLAGS) $(CFLAGS) $(ecs_pre_ccm)

if HAVE_PLUGINS

ecs_ccm.${SHARED_LIB_EXT}: ecs_pre_ccm.o
	@echo "Building ecs_ccm plugin..."
	rm -f $@
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py --echo --linker ${CC} --std-search-paths='no' -o $@ ${LDFLAGS_SHARED} ${CFLAGS_SHARED} ecs_pre_ccm.o ${LDFLAGS} $(CCM_LDFLAGS) $(CCM_LIBS) $(CCMRUNPATH)

local_libraries_add += ecs_ccm.${SHARED_LIB_EXT}

else

libprepostecs_a_LIBADD += ecs_pre_ccm.o

endif

endif

# Local targets

all-local: ${local_libraries_add}

clean-local:
	@list=; for p in ${local_libraries_add}; do \
	  if test -f $$p; then \
	    echo "rm -f $$p"; \
	    rm -f $$p; \
	  else :; fi; \
	done

install-exec-local: ${local_libraries_add}
	@list=; for p in ${local_libraries_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	test -z "$$list" -a -n "$(pkglibdir)" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(pkglibdir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(pkglibdir)" || exit 1; \
	  echo " cp -fp $$list '$(DESTDIR)$(pkglibdir)'"; \
	  cp -fpP $$list "$(DESTDIR)$(pkglibdir)"; \
	}

uninstall-local:
	@list=; for p in ${local_libraries_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	for p in $$list; do \
	  echo " rm -f '$(DESTDIR)$(pkglibdir)/$$p'"; \
	  rm -f "$(DESTDIR)$(pkglibdir)/$$p"; \
	done; \
	rmdir "$(DESTDIR)$(pkglibdir)" || exit 0
