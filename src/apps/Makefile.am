## Process this file with automake to produce Makefile.in

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

# Files not built through classical Automake rules need to be explicitely
# added to distribution.

EXTRA_DIST = cs_solver.c

# Build flat archives

# PLE dependency

if HAVE_INTERNAL_PLE
  CPPFLAGS_PLE = -I"$(top_srcdir)/libple/src" -I"$(top_builddir)/libple"
else
  CPPFLAGS_PLE = $(PLE_CPPFLAGS)
endif

# Relocatable installation

if HAVE_RELOCATABLE
  ORIGINRUNPATH = ${LDRPATH}'\$$ORIGIN/../../lib'
  PLUGINRUNPATH_R = ${LDRPATH}'\$$ORIGIN/../../lib/'$(PACKAGE)
else
  PLUGINRUNPATH_N = ${LDRPATH}$(pkglibdir)
endif

if HAVE_PLUGINS
  PLUGINRUNPATH = $(PLUGINRUNPATH_R)$(PLUGINRUNPATH_N)
endif

# Main part

if HAVE_BACKEND

AM_CPPFLAGS = \
-I$(top_srcdir)/src/atmo \
-I$(top_srcdir)/src/alge \
-I$(top_srcdir)/src/base \
-I$(top_srcdir)/src/bft \
-I$(top_srcdir)/src/cdo \
-I$(top_srcdir)/src/cogz \
-I$(top_srcdir)/src/comb \
-I$(top_srcdir)/src/ctwr \
-I$(top_srcdir)/src/elec \
-I$(top_srcdir)/src/fvm \
-I$(top_srcdir)/src/gui \
-I$(top_srcdir)/src/gwf \
-I$(top_srcdir)/src/lagr \
-I$(top_srcdir)/src/meg \
-I$(top_srcdir)/src/mesh \
-I$(top_srcdir)/src/pprt \
-I$(top_srcdir)/src/turb \
$(CPPFLAGS_PLE) \
$(MPI_CPPFLAGS) $(CPPFLAGS_BLAS) \
$(CPPFLAGS_MUMPS)

AM_CFLAGS = $(CFLAGS_DBG) $(CFLAGS_OPT) $(CFLAGS_SHARED)

# Global LDFLAGS must be used so as to handle the possibility of a CXX link.
# Per-executable flags will override this global definition.
AM_LDFLAGS = $(ORIGINRUNPATH) $(PLUGINRUNPATH)

# We use a local install rule for libraries here

lib_LIBRARIES =
noinst_LIBRARIES =

local_libraries_add = libsaturne.${LIB_EXT} libcs_solver.${LIB_EXT}
local_libraries_v_add =

libsaturne_a_libs = \
$(top_builddir)/src/fvm/libfvm_filters.a \
$(top_builddir)/src/base/libcsbase.a \
$(top_builddir)/src/base/libcscore.a \
$(top_builddir)/src/base/libcsrenumber.a \
$(top_builddir)/src/fvm/libfvm.a \
$(top_builddir)/src/pprt/libcspprt.a \
$(top_builddir)/src/cdo/libcscdo.a \
$(top_builddir)/src/alge/libcsalge.a \
$(top_builddir)/src/mesh/libcsmesh.a \
$(top_builddir)/src/mesh/libcspartition.a \
$(top_builddir)/src/turb/libcsturb.a \
$(top_builddir)/src/atmo/libcsatmo.a \
$(top_builddir)/src/cfbl/libcscfbl.a \
$(top_builddir)/src/cogz/libcscogz.a \
$(top_builddir)/src/comb/libcscomb.a \
$(top_builddir)/src/ctwr/libcsctwr.a \
$(top_builddir)/src/elec/libcselec.a \
$(top_builddir)/src/gui/libcsgui.a \
$(top_builddir)/src/gwf/libcsgwf.a \
$(top_builddir)/src/lagr/libcslagr.a \
$(top_builddir)/src/rayt/libcsrayt.a \
$(top_builddir)/src/user/libcsuser.a \
$(top_builddir)/src/meg/libcsmeg.a \
$(top_builddir)/src/libcsmodules.a \
$(top_builddir)/src/bft/libbft.a

cs_solver = $(top_srcdir)/src/apps/cs_solver.c
cs_solver.$(OBJEXT): $(cs_solver)
	$(CC) -c $(DEFS) $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $(AM_CFLAGS) $(CFLAGS) ${cs_solver}

if ENABLE_SHARED

libsaturne.${SHARED_LIB_EXT}: ${libsaturne_a_libs}
	rm -f $@ ; \
	ldadd_saturne=`PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	CS_TOP_BUILDDIR=$(top_builddir) \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	 --lib-flags-only` ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py -o $@ --echo --linker ${CS_LD} --version ${cs_version_lib} --whole-archive-start ${libsaturne_a_libs} --whole-archive-end ${LDFLAGS_SHARED} ${CXXFLAGS_SHARED} ${LDFLAGS} ${AM_LDFLAGS} $${ldadd_saturne}

libcs_solver.${SHARED_LIB_EXT}: libsaturne.${SHARED_LIB_EXT} cs_solver.$(OBJEXT)
	rm -f $@ ; \
	ldadd_saturne=`PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	CS_TOP_BUILDDIR=$(top_builddir) \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	 --lib-flags-only` ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py -o $@ --echo --linker ${CS_LD} --version ${cs_version_lib} ${LDFLAGS_SHARED} ${CXXFLAGS_SHARED} cs_solver.$(OBJEXT) ${LDFLAGS} -L. -lsaturne-${cs_version_lib} ${AM_LDFLAGS} $${ldadd_saturne}

local_libraries_v_add += \
libsaturne-${cs_version_lib}.${SHARED_LIB_EXT} \
libcs_solver-${cs_version_lib}.${SHARED_LIB_EXT}

else

libsaturne.a: ${libsaturne_a_libs}
	rm -f $@ ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py -o $@ -v --archive ${libsaturne_a_libs}

libcs_solver.a: cs_solver.$(OBJEXT)
	rm -f $@ ; \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py -o $@ -v --archive cs_solver.$(OBJEXT)

endif

endif

# Main program
#pkglibexecdir = $(libexecdir)/$(PACKAGE)
pkglibexec_PROGRAMS =

# code_saturne IO utility program

if HAVE_FRONTEND
pkglibexec_PROGRAMS += cs_io_dump
endif

if HAVE_BACKEND
cs_solver$(EXEEXT): libcs_solver.${LIB_EXT}
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	CS_TOP_BUILDDIR=$(top_builddir) \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	libcs_solver.${LIB_EXT}
endif

# code_saturne IO utility (minimal dependencies)

if HAVE_FRONTEND
cs_io_dump_CPPFLAGS = -DLOCALEDIR=\"'$(localedir)'\" -I$(top_srcdir)/src/base
cs_io_dump_SOURCES = cs_io_dump.c
cs_io_dump_LDADD =
endif

# cs_solver executable is built using the same Python script as when
# using user functions (rather than with libtool), so it must be
# handled separately

all-local: ${local_libraries_add} cs_solver$(EXEEXT)

# Remark:
# When installing shared libraries, remove library then copy it, rather than
# simply overwriting it, so that if a running instance is using that library,
# it will still access the old copy (not really removed until closed based
# on POSIX file semantics), instead of a mofified copy that can lead to
# SIGBUS errors due to an incomplete (in progress) copy or possibly other
# address changes.

if ENABLE_SHARED

install-exec-hook: ${local_libraries_add}
	@list=; for p in ${local_libraries_v_add} ${local_libraries_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	test -z "$$list" -a -n "$(libdir)" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(libdir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(libdir)" || exit 1; \
	  for p in $$list; do \
	    echo " rm -f '$(DESTDIR)$(libdir)/$$p'"; \
	    rm -f "$(DESTDIR)$(libdir)/$$p"; \
	    if test -h $$p; then \
	      echo " cp -fp $$p '$(DESTDIR)$(pkglibdir)'"; \
	      cp -fpP $$p "$(DESTDIR)$(libdir)"; \
	    fi; \
	  done; \
	  ldadd_saturne=`PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	  CS_TOP_BUILDDIR=$(top_builddir) \
	  $(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	  --lib-flags-only --mode=install --dest "$(DESTDIR)"` ; \
	  $(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py -o "$(DESTDIR)$(libdir)/libsaturne.${SHARED_LIB_EXT}" --echo --linker ${CS_LD} --version ${cs_version_lib} --whole-archive-start ${libsaturne_a_libs} --whole-archive-end ${LDFLAGS_SHARED} ${CXXFLAGS_SHARED} ${LDFLAGS} ${AM_LDFLAGS} $${ldadd_saturne}; \
	  $(PYTHON) -B $(top_srcdir)/build-aux/cs_link_library.py -o "$(DESTDIR)$(libdir)/libcs_solver.${SHARED_LIB_EXT}" --echo --linker ${CS_LD} --version ${cs_version_lib} ${LDFLAGS_SHARED} ${CXXFLAGS_SHARED} cs_solver.$(OBJEXT) ${LDFLAGS} -L"$(DESTDIR)$(libdir)" -lsaturne-${cs_version_lib} ${AM_LDFLAGS} $${ldadd_saturne}; \
	}; \
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	CS_TOP_BUILDDIR=$(top_builddir) \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	 --mode=install --dest "$(DESTDIR)"

else

install-exec-hook: ${local_libraries_add}
	@list=; for p in ${local_libraries_v_add} ${local_libraries_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	test -z "$$list" -a -n "$(libdir)" || { \
	  echo " $(MKDIR_P) '$(DESTDIR)$(libdir)'"; \
	  $(MKDIR_P) "$(DESTDIR)$(libdir)" || exit 1; \
	  for p in $$list; do \
	    echo " rm -f '$(DESTDIR)$(libdir)/$$p'"; \
	    rm -f "$(DESTDIR)$(libdir)/$$p"; \
	    echo " cp -fp $$p '$(DESTDIR)$(pkglibdir)'"; \
	    cp -fpP $$p "$(DESTDIR)$(libdir)"; \
	  done; \
	}; \
	PYTHONPATH=$(top_srcdir)/python/code_saturne/base \
	CS_TOP_BUILDDIR=$(top_builddir) \
	$(PYTHON) -B $(top_srcdir)/build-aux/cs_compile_build.py \
	 --mode=install --dest "$(DESTDIR)"

endif

clean-local:
	-rm -f cs_solver$(EXEEXT)
	@list=; for p in ${local_libraries_add} ${local_libraries_v_add}; do \
	  if test -f $$p; then \
	    echo "rm -f $$p"; \
	    rm -f $$p; \
	  else :; fi; \
	done

uninstall-local:
	rm $(pkglibexecdir)/cs_solver$(EXEEXT)
	@list=; for p in ${local_libraries_add} ${local_libraries_v_add}; do \
	  if test -f $$p; then \
	    list="$$list $$p"; \
	  else :; fi; \
	done; \
	for p in $$list; do \
	  echo " rm -f '$(DESTDIR)$(libdir)/$$p'"; \
	  rm -f "$(DESTDIR)$(libdir)/$$p"; \
	done; \
	rmdir "$(DESTDIR)$(libdir)" || exit 0
