dnl--------------------------------------------------------------------------------
dnl
dnl This file is part of Code_Saturne, a general-purpose CFD tool.
dnl
dnl Copyright (C) 1998-2019 EDF S.A.
dnl
dnl This program is free software; you can redistribute it and/or modify it under
dnl the terms of the GNU General Public License as published by the Free Software
dnl Foundation; either version 2 of the License, or (at your option) any later
dnl version.
dnl
dnl This program is distributed in the hope that it will be useful, but WITHOUT
dnl ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
dnl FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
dnl details.
dnl
dnl You should have received a copy of the GNU General Public License along with
dnl this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
dnl Street, Fifth Floor, Boston, MA 02110-1301, USA.
dnl
dnl--------------------------------------------------------------------------------

# CS_AC_TEST_DOCS
#----------------
# Macro function grouping the different tests for documentation generation

AC_DEFUN([CS_AC_TEST_DOCS], [

CS_AC_TEST_LATEX
CS_AC_TEST_DOXYGEN
CS_AC_TEST_SPHINX

])dnl

# CS_AC_TEST_LATEX
#-----------------
# modifies or sets cs_have_latex

AC_DEFUN([CS_AC_TEST_LATEX],[

cs_have_latex=yes

AC_ARG_VAR([PDFLATEX], [LaTeX PDF generator tool])

dnl where is pdflatex ?
AC_PATH_PROG(PDFLATEX, [pdflatex])
if test "x$PDFLATEX" = "x"; then
  AC_MSG_WARN(pdflatex not found)
  cs_have_latex=no
fi

AC_ARG_VAR([BIBTEX], [LaTeX references tool])

dnl where is bibtex ?
AC_PATH_PROG(BIBTEX, [bibtex])
if test "x$BIBTEX" = "x"; then
  AC_MSG_WARN(bibtex not found)
  cs_have_latex=no
fi

AC_ARG_VAR([MAKEINDEX], [LaTeX indexes tool])

dnl where is makeindex ?
AC_PATH_PROG(MAKEINDEX, [makeindex])
if test "x$MAKEINDEX" = "x"; then
  AC_MSG_WARN(makeindex not found)
  cs_have_latex=no
fi

AC_ARG_VAR([FIG2DEV], [Xfig translation tool])

dnl where is fig2dev ?
AC_PATH_PROG(FIG2DEV, [fig2dev])
if test "x$FIG2DEV" = "x"; then
  AC_MSG_WARN(fig2dev not found)
  cs_have_latex=no
fi

dnl So as to correctly set TEXINPUTS environment variable, one needs to use
dnl the system dependant path separator
if test "$host_os" = mingw32 ; then
  cs_tex_path_end=
  cs_tex_path_sep=';'
else
  cs_tex_path_end=':'
  cs_tex_path_sep=':'
fi
AC_SUBST([cs_tex_path_end])
AC_SUBST([cs_tex_path_sep])

AM_CONDITIONAL(HAVE_LATEX, [test $cs_have_latex = yes])
AC_SUBST(cs_have_latex)

])dnl

# CS_AC_TEST_DOXYGEN
#-------------------
# modifies or sets cs_have_doxygen

AC_DEFUN([CS_AC_TEST_DOXYGEN],[

cs_have_doxygen=yes

AC_ARG_VAR([DOXYGEN], [source code documentation generator])

dnl where is doxygen ?
AC_PATH_PROG(DOXYGEN, [doxygen])
if test "x$DOXYGEN" = "x"; then
  AC_MSG_WARN(doxygen not found)
  cs_have_doxygen=no
fi

AC_ARG_VAR([DOT], [graphs generator])

dnl where is dot ?
AC_PATH_PROG(DOT, [dot])
if test "x$DOT" = "x"; then
  AC_MSG_WARN(dot not found)
  cs_have_doxygen=no
fi

AM_CONDITIONAL(HAVE_DOXYGEN, [test $cs_have_doxygen = yes])
AC_SUBST(cs_have_doxygen)

])dnl

# CS_AC_TEST_SPHINX
#------------------
# modifies or sets cs_have_sphinx, SPHINX depending on libraries found

AC_DEFUN([CS_AC_TEST_SPHINX],[

cs_have_sphinx=yes

AC_ARG_VAR([SPHINXBUILD], [Sphinx documentation tool])

dnl where is sphinx ?
AC_PATH_PROG(SPHINXBUILD, sphinx-build)
if test "x$SPHINXBUILD" = "x"; then
  AC_MSG_WARN(sphinx-build not found)
  cs_have_sphinx=no
fi

AM_CONDITIONAL(HAVE_SPHINX, [test $cs_have_sphinx = yes])
AC_SUBST(cs_have_sphinx)

])dnl
