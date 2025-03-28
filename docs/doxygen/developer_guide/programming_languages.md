<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
-->

\page cs_dg_programming_languages Programming languages

[TOC]

Programming languages {#sec_prg_programming_languages}
=====================

The code\_saturne tool is written mostly in C++, but was only recently converted from C,
so most programming constructs encountered in the code will be Very
similar to C rather than "modern" C++.

Subsets of C++ used in code_saturne
-----------------------------------

- Functions with the same name but different arguments.
  This "syntaxic sugar" can allow handling variants of a function
  without requiring long and excessively verbose names.

- Templated functions.
  Using templated types allows a single function for various data types,
  which allows for better code factoring.

- Lambda functions.
  Lambda functions are an essential part of the `cs_dispatch` mechanism
  used to handle local parallel constructs, and enabling generation
  of CPU and GPU execution from the same source code.

Subsets of C++ not (or rarely) used in code_saturne.
----------------------------------------------------

- Class methods.
  Since most of code\_saturne was only recently converted to C++, although
  many structures use an object-oriented approach, they are still only defined as C structures and do not (yet)
  use a C++ method syntax.
  - For example, given a `cs_field_t` structure `*f`, we can call
    `a = cs_field_get_key_double(f, k)`, but not `a = f->get_key_double(k)`.
  - For the same reason, constructors and destructors are explicit.

- Inheritance.
  For the same reasons, code\_saturne does do not currently use class inheritance. Class inheritance may be used in the future, though
  caution must be used with virtual methods, to avoid performance
  issues (i.e. these must be used only for high-level contructs
  and operations).

- Exceptions.
  Actually, a few `try`/`catch` exceptions handling examples may be
  found in the SYCL portions of the code, but exceptions must be avoided
  in performance-ciritical portions of the code (and were not available
  in C).
