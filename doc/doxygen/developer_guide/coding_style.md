<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

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

\page cs_dg_coding_style Coding style guidelines

[TOC]

Master rule {#sec_prg_master_rule}
===========

**Keep the style consistent !**

This rule should be observed above all others. The coding style in code_saturne
has evolved over the years, but unless you are ready to update a whole
file to a more current style (in which case the other guidelines should be
followed), try to remain consistent with the style in the current file.

For new files, use recently updated or reference examples, such as
`src/base/cs_field.cpp` and  `src/base/cs_field.h`, `src/base/cs_array.cpp`,
and `src/base/cs_array.h` for C++.

If you consider your preferred coding style is better than the one used is
superior in some manner to the one used, suggestions are welcome, and
may lead to evolutions in the current style if the arguments are convincing,
but remember that an inconsistent style is more time-consuming and difficult
to read and understand.

Developers trying to follow a consistent style should expect help and
counseling when they have questions, but those whose time is too precious
to observe the existing code before providing their own or requesting assistance
should understand that this is also true of those whom they expect to read
their code.

General rules
=============

The following general rules are strongly recommended:

- Except for files in which they have a special meaning (such as
  Makefiles), use spaces, not tabs. *Absolutely* avoid tabs in
  Python code. Most importantly, use a decent text editor that does not
  randomly mix spaces and tabs.
  The `sbin/rmb` script can be used to remove trailing white-space and
  replace tabs with spaces, but this may appear to damage indentation when
  it is defined with an odd mix of spaces and tabs, so using a good editor
  in the first place is always preferred.
- 80 characters maximum line length; split lines longer than this to ensure
  readability on small screens, or (more commonly) when reading or comparing
  code side-by-side on wider screens. This rule is less important for LaTeX
  documentation, and may be difficult to follow for some Doxygen or Markdown
  constructs (one could argue that using one line per paragraph and relying
  on line wrapping would actually make revision merging simpler).

For new developments, prefer C++ to C, as the code base will be progressively
migrated to C++ (though C++ features will only be introduced by small
increments)..

Punctuation
-----------

Except when adding additional white space to align similar definitions
or arguments on several lines, standard English punctuation rules should be
followed:

- No white space before a punctuation mark (, ; . ), one white space
  after a punctuation mark.

- Use white space before an opening parenthesis, no white space after an opening
  parenthesis.

- No white space before a closing parenthesis, white-space after a closing
  parenthesis.

- Keep blank lines spacing symmetric: if there is a line after an opening
  brace, there should also be one after the closing brace, and vice-versa.

  For example, the following styles are valid:
  ```
  (
    for (cs_lnum_t i = 0; i < stride; i++) {
      w2[c_id][i] = (pvar[c_id][i] - p_mean[i]);
    }
  }
  ```

  and:
  ```
  (
    for (cs_lnum_t i = 0; i < stride; i++) {

      w2[c_id][i] = (pvar[c_id][i] - p_mean[i]);

    }
  }
  ```

  though in general, the denser style is preferred for short
  code blocks, and the one with more spacing for larger blocks.

- There should always be a blank line between code or comments and "full-width"
  separators, for example:
  ```
  /*=========================================================================*/

  /* Section starts here */
  ```

C++ coding style
================

The code_saturne coding style inherits from common conventions, with
a few specific additions.

General rules
-------------

The following presentation rules should be followed:

- Indentation step: 2 characters.

- Do not use tabulation characters, do not leave whitespace
  at the end of lines.

- Use [snake case](https://en.wikipedia.org/wiki/Snake_case) for
  class and variable names, with the following special cases:
  * Enumerations and macros must be in uppercase.
    - Some macros, such as those used for backwards compatibility
      may be in lowercase if this cannot be avoided.
  * A mix of lowercase and uppercase characters (for example CamelCase,
    often encountered in C++ libraries) is allowed in sections
    dealing specifically with external libraries using such coding styles.

- Header `.h` files should have a mechanism to prevent
  multiple inclusions.

The following coding rules are strongly recommended:

- All macro parameters must be enclosed inside parentheses.

- A function's return type must always be defined.

- Variables should be initialized before use
  (pointers are set to NULL). A good compiler should issue warnings when
  this is not the case, and those warnings must be acted upon;

- when a structure definition is only needed in a single file,
  it is preferred to define it directly in the C source file,
  so as to make as little visible as possible in the matching header file.
  structures only used through pointers may be made opaque in this
  manner, which ensures that their possible future modification should
  not require changes in other parts of the code or have unexpected side-effects.

- When a public function is defined in a C source file, a matching
  header file containing its prototype must be included.

- Usage of global variables must be kept to a minimum, though such
  variables may be useful to maintain state or references to mesh or
  variable structures.
  If a global variable is only needed inside a single file, it should
  be declared `static`. It it is needed in other files, then it must
  instead  be declared `extern`' in the matching header file.

- A `const` type should not be cast into a non-`const` type.
  * Exceptions are allowed when the loss of constedness is done for
    instrumentation purposes (for example adding performance or call
    counters to a class), or for synchronization purposes (e.g.
    updating an array's ghost cell values as a precaution where not
    changing the local values).

- Every `switch` construct should have a `default`
  clause (which may reduce to `assert(0)` to check code paths in
  debug mode, but at least this much must be ensured).

- A `const` attribute should be used when an array or structure
  is not  modified. Recall that for example `const cs_mesh_t *m`
  means that the contents of mesh structure `m` are not modified
  by the function, while `cs_mesh_t *const m` only means that
  the pointer to `m` is not modified; `const cs_mesh_t *const m`
  means both, but its usage in a function prototype gives no additional
  useful information on the function's side effects than the first form
  `const cs_mesh_t *m`, so that form is preferred, as it
  does not clutter the code.

- When an array is passed to a function, describing it as
  `array[]` is preferred to `*array`, as the array
  nature of the argument is made more visible.

- Where both a macro or an enumerated constant could be used,
  an enumeration is preferred, as values will appear with the
  enumerated value's name under a debugger, while only a macro's
  expanded value will appear. An additional advantage of enumerated
  values is that a compiler may issue a warning when a `switch`
  construct has no `case` for a given enumeration value.

Line spacing
------------

Try to be consistent regarding line spacing. Use 1 separation line to
highlight code blocks, but no more (i.e. avoid 2 successive blank lines).
Make the code symmetric also: when a blank line follows an opening brace
`{`, a similar line should precede the matching closing brace `}`.
When no such line is used after an opening brace, none should be used
either before the closing brace.

Separator lines
---------------

Double (`=====`) separator lines are used to separate main sections
of a file (e.g. first headers, than structure or class definitions,
then private, local functions, then public functions and classes).

Simple (`-----`) separators are used to separate function definitions.

Separator lines used to separate such sections must always end at column
79 when temrinated by a `*/` comment end, and column 77 otherwise, so
that the last `=` or `-` character is always at line 77.

Shorter separator lines (underlining a title) may be used inside
function bodies.

Separators define an implicit section hierarchy:
- File section level: full double (`=====`) line.
- Function or method separator: full single (`-----`) line.
- Main sections inside a function body: partial double line.
- Subsections inside a function: partial single line.

So a (`=====`) separated section should never be used
inside a section delimited by a single (`-----`) separator, as it
implies a higher-level section.

Language {#sec_prg_lang_cpp}
--------

For C++, the C++17 standard is required.

The choice of standard version is based on the availability of up to date
compilers on a wide range of machines, so upgrading to the C++ standard
will be considered only when C++20 compilers are available on all "old"
environments the code must run on (so probably not before 2028).
Note that building an up-to-date compiler could allow upgrading to a new
standard sooner, but this is not always easy with some organization's
IT rules.

For C code (e.g. in the preprocessor or PLE code), C11 or above is required.
C99 variable-length arrays should be avoided, as they are
an optional feature only.

Assertions
----------

Assertions are conditions which must always be verified. Several
expanded macro libraries may be available, but a standard C language
assertion has the following properties:

- It is only compiled in debug mode (and so incur no run-time
  performance penalty in production code, where the `NDEBUG`
  macro is defined).

- When its predicate is not verified, it causes a core dump;
  when running under a debugger, the code is stopped inside the
  assertion, but does not exit, which simplifies debugging.

Assertions are thus very useful to ensure that conditions
which are always expected (and not dependent on program input)
are met. They also make code more readable, in the sense that
it is made clear that conditions checked by an assertion
are always expected, and that not handling other cases is not an
programming error or omission.

If a condition may not be met for some program inputs,
and not just in case of programmer error, a more complete
test and call to an error handler (such as \ref bft_error)
is preferred.

Naming conventions
------------------

The following rules should be followed:

- Identifier names are in lowercase, except for macro or enumeration
  definitions, which are in uppercase; words in an identifier are
  separated by an underscore character (for example,
  `n_elt_groups_`.

- Global identifier names are prefixed by the matching library prefix,
  such as `cs_=` or `CS_`.

- Local (static) variable and function should be prefixed by an underscore
  character.

- Index arrays used with *0* to *n-1* (zero-based) numbering should
  be named using a `idx_` or `index_` prefix or suffix, while
  similar arrays using a *1* to *n* numbering (usually those that were
  also used in prior Fortran code) are named using a `pos_`
  prefix or suffix.

- In a similar manner, element identifiers should in general use
  a *0* to *n-1* (zero-based) numbering and be named using a `id_` prefix
  or suffix, while identifiers using a *1* to *n* numbering (usually those
  that were be also used in Fortran code) are named using a `num_`
  prefix or suffix.

Naming of enumerations
----------------------

The following form is preferred for enumerations:

```{.c}
typedef myclass { CS_MYCLASS_ENUM1,
                  CS_MYCLASS_ENUM2,
                  /* etc. */
                } cs_myclass_t;
```

Naming of structures and associated functions
---------------------------------------------

Macros and enumerations related to `myclass` structures
are prefixed by `CS_MYCLASS_`.

Public functions implementing methods are named
`cs_`*class_method*, while private functions are simply named:
*class_method* and are declared static.

Files containing these functions are named *class*`.c`.

Integer types {#sec_prg_lang_integer_types}
-------------

Several integer types are found in code_saturne:

- `cs_lnum_t` should be used for local entity (i.e. vertex, face,
   cell) numbers or connectivity. It is a signed integer, normally identical
   to `int`, but a larger size could be used in the future for very
   large meshes on shared memory machines.

- `cs_gnum_t` should be used for global entity numbers, usually
   necessary only for I/O aspects. It is an unsigned 32 or 64-bit integer,
   depending on whether the code was configured with the
   `--enable-long-gnum` option. Global numbers should always use
   this type, as for very large meshes, they may exceed the maximum size
   of a 32-bit integer (2 147 483 648). The choice of unsigned integers
   is two-fold: it doubles the range of available values, and good compilers
   will issue warnings when this type is mixed without precaution with
   the usual integer types. These warnings should be heeded, as they may
   avoid many hours of debugging.

- In all other cases, the standard C types `int` and `size_t`
  should be preferred (for example for loops over variables, probes, or
  any entity independent of mesh size.

Base functions and types
------------------------

In the code_saturne kernel, it is preferable to use base functions provided by the
BFT subsystem to the usual C functions, as those logging, exit and error-handling
functions will work correctly when running in parallel, and the memory
management macros ensure return value checking and allow additional logging.

The array below summarizes the replacements for usual functions:

Standard C function         | code_saturne macro or function | Header
----------------------------|--------------------------------|------------------
`exit`                      | \ref cs_exit                   | \ref cs_base.h
`exit` on error             | \ref bft_error                 | \ref bft_error.h
`printf`                    | \ref bft_printf                | \ref bft_printf.h
`printf` (to standard logs) | \ref cs_log_printf             | \ref cs_log.h
`malloc`                    | \ref CS_MALLOC                 | \ref cs_mem.h
` `                         | \ref CS_MALLOC_HD              | \ref cs_base_accel.h
`realloc`                   | \ref CS_REALLOC                | \ref cs_mem.h
` `                         | \ref CS_REALLOC_HD             | \ref cs_base_accel.h
`free`                      | \ref CS_FREE                   | \ref cs_mem.h
.h

Internationalization
--------------------

Messages for the user should always be defined in US English in the source
code (which avoids using extended characters and the accompanying
text encoding issues in source code).

To make future internationalization using a mechanism such as
`gettext()` possible, translatable strings should be
encased in a `_( )` macro (actually an abbreviation for a call
to `gettext()` if available, which reverts to an empty (identity)
macro is internationalization is unavailable or disabled).
Strings assigned to variables must be encased in a `N_( )`
macro (which is an ``empty'' macro, used by the `gettext` toolchain
to determine that those strings should appear in the translation dictionary),
and the variable to which such a string is assigned should be
encased in the `_( )` macro where used.

Note that with UTF-8 type character strings, accented or otherwise
extended characters are represented on multiple bytes. The `strlen`
C function will return the string's real size, which may be greater
than the number of output columns it uses. The \ref cs_log_strlen
function may be used to compute the printable width of a character
string, while and \ref cs_log_strpad and \ref cs_log_strpadl may be use
to pad a string.

C++ streams
-----------

Note that standard C++ stream output (using `std::cout` ) may be practical
for temporary debugging code, but should not be used instead of
`bft_printf` or `cs_log` for production code.

To use streams for production output, the following conditions must
be met first:
- Add a stream class whith configurable parallel behavior, similar to
  `bft_printf` and `cs_log(CS_LOG_DEFAULT, ...` so that output on MPI
  ranks > 0 can be either dropped or sent to separate files.
- Although internationalization is not activated since v7.0, the
  possibility of reactivating it is not currently ruled out. `printf`
  style functions can be internationalized with mechanisms such as
  `gettext`, which allows changing the order of some output arguments
  through the format string if language differences require it.
  Using streams, this is not possible. This may be possible using
  `std::format`, but this requires C++20 compatibility.

Defining a C++ style output stream for code_saturne will be possible
when the code switches to C++ 20 requirements.

C++ file template {#sec_prg_cpp_file_template}
-----------------

The following template illustrates a typical C++ file's preferred
organization, and should be followed so as to keep the code consistent.

```
/*============================================================================
 * Short description.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-<current_year> EDF S.A.

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
*/

/*----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>  // example
#include <string>    // example

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

// For example
#include "base/cs_math.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "section/current_file.h"

/*=============================================================================
 * Additional Doxygen documentation
 *============================================================================*/

/*!
  \file <file_name>
        Brief description.

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local macro definitions
 *============================================================================*/

/* Example local macro */

#define _CS_EXAMPLE_MACRO  16

/*============================================================================
 * Type definitions
 *============================================================================*/

// Class and structure definitions go here.

/* Field key definitions */

struct example_loca_struct {

  <type>  m1;        /*!< Description */
  <type>  m2;        /*!< Description */
  ...

};

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Static global variables (visible only from this file)
   go here, with `_` prefix. */

/*============================================================================
 * Global variables
 *============================================================================*/

// Global variables go here, with `cs_glob_` prefix.

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*
 * Private function description, legacy style.
 *
 * parameters:
 *   p1  <-- input parameter
 *   p2  <-> input/output parameter
 *   p3  --> output parameter
 *
 * returns  description of return output
 *----------------------------------------------------------------------------*/

static <return type>
_function_a(<type>  p1,
            <type>  p2,
            <type>  p3)
{
  ... // function body
}

/*----------------------------------------------------------------------------*/
/*!
 * Private function description, Doxygen style.
 *
 * \param[in]       p1  input parameter
 * \param[in, out]  p2  input/output parameter
 * \param[ouyt]     p3  output parameter
 *
 * \return  description of return output
 */
/*----------------------------------------------------------------------------*/

static <return type>
_function_b(<type>  p1,
            <type>  p2,
            <type>  p3)
{
  ... // function body
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * Include Doxygen documentation here
 */
/*----------------------------------------------------------------------------*/

return_type
function_a(...  // arguments
           )
{
  // function body
}

/*----------------------------------------------------------------------------*/
/*!
 * Include Doxygen documentation here
 */
/*----------------------------------------------------------------------------*/

return_type
function_b
(
  ...  // arguments
)
{
  // function body

  /* First section
     ============= */

  ...

  /* Subsection a
     ------------ */

  ...

  /* Second section
     ============= */

  ...
}

/*----------------------------------------------------------------------------*/
/*!
 * Include Doxygen documentation here
 */
/*----------------------------------------------------------------------------*/

return_type
class::member_function_a
(
  ...  // arguments
)
{
  // function body
}

/*----------------------------------------------------------------------------*/
```

Automatic formatting
--------------------

A `.clang-format` configuration file is available in the top-level
source directory, so the [clang-format](https://clang.llvm.org/docs/ClangFormat.html)
tool may be used to obtain an acceptable presentation.
Automated formatting does not exactly follow the above recommendations,
but is as close as possible using the available parameters, so can provide
a good starting point, though it is not activated automatically to allow for
manual fine-tuning.

No automatic formatter tested has been able to match the required
coding style without breaking consistency, because formatters
cannot handle the alignment of math expressions.
LLMs may be much better at following the style recommendations.

Python coding style
===================

The core python scripts use a coding style similar to PEP-8
[PEP-8](https://www.python.org/dev/peps/pep-0008/).
PEP-257 (docstring conventions) is also recommended.

Other parts of the code tend to use a camelCase naming, but should otherwise
adhere to the same standards. Moving the to PEP-8 style would be ideal,
though to avoid confusion, this should be done in an *atomic* step
for each module.

Documentation
=============

Documentation of the main code is based on the
[Doxygen](https://www.doxygen.nl/index.html/) tool, whose documentation
may be found on its web site.

Additional pages for the documentation may be found in the source tree,
under `docs/doxygen`. Files containing mostly examples may use the `.h` or
`.dox` extension (with `.dox` preferred for easier identification), and
pages which describe general aspects instead of code are preferably
written in Doxygen Markdown (`.md` extension), as this allows better
readability and interoperation with some editors (such as preview,
syntax highlighting, ...).

When building the code, remember to often use `make html` and check the
`docs/doxygen/doxygen_warn.log` file in the build directory for errors
and warnings.

When modifying arguments to a function or modifying structures, make sure
the special Doxygen comments are kept up to date. Using Doxygen
comments in the C code and simplified comments in the headers
(see \ref cs_medcoupling_postprocess.cxx and \ref cs_medcoupling_postprocess.h for example)
is recommended, but as this
adds to the coding effort, duplicating the headers from the reference
C code to the headers is not recommended since recent versions of Doxygen complain
about this, so avoiding duplicates is the preferred way.
For inlined functions or structures the documentation should be available in the
.h files since that is where the function or structure are defined.

For functions the following structure is advised (example extracted from
\ref cs_medcoupling_postprocess.cxx):
```
cs_real_t
cs_medcoupling_slice_scalar_mean_weighted
(
  const char        *name,     /*!<[in] name of the slice */
  const cs_real_t   *scalar,   /*!<[in] array of scalar values (size n_cells) */
  const cs_real_t   *weight_s, /*!<[in] scalar weight array (n_cells) */
  const cs_real_3_t *weight_v  /*!<[in] vectorial weight array (n_cells) */
)
{
  <Function code>
}
```

Private functions or structures should not appear in the documentation
(though their arguments should be documented in the source code),
so in most source files,
```
/*! \cond DOXYGEN_SHOULD_SKIP_THIS */
```
is used to mark the beginning of a section which should be ignored
by Doxygen, and
```
/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */
```
used to mark the end of that section. In most cases, this includes private
structures and functions in C code, but could be extended to public
function definitions if the `.h` file header already contains the same
Doxygen-formated comments.
