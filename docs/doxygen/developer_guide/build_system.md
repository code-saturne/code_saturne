<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

\page cs_dg_build_system Build system

[TOC]

In this section, some specific aspects of the build system are described.

Base tools {#sec_prg_build_base}
==========

The code_saturne build system is based mostly on the
[GNU Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html#Autotools-Introduction),
and provides a traditional `configure`, `make`, `make install` sequence).

Those tools have their sets of advantages and disadvantages, but are very
well documented. An very practical feature, not available in such a natural
form in some other build frameworks, is the possibility to provide a
quick and complete list of available options using `configure --help`.

Specific adaptations
====================

- To rebuild or update configure scripts and associated Makefile templates,
  run `./sbin/bootstap` from the top source directory. Some scripts may also
  be rebuilt automatically when the Autotools *maintainer mode* is enabled
  (which is the default when sources are pulled from Git), but to ensure
  a complete update, the above command is recommended.

- To allow detecting compilers used and provide default flags, a
  `config/cs_auto_flags.sh` script is run by `configure`.
  Modifying this file does not require running `./sbin/bootstap`, so
  it is possible to adapt it easily when porting to a new machine or set
  of compilers, even on cluster environments where the Autotools might
  not be available.

- Though Autotools allow both in-tree and out of tree builds, only
  out-of tree builds are recommended and supported for code_saturne, as
  specified in the install documentation. In-tree builds might work,
  but they will pollute the source tree with generated files, prohibit
  multiple builds (such as side-by-side debug and production builds).
  No precaution is taken in code_saturne to avoid name conflicts between
  provided files and generated files, since this is not an issue with
  out-of-tree builds.

- Since Automake has only partial support for modern Fortran, and does
  not handle module dependencies, all Fortran module dependencies must
  be handled explicitly. Using recursive make, they are build first
  based on rules in the `src/Makefile.am` Automake template
  (whereas most files to be compiled are listed in sub-directories).
  As Fortran is progressively being replaced by C in code_saturne, new
  modules are rarely added, limiting the inconvenience.

- Since code_saturne is managed through Python scripts, many settings
  detected at configure and build time are saved for use in Python through
  some specific files, output by `configure` (`bin/cs_config.py` and the
  intermediate `config/code_saturne_build.cfg.in` file) and `make install`
  (`bin/cs_package.py` and the final `config/code_saturne_build.cfg` file).
  Currently, most settings are used directly from the installed
  `cs_config.py` and `cs_package.py`, but migrating to import of settings
  from `config/code_saturne_build.cfg` is recommended.

- Libtool has some nice features relative to building and organizing
  libraries, but unfortunately decides that it is always right, even when
  it is wrong. There are no options to override some of its settings. So to
  avoid a long series of issues related to libtool, its use must be avoided
  when linking the final executable. This is detailed in a section below.

- In addition, `libtool` does not handle some compiler options, such
  as the CUDA compiler. To allow building CUDA code with dynamic libraries,

- Since code_saturne allows for user-defined functions and subroutines,
  a specific solution is required for those.

- All *m4* files specific to code_saturne are found in `m4`, and prefixed
  by `cs_` (so may be listed as `m4/cs_*.m4`). Other m4 files in the same
  directory are generated or copied by `automake`.

Compiling the main solver executable
------------------------------------

The main solver executable is usually `cs_solver`, but variants such
as `nc_solver`, or other executables for unit testing are handled in
a similar manner.

The main executable is compiled using the Python scripts. When installed or
installing, the appropriate code is found in `cs_compile.py` (under
`bin` in the code sources, and Python package install path in installed
code). When building, since the available path structure is different,
some methods are overloaded, using code from `build-aux/cs_compile_build.py`.

Code version determination
--------------------------

To avoid requiring extra steps when releasing a version, which have proved
to be error-prone in the past, the code version handling is semi-automatic.

Since it is good practice to provide some form of release notes, and such
notes require human intervention, the version number is based on entries
available in the `NEWS.md` file in the top source directory. This is
completed by info available from Git, if available. Based on this information,
a `build_aux/cs_version.py` script determines the major, minor, and
release numbers,n possibly with an additional Git hash, and a flag indicating
whether the version based on Git is further modified (i.e. if there are
uncommitted changes).

So when editing the `NEWS.md`, care must be taken to follow the existing
format and naming conventions.

The following logic is also applied in the `build_aux/cs_version.py` script:

- For each major code version *x*, there are 4 minor versions *y* (0, 1, 2, 3).
  When *y* would reach 4, *x* is incremented, and *y* set to 0.

- In non-release branches, such as the master branch, the version number
  is based on the version following the last release branch (so 6.2 is followed
  by 6.3, 6.3 is followed by 7.0, etc.), to which the *-alpha* qualifier
  is added.

- In release branches:
  - Prior to an *x.y.0* release, the version number is
    automatically suffixed by *beta*,  unless a tag containing *-rc* is
    found in Git, in which case this is used instead (for release candidates).
  - After a *x.y.0*, and while *x.y.z* is marked as *unreleased* in `NEWS.md`
    (rather than providing a release date), the version will appear as
    *x.y.(z-1)-patch*.

  Roadmap
  -------

  The GNU Autotools provide some nice features, and the code_saturne build
  system represents a large amount of work over multiple years, but this system
  is showing its age, and perhaps not evolving well enough to suit current
  needs.

  Worse, these tools have lacked consistent integration from the start. For
  example, Autoconf allows generating multiple `config.h` files (such as
  `cs_config.h` and `cs_config_priv.h`), but Automake can only handle one.
  This has required writing work-arounds. Libtool possibly causes more
  issues than it solves. And finally, its cross-platform aspects are very
  limited. So Windows ports are basically limited to the Windows Subsystem
  for Linux.

  Alternatives also have their advantages and disadvantages:

  - [CMake](https://cmake.org/) is very complete, cross-platform, and used by
    many software projects, but is also quite complex.
    - It does not have a nice equivalent to `configure --help`. Though
      `cmake -LAH` is similar, it immediately caches some results, and should
      thus be run in a specific directory, making it a bit more tricky to use.
    - It requires learning a specific language. Perhaps not more complicated
      than the *m4* scripts used in Autotools, but still one more than the
      main programming languages used in code_saturne.
    - We have not tested/evaluated the generation of configuration files
      (such as `config/code_saturne_build.cfg`) yet.
    - The caching behavior can lead to subtle issues.
      - With the GNU autotools, it is easy to obtain the full configuration
        command from `config.log` or `config.status`, adapt it if necessary,
        empty the build directory, and generate a new configuration. Doing the
        equivalent with CMake does not seem so straightforward.
    - It requires a CMake installation (though it is now almost ubiquitous).
    - On the advantages side, it does not require regenerating `configure`
      and `Makefile.in` files (using `sbin/bootstrap` in code_saturne) when
      CMake files are modified.

  - [Meson](https://mesonbuild.com/) seems like an interesting system,
    used by some large projects, and has the advantage of being based on
    Python, which is already familiar to most code_saturne developers.
    - A disadvantage is that its reliance on the
      [Ninja](https://ninja-build.org/) build system, adds an installation
      and build dependency.

  - Pure Python, such as used in the [PETSc](https://www.mcs.anl.gov/petsc/)
    build system, is interesting, but may require a significant coding
    investment.

For now, a safe solution seems to be to rewrite some parts of the build system
in Python, as is already done in quite a few areas (see files in `build-aux`),
and progressively reduce the reliance on the Autotools.
