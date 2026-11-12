General Information
===================

The PLE (Parallel Location and Exchange) library is designed to simplify
coupling of distributed parallel computational codes. It is maintained a
a part of [code_caturne](https://code-saturne.org) CFD tool,
EDF's general purpose Computational Fluid Dynamics (CFD) software,
but it may also be used with other tools, and is distributed under
a more permissive licence (LGPL instead of GPL).

PLE provides support for 2 categories of tasks: synchronizing
parallel codes at predifined points, and enabling parallel mapping
of points to meshes, and transfer of variables using this mapping.

Its target use is coupling of parallel tools based using
[MPI](https://www.mpi-forum.org/), though it could be adapted to other
runtimes.

Copying
-------

The PLE library is distrubuted under the GNU Lesser General Public Licence.
See the COPYING file for details.

Installation
------------

The installation system is based on the GNU autotools (without libtool)

If the code was obtained through a Git repository, a first step is required:
```
cd libple
./sbin/bootstrap
cd ..
```

To create a standalone installation of PLE, it is recommended
to create a specific build directory, and run PLE's configure script
from there:

```
mkdir ple_build
cd ple_build
<path_to_ple_sources>/configure
make
make install
```

The "configure" script takes several options, so running:
`<path_to_ple_sources>/configure --help`
first is recommended to list the available options.

To port PLE to a system on which it was not previously supported,
most of the adaptation should be possible by adapting the
`config/ple_auto_flags.sh` file, which does not require running
`./sbin/bootstrap` again.