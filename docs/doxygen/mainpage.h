/*============================================================================
 * Main Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

/*-----------------------------------------------------------------------------*/

/*!
  \mainpage Code_Saturne documentation

  \section mainpage_intro Introduction

  Code_Saturne is EDF's general purpose Computational Fluid Dynamics (CFD)
  software.

  The basic capabilities of Code_Saturne enable the handling of either
  incompressible or expandable flows with or without heat transfer and
  turbulence. Dedicated modules are available for specific physics such as
  radiative heat transfer, combustion (gas, coal, heavy fuel oil, ...),
  magneto-hydrodynamics, compressible flows, two-phase flows
  (Euler-Lagrange approach with two-way coupling), or atmospheric flows.

  \section mainpage_install Installation

  Code_Saturne may be configured and installed using the
  \c configure shell script and \c make. Please read the \c INSTALL
  file in the toplevel source directory if you are not familiar
  with configuration scripts generated through GNU autoconf and automake.

  \section mainpage_sublibs Sub-libraries

  See \ref mainpage_ple "PLE" (Parallel Location and Exchange),
  \ref bft_page and \ref fvm_page

  \section mainpage_user_sources User sources, functions, and subroutines

  Many user examples are available in the \ref cs_user_examples "user examples tab":
    - \ref cs_user_boundary_conditions_examples "User boundary conditions definitions",
    - \ref cs_user_extra_operations_examples "User extra operations",
    - \ref cs_user_initialization "Physical fields user initialization",
    - and so on ...

  \section additional_doc Additional documentation

  In addition to the Doxygen documentation, Code_Saturne is provided with six
  pdf documents:
    - a <a href="../../user.pdf"><b>user guide</b></a>,
    - a <a href="../../theory.pdf"><b>theory guide</b></a>,
    - a <a href="../../developer.pdf"><b>developer guide</b></a>,
    - a <a href="../../studymanager.pdf"><b>guide of studymanager</b></a>,
    - an <a href="../../install.pdf"><b>installation guide</b></a>,
    - a <a href="../../refcard.pdf"><b>refcard</b></a>.

  \page bft_page BFT

  \section bft BFT: Base Functions and Types

  The "Base Functions and Types" library is intended to simplify and enhance
  portability, memory and I/O use for scientific codes. It contains
  a number of system library wrappers for common services such as file
  I/O or memory management, ensuring portability of the calling code.

  \subsection bft_intro_portability Portability

  Using lower-level services in C or C++ often requires the definition
  of preprocessor symbols such as \c _POSIX_SOURCE (or \c OSF_SOURCE, or
  \c HPUX_SOURCE, ...). On certain systems, largefile support also requires
  additionnal preprocessor symbols such as \c _FILE_OFFSET_BITS or
  \c _LARGEFILE_SOURCE, and \c fseek/ftell replaced with \c fseeko/ftello.

  Authors of scientific code seeking portability should not have to worry
  about these issues unless they deliberately choose to use low-level functions.
  BFT tries to hide such portability issues while maintaining an API similar
  to that of the standard \c libc where applicable.

  \subsection bft_intro_retcode_error Return codes and error handling

  In most scientific codes destined to run in a batch environment, errors are
  usually fatal, especially when dealing with file access and memory allocation.
  The functions provided by the BFT library always check for return codes,
  and call the appropriate error handler when applicable. The default is
  to terminate the running program after printing the appropriate message,
  but the user may define and set other error handlers with different
  behavior.

  Error handling may also be modified by writing a specific error handler
  (see bft_error_handler_example).

  \subsection bft_intro_add_func Added functionnality

  BFT functions similar to \c libc functions add functionnality such
  as optional byte-swapping for conversion from internal to external
  data repressentation (or vice-versa), or optional memory-allocation
  logging and tracking of non-freed pointers.

  \subsection bft_intro_goals Goals and Limitations

  The BFT library tries to provide a set of utilitarian functions for
  common use, but does not seek to define a framework. As a general
  rule, functions provided by BFT should provide added portability
  or functionnality when compared to their \c libc or Posix counterparts
  (when such counterparts exist), as simple wrapping with no added
  functionnality only makes code less readable to an experienced developer
  and is to be avoided.

  Subsets of BFT may be used independently if desired, and are orthogonal,
  except as regards error handlers. With non-default error handlers, they
  can be made fully orthogonal. Only certain subsets may be used if
  preferred.

  The BFT library provides memory-usage measurement functions, whose
  implementations are system-dependent. If it has not yet been ported
  to a given type of environment, these functions should return 0.
  The user should thus check for the return values of such functions,
  but the API is guaranteed.

  \page fvm_page FVM

  \section fvm FVM: Finite Volume Mesh

  The "Finite Volume Mesh" library is intended to provide mesh and associated
  fields I/O and manipulation services for unstructure Finite Volume codes
  or other tools with similar requirements.

  \subsection fvm_intro_goals Goals and Limitations

  The FVM library is originaly intended for unstructured cell-centered
  finite-volume codes using almost arbitrary polyhedral cells. It
  may thus handle polygonal faces (convex or not) and polyhedral cells
  as well as more classical elements such as tetrahedra, prisms, pyramids,
  and hexahedra, but can currently only handle linear elements.
  There are also currently no optimizations for structured or
  block-structured meshes, wich are handled as unstructured meshes.

  \subsection fvm_nodal_writer Nodal Mesh and Writer Structures

  FVM is used to handle post-processing output in external formats
  (EnSight, CGNS, or MED are currently handled), possibly running in
  parallel using MPI. In the case of Code_Saturne, this implies reconstructing
  a nodal connectivity (cells -> vertices) from the faces -> cells
  connectivity (using an intermediate cells -> faces representation,
  passed to th FVM API). It is also possible to directly pass a nodal
  connectivity to FVM.

  So as to limit memory usage and avoid unecessary copies, the
  \c fvm_nodal_t structure associated to a mesh defined by its nodal
  connectivity is construed as a "view" on a mesh, and owns as little
  data as possible. As such, most main structures associated with
  this representation are defined by 2 arrays, one private, and one
  shared. For example, an fvm_nodal_t structure has 2 coordinate
  arrays:

  <tt>const *cs_coord_t  vertex_coords;</tt>

  <tt>*cs_coord_t       _vertex_coords;</tt>

  If coordinates are shared with the calling code (and owned by that
  code), we have  <tt>_vertex_coords = NULL</tt>, and \c vertex_coords
  points to the array passed from the calling code. If the
  coordinates array belongs to FVM (either having been "given" by
  the calling code, or generated by an operation wich invalidates
  coorindates sharing with the parent mesh, such as mesh extrusion),
  we have <tt>vertex_coords = _vertex_coords</tt>, which points to the
  private array.

  When an fvm_nodal_t object is destroyed, it destroys its private
  data and frees the corresponding memory, but obviously leaves its
  shared data untouched.

  If a \c fvm_nodal_t structure B is built from a structure A with
  which it shares its data, and a second \c fvm_nodal_t mesh C
  is a view on B (for example a restriction to a part of the domain
  that we whish to post-process more frequently), C must be destroyed
  before B, which must be destroyed before A.
  FVM does not use reference counters or a similar mechanism, so
  good managment of object lifecycles is of the responsiblity of
  the calling code. In practice, this logic is simple enough that
  no issues have been encountered with this model so far in the
  intended uses of the code.

  Another associated concept is that of "parent_number": if a mesh
  constitutes a "view" on another, sur un autre, a list of
  parent entity numbers allows accessing variables associated
  with the "parent" mesh, without needing to extract or duplicate
  values at the level of the calling code.

  Note that an FVM structure is global from an MPI point of
  view, as it may participate in collective parallel operations.
  Thus, if an FVM mesh represents a subset of a global domain,
  it may very well be empty for some ranks, but it must still exist
  on all ranks of the main communicator associated with FVM.

  For parallelism, another important concept is that of "global numbering",
  corresponding to entity numbers in a "serial" or "global" version
  of an object: two entities (vertices, elements, ...) on different
  processors having a same global number correspond to the same
  "absolute" object. Except when explicitely mentioned, all other
  data defining an object is based on local definitions and numberings.
  Parallelism is thus made quite "transparent" to the calling code.

  In parallel, post-processor output (using EnSight Gold, CGNS, or MED
  formats) uses the global numbering to serialize output, and generate
  a unique data set (independent of the number of processors used for
  the calculation). This operation is done by blocks so as to
  limit local memory consumption, possibly incurring a slight processing
  overhead. The aim is to never build a full global array on a single
  processor, so as to ensure scalability. This is fully implemented
  at least for the EnSight Gold format.

  As an option, generic polygonal or polyhedral elements may be
  split into simpler elements (triangles, tetrahedra, and pyramids)
  for output, this being done on the fly (so as to avoid a complete
  memory copy). This allows working around possible lack of support for
  these complex elements in certain tools or formats.

*/
