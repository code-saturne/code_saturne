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

\page cs_ug_mesh_prepare Mesh import and preprocessing

[TOC]

This page provides additional details on mesh import and preprocessing steps,
as well as how mesh regions can be managed in the code_saturne Solver.

Referring to the modules described in the [base architecture](@ref main_components)
presentation, a mesh can be imported using the code_saturne Preprocessor,
than further preprocessed in the Solver's initialization steps.

Importing meshes is done by the Preprocessor module, while most other
preprocessing operations are done mainly by the code Solver
(except for element orientation checking, which is done by the Preprocessor).

The Preprocessor module of code_saturne reads the
mesh file(s) (under any supported format) and translates the necessary
information into a Solver input file.

When multiple meshes are used, the Preprocessor is called once per mesh,
and each resulting output is added in a `mesh_input`
directory (instead of a single `mesh_input.csm`  file).

\subpage cs_ug_mesh_select_c

Throughout the use of code_saturne, including the selection of specific
mesh sections for preprocessing operations, it is also useful to understand
how high-level [mesh element selection](@ref cs_ug_mesh_select_c) can be handled.

Preprocessor module {#sec_preprocessor_module}
===================

The code_saturne run operation automatically calls the Preprocessor for every
mesh specified by the user with the GUI or in the `MESHES` list
defined in `cs_user_scripys.py:define_domain_parameters`.

Preprocessor options {#sec_prg_optappelecs}
--------------------

The executable of the Preprocessor module is named `cs_preprocess`, and
is normally called through the run script, so it is not fount
in the main executable program search paths.

Nonetheless, it may be useful to call the Preprocessor manually
in certain situations, especially for frequent verification when
building a mesh, so its use is described here. Verification
may also be done using the GUI or the mesh quality check mode
of the general run script.

The Preprocessor is controlled using command-line arguments.
A few [environment variables](@ref sec_env_var_pcs) allow advanced users
to modify its behavior or to obtain a trace of memory management.
If needed, it can be called directly on most systems using the
following syntax:
```
${install_prefix}/libexec/code_saturne/cs_preprocess [options] <mesh>
```
Note that on some systems (such as OpenSUSE), `libexec may be replaced by `lib`.

To obtain a description of all available option, type:
```
${install_prefix}/libexec/code_saturne/cs_preprocess --help
```
In the following sections, the `${install_prefix}/libexec/code_saturne`
prefix is omitted for conciseness.
Main choices are done using command-line options. For example:

`cs_preprocess --num 2 fluid.med`

means that we read the second mesh defined in the `fluid.med` file, while:

`cs_preprocess --no-write --post-volume med fluid.msh`

means that we read file `fluid.msh`, and do not produce a
`mesh_input.csm`  file, but do output `fluid.med` file
(effectively converting a Gmsh file to a MED file).

### Mesh selection {#sec_optpcs_mesh}
<!-- -->

Any use of the Preprocessor requires one mesh file (except for `cs_preprocess`
and `cs_preprocess -h` which respectively print the version number and
list of options). This file is selected as the last argument to `cs_preprocess`,
and its format is usually automatically determined based on its
[extension](@ref sec_mesh_viz_formats), but a `--format` option allows forcing
the format choice of the selected file.

For formats allowing multiple meshes in a single file, the
`--num` option followed by a strictly positive integer allows
selection of a specific mesh; by default, the first mesh is selected.

For meshes in [CGNS](@ref sec_fmtdesc_cgns) format, we may in addition use the
`--grp-cel` or `--grp-fac` options, followed by the `section`
or `zone` keywords, to define additional groups of cell or faces
based on the organization of the mesh in sections or zones. These sub-options
have no effect on meshes of other formats.

### Post-processing output {#sec_optpcs_post}
<!-- -->

By default, the Preprocessor does not generate any post-processor output.
By adding `--post-volume [format]`,
with the optional `format` argument being one of `ensight`,
`med`, or `cgns` to the command-line arguments, the output of the volume
mesh to the default or indicated format is provoked.

In case of errors, output of error visualization output is always
produced, and by adding `--post-error [format]`,
the format of that output may be selected (from one of `ensight`,
`med`, or `cgns`, assuming [MED](@ref sec_fmtdesc_med) and
[CGNS](@ref sec_fmtdesc_cgns) are available),

### Element orientation correction {#sec_optpcs_orient}
<!-- -->

Correction of element orientation is possible and can be activated using the
`--reorient` option.

Note that we cannot guarantee correction (or even detection) of a bad
orientation in all cases.
Not all local numbering possibilities of elements are tested,
as we focus on "common" numbering permutations. Moreover,
the algorithms used may produce false positives or fail to find
a correct renumbering in the case of highly non convex elements.
In this case, nothing may be done short of modifying the mesh, as
without a convexity hypothesis, it is not always possible to choose
between two possible definitions starting from a point set.

With a post-processing option such as `--post-error`
or, `--post-volume`, visualizable meshes of corrected elements as
well as remaining badly oriented elements are generated.

### Optional functionality {#sec_pcs_lib_opt}
<!-- -->

Some functions of the Preprocessor are based on external libraries,
which may not always be available. It is thus possible to configure
and compile the Preprocessor so as not to use these libraries.
When running the Preprocessor, the supported options are printed.
The following optional libraries may be used:

* CGNS library. In its absence, [CGNS](@ref sec_fmtdesc_cgns)
  format support is deactivated.

* med-file library. In its absence, [MED](@ref sec_fmtdesc_med)
  format support is simply deactivated.

* libCCMIO library. In its absence, [STAR-CCM+](@ref sec_fmtdesc_ccm)
  format support is simply deactivated.

* Read compressed files using Zlib. With this option, it is
  possible to directly read mesh files compressed with a
  *gzip* type algorithm and bearing a `.gz` extension.
  This is limited to formats not already based on an external
  library (i.e. it is not usable with CGNS, MED, or CCM files),
  and has memory and CPU time overhead, but may be practical.
  Without this library, files must be uncompressed before use.

General remarks
---------------

Note that the Preprocessor is in general capable of reading all "classical"
element types present in mesh files (triangles, quadrangles, tetrahedra,
pyramids, prisms, and hexahedra).
Quadratic or cubic elements are converted upon reading into their
linear counterparts. Vertices referenced by no element (isolated vertices
or centres of higher-degree elements) are discarded. Meshes are read
in the order defined by the user and are appended, vertex and element
indices being incremented appropriately.
Entity labels (numbers) are not maintained, as they would probably not be
unique when appending multiple meshes.

At this stage, volume elements are sorted by type, and the fluid domain
post-processing output is generated if required.

In general, groups assigned to vertices are ignored.
selections are thus based on faces or cells. with tools such
as Simail, faces of volume elements may be referenced directly, while
with I-deas or SALOME, a layer of surface elements bearing the required
colors and groups must be added. Internally, the Preprocessor always considers
that a layer of surface elements is added (i.e. when reading a Simail
mesh, additional faces are generated to bear cell face colors.
When building the \f$faces \rightarrow cells\f$ connectivity, all faces with the
same topology are merged: the initial presence of two layers of identical
surface elements belonging to different groups would thus lead to
a calculation mesh with faces belonging to two groups).

Files passed to the Solver {#sec_pcs_mode_comm}
--------------------------

Data passed to the Solver by the Preprocessor is transmitted using a
binary file, using a *big endian* data representation, named
`mesh_input.csm` (or contained in a `mesh_input` directory).

When using the Preprocessor for mesh verification, data for the Solver
is not always needed. In this case, the `--no-write` option may
avoid creating a Preprocessor output file.

Mesh preprocessing {#sec_prepro}
==================

Meshing remarks {#sec_prg_meshes}
---------------

\warning
Some turbulence models (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$, SSG, ...)
used in code_saturne are "High-Reynolds" models. Therefore the size of the cells
neighboring the wall must be greater than the thickness of the viscous
sub-layer (at the wall, \f$y^+>2.5\f$ is required, and \f$30<y^+<100\f$ is
preferable). If the mesh does not match this constraint, the results may
be false (particularly if thermal phenomena are involved). For more details
on these constraints, see the [iturb](@ref iturb) keyword.

Mesh joining {#sec_optpcs_join}
------------

Conforming joining of possibly non-conforming meshes may be done by the
solver, and defined either using the Graphical User Interface (GUI) or the
`cs_user_join` user function. In the GUI, the user must
add entries in the "Face joining" tab in the "Mesh/Proprocessing" page.
The user may specify faces to be joined, and can also modify basic joining
parameters.

\anchor fig_joining
\image html gui_mesh_join.png "Defining a mesh joining"

For a simple mesh, it is rarely useful to specify strict face selection
criteria, as joining is sufficiently automated to detect which faces
may actually be joined. For a more complex mesh, or a mesh with thin
walls which we want to avoid transforming into interior faces, it is
recommended to filter boundary faces that may be joined by using
face selection criteria. This has the additional advantage of reducing
the number of faces to test for in the intersection/overlap search,
and thus improving the performance of the joining algorithm.

One may also modify tolerance criteria using 2 options, explained
in the [figure](@ref fig_join_tolerance) below:

* `fraction` assigns value \f$r\f$ (where \f$0 < r < 0.49\f$) to
   the maximum intersection distance multiplier (\f$0.1\f$ by default).
   The maximum intersection distance for a given vertex is based on the
   length of the shortest incident edge, multiplied by *r*. The maximum
   intersection at a given point along an edge is interpolated from that
   at its vertices, as shown on the left of the figure
* `plane` assigns the maximum angle between normals for two faces
   to be considered coplanar (\f$25^{\circ}\f$ by default); this parameter is
   used in the second stage of the algorithm, to reconstruct conforming
   faces, as shown on the right of figure

\anchor fig_join_tolerance
\image html join_tolerance.svg "Maximum intersection tolerance and faces normal angle"

As shown above, verbosity and visualization levels can also be set,
to obtain more detailed logging of the joining operations, as
well as detailed visualizable output of the selected and joined faces.

In practice, we are sometimes led to increase the maximum intersection
distance multiplier to *0.2* or even *0.3* when joining curved surfaces,
so that all intersection are detected. As this influences merging
of vertices and thus simplification of reconstructed faces, but also
deformation of "lateral" faces, it is recommended only to modify it
if necessary. As for the `plane` parameter, its use has
only been necessary on a few meshes up to now, to
reduce the tolerance so that face reconstruction does not
try to generate faces from initial faces on different surfaces.

This operation can be specified either using the GUI or the
user-defined functions, as shown in various
[examples](@ref cs_user_mesh_h_cs_user_mesh_joining).

### Advanced parameters for mesh joining
<!-- -->

Advanced parameters may be set through user-defined functions by calling
\ref cs_join_set_advanced_param, as shown in the
[advanced mesh joining parameters](@ref cs_user_mesh_h_cs_user_mesh_add_advanced_joining) examples section.

Periodicity {#sec_optpcs_period}
-----------

Handling of periodicity is based on an extension of conforming joining,
as shown on the following [figure](@ref fig_join_periodic).
It is thus not necessary for the periodic faces to be conforming (though
it usually leads to better mesh quality). All options relative to conforming
joining of non-conforming faces also apply to periodicity. Note also that
once pre-processed, 2 periodic faces have the same orientation
(possibly adjusted by periodicity of rotation).

\anchor fig_join_periodic
\image html join_periodic.svg "Matching of periodic faces: base principle"

As with joining, it is recommended to filter boundary faces to process
using a selection criterion. As many periodicities may be built as desired,
as long as boundary faces are present. Once a periodicity is handled,
faces having periodic matches do not appear as boundary faces, but as
interior faces, and are thus not available anymore for other
periodicities. Translation, rotation, and mixed periodicities
can be defined.

This operation can be specified either using the GUI, as shown below, or
by user-defined functions, as shown in various
[examples](@ref cs_user_mesh_h_cs_user_mesh_periodicity).

\anchor fig_periodicities
\image html gui_mesh_periodicity.png "Defining a periodicity relation"

Modification of the mesh and geometry
-------------------------------------

### Mesh input modification or repetition
<!-- -->

Quite a few other preprocessing operations are available in the Solver:

The \ref cs_user_mesh_input function may be used
for advanced modification of the main \ref cs_mesh_t structure:

* Apply a geometric transformation or renaming groups upon reading a mesh

* reading a mesh files multiple times, in combination with application
  of individual geometric transformations and group renames for each instance

Examples with user-defined functions are provided in the
[mesh reading and modification](@ref cs_user_mesh_h_cs_user_mesh_input)
examples subsection.

### General mesh modification
<!-- -->

The \ref cs_user_mesh_modify function may be used
for advanced modification of the main \ref cs_mesh_t structure.

* extruding some faces

* inserting a refined boundary layer

* inserting a boundary between cells along selected faces

* applying a mesh refinement

The \ref cs_user_mesh_modify function may be used
for advanced modification of the main \ref cs_mesh_t structure.

Examples with user-defined functions are provided in the
[general mesh modification](@ref cs_user_mesh_h_cs_user_mesh_modifiy)
examples subsection.

\warning Caution must be exercised when modifying a mesh with
periodicity. Indeed, the periodicity parameters may not
be fully updated accordingly, meaning that the periodicity may not be valid
after mesh vertex coordinates have changed. It is particularly
true when one rescales the mesh. Most modifications should thus be done
in a separate run, before defining periodicity (joining can be done
in the same run without risk, as it is always applied first).

### Mesh smoothing utilities
<!-- -->

The principle of smoothers is to mitigate the local defects by averaging
the mesh quality. This procedure may help for calculation robustness or/and
results quality.

The \ref cs_user_mesh_smoothe user function allows to apply
smoothing operations.

An example with user-defined functions is provided in the
[mesh smoothing](@ref cs_user_mesh_h_cs_user_mesh_smoothing)
examples subsection.

\warning Caution must be exercised when using smoothing utilities
because the geometry may be modified. In order to preserve geometry,
the function \ref cs_mesh_smoother_fix_by_feature allows to
lock some boundary vertices by a feature angle criterion.
Fixing all boundary vertices ensures the geometry is preserved, but reduces
the smoothing algorithm's effectiveness.

#### Warped faces smoother
<!-- -->

The \ref cs_mesh_smoother_unwarp allows reducing face warping
in the calculation mesh.

Be aware that, in some cases, this algorithm may degrade other mesh quality
criteria.

<!-- ======================================================================= -->

\page cs_ug_mesh_select_c Face and cell mesh-defined properties and selection

Element groups {#sec_element_groups}
==============

The mesh entities may be referenced by the user during the mesh
creation. These references may then be used to mark out some mesh entities
according to the need (specification of boundary conditions, pressure
drop zones, ...). The references are generally of one of the two
following types:

* A **color** is an integer possibly associated with boundary faces and
  volume elements by the mesh generator. Depending on the tool,
  this concept may have different names, which code_saturne interprets
  as colors. Most tools allow only one color per face or element.
      - I-deas uses a color number with a default of
        7 (green) for elements, be they volume elements or boundary
          elements. Color 11 (red) is used for
          for vertices, but vertex properties are ignored by code_saturne.
      - Simail used the equivalent notions of "reference"
            for element faces, and "subdomain" for volume elements.
            By default, element faces were assigned no reference (0),
            and volume elements domain 1.
      - Gmsh uses "physical property" numbers.
      - EnSight has no similar notion, but if several parts
        are present in an EnSight 6 file, or several parts
        are present *and* vertex ids are given in an Ensight Gold file,
        the part number is interpreted as a color number by the Preprocessor.
      - The MED 2.3 model allowed integer "attributes" in addition to groups.

* Named **groups** of mesh entities may also be used with most current mesh
  generators or formats. In some cases, a given cell or face may belong
  to multiple groups (as some tools allow new groups to be defined
  by boolean operations on existing groups).
      - I-deas assigns a group number with each
        group, but by default, this number is just a counter.
        Only the group name is considered by code_saturne (so that elements
        belonging to two groups with identical names and different
        numbers are considered as belonging to the same group).
      - CGNS allows both for named boundary conditions and mesh
        sections. If present, boundary condition names are
        interpreted as group names, and groups may also be defined
        based on element section or zone names using additional
         Preprocessor [options](@ref sec_optpcs_mesh)
        (`-grp-cel` or `-grp-fac` followed by `section` or `zone`).
      - In the MED format and Salome platform, mesh groups are the main
        user-level concept to reference elements.
      - In current Gmsh versions, *physical entities* are interpreted
        as groups by code_saturne

Element selection criteria {#sec_selection_criteria}
==========================

Using the GUI or high-levels mesh element selections in user-defined
functions, *selection criteria* allow defining a selection of mesh
entities (usually cells, boundary faces, or interior faces) in a simple
and consistent manner.

Typically, a selection criteria is simply a string containing
the required group names (or color numbers for older formats),
possibly combined using boolean expressions. Simple geometric criteria are
also available.

A few examples are given below:

`ENTRY`                                     <br/>
`1 or 7`                                    <br/>
`all[]`                                     <br/>
`3.1 >= z >= -2 or not (15 or entry)`       <br/>
`range[04, 13, attribute]`                  <br/>
`sphere[0, 0, 0, 2] and (not no_group[])`   <br/>

Strings such as group names containing white-space
or having names similar to reserved operators may be protected
using "escape characters".

Note that for defining a string in Fortran, double quotes are easier to use,
as they do not conflict with Fortran's single quotes delimiting a string.
In C, the converse is true. Also, in C, to define a string such as
`\plane`, the string `\\plane` must be
used, as the first `\` character is used by the
compiler itself. Using the GUI, either notation is easy.
More complex examples of strings with protected strings are given here:

`"First entry" or Wall\ or\ sym`
`entry or \plane or "noone's output"`

The following operators and syntaxes are allowed (fully capitalized
versions of keywords are also allowed, but mixed upper-case/lower-case
versions are not):

<table>
<tr><th> Escape characters <th>
<tr><td> protect next character only <td> `\`
<tr><td> protect *string*            <td> `'string'` <br/>
                                          `"string"`
<tr><th> Basic operators <th>
<tr><td> priority        <td> `(    )`
<tr><td> not             <td>  `not    !    !=`
<tr><td> and             <td>  `and    &    &&`
<tr><td> or              <td>  `or    |    ||    ,    ;`
<tr><td> xor             <td>  `xor    ^`
<tr><th> General functions <th>
<tr><td> select all                          <td> `all[]`
<tr><td> entities having no group or color   <td> `no_group[]`
<tr><td> select a range of groups or colors  <td> `range[` *first*, *last*`]` <br/>
                                                  `range[` *first*, *last*`, group]` <br/>
                                                  `range[` *first*, *last*`, attribute]`
</table>

For the range operator, *first* and *last* values are inclusive.
For attribute (color) numbers, natural integer value ordering is used,
while for group names, alphabetical ordering is used. Note also that in
the bizarre (not recommended) case in which a mesh would contain for
example both a color number *15* and a group named "15", using
`range[15, 15, group]` or `range[15, 15, attribute]`
could be used to distinguish the two.

Geometric functions are also available. The coordinates considered are
those of the cell or face centres. Normals are of course
usable only for face selections, not cell selections.

<table>
<tr><th> Geometric functions <th>
<tr><td> face normals <td> `normal[` *x, y, z, epsilon*`]` <br/>
                           `normal[` *x, y, z*, `epsilon = ` *epsilon*`]`
<tr><td> plane, *ax + by + cz + d = 0* form
    <td> `plane[` *a, b, c, d, epsilon*`]` <br/>
         `plane[` *a, b, c, d*, ` epsilon = ` *epsilon*`]` <br/>
         `plane[` *a, b, c, d*, `inside]`  <br/>
         `plane[` *a, b, c, d*, `outside]`  <br/>
<tr><td> plane, normal + point in plane form
    <td> `plane[` *n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z, epsilon*`]` <br/>
 `plane[` *n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z,* `epsilon = ` *epsilon*`]` <br/>
 `plane[` *n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z,* `inside]` <br/>
 `plane[` *n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z,* `outside]`
<tr><td> box, extents (axis-aligned) form
    <td> `box[` *x<span style=" vertical-align:sub;">min</span>,
                 y<span style=" vertical-align:sub;">min</span>,
                 z<span style=" vertical-align:sub;">min</span>,
                 x<span style=" vertical-align:sub;">max</span>,
                 y<span style=" vertical-align:sub;">max</span>,
                 z<span style=" vertical-align:sub;">max</span>*`]`
<tr><td> box, origin + axes form
    <td> `box[` *x<span style=" vertical-align:sub;">0</span>,
                 y<span style=" vertical-align:sub;">0</span>,
                 z<span style=" vertical-align:sub;">0</span>,
                 dx<span style=" vertical-align:sub;">1</span>,
                 dy<span style=" vertical-align:sub;">1</span>,
                 dy<span style=" vertical-align:sub;">1</span>,
                 dx<span style=" vertical-align:sub;">2</span>,
                 dy<span style=" vertical-align:sub;">2</span>,
                 dz<span style=" vertical-align:sub;">2</span>,
                 dx<span style=" vertical-align:sub;">3</span>,
                 dy<span style=" vertical-align:sub;">3</span>,
                 dz<span style=" vertical-align:sub;">3</span>*`]`
<tr><td> cylinder
    <td> `cylinder[` *x<span style=" vertical-align:sub;">0</span>,
                      y<span style=" vertical-align:sub;">0</span>,
                      z<span style=" vertical-align:sub;">0</span>,
                      x<span style=" vertical-align:sub;">1</span>,
                      y<span style=" vertical-align:sub;">1</span>,
                      z<span style=" vertical-align:sub;">1</span>,
                      radius*`]`
<tr><td> sphere <td> `sphere[` *x, y, z, radius*`]`
<tr><td> inequalities <td> `>`, `<`, `>=`, `<=`
                            associated with `x`, `y`, `z`  or `X`, `Y`, `Z`
                            keywords and coordinate value <br/>
                           *x<span style=" vertical-align:sub;">min</span>* `<= x`
                           *x<span style=" vertical-align:sub;">max</span>*
                           type syntax is allowed
</table>

All selection criteria used are maintained in a list, so that
re-interpreting a criterion already encountered (such as at the previous
time step) is avoided. Lists of entities corresponding to a criteria
containing no geometric functions are also saved in a compact manner,
so re-using a previously used selection should be very fast.
For criteria containing geometric functions, the full list of
corresponding entities is not maintained, so each entity must be compared
to the criterion at each time step. Heavy use of many selection criteria
containing geometric functions may thus lead to reduced performance.

Using selection criteria in user code {#sec_fvm_selector}
------------------------

In order to use [selection criteria](@ref sec_selection_criteria) in C and
Fortran user subroutines, a collection of utility subroutines is provided.

for example:

* boundary conditions (c.f. `cs_user_boundary_conditions.f90`},
* volume initialization (c.f. \ref cs_user_initialization, ...),
* [zone](@ref sec_zones) definitions (cf. \ref cs_user_zones}),
* advanced post-processing (c.f. \ref cs_user_postprocess.c,
  \ref cs_user_extra_operations, ...),

### Selection criteria in Fortran

This section explains how to define surface or volume sections,
in the form of lists `lstelt` of `nlelt` elements
(internal faces, boundary faces or cells).
For each type of element, the user calls the appropriate Fortran subroutine:

* \ref getfbr for boundary faces
* \ref getfac for internal faces
* \ref getcel for cells.

Several examples of possible selections are given here:

*  `call getfbr("Face_1, Face_2", nlelt, lstelt)` selects
    boundary faces in groups *Face_1* or *Face_2*,

*  `call getfac("4", nlelt, lstelt)` selects internal
    faces of color *4*,

*  `call getfac("not(4)", nlelt, lstelt)` selects internal
    faces which have a different color than 4,

*  `call getfac("range[in_04, in_08]", nlelt, lstelt)` selects internal faces
    with group names between *in_04* and *in_08* (in lexicographical order),

*  `call getcel("1 or 2", nlelt, lstelt)` selects cells with colors 1 or 2,

*  `call getfbr("wall and y > 0", nlelt, lstelt)` selects boundary
    faces of group *wall* which have the coordinate *Y > 0*,

*  `call getfac("normal[1, 0, 0, 0.0001]", nlelt, lstelt)` selects
    internal faces which have a normal direction to the vector (1,0,0),

*  `call getcel("all[]", nlelt, lstelt)` selects all cells.

The user may then use a loop on the selected elements.
For instance, in the subroutine `cs_user_boundary_y_conditions` used to impose
boundary  conditions, let us consider the boundary faces of color
number 2 and which have the coordinate *X <= 0.01*
(so that `call getfbr('2 and x <= 0.01', nlelt,lstelt)`);
we can do a loop (`do ilelt = 1, nlelt`) and
obtain `ifac = lstelt(ilelt)`.

### Selection criteria in C

In C, the equivalent functions are:

* \ref cs_selector_get_b_face_list for boundary faces
* \ref cs_selector_get_i_face_list for internal faces
* \ref cs_selector_get_cell_list for cells.


More examples are available in the [User examples](@ref cs_user_examples) section.

Volume and boundary zones {#sec_zones}
=========================

Though selection criteria can and should be used directly in mesh preprocessing
operations and in some postprocessing subset extractions, for regions
with a specific meaning, it is preferable to build **zones**.

In most cases, zones can be mapped directly to mesh regions and groups,
but can be defined in a more complex manner using selection criteria.

Once built, boundary and volume zones can be used to quickly access
all matching elements, as they maintain lists of corresponding elements.

The GUI naturally builds and associates zones for boundary and volume
conditions. The \ref cs_user_zones user-defined functions (from
`cs_user_zones.c`) can be used to build such zones.
