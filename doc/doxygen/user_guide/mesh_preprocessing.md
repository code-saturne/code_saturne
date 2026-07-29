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

\page cs_ug_mesh_prepare Mesh import and preprocessing

[TOC]

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
`contains[X]`                               <br/>
`sphere[0, 0, 0, 2] and (not no_group[])`   <br/>

Strings such as group names containing white-space
or having names similar to reserved operators may be protected
using "escape characters".

In C++ (or C), to define a string such as
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
<tr><td> select a all groups or containing a substring  <td> `contains[` *substring* `]` <br/>
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
    <td> `plane[` <em>n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z, epsilon</em>`]` <br/>
 `plane[` <em>n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z,</em> `epsilon = ` *epsilon*`]` <br/>
 `plane[` <em>n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z,</em> `inside]` <br/>
 `plane[` <em>n<span style=" vertical-align:sub;">x</span>,
                   n<span style=" vertical-align:sub;">y</span>,
                   n<span style=" vertical-align:sub;">z</span>,
                   x, y, z,</em> `outside]`
<tr><td> box, extents (axis-aligned) form
    <td> `box[` <em>x<span style=" vertical-align:sub;">min</span>,
                 y<span style=" vertical-align:sub;">min</span>,
                 z<span style=" vertical-align:sub;">min</span>,
                 x<span style=" vertical-align:sub;">max</span>,
                 y<span style=" vertical-align:sub;">max</span>,
                 z<span style=" vertical-align:sub;">max</span></em>`]`
<tr><td> box, origin + axes form
    <td> `box[` <em>x<span style=" vertical-align:sub;">0</span>,
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
                 dz<span style=" vertical-align:sub;">3</span></em>`]`
<tr><td> cylinder
    <td> `cylinder[` <em>x<span style=" vertical-align:sub;">0</span>,
                      y<span style=" vertical-align:sub;">0</span>,
                      z<span style=" vertical-align:sub;">0</span>,
                      x<span style=" vertical-align:sub;">1</span>,
                      y<span style=" vertical-align:sub;">1</span>,
                      z<span style=" vertical-align:sub;">1</span>,
                      radius</em>`]`
<tr><td> sphere <td> `sphere[` *x, y, z, radius*`]`
<tr><td> inequalities <td> `>`, `<`, `>=`, `<=`
                            associated with `x`, `y`, `z`  or `X`, `Y`, `Z`
                            keywords and coordinate value <br/>
                           <em>x<span style=" vertical-align:sub;">min</span></em> `<= x`
                           <em>x<span style=" vertical-align:sub;">max</span></em>
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

In order to use [selection criteria](@ref sec_selection_criteria) in
user-defined functions, a collection of utility subroutines is provided.

for example:

* boundary conditions (c.f. `cs_user_boundary_conditions.f90`},
* volume initialization (c.f. \ref cs_user_initialization, ...),
* [zone](@ref sec_zones) definitions (cf. \ref cs_user_zones}),
* advanced post-processing (c.f. \ref cs_user_postprocess.c,
  \ref cs_user_extra_operations, ...),

### Selection criteria

This section explains how to define surface or volume sections,
in the form of lists `lstelt` of `nlelt` elements
(internal faces, boundary faces or cells).
For each type of element, the user calls the appropriate function:

* \ref cs_selector_get_b_face_list for boundary faces
* \ref cs_selector_get_i_face_list for internal faces
* \ref cs_selector_get_cell_list for cells.

Several examples of possible selections are given here:

*  `cs_selector_get_b_face_list("Face_1, Face_2", nlelt, lstelt)` selects
    boundary faces in groups *Face_1* or *Face_2*,

*  `cs_selector_get_i_face_list("4", nlelt, lstelt)` selects internal
    faces of color *4*,

*  `cs_selector_get_i_face_list("not(4)", nlelt, lstelt)` selects internal
    faces which have a different color than 4,

*  `cs_selector_get_i_face_list("range[in_04, in_08]", nlelt, lstelt)` selects internal faces
    with group names between *in_04* and *in_08* (in lexicographical order),

*  `cs_selector_get_b_face_list("contains[X]", nlelt, lstelt)` selects boundary faces
    with group names containing *X*,

*  `cs_selector_get_cell_list("1 or 2", nlelt, lstelt)` selects cells with colors 1 or 2,

*  `cs_selector_get_b_face_list("wall and y > 0", nlelt, lstelt)` selects boundary
    faces of group *wall* which have the coordinate *Y > 0*,

*  `cs_selector_get_b_face_list("normal[1, 0, 0, 0.0001]", nlelt, lstelt)` selects
    boundary faces which have a normal direction to the vector (1,0,0),

*  `cs_selector_get_cell_list("all[]", nlelt, lstelt)` selects all cells.

The user may then use a loop on the selected elements.
For instance, in the subroutine `cs_user_boundary_y_conditions` used to impose
boundary  conditions, let us consider the boundary faces of color
number 2 and which have the coordinate *X <= 0.01*
(so that `call getfbr('2 and x <= 0.01', nlelt,lstelt)`);
we can do a loop (`do ilelt = 1, nlelt`) and
obtain `ifac = lstelt(ilelt)`.

More examples are available in the [User examples](@ref cs_user_examples) section.
