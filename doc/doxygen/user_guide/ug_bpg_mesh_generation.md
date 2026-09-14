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

\page ug_bpg_mesh_generation Best Practice Guide for mesh generation

[TOC]

# BPG - Mesh Generation

## Introduction

The capability of *code_saturne* to import any mesh (any type of element,
non-conformal meshes, etc.) has been developed to help users during the mesh
generation stage.

However, the best suited meshes for the numerical schemes employed in
*code_saturne* are conformal meshes made of:

- **Cubes** with edges aligned with the streamlines or, if that is not
possible, bricks (orthogonal hexahedra) with a small aspect ratio and edges
aligned with the streamlines.
- **Equilateral tetrahedra** for flows without any privileged direction.

For RANS computations, at least, a good mesh consisting of tetrahedra (a fine
mesh of equilateral tetrahedra) is a better choice than a mesh of hexahedra
containing non-orthogonal cells.

\image html bpg_mesh/Fig1_bpg_mesh_generation.png "Figure 1: An orthogonal mesh" width=20%

## General Advice

When these simple rules cannot be applied to the whole mesh, the following
criteria should be respected as much as possible. They are practical
recommendations based on experience rather than strict mathematical
requirements.

### Extrusion at inlet and outlet boundaries for tetrahedral meshes

With tetrahedral meshes, it is necessary to apply an extrusion of at least 1
cell at inlet and outlet boundaries to prevent non physical behavior near
these boundaries. This will generate layers of orthogonal prisms, add them to
the rest of the mesh and preserve all group of boundary faces and cells.
The thickness of the cells should be similar to the size of tetrahedra near the
boundary face. The extrusion step can be defined either in the GUI or in
the `cs_user_mesh_modifiy` user source file.

### Mesh alignment

**Align the mesh with the streamlines or with the stratifications** (i.e. with
the expected isolines of the important quantities): in particular, it is
advisable to create 1 to 5 layers of hexahedra or of orthogonal prisms at the
walls (i.e. prisms with edges aligned with the normal to the wall).

### Mesh criterias

**Enforce the following geometrical criteria**, as much as possible (remembering
that difficulties may not necessarily appear on low quality cells if the
quantities are uniform in these regions):

#### Warping

Warping angles of 5° or more should be avoided (Figure 2).

- Dividing the warped faces into triangles (to eliminate the warping) does not
necessarily improve the precision (even if the local truncation error is a
priori reduced) and may deteriorate other quality criteria. However, dividing
the warped faces into triangles is compulsory for Lagrangian computations (to
avoid the risk of losing particles).

- Creating non-conformal meshes may lead to warped faces.

\image html bpg_mesh/Fig2_bpg_mesh_generation.png "Figure 2: a warped face,
colored in red" width=20%

#### Aspect ratio

The aspect ratio L/h (Figure 3 illustrates in 2D the 3D criterion which is
the ratio of the larger to the smaller characteristic length of the cell):
- The optimal value is 1.
- Maximal value is 10, and up to 100 for cells aligned with the streamlines (or
even more, depending on the physical phenomena considered).

\image html bpg_mesh/Fig3_bpg_mesh_generation.png "Figure 3: aspect ratio L/h" width=20%

#### Centering deviation

The centering deviation for a given face of a cell is the distance between
the center of the face F and the point O, defined as the intersection of the
plane of the face with the line defined by the centres of the neighboring cells
(Figure 4, Figure 5):
- optimal value: 0
- maximal value: try to keep the point O within the face
- creating non-conformal meshes and refining the mesh may modify the value of
the centering deviation (Figure 6)

\image html bpg_mesh/Fig4_bpg_mesh_generation.png "Figure 4: remarkable points" width=30%

\image html bpg_mesh/Fig5_bpg_mesh_generation.png "Figure 5: centering deviation (intersection outside of the face)" width=20%

\image html bpg_mesh/Fig6_bpg_mesh_generation.png "Figure 6: examples of modification of the centering deviation: improvement due to refinement (left), degradation due to refinement (centre), degradation on a non-conformal mesh (right)"

#### Non-orthogonality angle

The non-orthogonality angle is the angle between the normal to a face and, for internal
faces, the line joining the centre of the neighboring cells or, for boundary
faces, the line joining the centre of the face to the centre of the neighboring
cell (Figure 7):
- optimal value: 0
- usual values: lower than 15° for hexahedra and lower than 45° for tetrahedra
- maximal value: 80°; a large non-orthogonality angle is more easily handled by
the solver when it is located away from regions where the gradients of the
variables are large (for example non-orthogonality should be avoided near
walls). Creating non-conformal meshes may generate very large angles of
non-orthogonality. The maximum value that is supported depends on the physical
phenomena that are considered.

\image html bpg_mesh/Fig7_bpg_mesh_generation.png "Figure 7: non-orthogonality angle (left) and effect of non-conformal meshes (right)"

#### Maximal weighting

The maximal weighting (|distance FJ’/ distance I’J’|) for internal faces, with
I’ and J’ standing for the projection of I and J (centres of the neighboring
cells) on the line normal to their common face and containing the centre F of
the face (Figure 4); the accuracy my decrease if the weighting is too large:
- optimal value: 0.5

#### Growth rate

The growth rate is the ratio of the size of two successive cells in a given
direction; one should in particular avoid situations where this ratio oscillates
about 1 irregular mesh where cells are alternatively long and short as on Figure
8):
- optimal value: 1
- maximal value: 1.5 to 2

\image html bpg_mesh/Fig8_bpg_mesh_generation.png "Figure 8: example of an irregular mesh to avoid" width=40%

### Mesh joining

Multiple meshes may be joined and assembled into a single mesh, and faces lying
on a common surface may be joined, whether their interfaces are conformal
(i.e. their vertices are coincident) or not.

#### Non-conformal meshes and conformal joining

When non-conformal faces are joined,  they are subdivided into conforming
faces shared with adjacent cells, based on the subdivision of their edges along
their intersections, and addition of vertices if necessary. Very small edges
are avoided by merging their vertices, based on a geometric tolerance.

This operation thus modifies the underlying structure so as to obtain closed
polyhedra with conformal faces(Figure 9). Thus initially non-conforming
hexahedral meshes effectively are transformed inti conforming polyhedral meshes.

\image html bpg_mesh/Fig9_bpg_mesh_generation.png "Figure 9: examples of non-conformal meshes"

\image html bpg_mesh/Fig10_bpg_mesh_generation.png "Figure 10: examples of non-conformal meshes – coarsening ratio 2 cells / 3 cells (left), 1 cell / 2 cells (right)" width=50%

** recommendations **

- Avoid conformal joining for LES (as LES is sensitive to local jumps in cell
  size and modifications of the discretization which can lead to non-physical
  energy in the flow field).
- If several successive layers of conformal joining is used to coarsen, it is
advised to use a coarsening ratio of 2 cells / 3 cells rather than 1 cell / 2
cells (Figure 10). This is particularly important if the turbulence level is low
and if the main direction of the flow is normal to the conformal joining surface
(if the turbulence level is high, the mixing may help to eliminate the
perturbations that could appear because of the checkerboard structure of the
mesh).
- For the other cases of conformal joining, the coarsening ratio should be kept
of the order of 1 cell / 2 cell and below 1 cell / 5 cells (if needed, use several
steps/layers or redesign the mesh); otherwise, it may be necessary to
use an upwind convection scheme to stabilize the computation, bearing in mind
that the accuracy of the results maybe be reduced and LES turbulence modeling
will most probably be unreliable).
- If too many non-conformal joining interfaces are needed, it may be wise to
contemplate an unstructured mesh consisting of tetrahedra.
  * Joining non-conforming meshes may locally degrade their quality along
    joining interfaces, while using tetrahedra will lead to lower quality
    cells overall, but perhaps with a better "minimum cell quality".
- Where possible choose plane interfaces in which the joining meshes map onto
each other exactly. The meshes that should be joined should rest as exactly as
possible on the same geometrical surface, with as little overlapping or gaps as
possible.
- Place conformal joining interfaces that may produce non-orthogonal mesh cells
as far as possible from the regions of interest (and away from regions with
large gradients of the variables, in particular).
- Identify the sets of faces to be joined together. Indeed, it is possible to
let code_saturne decide which elements should be joined (on the basis of
geometry-based criteria). However, it is advised to use group names to
explicitly identify the faces that must be joined so as to speed up the process
and improve the robustness of the joining operation proper.
- As much as possible, differentiate between the faces associated with the each
conformal joining interface. In doing so avoid the use of the group name already
associated with the boundary conditions (this makes the completion of the
joining process much easier to check).

#### Mesh joining checks

- Check the logs and visualize the joined mesh to ensure that there is no face
or portion of face that was not joined correctly. Such faces will remain as
boundary faces, so will appear as such. Also, if all faces of a given group
should be joined, only interior faces of that group should appear after a
successful joining operation, so the log should not list boundary faces in
that group.
- In the presence of residual un-joined boundary faces or sub-faces, it is
possible to set a slip boundary condition on such residual portions of faces.
However, it remains necessary to visualize them to make sure that they may not
create any perturbation in the flow. For example, Figure 12 illustrates how
portions of joined faces may remain and produce small steps (two meshes of the
same circular-section pipe are considered; their refinement is different; they are
joined along a cross section perpendicular to the pipe 3 axis). More clearly,
one may think of a mesh approximating the circular section by an octagon and a
coarser mesh for which the refinement only allows to approximate the circular
section by an hexagon: the joining of these two sections creates residual
portions of faces that introduce irregularities of the surface if the code does
not manage to detect that the vertices shall be displaced locally to avoid this
artefact.

\image html bpg_mesh/Fig12_bpg_mesh_generation.png "Figure 12: example of conformal joining potentially leading to residual boundary faces" width=50%

### Predefined mesh patterns

- Use O-mesh for circular sections (the central pattern can be a square, a
hexagon, an octagon...) and around obstacles.
- Use pre-existing patterns for T-junctions and mixing grids (see existing
studies).

\image html bpg_mesh/Fig11_bpg_mesh_generation.png "Figure 11: O-meshes" width=40%

### Specify the cell size from physical considerations

- Evaluate the size of the cells at the wall from the boundary layer thickness
and from the constraints imposed by the selected turbulence model (see the
dedicated section).
- Evaluate the size of the cells in the core of the domain from the experience
acquired through previous studies on similar geometries and from the size of the
structures of the flow that shall be resolved.
- Use at least 5 cells between two facing walls (with less than 5 cells, the
fluid will go through, but the modeling will be too coarse to account for
anything but for the mass conservation; if this situation is not local
(associated with a singularity of the geometry), one should envisage to change
the wall boundary into a slip boundary and to add a head loss source term
accounting for the wall friction.

## Specific modelling

Some models require specific treatment:

### Second-Order RANS models (Reynolds Stress Models)

The structures that are resolved by second-Order RANS models are generally finer
than those captured by first-order models (k-epsilon or k-omega) and advected
further away during a longer periods of time. Indeed, second-order models are
generally less diffusive, they account for secondary motion (corners, flow
structures downstream a bend...), and the underlying system of equations has a
“more convective” nature than the system of equations that is associated with
first-order models (the latter is essentially based on the equilibrium between
production and dissipation source terms). Hence, the level of refinement that is
required to obtain a converged result with a first-order model is generally
lower than with that required for a second-order model. Moreover, for a given
coarse mesh, the results may be much worse with a second-order model than with a
first-order model.

### Large Eddy Simulation

For Large Eddy Simulation (LES), the mesh must be selected to resolve the large
anisotropic structures, the model being designed to represent the smaller
isotropic ones. To evaluate the cell size, one may - for want of something
better - carry out a preliminary RANS calculation to evaluate the size of the
turbulent structures:
- The integral scale \f$L_T = \alpha k^{3/2}/\varepsilon\f$, with α
approximately ranging from 0.1 to 0.3, provides the size of the large
structures.
- The Kolmogorov length-scale \f$\eta = (\nu^3/\varepsilon)^{1/4}\f$
provides the size of the smaller structures (the smaller vortices are
immediately dissipated by the fluid viscosity).

LES is theoretically applicable when \f$L_T\f$ >> η (i.e. when the “inertial
zone” of the turbulent spectrum is established: usually for a sufficiently
developed turbulence and with a turbulent Reynolds number
\f$Re_t = (L_T/\eta)\f$ large enough, typically superior to 1000). One should
select a mesh size of the order of (or smaller than) \f$L_T/10\f$ (with cells
smaller than the Kolmogorov length-scale, the simulation effectively becomes a
direct simulation). It is not always easy to honor this criterion, since the
integral scale may be very small (for example: in a channel flow, the integral
scale is proportional to the distance to the wall...).
