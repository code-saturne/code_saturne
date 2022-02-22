/*============================================================================
 * Intersect cells with CAD object.
 *============================================================================*/

/*
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
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C and C++ library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#include <iostream>
#include <iterator>
#include <map>
#include <vector>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"
#include "cs_base.h"
#include "cs_mesh.h"
#include "cs_mesh_connect.h"
#include "fvm_triangulate.h"

/*----------------------------------------------------------------------------
 * Open CASCADE headers
 *----------------------------------------------------------------------------*/

#include <Bnd_Box.hxx>
#include <gp_Pnt.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Solid.hxx>

#include <BRepLib_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>

#include <BRepPrimAPI_MakeBox.hxx>

#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Common.hxx>

#include <TopExp.hxx>
//#include <TopExp_Explorer.hxx>

#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <Interface_Static.hxx>
#include <STEPControl_Reader.hxx>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_cad_intersect.h"

/*=============================================================================
 * Additional doxygen documentation
 *============================================================================*/

/*!
  \file cs_cad_intersect.cxx

  Intersect selected cells with a CAD volume read from a STEP format file.

  The operation may be a cut (useful when the shape represents a solid
  portion to cut from the fluid), or a boolean "common" operation, when
  the CAD shape represents the actual fluid volume.

  As this feature is not integrated by default, it may require defining
  the following values in \c cs_user_scripts.py (adapting for
  actual install paths):

  \code{.py}
  occ_include_path = "/opt/occ/7.4.0.1/include/opencascade"
  occ_lib_paths = ("/opt/occ/7.4.0.1/lib", "/gl2ps/1.4.0.1/lib")

  occ_libs = "-lTKMesh -lTKernel -lTKG2d -lTKG3d -lTKMath -lTKIGES  -lTKXSBase -lTKBin -lTKBool -lTKBO -lTKCDF -lTKBRep -lTKTopAlgo -lTKGeomAlgo -lTKGeomBase -lTKOffset -lTKPrim -lTKSTEP -lTKSTEPBase -lTKSTEPAttr -lTKHLR -lTKFeat -lTKShHealing -lTKCAF -lTKBinL -lTKLCAF -lTKFillet -lTKSTEP209"

  domain.compile_cxxflags = "-std=c++11 -I" + occ_include_path;
  domain.compile_libs = ""
  for p in occ_lib_paths:
      domain.compile_libs += "-L" + p + " -Wl,-rpath -Wl," + p + " "
      domain.compile_libs += occ_libs
  \endcode

  \todo: replace usage of BRepAlgoAPI_Cut and BRepAlgoAPI_Common
  with BOPAlgo_PaveFiller or BOPAlgo_BOP, as this seems to allow
  finer control, and the documentation seems to flag the former
  tools as obsolete.
*/

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Private global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Update extents for one cell.
 *
 * \param[in]  m                       pointer to a cs_domain_t structure
 * \param[in]  cell_id                 caller cell id
 * \param[in]  n_cell_faces            number of faces for this cell
 * \param[in]  n_cell_face_num         cell face numbers (1-based, sign orients)
 * \param[in, out]  extents            extents
 */
/*----------------------------------------------------------------------------*/

static void
_update_extents(const cs_mesh_t    *m,
                cs_lnum_t           cell_id,
                const cs_lnum_t     n_cell_faces,
                const cs_lnum_t     cell_face_num[],
                cs_real_t           extents[6])
{
  cs_lnum_t n_sub_f = 0;

  const cs_lnum_t *f_vtx_idx, *f_vtx;

  for (cs_lnum_t j = 0; j < n_cell_faces; j++) {

    cs_lnum_t orient = cell_face_num[j] > 0 ? 1 : -1;
    cs_lnum_t f_id = CS_ABS(cell_face_num[j]) - 1;
    if (f_id < m->n_b_faces) {
      f_vtx_idx = m->b_face_vtx_idx;
      f_vtx = m->b_face_vtx_lst;
    }
    else {
      f_vtx_idx = m->i_face_vtx_idx;
      f_vtx = m->i_face_vtx_lst;
      f_id -= m->n_b_faces;
    }

    cs_lnum_t f_s_id = f_vtx_idx[f_id];
    cs_lnum_t f_e_id = f_vtx_idx[f_id+1];
    cs_lnum_t n_f_v = f_e_id - f_s_id;

    for (cs_lnum_t k = 0; k < n_f_v; k++) {
      cs_lnum_t v_id = (orient > 0) ? f_vtx[f_s_id + k] : f_vtx[f_e_id -1 -k];
      const cs_real_t *coo = m->vtx_coord + v_id*3;
      for (cs_lnum_t l = 0; l < 3; l++) {
        extents[l] = CS_MIN(extents[l], coo[l]);
        extents[l+3] = CS_MAX(extents[l+3], coo[l]);
      }
    }

  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Compute extents mesh selection.
 *
 * \param[in]   m                 pointer to a cs_domain_t structure
 * \param[in]   n_cells           number of selected cells
 * \param[in]   cell_ids          ids of selected cells
 * \param[in]   n_cell_faces      number of faces for this cell
 * \param[in]   n_cell_face_num   cell face numbers (1-based, sign orients)
 * \param[out]  extents           extents
 */
/*----------------------------------------------------------------------------*/

static void
_bounding_box(const cs_mesh_t    *m,
              cs_lnum_t           n_cells,
              const cs_lnum_t     cell_ids[],
              const cs_lnum_t     cell_face_idx[],
              const cs_lnum_t     cell_face_num[],
              TopoDS_Solid       &box)
{
  cs_real_t extents[6] = {HUGE_VAL, HUGE_VAL, HUGE_VAL,
                          -HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

  for (cs_lnum_t i = 0; i < m->n_cells; i++) {

    cs_lnum_t c_id = (cell_ids != NULL) ? cell_ids[i] : i;
    cs_lnum_t s_id = cell_face_idx[i] -1;
    cs_lnum_t e_id = cell_face_idx[i+1] -1;

    _update_extents(m,
                    c_id,
                    e_id - s_id,
                    cell_face_num + s_id,
                    extents);

  }

  for (int i = 0; i < 3; i++) {
    cs_real_t d = extents[i+3] - extents[i];
    extents[i] -= d*0.01;
    extents[i] += d*0.01;
  }

  gp_Pnt pnt1(extents[0], extents[1], extents[2]);
  gp_Pnt pnt2(extents[3], extents[4], extents[5]);

  BRepPrimAPI_MakeBox mb (pnt1, pnt2);
  mb.Build();
  box = mb.Solid();
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build the shape matching one cell.
 *
 * \param[in]  m                       pointer to a cs_domain_t structure
 * \param[in]  cell_id                 caller cell id
 * \param[in]  n_cell_faces            number of faces for this cell
 * \param[in]  n_cell_face_num         cell face numbers (1-based, sign orients)
 * \param[in]  c_vtx                   cell vertices point map
 * \param[out] solid                   solid shape
 * \param[in, out]  n_max_triangles    max number of triangles
 * \param[in, out]  triangle_vertices  work array
 * \param[in, out]  triangulate_state  polygon triangulation helper
 *
 * \return 0 in case of success, error code otherwise
 */
/*----------------------------------------------------------------------------*/

static int
_build_cell(const cs_mesh_t              *m,
            cs_lnum_t                     cell_id,
            const cs_lnum_t               n_cell_faces,
            const cs_lnum_t               cell_face_num[],
            std::map<cs_lnum_t, gp_Pnt>   c_vtx,
            TopoDS_Solid                 &solid,
            std::vector<TopoDS_Shape>    &surfaces,
            cs_lnum_t                    *n_max_triangles,
            cs_lnum_t                   **triangle_vertices,
            fvm_triangulate_state_t      *triangulate_state)
{
  cs_lnum_t n_sub_f = 0;

  const cs_lnum_t *f_vtx_idx, *f_vtx;

  /* First, map cell vertices */

  cs_lnum_t n_v = 0;
  std::map<cs_lnum_t, gp_Pnt>::iterator c_vtx_it;
  BRepBuilderAPI_Sewing sew(1e-6); // default tolerance: 1e-6

  for (cs_lnum_t j = 0; j < n_cell_faces; j++) {

    cs_lnum_t orient = cell_face_num[j] > 0 ? 1 : -1;
    cs_lnum_t f_id = CS_ABS(cell_face_num[j]) - 1;
    if (f_id < m->n_b_faces) {
      f_vtx_idx = m->b_face_vtx_idx;
      f_vtx = m->b_face_vtx_lst;
    }
    else {
      f_vtx_idx = m->i_face_vtx_idx;
      f_vtx = m->i_face_vtx_lst;
      f_id -= m->n_b_faces;
    }

    cs_lnum_t f_s_id = f_vtx_idx[f_id];
    cs_lnum_t f_e_id = f_vtx_idx[f_id+1];
    cs_lnum_t n_f_v = f_e_id - f_s_id;

    // insert vertices if needed

    for (cs_lnum_t k = 0; k < n_f_v; k++) {
      cs_lnum_t v_id = (orient > 0) ? f_vtx[f_s_id + k] : f_vtx[f_e_id -1 -k];
      c_vtx_it = c_vtx.find(v_id);
      if (c_vtx_it == c_vtx.end()) {
        double x = m->vtx_coord[v_id*3];
        double y = m->vtx_coord[v_id*3+1];
        double z = m->vtx_coord[v_id*3+2];
        gp_Pnt pnt(x, y, z);
        c_vtx[v_id] = pnt;
        n_v++;
      }
    }

    // Add face

    {
      BRepLib_MakePolygon poly;
      for (cs_lnum_t k = 0; k < n_f_v; k++) {
        cs_lnum_t v_id = (orient > 0) ? f_vtx[f_s_id + k] : f_vtx[f_e_id -1 -k];
        poly.Add(c_vtx[v_id]);
      }
      poly.Close();
      TopoDS_Wire wire = poly.Wire();
      BRepBuilderAPI_MakeFace FB(wire);

      // If wire is planar, the face is built and added
      if (false && FB.IsDone()) {
        TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);
        sew.Add(face);
        surfaces.push_back(face);
        n_sub_f += 1;
      }

      // Otherwise, we need to subdivide it
      else {

        cs_lnum_t n_triangles = n_f_v - 2;
        if (n_triangles > *n_max_triangles) {
          BFT_REALLOC(*triangle_vertices, 3*n_triangles, cs_lnum_t);
          *n_max_triangles = n_triangles;
        }

        cs_lnum_t *_triangle_vertices = *triangle_vertices;

        if (n_f_v == 4)
          fvm_triangulate_quadrangle(3, // dimension
                                     0, // base
                                     m->vtx_coord,
                                     NULL,
                                     f_vtx + f_s_id,
                                     _triangle_vertices);
        else if (n_f_v > 4) {
          if (triangulate_state == NULL)
            triangulate_state = fvm_triangulate_state_create(n_f_v);
          n_triangles = fvm_triangulate_polygon(3, // dimension
                                                0, // base
                                                n_f_v,
                                                m->vtx_coord,
                                                NULL,
                                                f_vtx + f_s_id,
                                                FVM_TRIANGULATE_MESH_DEF,
                                                _triangle_vertices,
                                                triangulate_state);
          assert(n_triangles == n_f_v -2); // unless face is strongly warped.
        }

        TopoDS_Face face;
        BRepBuilderAPI_Sewing fsew; // default tolerance: 1e-6

        for (cs_lnum_t ti = 0; ti < n_triangles; ti++) {

          cs_lnum_t i0, i1, i2;
          if (orient == 1) {
            i0 = _triangle_vertices[ti*3];
            i1 = _triangle_vertices[ti*3 + 1];
            i2 = _triangle_vertices[ti*3 + 2];
          }
          else {
            i0 = _triangle_vertices[ti*3 + 2];
            i1 = _triangle_vertices[ti*3 + 1];
            i2 = _triangle_vertices[ti*3];
          }

          BRepLib_MakePolygon tria(c_vtx[i0], c_vtx[i1], c_vtx[i2]);
          tria.Close();
          TopoDS_Wire tw = tria.Wire();
          BRepBuilderAPI_MakeFace TB(tw);
          if (TB.IsDone()) {
            TopoDS_Face subface = BRepBuilderAPI_MakeFace(tw);
            fsew.Add(subface);
            n_sub_f += 1;
          }
          else
            bft_error(__FILE__, __LINE__, 0,
                      _("Error building face from vertices with coordinates\n"
                        "  [%g %g %g]\n"
                        "  [%g %g %g]\n"
                        "  [%g %g %g]"),
                      m->vtx_coord[_triangle_vertices[i0*3]],
                      m->vtx_coord[_triangle_vertices[i0*3 + 1]],
                      m->vtx_coord[_triangle_vertices[i0*3 + 2]],
                      m->vtx_coord[_triangle_vertices[i1*3]],
                      m->vtx_coord[_triangle_vertices[i1*3 + 1]],
                      m->vtx_coord[_triangle_vertices[i1*3 + 2]],
                      m->vtx_coord[_triangle_vertices[i2*3]],
                      m->vtx_coord[_triangle_vertices[i2*3 + 1]],
                      m->vtx_coord[_triangle_vertices[i2*3 + 2]]);
        }

        fsew.Perform();
        TopoDS_Shape sewedFace = fsew.SewedShape();

        sew.Add(sewedFace);
        surfaces.push_back(sewedFace);

      }
    }

  }

  sew.Perform();
  TopoDS_Shape sewedShape = sew.SewedShape();

  // Check for sewing error

  int sew_error = 0;
  if (sewedShape.IsNull())
    sew_error = 1;
  else if (! sewedShape.Closed())
    sew_error = 2;
  if (sew_error) {
    cs_gnum_t g_cell_num = cell_id+1;
    if (m->global_cell_num != NULL)
      g_cell_num = m->global_cell_num[cell_id];
    const char sew_err_str[]
      = N_("Cell %llu: error sewing shell with %d faces and %d sub-faces");
    const char not_closed_str[]
      = N_("Cell %llu: sewed shell with %d faces and %d sub-faces not closed");
    const char *err_str = sew_err_str;
    if (sew_error == 2)
      err_str = not_closed_str;
    bft_error(__FILE__, __LINE__, 0,
              _(err_str),
              (unsigned long long)g_cell_num,
              (int)n_cell_faces, (int)n_sub_f);
    return 1;
  }

  BRepBuilderAPI_MakeSolid aMkSolid;

  // Build solid

#if 0
  TopTools_IndexedMapOfShape shellMap;
  TopExp::MapShapes(sewedShape, TopAbs_SHELL, shellMap);
  for(int ishell = 1; ishell <= shellMap.Extent(); ishell++) {
    const TopoDS_Shell& shell = TopoDS::Shell(shellMap(ishell));
    aMkSolid.Add(shell);
  }
#else
  aMkSolid.Add(TopoDS::Shell(sewedShape));
#endif

  solid = aMkSolid.Solid();

  return 0;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Intersect selected cells with CAD shape.
 *
 * \param[in]       m                pointer to a mesh structure
 * \param[in]       path             path to CAD file
 * \param[in]       op               common if CAD represents fluid domain
 *                                   cut if CAD represents solid complement
 * \param[in]       n_cells          number of selected cells
 * \param[in]       cell_ids         ids of selected cells
 * \param[in, out]  cell_porosity    cell porosity
 * \param[in, out]  cell_f_center    cell fluid center, or NULL
 * \param[in, out]  i_face_porosity  interior face porosity, or NULL
 * \param[in, out]  i_face_f_center  interior face fluid center, or NULL
 * \param[in, out]  b_face_porosity  boundary face porosity, or NULL
 * \param[in, out]  b_face_f_center  boundary face fluid center, or NULL
 */
/*----------------------------------------------------------------------------*/

static void
_cad_intersect(const cs_mesh_t        *m,
               const char             *path,
               cs_cad_intersect_op_t   op,
               cs_lnum_t               n_cells,
               const cs_lnum_t         cell_ids[],
               cs_real_t               cell_porosity[],
               cs_real_t               cell_f_center[][3],
               cs_real_t               i_face_porosity[],
               cs_real_t               i_face_f_center[][3],
               cs_real_t               b_face_porosity[],
               cs_real_t               b_face_f_center[][3])
{
  if (n_cells < 1)
    return;

  const cs_lnum_t n_cells_ext = m->n_cells_with_ghosts;

  cs_lnum_t  *cell_face_idx = NULL, *cell_face_num = NULL;

  /* Initialize flags */

  bool compute_face_quantities = false;
  bool compute_face_centers = true;

  if (i_face_f_center != NULL || b_face_f_center != NULL)
    compute_face_centers = true;
  if (compute_face_centers || i_face_porosity != NULL || b_face_porosity != NULL)
    compute_face_quantities = true;

  /* Read CAD shape */

  STEPControl_Reader reader;
  IFSelect_ReturnStatus status = reader.ReadFile(path);
  Interface_Static::SetCVal ("xstep.cascade.unit", "M");
  Interface_Static::SetIVal("read.step.ideas", 1);
  Interface_Static::SetIVal("read.step.nonmanifold", 1);

  if (status == IFSelect_RetDone) {
    bft_printf(_("\n"
                 "  Loaded CAD file: %s\n"), path);
    // check UnitFlag to units from file
    TColStd_SequenceOfAsciiString anUnitLengthNames;
    TColStd_SequenceOfAsciiString anUnitAngleNames;
    TColStd_SequenceOfAsciiString anUnitSolidAngleNames;
    reader.FileUnits(anUnitLengthNames, anUnitAngleNames, anUnitSolidAngleNames);
    if (anUnitLengthNames.Length() > 0) {
      TCollection_AsciiString unit_name = anUnitLengthNames.First();
      const char *c_unit_name = unit_name.ToCString();
      bft_printf(_("    length units: %s\n"), c_unit_name);
    }

    Standard_Boolean failsonly = Standard_False;
    // reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);

    // Root transfers
    Standard_Integer aNbRoots = reader.NbRootsForTransfer();
    Standard_Integer i;

    // reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);

    for (i = 1; i <= aNbRoots; i++) {
      reader.TransferRoot(i);
    }
  }

  if (status == IFSelect_RetVoid)
    bft_error(__FILE__, __LINE__, 0,
              _("Empty result when reading file %s"),
              path);
  else if (status != IFSelect_RetDone)
    bft_error(__FILE__, __LINE__, 0,
              _("Error type %d when reading file %s"),
              (int)status, path);

  IFSelect_PrintCount mode = IFSelect_CountByItem;

  //Standard_Integer num = reader.TransferRoots();
  TopoDS_Shape cad_shape = reader.OneShape();

  {
    double x_min, y_min, z_min, x_max, y_max, z_max;
    Bnd_Box bb;
    BRepBndLib::Add(cad_shape, bb);
    bb.Get(x_min, y_min, z_min, x_max, y_max, z_max);
    GProp_GProps VProps;
    BRepGProp::VolumeProperties(cad_shape, VProps, Standard_True);
    cs_real_t cad_volume = VProps.Mass();
    bft_printf(_("    bounding box: [%g %g %g]\n"
                 "                  [%g %g %g]\n"
                 "    volume:       %g\n"),
               x_min, y_min, z_min, x_max, y_max, z_max, cad_volume);
  }

  /* Extract cell->faces index for selected subset */

  cs_lnum_t  *c_restrict_id;
  BFT_MALLOC(c_restrict_id, n_cells_ext, cs_lnum_t);
  if (cell_ids != NULL) {
    for (cs_lnum_t i = 0; i < n_cells_ext; i++)
      c_restrict_id[i] = -1;
    for (cs_lnum_t i = 0; i < n_cells; i++)
      c_restrict_id[cell_ids[i]] = i;
  }
  else {
    for (cs_lnum_t i = 0; i < n_cells_ext; i++)
      c_restrict_id[i] = i;
  }

  cs_mesh_connect_get_cell_faces(m,
                                 n_cells,
                                 c_restrict_id,
                                 &cell_face_idx,
                                 &cell_face_num);

  /* Option (optimization):
     use only part of shape within bounding box */

  if (true) {
    TopoDS_Solid   box;

    _bounding_box(m,
                  n_cells,
                  cell_ids,
                  cell_face_idx,
                  cell_face_num,
                  box);

    BRepAlgoAPI_Common c = BRepAlgoAPI_Common(box, cad_shape);
    c.Build();
    if (c.IsDone()) {
      cad_shape = c.Shape();
    }
  }

  /* Loop on identified cells */

  std::map<cs_lnum_t, gp_Pnt> c_vtx;
  std::map<cs_lnum_t, gp_Pnt>::iterator c_vtx_it;
  std::vector<TopoDS_Shape> surfaces;

  cs_lnum_t    n_max_triangles = 0;
  cs_lnum_t   *triangle_vertices = NULL;
  fvm_triangulate_state_t  *triangulate_state = NULL;

  for (cs_lnum_t i = 0; i < n_cells; i++) {

    cs_lnum_t c_id = (cell_ids != NULL) ? cell_ids[i] : i;
    cs_lnum_t s_id = cell_face_idx[i] -1;
    cs_lnum_t e_id = cell_face_idx[i+1] -1;

    TopoDS_Solid solid;

    _build_cell(m,
                c_id,
                e_id - s_id,
                cell_face_num + s_id,
                c_vtx,
                solid,
                surfaces,
                &n_max_triangles,
                &triangle_vertices,
                triangulate_state);

    GProp_GProps VProps;
    BRepGProp::VolumeProperties(solid, VProps, Standard_True);
    cs_real_t cell_volume = VProps.Mass();
    cs_lnum_t n_cell_faces = surfaces.size();
    assert(n_cell_faces == e_id - s_id);

    bool bool_ok = false;
    TopoDS_Shape  fluid_cell;

    if (op == CS_CAD_INTERSECT_CUT) {
      BRepAlgoAPI_Cut c = BRepAlgoAPI_Cut(solid, cad_shape);
      c.Build();
      if (c.IsDone()) {
        bool_ok = true;
        fluid_cell = c.Shape();
      }
    }
    else {
      BRepAlgoAPI_Common c = BRepAlgoAPI_Common(solid, cad_shape);
      c.Build();
      if (c.IsDone()) {
        bool_ok = true;
        fluid_cell = c.Shape();
      }
    }

    if (bool_ok) {
      BRepGProp::VolumeProperties(fluid_cell, VProps, Standard_True);
      cs_real_t fluid_volume = VProps.Mass();
      cell_porosity[c_id] = fluid_volume / cell_volume;
      if (cell_f_center != NULL) {
        gp_Pnt c_cen = VProps.CentreOfMass();
        cell_f_center[c_id][0] = c_cen.X();
        cell_f_center[c_id][1] = c_cen.Y();
        cell_f_center[c_id][2] = c_cen.Z();
      }
      if (compute_face_quantities) {
        for (cs_lnum_t j = 0; j < n_cell_faces; j++) {

          cs_real_t *f_porosity = NULL;
          cs_real_t *f_f_center = NULL;

          cs_lnum_t face_num = cell_face_num[s_id + j];
          cs_lnum_t f_id = CS_ABS(face_num) - 1;
          if (f_id < m->n_b_faces) {
            f_porosity = b_face_porosity + f_id;
            if (b_face_f_center != NULL)
              f_f_center = b_face_f_center[f_id];
          }
          else {
            f_id -= m->n_b_faces;
            f_porosity = i_face_porosity + f_id;
            if (i_face_f_center != NULL)
              f_f_center = i_face_f_center[f_id];
          }

          GProp_GProps SProps;
          BRepGProp::SurfaceProperties(surfaces[j], SProps, Standard_True);
          cs_real_t face_surface = SProps.Mass();
          cs_real_t face_porosity = 1;

          TopoDS_Shape  fluid_face;

          if (op == CS_CAD_INTERSECT_CUT) {
            BRepAlgoAPI_Cut c = BRepAlgoAPI_Cut(surfaces[j], cad_shape);
            c.Build();
            if (c.IsDone()) {
              bool_ok = true;
              fluid_face = c.Shape();
            }
          }
          else {
            BRepAlgoAPI_Common c = BRepAlgoAPI_Common(surfaces[j], cad_shape);
            c.Build();
            if (c.IsDone()) {
              bool_ok = true;
              fluid_face = c.Shape();
            }
          }

          if (bool_ok) {
            BRepGProp::SurfaceProperties(fluid_face, SProps, Standard_True);
            face_porosity = SProps.Mass() / face_surface;

            if (f_porosity != NULL)
              f_porosity[0] = face_porosity;

            if (f_f_center != NULL) {
              gp_Pnt f_cen = SProps.CentreOfMass();
              f_f_center[0] = f_cen.X();
              f_f_center[1] = f_cen.Y();
              f_f_center[2] = f_cen.Z();
            }
          }
        }
      }
    }
    else {
      cell_porosity[c_id] = -1;
    }

    // clear local vertices map and surfaces vector
    c_vtx.clear();
    surfaces.clear();

  }

  // Cleanup

  BFT_FREE(triangle_vertices);
  if (triangulate_state != NULL)
    triangulate_state = fvm_triangulate_state_destroy(triangulate_state);

  BFT_FREE(c_restrict_id);
}

/*----------------------------------------------------------------------------*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Intersect selected cells with CAD shape.
 *
 * \param[in]       m                pointer to a mesh structure
 * \param[in]       path             path to CAD file
 * \param[in]       op               common if CAD represents fluid domain
 *                                   cut if CAD represents solid complement
 * \param[in]       n_cells          number of selected cells
 * \param[in]       cell_ids         ids of selected cells
 * \param[in, out]  cell_porosity    cell porosity
 * \param[in, out]  cell_f_center    cell fluid center, or NULL
 * \param[in, out]  i_face_porosity  interior face porosity, or NULL
 * \param[in, out]  i_face_f_center  interior face fluid center, or NULL
 * \param[in, out]  b_face_porosity  boundary face porosity, or NULL
 * \param[in, out]  b_face_f_center  boundary face fluid center, or NULL
 */
/*----------------------------------------------------------------------------*/

void
cs_cad_intersect(const cs_mesh_t        *m,
                 const char             *path,
                 cs_cad_intersect_op_t   op,
                 cs_lnum_t               n_cells,
                 const cs_lnum_t         cell_ids[],
                 cs_real_t               cell_porosity[],
                 cs_real_t               cell_f_center[][3],
                 cs_real_t               i_face_porosity[],
                 cs_real_t               i_face_f_center[][3],
                 cs_real_t               b_face_porosity[],
                 cs_real_t               b_face_f_center[][3])
{
  _cad_intersect(m,
                 path,
                 op,
                 n_cells,
                 cell_ids,
                 cell_porosity,
                 cell_f_center,
                 i_face_porosity,
                 i_face_f_center,
                 b_face_porosity,
                 b_face_f_center);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
