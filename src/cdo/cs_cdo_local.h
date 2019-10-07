#ifndef __CS_CDO_LOCAL_H__
#define __CS_CDO_LOCAL_H__

/*============================================================================
 * Routines to handle low-level routines related to CDO local quantities:
 * - local matrices (stored in dense format),
 * - local quantities related to a cell.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_flag.h"
#include "cs_param_cdo.h"
#include "cs_sdm.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*! \struct cs_cell_builder_t
 *  \brief Set of local and temporary buffers useful for building the algebraic
 *  system with a cellwise process. This structure belongs to one thread.
 */
typedef struct {

  /* Specific members for the weakly enforcement of Dirichlet BCs (diffusion) */
  double     eig_ratio; /*!< ratio of the eigenvalues of the diffusion tensor */
  double     eig_max;   /*!< max. value among eigenvalues */

  /* Store the cellwise value for the diffusion, curl-curl, grad-div, time
     and reaction properties */
  cs_real_33_t  dpty_mat; /*!< Property tensor if not isotropic for diffusion */
  double        dpty_val; /*!< Property value if isotropic for diffusion */

  cs_real_33_t  cpty_mat; /*!< Property tensor if not isotropic for curl-curl */
  double        cpty_val; /*!< Property value if isotropic for curl-curl */

  double        gpty_val; /*!< Property value if isotropic for grad-div */

  double        tpty_val; /*!< Property value for time operator */

  /*! Property values for the reaction operator */
  double        rpty_vals[CS_CDO_N_MAX_REACTIONS];
  double        rpty_val; /*!< Sum of all reaction property values  */

  /* Advection-related values */
  double       *adv_fluxes;

  /* Temporary buffers (erase and updated several times during the system
     build */
  int          *ids;     /*!< local ids */
  double       *values;  /*!< local values */
  cs_real_3_t  *vectors; /*!< local 3-dimensional vectors */

  /* Structures used to build specific terms composing the algebraic system */
  cs_sdm_t     *hdg;   /*!< local hodge matrix for diffusion (may be NULL) */
  cs_sdm_t     *loc;   /*!< local square matrix of size n_cell_dofs */
  cs_sdm_t     *aux;   /*!< auxiliary local square matrix of size n_cell_dofs */

} cs_cell_builder_t;

/*! \struct cs_cell_sys_t
 *  \brief Set of arrays and local (small) dense matrices related to a mesh cell
 *  This is a key structure for building the local algebraic system.
 *  This structure belongs to one thread and only.
 */
typedef struct {

  cs_lnum_t   c_id;       /*!< cell id  */
  cs_flag_t   cell_flag;  /*!< matadata related to the cell */

  int         n_dofs;   /*!< Number of Degrees of Freedom (DoFs) in this cell */

  cs_lnum_t  *dof_ids;  /*!< DoF ids */
  cs_flag_t  *dof_flag; /*!< size = number of DoFs */

  cs_sdm_t   *mat;      /*!< cellwise view of the system matrix */
  double     *rhs;      /*!< cellwise view of the right-hand side */
  double     *source;   /*!< cellwise view of the source term array */
  double     *val_n;    /*!< values of the unkown at previous time t_n */

  /* Boundary conditions for the local system */
  short int   n_bc_faces;  /*!< Number of border faces associated to a cell */
  short int  *_f_ids;      /*!< List of face ids in the cell numbering */
  cs_lnum_t  *bf_ids;      /*!< List of face ids in the border face numbering */
  cs_flag_t  *bf_flag;     /*!< Boundary face flag; size n_bc_faces */

  /* Dirichlet BCs */
  bool        has_dirichlet; /*!< Dirichlet BCs ?*/
  double     *dir_values;    /*!< Values of the Dirichlet BCs (size = n_dofs) */

  /* Neumann BCs */
  bool        has_nhmg_neumann; /*!< Non-homogeneous Neumann BCs ? */
  double     *neu_values;       /*!< Neumann BCs values; size = n_dofs */

  /* Robin BCs */
  bool        has_robin;     /*!< Robin BCs ? */
  double     *rob_values;    /*!< Robin BCs values; size = 3*n_dofs */

  /* Sliding BCs */
  bool        has_sliding;   /*!< Sliding BCs ? */

  /* Internal enforcement of DoFs */
  bool        has_internal_enforcement;  /*!< Internal enforcement ? */
  cs_lnum_t  *intern_forced_ids;         /*!< Id in the enforcement array */

} cs_cell_sys_t;

/*! \struct cs_cell_mesh_t
 *  \brief Set of local quantities and connectivities related to a mesh cell
 *  This is a key structure for all cellwise processes. This structure belongs
 *  to one thread and only.
 *  This structure used allows one to get a better memory locality.
 *  It maps the existing global mesh and other related structures into a more
 *  compact one dedicated to a cell. Arrays are allocated to the max number of
 *  related entities (e.g. n_max_vbyc or to n_max_ebyc).
 *  The cell-wise numbering is based on the c2v, c2e and c2f connectivity.
*/

typedef struct {

  cs_eflag_t     flag;    /*!< indicate which quantities have to be computed */
  fvm_element_t  type;    /*!< type of element related to this cell */

  /* Sizes used to allocate buffers */
  short int      n_max_vbyc;
  short int      n_max_ebyc;
  short int      n_max_fbyc;

  /* Cell information */
  cs_lnum_t      c_id;    /*!< id of related cell */
  cs_real_3_t    xc;      /*!< coordinates of the cell center */
  double         vol_c;   /*!< volume of the current cell */
  double         diam_c;  /*!< diameter of the current cell */

  /* Vertex information */
  short int    n_vc;  /*!< number of vertices in a cell */
  cs_lnum_t   *v_ids; /*!< vertex ids on this rank */
  double      *xv;    /*!< local vertex coordinates (copy) */
  double      *wvc;   /*!< weight |dualvol(v) cap vol_c|/|vol_c|, size: n_vc */

  /* Edge information */
  short int    n_ec;   /*!< number of edges in a cell */
  cs_lnum_t   *e_ids;  /*!< edge ids on this rank */
  cs_quant_t  *edge;   /*!< edge quantities (xe, length and unit vector) */
  cs_nvec3_t  *dface;  /*!< dual face quantities (area and unit normal) */
  cs_real_t   *pvol_e; /*!< volume associated to an edge in the cell */

  /* Face information */
  short int    n_fc;        /*!< number of faces in a cell */
  cs_lnum_t    bface_shift; /*!< shift to get the boundary face numbering */
  cs_lnum_t   *f_ids;       /*!< face ids on this rank */
  short int   *f_sgn;       /*!< incidence number between f and c */
  double      *f_diam;      /*!< diameters of local faces */
  double      *hfc;         /*!< height of the pyramid of basis f and apex c */
  cs_quant_t  *face;        /*!< face quantities (xf, area and unit normal) */
  cs_nvec3_t  *dedge;      /*!< dual edge quantities (length and unit vector) */
  cs_real_t   *pvol_f;      /*!< volume associated to a face in the cell */

  /* Local e2v connectivity: size 2*n_ec (allocated to 2*n_max_ebyc) */
  short int   *e2v_ids;  /*!< cell-wise edge->vertices connectivity */
  short int   *e2v_sgn;  /*!< cell-wise edge->vertices orientation (-1 or +1) */

  /* Local f2v connectivity: size = 2*n_max_ebyc */
  short int   *f2v_idx;  /*!< size n_fc + 1 */
  short int   *f2v_ids;  /*!< size 2*n_max_ebyc */

  /* Local f2e connectivity: size = 2*n_max_ebyc */
  short int   *f2e_idx;  /*!< cellwise face->edges connectivity (size n_fc+1) */
  short int   *f2e_ids;  /*!< cellwise face->edges ids (size 2*n_max_ebyc) */
  short int   *f2e_sgn;  /*!< cellwise face->edges orientation (-1 or +1) */
  double      *tef;      /*!< area of the triangle of base |e| and apex xf */

  /* Local e2f connectivity: size 2*n_ec (allocated to 2*n_max_ebyc) */
  short int   *e2f_ids;  /*!< cell-wise edge -> faces connectivity */
  cs_nvec3_t  *sefc;     /*!< portion of dual faces (2 triangles by edge) */

} cs_cell_mesh_t;

/*! \struct cs_face_mesh_t
    \brief Set of local quantities and connectivities related to a mesh face
    Structure used to get a better memory locality. Map existing structure
    into a more compact one dedicated to a face.
    Arrays are allocated to n_max_vbyf (= n_max_ebyf).
    Face-wise numbering is based on the f2e connectivity.
*/

typedef struct {

  short int    n_max_vbyf; /*!< = n_max_ebyf */

  cs_lnum_t    c_id;    /*!< id of related cell */
  cs_real_3_t  xc;      /*!< pointer to the coordinates of the cell center */

  /* Face information */
  cs_lnum_t    f_id;   /*!< local mesh face id */
  short int    f_sgn;  /*!< incidence number between f and c */
  cs_quant_t   face;   /*!< face quantities (xf, area and unit normal) */
  cs_nvec3_t   dedge;  /*!< its dual edge quantities (length and unit vector) */

  /* Vertex information */
  short int    n_vf;    /*!< local number of vertices on this face */
  cs_lnum_t   *v_ids;   /*!< vertex ids (in the mesh or cellwise numbering) */
  double      *xv;      /*!< local vertex coordinates (copy) */
  double      *wvf;     /*!< weight related to each vertex */

  /* Edge information */
  short int    n_ef;    /*!< local number of edges in on this face (= n_vf) */
  cs_lnum_t   *e_ids;   /*!< edge ids (in the mesh or cellwise numbering) */
  cs_quant_t  *edge;    /*!< edge quantities (xe, length and unit vector) */
  double      *tef;     /*!< area of the triangle of base e and apex xf */

  /* Local e2v connectivity: size 2*n_ec (allocated to 2*n_max_ebyf) */
  short int   *e2v_ids;  /*!< face-wise edge -> vertices connectivity */

} cs_face_mesh_t;

/*
   A cs_face_mesh_light_t structure is close to a cs_face_mesh_t structure
   There are less members to be buildt quicker.
   Such structure is always associated to a cs_cell_mesh_t structure
*/

typedef struct {

  short int    n_max_vbyf;  /* Max number of vertices belonging to a face
                               (= n_max_ebyf) */

  cs_lnum_t    c_id;    /* id of related cell in the mesh numbering */
  short int    f;       /* id of the face in the cell mesh numbering */

  /* Vertex information */
  short int    n_vf;    /* local number of vertices on this face */
  short int   *v_ids;   /* vertex ids in the cellwise numbering */
  double      *wvf;     /* weights related to each vertex */

  /* Edge information */
  short int    n_ef;    /* local number of edges on this face (= n_vf) */
  short int   *e_ids;   /* edge ids in the cellwise numbering */
  double      *tef;     /* area of the triangle of base e and apex xf */

} cs_face_mesh_light_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

extern cs_cell_mesh_t        **cs_cdo_local_cell_meshes;
extern cs_face_mesh_t        **cs_cdo_local_face_meshes;
extern cs_face_mesh_light_t  **cs_cdo_local_face_meshes_light;

/*============================================================================
 * Staic inline function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the vertex id in the cellwise numbering associated to the
 *          given vertex id in the mesh numbering
 *
 * \param[in]       v_id    vertex id in the mesh numbering
 * \param[in]       cm      pointer to a cs_cell_mesh_t structure
 *
 * \return the vertex id in the cell numbering or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static inline short int
cs_cell_mesh_get_v(const cs_lnum_t              v_id,
                   const cs_cell_mesh_t  *const cm)
{
  if (cm == NULL)
    return -1;
  for (short int v = 0; v < cm->n_vc; v++)
    if (cm->v_ids[v] == v_id)
      return v;
  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the edge id in the cellwise numbering associated to the
 *          given edge id in the mesh numbering
 *
 * \param[in]       e_id    vertex id in the mesh numbering
 * \param[in]       cm      pointer to a cs_cell_mesh_t structure
 *
 * \return the edge id in the cell numbering or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static inline short int
cs_cell_mesh_get_e(const cs_lnum_t              e_id,
                   const cs_cell_mesh_t  *const cm)
{
  if (cm == NULL)
    return -1;
  for (short int e = 0; e < cm->n_ec; e++)
    if (cm->e_ids[e] == e_id)
      return e;
  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the face id in the cellwise numbering associated to the
 *          given face id in the mesh numbering
 *
 * \param[in]       f_id    face id in the mesh numbering
 * \param[in]       cm      pointer to a cs_cell_mesh_t structure
 *
 * \return the face id in the cell numbering or -1 if not found
 */
/*----------------------------------------------------------------------------*/

static inline short int
cs_cell_mesh_get_f(const cs_lnum_t              f_id,
                   const cs_cell_mesh_t  *const cm)
{
  if (cm == NULL)
    return -1;
  for (short int f = 0; f < cm->n_fc; f++)
    if (cm->f_ids[f] == f_id)
      return f;
  return -1;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Retrieve the list of vertices attached to a face
 *
 * \param[in]       f       face id in the cell numbering
 * \param[in]       cm      pointer to a cs_cell_mesh_t structure
 * \param[in, out]  n_vf    pointer of pointer to a cellwise view of the mesh
 * \param[in, out]  v_ids   list of vertex ids in the cell numbering
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_cell_mesh_get_f2v(short int                    f,
                     const cs_cell_mesh_t        *cm,
                     short int                   *n_vf,
                     short int                   *v_ids)
{
  /* Reset */
  *n_vf = 0;
  for (short int v = 0; v < cm->n_vc; v++) v_ids[v] = -1;

  /* Tag vertices belonging to the current face f */
  for (short int i = cm->f2e_idx[f]; i < cm->f2e_idx[f+1]; i++) {

    const int  shift_e = 2*cm->f2e_ids[i];
    v_ids[cm->e2v_ids[shift_e]] = 1;
    v_ids[cm->e2v_ids[shift_e+1]] = 1;

  } /* Loop on face edges */

  for (short int v = 0; v < cm->n_vc; v++) {
    if (v_ids[v] > 0)
      v_ids[*n_vf] = v, *n_vf += 1;
  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the next three vertices in a row from a face to edge connectivity
 *         and a edge to vertex connectivity
 *
 * \param[in]       f2e_ids     face-edge connectivity
 * \param[in]       e2v_ids     edge-vertex connectivity
 * \param[in, out]  v0          id of the first vertex
 * \param[in, out]  v1          id of the second vertex
 * \param[in, out]  v2          id of the third vertex
 */
/*----------------------------------------------------------------------------*/

static inline void
cs_cell_mesh_get_next_3_vertices(const short int   *f2e_ids,
                                 const short int   *e2v_ids,
                                 short int         *v0,
                                 short int         *v1,
                                 short int         *v2)
{
  const short int e0  = f2e_ids[0];
  const short int e1  = f2e_ids[1];
  const short int tmp = e2v_ids[2*e1];

  *v0 = e2v_ids[2*e0];
  *v1 = e2v_ids[2*e0+1];
  *v2 = ((tmp != *v0) && (tmp != *v1)) ? tmp : e2v_ids[2*e1+1];
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Is the face a boundary one ?
 *
 * \param[in]  cm     pointer to a \ref cs_cell_mesh_t structure
 * \param[in]  f      id of the face in the cellwise numbering
 *
 * \return true if this is a boundary face otherwise false
 */
/*----------------------------------------------------------------------------*/

static inline bool
cs_cell_mesh_is_boundary_face(const cs_cell_mesh_t    *cm,
                              const short int          f)
{
  if (cm->f_ids[f] - cm->bface_shift > -1)
    return true;
  else
    return false;
}

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate global structures related to a cs_cell_mesh_t and
 *         cs_face_mesh_t structures
 *
 * \param[in]   connect   pointer to a \ref cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_initialize(const cs_cdo_connect_t     *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to \ref cs_cell_mesh_t and
 *         \ref cs_face_mesh_t structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a \ref cs_cell_sys_t structure
 *
 * \param[in]   n_max_dofbyc    max number of entries
 * \param[in]   n_max_fbyc      max number of faces in a cell
 * \param[in]   n_blocks        number of blocks in a row/column
 * \param[in]   block_sizes     size of each block or NULL if n_blocks = 1
 *                              Specific treatment n_blocks = 1.
 *
 * \return a pointer to a new allocated \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_sys_t *
cs_cell_sys_create(int    n_max_dofbyc,
                   int    n_max_fbyc,
                   int    n_blocks,
                   int   *block_sizes);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Reset all members related to BC and some other ones in a
 *         \ref cs_cell_sys_t structure
 *
 * \param[in]      n_fbyc     number of faces in a cell
 * \param[in, out] csys       pointer to the \ref cs_cell_sys_t struct to reset
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_reset(int              n_fbyc,
                  cs_cell_sys_t   *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_cell_sys_t structure
 *
 * \param[in, out]  p_csys   pointer of pointer to a \ref cs_cell_sys_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_free(cs_cell_sys_t     **p_csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local system for debugging purpose
 *
 * \param[in]       msg     associated message to print
 * \param[in]       csys    pointer to a \ref cs_cell_sys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_sys_dump(const char              msg[],
                 const cs_cell_sys_t    *csys);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate \ref cs_cell_builder_t structure
 *
 * \return a pointer to the new allocated \ref cs_cell_builder_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_builder_t *
cs_cell_builder_create(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_cell_builder_t structure
 *
 * \param[in, out]  p_cb  pointer of pointer to a \ref cs_cell_builder_t struct
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_builder_free(cs_cell_builder_t     **p_cb);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate and initialize a cs_cell_mesh_t structure
 *
 * \param[in]  connect        pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cell_mesh_create(const cs_cdo_connect_t   *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_cell_mesh_t structure corresponding to mesh id
 *
 * \param[in]   mesh_id   id in the array of pointer to cs_cell_mesh_t struct.
 *
 * \return a pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cdo_local_get_cell_mesh(int    mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize to invalid values a cs_cell_mesh_t structure
 *
 * \param[in]  cm         pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_reset(cs_cell_mesh_t   *cm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Dump a cs_cell_mesh_t structure
 *
 * \param[in]    cm    pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_dump(const cs_cell_mesh_t     *cm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cell_mesh_t structure
 *
 * \param[in, out]  p_cm   pointer of pointer to a cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_free(cs_cell_mesh_t     **p_cm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_cell_mesh_t structure for a given cell id. According
 *         to the requested level, some quantities may not be defined;
 *
 * \param[in]       c_id        cell id
 * \param[in]       build_flag  indicate which members are really built
 * \param[in]       connect     pointer to a cs_cdo_connect_t structure
 * \param[in]       quant       pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  cm          pointer to a cs_cell_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_build(cs_lnum_t                    c_id,
                   cs_eflag_t                   build_flag,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_cell_mesh_t              *cm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_face_mesh_t structure
 *
 * \param[in]  n_max_vbyf    max. number of vertices for a face
 *
 * \return a pointer to a new allocated cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_face_mesh_create(short int   n_max_vbyf);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_face_mesh_t structure corresponding to mesh id
 *
 * \param[in]   mesh_id   id in the array of pointer to cs_face_mesh_t struct.
 *
 * \return a pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_cdo_local_get_face_mesh(int    mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_face_mesh_t structure
 *
 * \param[in, out]  p_fm   pointer of pointer to a cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_free(cs_face_mesh_t     **p_fm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given face/cell id.
 *
 * \param[in]       c_id      cell id
 * \param[in]       f_id      face id in the mesh structure
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_build(cs_lnum_t                    c_id,
                   cs_lnum_t                    f_id,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_face_mesh_t              *fm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_t structure for a given cell from a
 *         cs_cell_mesh_t structure.
 *         v_ids and e_ids are defined in the cell numbering given by cm
 *
 * \param[in]       cm        pointer to the reference cs_cell_mesh_t structure
 * \param[in]       f         face id in the cs_cell_mesh_t structure
 * \param[in, out]  fm        pointer to a cs_face_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_build_from_cell_mesh(const cs_cell_mesh_t    *cm,
                                  short int                f,
                                  cs_face_mesh_t          *fm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_face_mesh_light_t structure
 *
 * \param[in]  n_max_vbyf    max. number of vertices for a face
 * \param[in]  n_max_vbyc    max. number of vertices for a cell
 *
 * \return a pointer to a new allocated cs_face_mesh_light_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_light_t *
cs_face_mesh_light_create(short int   n_max_vbyf,
                          short int   n_max_vbyc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get a pointer to a cs_face_mesh_light_t structure corresponding to
 *         mesh id
 *
 * \param[in]   mesh_id   id in the cs_face_mesh_light_t array
 *
 * \return a pointer to a cs_face_mesh_light_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_light_t *
cs_cdo_local_get_face_mesh_light(int    mesh_id);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_face_mesh_light_t structure
 *
 * \param[in, out]  p_fm   pointer of pointer to a cs_face_mesh_light_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_light_free(cs_face_mesh_light_t     **p_fm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a cs_face_mesh_light_t structure starting from a
 *         cs_cell_mesh_t structure.
 *
 * \param[in]       cm     pointer to the reference cs_cell_mesh_t structure
 * \param[in]       f      face id in the cs_cell_mesh_t structure
 * \param[in, out]  fm     pointer to a cs_face_mesh_light_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_face_mesh_light_build(const cs_cell_mesh_t    *cm,
                         short int                f,
                         cs_face_mesh_light_t    *fm);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_LOCAL_H__ */
