#ifndef __CS_CDO_LOCAL_H__
#define __CS_CDO_LOCAL_H__

/*============================================================================
 * Routines to handle low-level routines related to CDO local quantities:
 * - local matrices (stored in dense format),
 * - local quantities related to a cell.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2016 EDF S.A.

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

#include "cs_cdo.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_toolbox.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_CDO_LOCAL_V   (1 << 0) //  1: local information related to vertices
#define CS_CDO_LOCAL_E   (1 << 1) //  2: local information related to edges
#define CS_CDO_LOCAL_F   (1 << 2) //  4: local information related to faces
#define CS_CDO_LOCAL_EV  (1 << 3) //  8: cell-wise edge --> vertices connect.
#define CS_CDO_LOCAL_FE  (1 << 4) // 16: cell-wise face --> edges connect.
#define CS_CDO_LOCAL_EF  (1 << 5) // 32: cell-wise edge --> faces connect.

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure used to store a local system (cell-wise for instance) */
typedef struct {

  cs_locmat_t  *mat;      // cell-wise view of the system matrix
  double       *rhs;      // cell-wise view of the right-hand side
  double       *dir_bc;   // Dirichlet values for boundary degrees of freedom

} cs_cdo_locsys_t;

/* Structure used to get a better memory locality. Map existing structure
   into a more compact one dedicated to a cell.
   Arrays are allocated to n_max_vbyc or to n_max_ebyc.
   Cell-wise numbering is based on the c2e and c2v connectivity.
*/

typedef struct {

  short int    n_max_vbyc;
  short int    n_max_ebyc;
  short int    n_max_fbyc;

  cs_flag_t    flag;   // indicate which quantities are defined
  cs_lnum_t    c_id;   // id of related cell
  cs_real_t   *xc;     // pointer to the coordinates of the cell center
  double       vol_c;  // volume of the current cell

  /* Vertex information */
  short int    n_vc;   // local number of vertices in a cell
  cs_lnum_t   *v_ids;  // vertex ids on this rank
  double      *xv;     // local vertex coordinates (copy)
  short int   *vtag;   // link between mesh and cell-wise numbering (-1 not set)
  double      *wvc;    // weight |vol_dc(v) cap vol_c|/|vol_c for each cell vtx

  /* Edge information */
  short int    n_ec;   // local number of edges in a cell
  cs_lnum_t   *e_ids;  // edge ids on this rank
  short int   *etag;   // link between mesh and cell-wise numbering (-1 not set)
  cs_quant_t  *edge;   // local edge quantities (xe, length and unit vector)
  cs_nvec3_t  *dface;  // local dual face quantities (area and unit normal)

  /* Face information */
  short int    n_fc;   // local number of faces in a cell
  cs_lnum_t   *f_ids;  // face ids on this rank
  short int   *f_sgn;  // incidence number between f and c
  cs_quant_t  *face;   // local face quantities (xf, area and unit normal)
  cs_nvec3_t  *dedge;  // local dual edge quantities (length and unit vector)

  /* Local e2v connectivity: size 2*n_ec (allocated to 2*n_max_ebyc) */
  short int   *e2v_ids;  // cell-wise edge -> vertices connectivity
  short int   *e2v_sgn;  // cell-wise edge -> vertices orientation (-1 or +1)

  /* Local f2e connectivity: size = n_fc*n_max_ebyf */
  short int   *f2e_idx;  // size n_fc + 1
  short int   *f2e_ids;  // size f2e_idx[n_fc]

  /* Local e2f connectivity: size 2*n_ec (allocated to 2*n_max_ebyc) */
  short int   *e2f_ids;  // cell-wise edge -> faces connectivity
  short int   *e2f_sgn;  // cell-wise edge -> faces orientation (-1 or +1)

} cs_cell_mesh_t;

/* Structure used to get a better memory locality. Map existing structure
   into a more compact one dedicated to a face.
   Arrays are allocated to n_max_vbyf (= n_max_ebyf).
   Face-wise numbering is based on the f2e connectivity.
*/

typedef struct {

  short int    n_max_vbyf; // = n_max_ebyf

  cs_lnum_t    c_id;    // id of related cell
  cs_real_t   *xc;      // pointer to the coordinates of the cell center

  /* Face information */
  cs_lnum_t    f_id;    // local mesh face id
  short int    f_sgn;   // incidence number between f and c
  cs_quant_t   face;    // local face quantities (xf, area and unit normal)
  cs_nvec3_t   dedge;   // local dual edge quantities (length and unit vector)

  /* Vertex information */
  short int    n_vf;    // local number of vertices on this face
  cs_lnum_t   *v_ids;   // vertex ids on this rank or in the cellwise numbering
  double      *xv;      // local vertex coordinates (copy)
  double      *wvf;     // weight related to each vertex

  /* Edge information */
  short int    n_ef;    // local number of edges in on this face (= n_vf)
  cs_lnum_t   *e_ids;   // edge ids on this rank or in the cellwise numbering
  cs_quant_t  *edge;    // local edge quantities (xe, length and unit vector)
  double      *tef;     // area of the triangle of base e and apex xf

  /* Local e2v connectivity: size 2*n_ec (allocated to 2*n_max_ebyf) */
  short int   *e2v_ids;  // face-wise edge -> vertices connectivity

} cs_face_mesh_t;

/*============================================================================
 *  Global variables
 *============================================================================*/

extern cs_cell_mesh_t  **cs_cdo_local_cell_meshes       ;
extern cs_face_mesh_t  **cs_cdo_local_face_meshes;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_cdo_locsys_t structure
 *
 * \param[in]   n_max_ent    max number of entries
 *
 * \return a pointer to a new allocated cs_cdo_locsys_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cdo_locsys_t *
cs_cdo_locsys_create(int    n_max_ent);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a cs_cdo_locsys_t structure
 *
 * \param[in, out]  p_ls   pointer of pointer to a cs_cdo_locsys_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_locsys_free(cs_cdo_locsys_t     **p_ls);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate global structures related to a cs_cell_mesh_t and
 *         cs_face_mesh_t structures
 *
 * \param[in]   connect   pointer to a cs_cdo_connect_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_initialize(const cs_cdo_connect_t     *connect);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free global structures related to cs_cell_mesh_t and cs_face_mesh_t
 *         structures
 */
/*----------------------------------------------------------------------------*/

void
cs_cdo_local_finalize(void);

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
 * \brief  Allocate a cs_cell_mesh_t structure
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_cell_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_cell_mesh_t *
cs_cell_mesh_create(const cs_cdo_connect_t     *connect);

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
 * \param[in]       c_id      cell id
 * \param[in]       level     indicate which members are really defined
 * \param[in]       connect   pointer to a cs_cdo_connect_t structure
 * \param[in]       quant     pointer to a cs_cdo_quantities_t structure
 * \param[in, out]  cm        pointer to a cs_cell_mesh_t structure to set
 */
/*----------------------------------------------------------------------------*/

void
cs_cell_mesh_build(cs_lnum_t                    c_id,
                   cs_flag_t                    level,
                   const cs_cdo_connect_t      *connect,
                   const cs_cdo_quantities_t   *quant,
                   cs_cell_mesh_t              *cm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Allocate a cs_face_mesh_t structure
 *
 * \param[in]      connect   pointer to a cs_cdo_connect_t structure
 *
 * \return a pointer to a new allocated cs_face_mesh_t structure
 */
/*----------------------------------------------------------------------------*/

cs_face_mesh_t *
cs_face_mesh_create(const cs_cdo_connect_t     *connect);

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

END_C_DECLS

#endif /* __CS_CDO_LOCAL_H__ */
