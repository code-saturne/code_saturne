/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2009 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_POST_H__
#define __CS_POST_H__

/*============================================================================
 * Post-processing management
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * FVM library headers
 *----------------------------------------------------------------------------*/

#include <fvm_nodal.h>
#include <fvm_writer.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Local type definitions
 *============================================================================*/

/* Datatype enumeration */

typedef enum {
  CS_POST_TYPE_cs_int_t,
  CS_POST_TYPE_cs_real_t,
  CS_POST_TYPE_int,
  CS_POST_TYPE_float,
  CS_POST_TYPE_double
} cs_post_type_t;

/* Function pointer associated with a specific post-processing variables
   output: such functions are registered using the cs_post_add_time_dep_var(),
   and all registered functions are automatically called by PSTVAR. */

typedef void
(cs_post_time_dep_var_t) (int        instance_id,
                          int        nt_cur_abs,
                          cs_real_t  t_cur_abs);

/*=============================================================================
 * Global variables
 *============================================================================*/

/*============================================================================
 * Public Fortran function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a writer based on Fortran data; this object is based on a choice
 * of a case, directory, and format, as well as indicator for associated
 * meshe's time dependency, and the default output frequency for associated
 * variables.
 *
 * Fortran Interface: use PSTCWR (see cs_post_util.F)
 *
 * SUBROUTINE PSTCW1 (NUMGEP, NOMCAS, NOMREP, NOMFMT, OPTFMT,
 * *****************
 *                    LNMCAS, LNMFMT, LNMREP, LOPFMT,
 *                    INDMOD, NTCHR)
 *
 * INTEGER          NUMWRI      : <-- : Number of writer to create (< 0 for
 *                              :     : standard writer, > 0 for user writer)
 * CHARACTER        NOMCAS      : <-- : Name of associated case
 * CHARACTER        NOMREP      : <-- : Name of associated directory
 * INTEGER          NOMFMT      : <-- : Name of associated format
 * INTEGER          OPTFMT      : <-- : Additional format options
 * INTEGER          LNMCAS      : <-- : Case name length
 * INTEGER          LNMREP      : <-- : Directory name length
 * INTEGER          LNMFMT      : <-- : Format name length
 * INTEGER          LOPFMT      : <-- : Format options string length
 * INTEGER          INDMOD      : <-- : 0 if fixed, 1 if deformable,
 *                              :     : 2 if topology changes
 * INTEGER          NTCHR       : <-- : Default output frequency
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcw1, PSTCW1)
(
 const cs_int_t  *numwri,
 const char      *nomcas,
 const char      *nomrep,
 const char      *nomfmt,
 const char      *optfmt,
 const cs_int_t  *lnmcas,
 const cs_int_t  *lnmrep,
 const cs_int_t  *lnmfmt,
 const cs_int_t  *lopfmt,
 const cs_int_t  *indmod,
 const cs_int_t  *ntchr
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Create a post-processing mesh; lists of cells or faces to extract are
 * sorted upon exit, whether they were sorted upon calling or not.
 *
 * The list of associated cells is only necessary if the number of cells
 * to extract is strictly greater than 0 and less than the number of cells
 * of the computational mesh.
 *
 * Lists of faces are ignored if the number of extracted cells is nonzero;
 * otherwise, if the number of boundary faces to extract is equal to the
 * number of boundary faces in the computational mesh, and the number of
 * interior faces to extract is zero, than we extrac by default the boundary
 * mesh, and the list of associated boundary faces is thus not necessary.
 *
 * Fortran interface: use PSTCMA (see cs_post_util.F)
 *
 * SUBROUTINE PSTCM1 (NUMMAI, NOMMAI, LNMMAI,
 * *****************
 *                    NBRCEL, NBRFAC, NBRFBR, LSTCEL, LSTFAC, LSTFBR)
 *
 * INTEGER          NUMMAI      : <-- : Number of output mesh to create
 *                              :     : (< 0 for standard mesh,
 *                              :     : > 0 for user mesh)
 * CHARACTER        NOMMAI      : <-- : Name of associated output mesh
 * INTEGER          LNMMAI      : <-- : Mesh name length
 * INTEGER          NBRCEL      : <-- : Number of associated cells
 * INTEGER          NBRFAC      : <-- : Number of associated interior faces
 * INTEGER          NBRFBR      : <-- : Nulber of associated boundary faces
 * INTEGER          LSTCEL      : <-- : List of associated cells
 * INTEGER          LSTFAC      : <-- : List of associated interior faces
 * INTEGER          LSTFBR      : <-- : List of associated boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcm1, PSTCM1)
(
 const cs_int_t  *nummai,
 const char      *nommai,
 const cs_int_t  *lnmmai,
 const cs_int_t  *nbrcel,
 const cs_int_t  *nbrfac,
 const cs_int_t  *nbrfbr,
       cs_int_t   lstcel[],
       cs_int_t   lstfac[],
       cs_int_t   lstfbr[]
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*----------------------------------------------------------------------------
 * Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTEDG (NUMMAI, NUMREF)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the edges mesh to create
 * INTEGER          NUMREF      : <-- : Number of the existing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstedg, PSTEDG)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numref
);

/*----------------------------------------------------------------------------
 * Create an alias to a post-processing mesh.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTALM (NUMMAI, NUMREF)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the alias to create
 * INTEGER          NUMREF      : <-- : Number of the associated output mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstalm, PSTALM)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numref
);

/*----------------------------------------------------------------------------
 * Associate a writer to a post-processing mesh.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTASS (NUMMAI, NUMWRI)
 * *****************
 *
 * INTEGER          NUMMAI      : <-- : Number of the associated output mesh
 * INTEGER          NUMWRI      : <-- : Number of the writer to associate
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstass, PSTASS)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numwri
);

/*----------------------------------------------------------------------------
 * Update the "active" or "inactive" flag for writers based on the current
 * time step and their default output frequency.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTNTC (NTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : <-- : Current time step number
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t  *ntcabs
);

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTNTC (NUMWRI, INDACT)
 * *****************
 *
 * INTEGER          NUMWRI      : <-- : Writer number, or 0 for all writers
 * INTEGER          INDACT      : <-- : 0 to deactivate, 1 to activate
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstact, PSTACT)
(
 const cs_int_t  *numwri,
 const cs_int_t  *indact
);

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTEMA (NTCABS, TTCABS)
 * *****************
 *
 * INTEGER          NTCABS      : <-- : Current time step number
 * DOUBLE PRECISION TTCABS      : <-- : Current physical time
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstema, PSTEMA)
(
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs
);

/*----------------------------------------------------------------------------
 * Loop on post-processing meshes to output variables
 *
 * Fortran interface:
 *
 * SUBROUTINE PSTVAR (IDBIA0, IDBRA0,
 * *****************
 *                    NDIM,   NTCABS, NCELET, NCEL,   NFAC,   NFABOR,
 *                    NFML,   NPRFML, NNOD,   LNDFAC, LNDFBR, NCELBR,
 *                    NVAR,   NSCAL,  NPHAS,  NVLSTA, NVISBR,
 *                    NIDEVE, NRDEVE, NITUSE, NRTUSE,
 *                    IFACEL, IFABOR, IFMFBR, IFMCEL, IPRFML,
 *                    IPNFAC, NODFAC, IPNFBR, NODFBR,
 *                    IDEVEL, ITUSER, IA,
 *                    TTCABS, XYZCEN, SURFAC, SURFBO, CDGFAC, CDGFBO,
 *                    XYZNOD, VOLUME,
 *                    DT,     RTPA,   RTP,    PROPCE, PROPFA, PROPFB,
 *                    COEFA,  COEFB,
 *                    STATCE, STATIV, STATFB,
 *                    RDEVEL, RTUSER, RA)
 *
 * INTEGER          IDBIA0      : <-- : Number of first free position in IA
 * INTEGER          IDBRA0      : <-- : Number of first free position in RA
 * INTEGER          NDIM        : <-- : Spatial dimension
 * INTEGER          NTCABS      : --> : Current time step number
 * INTEGER          NCELET      : <-- : Number of extended (real + ghost) cells
 * INTEGER          NFAC        : <-- : Number of interior faces
 * INTEGER          NFABOR      : <-- : Number of boundary faces
 * INTEGER          NFML        : <-- : Number of families (group classes)
 * INTEGER          NPRFML      : <-- : Number of family properties
 * INTEGER          NNOD        : <-- : Number of vertices
 * INTEGER          LNDFAC      : <-- : Size of nodfac
 * INTEGER          LNDFBR      : <-- : Size of nodfbr
 * INTEGER          NCELBR      : <-- : Number of cells on boundary
 * INTEGER          NVAR        : <-- : Number of variables
 * INTEGER          NSCAL       : <-- : Number of scalars
 * INTEGER          NPHAS       : <-- : Number of phases
 * INTEGER          NVLSTA      : <-- : Number of statistical variables (lagr)
 * INTEGER          NVISBR      : <-- : Number of boundary stat. variables (lagr)
 * INTEGER          NIDEVE      : <-- : Size of IDEVEL integer array
 * INTEGER          NRDEVE      : <-- : Size of RDEVEL floating-point array
 * INTEGER          NITUSE      : <-- : Size of ITUSER integer array
 * INTEGER          NRTUSE      : <-- : Size of RTUSER floating-point array
 * INTEGER          IFACEL      : <-- : Interior faces -> cells connectivity
 * INTEGER          IFABOR      : <-- : Boundary faces -> cell connectivity
 * INTEGER          IFMFBR      : <-- : Boundary face families
 * INTEGER          IFMCEL      : <-- : Cell families
 * INTEGER          IPRFML      : <-- : List of family properties
 * INTEGER          IPNFAC      : <-- : Interior faces -> vertices connect. idx.
 * INTEGER          NODFAC      : <-- : Interior faces -> vertices connectivity
 * INTEGER          IPNFBR      : <-- : Boundary faces -> vertices connect. idx.
 * INTEGER          NODFBR      : <-- : Boundary faces -> vertices connectivity
 * INTEGER          IDEVEL      : <-- : IDEVEL integer array
 * INTEGER          ITUSER      : <-- : ITUSER integer array
 * INTEGER          IA          : <-- : IA integer array
 * DOUBLE PRECISION TTCABS      : <-- : Current physical time
 * DOUBLE PRECISION XYZCEN      : <-- : Points associated with cell centers
 * DOUBLE PRECISION SURFAC      : <-- : Interior face surface vectors
 * DOUBLE PRECISION SURFBO      : <-- : Boundary face surface vectors
 * DOUBLE PRECISION CDGFAC      : <-- : Interior face centers
 * DOUBLE PRECISION CDGFBO      : <-- : Boundary face vectors
 * DOUBLE PRECISION XYZNOD      : <-- : Vertex coordinates (optional)
 * DOUBLE PRECISION VOLUME      : <-- : Cell volumes
 * DOUBLE PRECISION DT          : <-- : Local time step
 * DOUBLE PRECISION RTPA        : <-- : Cell variables at previous time step
 * DOUBLE PRECISION RTP         : <-- : Cell variables
 * DOUBLE PRECISION PROPCE      : <-- : Cell physical properties
 * DOUBLE PRECISION PROPFA      : <-- : Interior face physical properties
 * DOUBLE PRECISION PROPFB      : <-- : Boundary face physical properties
 * DOUBLE PRECISION COEFA       : <-- : Boundary conditions array
 * DOUBLE PRECISION COEFB       : <-- : Boundary conditions array
 * DOUBLE PRECISION STATCE      : <-- : Cell statistics (Lagrangian)
 * DOUBLE PRECISION STATIV      : <-- : Cell variance statistics (Lagrangian)
 * DOUBLE PRECISION STATFB      : <-- : Boundary face statistics (Lagrangian)
 * DOUBLE PRECISION RDEVEL      : <-- : RDEVEL floating-point array
 * DOUBLE PRECISION RTUSER      : <-- : RTUSER floating-point array
 * DOUBLE PRECISION RA          : <-- : RA floating-point array
 *
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *idbia0,
 const cs_int_t   *idbra0,
 const cs_int_t   *ndim,
 const cs_int_t   *ntcabs,
 const cs_int_t   *ncelet,
 const cs_int_t   *ncel,
 const cs_int_t   *nfac,
 const cs_int_t   *nfabor,
 const cs_int_t   *nfml,
 const cs_int_t   *nprfml,
 const cs_int_t   *nnod,
 const cs_int_t   *lndfac,
 const cs_int_t   *lndfbr,
 const cs_int_t   *ncelbr,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nphas,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr,
 const cs_int_t   *nideve,
 const cs_int_t   *nrdeve,
 const cs_int_t   *nituse,
 const cs_int_t   *nrtuse,
 const cs_int_t    ifacel[],
 const cs_int_t    ifabor[],
 const cs_int_t    ifmfbr[],
 const cs_int_t    ifmcel[],
 const cs_int_t    iprfml[],
 const cs_int_t    ipnfac[],
 const cs_int_t    nodfac[],
 const cs_int_t    ipnfbr[],
 const cs_int_t    nodfbr[],
 const cs_int_t    idevel[],
       cs_int_t    ituser[],
       cs_int_t    ia[],
 const cs_real_t  *ttcabs,
 const cs_real_t   xyzcen[],
 const cs_real_t   surfac[],
 const cs_real_t   surfbo[],
 const cs_real_t   cdgfac[],
 const cs_real_t   cdgfbo[],
 const cs_real_t   xyznod[],
 const cs_real_t   volume[],
 const cs_real_t   dt[],
 const cs_real_t   rtpa[],
 const cs_real_t   rtp[],
 const cs_real_t   propce[],
 const cs_real_t   propfa[],
 const cs_real_t   propfb[],
 const cs_real_t   coefa[],
 const cs_real_t   coefb[],
 const cs_real_t   statce[],
 const cs_real_t   stativ[],
 const cs_real_t   statfb[],
 const cs_real_t   rdevel[],
       cs_real_t   rtuser[],
       cs_real_t   ra[]
);

/*----------------------------------------------------------------------------
 * Post-processing output of a variable defined on cells or faces of a mesh
 * using associated writers.
 *
 * Fortran interface; use PSTEVA (see cs_post_util.F)
 *
 * SUBROUTINE PSTEV1 (NUMMAI, NOMVAR, LNMVAR, IDIMT,  IENTLA, IVARPR,
 * *****************
 *                    NTCABS, TTCABS, VARCEL, VARFAC, VARFBR)
 *
 * INTEGER          NUMMAI      : <-- : Number of associated output mesh
 * CHARACTER        NOMVAR      : <-- : Name of associated variable
 * INTEGER          LNMVAR      : <-- : Variable name length
 * INTEGER          IDIMT       : <-- : 1 for scalar, 3 for vector
 * INTEGER          IENTLA      : <-- : If a vector, 1 for interlaced values
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 otherwise (x1, x2, ...xn, y1, y2, ...)
 * INTEGER          IVARPR      : <-- : 1 if variable is defined on "parent"
 *                              :     : mesh, 2 if defined on output mesh
 * INTEGER          NTCABS      : <-- : Current time step number
 * DOUBLE PRECISION TTCABS      : <-- : Current physical time
 * DOUBLE PRECISION VARCEL(*)   : <-- : Cell values
 * DOUBLE PRECISION VARFAC(*)   : <-- : Interior face values
 * DOUBLE PRECISION VARFBO(*)   : <-- : Boundary face values
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstev1, PSTEV1)
(
 const cs_int_t   *nummai,
 const char       *nomvar,
 const cs_int_t   *lnmvar,
 const cs_int_t   *idimt,
 const cs_int_t   *ientla,
 const cs_int_t   *ivarpr,
 const cs_int_t   *ntcabs,
 const cs_real_t  *ttcabs,
 const cs_real_t   varcel[],
 const cs_real_t   varfac[],
 const cs_real_t   varfbr[]
 CS_ARGF_SUPP_CHAINE              /*     (possible 'length' arguments added
                                         by many Fortran compilers) */
);

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a writer; this objects manages a case's name, directory, and format,
 * as well as associated meshe's time dependency, and the default output
 * frequency for associated variables.
 *
 * parameters:
 *   writer_id <-- id of writer to create (< 0 reserved, > 0 for user)
 *   case_name <-- associated case name
 *   dir_name  <-- associated directory name
 *   fmt_name  <-- associated format name
 *   fmt_opts  <-- associated format options
 *   mod_flag  <-- 0 if fixed, 1 if deformable, 2 if topolygy changes,
 *                 +10 add a displacement field
 *   frequency <-- default output frequency
 *----------------------------------------------------------------------------*/

void
cs_post_add_writer(int          writer_id,
                   const char  *case_name,
                   const char  *dir_name,
                   const char  *fmt_name,
                   const char  *fmt_opts,
                   cs_int_t     mod_flag,
                   cs_int_t     frequency);

/*----------------------------------------------------------------------------
 * Create a post-processing mesh; lists of cells or faces to extract are
 * sorted upon exit, whether they were sorted upon calling or not.
 *
 * The list of associated cells is only necessary if the number of cells
 * to extract is strictly greater than 0 and less than the number of cells
 * of the computational mesh.
 *
 * Lists of faces are ignored if the number of extracted cells is nonzero;
 * otherwise, if the number of boundary faces to extract is equal to the
 * number of boundary faces in the computational mesh, and the number of
 * interior faces to extract is zero, than we extrac by default the boundary
 * mesh, and the list of associated boundary faces is thus not necessary.
 *
 * parameters:
 *   mesh_id     <-- id of mesh to create (< 0 reserved, > 0 for user)
 *   mesh_name   <-- associated mesh name
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-> list of associated cells
 *   i_face_list <-> list of associated interior faces
 *   b_face_list <-> list of associated boundary faces
 *----------------------------------------------------------------------------*/

void
cs_post_add_mesh(int          mesh_id,
                 const char  *mesh_name,
                 cs_int_t     n_cells,
                 cs_int_t     n_i_faces,
                 cs_int_t     n_b_faces,
                 cs_int_t     cell_list[],
                 cs_int_t     i_face_list[],
                 cs_int_t     b_face_list[]);

/*----------------------------------------------------------------------------
 * Create a post-processing mesh associated with an existing exportable mesh
 * representation.
 *
 * If the exportable mesh is not intended to be used elsewhere, one can choose
 * to transfer its property to the post-processing mesh, which will then
 * manage its lifecycle based on its own requirements.
 *
 * If the exportable mesh must still be shared, one must be careful to
 * maintain consistency between this mesh and the post-processing output.
 *
 * parameters:
 *   mesh_id   <-- number of mesh to create (< 0 reserved, > 0 for user)
 *   exp_mesh  <-- mesh in exportable representation (i.e. fvm_nodal_t)
 *   transfer  <-- if true, ownership of exp_mesh is transferred to the
 *                 post-processing mesh
 *----------------------------------------------------------------------------*/

void
cs_post_add_existing_mesh(int           mesh_id,
                          fvm_nodal_t  *exp_mesh,
                          cs_bool_t     transfer);

/*----------------------------------------------------------------------------
 * Create a mesh based upon the extraction of edges from an existing mesh.
 *
 * The newly created edges have no link to their parent elements, so
 * no variable referencing parent elements may be output to this mesh,
 * whose main use is to visualize "true" face edges when polygonal faces
 * are subdivided by the writer. In this way, even highly non-convex
 * faces may be visualized correctly if their edges are overlaid on
 * the surface mesh with subdivided polygons.
 *
 * parameters:
 *   edges_id <-- id of edges mesh to create (< 0 reserved, > 0 for user)
 *   base_id  <-- id of existing mesh (< 0 reserved, > 0 for user)
 *----------------------------------------------------------------------------*/

void
cs_post_add_mesh_edges(int  edges_id,
                       int  base_id);

/*----------------------------------------------------------------------------
 * Create an alias to a post-processing mesh.
 *
 * An alias allows association of an extra identifier (number) to an
 * existing post-processing mesh, and thus to associate different writers
 * than those associated with the existing mesh. For example, this allows
 * outputting a set of main variables every n1 time steps with one writer,
 * and outputting a specific set of variables every n2 time time steps to
 * another post-processing set using another writer, without the overhead
 * that would be incurred by duplication of the post-processing mesh.
 *
 * An alias is thus treated in all points like its associated mesh;
 * if the definition of either one is modified, that of the other is
 * modified also.
 *
 * It is forbidden to associate an alias to another alias (as there is no
 * identified use for this, and it would make consistency checking more
 * difficult), but multiple aliases may be associated with a given mesh.
 *
 * parameters:
 *   alias_id <-- id of alias to create (< 0 reserved, > 0 for user)
 *   mesh_id  <-- id of associated mesh
 *----------------------------------------------------------------------------*/

void
cs_post_alias_mesh(int  alias_id,
                   int  mesh_id);

/*----------------------------------------------------------------------------
 * Check for the existence of a writer of the given id.
 *
 * parameters:
 *   writer_id <-- writer id to check
 *
 * returns:
 *   true if writer with this id exists, false otherwise
 *----------------------------------------------------------------------------*/

cs_bool_t
cs_post_writer_exists(int  writer_id);

/*----------------------------------------------------------------------------
 * Return a pointer to the FVM library writer associated to a writer_id.
 *
 * parameters:
 *   writer_id <-- associated writer id
 *
 * Returns:
 *  a pointer to a fvm_writer_t structure
 *----------------------------------------------------------------------------*/

fvm_writer_t *
cs_post_get_writer(cs_int_t  writer_id);

/*----------------------------------------------------------------------------
 * Check for the existence of a post-processing mesh of the given id.
 *
 * parameters:
 *   mesh_id <-- mesh id to check
 *
 * returns:
 *   true if mesh with this id exists, false otherwise
 *----------------------------------------------------------------------------*/

cs_bool_t
cs_post_mesh_exists(int  mesh_id);

/*----------------------------------------------------------------------------
 * Modify an existing post-processing mesh.
 *
 * The lists of cells or faces are redefined, for example to update an
 * extracted mesh based in "interesting" zones.
 *
 * It is not necessary to use this function if a mesh is simply deformed.
 *
 * parameters:
 *   mesh_id     <-- id of mesh to modify (< 0 reserved, > 0 for user)
 *   n_cells     <-- number of associated cells
 *   n_i_faces   <-- number of associated interior faces
 *   n_b_faces   <-- number of associated boundary faces
 *   cell_list   <-> list of associated cells
 *   i_face_list <-> list of associated interior faces
 *   b_face_list <-> list of associated boundary faces
 *
 *----------------------------------------------------------------------------*/

void
cs_post_modify_mesh(int       mesh_id,
                    cs_int_t  n_cells,
                    cs_int_t  n_i_faces,
                    cs_int_t  n_b_faces,
                    cs_int_t  cell_list[],
                    cs_int_t  i_face_list[],
                    cs_int_t  b_face_list[]);

/*----------------------------------------------------------------------------
 * Return the next "reservable" (i.e. non-user) mesh id available.
 *
 * Returns:
 *   the smallest negative integer present, -1
 *----------------------------------------------------------------------------*/

int
cs_post_get_free_mesh_id(void);

/*----------------------------------------------------------------------------
 * Associate a writer with a post-processing mesh.
 *
 * parameters:
 *   mesh_id   <-- id of associated mesh
 *   writer_id <-- id of associated writer
 *----------------------------------------------------------------------------*/

void
cs_post_associate(int  mesh_id,
                  int  writer_id);

/*----------------------------------------------------------------------------
 * Update "active" or "inactive" flag of writers whose output frequency
 * is a divisor of the current time step number.
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *----------------------------------------------------------------------------*/

void
cs_post_activate_if_default(int  nt_cur_abs);

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * parameters:
 *   writer_id <-- writer id, or 0 for all writers
 *   activate  <-- 0 to deactivate, 1 to activate
 *----------------------------------------------------------------------------*/

void
cs_post_activate_writer(int  writer_id,
                        int  activate);

/*----------------------------------------------------------------------------
 * Output post-processing meshes using associated writers.
 *
 * parameters:
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *----------------------------------------------------------------------------*/

void
cs_post_write_meshes(int     nt_cur_abs,
                     double  t_cur_abs);

/*----------------------------------------------------------------------------
 * Output a variable defined at cells or faces of a post-processing mesh
 * using associated writers.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   var_name    <-- name of variable to output
 *   var_dim     <-- 1 for scalar, 3 for vector
 *   interlace   <-- if a vector, true for interlaced values, false otherwise
 *   use_parent  <-- true if values are defined on "parent" mesh,
 *                   false if values are defined on post-processing mesh
 *   var_type    <-- variable's data type
 *   nt_cur_abs  <-- current time step number
 *   t_cur_abs   <-- current physical time
 *   cel_vals    <-- cell values
 *   i_face_vals <-- interior face values
 *   b_face_vals <-- boundary face values
 *----------------------------------------------------------------------------*/

void
cs_post_write_var(int              mesh_id,
                  const char      *var_name,
                  cs_int_t         var_dim,
                  cs_bool_t        interlace,
                  cs_bool_t        use_parent,
                  cs_post_type_t   var_type,
                  cs_int_t         nt_cur_abs,
                  cs_real_t        t_cur_abs,
                  const void      *cel_vals,
                  const void      *i_face_vals,
                  const void      *b_face_vals);

/*----------------------------------------------------------------------------
 * Output a variable defined at vertices of a post-processing mesh using
 * associated writers.
 *
 * parameters:
 *   mesh_id    <-- id of associated mesh
 *   var_name   <-- name of variable to output
 *   var_dim    <-- 1 for scalar, 3 for vector
 *   interlace  <-- if a vector, true for interlaced values, false otherwise
 *   use_parent <-- true if values are defined on "parent" mesh,
 *                  false if values are defined on post-processing mesh
 *   var_type   <-- variable's data type
 *   nt_cur_abs <-- current time step number
 *   t_cur_abs  <-- current physical time
 *   vtx_vals   <-- vertex values
 *----------------------------------------------------------------------------*/

void
cs_post_write_vertex_var(int              mesh_id,
                         const char      *var_name,
                         cs_int_t         var_dim,
                         cs_bool_t        interlace,
                         cs_bool_t        use_parent,
                         cs_post_type_t   var_type,
                         cs_int_t         nt_cur_abs,
                         cs_real_t        t_cur_abs,
                         const void      *vtx_vals);

/*----------------------------------------------------------------------------
 * Update references to parent mesh of post-processing meshes in case of
 * computational mesh cell renumbering.
 *
 * This function may be called only once, after possible renumbering of cells,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * parameters:
 *   init_cell_num <-- initial cell numbering (1 to n, new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_cells(const cs_int_t  init_cell_num[]);

/*----------------------------------------------------------------------------
 * Update references to parent mesh of post-processing meshes in case of
 * computational mesh interior and/or boundary faces renumbering.
 *
 * This function may be called only once, after possible renumbering of faces,
 * to update existing post-processing meshes. Post-processing meshes defined
 * after renumbering will automatically be based upon the new numbering,
 * so this function will not need to be called again.
 *
 * parameters:
 *   init_i_face_num <-- initial interior numbering (1 to n, new -> old)
 *   init_b_face_num <-- initial boundary numbering (1 to n, new -> old)
 *----------------------------------------------------------------------------*/

void
cs_post_renum_faces(const cs_int_t  init_i_face_num[],
                    const cs_int_t  init_b_face_num[]);

/*----------------------------------------------------------------------------
 * Destroy all structures associated with post-processing
 *----------------------------------------------------------------------------*/

void
cs_post_finalize(void);

/*----------------------------------------------------------------------------
 * Initialize main post-processing writer
 *----------------------------------------------------------------------------*/

void
cs_post_init_main_writer(void);

/*----------------------------------------------------------------------------
 * Initialize main post-processing meshes
 *----------------------------------------------------------------------------*/

void
cs_post_init_main_meshes(void);

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * troubleshoot errors.
 *----------------------------------------------------------------------------*/

void
cs_post_init_error_writer(void);

/*----------------------------------------------------------------------------
 * Register a processing of a time-dependent variable to the call to PSTVAR.
 *
 * The instance identifier associated with the function allows registering
 * the same function several times, with a diferent identifier allowing the
 * function to select a specific operation or data.
 *
 * parameters:
 *   function    <-- function to register
 *   instance_id <-- instance id associated with this registration
 *----------------------------------------------------------------------------*/

void
cs_post_add_time_dep_var(cs_post_time_dep_var_t  *function,
                         int                      instance_id);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_POST_H__ */
