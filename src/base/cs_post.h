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
 * mesh's time dependency, and the default output frequency for associated
 * variables.
 *
 * Fortran Interface: use pstcwr (see cs_post_f2c.f90)
 *
 * subroutine pstcw1 (numgep, nomcas, nomrep, nomfmt, optfmt,
 * *****************
 *                    lnmcas, lnmfmt, lnmrep, lopfmt, indmod, ntchr)
 *
 * integer          numwri      : <-- : number of writer to create (< 0 for
 *                              :     : standard writer, > 0 for user writer)
 * character        nomcas      : <-- : name of associated case
 * character        nomrep      : <-- : name of associated directory
 * integer          nomfmt      : <-- : name of associated format
 * integer          optfmt      : <-- : additional format options
 * integer          lnmcas      : <-- : case name length
 * integer          lnmrep      : <-- : directory name length
 * integer          lnmfmt      : <-- : format name length
 * integer          lopfmt      : <-- : format options string length
 * integer          indmod      : <-- : 0 if fixed, 1 if deformable,
 *                              :     : 2 if topology changes
 * integer          ntchr       : <-- : default output frequency in time-steps
 * double precision frchr       : <-- : default output frequency in seconds
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
 const cs_int_t  *ntchr,
 const cs_real_t *frchr
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
 * Fortran interface: use pstcma (see cs_post_f2c.f90)
 *
 * subroutine pstcm1 (nummai, nommai, lnmmai,
 * *****************
 *                    nbrcel, nbrfac, nbrfbr, lstcel, lstfac, lstfbr)
 *
 * integer          nummai      : <-- : number of output mesh to create
 *                              :     : (< 0 for standard mesh,
 *                              :     : > 0 for user mesh)
 * character        nommai      : <-- : name of associated output mesh
 * integer          lnmmai      : <-- : mesh name length
 * integer          indgrp      : <-- : 1 to add group information, or O
 * integer          nbrcel      : <-- : number of associated cells
 * integer          nbrfac      : <-- : number of associated interior faces
 * integer          nbrfbr      : <-- : nulber of associated boundary faces
 * integer          lstcel      : <-- : list of associated cells
 * integer          lstfac      : <-- : list of associated interior faces
 * integer          lstfbr      : <-- : list of associated boundary faces
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcm1, PSTCM1)
(
 const cs_int_t  *nummai,
 const char      *nommai,
 const cs_int_t  *indgrp,
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
 * subroutine pstedg (nummai, numref)
 * *****************
 *
 * integer          nummai      : <-- : number of the edges mesh to create
 * integer          numref      : <-- : number of the existing mesh
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstedg, PSTEDG)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numref
);

/*----------------------------------------------------------------------------
 * Assign a category to a post-processing mesh.
 *
 * By default, each mesh is assigned a category id identical to its id.
 * The automatic variables output associated with the main volume and
 * boundary meshes will also be applied to meshes of the same categories
 * (i.e. -1 and -2 respectively, whether meshes -1 and -2 are actually
 * defined or not), so setting a user mesh's category to one of these
 * values will automatically provide the same automatic variable output to
 * the user mesh.
 *
 * Fortran interface:
 *
 * subroutine pstcat (nummai, numwri)
 * *****************
 *
 * integer          nummai      : <-- : number of the alias to create
 * integer          numcat      : <-- : number of the assigned category
 *                                      (-1: as volume, -2: as boundary)
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstcat, PSTCAT)
(
 const cs_int_t  *nummai,
 const cs_int_t  *numcat
);

/*----------------------------------------------------------------------------
 * Create an alias to a post-processing mesh.
 *
 * Fortran interface:
 *
 * subroutine pstalm (nummai, numref)
 * *****************
 *
 * integer          nummai      : <-- : number of the alias to create
 * integer          numref      : <-- : number of the associated output mesh
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
 * subroutine pstass (nummai, numwri)
 * *****************
 *
 * integer          nummai      : <-- : number of the associated output mesh
 * integer          numwri      : <-- : number of the writer to associate
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
 * subroutine pstntc (ntcabs, ttcabs)
 * *****************
 *
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : absolute time at the current time step
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstntc, PSTNTC)
(
 const cs_int_t  *ntcabs,
 const cs_real_t *ttcabs
);

/*----------------------------------------------------------------------------
 * Force the "active" or "inactive" flag for a specific writer or for all
 * writers for the current time step.
 *
 * Fortran interface:
 *
 * subroutine pstact (numwri, indact)
 * *****************
 *
 * integer          numwri      : <-- : writer number, or 0 for all writers
 * integer          indact      : <-- : 0 to deactivate, 1 to activate
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
 * subroutine pstema (ntcabs, ttcabs)
 * *****************
 *
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : current physical time
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
 * subroutine pstvar (idbia0, idbra0,
 * *****************
 *                    ntcabs,
 *                    nvar,   nscal,  nphas,  nvlsta, nvisbr,
 *                    nideve, nrdeve, nituse, nrtuse,
 *                    idevel, ituser, ia,
 *                    ttcabs,
 *                    dt,     rtpa,   rtp,    propce, propfa, propfb,
 *                    coefa,  coefb,
 *                    statce, stativ, statfb,
 *                    rdevel, rtuser, ra)
 *
 * integer          idbia0      : <-- : number of first free position in ia
 * integer          idbra0      : <-- : number of first free position in ra
 * integer          ntcabs      : --> : current time step number
 * integer          nvar        : <-- : number of variables
 * integer          nscal       : <-- : number of scalars
 * integer          nphas       : <-- : number of phases
 * integer          nvlsta      : <-- : number of statistical variables (lagr)
 * integer          nvisbr      : <-- : number of boundary stat. variables (lagr)
 * integer          nideve      : <-- : size of idevel integer array
 * integer          nrdeve      : <-- : size of rdevel floating-point array
 * integer          nituse      : <-- : size of ituser integer array
 * integer          nrtuse      : <-- : size of rtuser floating-point array
 * integer          idevel      : <-- : idevel integer array
 * integer          ituser      : <-- : ituser integer array
 * integer          ia          : <-- : ia integer array
 * double precision ttcabs      : <-- : current physical time
 * double precision dt          : <-- : local time step
 * double precision rtpa        : <-- : cell variables at previous time step
 * double precision rtp         : <-- : cell variables
 * double precision propce      : <-- : cell physical properties
 * double precision propfa      : <-- : interior face physical properties
 * double precision propfb      : <-- : boundary face physical properties
 * double precision coefa       : <-- : boundary conditions array
 * double precision coefb       : <-- : boundary conditions array
 * double precision statce      : <-- : cell statistics (lagrangian)
 * double precision stativ      : <-- : cell variance statistics (lagrangian)
 * double precision statfb      : <-- : boundary face statistics (lagrangian)
 * double precision rdevel      : <-- : rdevel floating-point array
 * double precision rtuser      : <-- : rtuser floating-point array
 * double precision ra          : <-- : ra floating-point array
 *----------------------------------------------------------------------------*/

void CS_PROCF (pstvar, PSTVAR)
(
 const cs_int_t   *idbia0,
 const cs_int_t   *idbra0,
 const cs_int_t   *ntcabs,
 const cs_int_t   *nvar,
 const cs_int_t   *nscal,
 const cs_int_t   *nphas,
 const cs_int_t   *nvlsta,
 const cs_int_t   *nvisbr,
 const cs_int_t   *nideve,
 const cs_int_t   *nrdeve,
 const cs_int_t   *nituse,
 const cs_int_t   *nrtuse,
 const cs_int_t    idevel[],
       cs_int_t    ituser[],
       cs_int_t    ia[],
 const cs_real_t  *ttcabs,
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
 * fortran interface; use psteva (see cs_post_f2c.f90)
 *
 * subroutine pstev1 (nummai, nomvar, lnmvar, idimt,  ientla, ivarpr,
 * *****************
 *                    ntcabs, ttcabs, varcel, varfac, varfbr)
 *
 * integer          nummai      : <-- : number of associated output mesh
 * character        nomvar      : <-- : name of associated variable
 * integer          lnmvar      : <-- : variable name length
 * integer          idimt       : <-- : 1 for scalar, 3 for vector
 * integer          ientla      : <-- : if a vector, 1 for interlaced values
 *                              :     : (x1, y1, z1, x2, y2, ..., yn, zn),
 *                              :     : 0 otherwise (x1, x2, ...xn, y1, y2, ...)
 * integer          ivarpr      : <-- : 1 if variable is defined on "parent"
 *                              :     : mesh, 2 if defined on output mesh
 * integer          ntcabs      : <-- : current time step number
 * double precision ttcabs      : <-- : current physical time
 * double precision varcel(*)   : <-- : cell values
 * double precision varfac(*)   : <-- : interior face values
 * double precision varfbo(*)   : <-- : boundary face values
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
 * as well as associated mesh's time dependency, and the default output
 * frequency for associated variables.
 *
 * parameters:
 *   writer_id   <-- number of writer to create (< 0 reserved, > 0 for user)
 *   case_name   <-- associated case name
 *   dir_name    <-- associated directory name
 *   fmt_name    <-- associated format name
 *   fmt_opts    <-- associated format options
 *   mod_flag    <-- 0 if fixed, 1 if deformable, 2 if topolygy changes,
 *                   +10 add a displacement field
 *   frequency_n <-- default output frequency in time-steps
 *   frequency_t <-- default output frequency in seconds
 *----------------------------------------------------------------------------*/

void
cs_post_add_writer(int          writer_id,
                   const char  *case_name,
                   const char  *dir_name,
                   const char  *fmt_name,
                   const char  *fmt_opts,
                   cs_int_t     mod_flag,
                   cs_int_t     frequency_n,
                   cs_real_t    frequency_s);

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
 *   mesh_id      <-- number of mesh to create (< 0 reserved, > 0 for user)
 *   mesh_name    <-- associated mesh name
 *   add_families <-- add family information if possible
 *   n_cells      <-- number of associated cells
 *   n_i_faces    <-- number of associated interior faces
 *   n_b_faces    <-- number of associated boundary faces
 *   cell_list    <-- list of associated cells
 *   i_face_list  <-- list of associated interior faces
 *   b_face_list  <-- list of associated boundary faces
 *----------------------------------------------------------------------------*/

void
cs_post_add_mesh(int          mesh_id,
                 const char  *mesh_name,
                 cs_bool_t    add_families,
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
 * The mesh in exportable dimension may be of a lower dimension than
 * its parent mesh, if it has been projected. In this case, a
 * dim_shift value of 1 indicates that parent cells are mapped to
 * exportable faces, and faces to edges, while a dim_shift value of 2
 * would indicate that parent cells are mapped to edges.
 * This is important when variables values are exported.
 *
 * parameters:
 *   mesh_id   <-- number of mesh to create (< 0 reserved, > 0 for user)
 *   exp_mesh  <-- mesh in exportable representation (i.e. fvm_nodal_t)
 *   dim_shift <-- nonzero if exp_mesh has been projected
 *   transfer  <-- if true, ownership of exp_mesh is transferred to the
 *                 post-processing mesh
 *----------------------------------------------------------------------------*/

void
cs_post_add_existing_mesh(int           mesh_id,
                          fvm_nodal_t  *exp_mesh,
                          int           dim_shift,
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
 * Remove a post-processing mesh.
 *
 * No further post-processing output will be allowed on this mesh,
 * so the associated structures may be freed.
 *
 * A post-processing mesh that has been associated with a time-varying
 * writer or that is referenced by an alias may not be removed.
 *
 * parameters:
 *   mesh_id <-- id of mesh to remove
 *----------------------------------------------------------------------------*/

void
cs_post_free_mesh(int  mesh_id);

/*----------------------------------------------------------------------------
 * Assign a category to a post-processing mesh.
 *
 * By default, each mesh is assigned a category id identical to its id.
 * The automatic variables output associated with the main volume and
 * boundary meshes will also be applied to meshes of the same categories
 * (i.e. -1 and -2 respectively, whether meshes -1 and -2 are actually
 * defined or not), so setting a user mesh's category to one of these
 * values will automatically provide the same automatic variable output to
 * the user mesh.
 *
 * parameters:
 *   mesh_id     <-- id of associated mesh
 *   category_id <-- id of mesh category (-1: as volume, -2: as boundary)
 *----------------------------------------------------------------------------*/

void
cs_post_set_mesh_category(int  mesh_id,
                          int  category_id);

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
 * Return the next "reservable" (i.e. non-user) writer id available.
 *
 * Returns:
 *   the smallest negative integer present, -1
 *----------------------------------------------------------------------------*/

int
cs_post_get_free_writer_id(void);

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
 *   t_cur_abs  <-- absolute time at the current time step
 *----------------------------------------------------------------------------*/

void
cs_post_activate_if_default(int     nt_cur_abs,
                            double  t_cur_abs);

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
 *
 * The check_flag variable is a mask, used for additionnal post-processing:
 *
 *  - If (check_flag & 1), volume submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 *  - If (check_flag & 2), boundary submeshes are output by groups if more
 *    than one group is present and the default writer uses the EnSight format.
 *
 * parameters:
 *   check_flag <-- mask used for additional output
 *----------------------------------------------------------------------------*/

void
cs_post_init_main_meshes(int check_mask);

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, intended to
 * troubleshoot errors.
 *----------------------------------------------------------------------------*/

void
cs_post_init_error_writer(void);

/*----------------------------------------------------------------------------
 * Initialize post-processing writer with same format and associated
 * options as default writer, but no time dependency, and associate
 * and output global volume mesh.
 *
 * This is intended to help troubleshoot errors using fields based
 * on cells.
 *
 * returns:
 *   id of error output mesh (< 0), or 0 if all writers are deactivated
 *----------------------------------------------------------------------------*/

int
cs_post_init_error_writer_cells(void);

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
