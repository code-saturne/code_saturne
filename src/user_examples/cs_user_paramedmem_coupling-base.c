/*============================================================================
 * User functions for input of ParaMEDMEM coupling parameters
 *============================================================================*/

/* VERS */

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
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_paramedmem_coupling.c
 *
 * \brief User functions for input of ParaMEDMEM coupling parameters
 *
 * \brief User functions for input of calculation parameters.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define ParaMEDMEM coupling(s)
 */
/*----------------------------------------------------------------------------*/

void
cs_user_paramedmem_define_couplings(void)
{

  /* Define a coupling using ParaMEDMEM.
   * Coupled instances (apps) are "SAT" and "NCFD".
   * Coupling name is "CPL1".
   */
  /*! [paramedmem_coupling_define1] */
  {
    cs_paramedmem_coupling_t *c = cs_paramedmem_coupling_create("SAT",
                                                                "NCFD",
                                                                "CPL1");
  }
  /*! [paramedmem_coupling_define1] */

  /* Define a coupling using ParaMEDMEM.
   * First parameter is set to NULL, hence app #1 is this instance.
   * Second instance is the app named "PARTNER".
   * Third parameter is NULL, hence the coupling name will be set automatically
   * to '<app1_name>_<app2_name>_cpl'
   */
  /*! [paramedmem_coupling_define2] */
  {
    cs_paramedmem_coupling_t *c = cs_paramedmem_coupling_create(NULL,
                                                                "PARTNER",
                                                                NULL);
  }
  /*! [paramedmem_coupling_define2] */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define coupled meshes
 */
/*----------------------------------------------------------------------------*/

void
cs_user_paramedmem_define_meshes(void)
{

  /* Define the coupled mesh for a given coupling based on selection criteria.
   *
   * Input paramters are :
   *    (1) cs_paramedmem_coupling_t pointer
   *    (2) selection criteria (const char)
   *    (3) dimension of mesh elements (3 for cells, 2 for faces)
   *
   * Here we work with the coupling named "CPL1".
   * We define a volume mesh , based on a selection criteria "x < 0.5".
   */
  /*! [paramedmem_coupling_define_mesh1] */
  {
    cs_paramedmem_coupling_t *c = cs_paramedmem_coupling_by_name("CPL1");
    cs_paramedmem_add_mesh_from_criteria(c, "x < 0.5", 3);
  }
  /*! [paramedmem_coupling_define_mesh1] */

  /* Define the coupled mesh for a given coupling based on a cs_zone_t pointer.
   *
   * Input paramters are :
   *    (1) cs_paramedmem_coupling_t pointer
   *    (2) cs_zone_t pointer
   */
  /*! [paramedmem_coupling_define_mesh2] */
  {
    cs_paramedmem_coupling_t *c = cs_paramedmem_coupling_by_name("CPL1");
    const cs_zone_t *z = cs_volume_zone_by_name("zone_pmm1");
    cs_paramedmem_add_mesh_from_zone(c, z);
  }
  /*! [paramedmem_coupling_define_mesh2] */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define fields to couple with ParaMEDMEM
 */
/*----------------------------------------------------------------------------*/

void
cs_user_paramedmem_define_fields(void)
{

  /* Define coupled arrays using ParaMEDMEM */
  /*! [paramedmem_coupling_define_field1] */
  {
    cs_paramedmem_coupling_t *c = cs_paramedmem_coupling_by_name("CPL1");

    /* Define an array.
     *
     * Input parameters are:
     *
     *  (1) cs_paramedmeme_coupling_t pointer
     *  (2) field array associated name ("f1" here)
     *  (3) field array number of components (1 here)
     *  (4) Interpolation method for the field. Choice indicates whether
     *      field is intensive or extensive, and whether interpolation shoud
     *      ensure the maximum principle or volumic integral.
     *      !! WARNING !!
     *      This choice should be the same for both codes communicating
     *      !! WARNING !!
     *      Options are:
     *        CS_MEDCPL_FIELD_EXT_CONSERVATION
     *        CS_MEDCPL_FIELD_EXT_MAXIMUM
     *        CS_MEDCPL_FIELD_INT_CONSERVATION
     *        CS_MEDCPL_FIELD_INT_MAXIMUM
     *  (5) Localisation of values:
     *        CS_MEDCPL_ON_CELLS -> Cell centers
     *        CS_MEDCPL_ON_NODES -> Vertices
     *  (6) Time discretisation for coupling. options are:
     *        CS_MEDCPL_NO_TIME
     *        CS_MEDCPL_ONE_TIME
     *        CS_MEDCPL_LINEAR_TIME
     *
     * function returns the index for the field array.
     */
    int f_id1 = cs_paramedmem_def_coupled_field(c,
                                                "f1",
                                                1,
                                                CS_MEDCPL_FIELD_INT_CONSERVATION,
                                                CS_MEDCPL_ON_CELLS,
                                                CS_MEDCPL_NO_TIME);

    /* Define a coupled array based on a cs_field_t pointer.
     *
     * Input paramters are:
     *
     *  (1) cs_paramedmeme_coupling_t pointer (c here)
     *  (2) cs_field_t pointer (f here)
     *  (3) Interpolation method for the field. Choice indicates whether
     *      field is intensive or extensive, and whether interpolation shoud
     *      ensure the maximum principle or volumic integral.
     *      !! WARNING !!
     *      This choice should be the same for both codes communicating
     *      !! WARNING !!
     *      Options are:
     *        CS_MEDCPL_FIELD_EXT_CONSERVATION
     *        CS_MEDCPL_FIELD_EXT_MAXIMUM
     *        CS_MEDCPL_FIELD_INT_CONSERVATION
     *        CS_MEDCPL_FIELD_INT_MAXIMUM
     *  (4) Time discretisation for coupling. options are:
     *        CS_MEDCPL_NO_TIME
     *        CS_MEDCPL_ONE_TIME
     *        CS_MEDCPL_LINEAR_TIME
     */
    cs_field_t *f = cs_field_by_name("coupling_field1");
    int f_id2 =
      cs_paramedmem_def_coupled_field_from_cs_field(c, f,
                                                    CS_MEDCPL_FIELD_INT_CONSERVATION,
                                                    CS_MEDCPL_NO_TIME);
  }
  /*! [paramedmem_coupling_define_field1] */

}

END_C_DECLS
