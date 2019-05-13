#ifndef __CS_PARAMEDMEM_REMAPPER_HXX__
#define __CS_PARAMEDMEM_REMAPPER_HXX__

/*============================================================================
 * Parallel interpolation using ParaMEDMEM OverlapDEC
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

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

#if defined(HAVE_MEDCOUPLING_LOADER)

typedef struct _cs_paramedmem_remapper_t cs_paramedmem_remapper_t;

/* -------------------------------------------------------------------------- */
/*!
 * \brief   Creates a new cs_paramedmem_remapper_t instance
 *
 * \param[in] name          name of the remapper
 * \param[in] sel_criteria  cells selection criteria
 * \param[in] fileName      med file name
 * \param[in] meshName      name of the mesh in the med file
 *
 * \return  cs_paramedmem_remapper_t struct
 */
/* -------------------------------------------------------------------------- */

cs_paramedmem_remapper_t *
cs_paramedmem_remapper_create(char       *name,
                              const char *sel_criteria,
                              char        *fileName,
                              char        *meshName);

/* -------------------------------------------------------------------------- */
/*!
 * \brief Interpolate a given field on the local mesh for a given time
 *
 * \param[in] r             pointer to cs_paramedmem_remapper_t struct
 * \param[in] fieldName     name of the field to remap from the file
 * \param[in] time_choice   Choice of the time interpolation.
 *                          0: Value of field interpolated at t=tval from the
 *                          med file.
 *                          1: Returns field values for the first time step in
 *                          the file. tval is then ignored.
 *                          2: Returns field values for the last time step in
 *                          the file. tval is then ignored.
 * \param[in] tval          requested time instant. If time choice is 0 and
 *                          tval outside of the file time bounds, return value
 *                          will be at the the first time step (if tval < tmin)
 *                          or last time step (if tval > tmax)
 *
 * \return  cs_real_t pointer containing the new values on target mesh
 */
/* -------------------------------------------------------------------------- */

cs_real_t *
cs_paramedmem_remap_field(cs_paramedmem_remapper_t *r,
                          char                     *fieldName,
                          int                       time_choice,
                          double                    tval);

END_C_DECLS
#endif
#endif /* __CS_PARAMEDMEM_REMAPPER_HXX__ */
