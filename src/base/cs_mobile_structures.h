#ifndef __CS_MOBILE_STRUCTURES_H__
#define __CS_MOBILE_STRUCTURES_H__

/*============================================================================
 * Mobile structures management.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2025 EDF S.A.

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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

#include "base/cs_restart.h"
#include "base/cs_time_control.h"
#include "base/cs_time_plot.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS


/*! Mobile_structures type */
/*-------------------------*/

typedef struct {

  /* Base structure definitions and input */

  int n_int_structs; /*!< number of internal structures */

  bool has_ext_structs; /*!< has external structures ? */

  cs_real_t      aexxst;     /*!< coefficient for the predicted displacement */
  cs_real_t      bexxst;     /*!< coefficient for the predicted displacement */

  cs_real_t      cfopre;     /*!< coefficient for the predicted force */

  cs_real_t      alpnmk;     /*!< alpha coefficient for the Newmark hht method */
  cs_real_t      betnmk;     /*!< beta coefficient for the Newmark hht method */
  cs_real_t      gamnmk;     /*!< gamma coefficient for the Newmark hht method */

  cs_real_33_t  *xmstru;     /*!< mass matrices (kg) */
  cs_real_33_t  *xcstru;     /*!< damping matrix coefficients (kg/s) */
  cs_real_33_t  *xkstru;     /*!< spring matrix constants (kg/s2 = N/m) */

  /* Output (plotting) control */

  int                plot;               /*!< monitoring format mask
                                           0: no plot
                                           1: plot to text (.dat) format
                                           2: plot to .csv format *
                                           3: plot to both formats */

  cs_time_control_t  plot_time_control;  /*!< time control for plotting */
  char              *plot_dir_name;      /*!< monitoring output directory */

  /* Computed structure values */

  cs_real_3_t   *xstr;       /*!< displacement vectors compared to structure
                              *   positions in the initial mesh (m) */
  cs_real_3_t  *xsta;        /*!< values of xstr at the previous time step */
  cs_real_3_t  *xstp;        /*!< predicted values of xstr */
  cs_real_3_t  *xstreq;      /*!< equilibrum positions of a structure (m) */

  cs_real_3_t  *xpstr;       /*!< velocity vectors (m/s) */
  cs_real_3_t  *xpsta;       /*!< xpstr at previous time step */

  cs_real_3_t  *xppstr;      /*!< acceleration vectors (m/s2) */
  cs_real_3_t  *xppsta;      /*!< acceleration vectors at previous
                              *   time step (m/s2) */

  cs_real_3_t  *forstr;      /*!< force vectors acting on the structure (N) */
  cs_real_3_t  *forsta;      /*!< forstr at previous time step (N) */
  cs_real_3_t  *forstp;      /*!< predicted force vectors (N) */

  cs_real_t  *dtstr;         /*!< time step used to solve structure movements */
  cs_real_t *dtsta; /*!< previous time step used to solve structure movements */

  /* Association with mesh */

  int        *idfstr;        /*!< structure number associated to each
                              *   boundary face:
                              *   - 0 if face is not coupled to a structure
                              *   - if > 0, internal structure id + 1
                              *   - if < 0, - code_aster instance id  - 1 */

  /* Plotting */

  int            n_plots;    /*!< number of plots for format */

  cs_time_plot_t  **plot_files[2];  /*!< Associated plot files */

} cs_mobile_structures_t;

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

/*! Maximum number of implicitation iterations of the structure displacement */
extern int cs_glob_mobile_structures_n_iter_max;

/*! Relative precision of implicitation of the structure displacement */
extern double cs_glob_mobile_structures_i_eps;

extern cs_mobile_structures_t  *_mobile_structures;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Initialize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_initialize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Finalize mobile structures with ALE for internal coupling.
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_finalize(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Log structures and coupling information
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_log_setup(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query number of internal mobile structures defined.
 *
 * \return  number of internal mobile structures
 */
/*----------------------------------------------------------------------------*/

int
cs_mobile_structures_get_n_int_structures(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Query number of external mobile structures defined.
 *
 * \return  number of external mobile structures
 */
/*----------------------------------------------------------------------------*/

int
cs_mobile_structures_get_n_ext_structures(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add internal mobile structures.
 *
 * This function may be called multiple time to change the number of
 * mobile structures.
 *
 * \param[in]   n_structures  number of internal mobile structures
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_add_n_structures(int  n_structures);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add external mobile structures.
 *
 * This function may be called multiple time to change the number of
 * mobile structures.
 *
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_add_external_structures(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set Newmark coefficients for internal mobile structures.
 *
 * \param[in]   alpha  alpha coefficient for Newmark algorithm
 * \param[in]   beta   beta coefficient for Newmark algorithm
 * \param[in]   gamma  gamma coefficient for Newmark algorithm
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_set_newmark_coefficients(cs_real_t  alpha,
                                              cs_real_t  beta,
                                              cs_real_t  gamma);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Predict displacement of mobile structures with ALE.
 *
 * \param[in]   bc_coeffs_vel   velocity boundary condition structure
 * \param[in]   itrale          ALE iteration number
 * \param[in]   italim          implicit coupling iteration number
 * \param[in]   ineefl          indicate whether fluxes should be saved
 * \param[out]  impale          imposed displacement indicator
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_prediction(cs_field_bc_coeffs_t *bc_coeffs_vel,
                                int                   itrale,
                                int                   italim,
                                int                   ineefl,
                                int                   impale[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Displacement of mobile structures with ALE for internal coupling.
 *
 * \param[in]       itrale   ALE iteration number
 * \param[in]       italim   implicit coupling iteration number
 * \param[in, out]  itrfin   indicator for last iteration of implicit coupling
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_displacement(int   itrale,
                                  int   italim,
                                  int  *itrfin);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Read mobile structures data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_restart_read(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Write mobile structures data to checkpoint.
 *
 * \param[in, out]  r   associated restart file pointer
 */
/*----------------------------------------------------------------------------*/

void
cs_mobile_structures_restart_write(cs_restart_t  *r);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_MOBILE_STRUCTURES_H__ */
