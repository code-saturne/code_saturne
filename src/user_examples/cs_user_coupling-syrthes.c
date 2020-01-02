/*============================================================================
 * Code couplings definition with SYRTHES and Code_Saturne.
 *
 * 1) Define conjuguate heat transfer couplings with the SYRTHES code
 * 2) Define couplings with other instances of Code_Saturne
 *============================================================================*/

/* VERS */

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

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
 * \file cs_user_coupling-syrthes.c
 *
 * \brief Code couplings definition with SYRTHES and Code_Saturne.
 *
 * See \subpage user_coupling for examples.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define couplings with SYRTHES code.
 *
 * This is done by calling the \ref cs_syr_coupling_define function for each
 * coupling to add.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_syrthes_coupling(void)
{
  /*! [coupling_syrthes_1] */
  {
    int  verbosity = 1, plot = 1;
    float tolerance = 0.1;
    bool allow_nonmatching = false;

    /*-------------------------------------------------------------------------
     * Example 1:
     *
     * Boundary faces of group '3' coupled with instance named 'SYRTHES_01'.
     *-------------------------------------------------------------------------*/

    cs_syr_coupling_define("SYRTHES_01",
                           "3",               /* boundary criteria */
                           NULL,              /* volume_criteria */
                           ' ',               /* projection_axis */
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           plot);

  }
  /*! [coupling_syrthes_1] */

  /*! [coupling_syrthes_2] */
  {
    int  verbosity = 1, plot = 1;
    float tolerance = 0.1;
    bool allow_nonmatching = false;

    /*-------------------------------------------------------------------------
     * Example 2:
     *
     * Boundary faces of group 'Wall' coupled with 2D SYRTHES instance
     * named 'SYRTHES_02'.
     *-------------------------------------------------------------------------*/

    cs_syr_coupling_define("SYRTHES_02",
                           "Wall",            /* boundary criteria */
                           NULL,              /* volume_criteria */
                           'z',               /* projection_axis */
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           plot);

  }
  /*! [coupling_syrthes_2] */

  /*! [coupling_syrthes_3] */
  {
    int  verbosity = 1, plot = 1;
    float tolerance = 0.1;
    bool allow_nonmatching = false;

    /*-------------------------------------------------------------------------
     * Example 3:
     *
     * Cells in box with corners (0, 0, 0) and (1, 1, 1) coupled with
     * SYRTHES instance named 'Solid' (volume coupling).
     *-------------------------------------------------------------------------*/

    cs_syr_coupling_define("Solid",
                           NULL,                          /* boundary */
                           "box[0., 0., 0., 1., 1., 1.]", /* volume */
                           ' ',                           /* projection */
                           allow_nonmatching,
                           tolerance,
                           verbosity,
                           plot);
  }
  /*! [coupling_syrthes_3] */

  /* By default, conservativity forcing flag is switched off (value 0)
     If one wants to switch on the conservativity forcing flag:

     cs_syr_coupling_set_conservativity(1);
  */

  /* Only for a volume coupling:
      By default, implicit treatment is done. You can switch to
      an explicit treatment by using the following function:

     cs_syr_coupling_set_explicit_treatment();
  */

}

END_C_DECLS
