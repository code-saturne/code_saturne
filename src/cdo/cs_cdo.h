#ifndef __CS_CDO_H__
#define __CS_CDO_H__

/*============================================================================
 * General functions or variables for the INNOV module
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

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

#define CS_CDO_LEN_NAME 64

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef unsigned short int cs_flag_t;

/* Type of numerical scheme */
typedef enum {

  CS_SPACE_SCHEME_NONE,
  CS_SPACE_SCHEME_CDOVB,   /* CDO scheme with vertex-based positionning */
  CS_SPACE_SCHEME_CDOFB,   /* CDO cell-based scheme with hybridization */
  CS_SPACE_N_SCHEMES

} cs_space_scheme_t;

/* Values associated to the different ways to retrieve data */
typedef union {

  cs_flag_t              flag;      // flag
  char                  *name;      // file name for instance
  int                    id;        // identification number
  cs_lnum_t              num;       // local number
  cs_real_t              val;       // value
  cs_real_2_t            couple;    // two values
  cs_real_3_t            vect;      // vector: 3 values
  cs_real_6_t            twovects;  // two vectors
  cs_real_33_t           tens;      // tensor: 9 values

} cs_get_t;

/* Analytic definition through a function */
typedef void
(cs_analytic_func_t)(cs_real_t     time,
                     cs_real_3_t   xyz,
                     cs_get_t     *retval);

typedef void
(cs_user_func_t) (const void         *input1,
                  const void         *input2,
                  cs_real_t           cur_time,
                  cs_real_3_t         xyz,
                  cs_get_t           *output);

/*============================================================================
 * Global variables
 *============================================================================*/

/* Zero threshold */
extern double  cs_base_zthreshold;

/* Separation lines: long, medium, short */
extern const char lsepline[];
extern const char msepline[];
extern const char ssepline[];

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Return a string "true" or "false" according to the boolean
 *
 * \param[in]  boolean     bool  type
 *
 * \return a string "true" or "false"
 */
/*----------------------------------------------------------------------------*/

const char *
cs_base_strtf(bool  boolean);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute epsilon which is the machine precision
 */
/*----------------------------------------------------------------------------*/

void
cs_set_eps_machine(void);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Get the machine precision
 */
/*----------------------------------------------------------------------------*/

double
cs_get_eps_machine(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CDO_H__ */
