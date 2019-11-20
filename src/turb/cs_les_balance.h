#ifndef __CS_LES_BALANCE_H__
#define __CS_LES_BALANCE_H__

/*============================================================================
 * Filters for dynamic models.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2018 EDF S.A.

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

#include "cs_base.h"
#include "cs_defs.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/* LES balance type masks */

/*! Active LES Rij balance */
#define CS_LES_BALANCE_RIJ               (1 << 0)

/*! basic LES Rij balance */
#define CS_LES_BALANCE_RIJ_BASE          (1 << 1)

/*! full LES Rij balance */
#define CS_LES_BALANCE_RIJ_FULL          (1 << 2)

/*! Active LES Tui balance */
#define CS_LES_BALANCE_TUI               (1 << 3)

/*! basic LES Tui balance */
#define CS_LES_BALANCE_TUI_BASE          (1 << 4)

/*! full LES Tui balance */
#define CS_LES_BALANCE_TUI_FULL          (1 << 5)

/*============================================================================
 * Type definition
 *============================================================================*/

/* Additional types */

/* Rij LES balance structure */
/*---------------------------*/

typedef struct {

  /* Working arrays */
  cs_real_t      *pp2;
  cs_real_t      *smagp2;
  cs_real_6_t    *prodij;
  cs_real_6_t    *phiij;
  cs_real_6_t    *epsij;
  cs_real_6_t    *difftij;
  cs_real_6_t    *difftpij;
  cs_real_6_t    *unstij;
  cs_real_6_t    *convij;
  cs_real_6_t    *difflamij;
  cs_real_6_t    *budsgsij;
  cs_real_96_t   *budsgsfullij;

} cs_les_balance_rij_t;

/* Tui LES balance structure */
/*---------------------------*/

typedef struct {

  /* Field id */
  int              f_id;

  /* Working arrays */
  cs_real_t        *unstvar;
  cs_real_t        *tptp;
  cs_real_t        *prodvar;
  cs_real_t        *epsvar;
  cs_real_t        *difftvar;
  cs_real_t        *convvar;
  cs_real_t        *difflamvar;
  cs_real_t        *budsgsvar;
  cs_real_3_t      *tpuip;
  cs_real_3_t      *unstti;
  cs_real_3_t      *prodtUi;
  cs_real_3_t      *prodtTi;
  cs_real_3_t      *phiti;
  cs_real_3_t      *epsti;
  cs_real_3_t      *difftti;
  cs_real_3_t      *diffttpi;
  cs_real_3_t      *convti;
  cs_real_3_t      *difflamti;
  cs_real_3_t      *budsgstui;
  cs_real_6_t      *budsgsvarfull;
  cs_real_3_t     **budsgstuifull;

} cs_les_balance_tui_t;

/* LES balance structure */
/*-----------------------*/

typedef struct {

  cs_les_balance_rij_t    *brij;

  cs_les_balance_tui_t   **btui;

  int                      i_les_balance;

  int                      type;

  int                      frequency_n;

} cs_les_balance_t;

/* Pointer to main LES balance structure */
extern cs_les_balance_t *cs_glob_les_balance;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create fields used in LES balance computation
 *----------------------------------------------------------------------------*/

void
cs_les_balance_create_fields(void);

/*----------------------------------------------------------------------------
 * Provide access to cs_glob_les_balance.
 *
 * returns:
 *  pointer to LES balance global structure
 *----------------------------------------------------------------------------*/

cs_les_balance_t *
cs_get_glob_les_balance(void);

/*----------------------------------------------------------------------------
 * Create a LES balance descriptor.
 *----------------------------------------------------------------------------*/

void
cs_les_balance_create(void);

/*----------------------------------------------------------------------------
 * Update gradients needed in LES balance
 *----------------------------------------------------------------------------*/

void
cs_les_balance_update_gradients(void);

/*----------------------------------------------------------------------------
 * Compute Rij LES balance.
 *----------------------------------------------------------------------------*/

void
cs_les_balance_compute_rij(void);

/*----------------------------------------------------------------------------
 * Compute Tui LES balance.
 *----------------------------------------------------------------------------*/

void
cs_les_balance_compute_tui(void);

/*----------------------------------------------------------------------------
 * Write the LES balance structure in the auxiliary restart file.
 *----------------------------------------------------------------------------*/

void
cs_les_balance_write_restart(void);

/*----------------------------------------------------------------------------
 * Destroy the LES balance structure.
 *----------------------------------------------------------------------------*/

void
cs_les_balance_finalize(void);

/*----------------------------------------------------------------------------
 * Active the LES balance module.
 *
 * parameters:
 *  type_flag   --> mask of LES balance type
 *  frequency_n --> balance computing frequency in time-steps
 *----------------------------------------------------------------------------*/

void
cs_les_balance_activate(int     type_flag,
                        int     frequency_n);

/*----------------------------------------------------------------------------
 * Compute the LES balance
 *----------------------------------------------------------------------------*/

void
cs_les_balance_compute(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_LES_BALANCE_H__ */
