/*============================================================================
 * code_saturne documentation page
 *============================================================================*/

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

/*-----------------------------------------------------------------------------*/

/*!
  \page us_pt1d Data setting for the 1D-wall thermal module (cs_user_1d_wall_thermal.c)

  The \ref cs_user_1d_wall_thermal subroutine is used to set the 1D-wall thermal
  module parameters.

  \section uspt1d_h_loc_var Local variables declaration

  \snippet cs_user_1d_wall_thermal.c loc_var_dec

  \section uspt1d_h_allocate Allocation

  \snippet cs_user_1d_wall_thermal.c allocate

  \section uspt1d_h_restart Rereading of the restart file

  \snippet cs_user_1d_wall_thermal.c restart

  \section uspt1d_h_iappel_12 iappel = 1 or 2

  \snippet cs_user_1d_wall_thermal.c iappel_12

  \section uspt1d_h_iappel_2 iappel = 2

  \snippet cs_user_1d_wall_thermal.c iappel_2

  \section uspt1d_h_iappel_3 iappel = 3

  \snippet cs_user_1d_wall_thermal.c iappel_3

  \section uspt1d_h_deallocate Deallocation

  \snippet cs_user_1d_wall_thermal.c deallocate

*/
