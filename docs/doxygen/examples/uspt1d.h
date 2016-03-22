/*============================================================================
 * Code_Saturne documentation page
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

/*-----------------------------------------------------------------------------*/

/*!
  \page us_pt1d Data setting for the 1D-wall thermal module (uspt1d.f90)

  The \ref uspt1d subroutine is used to set the 1D-wall thermal
  module parameters.

  \section uspt1d_h_arg Arguments of uspt1d

  \snippet uspt1d.f90 arg

  \section uspt1d_h_loc_var Local variables declaration

  \snippet uspt1d.f90 loc_var_dec

  \section uspt1d_h_allocate Allocation

  \snippet uspt1d.f90 allocate

  \section uspt1d_h_restart Rereading of the restart file

  \snippet uspt1d.f90 restart

  \section uspt1d_h_iappel_12 iappel = 1 or 2

  \snippet uspt1d.f90 iappel_12

  \section uspt1d_h_iappel_2 iappel = 2

  \snippet uspt1d.f90 iappel_2

  \section uspt1d_h_iappel_3 iappel = 3

  \snippet uspt1d.f90 iappel_3

  \section uspt1d_h_deallocate Deallocation

  \snippet uspt1d.f90 deallocate

*/
