/*============================================================================
 * Code_Saturne documentation page
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

/*-----------------------------------------------------------------------------*/

/*!
  \page cs_user_modules Creating user arrays (cs_user_modules.f90)

  \section cs_user_module_h_intro Introduction

  User subroutine for the definition of a user module. The following example
  illustrates the creation of user arrays called \c rwork and \c iwork.

  \section cs_user_module_h_variable Arrays declaration

  The first step is to declare these two arrays at the very
  beginning of the module. The first one \c iwork is an allocatable
  array of one dimension whereas \c rwork is a pointer to a bidimensionnal
  array.

  \snippet cs_user_modules-user-arrays.f90 variables

  \section cs_user_module_h_allocate Arrays allocation

  The \ref init_user_module subroutine allocates \c rwork and \c iwork if they
  are not already allocated or associated (for \c rwork).

  \snippet cs_user_modules-user-arrays.f90 allocate

  \section cs_user_module_h_c_pointer Access to arrays in C

  It is possible to access the \c rwork array in the C part of the code.
  This can be done by using the \ref get_user_module_rwork subroutine.

  \snippet cs_user_modules-user-arrays.f90 c_pointer

  \section cs_user_module_h_free Arrays freeing

  Eventually, the \ref finalize_user_module subroutine is used to free
  \c iwork and \c rwork.

  \snippet cs_user_modules-user-arrays.f90 free

*/
