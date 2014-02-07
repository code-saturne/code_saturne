/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
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

#ifndef __CS_LAGR_DLVO_H__
#define __CS_LAGR_DLVO_H__

/*============================================================================
 * Functions and types for the clogging modeling
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"
#include "cs_lagr_tracking.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

cs_real_t
van_der_waals_sphere_plane(                   cs_real_t distp,
                                              cs_real_t rpart,
                                              cs_real_t lambwl,
                                              cs_real_t cstham
  );

cs_real_t
EDL_sphere_plane (   cs_real_t distp,
                     cs_real_t rpart,
                     cs_real_t phi1,
                     cs_real_t phi2,
                     cs_real_t kboltz,
                     cs_real_t temp,
                     cs_real_t debye_length,
                     cs_real_t free_space_permit,
                     cs_real_t water_permit
  );

cs_real_t
van_der_waals_sphere_sphere(           cs_real_t              distcc,
                                       cs_real_t              rpart1,
                                       cs_real_t              rpart2,
                                       cs_real_t              lambwl,
                                       cs_real_t              cstham
  );


cs_real_t
EDL_sphere_sphere(            cs_real_t              distcc,
                              cs_real_t              rpart1,
                              cs_real_t              rpart2,
                              cs_real_t              phi1,
                              cs_real_t              phi2,
                              cs_real_t              kboltz,
                              cs_real_t              temp,
                              cs_real_t              debye_length,
                              cs_real_t              free_space_permit,
                              cs_real_t              water_permit
  );


END_C_DECLS

#endif /* __CS_LAGR_DLVO_H__ */

