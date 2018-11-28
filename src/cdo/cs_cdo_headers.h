#ifndef __CS_CDO_HEADERS_H__
#define __CS_CDO_HEADERS_H__

/*============================================================================
 * Global Code_Saturne headers file for easier include
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


/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_advection_field.h"
#include "cs_basis_func.h"
#include "cs_cdo_advection.h"
#include "cs_cdo_bc.h"
#include "cs_cdo_connect.h"
#include "cs_cdo_diffusion.h"
#include "cs_cdo_headers.h"
#include "cs_cdo_local.h"
#include "cs_cdo_main.h"
#include "cs_cdo_quantities.h"
#include "cs_cdo_time.h"
#include "cs_cdofb_ac.h"
#include "cs_cdofb_navsto.h"
#include "cs_cdofb_priv.h"
#include "cs_cdofb_scaleq.h"
#include "cs_cdofb_uzawa.h"
#include "cs_cdofb_vecteq.h"
#include "cs_cdovb_priv.h"
#include "cs_cdovb_scaleq.h"
#include "cs_cdovb_vecteq.h"
#include "cs_cdovcb_scaleq.h"
#include "cs_dbg.h"
#include "cs_domain.h"
#include "cs_domain_op.h"
#include "cs_domain_setup.h"
#include "cs_equation.h"
#include "cs_equation_bc.h"
#include "cs_equation_common.h"
#include "cs_equation_param.h"
#include "cs_equation_priv.h"
#include "cs_evaluate.h"
#include "cs_flag.h"
#include "cs_gwf.h"
#include "cs_gwf_soil.h"
#include "cs_gwf_tracer.h"
#include "cs_hho_builder.h"
#include "cs_hho_scaleq.h"
#include "cs_hho_stokes.h"
#include "cs_hho_vecteq.h"
#include "cs_hodge.h"
#include "cs_mesh_deform.h"
#include "cs_navsto_coupling.h"
#include "cs_navsto_param.h"
#include "cs_navsto_system.h"
#include "cs_param.h"
#include "cs_param_cdo.h"
#include "cs_property.h"
#include "cs_quadrature.h"
#include "cs_reco.h"
#include "cs_scheme_geometry.h"
#include "cs_sdm.h"
#include "cs_sla.h"
#include "cs_source_term.h"
#include "cs_static_condensation.h"
#include "cs_walldistance.h"
#include "cs_xdef.h"
#include "cs_xdef_cw_eval.h"
#include "cs_xdef_eval.h"

#endif /* __CS_CDO_HEADERS_H__ */
