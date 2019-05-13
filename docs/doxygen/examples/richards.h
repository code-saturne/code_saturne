/*============================================================================
 * Code_Saturne documentation page
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

/*-----------------------------------------------------------------------------*/

/*!
  \page richards Data setting for the groundwater flow module

  \section richards_intro Introduction

  The Hydrogeology module of \CS is a numerical model for water flow and solute
  transport in continuous porous media. The flow part is based on the Richards
  equation, derived from the Darcy law and the conservation of mass. The
  transport part is based on the the classical advection-diffusion equation,
  slightly modified to account the specificities of underground transport.

  This module can be used to simulate transfers of water and solutes in several
  saturated and/or unsaturated porous media. The flow part can be steady or
  unsteady, with scalar or tensorial permeabilities and allows any type of soil
  water retention model, such as the van Genuchten model. The transport part
  considers dispersion, sorption and radioactive decay. The partition between
  soil and water phases can be modeled by a classical Kd approach or an
  alternative EK (equilibrium-kinetic) model. Additionaly solutes
  precipitation/dissolution phenomena can also be taken into account by an
  instantaneous model.

  Physical concepts and equations are presented in the \ref theory guide.

  The groundwater flow module is recent, and thus, has few limitations:
  - The weighted gradient computation, required for tetrahedral meshes and high
    permeability ratio, can only be used for isotropic soils.
  - Only one solute with anisotropic dispersion can be treated.

  \section richards_activ Activation of the module

  The module can be activated in the \ref usppmo routine in
  \ref cs_user_parameters.f90. The corresponding keyword is idarcy in the
  \ref optcal module:

  \snippet cs_user_parameters-richards.f90 richards_activ

  Note that the activation of the module requires to desactivation the turbulent
  model in \ref usipph routine in \ref cs_user_parameters.f90 file:

  \snippet cs_user_parameters-richards.f90 richards_warning

  \section richards_parameters Specific parameters

  When the module is activated, its specific input parameters should be set in
  the \ref user_darcy_ini1 routine of \ref cs_user_parameters.f90 file. An example
  is given in cs_user_parameters-richards.f90.

  \subsection richards_parameters_flow Flow part

  The permeability can be isotropic (scalar) or anisotropic (tensor) but all
  soils will be treated in the same way (isotropic or anisotropic):

  \snippet cs_user_parameters-richards.f90 richards_perm

  The primary variable of the groundwater flow module is the hydraulic head H=h+z.
  In order to switch easily to the pressure head, the keyword \ref darcy_gravity
  can be used and the value \ref darcy_gravity_x/y/z will defined the direction:

  \snippet cs_user_parameters-richards.f90 richards_grav

  The convergence criteron of the Newton scheme can be set over pressure or over
  velocity. It is recommended to keep the criteron over pressure:

  \snippet cs_user_parameters-richards.f90 richards_conv

  \subsection richards_parameters_trpt Transport part

  The dispersion can be isotropic (scalar) or anisotropic (tensor) but all
  solutes will be treated in the same way (isotropic or anisotropic):

  \snippet cs_user_parameters-richards.f90 richards_disp

  The transient transport can be based on a steady or unsteasy darcian velocity
  field:

  \snippet cs_user_parameters-richards.f90 richards_steady

  The partition between solid and liquid phases can be modelled by a classical
  Kd approach or an alternative EK (equilibrium-kinetic) model.
  Additionally solutes precipitation/dissolution phenomena can also be taken
  into account by an instantaneous model.

  \snippet cs_user_parameters-richards.f90 richards_partition

  The radioactive decay is treated as a source term in the transport equation
  (in \ref covofi.f90). It is set in \ref cs_user_parameters.f90 as follows:

  \snippet cs_user_parameters-richards.f90 richards_decay

  \section richards_numerics Numerical parameters

  Specific numerical parameters can be set in \ref usipsu routine of \ref
  cs_user_parameters.f90 file. An example is given in \ref
  cs_user_parameters-richards.f90:

  \subsection richards_numerics_flow Flow part

  In the case of soils of very different permeabilities, the resolution of the
  Richards equation requires a weighted gradient computation for tetrahedral
  meshes.
  This option is only available for soils with isotropic permeabilities
  (\ref darcy_anisotropic_permeability = 0) for now.

  \snippet cs_user_parameters-richards.f90 richards_iwgrec

  It is recommended to choose low criteria for gradient reconstruction in order
  to obtain a smooth darcian velocity field for the transport part. For
  instance:

  \snippet cs_user_parameters-richards.f90 richards_num_flow

  \subsection richards_numerics_trpt Transport part

  In the case of soils of very different diffusion (dispersion or molecular
  diffusion), the resolution of the transport equation requires a weighted
  gradient computation for tetrahedral meshes.

  \snippet cs_user_parameters-richards.f90 richards_num_trpt

  \subsection richards_numerics_time Time parameters

  The total number of iterations and the reference time step are also set in this
  routine.

  \snippet cs_user_parameters-richards.f90 richards_num_time

  However, the time step can be modified in \ref cs_user_extra_operations.f90
  (see \ref cs_user_extra_operations-flux.f90) in order to modif the time step
  with time:

  \snippet cs_user_extra_operations-richards_flux.f90 richards_time_modif

  \section richards_phys_prop Physical parameters

  Physical parameters can be set in \ref usphyv routine of \ref
  cs_user_physical_properties.f90 file. This section presents two examples that
  can be found in in \ref cs_user_physical_properties-richards_sat.f90 and \ref
  cs_user_physical_properties-richards_unsat.f90.

  Note that, in the flow part, depending on the variable \ref
  darcy_anisotropic_permeability, the permeability storage table is \ref
  permeability (isotropic case) or \ref tensor_permeability (anisotropic case).
  For the transport part, the isotropic part of the diffusion (molecular
  diffusion and isotropic dispersion) is always stored in \ref cpro_vscalt.
  Only anisotropic dispersion (i.e. \ref darcy_anisotropic_diffusion = 1)
  is stored in the tensor \ref visten.

  \subsection richards_phys_prop_sat Case of two saturated soils

  This example shows how to set physical parameters for two fully saturated
  soils (isotropic or anisotropic permeability) and several solutes with
  isotropic dispersion:

  \subsubsection richards_phys_prop_sat_flow Flow part

  \snippet cs_user_physical_properties-richards_sat.f90 richards_flow_soils

  \subsubsection richards_phys_prop_sat_trpt Transport part

  For every solute, the isotropic dispersion should be set in all soils.
  For instance:

  \snippet cs_user_physical_properties-richards_sat.f90 richards_flow_solut

  \subsection richards_phys_prop_unsat Unsaturated media

  This example shows how to set physical parameters for a single variably
  saturated soil (isotropic or anisotropic permeability) and a single solute
  with molecular diffusion, anisotropic dispersivity and sorption. The van
  Genuchten model, coupled with the Mualem condition, is used to determine
  the relation between the moisture content and the presure head (h).

  \subsubsection richards_phys_prop_unsat_flow Flow part

  First the permeability and the van Genuchten parameters are set:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_set_genuch

  As the van Genuchten law is based on the pressure head (h), the gravity term
  is added if necessary:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_set_press

  In order to compute the capacity, the saturation and the permeability, the
  saturated and the unsaturated parts are treated differently.

  In the saturated part (h>=0), the water content is equal to the saturated
  water content, the permeability is equal to the saturated permeability and
  the capacity is equal to zero:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_sat_part

  In the unsaturated part (h<0), the water content, the capacity and the
  permeability are function of the pressure head. They are determined thanks
  to the van Genuchten law:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_unsat_part

  \subsubsection richards_phys_prop_unsat_trpt Transport part

  First, the values of the longitudinal and transversal dispersivity as well
  as the molecular diffusion of the single solute are set as follow:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_unsat_trpt_init

  The molecular diffusion (isotropic term) is stored in \ref cpro_vscalt and
  computed as:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_unsat_mol_diff

  The anisotropic dispersion is stored in \ref visten and computed as:

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_unsat_aniso_disp

  Sorption parameters like Kd, kplus or kminus (in case of EK model) can be set
  as follows. Soil density is also needed to compute the transport delay.
  If precipitation is taken into account, the solubility index needs to be set as well.

  \snippet cs_user_physical_properties-richards_unsat.f90 richards_unsat_soilwater_partition

  \section richards_source Source terms

  Source terms can be added to scalar transport equation in \ref ustsns routine
  of \ref cs_user_source_terms.f90 file.

  \subsection richards_phys_prop_chem_rel Chemicals release

  Substances can be gradually released within the soil.

  \snippet cs_user_source_terms-richards.f90 richards_leaching

  \section richards_init Initialisation

  The initialisation of the variables required for the flow part (hydraulic
  head H) and transport part (concentration c) can be done globally:

  \snippet cs_user_initialization-richards.f90 richards_init_cell

  or by selecting a precise soil:

  \snippet cs_user_initialization-richards.f90 richards_init_grp

  If EK model is considered, sorbed concentration must be initialized.

  \snippet cs_user_initialization-richards.f90 richards_init_sorb

  If precipitation phenomenon is taken into account, precipitated concentration
  has to be initialized.

  \snippet cs_user_initialization-richards.f90 richards_init_precip

  \section richards_bound_cond Boundary conditions

  For groundwater flows of water and solutes, the undefined type face \ref
  iindef is used to impose Dirichlet, Neumann and mixte boundary conditions
  on hydraulic head H (here pressure) and solutes. Several examples can be
  found in \ref cs_user_boundary_conditions-richards.f90.

  \subsection richards_bound_cond_diri Dirichlet boundary conditions

  Dirichlet boundary conditions can be used to impose a value for the hydraulic
  head H and the concentration c at a given boundary face:

  \snippet cs_user_boundary_conditions-richards.f90 richards_BC_ex1

  It can also be used to impose a hydraulic head profile at another face:

  \snippet cs_user_boundary_conditions-richards.f90 richards_BC_ex2

  \subsection richards_bound_cond_neum Neumann boundary conditions

  Neumann boundary conditions can be used to impose fluxes at boundaries:

  \snippet cs_user_boundary_conditions-richards.f90 richards_BC_ex3

  Note that, for the transport part, Neumann boundary conditions can
  only be used for boundary surface with outward or null normal flow.
  In both cases, the prescribed flux is the diffusive flux.

  \subsection richards_bound_cond_mixte Mixte boundary conditions

  The mixte boundary conditions (Robin) can be used to impose a concentration
  flux at an entrance (inward normal flow at a boundary face). The following
  example explains how to determine the two parameters of the mixte boundary
  in order to impose a total flux:

  \snippet cs_user_boundary_conditions-richards.f90 richards_BC_ex4

  \section richards_post_proc Flux computations

  In order to compute fluxes at innner or boundary surfaces, the file \ref
  cs_user_extra_operations.f90 can be used. An example can be found in \ref
  cs_user_extra_operations-richards_flux.f90. It shows how to compute the
  total fluxes at to different surface and how to write the evolution with
  time:

  \snippet cs_user_extra_operations-richards_flux.f90 richards_flux_comp

*/
