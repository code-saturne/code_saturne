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

/*----------------------------------------------------------------------------*/

/*!
  \page cs_var_dico Variables and structures reference (C and Fortran)

  \section cs_var_dico_intro Introduction

  This page is meant to help users find their way through when implementing
  user C functions or Fortran subroutines or even developing inside the code
  kernel. It provides cross-reference tables containing the names of Fortran
  variables and their C counterparts as well as guidelines about
  how to manage the mesh entities and the fields (variables, properties, ...).
  In addition, some naming conventions are described.

  Note: variables with the same name in Fortran and in C are most of the time not listed.

  \section cs_var_dico_namingcontent Contents

  Cross-reference tables and guidelines are organized as follows:

   - \subpage mesh
   - \subpage field
   - \subpage local

*/
// _____________________________________________________________________________
/*!

  \page mesh How to access and manage the mesh entities and mesh quantities ?

  \section cs_var_dico_mesh_vars Mesh variables

  These variables are defined in the files \ref mesh.f90 and \ref cs_mesh.h.

  - The Fortran variables are global in the code, accessible directly by their
    name.
  - Members of the global C structure \c cs_glob_mesh are accessed as: \n
    <tt>cs_glob_mesh->name</tt>, \n
    \e i.e. adding <tt> cs_glob_mesh-> </tt> in front of the name of the
    variable.

  Fortran code       | C code                                | Description
  ------------------ | ------------------------------------- | ------------
  <tt> \ref ndim     | cs_glob_mesh->dim                     | Space dimension
  \ref ncelet        | cs_glob_mesh->n_cells_with_ghosts     | Total number of cells on the local rank \n (n_cells + n_ghost_cells)
  \ref ncel          | cs_glob_mesh->n_cells                 | Number of cells
  \ref nfac          | cs_glob_mesh->n_i_faces               | Number of interior faces
  \ref nfabor        | cs_glob_mesh->n_b_faces               | Number of boundary faces
  \ref nnod          | cs_glob_mesh->n_vertices              | Number of vertices
  -                  | cs_glob_mesh->n_b_cells               | Number of boundary cells
  \ref lndfac        | cs_glob_mesh->i_face_vtx_connect_size | Size of the connectivity \n interior faces -> vertices
  \ref lndfbr        | cs_glob_mesh->b_face_vtx_connect_size | Size of the connectivity \n boundary faces -> vertices
  \ref nfml          | cs_glob_mesh->n_families              | Number of families
  \ref ifacel        | cs_glob_mesh->i_face_cells            | Interior faces -> cells connectivity
  \ref ifabor        | cs_glob_mesh->b_face_cells            | Boundary faces -> cells connectivity
  \ref ipnfac        | cs_glob_mesh->i_face_vtx_idx          | Interior faces -> vertices index
  \ref nodfac        | cs_glob_mesh->i_face_vtx_lst          | Interior faces -> vertices connectivity
  \ref ipnfbr        | cs_glob_mesh->b_face_vtx_idx          | Boundary faces -> vertices index
  \ref nodfbr        | cs_glob_mesh->b_face_vtx_lst          | Boundary faces -> vertices connectivity
  \ref ifmfbr        | cs_glob_mesh->b_face_family           | Boundary face family
  \ref ifmcel        | cs_glob_mesh->cell_family             | Cell family
  -                  | cs_glob_mesh->b_cells                 | Boundary cell list
  \ref xyznod        | cs_glob_mesh->vtx_coord               | Vertex coordinates </tt>

  \section cs_var_dico_mesh_q_vars Mesh quantity variables

  These variables are defined in the files \ref mesh.f90 and
  \ref cs_mesh_quantities.h.

  - The Fortran variables are global in the code, accessible directly by their
    name.
  - Members of the global C structure \c cs_glob_mesh_quantities are accessed
    as: \n
    <tt>cs_glob_mesh_quantities->name</tt>, \n
    \e i.e. adding <tt> cs_glob_mesh_quantities-> </tt> in front of the name of
    the variable.

  Fortran code       | C code                                  | Description
  ------------------ | --------------------------------------- |-------------
  <tt> \ref isympa   |  cs_glob_mesh_quantities->b_sym_flag    | Symmetry flag for boundary faces
  \ref xyzcen        |  cs_glob_mesh_quantities->cell_cen      | Cell center coordinates
  \ref surfac        |  cs_glob_mesh_quantities->i_face_normal | Surface normal of interior faces
  \ref surfbo        |  cs_glob_mesh_quantities->b_face_normal | Surface normal of border faces
  \ref cdgfac        |  cs_glob_mesh_quantities->i_face_cog    | Center of gravity of interior faces
  \ref cdgfbo        |  cs_glob_mesh_quantities->b_face_cog    | Center of gravity of border faces
  \ref volume        |  cs_glob_mesh_quantities->cell_vol      | Cell volume
  \ref surfan        |  cs_glob_mesh_quantities->i_face_surf   | Surface of interior faces
  \ref surfbn        |  cs_glob_mesh_quantities->b_face_surf   | Surface of boundary faces
  \ref dist          |  cs_glob_mesh_quantities->i_dist        | Distance between the cell center and \n the center of gravity of interior faces
  \ref distb         |  cs_glob_mesh_quantities->b_dist        | Distance between the cell center and \n the center of gravity of border faces
  \ref pond          |  cs_glob_mesh_quantities->weight        | Interior faces weighting factor </tt>


*/
// _____________________________________________________________________________
/*!

  \page field How to access and manage variables and properties using the cs_field API and the deprecated array \c propce ?

  \ref cs_var_dico_vars "Variables" and \ref cs_var_dico_props "properties" can be accessed both in Fortran and in C using the \ref field.f90 "cs_field" API. Some Fortran
  properties can also still be accessed through the deprecated array \c propce.

  \par Accessing variables and properties in Fortran:

    - Both \ref cs_var_dico_vars "variables" and \ref cs_var_dico_props "properties" can be accessed via the
      \ref field.f90 "cs_field" API, as in the following examples: \n
      - For one-dimensional arrays :\n\n
      <tt>call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref numvar::ipr "ipr"), cvar_pr)
          \n pres =  cvar_pr(iel)</tt>, \n\n
      <tt>call \ref field::field_get_val_s_by_name "field_get_val_s_by_name"("pressure", cvar_pr)
          \n pres =  cvar_pr(iel)</tt>, \n\n
      <tt>call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref cstphy::icp "icp"), cpro_cp) \n
          cp = cpro_cp(iel)</tt>, \n\n
          The scalar values are accessed as follows:\n\n
       <tt>call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref isca "isca"(iscalt)), cvar_scalt) \n
          temp = cvar_scalt(iel)</tt>, \n\n
      - For interleaved multidimensional arrays:\n\n
            <tt>call \ref field::field_get_val_v "field_get_val_v"(ivarfl(\ref numvar::iu "iu"), cvar_vel)
          \n ux = cvar_vel(1,iel)</tt>, \n\n
      .
       where \ref numvar::ipr "ipr", \ref numvar::iu "iu" are variable indexes,
       \ref cstphy::icp "icp" is a property index and \ref isca "isca"(iscalt) is a scalar index.\n\n
    - \ref cs_var_dico_props "Properties" can alternatively be accessed through the
       deprecated array \c propce as: \n
      <tt>propce(iel,index)</tt>.

  \par Accessing variables and properties in C:

    - Almost all \ref cs_var_dico_vars "variables" and \ref cs_var_dico_props "properties" can be accessed using the \ref CS_F_ macro: \n
     - For one-dimensional arrays :\n\n
      <tt>press = CS_F_(p)->val[cell_id]</tt>, \n\n
      <tt>cp = CS_F_(cp)->val[cell_id]</tt>, \n\n
      <tt>temp = CS_F_(t)->val[cell_id]</tt>, \n\n
     - For multidimensional arrays:\n\n
      <tt>uz = CS_F_(u)->val[3*cell_id + 2]</tt>\n\n
      These arrays can also be casted as follows (for a 3-D array):\n\n
      <tt>\ref cs_real_3_t *cvar_vel = (\ref cs_real_3_t *)CS_F_(u)->val</tt> \n\n
      The cell value can then be accessed as : \n\n
      <tt>ux = cvar_vel[cell_id][0]</tt>\n\n
     .
      <b>\c u, \c p, or \c cp </b> are defined in \ref cs_field_pointer.h. \n
    - Indexed variables (such as user scalars) and indexed properties
      are accessed as: \n
      <tt>\ref CS_FI_(name,ii-1)->val[cell_id]</tt>.\n
    - In any case, all variables can be accessed using the function \ref cs_field_by_name :\n
      <tt>\ref cs_real_t *cvar_pr = \ref cs_field_by_name("pressure")->val</tt>

  \remark Note that indexes in C begin at 0, while indexes in Fortran begin at 1.
          Thus, Fortran and C loop counters are related in the following as:\n
          <tt>cell_id = iel-1</tt>.

  Cross-reference tables are available for the variables and properties of the
  standard solver and the specific physics features:
    - \ref cs_var_dico_vars
    - \ref cs_var_dico_props
    - \ref cs_var_dico_part
    - \ref cs_var_dico_atmo
    - \ref cs_var_dico_comb
    - \ref cs_var_dico_cfbl
    - \ref cs_var_dico_elec
    - \ref cs_var_dico_cogz
    - \ref cs_var_dico_rayt

  \section cs_var_dico_vars Variables

  The Fortran variables indexes are defined in the files \ref numvar.f90 (with
  the exception of \c ihm and \c iscal, which are respectively defined in \ref
  ppincl.f90 and \ref optcal.f90) and the C variables names are defined in
  \ref cs_field_pointer.h. \n Note that \c dt is just an \c allocatable array in
  Fortran while it is mapped as a field in C.

  Fortran code                             | C code                       | Description
  ------------------------------------------------ | ---------------------------- | ------------
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ipr), cvar_pr)       | CS_F_(p)->val       | Pressure
  call \ref field::field_get_val_v "field_get_val_v"(ivarfl(\ref iu), cvar_vel)       | CS_F_(u)->val       | Velocity
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ivoidf), cvar_voidf) | CS_F_(void_f)->val  | Void fraction for cavitation modelling
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ik  ), cvar_k  )     | CS_F_(k)->val       | Turbulent kinetic energy \f$ k \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref iep ), cvar_eps)     | CS_F_(eps)->val     | Turbulent dissipation \f$ \varepsilon \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ir11), cvar_r11)     | CS_F_(r11)->val     | Reynolds stress component \f$ R_{xx} \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ir22), cvar_r22)     | CS_F_(r22)->val     | Reynolds stress component \f$ R_{yy} \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ir33), cvar_r33)     | CS_F_(r33)->val     | Reynolds stress component \f$ R_{zz} \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ir12), cvar_r12)     | CS_F_(r12)->val     | Reynolds stress component \f$ R_{xy} \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ir23), cvar_r23)     | CS_F_(r23)->val     | Reynolds stress component \f$ R_{yz} \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ir13), cvar_r13)     | CS_F_(r13)->val     | Reynolds stress component \f$ R_{xz} \f$
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref iphi), cvar_phi)     | CS_F_(phi)->val     | \f$ \phi \f$ for \f$ \phi-f_b \f$ model
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ifb ), cvar_fb )     | CS_F_(f_bar)->val   | \f$ f_b \f$ for \f$ \phi-f_b \f$ model
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref ial ), cvar_al )     | CS_F_(alpha)->val   | \f$ \alpha \f$ for \f$ Bl-v^2-k \f$ \n or EBRSM model
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref iomg), cvar_omg)     | CS_F_(omg)->val     | \f$ \omega \f$ for \f$ k-\omega \f$ SST model
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref inusa), cvar_nusa)   | CS_F_(nusa)->val    | \f$ \widetilde{\nu}_T \f$ for Spalart-Allmaras
  call \ref field::field_get_val_v "field_get_val_v"(ivarfl(\ref iuma), cvar_mesh_v)  | CS_F_(mesh_u)->val  | Mesh velocity
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(isca(\ref ppincl::ihm "ihm")), cvar_hm) | CS_F_(h)->val | Enthalpy
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(isca(\ref optcal::iscalt "iscalt")), cvar_scalt) | CS_F_(t)->val | Temperature </tt>


  \section cs_var_dico_props Properties

  These properties are defined in the files \ref numvar.f90 and
  \ref cs_field_pointer.h.

  Deprecated Fortran code  | Current Fortran code                                                      |C code                             | Description
  ------------------------ | ------------------------------------------------------------------------- |---------------------------------- | ------------
  <tt> dt                  | dt                                                                        |CS_F_(dt)->val            | Local time step
  propce(iel,irom)         | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::irom "irom"), cpro_rho  )        |CS_F_(rho)->val           | Density at the current time step
  propce(iel,irom)         | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::iroma "iroma"), cpro_rhoa )      |CS_F_(rho)->val_pre       | Density at the previous time step
  propce(iel,iromaa)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::iromaa "iromaa"), cpro_rhoaa)    |\ref cs_real_t *cpro_rhoaa = \ref cs_field_by_name("density_old")->val    | Density at the second previous time
  propce(iel,iviscl)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::iviscl "iviscl"), cpro_viscl)    |CS_F_(mu)->val            | Molecular viscosity
  propce(iel,ivisct)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ivisct "ivisct"), cpro_visct)    |CS_F_(mu_t)->val          | Turbulent dynamic viscosity
  propce(iel,ivisla)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ivisla "ivisla"), cpro_romaa)    |CS_F_(mu)->val_pre        | Dynamic molecular viscosity (in kg/(m.s)) \n at the previous time-step
  propce(iel,ivista)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ivista "ivista"), cpro_viscta)   |CS_F_(mu_t)->val_pre      | Dynamic turbulent viscosity \n at the previous time-step
  propce(iel,icp)          | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref cstphy::icp "icp"), cpro_cp)             |CS_F_(cp)->val            | Specific heat
  propce(iel,icpa)         | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::icpa "icpa"), cpro_cpa)          |CS_F_(cp)->val_pre        | specific heat at the previous time-step
  propce(iel,icrom)        | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::icrom "icrom"), cpro_crom)       |CS_F_(rho)->val           | Density (at cells)
  propfb(ifac,ibrom)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ibrom "ibrom"), bpro_rho)        |CS_F_(rho_b)->val[face_id]         | Density (at boundary faces)
  propce(iel,ismago)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ismago "ismago"), cpro_smago)    |\ref cs_real_t *cpro_smago = \ref cs_field_by_name("smagorinsky_constant^2")->val| Field id of the anisotropic turbulent viscosity
  propce(iel,icour)        | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::icour "icour"), cpro_cour)       |\ref cs_real_t *cpro_cour = \ref cs_field_by_name("courant_number")->val | Courant number
  propce(iel,ifour)        | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ifour "ifour"), cpro_four)       |\ref cs_real_t *cpro_four = \ref cs_field_by_name("fourier_number")->val | Fourier number
  propce(iel,iprtot)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::iprtot "iprtot"), cpro_prtot)    |\ref cs_real_t *cpro_prtot = \ref cs_field_by_name("total_pressure")->val | Total pressure at cell centers
  propce(iel,ivisma(1)) \n propce(iel,ivisma(2)) \n propce(iel,ivisma(3)) | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ivisma "ivisma"(1)), cpro_vism1) \n call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ivisma "ivisma"(2)), cpro_vism2) \n call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ivisma "ivisma"(3)), cpro_vism3) | \ref cs_real_t *cpro_vism1  = \ref cs_field_by_name("mesh_viscosity_1")->val \n \ref cs_real_t *cpro_vism2 = \ref cs_field_by_name("mesh_viscosity_2")->val \n \ref cs_real_t *cpro_vism3  = \ref cs_field_by_name("mesh_viscosity_3")->val | Mesh velocity viscosity for the ALE module
  propce(iel,itsrho)       | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::itsrho "itsrho"), cpro_tsrho )   |\ref cs_real_t *cpro_tsrho = \ref cs_field_by_name("dila_st")->val        | Global dilatation source terms
  propce(iel,ibeta)        | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref numvar::ibeta "ibeta"), cpro_beta  )     |\ref cs_real_t *cpro_beta = \ref cs_field_by_name("thermal_expansion")->val | Thermal expansion coefficient
  -                        | call \ref field::field_get_val_s "field_get_val_s"(\ref numvar::ipori "ipori", cpro_ipori)              |CS_F_(poro)->val                   | Porosity
  -                        | call \ref field::field_get_val_v "field_get_val_v"(\ref numvar::iporf "iporf", cpro_iporf)              |CS_F_(t_poro)->val                 | Tensorial porosity
  -                        | call \ref field::field_get_val_v "field_get_val_v"(\ref numvar::iforbr "iforbr", bpro_forbr)            |\ref cs_real_t *bpro_forbr = \ref cs_field_by_name("boundary_forces")->val| Field id of the stresses at boundary
  -                        | call \ref field::field_get_val_s "field_get_val_s"(\ref numvar::iyplbr "iyplbr", bpro_yplus)            |\ref cs_real_t *bpro_yplus = \ref cs_field_by_name("yplus")->val          | Field id of \f$y^+\f$ at boundary
  -                        | call \ref field::field_get_val_v "field_get_val_v"(\ref numvar::idtten "idtten", dttens)                |\ref cs_real_t *dttens = \ref cs_field_by_name("dttens")->val         | Field id for the dttens tensor
  -                        | call \ref field::field_get_val_s "field_get_val_s"(\ref numvar::itempb "itempb", t_b)                   |CS_F_(t_b)->val                   | Boundary temperature </tt>


  \section cs_var_dico_part Specific physics

  \subsection cs_var_dico_atmo Atmospheric

  Defined in \ref optcal.f90, \ref atincl.f90, \ref atvarp.f90 and
  \ref cs_field_pointer.h.

  Fortran code                                                                               | C code                                 | Description
  ------------------------------------------------------------------------------------------ | -------------------------------------- | ------------
  <tt> call \ref field::field_get_val_s "field_get_val_s"(ivarfl(isca(iscalt)), cvar_scalt)                                | CS_F_(pot_t)->val             | Potential temperature
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref isca(\ref atincl::itotwt "itotwt")), cvar_totwt)          | CS_F_(totwt)->val             | Total water content
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref isca(\ref atincl::intdrp "intdrp")), cvar_intdrp)         | CS_F_(ntdrp)->val             | Total number of droplets
  call \ref field::field_get_val_s "field_get_val_s"(ivarfl(\ref isca(\ref atchem::isca_chem "isca_chem"(iesp))), cvar_sc) | \ref CS_FI_(chemistry,iesp-1)->val | Chemistry species (indexed) </tt>


  \subsection cs_var_dico_comb Coal combustion

  Defined in \ref ppincl.f90, \ref ppcpfu.f90 and \ref cs_field_pointer.h.

  Fortran code                                                                    | C code                           | Description
  ------------------------------------------------------------------------------- | -------------------------------- | ------------
  <tt> call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::inp "inp"(iesp)), cvar_inpcl)  | \ref CS_FI_(np,iesp-1)->val  | Particles per kg for coal class
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ixch "ixch"(iesp)), cvar_xchcl)     | \ref CS_FI_(xch,iesp-1)->val | Reactive coal mass fraction for coal class
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ixck "ixck"(iesp)), cvar_xckcl)     | \ref CS_FI_(xck,iesp-1)->val | Coke mass fraction for coal class
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ixwt "ixwt"(iesp)), cvar_xwtcl)     | \ref CS_FI_(xwt,iesp-1)->val | Water mass fraction for coal class
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ih2 "ih2"(iesp)), cvar_h2cl)        | \ref CS_FI_(h2,iesp-1)->val  | Mass enthalpy for coal class (permeatic case)
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if1m "if1m"(iesp)), cvar_f1mcl)     | \ref CS_FI_(f1m,iesp-1)->val | Mean value light volatiles for coal class
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if2m "if2m"(iesp)), cvar_f2mcl)     | \ref CS_FI_(f2m,iesp-1)->val | Mean value heavy volatiles for coal class
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if4m "if4m"), cvar_f4m)             | CS_F_(f4m)->val         | Oxydant 2 mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if5m "if5m"), cvar_f5m))            | CS_F_(f5m)->val         | Oxydant 3 mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if6m "if6m"), cvar_f6m))            | CS_F_(f6m)->val         | Water from coal drying mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if7m "if7m"), cvar_f7m))            | CS_F_(f7m)->val         | Carbon from coal oxidyzed by O2 mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if8m "if8m"), cvar_f8m))            | CS_F_(f8m)->val         | Carbon from coal gasified by CO2 mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::if9m "if9m"), cvar_f9m))            | CS_F_(f9m)->val         | Carbon from coal gasified by H2O mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ifvp2m "ifvp2m"), cvar_fvp2m)       | CS_F_(fvp2m)->val       | f1f2 variance
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppcpfu::iyco2 "iyco2"), cvar_yco2)          | CS_F_(yco2)->val        | CO2 fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppcpfu::iyhcn "iyhcn"), cvar_yhnc)          | CS_F_(yhcn)->val        | HCN fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppcpfu::iyno "iyno"), cvar, yno)            | CS_F_(yno)->val         | NO fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppcpfu::iynh3 "iynh3"), cvar_ynh3)          | CS_F_(ynh3)->val        | NH3 enthalpy
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppcpfu::ihox "ihox"), cvar_hox)             | CS_F_(hox)->val         | Ox enthalpy </tt>


  \subsection cs_var_dico_cfbl Compressible

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                                                                   | C code                        | Description
  ------------------------------------------------------------------------------ | ----------------------------- | ------------
  <tt> call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ienerg "ienerg"), cvar_energ) | CS_F_(energy)->val   | Total energy
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::itempk "itempk"), cvar_tempk)      | CS_F_(t_kelvin)->val | Temperature, in Kelvin </tt>


  \subsection cs_var_dico_elec Electric arcs

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                                                                  | C code                             | Description
  ----------------------------------------------------------------------------- | ---------------------------------- | ------------
  <tt> call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ipotr "ipotr"), cvar_potr)   | CS_F_(potr)->val          | Electric potential, real part
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ipoti "ipoti"), cvar_poti)        | CS_F_(poti)->val          | Electric potential, imaginary part
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ipotva "ipotva(1)"), cvar_potva1) \n call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ipotva "ipotva(2")), cvar_potva2) \n call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ipotva "ipotva(3)"), cvar_potva3) | CS_F_(potva)->val | Vector potential
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::iycoel "iycoel"(iesp)), cvar_ycoel(iesp)) | \ref CS_FI_(ycoel,iesp-1)->val | Constituent mass fraction </tt>


  \subsection cs_var_dico_cogz Gas combustion

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                                                              | C code                     | Description
  ------------------------------------------------------------------------- | -------------------------- | ------------
  <tt> call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ifm "ifm"), cvar_fm)     | CS_F_(fm)->val    | Mixture fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ifp2m "ifp2m"), cvar_fp2m)    | CS_F_(fp2m)->val  | Mixture fraction variance
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::ifsm "ifsm"), cvar_fsm)       | CS_F_(fsm)->val   | Soot mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::inpm "inpm"), cvar_npm)       | CS_F_(npm)->val   | Soot precursor number
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::iygfm "iygfm"), cvar_ygfm)    | CS_F_(ygfm)->val  | Fresh gas fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::iyfm "iyfm"), cvar_yfm)       | CS_F_(yfm)->val   | Mass fraction
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::iyfp2m "iyfp2m"), cvar_yfp2m) | CS_F_(yfp2m)->val | Mass fraction variance
  call \ref field::field_get_val_s "field_get_val_s"(\ref isca(\ref ppincl::icoyfp "icoyfp"), cvar_coyfp) | CS_F_(coyfp)->val | Mass fraction covariance </tt>


  \subsection cs_var_dico_rayt Radiative transfer

  Defined in \ref radiat.f90 and \ref cs_field_pointer.h.

  Deprecated Fortran code                       | Fortran code                                                              | C code                       | Description
  --------------------------------------------- | ------------------------------------------------------------------------- | ---------------------------- | ------------
  <tt> propce(iel,\ref radiat::ilumin "ilumin") |  call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::ilumin "ilumin"), cpro_lumin)       | CS_F_(rad_lumin)->val        | Radiative luminance
  propce(iel,\ref radiat::iqx "iqx") \n propce(iel,\ref radiat::iqy "iqy") \n propce(iel,\ref radiat::iqz "iqz")            | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::iqx "iqx"), cpro_qx) \n call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::iqy "iqy"), cpro_qy) \n call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::iqz "iqz"), cpro_qz) | CS_F_(rad_q)->val | Radiative flux
  propce(iel,\ref radiat::itsre "itsre"(iesp))  | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::itsre "itsre"(iesp)), cpro_tsre) | \ref CS_FI_(rad_ets,iesp-1)->val  | Radiative flux explicit source term
  propce(iel,\ref radiat::itsri "itsri"(iesp))  | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::itsri "itsri"(iesp)), cpro_tsri) | \ref CS_FI_(rad_its,iesp-1)->val  | Radiative flux implicit source term
  propce(iel,\ref radiat::iabso "iabso"(iesp))  | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::iabso "iabso"(iesp)), cpro_abso) | \ref CS_FI_(rad_abs,iesp-1)->val  | Radiative absorption
  propce(iel,\ref radiat::iemi "iemi"(iesp))    | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::iemi "iemi"(iesp)), cpro_emi)    | \ref CS_FI_(rad_emi,iesp-1)->val  | Radiative emission
  propce(iel,\ref radiat::icak "icak"(iesp))    | call \ref field::field_get_val_s "field_get_val_s"(iprpfl(\ref radiat::icak "icak"(iesp)), cpro_cak)    | \ref CS_FI_(rad_cak,iesp-1)->val  | Radiative absorption coefficient
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::iqinci "iqinci", bqinci)                | CS_F_(qinci)->val            | Radiative incident radiative flux density
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::ixlam "ixlam", bxlam)                   | CS_F_(xlam)->val             | Wall thermal conductivity
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::iepa "iepa", bepa)                      | CS_F_(epa)->val              | Wall thickness
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::ieps "ieps", beps)                      | CS_F_(emissivity)->val       | Wall emissivity
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::ifnet "ifnet", bfnet)                   | CS_F_(fnet)->val             | Boundary radiative flux
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::ifconv "ifconv", bfconv)                | CS_F_(fconv)->val            | Boundary radiative convective flux
  -                                             | call \ref field::field_get_val_s "field_get_val_s"(\ref radiat::ihconv "ihconv", bhconv)                | CS_F_(hconv)->val            | Radiative exchange coefficient </tt>

*/
// _____________________________________________________________________________
/*!

  \page local How to name common local variables ?

  The following table provides a non-exhaustive list of local variables which
  are used in the code in a recurring manner.

  Fortran code  | C code     | Description
  ------------- | ---------- | ------------
  <tt> iel      | cell_id    | Cell index
  ifac          | face_id    | Face index
  ig            | g_id       | Interior face number of associated groups (OpenMP)
  it            | t_id       | Interior face number of threads (OpenMP)
  idimtr        | tr_dim     | Indicator for tensor perodicity of rotation
  flumas        | i_massflux | Mass flux at interior faces
  flumab        | b_massflux | Mass flux at boundary faces
  viscf         | i_visc     | \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$ \n  at interior faces for the r.h.s.
  viscb         | b_visc     | \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$ \n  at border faces for the r.h.s.
  smbrp         | rhs        | Right hand side \f$ \vect{Rhs} \f$ </tt>

  \section cs_var_dico_conv Local naming convention for fields (Fotran and C)

  Rules have been stablished for local names denoting fields, depending on their nature. The convention, applying both in Fortran and in C, is as follows:

  - The first letter of the name indicates the location at which the field values are defined:
    - \b c for values at the cell centers.
    - \b i for values at the interior faces.
    - \b b for values at the boundary faces.
  - The next three letters indicate if the field is a variable (at the current time step or the previous time step) or a property:
    - \b var for variables at the current time step.
    - \b vara for variables at the previous time step.
    - \b pro for properties.
  - An underscore \b _ follows.
  - Finally, the <b> short name </b> of the variable/property is specified. This short name is built from the variable/property Fortran index, removing the \c i at the beginning of the word.

  The following examples ilustrate this convention:

  \c cvar_pr: Values of the variable pressure field defined at the cell centers, at the current time step. \n
  \c cvara_pr: Values of the variable pressure field defined at the cell centers, at the previous time step. \n
  \c cpro_cp: Values of the property specific heat defined field at the cell centers. \n

*/
