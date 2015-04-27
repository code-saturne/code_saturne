/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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
  \page cs_var_dico Fortran-C naming reference

  \section cs_var_dico_intro

  This page provides cross-reference tables containing the names of Fortran
  variables and their C counterparts.

  Note: Variables with the same name in Fortran and in C are not listed.

  \section cs_var_dico_namingcontent Contents

  Cross-reference tables are available for the following categories of
  variables:

   - \subpage mesh
   - \subpage field
   - \subpage local

*/
// _____________________________________________________________________________
/*!

  \page mesh Mesh and mesh quantities variables

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

  \page field Variables, properties and fields (deprecated \c propce)

  - In Fortran:
    - \ref cs_var_dico_props "Properties" are accessed through the
       deprecated array \c propce as: \n
      <tt>propce(iel,index)</tt>.
    - The scalar indexes are accessed as \ref isca "isca"(j), and the scalar values as
      follows: \n
      <tt>isca(iscalt)</tt>.
    - Both variables and properties can be accessed via the
      \ref field.f90 "cs_field" API, as in the following examples: \n
      <tt>call \ref field::field_get_val_s "field_get_val_s"(ivarfl(ipr), cvar_pr)
          \n cvar_pr(iel)</tt>, \n\n
      <tt>call \ref field::field_get_val_v "field_get_val_v"(ivarfl(iu), cvar_vel)
          \n cvar_vel(isou,iel)</tt>, \n\n
      <tt>call \ref field::field_get_val_s "field_get_val_s"(iprpfl(icp), cpro_cp) \n
          cpro_cp(iel)</tt>, \n\n
       where \ref numvar::ipr "ipr", \ref numvar::iu "iu" are variable indexes
       and \ref cstphy::icp "icp" is a property index.
  - In C:
    - Both variables and properties are accessed as: \n
      <tt>CS_F_(name)->val[cell_id]</tt>, \n
      where the \c name is defined in \ref cs_field_pointer.h. \n
    - Indexed variables (such as user scalars) and indexed properties
      are accessed as: \n
      <tt>CS_FI_(name,ii-1)->val[cell_id]</tt>.

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

  V3.0 and older Fortran code     | Current Fortran code                             | C code                       | Description
  ------------------------------- | ------------------------------------------------ | ---------------------------- | ------------
  <tt> rtp(iel,ipr)               | call field_get_val_s(ivarfl(\ref ipr), cvar_pr)       | CS_F_(p)->val[cell_id]       | Pressure
  <EM>(WARNING: removed)</EM> \n rtp(iel,iu) \n rtp(iel,iv) \n rtp(iel,iw) | call field_get_val_v(ivarfl(\ref iu), cvar_vel) \n cvar_vel(1, iel) \n cvar_vel(2, iel) \n cvar_vel(3, iel) | \n CS_F_(u)->val[cell_id][0] \n CS_F_(u)->val[cell_id][1] \n CS_F_(u)->val[cell_id][2] | Velocity
  rtp(iel,ivoidf)                 | call field_get_val_s(ivarfl(\ref ivoidf), cvar_voidf) | CS_F_(void_f)->val[cell_id]  | Void fraction for cavitation modelling
  rtp(iel,ik)                     | call field_get_val_s(ivarfl(\ref ik  ), cvar_k  )     | CS_F_(k)->val[cell_id]       | Turbulent kinetic energy \f$ k \f$
  rtp(iel,iep)                    | call field_get_val_s(ivarfl(\ref iep ), cvar_eps)     | CS_F_(eps)->val[cell_id]     | Turbulent dissipation \f$ \varepsilon \f$
  rtp(iel,ir11)                   | call field_get_val_s(ivarfl(\ref ir11), cvar_r11)     | CS_F_(r11)->val[cell_id]     | Reynolds stress component \f$ R_{xx} \f$
  rtp(iel,ir22)                   | call field_get_val_s(ivarfl(\ref ir22), cvar_r22)     | CS_F_(r22)->val[cell_id]     | Reynolds stress component \f$ R_{yy} \f$
  rtp(iel,ir33)                   | call field_get_val_s(ivarfl(\ref ir33), cvar_r33)     | CS_F_(r33)->val[cell_id]     | Reynolds stress component \f$ R_{zz} \f$
  rtp(iel,ir12)                   | call field_get_val_s(ivarfl(\ref ir12), cvar_r12)     | CS_F_(r12)->val[cell_id]     | Reynolds stress component \f$ R_{xy} \f$
  rtp(iel,ir23)                   | call field_get_val_s(ivarfl(\ref ir23), cvar_r23)     | CS_F_(r23)->val[cell_id]     | Reynolds stress component \f$ R_{yz} \f$
  rtp(iel,ir13)                   | call field_get_val_s(ivarfl(\ref ir13), cvar_r13)     | CS_F_(r13)->val[cell_id]     | Reynolds stress component \f$ R_{xz} \f$
  rtp(iel,iphi)                   | call field_get_val_s(ivarfl(\ref iphi), cvar_phi)     | CS_F_(phi)->val[cell_id]     | \f$ \phi \f$ for \f$ \phi-f_b \f$ model
  rtp(iel,ifb)                    | call field_get_val_s(ivarfl(\ref ifb ), cvar_fb )     | CS_F_(f_bar)->val[cell_id]   | \f$ f_b \f$ for \f$ \phi-f_b \f$ model
  rtp(iel,ial)                    | call field_get_val_s(ivarfl(\ref ial ), cvar_al )     | CS_F_(alpha)->val[cell_id]   | \f$ \alpha \f$ for \f$ Bl-v^2-k \f$ \n or EBRSM model
  rtp(iel,iomg)                   | call field_get_val_s(ivarfl(\ref iomg), cvar_omg)     | CS_F_(omg)->val[cell_id]     | \f$ \omega \f$ for \f$ k-\omega \f$ SST model
  rtp(iel,inusa)                  | call field_get_val_s(ivarfl(\ref inusa), cvar_nusa)   | CS_F_(nusa)->val[cell_id]    | \f$ \widetilde{\nu}_T \f$ for Spalart-Allmaras
  rtp(iel,iuma) \n rtp(iel,ivma) \n rtp(iel,iwma) | call field_get_val_v(ivarfl(\ref iuma), cvar_mesh_v)  | CS_F_(mesh_u)->val[3*cell_id] \n CS_F_(mesh_u)->val[3*cell_id+1] \n CS_F_(mesh_u)->val[3*cell_id+2] | Mesh velocity
  rtp(iel,ihm)                    | call field_get_val_s(ivarfl(isca(\ref ppincl::ihm "ihm")), cvar_hm) | CS_F_(h)->val[cell_id] | Enthalpy
  rtp(iel,isca(iscalt)) | call field_get_val_s(ivarfl(isca(\ref optcal::iscalt "iscalt")), cvar_scalt) | CS_F_(t)->val[cell_id] | Temperature </tt>


  \section cs_var_dico_props Properties

  These properties are defined in the files \ref numvar.f90 and
  \ref cs_field_pointer.h.

  Deprecated Fortran code  | Current Fortran code                                                      |C code                             | Description
  ------------------------ | ------------------------------------------------------------------------- |---------------------------------- | ------------
  <tt> dt                  | dt                                                                        |CS_F_(dt)->val[cell_id]            | Local time step
  propce(iel,irom)         | call field_get_val_s(iprpfl(\ref numvar::irom "irom"), cpro_rho  )        |CS_F_(rho)->val[cell_id]           | Density at the current time step
  propce(iel,irom)         | call field_get_val_s(iprpfl(\ref numvar::iroma "iroma"), cpro_rhoa )      |CS_F_(rho)->val_pre[cell_id]       | Density at the previous time step
  propce(iel,iromaa)       | call field_get_val_s(iprpfl(\ref numvar::iromaa "iromaa"), cpro_rhoaa)    |cs_field_by_name("density_old")    | Density at the second previous time
  propce(iel,iviscl)       | call field_get_val_s(iprpfl(\ref numvar::iviscl "iviscl"), cpro_viscl)    |CS_F_(mu)->val[cell_id]            | Molecular viscosity
  propce(iel,ivisct)       | call field_get_val_s(iprpfl(\ref numvar::ivisct "ivisct"), cpro_visct)    |CS_F_(mu_t)->val[cell_id]          | Turbulent dynamic viscosity
  propce(iel,ivisla)       | call field_get_val_s(iprpfl(\ref numvar::ivisla "ivisla"), cpro_romaa)    |CS_F_(mu)->val_pre[cell_id]        | Dynamic molecular viscosity (in kg/(m.s)) \n at the previous time-step
  propce(iel,ivista)       | call field_get_val_s(iprpfl(\ref numvar::ivista "ivista"), cpro_viscta)   |CS_F_(mu_t)->val_pre[cell_id]      | Dynamic turbulent viscosity \n at the previous time-step
  propce(iel,icp)          | call field_get_val_s(iprpfl(\ref cstphy::icp "icp"), cpro_cp)             |CS_F_(cp)->val[cell_id]            | Specific heat
  propce(iel,icpa)         | call field_get_val_s(iprpfl(\ref numvar::icpa "icpa"), cpro_cpa)          |CS_F_(cp)->val_pre[cell_id]        | specific heat at the previous time-step
  propce(iel,icrom)        | call field_get_val_s(iprpfl(\ref numvar::icrom "icrom"), cpro_crom)       |CS_F_(rho)->val[cell_id]           | Density (at cells)
  propce(iel,ibrom)        | call field_get_val_s(iprpfl(\ref numvar::ibrom "ibrom"), bpro_rho)        |CS_F_(rho_b)->val[cell_id]         | Density (at boundary faces)
  propce(iel,ismago)       | call field_get_val_s(iprpfl(\ref numvar::ismago "ismago"), cpro_smago)    |cs_field_by_name("smagorinsky_constant^2")| Field id of the anisotropic turbulent viscosity
  propce(iel,icour)        | call field_get_val_s(iprpfl(\ref numvar::icour "icour"), cpro_cour)       |cs_field_by_name("courant_number") | Courant number
  propce(iel,ifour)        | call field_get_val_s(iprpfl(\ref numvar::ifour "ifour"), cpro_four)       |cs_field_by_name("fourier_number") | Fourier number
  propce(iel,iprtot)       | call field_get_val_s(iprpfl(\ref numvar::iprtot "iprtot"), cpro_prtot)    |cs_field_by_name("total_pressure") | Total pressure at cell centers
  propce(iel,ivisma(1)) \n propce(iel,ivisma(2)) \n propce(iel,ivisma(3)) | call field_get_val_s(iprpfl(\ref numvar::ivisma "ivisma"(1)), cpro_vism1) \n call field_get_val_s(iprpfl(\ref numvar::ivisma "ivisma"(2)), cpro_vism2) \n call field_get_val_s(iprpfl(\ref numvar::ivisma "ivisma"(3)), cpro_vism3) | cs_field_by_name("mesh_viscosity_1") \n cs_field_by_name("mesh_viscosity_2") \n cs_field_by_name("mesh_viscosity_3") | Mesh velocity viscosity for the ALE module
  propce(iel,itsrho)       | call field_get_val_s(iprpfl(\ref numvar::itsrho "itsrho"), cpro_tsrho )   |cs_field_by_name("dila_st")        | Global dilatation source terms
  propce(iel,ibeta)        | call field_get_val_s(iprpfl(\ref numvar::ibeta "ibeta"), cpro_beta  )     |cs_field_by_name("thermal_expansion") | Thermal expansion coefficient
  -                        | call field_get_val_s(\ref numvar::ipori "ipori", cpro_ipori)              |CS_F_(poro)->val[cell_id]          | Porosity
  -                        | call field_get_val_v(\ref numvar::iporf "iporf", cpro_iporf)              |CS_F_(t_poro)->val[cell_id][ii]    | Tensorial porosity
  -                        | call field_get_val_v(\ref numvar::iforbr "iforbr", bpro_forbr)            |cs_field_by_name("boundary_forces")| Field id of the stresses at boundary
  -                        | call field_get_val_s(\ref numvar::iyplbr "iyplbr", bpro_yplus)            |cs_field_by_name("yplus")          | Field id of \f$y^+\f$ at boundary
  -                        | call field_get_val_v(\ref numvar::idtten "idtten", dttens)                |cs_field_by_name("dttens")         | Field id for the dttens tensor</tt>


  \section cs_var_dico_part Specific physics

  \subsection cs_var_dico_atmo Atmospheric

  Defined in \ref optcal.f90, \ref atincl.f90, \ref atvarp.f90 and
  \ref cs_field_pointer.h.

  Fortran code                                                                               | C code                                 | Description
  ------------------------------------------------------------------------------------------ | -------------------------------------- | ------------
  <tt> call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)                                | CS_F_(pot_t)->val[cell_id]             | Potential temperature
  call field_get_val_s(ivarfl(\ref isca(\ref atincl::itotwt "itotwt")), cvar_totwt)          | CS_F_(totwt)->val[cell_id]             | Total water content
  call field_get_val_s(ivarfl(\ref isca(\ref atincl::intdrp "intdrp")), cvar_intdrp)         | CS_F_(ntdrp)->val[cell_id]             | Total number of droplets
  call field_get_val_s(ivarfl(\ref isca(\ref atchem::isca_chem "isca_chem"(iesp))), cvar_sc) | CS_FI_(chemistry,iesp-1)->val[cell_id] | Chemistry species (indexed) </tt>


  \subsection cs_var_dico_comb Coal combustion

  Defined in \ref ppincl.f90, \ref ppcpfu.f90 and \ref cs_field_pointer.h.

  Fortran code                                                                    | C code                           | Description
  ------------------------------------------------------------------------------- | -------------------------------- | ------------
  <tt> call field_get_val_s(\ref isca(\ref ppincl::inp "inp"(iesp)), cvar_inpcl)  | CS_FI_(np,iesp-1)->val[cell_id]  | Particles per kg for coal class
  call field_get_val_s(\ref isca(\ref ppincl::ixch "ixch"(iesp)), cvar_xchcl)     | CS_FI_(xch,iesp-1)->val[cell_id] | Reactive coal mass fraction for coal class
  call field_get_val_s(\ref isca(\ref ppincl::ixck "ixck"(iesp)), cvar_xckcl)     | CS_FI_(xck,iesp-1)->val[cell_id] | Coke mass fraction for coal class
  call field_get_val_s(\ref isca(\ref ppincl::ixwt "ixwt"(iesp)), cvar_xwtcl)     | CS_FI_(xwt,iesp-1)->val[cell_id] | Water mass fraction for coal class
  call field_get_val_s(\ref isca(\ref ppincl::ih2 "ih2"(iesp)), cvar_h2cl)        | CS_FI_(h2,iesp-1)->val[cell_id]  | Mass enthalpy for coal class (permeatic case)
  call field_get_val_s(\ref isca(\ref ppincl::if1m "if1m"(iesp)), cvar_f1mcl)     | CS_FI_(f1m,iesp-1)->val[cell_id] | Mean value light volatiles for coal class
  call field_get_val_s(\ref isca(\ref ppincl::if2m "if2m"(iesp)), cvar_f2mcl)     | CS_FI_(f2m,iesp-1)->val[cell_id] | Mean value heavy volatiles for coal class
  call field_get_val_s(\ref isca(\ref ppincl::if4m "if4m"), cvar_f4m)             | CS_F_(f4m)->val[cell_id]         | Oxydant 2 mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::if5m "if5m"), cvar_f5m))            | CS_F_(f5m)->val[cell_id]         | Oxydant 3 mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::if6m "if6m"), cvar_f6m))            | CS_F_(f6m)->val[cell_id]         | Water from coal drying mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::if7m "if7m"), cvar_f7m))            | CS_F_(f7m)->val[cell_id]         | Carbon from coal oxidyzed by O2 mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::if8m "if8m"), cvar_f8m))            | CS_F_(f8m)->val[cell_id]         | Carbon from coal gasified by CO2 mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::if9m "if9m"), cvar_f9m))            | CS_F_(f9m)->val[cell_id]         | Carbon from coal gasified by H2O mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::ifvp2m "ifvp2m"), cvar_fvp2m)       | CS_F_(fvp2m)->val[cell_id]       | f1f2 variance
  call field_get_val_s(\ref isca(\ref ppcpfu::iyco2 "iyco2"), cvar_yco2)          | CS_F_(yco2)->val[cell_id]        | CO2 fraction
  call field_get_val_s(\ref isca(\ref ppcpfu::iyhcn "iyhcn"), cvar_yhnc)          | CS_F_(yhcn)->val[cell_id]        | HCN fraction
  call field_get_val_s(\ref isca(\ref ppcpfu::iyno "iyno"), cvar, yno)            | CS_F_(yno)->val[cell_id]         | NO fraction
  call field_get_val_s(\ref isca(\ref ppcpfu::iynh3 "iynh3"), cvar_ynh3)          | CS_F_(ynh3)->val[cell_id]        | NH3 enthalpy
  call field_get_val_s(\ref isca(\ref ppcpfu::ihox "ihox"), cvar_hox)             | CS_F_(hox)->val[cell_id]         | Ox enthalpy </tt>


  \subsection cs_var_dico_cfbl Compressible

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                                                                   | C code                        | Description
  ------------------------------------------------------------------------------ | ----------------------------- | ------------
  <tt> call field_get_val_s(\ref isca(\ref ppincl::ienerg "ienerg"), cvar_energ) | CS_F_(energy)->val[cell_id]   | Total energy
  call field_get_val_s(\ref isca(\ref ppincl::itempk "itempk"), cvar_tempk)      | CS_F_(t_kelvin)->val[cell_id] | Temperature, in Kelvin </tt>


  \subsection cs_var_dico_elec Electric arcs

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                                                                  | C code                             | Description
  ----------------------------------------------------------------------------- | ---------------------------------- | ------------
  <tt> call field_get_val_s(\ref isca(\ref ppincl::ipotr "ipotr"), cvar_potr)   | CS_F_(potr)->val[cell_id]          | Electric potential, real part
  call field_get_val_s(\ref isca(\ref ppincl::ipoti "ipoti"), cvar_poti)        | CS_F_(poti)->val[cell_id]          | Electric potential, imaginary part
  call field_get_val_s(\ref isca(\ref ppincl::ipotva "ipotva(1)"), cvar_potva1) \n call field_get_val_s(\ref isca(\ref ppincl::ipotva "ipotva(2")), cvar_potva2) \n call field_get_val_s(\ref isca(\ref ppincl::ipotva "ipotva(3)"), cvar_potva3) | CS_F_(potva)->val[cell_id][0] \n CS_F_(potva)->val[cell_id][1] \n CS_F_(potva)->val[cell_id][2] | Vector potential
  call field_get_val_s(\ref isca(\ref ppincl::iycoel "iycoel"(iesp)), cvar_ycoel(iesp)) | CS_FI_(ycoel,iesp-1)->val[cell_id] | Constituent mass fraction </tt>


  \subsection cs_var_dico_cogz Gas combustion

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                                                              | C code                     | Description
  ------------------------------------------------------------------------- | -------------------------- | ------------
  <tt> call field_get_val_s(\ref isca(\ref ppincl::ifm), cvar_fm)           | CS_F_(fm)->val[cell_id]    | Mixture fraction
  call field_get_val_s(\ref isca(\ref ppincl::ifp2m "ifp2m"), cvar_fp2m)    | CS_F_(fp2m)->val[cell_id]  | Mixture fraction variance
  call field_get_val_s(\ref isca(\ref ppincl::ifsm "ifsm"), cvar_fsm)       | CS_F_(fsm)->val[cell_id]   | Soot mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::inpm "inpm"), cvar_npm)       | CS_F_(npm)->val[cell_id]   | Soot precursor number
  call field_get_val_s(\ref isca(\ref ppincl::iygfm "iygfm"), cvar_ygfm)    | CS_F_(ygfm)->val[cell_id]  | Fresh gas fraction
  call field_get_val_s(\ref isca(\ref ppincl::iyfm "iyfm"), cvar_yfm)       | CS_F_(yfm)->val[cell_id]   | Mass fraction
  call field_get_val_s(\ref isca(\ref ppincl::iyfp2m "iyfp2m"), cvar_yfp2m) | CS_F_(yfp2m)->val[cell_id] | Mass fraction variance
  call field_get_val_s(\ref isca(\ref ppincl::icoyfp "icoyfp"), cvar_coyfp) | CS_F_(coyfp)->val[cell_id] | Mass fraction covariance </tt>


  \subsection cs_var_dico_rayt Radiative transfer

  Defined in \ref radiat.f90 and \ref cs_field_pointer.h.

  Fortran code                                                              | C code                               | Description
  ------------------------------------------------------------------------- | ------------------------------------ | ------------
  <tt> propce(iel,\ref radiat::ilumin "ilumin")                             | CS_F_(rad_lumin)->val[cell_id]       | Radiative luminance
  propce(iel,\ref radiat::iqx "iqx") \n propce(iel,\ref radiat::iqy "iqy") \n propce(iel,\ref radiat::iqz "iqz") | CS_F_(rad_q)->val[cell_id][0] \n CS_F_(rad_q)->val[cell_id][1] \n CS_F_(rad_q)->val[cell_id][2] | Radiative flux
  propce(iel,\ref radiat::itsre "itsre"(iesp))                              | CS_FI_(rad_ets,iesp-1)->val[cell_id] | Radiative flux explicit source term
  propce(iel,\ref radiat::itsri "itsri"(iesp))                              | CS_FI_(rad_its,iesp-1)->val[cell_id] | Radiative flux implicit source term
  propce(iel,\ref radiat::iabso "iabso"(iesp))                              | CS_FI_(rad_abs,iesp-1)->val[cell_id] | Radiative absorption
  propce(iel,\ref radiat::iemi "iemi"(iesp))                                | CS_FI_(rad_emi,iesp-1)->val[cell_id] | Radiative emission
  propce(iel,\ref radiat::icak "icak"(iesp))                                | CS_FI_(rad_cak,iesp-1)->val[cell_id] | Radiative absorption coefficient
  call field_get_val_s(\ref radiat::itparo "itparo", btparo) \n btparo(iel) | CS_F_(tparo)->val[cell_id]           | Wall temperature
  call field_get_val_s(\ref radiat::iqinci "iqinci", bqinci) \n bqinci(iel) | CS_F_(qinci)->val[cell_id]           | Radiative incident radiative flux density
  call field_get_val_s(\ref radiat::ixlam "ixlam", bxlam) \n bxlam(iel)     | CS_F_(xlam)->val[cell_id]            | Wall thermal conductivity
  call field_get_val_s(\ref radiat::iepa "iepa", bepa) \n bepa(iel)         | CS_F_(epa)->val[cell_id]             | Wall thickness
  call field_get_val_s(\ref radiat::ieps "ieps", beps) \n beps(iel)         | CS_F_(emissivity)->val[cell_id]      | Wall emissivity
  call field_get_val_s(\ref radiat::ifnet "ifnet", bfnet) \n bfnet(iel)     | CS_F_(fnet)->val[cell_id]            | Boundary radiative flux
  call field_get_val_s(\ref radiat::ifconv "ifconv", bfconv) \n bfconv(iel) | CS_F_(fconv)->val[cell_id]           | Boundary radiative convective flux
  call field_get_val_s(\ref radiat::ihconv "ihconv", bhconv) \n bhconv(iel) | CS_F_(hconv)->val[cell_id]           | Radiative exchange coefficient </tt>

*/
// _____________________________________________________________________________
/*!

  \page local Common local variables

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
